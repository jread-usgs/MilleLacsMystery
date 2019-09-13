library(glmtools)

generate_secchi = function(min, years){
  max = 1.569*min
  sdt = data.frame(yday = seq(100, 320), Secchi = customLogitModel(x=seq(100,320), max=max, min=min))
  sdt = rbind(data.frame(yday=c(1, 365), Secchi=c(min, min)), sdt)
  sdt = as.data.frame(approx(sdt$yday, sdt$Secchi, 1:365))
  names(sdt) = c('yday', 'Secchi')
  
  sdt = do.call(rbind, lapply(years, function(x){sdt$year = x;return(sdt)}))
  sdt$time = as.Date(ISOdate(sdt$year, 1, 1)) - 1 + sdt$yday
  sdt$Kd = 1.7/sdt$Secchi
  return(sdt[, c('time', 'Kd')])
}

run_mod = function(min){
  nml = read_nml('glm2.nml')
  
  #' Generate and save kd
  runpath = tempdir()
  tmpkd = file.path(runpath, 'tmpkd.csv')
  gs = generate_secchi(min, 1997:1998)
  write.table(gs, tmpkd, sep=',', quote=FALSE, row.names=FALSE)
  
  
  #' copy driver file over and run against same year (1997-1998)
  tmpdvr = file.path(runpath, 'tmpdvr.csv')
  file.copy('inst/extdata/ml_nldas_drivers.csv', tmpdvr, overwrite = TRUE)

  nml = set_nml(nml, 'meteo_fl', 'tmpdvr.csv')
  nml = set_nml(nml, 'start', '1997-01-01')
  nml = set_nml(nml, 'stop', '1998-12-31')
  nml = set_nml(nml, 'Kw_file', tmpkd)
  
  glmtools::write_nml(nml, file.path(runpath, 'glm2.nml'))
  
  origin = getwd()
  glmpath = file.path(getwd(), 'inst/extbin/glm')
  setwd(runpath)
  system2(glmpath)
  setwd(origin)
  
  return(list(ncfile=file.path(runpath, 'output.nc'), kdfile=tmpkd))
}

calc_toha = function(ncfile, kdfile, myear=1998){
  
  bathy   = read.table('inst/extdata/ML_hypso.tsv', sep='\t', header=TRUE)
  kd      = read.csv(kdfile)
  kd$time = as.POSIXct(kd$time, tz='Etc')
  
  names(bathy) = c('depths', 'areas')
  io = get_var(ncfile, 'I_0')
  io = subset(io, lubridate::year(DateTime) == myear)
  opt_wtr = get_temp(ncfile, z_out=bathy$depths, reference='surface')
  opt_wtr = subset(opt_wtr, lubridate::year(DateTime) == myear)
  
  opt_wtr$DateTime = as.POSIXct(opt_wtr$DateTime, tz='Etc')
  kd      = subset(kd, time %in% opt_wtr$DateTime)
  
  uyears = unique(lubridate::year(opt_wtr$DateTime))
  out = data.frame()
  
  for(i in 1:(nrow(opt_wtr)-1)){
    #   opti = mda.lakes::opti_thermal_habitat(subset(opt_wtr, lubridate::year(DateTime) == uyears[i] & month(DateTime) %in% season), 
    #                                          subset(io, lubridate::year(DateTime) == uyears[i] & month(DateTime) %in% season), 
    #                                          subset(kd, lubridate::year(time) == uyears[i] & month(time) %in% season)$Kd, ll[1], ll[2], bathy, irr_thresh = c(0.0762, 0.6476), 
    #                                   wtr_thresh=c(11,25), interp_hour=TRUE, area_type="benthic")
    #   opti$year = uyears[i]
    opti = mda.lakes::opti_thermal_habitat(opt_wtr[i:(i+1),],
                                           io[i:(i+1),],  
                                           kd[i:(i+1),]$Kd, ll[1], ll[2], bathy, irr_thresh = c(0.0762, 0.6476), 
                                           wtr_thresh=c(11,25), interp_hour=TRUE, area_type="benthic")
    
    opti$DateTime = opt_wtr$DateTime[i]
    out = rbind(opti, out)
  }
  return(out)
}

ll = c(46.244647, -93.652001)
msecchis = seq(0.5, 8, by=0.1)
out = list()
for(i in seq_along(msecchis)){
  run_info = run_mod(msecchis[i])
  out [[i]] = calc_toha(run_info$ncfile, run_info$kdfile)
  out[[i]]$min_secchi = msecchis[i]
}

tohas = do.call(rbind, out)
write.table(tohas, 'out/toha_gradient.tsv', sep='\t', row.names=FALSE)
