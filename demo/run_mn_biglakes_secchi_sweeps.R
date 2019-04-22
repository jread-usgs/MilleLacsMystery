
#first crack at running mn biglakes GLM
library(mda.lakes)
library(glmtools)
library(lakemodeltools)
library(lubridate)

run_ml = function(lk_name, custom_nml=NULL, kw_factor=1, Kw_file, ll, bathy, wtr, dvr_file){
  
  run_dir = file.path(tempdir(), lk_name)
  dir.create(run_dir)
  dvr = read.csv(dvr_file, as.is=TRUE)
  tmpdvr = file.path(run_dir, 'newdvr.csv')
  
  dvr_rain = mda.lakes::driver_add_rain(dvr, months = 7:9, rain_add = 1)
  write.csv(dvr_rain, tmpdvr, quote=FALSE, row.names=FALSE)
  
  
  #some default parameters from kevin's run
  nml_args=list(
    dt=3600, subdaily=FALSE, nsave=24,
    timezone=-6,
    csv_point_nlevs=0, 
    meteo_fl=gsub("\\\\", '/', tmpdvr),
    'min_layer_thick'=0.01)
  
  #bathy = read.table(paste0(getwd(), '/inst/extdata/hypsos/', lk_name, '.csv'), sep=',', header=TRUE)
  #names(bathy) = c('depths', 'areas')
  
  nml_obj = populate_base_lake_nml(site_id='ml', kd=1, driver = 'dummy.csv', zmax = max(bathy$depth),
                                   cd=coef_drag(wind_sheltering(max(bathy$area))), bathy = bathy, elev = 700, latlon = ll, 
                                   area = max(bathy$area))
  
  nml_obj = set_nml(nml_obj, arg_list=c(nml_args, custom_nml))
  
  #if we're using variable Kw file, add it to NML
  nml_obj[['glm_setup']]$Kw_file = Kw_file
  
  
  #@jread edits
  # nml_obj = set_nml(nml_obj, 'max_layer_thick', 3)
  # nml_obj = set_nml(nml_obj, 'sw_factor',1.08)
  # nml_obj = set_nml(nml_obj, 'coef_wind_stir',0.376)
  # nml_obj = set_nml(nml_obj, 'coef_mix_hyp',0.22)
  # nml_obj = set_nml(nml_obj, 'coef_mix_KH', 0.074)
  # nml_obj = set_nml(nml_obj, 'cd', 0.0011)
  # nml_obj = set_nml(nml_obj, 'ce', 0.00095)
  nml_obj = set_nml(nml_obj, 'sed_temp_mean',4.5)
  nml_obj = set_nml(nml_obj, 'sed_temp_amplitude',0.25)
  # nml_obj = set_nml(nml_obj, 'min_layer_thick', 0.1)
  # nml_obj = set_nml(nml_obj, 'coef_mix_conv', 0.25) # was 0.23, from Bruce et al is 0.20
  
  nml_obj = set_nml(nml_obj, 'start', '2008-01-01')
  nml_obj = set_nml(nml_obj, 'stop', '2010-02-01') #2009 for tohas
  nml_obj[['sed_heat']] = NULL
  
  write_nml(nml_obj, file = file.path(run_dir, 'glm2.nml'))
  
  run_glm(run_dir)
  # origin = getwd()
  # glmpath = file.path(getwd(), 'inst/extbin/glm')
  # setwd(run_dir)
  # system2(glmpath)
  # setwd(origin)
  
  if(!missing(wtr)){
    #going to assume passed df is correctly formatted
    fname = tempfile(fileext = '.csv')
    write.csv(wtr, fname, quote=FALSE, row.names=FALSE)
    output = resample_to_field(file.path(run_dir, 'output.nc'), 
                               field_file = fname)
  }else{
    output = data.frame()
  }
  
  attr(output, 'run_dir') = run_dir
  
  return(output)
}

calc_toha = function(run_dir, bathy, kd, ll){
  
  #bathy   = read.table(paste0(getwd(), '/inst/extdata/hypsos/', lk_name, '.csv'), sep=',', header=TRUE)
  #kd      = read.csv(kdfile, as.is=TRUE)
  #kd$time = as.POSIXct(kd$time, tz='Etc')
  
  names(bathy) = c('depths', 'areas')
  ncfile  = file.path(run_dir, 'output.nc')
  nmlfile  = file.path(run_dir, 'glm2.nml')
  io = get_var(ncfile, 'I_0')
  
  if(length(kd) ==1){
    kd = data.frame(time=io$DateTime, Kd=kd)
  }
  
  opt_wtr = get_temp(ncfile, z_out=bathy$depths, reference='surface')
  
  #### read this from file jordan sent
  #opt_wtr = read.csv('debiased_wtr_for_toha.csv', sep=',', header=TRUE)
  opt_wtr$DateTime = as.POSIXct(opt_wtr$DateTime, tz='Etc')
  
  kd      = subset(kd, time %in% opt_wtr$DateTime)
  
  uyears = unique(lubridate::year(opt_wtr$DateTime))
  out = list()
  #  season = 4:6
  
  library(parallel)
  c1 = parallel::makePSOCKcluster(rep('localhost', 5))
  clusterExport(c1, varlist = c('opt_wtr', 'io', 'kd', 'll', 'bathy'), envir =  environment())
  
  #for(i in 1:(nrow(opt_wtr))){
  out = clusterApplyLB(c1, 1:nrow(opt_wtr), 
                       function(i){
    #   opti = mda.lakes::opti_thermal_habitat(subset(opt_wtr, lubridate::year(DateTime) == uyears[i] & month(DateTime) %in% season), 
    #                                          subset(io, lubridate::year(DateTime) == uyears[i] & month(DateTime) %in% season), 
    #                                          subset(kd, lubridate::year(time) == uyears[i] & month(time) %in% season)$Kd, ll[1], ll[2], bathy, irr_thresh = c(0.0762, 0.6476), 
    #                                   wtr_thresh=c(11,25), interp_hour=TRUE, area_type="benthic")
    #   opti$year = uyears[i]
    naomitwtr = opt_wtr[i,]
    tmpwtr = naomitwtr[,!sapply(naomitwtr, is.na)]
    opti = mda.lakes::opti_thermal_habitat(tmpwtr,
                                           io[i,],  
                                           kd[i,]$Kd, ll[1], ll[2], bathy, irr_thresh = c(0.0762, 0.6476), 
                                           wtr_thresh=c(11,25), interp_hour=TRUE, area_type="benthic", approx_method='constant')
    
    opti$DateTime = opt_wtr$DateTime[i]
    #out[[i]] = opti
    return(opti)
  })
  
  stopCluster(c1)
  
  out = dplyr::bind_rows(out)
  
  return(out)
  #write.csv(opt_wtr, 'out/2018-01-09_wtr_for_toha.csv', row.names=FALSE)
  #write.csv(out, 'out/2018-01-09_daily_toha_estimate.csv', row.names=FALSE)
}


lakell = read.csv('inst/extdata/mn_big_lake_table.csv', as.is=TRUE, colClasses=c(DOW='character'))
lakes = read.csv('inst/extdata/Large_lakes_DOW_all.csv', as.is=TRUE, colClasses=c(DOW='character'))
lakes = merge(lakes, lakell[,c('Lake', 'Lat', 'Lon')], by='Lake')

bthfiles = data.frame(bathyfile = file.path(getwd(), Sys.glob('inst/extdata/hypsos/*')), stringsAsFactors = FALSE)
bthfiles$Lake = stringr::str_match(basename(bthfiles$bathyfile), '(.*)_([0-9]*)')[,2]
lakes = merge(lakes, bthfiles, by='Lake')

kdfiles = data.frame(kdfile = file.path(getwd(), Sys.glob('inst/extdata/kds/*')), stringsAsFactors = FALSE)
kdfiles$Lake = stringr::str_match(basename(kdfiles$kdfile), '(.*)_([0-9]*)')[,2]
lakes = merge(lakes, kdfiles, by='Lake')

dvrfiles = data.frame(dvrfile = file.path(getwd(), Sys.glob('inst/extdata/drivers/*')), stringsAsFactors = FALSE)
dvrfiles$Lake = stringr::str_match(basename(dvrfiles$dvrfile), '(.*)_([0-9]*)')[,2]
lakes = merge(lakes, dvrfiles, by='Lake')

allwtr = readRDS('inst/extdata/temps/mn_biglakes_all_temp.rds')

ulakes = unique(lakes$Lake)

tohas = list()
for(i in 1:length(ulakes)){
  #tryCatch({
  lakemeta = subset(lakes, Lake == ulakes[i])
  bathy = read.table(lakemeta$bathyfile[1], sep=',', header=TRUE)
  wtr = subset(allwtr, DOW.x %in% lakemeta$DOW) 
  wtr = wtr[, c('DateTime', 'Depth', 'temp')]
  
  secchi = seq(0.5, 10, by=0.2)
  kds = 1.7/secchi
  
  for(j in 1:length(kds)){
    out = run_ml(lakemeta$Lake[1], Kw_file=NULL, ll = as.numeric(lakemeta[1, c('Lat', 'Lon')]), 
                 bathy=bathy, dvr_file = lakemeta$dvrfile[1], custom_nml = c('Kw'=kds[j]))
    
    #plot(glmtools::get_surface_height(file.path(attr(out, 'run_dir'), 'output.nc')))
    
    cat('Running ', lakemeta$Lake[1], ' kd:', kds[j], '...\n')
    
    toha = calc_toha(attr(out, 'run_dir'), ll = as.numeric(lakemeta[1, c('Lat', 'Lon')]), 
                     bathy=bathy, kd=kds[j])
    
    toha$Lake = lakemeta$Lake[1]
    toha$secchi = secchi[j]
    
    tohas[[length(tohas)+1]] = toha
  }
  save.image('b:/test/secchi_sweep_backup.Rdata')
  
  unlink(attr(out, 'run_dir'), recursive=TRUE)
  #}, error=function(e){unlink(attr(out, 'run_dir'), recursive=TRUE)})
}

alltoha = dplyr::bind_rows(tohas)

alltoha = subset(alltoha, DateTime >= as.POSIXct('2009-01-01') & as.POSIXct(DateTime) <= '2009-12-31')
write.csv(alltoha, 'b:/test/2018-05-25_all_toha_secchi_sweeps.csv', row.names=FALSE, quote=FALSE)

