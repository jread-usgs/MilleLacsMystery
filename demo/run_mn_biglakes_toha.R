
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
    meteo_fl=tmpdvr,
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
  
  nml_obj = set_nml(nml_obj, 'start', '1980-04-01')
  nml_obj = set_nml(nml_obj, 'stop', '2016-10-01')
  nml_obj[['sed_heat']] = NULL
  
  write_nml(nml_obj, file = file.path(run_dir, 'glm2.nml'))
  
  #/run_glm(run_dir)
  origin = getwd()
  glmpath = file.path(getwd(), 'inst/extbin/glm')
  setwd(run_dir)
  system2(glmpath)
  setwd(origin)
  
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

calc_toha = function(run_dir, bathy, kdfile, ll){
  
  #bathy   = read.table(paste0(getwd(), '/inst/extdata/hypsos/', lk_name, '.csv'), sep=',', header=TRUE)
  kd      = read.csv(kdfile, as.is=TRUE)
  kd$time = as.POSIXct(kd$time, tz='Etc')
  
  names(bathy) = c('depths', 'areas')
  ncfile  = file.path(run_dir, 'output.nc')
  nmlfile  = file.path(run_dir, 'glm2.nml')
  io = get_var(ncfile, 'I_0')
  
  
  opt_wtr = get_temp(ncfile, z_out=bathy$depths, reference='surface')
  
  #### read this from file jordan sent
  #opt_wtr = read.csv('debiased_wtr_for_toha.csv', sep=',', header=TRUE)
  opt_wtr$DateTime = as.POSIXct(opt_wtr$DateTime, tz='Etc')
  
  kd      = subset(kd, time %in% opt_wtr$DateTime)
  
  uyears = unique(lubridate::year(opt_wtr$DateTime))
  out = list()
#  season = 4:6
  
  for(i in 1:(nrow(opt_wtr))){
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
    out[[i]] = opti
    if(i %% 100 == 0){
      cat(sprintf('%0.2g%% done\n', 100*i/nrow(opt_wtr)))
    }
  }
  
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
  
  out = run_ml(lakemeta$Lake[1], Kw_file=lakemeta$kdfile[1], ll = as.numeric(lakemeta[1, c('Lat', 'Lon')]), 
               bathy=bathy, dvr_file = lakemeta$dvrfile[1])
  
  plot(glmtools::get_surface_height(file.path(attr(out, 'run_dir'), 'output.nc')))
  
  cat('Running ', lakemeta$Lake[1], '...\n')
  
  toha = calc_toha(attr(out, 'run_dir'), ll = as.numeric(lakemeta[1, c('Lat', 'Lon')]), 
                   bathy=bathy, kdfile=lakemeta$kdfile[1])
  
  toha$Lake = lakemeta$Lake[1]
  
  tohas[[i]] = toha
  unlink(attr(out, 'run_dir'), recursive=TRUE)
  #}, error=function(e){unlink(attr(out, 'run_dir'), recursive=TRUE)})
}

alltoha = dplyr::bind_rows(tohas)


tomod = function(param, Kw_file=NULL){
  
  if(any('kw_factor' %in% names(param))){
    kw_factor = param['kw_factor']
    param = param[!names(param) %in% 'kw_factor'] 
  }else{
    kw_factor = 1
  }
  
  
  custom_nml=as.list(param) #list('coef_wind_stir'=param[1], coef_mix_hyp=param[2], coef_mix_conv=param[3], 'cd'=param[4], 'ce'=param[5], 'coef_mix_KH'=param[6])
  fixed_nml = as.list(c(coef_mix_shear=0.2, coef_mix_turb=0.51, 'coef_wind_stir'=0.376, 'coef_mix_hyp'=0.22, 'coef_mix_KH'=0.074, 
                        'coef_mix_conv'=0.25))
  
  #tryCatch({
    
    datasp = run_ml(custom_nml=c(custom_nml, fixed_nml), kw_factor=kw_factor, Kw_file=Kw_file)
    
    
    datasp$error = datasp$Modeled_temp - datasp$Observed_temp
    
    par(cex=0.8, mfrow=c(2,1))
    plot(datasp$Observed_temp, datasp$Modeled_temp, xlab='Observed', ylab='Modeled')
    boxplot(error~floor(Depth), datasp, main=paste(names(param), param, collapse=' ', sep=':'), xlab=paste0('RMSE:', RMSE(datasp[,3:4])))
    abline(h=0)
    browser()
    return(RMSE(datasp[,3:4]))
  #}, error = function(e){return(NA)})
  
}

##### Run and calc TOHA #####

initial_params = c('cd'=0.0014, 'ce'=0.0014, ch=0.0014, coef_wind_stir=0.376, 
                   coef_mix_hyp=0.22, coef_mix_conv=0.2, coef_mix_KH=0.074)



tomod(initial_params, Kw_file=file.path(getwd(), 'inst/extdata/ml_model_kd_2017.csv'))



plot(out$DateTime, out$opti_therm_hab, type='l', col='green')
abline(v=as.POSIXct('1988-01-01'));abline(v=as.POSIXct('1993-01-01'));abline(v=as.POSIXct('2013-01-01'))

out$year = year(out$DateTime)
out = plyr::ddply(out, 'year', plyr::summarize, opti_therm_hab=sum(opti_therm_hab), therm_hab=sum(therm_hab), opti_hab=sum(opti_hab))


png('out/figures/toha_annual.png', res=450, width=2100, height=3000)
par(mfrow = c(3,1), mar=c(1,4,1,1), oma=c(3,0,0,0))
plot(out$year, out$opti_therm_hab, ylab='Thermo-Optical Habitat', xaxt='n', type='o', pch=16)
abline(v=1988);abline(v=1993);abline(v=2013)
plot(out$year, out$therm_hab, ylab='Thermal Habitat', xaxt='n', type='o', pch=16)
abline(v=1988);abline(v=1993);abline(v=2013)
plot(out$year, out$opti_hab, ylab='Optical Habitat', xlab='year', type='o', pch=16)
abline(v=1988);abline(v=1993);abline(v=2013)
dev.off()




##### Run Optim on Model #####
lower = initial_params * 0.6
upper = initial_params * 1.4

initial_params = c('cd'=0.0014, 'ce'=0.0014, ch=0.0014, coef_wind_stir=0.376, 
                   coef_mix_conv=0.2, sw_factor=1, at_factor=1, rh_factor=1)

res = optim(fn = tomod, par=initial_params, control=list(parscale=c(initial_params*0.1)), 
            Kw_file=file.path(getwd(), 'inst/extdata/ml_model_kd.csv'))

