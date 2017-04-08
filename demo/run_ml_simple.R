#first crack at running ML GLM
library(mda.lakes)
library(glmtools)
library(lakemodeltools)

ll = c(46.244647, -93.652001)
# dvr = lakemodeltools::nldas_to_glm_drivers(nldastools::get_primary_forcing_local(ll[2], ll[1]))
# dvr = lakemodeltools::driver_aggregate(dvr)
# dvr = nldastools::driver_nldas_wind_debias(dvr)
# 
# write.csv(dvr, 'inst/extdata/ml_nldas_drivers.csv', quote=FALSE, row.names=FALSE)


run_ml = function(custom_nml=NULL, fixed_kd, kw_factor=1){
  
  run_dir = file.path(tempdir(), 'ml_run')
  dir.create(run_dir)
  
  #some default parameters from kevin's run
  nml_args=list(
    dt=3600, subdaily=FALSE, nsave=24,
    timezone=-6,
    csv_point_nlevs=0, 
    meteo_fl=file.path(getwd(), 'inst/extdata/ml_nldas_drivers.csv'),
    'min_layer_thick'=0.01)
  
  bathy = read.table('inst/extdata/ML_hypso.tsv', sep='\t', header=TRUE)
  names(bathy) = c('depths', 'areas')
    
  nml_obj = populate_base_lake_nml(site_id='ml', kd=1, driver = 'dummy.csv', zmax = max(bathy$depth),
                                  cd=coef_drag(wind_sheltering(max(bathy$area))), bathy = bathy, elev = 700, latlon = ll, 
                                  area = max(bathy$area))
  
  nml_obj = set_nml(nml_obj, arg_list=c(nml_args, custom_nml))
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
  
  wtemp = read.table('inst/extdata/ML_observed_temperatures.txt', sep='\t', header=TRUE)
  wtemp = subset(wtemp, temp.f > 32)
  wtemp$Depth = wtemp$depth.ft * 0.3048
  wtemp$temp  = (wtemp$temp.f - 32)*5/9
  wtemp$DateTime = lubridate::mdy(wtemp$Date)
  
  ## write wtemp obs file
  #having a weird issue with resample_to_field, make unique
  wtemp = wtemp[!duplicated(wtemp[,c('DateTime', 'Depth')]), ]
  
  write.table(wtemp[, c('DateTime', 'Depth', 'temp')], file.path(run_dir, 'obs.tsv'), sep='\t', row.names=FALSE)
  
  nml_obj = set_nml(nml_obj, 'start', '1980-04-01')
  nml_obj = set_nml(nml_obj, 'stop', '2016-10-01')
  nml_obj[['sed_heat']] = NULL
  
  if(missing(fixed_kd)){
    nml_obj = set_nml(nml_obj, 'Kw', median(0.54, na.rm=TRUE)*kw_factor)
  }else{
    nml_obj = set_nml(nml_obj, 'Kw', fixed_kd)
  }
  
  write_nml(nml_obj, file = file.path(run_dir, 'glm2.nml'))
  
  run_glm(run_dir)
  
  output = resample_to_field(file.path(run_dir, 'output.nc'), 
                                   field_file = file.path(run_dir, 'obs.tsv'))
  
  return(output)
  
}



tomod = function(param){
  
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
    
    datasp = run_ml(custom_nml=c(custom_nml, fixed_nml), kw_factor=kw_factor)
    
    datasp$error = datasp$Modeled_temp - datasp$Observed_temp
    
    par(cex=0.8, mfrow=c(2,1))
    plot(datasp$Observed_temp, datasp$Modeled_temp, xlab='Observed', ylab='Modeled')
    boxplot(error~floor(Depth), datasp, main=paste(names(param), param, collapse=' ', sep=':'), xlab=paste0('RMSE:', RMSE(datasp[,3:4])))
    abline(h=0)
    return(RMSE(datasp[,3:4]))
  #}, error = function(e){return(NA)})
  
}


initial_params = c('cd'=0.0014, 'ce'=0.0014, ch=0.0014, coef_wind_stir=0.376, coef_mix_hyp=0.22, coef_mix_conv=0.2, coef_mix_KH=0.074)

tomod(initial_params)

lower = initial_params * 0.6
upper = initial_params * 1.4
res = optim(fn = tomod, par=initial_params, control=list(parscale=c(initial_params*0.1)))

