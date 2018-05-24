
#first crack at running ML GLM
library(mda.lakes)
library(glmtools)
library(lakemodeltools)
library(lubridate)

ll = c(46.244647, -93.652001)
# dvr = lakemodeltools::nldas_to_glm_drivers(nldastools::get_primary_forcing_local(ll[2], ll[1]))
# dvr = lakemodeltools::driver_aggregate(dvr)
# dvr = nldastools::driver_nldas_wind_debias(dvr)
# 
# write.csv(dvr, 'inst/extdata/ml_nldas_drivers.csv', quote=FALSE, row.names=FALSE)


run_ml = function(lk_name, custom_nml=NULL, fixed_kd, kw_factor=1, Kw_file=NULL){
  
  run_dir = file.path(tempdir(), lk_name)
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
  
  #if we're using variable Kw file, add it to NML
  if(!is.null(Kw_file)){
    nml_obj[['glm_setup']]$Kw_file = Kw_file
  }
  
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
  
  #/run_glm(run_dir)
  origin = getwd()
  glmpath = file.path(getwd(), 'inst/extbin/glm')
  setwd(run_dir)
  system2(glmpath)
  setwd(origin)
  
  output = resample_to_field(file.path(run_dir, 'output.nc'), 
                                   field_file = file.path(run_dir, 'obs.tsv'))
  
  return(output)
  
}



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



bathy   = read.table('inst/extdata/ML_hypso.tsv', sep='\t', header=TRUE)
kd      = read.csv("inst/extdata/ml_model_kd_2017.csv")
kd$time = as.POSIXct(kd$time, tz='Etc')

names(bathy) = c('depths', 'areas')
ncfile  = './output.nc'
nmlfile  = './glm2.nml'
io = get_var(ncfile, 'I_0')


opt_wtr = get_temp(ncfile, z_out=bathy$depths, reference='surface')

#### read this from file jordan sent
#opt_wtr = read.csv('debiased_wtr_for_toha.csv', sep=',', header=TRUE)
opt_wtr$DateTime = as.POSIXct(opt_wtr$DateTime, tz='Etc')

kd      = subset(kd, time %in% opt_wtr$DateTime)

uyears = unique(lubridate::year(opt_wtr$DateTime))
out = list()
season = 4:6

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
}

out = dplyr::bind_rows(out)

write.csv(opt_wtr, 'out/2018-01-09_wtr_for_toha.csv', row.names=FALSE)
write.csv(out, 'out/2018-01-09_daily_toha_estimate.csv', row.names=FALSE)


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

