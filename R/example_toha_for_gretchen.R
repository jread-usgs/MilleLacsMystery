library(glmtools)
library(mda.lakes)


bathy   = read.table('inst/extdata/ML_hypso.tsv', sep='\t', header=TRUE)
kd      = read.csv("inst/extdata/ml_model_kd.csv")
kd$time = as.POSIXct(kd$time, tz='Etc')

names(bathy) = c('depths', 'areas')
ncfile  = '/tmp/RtmpKTMTDT/ml_run/output.nc'
nmlfile  = '/tmp/RtmpKTMTDT/ml_run/glm2.nml'
io = get_var(ncfile, 'I_0')
#opt_wtr = get_temp(ncfile, z_out=bathy$depths, reference='surface')
opt_wtr = read.csv('debiased_wtr_for_toha.csv', sep=',', header=TRUE)
opt_wtr$DateTime = as.POSIXct(opt_wtr$DateTime, tz='Etc')
kd      = subset(kd, time %in% opt_wtr$DateTime)

uyears = unique(lubridate::year(opt_wtr$DateTime))
out = data.frame()
season = 4:6

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
