

estimate_toha <- function(file_out, ncfile){
  
  sim_dir <- dirname(ncfile) %>% dirname
  io <- get_var(ncfile, 'I_0')
  
  nml <- read_nml(nml_file = file.path(sim_dir, 'glm3.nml'))
  bathy <- data.frame(H = nml$morphometry$H, A = nml$morphometry$A)
  
  hypsos <- bathy %>% mutate(depths = max(H) - H, areas = A) %>% 
    arrange(depths) %>% select(depths, areas)
  
  opt_wtr <- get_temp(ncfile, z_out=hypsos$depths, reference='surface') %>% 
    mutate(DateTime = as.POSIXct(as.character(DateTime), tz='Etc')) # force GMT...
  
  kd <- readr::read_csv(file.path(sim_dir, nml$light$Kw_file), col_types = 'cd') %>% 
    mutate(time = as.POSIXct(time, tz='Etc')) %>% filter(time %in% opt_wtr$DateTime)
  
  uyears = unique(lubridate::year(opt_wtr$DateTime))
  out = list()
  
  for(i in 1:(nrow(opt_wtr))){
    
    naomitwtr = opt_wtr[i,]
    tmpwtr = naomitwtr[,!sapply(naomitwtr, is.na)]
    opti = mda.lakes::opti_thermal_habitat(tmpwtr,
                                           io[i,],  
                                           kd[i,]$Kd, lat = lat, lon = lon, hypsos, irr_thresh = c(0.0762, 0.6476), 
                                           wtr_thresh=c(11,25), interp_hour=TRUE, area_type="benthic", approx_method='constant')
    
    opti$DateTime = opt_wtr$DateTime[i]
    out[[i]] = opti
  }
  
  out = dplyr::bind_rows(out)
  
  readr::write_csv(out, path = file_out)
}