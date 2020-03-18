build_nml_list <- function(nml_ids, H_A_file, cd_file, lat_lon_file, len_wid_file, lake_depth_file, layer_thick_file, meteo_fl_file){
  
  
  nml_df_data <- readRDS(cd_file) %>% 
    inner_join(readRDS(lat_lon_file), by = 'site_id') %>% 
    inner_join(readRDS(len_wid_file), by = 'site_id') %>% 
    inner_join(readRDS(lake_depth_file), by = 'site_id') %>% 
    inner_join(readRDS(layer_thick_file), by = 'site_id') %>% 
    inner_join(readRDS(meteo_fl_file), by = 'site_id')
  
  # now note the H_A_file is a list()
  nml_list <- split(nml_df_data, seq(nrow(nml_df_data))) %>% setNames(nml_df_data$site_id)
  
  H_A_list <- readRDS(H_A_file)
  
  H_A_site_ids <- names(H_A_list)
  
  for (id in H_A_site_ids[H_A_site_ids %in% names(nml_list)]){
    nml_list[[id]] = as.list(nml_list[[id]])
    nml_list[[id]]$A = H_A_list[[id]]$A
    nml_list[[id]]$H = H_A_list[[id]]$H
  }
  
  return(nml_list)
}

Kw_file_dir <- '~/Downloads/USGS_Kd_Predictions_byNHDHR'

get_ids <- function(H_A_file, xwalk_file, meteo_file){
  
  Kw_ids <- dir(Kw_file_dir) %>% stringr::str_remove('_kw.csv')
    
  
  meteo_ids <- readRDS(meteo_file) %>% mutate(is_file = file.exists(file.path('../lake-temperature-model-prep/7_drivers_munge/out', basename(meteo_fl)))) %>% 
    filter(is_file) %>% pull(site_id)
  HA_values <- readRDS(H_A_file) %>% sapply(FUN = function(x)length(x$H) > 2)
  HA_ids <- readRDS(H_A_file)[HA_values] %>% names
  nhd_ids <- readRDS(xwalk_file) %>% filter(site_id %in% Kw_ids) %>% 
    filter(site_id %in% HA_ids) %>% 
    filter(site_id %in% meteo_ids) %>% 
    pull(site_id) %>% unique()
  
  return(nhd_ids)
    
}

prep_model <- function(glm_target, site_id, nml_list, nml_obj, kw_file){
  
  dir.create(dirname(glm_target))
  nml_obj <- split_list(nml_list, site_id) %>% 
    build_nml()
  
  write_nml(nml_obj, file = glm_target)
  dir.create(file.path(dirname(glm_target), 'drivers'))

  meteo_fl <- glmtools::get_nml_value(nml_obj, 'meteo_fl') %>% basename() %>% 
    file.path('../lake-temperature-model-prep/7_drivers_munge/out', .)
  
  readr::read_csv(meteo_fl, col_types = 'Dddddddd') %>% mutate(Rain = Rain *0) %>% readr::write_csv(file.path(dirname(glm_target), nml_obj$meteorology$meteo_fl))
  
  read_csv(kw_file, col_types = 'Dd') %>% write_csv(file.path(dirname(glm_target), nml_obj$light$Kw_file))
  
}

build_nml <- function(nml_args){
  nml_args$sim_name <- nml_args$site_id
  nml_args$site_id <- NULL
  
  nml_args$sw_factor = 1
  nml_args$start = '1980-04-01'
  nml_args$stop = '2018-12-31'
  nml_args$dt=3600
  
  nml_args$bsn_vals = length(nml_args$H)
  nml_args$the_depths = c(0, floor(nml_args$lake_depth*100)/100)
  nml_args$timefmt = 2
  nml <- read_nml(GLM3r::nml_template_path()) %>% 
    set_nml(arg_list = nml_args)
  nml$debugging <- list(disable_evap = TRUE)
  nml$sediment <- NULL
  # $sed_heat_Ksoil = 0.01
  # nml$sediment$n_zones = 3
  # nml$sediment$sed_temp_mean = c(7, 7.5, 8)
  nml$light$Kw_file <- 'Kw_file.csv'
  return(nml)
}

split_list <- function(list = nml_list, name){
  return(list[[name]])
}


run_glm3 <- function(ind_path, nml_filename, nc_filename){
  # nml <- read_nml(nml_file)
  # nml$light$Kw_file <- NULL
  # write_nml(nml, nml_file)
  GLM3r::run_glm(sim_folder = dirname(nml_filename), verbose = FALSE)
  nml <- read_nml(nml_filename)
  sim_output <- paste0(dirname(nml_filename), '/', get_nml_value(nml, 'out_dir'), '/', get_nml_value(nml, 'out_fn'), '.nc')
  file.copy(from = sim_output, to = nc_filename, overwrite = TRUE)
  unlink(sim_output)
  yaml::write_yaml(file = ind_path, x = as.list(tools::md5sum(nc_filename)))
}
  
