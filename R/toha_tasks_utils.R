create_toha_model_plan <- function(ids, nml_list){
  
  

  write_nml_step <- create_task_step(  
    step_name = 'write_nml',
    target_name = function(task_name, step_name, ...) {
      sprintf("toha_models/%s/glm3.nml", task_name)
    },
    command = function(task_name, step_name, ...) {
      sprintf("prep_model(target_name, 
      site_id = I('%s'), nml_list = nml_list,
      kw_file = '~/Downloads/USGS_Kd_Predictions_byNHDHR/%s_kw.csv')", task_name, task_name)
    }
  )
  
  
  run_model_step <- create_task_step(
    step_name = 'run_model',
    target_name = function(task_name, step_name, ...) {
      sprintf("toha_models/%s.ind", task_name)
    },
    command = function(task_name, step_name, ...) {
      sprintf("run_glm3(target_name, 
      nml_filename = 'toha_models/%s/glm3.nml', 
      nc_filename = I('/Volumes/ThunderBlade/TOHA_GLM3/%s/output.nc'))",task_name, task_name)
    }
  )
  
  export_temp_step <- create_task_step(
    step_name = 'export_temp',
    target_name = function(task_name, step_name, ...) {
      sprintf("toha_models/%s/%s_temperature.feather", task_name, task_name)
    },
    command = function(task_name, step_name, ...) {
      sprintf("export_temp_glm3(target_name, 'toha_models/%s.ind', 'toha_models/%s/glm3.nml')", task_name, task_name)
    }
  )
  
  calc_toha_step <- create_task_step(
    step_name = 'calc_toha',
    target_name = function(task_name, step_name, ...) {
      sprintf("'toha_models/shared/%s_toha.csv'", task_name)
    },
    command = function(task_name, step_name, ...) {
      sprintf("estimate_toha(target_name, 'toha_models/%s.ind')", task_name)
    }
  )

  create_task_plan(names(nml_list)[names(nml_list) %in% ids], list(write_nml_step, run_model_step, export_temp_step, calc_toha_step), final_steps='run_model', add_complete = FALSE)# 'calc_toha', add_complete = FALSE)
  
}


create_toha_model_makefile <- function(makefile, task_plan, final_targets){
  include <- "remake.yml"
  packages <- c('GLM3r', 'readr', 'glmtools', 'mda.lakes')
  sources <- c('R/glm_setup_utils.R', 'R/toha_calc_utils.R')
  
  create_task_makefile(task_plan, makefile, include = include, packages = packages, sources = sources, final_targets = final_targets)
}





