###########################################################
# Filename: calibration_compile.R
# Author: Nicholas von Turkovich
# Date: 11/23/2021
# Note(s): 
###########################################################

# environment setup -------------------------------------------------------
rm(list = ls())

library(tidyverse)
library(here)

setwd(paste0(here(), "/Output/Logs/Calibration_Deterministic/"))
print(getwd())

# Calibration type
args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0){
  calib_type = args[2]
} else {
  calib_type = "calibration"
}

algs = c("NLOPT_GN_DIRECT", "NLOPT_GN_DIRECT_L_RAND", "NLOPT_G_MLSL_LDS", "NLOPT_GN_CRS2_LM", "NLOPT_GN_ESCH")

# Compile for particular algorithm ----------------------------------------------------------------

for (i in 1:length(algs)){

  files = list.files()
  filename = paste0(calib_type, "\\d+", "_", algs[i])
  files = files[regexpr(paste0(filename, ".csv"), files) == 1]
  files = sort(files)

  if (length(files) == 0){
    next
  }
  calibration = files %>%
    map(.f = read.csv)

  for (j in 1:length(files)){

    if (substr(files[j],nchar(calib_type)+2,nchar(calib_type)+2) == "_"){
      calibration[[j]] = calibration[[j]] %>%
        mutate(calibration = substr(files[j],nchar(calib_type)+1,nchar(calib_type)+1))
    } else {
      calibration[[j]] = calibration[[j]] %>%
        mutate(calibration = substr(files[j],nchar(calib_type)+1,nchar(calib_type)+2))
    }

  }

  if (calib_type != "calibration_entry") {
    calibration = calibration %>%
      map(.f = function(X){
        X = X %>% mutate(tfp_growth_target = X[3,"tfp_growth"],
                         rd_gdp_target = X[3,"rd_gdp"],
                         mean_markup_target = X[3,"mean_markup"],
                         median_markup_target = X[3,"median_markup"],
                         pct75_markup_target = X[3,"pct75_markup"],
                         pct90_markup_target = X[3,"pct90_markup"],
                         pct95_markup_target = X[3,"pct95_markup"],
                         pct99_markup_target = X[3,"pct99_markup"],
                         std_markup_target = X[3,"std_markup"],
                         FHK_WITHIN_adj_sh_5y_target = X[3,"FHK_WITHIN_adj_sh_5y"],
                         mean_innov_target = X[3,"mean_innov"],
                         pct90_innov_target = X[3,"pct90_innov"],
                         rd_leader_share_target = X[3,"rd_leader_share"],
                         mu_end_target = X[3,"mu_end"],
                         B_lower = X[1,"B"],
                         eta_lower = X[1,"eta"],
                         lambda_lower = X[1,"lambda"],
                         phi_lower = X[1,"phi"],
                         B_upper = X[2,"B"],
                         eta_upper = X[2,"eta"],
                         lambda_upper = X[2,"lambda"],
                         phi_upper = X[2,"phi"])})
  } else {
    calibration = calibration %>%
      map(.f = function(X){
        X = X %>% mutate(tfp_growth_target = X[3,"tfp_growth"],
                         rd_gdp_target = X[3,"rd_gdp"],
                         mean_markup_target = X[3,"mean_markup"],
                         median_markup_target = X[3,"median_markup"],
                         pct75_markup_target = X[3,"pct75_markup"],
                         pct90_markup_target = X[3,"pct90_markup"],
                         pct95_markup_target = X[3,"pct95_markup"],
                         pct99_markup_target = X[3,"pct99_markup"],
                         std_markup_target = X[3,"std_markup"],
                         FHK_ENTRY_sh_5y_target = X[3,"FHK_ENTRY_sh_5y"],
                         mean_innov_target = X[3,"mean_innov"],
                         pct90_innov_target = X[3,"pct90_innov"],
                         rd_leader_share_target = X[3,"rd_leader_share"],
                         mu_end_target = X[3,"mu_end"],
                         emp5_target = X[3,"emp5"],
                         emp10_target = X[3,"emp10"],
                         B_lower = X[1,"B"],
                         eta_lower = X[1,"eta"],
                         lambda_lower = X[1,"lambda"],
                         phi_lower = X[1,"phi"],
                         B_e_lower = X[1,"B_e"],
                         phi_e_lower = X[1,"phi_e"],
                         B_upper = X[2,"B"],
                         eta_upper = X[2,"eta"],
                         lambda_upper = X[2,"lambda"],
                         phi_upper = X[2,"phi"],
                         B_e_upper = X[2,"B_e"],
                         phi_e_upper = X[2,"phi_e"])})

  }

  calibration = calibration %>%
    map(.f = function(X){
      X = tail(X, -4)
      X = X %>%
        mutate(run = row_number(),
               iter = dim(X)[1])
      return(X)
    })

  calibration = calibration %>%
    map(.f = function(X){
      X = X %>%
        arrange(loss) %>%
        slice(1)
      return(X)
    }) %>% reduce(.f = bind_rows); calibration


  write.csv(calibration, paste0("./Compiled/compiled_", calib_type, "_", algs[i], ".csv"), row.names = F)

}

# Compile broadly ---------------------------------

files = list.files("./Compiled/")
files = files[startsWith(files, paste0("compiled_", calib_type, "_NLOPT"))]
files = paste0("./Compiled/", files)
files = files[grepl(".csv", files)]
algorithms = gsub(paste0("./Compiled/compiled_", calib_type, "_"), "", files)
algorithms = gsub(".csv", "", algorithms)

# Create master csv
compiled = files %>%
  map(.f = read.csv)

for (i in 1:length(algorithms)){
  compiled[[i]] = compiled[[i]] %>% mutate(algorithm = algorithms[i])
}

compiled = compiled %>%
  reduce(.f = bind_rows); head(compiled %>% arrange(loss))

write.csv(compiled, paste0("./Compiled/master_compiled_", calib_type, ".csv"), row.names = FALSE)
  
