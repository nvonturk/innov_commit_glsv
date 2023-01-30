###########################################################
# Filename: tripcons_compile.R
# Author: Nicholas von Turkovich
# Date: 4/18/2021
# Note(s): 
###########################################################

# environment setup -------------------------------------------------------
rm(list = ls())

library(tidyverse)
library(here)

setwd(paste0(here(), "/Output/Logs/Main/"))
print(getwd())

args = commandArgs(trailingOnly = TRUE)
if (length(args) == 0){
  stop("ERROR: No filename supplied")
} else if (length(args) > 1) {
  stop("ERROR: Too many arguments supplied")
} else {
  filename = args[1]
}

# Pull in data from different eta1s ----------------------------------------------------------------

# Keep only the log for the minimum loss search
files = list.files()
files = files[grepl(filename, files)]
files = files[grepl(".csv", files)]
files = files[!grepl("curve", files) & !grepl("compiled", files)]

print(files)

if (length(files) > 0){
  
  # Read all files and take the entry with minimum loss
  merged = files %>%
    map(.f = function(X){
      X = read.csv(X) %>%
        arrange(loss) %>%
        slice(1)
    }) %>% reduce(.f = bind_rows) %>%
    filter(loss < 1e-3) %>%
    arrange(heta1)
  
  # Keep all the guesses in one long format file
  merged_long = files %>%
    map(.f = function(X){
      X = read.csv(X)
    }) %>% reduce(.f = bind_rows)
  
  write.csv(merged, paste0(filename, "_compiled.csv"), row.names = F)
  write.csv(merged_long, paste0(filename, "_compiled_long.csv"), row.names = F)
  
} else {
  write.csv(data.frame(matrix(ncol = 1, nrow = 0)), paste0(filename, "_compiled.csv"), row.names = F)
}

