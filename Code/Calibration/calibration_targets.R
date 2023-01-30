###########################################################
# Filename: calibration_targets.R
# Author: Nicholas von Turkovich
# Date: 11/22/2021
# Note(s): Code compiles moment targets for calibration_deterministic.m 
###########################################################

# Environment setup -------------------------------------------------------
rm(list = ls())

library(tidyverse)
library(here)

dir = dirname(sys.frame(1)$ofile)
dir = gsub("Code/Calibration", "", dir)
setwd(dir)

# Settings ----------------------------------------------------------------

set.seed(2019)

sample_period_long <- seq(1960, 2019, by = 1)
sample_period_medium <- seq(1970, 1999, by = 1)
sample_period_short <- seq(2000, 2019, by = 1)
sample_period_vshort <- seq(2004,2019,by = 1)

sample_periods <- list(sample_period_long,
                       sample_period_medium,
                       sample_period_short,
                       sample_period_vshort)

sample_periods_names <- list("1960 - 2019", "1970 - 1999", "2000 - 2019", "2004 - 2019")

# TFP data --------------------------------------------------------

# Read in data from FRBSF: https://www.frbsf.org/economic-research/indicators-data/total-factor-productivity-tfp/
tfp <- readxl::read_xlsx(path = "Output/Tracked/Data/quarterly_tfp.xlsx", 
                         sheet = "annual")

# Compute the average growth rate for each sample period
tfp_stats <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(tfp_stats) = c("n", "dtfp_mean", "dtfp_util_mean")

for (i in 1:length(sample_periods)){
  
  tfp_temp <- tfp %>%
    filter(date %in% sample_periods[[i]]) %>%
    summarise(n = n(),
              dtfp_mean = mean(dtfp),
              dtfp_util_mean = mean(dtfp_util))
  
  tfp_stats <- bind_rows(tfp_stats, tfp_temp)
  
}

row.names(tfp_stats) <- sample_periods_names


# Hall's Lerner data ------------------------------------------------

# Read in data from https://web.stanford.edu/~rehall/Recent_Unpublished_Papers.html
lerner <- readxl::read_xlsx(path = "Output/Tracked/Data/MPresults.xlsx", 
                         sheet = "Trend by ind",
                         range = "U5:X33")
lerner <- cbind(lerner[,1], lerner[,4]) %>% data.frame
colnames(lerner) <- c("date", "wgt_avg_lerner")

lerner_means <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(lerner_means) <- c("n", "mean_lerner")

# Take the average Lerner value for each different time span
for (i in 1:length(sample_periods)){
  
  lerner_temp <- lerner %>%
    filter(date %in% sample_periods[[i]]) %>%
    summarise(n = n(),
              mean_lerner = mean(wgt_avg_lerner))
  
  lerner_means <- bind_rows(lerner_means, lerner_temp)
  
}

# Take series of draws from Lerner distribution in paper for comparison
lerner_draws <- rbeta(1e6, 1.36, 8)
lerner_means <- rbind(c(NA, mean(lerner_draws)), lerner_means)
row.names(lerner_means) <- c("Paper", sample_periods_names)

# Given the two means, compute the implied alpha parameter for the beta distribution (note for the Paper row this should return ~1.36)
lerner_means$beta = 8
lerner_means$alpha = (lerner_means$beta*lerner_means$mean_lerner)/(1 - lerner_means$mean_lerner)

lerner_stats = as.data.frame(matrix(NA, length(sample_periods)+1, 3))
colnames(lerner_stats) = c("mean_markup", "median_markup", "pct90_markup")

# Generate draws from beta distribution and compute stats on markups
for (i in 1:(length(sample_periods)+1)){
  
  lerner_draws <-  rbeta(1e6, lerner_means$alpha[i], lerner_means$beta[i])
  lerner_markups <- (1/(1-lerner_draws) - 1)*100
  
  markup_dist = c(mean(lerner_markups), median(lerner_markups), quantile(lerner_markups, probs = 0.9))

  lerner_stats[i,] = markup_dist
  
}

hall_stats <- round(cbind(lerner_means, lerner_stats), 2); hall_stats

# R&D to GDP data from FRED -----------------------------------------------

clean_fred <- function(X){
  quantmod::getSymbols(X, src = "FRED")
  temp <- data.frame(date = zoo::index(eval(sym(X))),
                     zoo::coredata(eval(sym(X))))
  colnames(temp) = tolower(colnames(temp))
  
  return(temp)
}

fred_symbols <- c("B985RX1Q020SBEA", "Y694RX1Q020SBEA", "GDPC1")

fred_series <- purrr::map(fred_symbols, clean_fred) %>%
  purrr::reduce(full_join, "date") %>%
  arrange(date)

fred_series <- fred_series %>%
  mutate(rd_gdp = y694rx1q020sbea/gdpc1*100,
         rd_gdp_wsoftware = (y694rx1q020sbea + b985rx1q020sbea) / gdpc1 * 100) %>%
  mutate(year = lubridate::year(date))

# Compute the average ratio for each sample period

rdgdp_stats <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(rdgdp_stats) = c("n", "rd_gdp_mean", "rd_gdp_wsoftware_mean")

for (i in 1:length(sample_periods)){
  
  rdgdp_temp <- fred_series %>%
    filter(year %in% sample_periods[[i]]) %>%
    summarise(n = n(),
              rd_gdp_mean = mean(rd_gdp, na.rm = T),
              rd_gdp_wsoftware_mean = mean(rd_gdp_wsoftware, na.rm = T))
  
  rdgdp_stats <- bind_rows(rdgdp_stats, rdgdp_temp)
  
}

row.names(rdgdp_stats) <- sample_periods_names

print("---------------------------------------------------")
print(tfp_stats)

print("---------------------------------------------------")
print(hall_stats)

print("---------------------------------------------------")
print(rdgdp_stats)

