#### Preamble ####
# Purpose: Tests
# Author: Yingzhi Zhang
# Date: 10 Feburary 2024
# Contact: yingzhi.zhang@mail.utoronto.ca
# License: MIT

library(tidyverse)
library(janitor)
library(readxl)
library(haven)
library(broom)
library(dotwhisker)

raw_meta_data <- read_xlsx(paste0(gsub("/script", "", getwd()), "/output/data/raw_meta_data.xlsx"))

sum(raw_meta_data$ESS, na.rm = TRUE) ==1987495

#outgroup hostility data
dat_outgroup <- read_xlsx(paste0(gsub("/script", "", getwd()), "/output/data/dat_outgroup.xlsx"))
length(unique(dat_outgroup$ID_R)) == 126

#conservative shift data
dat_conservatism <- read_xlsx(paste0(gsub("/script", "", getwd()), "/output/data/dat_conservatism.xlsx"))
length(unique(dat_conservatism$ID_R)) == 144

#rally data
dat_rally <- read_xlsx(paste0(gsub("/script", "", getwd()), "/output/data/dat_rally.xlsx"))
length(unique(dat_rally$ID_R)) == 72

#check total number
count(dat_outgroup) + count(dat_conservatism) + count(dat_rally) == 1733
