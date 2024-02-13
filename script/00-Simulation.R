#### Preamble ####
# Purpose: Simulation for data
# Author: Yingzhi Zhang
# Date: 10 Feburary 2024
# Contact: yingzhi.zhang@mail.utoronto.ca
# License: MIT


#### Workspace setup ####
library(tidyverse)


#### Simulate data ####
set.seed(90)

category_name <- c("Acts of Violence", "Threat", "Other", "Reported Exposure", "Fear", "Anger")

simulated_data <-
  tibble(
    "year" = rep(1985:2023, each = 6),
    "category" = rep(category_name, times = 39),
    "number" = round (runif(n = 39*6, min = 0, max = 1000), 0)
  )

head(simulated_data)





