#### Preamble ####
# Purpose: Data cleaning and processing
# Author: Yingzhi Zhang
# Date: 10 Feburary 2024
# Contact: yingzhi.zhang@mail.utoronto.ca
# License: MIT


#### Workspace setup ####
library(tidyverse)
library(janitor)
library(readxl)
library(haven)

###### Import and split data ###### 
raw_meta_data <- read_excel("input/replication_package/data/02-metaanalysis-data.xlsx", sheet = "Data", 
                  col_names = TRUE, skip = 3)

writexl::write_xlsx(raw_meta_data, "output/data/raw_meta_data.xlsx")

# dataframe of studies
dat_studies <- raw_meta_data[c(1:87)] #studies.
dat_studies <- group_by(dat_studies, ID_RS) %>% slice(1)
sum(raw_meta_data$ESS, na.rm = TRUE) #number of unique respondents

# dataframe of reports
dat_reports <- raw_meta_data[c(1:87)] #reports/publications.
dat_reports <- group_by(dat_reports, ID_R) %>% slice(1)

# set wished layout settings
mytheme <- theme(plot.title = element_text(face = "bold", size = (22), colour = "black"), 
                 axis.text = element_text(size = (18), colour = "black"),
                 axis.ticks.length = unit(0.5, "cm"),
                 axis.title.y = element_text(size = (22), colour = "black"),
                 axis.title.x = element_text(size = (22), colour = "white"))

# Save the data for US
US_meta_data <-
  raw_meta_data %>% 
  filter(Country == "US")
writexl::write_xlsx(US_meta_data, "output/data/US_meta_data.xlsx")

# Save the data for Israel
Israel_meta_data <-
  raw_meta_data %>% 
  filter(Country == "Israel"| Country == "Israel + NI")
writexl::write_xlsx(Israel_meta_data, "output/data/Israel_meta_data.xlsx")
