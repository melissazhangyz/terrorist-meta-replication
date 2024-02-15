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
library(broom)
library(dotwhisker)


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

### Rregression result
# Data preparation: Split dataset in 3 based on outcome type
# 1) outgroup dataset
dat_outgroup <- raw_meta_data[raw_meta_data$PA_Category == 10 |
                      raw_meta_data$PA_Category == 99, ]
length(unique(dat_outgroup$ID_R))#number of manuscripts (N_j)
writexl::write_xlsx(dat_outgroup, "output/data/dat_outgroup.xlsx")

# 2) conservative shift dataset
dat_conservatism <- raw_meta_data[raw_meta_data$PA_Category == 5 |
                          raw_meta_data$PA_Category == 6 |
                          raw_meta_data$PA_Category == 7 |
                          raw_meta_data$PA_Category == 8 |
                          raw_meta_data$PA_Category == 9 |
                          raw_meta_data$PA_Category == 11, ]
length(unique(dat_conservatism$ID_R))#number of manuscripts (N_j)
writexl::write_xlsx(dat_conservatism, "output/data/dat_conservatism.xlsx")

# 3) rally-around-the-flag dataset
dat_rally <- raw_meta_data[raw_meta_data$PA_Category == 1 |
                   raw_meta_data$PA_Category == 2 |
                   raw_meta_data$PA_Category == 3 |
                   raw_meta_data$PA_Category == 4, ]
length(unique(dat_rally$ID_R))#number of manuscripts (N_j)
writexl::write_xlsx(dat_rally, "output/data/dat_rally.xlsx")

# double check (total number of effect sizes; N_k)
count(dat_outgroup) + count(dat_conservatism) + count(dat_rally)

#Outgroup Location
#check data
summary(as.factor(dat_outgroup$Country))
#recode data
dat_outgroup$Country2[dat_outgroup$Country=="US"] <- "US"
dat_outgroup$Country2[dat_outgroup$Country=="Israel"] <- "Israel"
dat_outgroup$Country2[dat_outgroup$Country=="Spain"] <- "Spain"
dat_outgroup$Country2[dat_outgroup$Country=="France"] <- "France"
dat_outgroup$Country2[dat_outgroup$Country=="UK"] <- "UK"
dat_outgroup$Country2[dat_outgroup$Country=="Canada"] <- "Canada"
dat_outgroup$Country2[is.na(dat_outgroup$Country2)] <- "Other"
dat_outgroup$Country2 <- as.factor(dat_outgroup$Country2)
table(dat_outgroup$Country2)
#convert characters into dummy variables
dat_outgroup$US <- ifelse(dat_outgroup$Country=="US", yes=1, no=0)
dat_outgroup$Israel <- ifelse(dat_outgroup$Country=="Israel", yes=1, no=0)
dat_outgroup$Spain <- ifelse(dat_outgroup$Country=="Spain", yes=1, no=0)
dat_outgroup$France <- ifelse(dat_outgroup$Country=="France", yes=1, no=0)
dat_outgroup$UK <- ifelse(dat_outgroup$Country=="UK", yes=1, no=0)
dat_outgroup$Canada <- ifelse(dat_outgroup$Country=="Canada", yes=1, no=0)
dat_outgroup$Other_c <- ifelse(dat_outgroup$US==0 
                               & dat_outgroup$Israel==0, yes=1, no=0)

#fit model and plot results
set.seed(1234)
og_location <- rma.mv(y=Fisher, V=Variance_F, mods = cbind(US, Israel, Spain, France, UK, Canada),
                      intercept=F, 
                      random = ~ 1 | ID_R/ID_ES_Unique, data=dat_outgroup,
                      method = "ML")
summary(og_location)
og_location_df <- tidy(og_location) # create data.frame of regression results
og_location_df # a tidy data.frame available for dwplot
og_location_df <- og_location_df %>%  #re-label variable
  relabel_predictors(c(US = "United States",
                       Israel  = "Israel",
                       Spain = "Spain",
                       France = "France",
                       UK = "UK",
                       Canada = "Canada"))
write_csv(og_location_df, "output/data/og_location_df.csv")

