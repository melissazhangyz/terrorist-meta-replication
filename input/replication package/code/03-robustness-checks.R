# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
# PROJECT:    How Terrorism Does (and Does Not) Affect Citizens' 
#             Political Attitudes: A Meta-Analysis. 
# AUTHOR:     Amelie Godefroidt
# CONTACT:    amelie.godefroidt@ntnu.no
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
# This R file contains the code necessary to replicate the **robustness checks**,
# including the sensitivity analyses and publication bias analyses (Supp. Info.).
# 
# Last successful replication: 01.11.2021
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 

# first things first: set your settings.

rm(list=ls())
setwd(".../Godefroidt_2021_ReplicationFiles/") #set your directory

### layout settings
mytheme <- theme(plot.title = element_text(face = "bold", size = (20), colour = "black"), 
                 axis.text = element_text(size = (17), colour = "black"),
                 axis.ticks.length = unit(0.4, "cm"),
                 axis.title = element_text(size = (22), colour = "black"))


###### Install and load packages ###### 

# ipak function: install and load multiple R packages.
# check to see if packages are installed. Install them if they are not, then load them into the R session.

ipak <- function(pkg){  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
if(length(new.pkg)) install.packages(new.pkg, dependencies = TRUE)
sapply(pkg, require, character.only = TRUE)
}

packages <- c("readxl", "psych", "dplyr",
              "metaSEM", "metafor", "metaviz", "PerformanceAnalytics", "weightr", 
              "ggplot2", "ggpubr", "wesanderson", "stringr", "gridExtra", "grid") 
ipak(packages)


###### Import and split data ###### 

### import data file
dat <- read_excel("data/02-metaanalysis-data.xlsx", sheet = "Data", 
                  col_names = TRUE, skip = 3)

### dataframe of studies and reports
dat_studies <- dat[c(1:87)] #studies.
dat_studies <- group_by(dat_studies, ID_RS) %>% slice(1)
sum(dat_studies$ESS, na.rm = TRUE) #number of unique respondents
dat_reports <- dat[c(1:87)] #reports/publications.
dat_reports <- group_by(dat_reports, ID_R) %>% slice(1)

### split dataset in 3
# 1) outgroup dataset
dat_outgroup <- dat[dat$PA_Category == 10 |
                      dat$PA_Category == 99, ]
# 2) conservative shift dataset
dat_conservatism <- dat[dat$PA_Category == 5 |
                          dat$PA_Category == 6 |
                          dat$PA_Category == 7 |
                          dat$PA_Category == 8 |
                          dat$PA_Category == 9 |
                          dat$PA_Category == 11, ]
# 3) rally-around-the-flag dataset
dat_rally <- dat[dat$PA_Category == 1 |
                   dat$PA_Category == 2 |
                   dat$PA_Category == 3 |
                   dat$PA_Category == 4, ]
# double check (n)
count(dat_outgroup) + count(dat_conservatism) + count(dat_rally)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
##### Main 3-level models reported in the manuscript (for comparison)
##### Table C.4. in the Supplementary Information
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
#outgroup hostility
model_og_LBCI <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                       data = dat_outgroup, 
                       I2=c("I2q", "ICC"),
                       intervals.type = "LB",
                       model.name = "3 level outgroup LB")
summary(model_og_LBCI)

#conservative shift
model_cs_LBCI <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                       data = dat_conservatism, 
                       I2=c("I2q", "ICC"),
                       intervals.type = "LB",
                       model.name = "3 level conservatism LB")
summary(model_cs_LBCI)

#rally-around-the-flag
model_rf_LBCI <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                       data = dat_rally, 
                       I2=c("I2q", "ICC"),
                       intervals.type = "LB",
                       model.name = "3 level rally LB")
summary(model_rf_LBCI)


##### Same model but fitted with rma.mv instead of meta3 #####

#outgroup hostility
rma_mv_og <- rma.mv(y=Fisher, V=Variance_F, 
                          random = ~ 1 | ID_R/ID_ES_Unique, data=dat_outgroup,
                          method = "ML") #fit 3-Level with RMA to check consistency
summary(rma_mv_og)

#conservative shift
rma_mv_cs <- rma.mv(y=Fisher, V=Variance_F, 
                          random = ~ 1 | ID_R/ID_ES_Unique, data=dat_conservatism,
                          method = "ML") #fit 3-Level with RMA to check consistency
summary(rma_mv_cs)

#rally effects
rma_mv_rf <- rma.mv(y=Fisher, V=Variance_F, 
                          random = ~ 1 | ID_R/ID_ES_Unique, data=dat_rally,
                          method = "ML") #fit 3-Level with RMA to check consistency
summary(rma_mv_rf)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
###### A. SENSITIVITY ANALYSES ######
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
##### 1. Sensitivity to additional clustering (Table C.5) ##### 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 

##### 1.1. Outgroup Hostility ##### 
#sensitivity analysis 1: Changing the Level-3
model_og_LBCI_RS <- meta3(y=Fisher, v=Variance_F, cluster=ID_RS,
                          data = dat_outgroup, 
                          I2=c("I2q", "ICC"),
                          intervals.type = "LB",
                          model.name = "3 level outgroup LB_RS")
summary(model_og_LBCI_RS)
#sensitivity analysis 2: Adding a Level-4
four_levels_og <- rma.mv(y=Fisher, V=Variance_F, 
                         random = ~ 1 | ID_R/ID_RS/ID_ES_Unique, 
                         data=dat_outgroup,
                         method = "ML")
summary(four_levels_og)

##### 1.2. Conservatism ##### 
#sensitivity analysis 1: Changing the Level-3
model_cs_LBCI_RS <- meta3(y=Fisher, v=Variance_F, cluster=ID_RS,
                          data = dat_conservatism, 
                          I2=c("I2q", "ICC"),
                          intervals.type = "LB",
                          model.name = "3 level conservatism LB_RS")
summary(model_cs_LBCI_RS)
#sensitivity analysis 2: Adding a Level-4
four_levels_cs <- rma.mv(y=Fisher, V=Variance_F, 
                         random = ~ 1 | ID_R/ID_RS/ID_ES_Unique, 
                         data=dat_conservatism)
summary(four_levels_cs)

##### 1.3. Rally Effects ##### 
#sensitivity analysis 1: Changing the Level-3
model_rf_LBCI_RS <- meta3(y=Fisher, v=Variance_F, cluster=ID_RS,
                          data = dat_rally, 
                          I2=c("I2q", "ICC"),
                          intervals.type = "LB",
                          model.name = "3 level rally LB_RS")
summary(model_rf_LBCI_RS)
#sensitivity analysis 2: Adding a Level-4
four_levels_rf <- rma.mv(y=Fisher, V=Variance_F, 
                         random = ~ 1 | ID_R/ID_RS/ID_ES_Unique, 
                         data=dat_rally)
summary(four_levels_rf)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
##### 2. Sensitivity to outliers (Table C.6) #####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 

##### 2.1. Outgroup Hostility ##### 

##### 2.1.1. Extreme Effect Sizes
confint(model_og_LBCI)[1] #lowerbound
confint(model_og_LBCI)[4] #upperbound
# Calculate CI for all observed effect sizes
dat_outgroup$upperci <- dat_outgroup$Fisher + 1.96 * sqrt(dat_outgroup$Variance_F)
dat_outgroup$lowerci <- dat_outgroup$Fisher - 1.96 * sqrt(dat_outgroup$Variance_F)
# Create filter variable
dat_outgroup$outlier <- dat_outgroup$upperci < confint(model_og_LBCI)[1] | dat_outgroup$lowerci > confint(model_og_LBCI)[4]
# Count number of outliers:
sum(dat_outgroup$outlier)
# Make a basic plot: specify that the x-variable is the effect size, and 
# the color and fill of the histogram bars are based on the value of 'outlier'
ggplot(data = dat_outgroup, aes(x = Fisher, colour = outlier, fill = outlier)) +
  # Add a histogram with transparent bars (alpha = .2)
  geom_histogram(alpha = .2, bins = 50) +
  # Add a vertical line at the pooled effect value (m_re$b[1])
  geom_vline(xintercept = coef(model_og_LBCI)[1], colour = "grey10", linetype = 2, size = 1) +
  # Apply a black and white theme
  theme_classic2() +
  scale_fill_grey(start = 0.2, end = 0.6) +
  scale_color_grey(start = 0.2, end = 0.6) +
  labs(x = "Fisher's Z Correlation Coefficient", 
       y = "", 
       colour="Outlier", fill="Outlier") 
#sensitivity analysis 1: outliers excluded 
dat_outgroup_no <- dat_outgroup[ which(dat_outgroup$outlier==FALSE), ]
model_og_no1 <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                     data = dat_outgroup_no, 
                     I2=c("I2q", "ICC"),
                     intervals.type = "LB",
                     model.name = "3 level og no1")
summary(model_og_no1)
model_og_no1 <- rerun(model_og_no1)
summary(model_og_no1)
#sensitivity analysis 2: outliers winsorised 
#set extremely low ES's to lower bound
#set extremely high ES's to higher bound
dat_outgroup <- dat_outgroup %>% 
  mutate(
    F_1 = if_else(dat_outgroup$upperci < confint(model_og_LBCI)[1], 
                  confint(model_og_LBCI)[1],
                  Fisher),
    F_2 = if_else(dat_outgroup$lowerci > confint(model_og_LBCI)[4], 
                  confint(model_og_LBCI)[4],
                  F_1))
model_og_no1w <- meta3(y=F_2, v=Variance_F, cluster=ID_R,
                      data = dat_outgroup, 
                      I2=c("I2q", "ICC"),
                      intervals.type = "LB",
                      model.name = "3 level og no1w")
summary(model_og_no1w)
model_og_no1w <- rerun(model_og_no1w)
summary(model_og_no1w)

##### 2.1.2. Extreme Variances
dat_outgroup <- dat_outgroup %>% 
  mutate(
    small_var = Variance_F <= quantile(Variance_F, .1, na.rm = TRUE),
    Variance_F_2 = if_else(Variance_F <= quantile(Variance_F, .1, na.rm = TRUE), 
                    quantile(Variance_F, .1, na.rm = TRUE),
                    Variance_F))
#sensitivity analysis 1: outliers excluded 
dat_outgroup_no2 <- dat_outgroup[ which(dat_outgroup$small_var==FALSE), ]
model_og_no2 <- meta3(y=Fisher, v=Variance_F_2, cluster=ID_R,
                      data = dat_outgroup_no2, 
                      I2=c("I2q", "ICC"),
                      intervals.type = "LB",
                      model.name = "3 level og no2")
summary(model_og_no2)
#sensitivity analysis 2: outliers winsorised 
model_og_no2w <- meta3(y=Fisher, v=Variance_F_2, cluster=ID_R,
                     data = dat_outgroup, 
                     I2=c("I2q", "ICC"),
                     intervals.type = "LB",
                     model.name = "3 level og no2w")
summary(model_og_no2w)


##### 2.2. Conservatism #####  

##### 2.2.1. Extreme Effect Sizes
confint(model_cs_LBCI)[1] #lowerbound
confint(model_cs_LBCI)[4] #upperbound
# Calculate CI for all observed effect sizes
dat_conservatism$upperci <- dat_conservatism$Fisher + 1.96 * sqrt(dat_conservatism$Variance_F)
dat_conservatism$lowerci <- dat_conservatism$Fisher - 1.96 * sqrt(dat_conservatism$Variance_F)
# Create filter variable
dat_conservatism$outlier <- dat_conservatism$upperci < confint(model_cs_LBCI)[1] | dat_conservatism$lowerci > confint(model_cs_LBCI)[4]
# Count number of outliers:
sum(dat_conservatism$outlier)
# Make a basic plot: specify that the x-variable is the effect size, and 
# the color and fill of the histogram bars are based on the value of 'outlier'
ggplot(data = dat_conservatism, aes(x = Fisher, colour = outlier, fill = outlier)) +
  # Add a histogram with transparent bars (alpha = .2)
  geom_histogram(alpha = .2, bins = 50) +
  # Add a vertical line at the pooled effect value (m_re$b[1])
  geom_vline(xintercept = coef(model_cs_LBCI)[1], colour = "grey10", linetype = 2, size = 1) +
  # Apply a black and white theme
  theme_classic2() +
  scale_fill_grey(start = 0.2, end = 0.6) +
  scale_color_grey(start = 0.2, end = 0.6) +
  labs(x = "Fisher's Z Correlation Coefficient", 
       y = "", 
       colour="Outlier", fill="Outlier") 
#sensitivity analysis 1: outliers excluded 
dat_conservatism_no <- dat_conservatism[ which(dat_conservatism$outlier==FALSE), ]
model_cs_no1 <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                     data = dat_conservatism_no, 
                     I2=c("I2q", "ICC"),
                     intervals.type = "LB",
                     model.name = "3 level cs no1")
summary(model_cs_no1)
#sensitivity analysis 2: outliers winsorised 
#set extremely low ES's to lower bound
#set extremely high ES's to higher bound
dat_conservatism <- dat_conservatism %>% 
  mutate(
    F_1 = if_else(dat_conservatism$upperci < confint(model_cs_LBCI)[1], 
                  confint(model_cs_LBCI)[1],
                  Fisher),
    F_2 = if_else(dat_conservatism$lowerci > confint(model_cs_LBCI)[4], 
                  confint(model_cs_LBCI)[4],
                  F_1))
model_cs_no1w <- meta3(y=F_2, v=Variance_F, cluster=ID_R,
                       data = dat_conservatism, 
                       I2=c("I2q", "ICC"),
                       intervals.type = "LB",
                       model.name = "3 level cs no1w")
summary(model_cs_no1w)
#summary(rerun(model_cs_no1w))

##### 2.2.2. Extreme Variances
dat_conservatism <- dat_conservatism %>% 
  mutate(
    small_var = Variance_F <= quantile(Variance_F, .1, na.rm = TRUE),
    Variance_F_2 = if_else(Variance_F <= quantile(Variance_F, .1, na.rm = TRUE), 
                           quantile(Variance_F, .1, na.rm = TRUE),
                           Variance_F))
#sensitivity analysis 1: outliers excluded 
dat_conservatism_no2 <- dat_conservatism[ which(dat_conservatism$small_var==FALSE), ]
model_cs_no2 <- meta3(y=Fisher, v=Variance_F_2, cluster=ID_R,
                      data = dat_conservatism_no2, 
                      I2=c("I2q", "ICC"),
                      intervals.type = "LB",
                      model.name = "3 level cs no2")
summary(model_cs_no2)
#sensitivity analysis 2: outliers winsorised 
model_cs_no2w <- meta3(y=Fisher, v=Variance_F_2, cluster=ID_R,
                       data = dat_conservatism, 
                       I2=c("I2q", "ICC"),
                       intervals.type = "LB",
                       model.name = "3 level cs no2w")
summary(model_cs_no2w)


##### 2.3. Rally Effects ##### 

##### 2.3.1 Extreme Effect Sizes
confint(model_rf_LBCI)[1] #lowerbound
confint(model_rf_LBCI)[4] #upperbound
# Calculate CI for all observed effect sizes
dat_rally$upperci <- dat_rally$Fisher + 1.96 * sqrt(dat_rally$Variance_F)
dat_rally$lowerci <- dat_rally$Fisher - 1.96 * sqrt(dat_rally$Variance_F)
# Create filter variable
dat_rally$outlier <- dat_rally$upperci < confint(model_rf_LBCI)[1] | dat_rally$lowerci > confint(model_rf_LBCI)[4]
# Count number of outliers:
sum(dat_rally$outlier)
# Make a basic plot: specify that the x-variable is the effect size, and 
# the color and fill of the histogram bars are based on the value of 'outlier'
ggplot(data = dat_rally, aes(x = Fisher, colour = outlier, fill = outlier)) +
  # Add a histogram with transparent bars (alpha = .2)
  geom_histogram(alpha = .2, bins = 50) +
  # Add a vertical line at the pooled effect value (m_re$b[1])
  geom_vline(xintercept = coef(model_rf_LBCI)[1], colour = "grey10", linetype = 2, size = 1) +
  # Apply a black and white theme
  theme_classic2() +
  scale_fill_grey(start = 0.2, end = 0.6) +
  scale_color_grey(start = 0.2, end = 0.6) +
  labs(x = "Fisher's Z Correlation Coefficient", 
       y = "", 
       colour="Outlier", fill="Outlier") 
#sensitivity analysis 1: outliers excluded 
dat_rally_no <- dat_rally[ which(dat_rally$outlier==FALSE), ]
model_rf_no <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                     data = dat_rally_no, 
                     I2=c("I2q", "ICC"),
                     intervals.type = "LB",
                     model.name = "3 level rf no1")
summary(model_rf_no)
#sensitivity analysis 2: outliers winsorised 
#set extremely low ES's to lower bound
#set extremely high ES's to higher bound
dat_rally <- dat_rally %>% 
  mutate(
    F_1 = if_else(dat_rally$upperci < confint(model_rf_LBCI)[1], 
                  confint(model_rf_LBCI)[1],
                  Fisher),
    F_2 = if_else(dat_rally$lowerci > confint(model_rf_LBCI)[4], 
                  confint(model_rf_LBCI)[4],
                  F_1))
model_rf_no1w <- meta3(y=F_2, v=Variance_F, cluster=ID_R,
                       data = dat_rally, 
                       I2=c("I2q", "ICC"),
                       intervals.type = "LB",
                       model.name = "3 level rf no1w")
summary(model_rf_no1w)

##### 1.3.2. Extreme Variances
dat_rally <- dat_rally %>% 
  mutate(
    small_var = Variance_F <= quantile(Variance_F, .1, na.rm = TRUE),
    Variance_F_2 = if_else(Variance_F <= quantile(Variance_F, .1, na.rm = TRUE), 
                           quantile(Variance_F, .1, na.rm = TRUE),
                           Variance_F))
#sensitivity analysis 1: outliers excluded 
dat_rally_no2 <- dat_rally[ which(dat_rally$small_var==FALSE), ]
model_rf_no2 <- meta3(y=Fisher, v=Variance_F_2, cluster=ID_R,
                      data = dat_rally_no2, 
                      I2=c("I2q", "ICC"),
                      intervals.type = "LB",
                      model.name = "3 level rf no2")
summary(model_rf_no2)
#sensitivity analysis 2: outliers winsorised 
model_rf_no2w <- meta3(y=Fisher, v=Variance_F_2, cluster=ID_R,
                       data = dat_rally, 
                       I2=c("I2q", "ICC"),
                       intervals.type = "LB",
                       model.name = "3 level rf no2w")
summary(model_rf_no2w)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
##### 3. Sensitivity to study quality (Table C.7) #####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 

##### 3.1. Outgroup Hostility ##### 
dat_outgroup <- dat_outgroup %>% 
  mutate(quality_1 = if_else(DV_Quality == 3, 1, 0),
         quality_2 = if_else(IV_Quality == 1, 1, 0),
         quality_3 = if_else(`NR-AT` == 1, 1, 0),
         quality_4 = if_else(Hypotheses == 0, 0, 1),
         quality_5 = if_else(StudentPop == 1, 0, 1),
         quality_6 = if_else(Preregistration == 1, 1, 0),
         quality_7 = if_else(TypeReport == 1, 1, 0)) %>% 
  mutate(quality = select(., quality_1:quality_7) %>% rowSums(na.rm = TRUE))
format(round(psych::describe(dat_outgroup$quality), 3), nsmall=3)
#sensitivity analysis 1: impact of study quality
og_quality <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                    x=scale(quality, center = T, scale = T),
                    data = dat_outgroup, 
                    model.name = "og quality1")
summary(og_quality)
og_quality_1 <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                    x=quality_1,
                    data = dat_outgroup, 
                    model.name = "og quality1_1")
summary(og_quality_1)
og_quality_2 <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                      x=quality_2,
                      data = dat_outgroup, 
                      model.name = "og quality1_2")
summary(og_quality_2)
og_quality_3 <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                      x=quality_3,
                      data = dat_outgroup, 
                      model.name = "og quality1_3")
summary(og_quality_3)
og_quality_4 <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                      x=quality_4,
                      data = dat_outgroup, 
                      model.name = "og quality1_4")
summary(og_quality_4)
og_quality_5 <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                      x=quality_5,
                      data = dat_outgroup, 
                      model.name = "og quality1_5")
summary(og_quality_5) #significant indicator
og_quality_6 <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                      x=quality_6,
                      data = dat_outgroup, 
                      model.name = "og quality1_6")
summary(og_quality_6) #significance .051
og_quality_7 <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                      x=quality_7,
                      data = dat_outgroup, 
                      model.name = "og quality1_7")
summary(og_quality_7)
#sensitivity analysis 2: impact of effective sample size used
og_samplesize <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                       x=scale(ESS, center = T, scale = T),
                       data = dat_outgroup, 
                       model.name = "og quality2")
summary(og_samplesize)
#sensitivity analysis 3: impact of impact factor
og_if <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
               x=scale(ImpactFactor, center = T, scale = T),
               data = dat_outgroup, 
               model.name = "og quality3")
summary(og_if)

##### 3.2. Conservatism ##### 
dat_conservatism <- dat_conservatism %>% 
  mutate(quality_1 = if_else(DV_Quality == 3, 1, 0),
         quality_2 = if_else(IV_Quality == 1, 1, 0),
         quality_3 = if_else(`NR-AT` == 1, 1, 0),
         quality_4 = if_else(Hypotheses == 0, 0, 1),
         quality_5 = if_else(StudentPop == 1, 0, 1),
         quality_6 = if_else(Preregistration == 1, 1, 0),
         quality_7 = if_else(TypeReport == 1, 1, 0)) %>% 
  mutate(quality = select(., quality_1:quality_7) %>% rowSums(na.rm = TRUE))
format(round(psych::describe(dat_conservatism$quality), 3), nsmall=3)
#sensitivity analysis 1: impact of study quality
cs_quality <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                    x=scale(quality, center = T, scale = T),
                    data = dat_conservatism, 
                    model.name = "cs quality1")
summary(cs_quality)
cs_quality_1 <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                      x=quality_1,
                      data = dat_conservatism, 
                      model.name = "cs quality1_1")
summary(cs_quality_1)
cs_quality_2 <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                      x=quality_2,
                      data = dat_conservatism, 
                      model.name = "cs quality1_2")
summary(cs_quality_2) #significant
cs_quality_3 <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                      x=quality_3,
                      data = dat_conservatism, 
                      model.name = "cs quality1_3")
summary(cs_quality_3)
cs_quality_4 <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                      x=quality_4,
                      data = dat_conservatism, 
                      model.name = "cs quality1_4")
summary(cs_quality_4)
cs_quality_5 <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                      x=quality_5,
                      data = dat_conservatism, 
                      model.name = "cs quality1_5")
summary(cs_quality_5)
cs_quality_6 <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                      x=quality_6,
                      data = dat_conservatism, 
                      model.name = "cs quality1_6")
summary(cs_quality_6)
cs_quality_7 <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                      x=quality_7,
                      data = dat_conservatism, 
                      model.name = "cs quality1_7")
summary(cs_quality_7)
#sensitivity analysis 2: impact of effective sample size used
cs_samplesize <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                       x=scale(ESS, center = T, scale = T),
                       data = dat_conservatism, 
                       model.name = "cs quality2")
summary(cs_samplesize)
#sensitivity analysis 3: impact of impact factor
cs_if <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
               x=scale(ImpactFactor, center = T, scale = T),
               data = dat_conservatism, 
               model.name = "cs quality3")
summary(cs_if)

##### 3.3. Rally Effects ##### 
dat_rally <- dat_rally %>% 
  mutate(quality_1 = if_else(DV_Quality == 3, 1, 0),
         quality_2 = if_else(IV_Quality == 1, 1, 0),
         quality_3 = if_else(`NR-AT` == 1, 1, 0),
         quality_4 = if_else(Hypotheses == 0, 0, 1),
         quality_5 = if_else(StudentPop == 1, 0, 1),
         quality_6 = if_else(Preregistration == 1, 1, 0),
         quality_7 = if_else(TypeReport == 1, 1, 0)) %>% 
  mutate(quality = select(., quality_1:quality_7) %>% rowSums(na.rm = TRUE))
format(round(psych::describe(dat_rally$quality), 3), nsmall = 3)
#sensitivity analysis 1: impact of study quality
rf_quality <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                    x=scale(quality, center = T, scale = T),
                    data = dat_rally, 
                    model.name = "rf quality1")
summary(rf_quality)
rf_quality_1 <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                      x=quality_1,
                      data = dat_rally, 
                      model.name = "rf quality1_1")
summary(rf_quality_1) #significant
rf_quality_2 <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                      x=quality_2,
                      data = dat_rally, 
                      model.name = "rf quality1_2")
summary(rf_quality_2)
rf_quality_3 <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                      x=quality_3,
                      data = dat_rally, 
                      model.name = "rf quality1_3")
summary(rf_quality_3)
rf_quality_4 <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                      x=quality_4,
                      data = dat_rally, 
                      model.name = "rf quality1_4")
summary(rf_quality_4)
rf_quality_5 <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                      x=quality_5,
                      data = dat_rally, 
                      model.name = "rf quality1_5")
summary(rf_quality_5)
rf_quality_6 <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                      x=quality_6,
                      data = dat_rally, 
                      model.name = "rf quality1_6")
summary(rf_quality_6)
rf_quality_7 <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                      x=quality_7,
                      data = dat_rally, 
                      model.name = "rf quality1_7")
summary(rf_quality_7)
#sensitivity analysis 2: impact of effective sample size used
rf_samplesize <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                       x=scale(ESS, center = T, scale = T),
                       data = dat_rally, 
                       model.name = "rf quality2")
summary(rf_samplesize)
#sensitivity analysis 3: impact of impact factor
rf_if <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
               x=scale(ImpactFactor, center = T, scale = T),
               data = dat_rally, 
               model.name = "rf quality3")
summary(rf_if)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
##### 4. Sensitivity to coefficients from multivariate regression (Table C.8) #####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 

##### 4.1. Outgroup Hostility ##### 
dat_outgroup_nr <- dat_outgroup[ which(dat_outgroup$Regression==0), ]
model_og_nr <- meta3(y=Fisher, v=Variance_F_2, cluster=ID_R,
                     data = dat_outgroup_nr, 
                     I2=c("I2q", "ICC"),
                     intervals.type = "LB",
                     model.name = "3 level og no reg coefs")
summary(model_og_nr)
##### 4.2. Conservatism ##### 
dat_conservatism_nr <- dat_conservatism[ which(dat_conservatism$Regression==0), ]
model_cs_nr <- meta3(y=Fisher, v=Variance_F_2, cluster=ID_R,
                     data = dat_conservatism_nr, 
                     I2=c("I2q", "ICC"),
                     intervals.type = "LB",
                     model.name = "3 level cs no reg coefs")
summary(model_cs_nr)
##### 4.3. Rally Effects ##### 
dat_rally_nr <- dat_rally[ which(dat_rally$Regression==0), ]
model_rf_nr <- meta3(y=Fisher, v=Variance_F_2, cluster=ID_R,
                     data = dat_rally_nr, 
                     I2=c("I2q", "ICC"),
                     intervals.type = "LB",
                     model.name = "3 level cs no reg coefs")
summary(model_rf_nr)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
###### B. PUBLICATION BIAS ANALYSES ######
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
##### 1. Detecting Publication Bias (Table C.9) #####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 

##### 1.1. Outgroup Hostility ##### 
summary(as.factor(dat_outgroup$TypeReport))
dat_outgroup$Published <- ifelse(dat_outgroup$TypeReport==1, yes=1, no=0)
summary(as.factor(dat_outgroup$Published))
og_pub<- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                      x=scale(Published),
                      data = dat_outgroup, 
                      model.name = "og pub")
summary(og_pub)
og_ess <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                       x=scale(ESS, center = T, scale = T),
                       data = dat_outgroup, 
                       model.name = "eggers ess")
summary(og_ess)
og_sei<- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                      x=scale((1/SE_F), center = T, scale = T),
                      data = dat_outgroup, 
                      model.name = "eggers sei")
summary(og_sei)
og_var<- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                      x=scale((1/Variance_F), center = T, scale = T),
                      data = dat_outgroup, 
                      model.name = "eggers var")
summary(og_var)


##### 1.2. Conservatism ##### 
summary(as.factor(dat_conservatism$TypeReport))
dat_conservatism$Published <- ifelse(dat_conservatism$TypeReport==1, yes=1, no=0)
summary(as.factor(dat_conservatism$Published))
cs_pub<- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                      x=scale(Published),
                      data = dat_conservatism, 
                      model.name = "cs pub")
summary(cs_pub)
cs_ess <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                       x=scale(ESS, center = T, scale = T),
                       data = dat_conservatism, 
                       model.name = "eggers ess")
summary(cs_ess)
cs_sei<- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                      x=scale((1/SE_F), center = T, scale = T),
                      data = dat_conservatism, 
                      model.name = "eggers sei")
summary(cs_sei)
cs_var<- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                      x=scale((1/Variance_F), center = T, scale = T),
                      data = dat_conservatism, 
                      model.name = "eggers var")
summary(cs_var)


##### 1.3. Rally Effects ##### 
summary(as.factor(dat_rally$TypeReport))
dat_rally$Published <- ifelse(dat_rally$TypeReport==1, yes=1, no=0)
summary(as.factor(dat_rally$Published))
rf_pub<- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                      x=scale(Published),
                      data = dat_rally, 
                      model.name = "cs pub")
summary(rf_pub)
rf_ess <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                       x=scale(ESS, center = T, scale = T),
                       data = dat_rally, 
                       model.name = "eggers ess")
summary(rf_ess)
rf_sei<- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                      x=scale((1/SE_F), center = T, scale = T),
                      data = dat_rally, 
                      model.name = "eggers sei")
summary(rf_sei)
rf_var<- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                      x=scale((1/Variance_F), center = T, scale = T),
                      data = dat_rally, 
                      model.name = "eggers var")
summary(rf_var)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
##### 2. Assessing Publication Bias (Figure C.1 and Table C.10) #####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 

#### 2.1. Outgroup Hostility #### 

#trim and fill
og_funnel = data.frame(dat_outgroup[, c("Fisher", "SE_F")])
og_funnel <- viz_funnel(og_funnel, y_axis = "se", method = "DL", 
                        contours = T, sig_contours = F, addev_contours = F,
                        trim_and_fill = T, trim_and_fill_side = "left",
                        egger = F,
                        xlab = "Fisher's Z Correlation Coefficient \n",
                        x_breaks = (c(-1, -0.5, 0, 0.5, 1)))
og_funnel <- og_funnel + mytheme + ggtitle('A. Outgroup Hostility')
#PET
og_PET <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                x=scale(SE_F),
                data = dat_outgroup,
                #intervals.type = "LB",
                model.name = "cs PET")
summary(og_PET)
#PEESE
og_PEESE <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                x=scale(Variance_F),
                data = dat_outgroup,
                #intervals.type = "LB",
                model.name = "cs PEESE")
summary(og_PEESE)


#### 2.2. Conservative Shift #### 

#trim and fill
cs_funnel = data.frame(dat_conservatism[, c("Fisher", "SE_F")])
cs_funnel <- viz_funnel(cs_funnel, y_axis = "se", method = "DL", 
                        contours = T, sig_contours = F, addev_contours = F,
                        trim_and_fill = T, trim_and_fill_side = "left",
                        egger = F,
                        xlab = "Fisher's Z Correlation Coefficient \n",
                        x_breaks = (c(-1, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1)))
cs_funnel <- cs_funnel + mytheme + ggtitle('B. Conservative Shift')
#PET
cs_PET <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                x=scale(SE_F),
                data = dat_conservatism,
                #intervals.type = "LB",
                model.name = "cs PET")
summary(cs_PET)
#PEESE
cs_PEESE <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                  x=scale(Variance_F),
                  data = dat_conservatism,
                  #intervals.type = "LB",
                  model.name = "cs PEESE")
summary(cs_PEESE)


#### 2.3. Rally Effects #### 

#trim and fill
rf_funnel = data.frame(dat_rally[, c("Fisher", "SE_F")])
rf_funnel <- viz_funnel(rf_funnel, y_axis = "se", method = "DL", 
                        contours = T, sig_contours = F, addev_contours = F,
                        trim_and_fill = T, trim_and_fill_side = "left",
                        egger = F,
                        xlab = "Fisher's Z Correlation Coefficient",
                        x_breaks = (c(-1, -0.5, 0, 0.5, 1)))
rf_funnel <- rf_funnel + mytheme + ggtitle('C. Rally-Around-The-Flag') 
#PET
rf_PET <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                x=scale(SE_F),
                data = dat_rally,
                #intervals.type = "LB",
                model.name = "cs PET")
summary(rf_PET)
#PEESE
rf_PEESE <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                  x=scale(Variance_F),
                  data = dat_rally,
                  #intervals.type = "LB",
                  model.name = "cs PET")
summary(rf_PEESE)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
################################ FIGURE S1 ###############################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
layout_matrix <- matrix(c(1, 1, 2, 2, 4, 3, 3, 4), nrow = 2, byrow = TRUE)
grid.arrange(og_funnel, cs_funnel, rf_funnel, layout_matrix = layout_matrix)
jpeg("Figures/FigS1.jpeg", width = 12, height = 11, units = 'in', res = 500)
grid.arrange(og_funnel, cs_funnel, rf_funnel, 
             layout_matrix = layout_matrix)
dev.off() #safe picture in high resolution.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
