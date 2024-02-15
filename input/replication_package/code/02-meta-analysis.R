# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
# PROJECT:    How Terrorism Does (and Does Not) Affect Citizens' 
#             Political Attitudes: A Meta-Analysis. 
# AUTHOR:     Amelie Godefroidt
# CONTACT:    amelie.godefroidt@ntnu.no
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
# This R file contains the code necessary to replicate all main results,
# including the descriptive statistics, meta-analysis, and meta-regressions.
# 
# Last successful replication: 01.11.2021
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 

# first things first: set your settings.

rm(list=ls())
setwd(".../Godefroidt_2021_ReplicationFiles/") #set your directory
options("scipen"=100, "digits"=8) #set digits
getwd()


###### Install and load packages ###### 

# ipak function: install and load multiple R packages.
# check to see if packages are installed. Install them if they are not, then load them into the R session.

ipak <- function(pkg){  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
if(length(new.pkg)) install.packages(new.pkg, dependencies = TRUE)
sapply(pkg, require, character.only = TRUE)
}

packages <- c("readxl", "haven", "psych", "dplyr", "DescTools", "broom", "cowplot",
              "pwr", "esc", "metaSEM", "metafor", "metaviz", "PerformanceAnalytics", 
              "ggplot2", "ggpubr", "wesanderson", "stringr", "dotwhisker") 
ipak(packages)


###### Import and split data ###### 

dat <- read_excel("data/02-metaanalysis-data.xlsx", sheet = "Data", 
                  col_names = TRUE, skip = 3)

### dataframe of studies
dat_studies <- dat[c(1:87)] #studies.
dat_studies <- group_by(dat_studies, ID_RS) %>% slice(1)
sum(dat_studies$ESS, na.rm = TRUE) #number of unique respondents

### dataframe of reports
dat_reports <- dat[c(1:87)] #reports/publications.
dat_reports <- group_by(dat_reports, ID_R) %>% slice(1)


### set wished layout settings
mytheme <- theme(plot.title = element_text(face = "bold", size = (22), colour = "black"), 
                 axis.text = element_text(size = (18), colour = "black"),
                 axis.ticks.length = unit(0.5, "cm"),
                 axis.title.y = element_text(size = (22), colour = "black"),
                 axis.title.x = element_text(size = (22), colour = "white"))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
###### A. DESCRIPTIVES: THE FIELD OF EMPIRICAL TERRORISM STUDIES ######
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 

## mean ES's per report (reported in footnote 6)
format(round(psych::describe(aggregate(Correlation ~ ID_R, data = dat, FUN = length)),3), nsmall=3)


####### Year of study and publication, and country of study ######

# year of study

psych::describe(dat_studies$StudyYear)
# Figure 2A
StudyYear <- ggplot(dat_studies, aes(x=StudyYear)) + 
  geom_histogram(binwidth=1,
                 color="#000000", fill="#939598") +
  labs(x="to fill some space", y = "", title = "Year of Data Collection") +
  scale_y_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45)) +
  scale_x_continuous(breaks = c(1985,1990,1995,2000,2005,2010, 2015,2020)) +
  theme_classic() #histogram year of study 
StudyYear <- StudyYear + geom_vline(xintercept=2001, linetype="dashed", size = 1) #add 9/11 line
StudyYear <- StudyYear + geom_vline(xintercept=2004, linetype="dashed", size = 1) #add line
StudyYear <- StudyYear + geom_vline(xintercept=2015, linetype="dashed", size = 1) #add ISIS line
StudyYear <- StudyYear + geom_label(data = dat_studies, 
                                    aes(x=2001, y=27, label = "9/11"),
                                    size=8) #add 9/11 label
StudyYear <- StudyYear + geom_label(data = dat_studies, 
                                    aes(x=2007, y=35, label = "Madrid & \n Israel-Palestine"),
                                    size=6) #add label
StudyYear <- StudyYear + geom_label(data = dat_studies, 
                                    aes(x=2014, y=44, label = "IS Attacks"),
                                    size=8) #add ISIS label
StudyYear <- StudyYear + mytheme


# year of publication

describe(dat_reports$Year) 
# Figure 2B
PubYear <- ggplot(dat_reports, aes(x=Year)) + 
  geom_histogram(binwidth=1, 
                 color="#000000", fill="#939598") +
  labs(title ="Year of Publication", x = "to fill some space", y = "") +
  scale_y_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45)) +
  scale_x_continuous(breaks = c(1985,1990,1995,2000,2005,2010, 2015,2020)) +
  theme_classic() #histogram year of study 
PubYear <- PubYear + geom_vline(xintercept=2001, linetype="dashed", size = 1) #add 9/11 line
PubYear <- PubYear + geom_vline(xintercept=2015, linetype="dashed", size = 1) #add ISIS France line
PubYear <- PubYear + geom_label(data = dat_reports,
                                aes(x=2001, y=27, label = "9/11"),
                                size=8) #add 9/11 label
PubYear <- PubYear + geom_label(data = dat_reports, 
                                aes(x=2014, y=44, label = "IS Attacks"),
                                size=8) #add ISIS label
PubYear <- PubYear + mytheme


####### Terrorism measures, incl. cognitive and affective reactions ######

# IV measurements used 
summary(as.factor(dat$IV_Code))

#wrangling
dat$TerrorMeasure[dat$IV_Code==1] <- "Acts of violence"
dat$TerrorMeasure[dat$IV_Code==2] <- "Acts of violence"
dat$TerrorMeasure[dat$IV_Code==3] <- "Threat"
dat$TerrorMeasure[dat$IV_Code==4] <- "Threat"
dat$TerrorMeasure[dat$IV_Code==5] <- "Other"
dat$TerrorMeasure[dat$IV_Code==6] <- "Reported exposure"
dat$TerrorMeasure[dat$IV_Code==7] <- "Reported exposure"
dat$TerrorMeasure[dat$IV_Code==8] <- "Reported exposure"
dat$TerrorMeasure[dat$IV_Code==9] <- "Threat"
dat$TerrorMeasure[dat$IV_Code==10] <- "Threat"
dat$TerrorMeasure[dat$IV_Code==11] <- "Fear"
dat$TerrorMeasure[dat$IV_Code==12] <- "Anger"
dat$TerrorMeasure[dat$IV_Code==13] <- "Other"
dat$TerrorMeasure[dat$IV_Code==14] <- "Other"
dat$TerrorMeasure[dat$IV_Code==15] <- "Other"
dat$TerrorMeasure[dat$IV_Code==16] <- "Other"

dat$TerrorMeasure <- as.factor(dat$TerrorMeasure)
table(dat$TerrorMeasure)
prop.table(table(dat$TerrorMeasure))
table(dat$ExactAttack) #which attacks are studied?
prop.table(table(dat$ExactAttack)) #which attacks are studied?

#visualization: Figure 2C
tmdat <- data.frame(prop.table(table(dat$TerrorMeasure)))
TerrorMeasures <- ggplot(tmdat, aes(x = reorder(Var1, -Freq), y = Freq)) +
  geom_bar(stat="identity", width=0.7,
           color="#000000", fill="#939598") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  scale_y_continuous(labels = scales::percent,
                     limits = c(0,.6),
                     breaks = c(0,.1,.2,.3,.4,.5,.6)) +
  labs(title = "Terrorism Measures", x="", y = "") +
  theme_classic() 
TerrorMeasures <- TerrorMeasures + mytheme


# type of terrorism
summary(as.factor(dat$Type))

#wrangling
dat$TerrorType[dat$Type==0] <- "No ideology"
dat$TerrorType[dat$Type==1] <- "Islamist"
dat$TerrorType[dat$Type==2] <- "Extreme right"
dat$TerrorType[dat$Type==3] <- "Other"
dat$TerrorType[dat$Type==4] <- "Other"
dat$TerrorType[dat$Type==5] <- "State terror"
dat$TerrorType[dat$Type==6] <- "Other"
dat$TerrorType <- as.factor(dat$TerrorType)
summary(dat$TerrorType)
prop.table(table(dat$TerrorType))

#visualization: Figure 2D
ttdat <- data.frame(prop.table(table(dat$TerrorType)))
TerrorType <- ggplot(ttdat, aes(x = reorder(Var1, -Freq), y = Freq)) +
  geom_bar(stat="identity", width=0.7,
           color="#000000", fill="#939598") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  scale_y_continuous(labels = scales::percent,
                     breaks = c(0,.1,.2,.3,.4,.5,.6)) +
  labs(title = "Type of Terrorism", x = "", y = "") +
  theme_classic() 
TerrorType <- TerrorType + mytheme


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
################################# FIGURE 2 ###############################
jpeg("Figures/Fig2.jpeg", width = 15, height = 16, units = 'in', res = 300)
ggarrange(StudyYear, PubYear, 
          TerrorMeasures, TerrorType,
          labels = c("(A)", "(B)", "(C)", "(D)"),
          font.label = list(size = 22, color = "black"),
          ncol = 2, nrow = 2) #arrange all plots
dev.off() #safe picture in high resolution.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 


####### Other study characteristics ######

# country of study
dat_studies$Country <- as.factor(dat_studies$Country)
table(dat_studies$Country)
prop.table(table(dat_studies$Country))
length(unique(dat_studies$Country))
#Europe: All European regions, based on UN geoscheme
Europe = 1+ #Austria
  15+ #Belgium 1
  1+ #Bosnia 2
  1+ #Crimea 3
  42+ #Denmark 4
  7+ #EU
  1+ #Finland 5
  14+ #France 6
  22+ #Germany 7
  1+ #Northern Ireland (study mixed with Israel, no info on separate samples)
  2+ #Italy 8
  9+ #Netherlands 9
  1+ #Northern Ireland
  7+ #Norway 10
  5+ #Poland 11
  14+ #Spain 12
  1+ #Sweden 13
  1+ #Switzerland 14
  17+ #UK 15
  1 #US, UK, and Australia (no info on separate samples))
Europe
Europe/326
#Non-western
1+ # Colombia
  1+ #Egypt
  1+ #Morocco
  2+ #Nigeria
  1+ #Pakistan
  1 #South-Africa
7/326

# research design
summary(as.factor(dat_studies$TypeStudy))
dat_studies$Design[dat_studies$TypeStudy==1] <- "Experiment"
dat_studies$Design[dat_studies$TypeStudy==2] <- "QuasiExperiment"
dat_studies$Design[dat_studies$TypeStudy==3] <- "Correlation"
dat_studies$Design[dat_studies$TypeStudy==4] <- "Longitudinal"
dat_studies$Design <- as.factor(dat_studies$Design)
table(dat_studies$Design)
prop.table(table(dat_studies$Design))

# sampling
dat_studies$Sample[dat_studies$GeneralPop==1] <- "General population"
dat_studies$Sample[dat_studies$StudentPop==1] <- "Student population"
dat_studies$Sample[is.na(dat_studies$Sample)] <- "Other"
dat_studies$Sample <- as.factor(dat_studies$Sample)
table(dat_studies$Sample)
prop.table(table(dat_studies$Sample))

# minority/majority distribution
table(as.factor(dat_studies$EthCat))
dat_studies$Majority[dat_studies$EthCat==1] <- "Minority"
dat_studies$Majority[dat_studies$EthCat==2] <- "Majority"
dat_studies$Majority[dat_studies$EthCat==3] <- "Mixed"
dat_studies$Majority[dat_studies$EthCat==4] <- "Unknown"
dat_studies$Majority <- as.factor(dat_studies$Majority)
table(dat_studies$Majority)
prop.table(table(dat_studies$Majority))

# mean age
format(round(psych::describe(dat_studies$AgeMean), 3), nsmall = 3)

# gender distribution
format(round(psych::describe(dat_studies$FemaleP), 3), nsmall = 3)

# effective sample size
format(round(psych::describe(dat_studies$ESS), 3), nsmall = 3)

# preregistration
summary(as.factor(dat_studies$Preregistration))
11/(315+11)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
######## B. META-ANALYSIS: for each of the hypotheses ######## 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 

# Data preparation: Split dataset in 3 based on outcome type
# 1) outgroup dataset
dat_outgroup <- dat[dat$PA_Category == 10 |
                      dat$PA_Category == 99, ]
length(unique(dat_outgroup$ID_R))#number of manuscripts (N_j)
# 2) conservative shift dataset
dat_conservatism <- dat[dat$PA_Category == 5 |
                          dat$PA_Category == 6 |
                          dat$PA_Category == 7 |
                          dat$PA_Category == 8 |
                          dat$PA_Category == 9 |
                          dat$PA_Category == 11, ]
length(unique(dat_conservatism$ID_R))#number of manuscripts (N_j)
# 3) rally-around-the-flag dataset
dat_rally <- dat[dat$PA_Category == 1 |
                   dat$PA_Category == 2 |
                   dat$PA_Category == 3 |
                   dat$PA_Category == 4, ]
length(unique(dat_rally$ID_R))#number of manuscripts (N_j)
# double check (total number of effect sizes; N_k)
count(dat_outgroup) + count(dat_conservatism) + count(dat_rally)


##### 1. Fit empty 3-level model using the meta3 function #####

### To replicate Table 2: call summary() function for each of outcome cluster
      ## coef(model...)[1] = Fishers' Z correlation coefficients (Zr)
### To replicate the correlation coefficients mentioned in the text, 
### call the FisherZinv() function for each outcome cluster.

# Outgroup hostility (Table 2)
set.seed(1234)
model_og_LBCI <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                       data = dat_outgroup, 
                       I2=c("I2q", "ICC"),
                       intervals.type = "LB",
                       model.name = "3 level outgroup LB")
summary(model_og_LBCI)#results presented in Table 2, as well as
#throughout the Supplementary Information (esp. Table B.3).
FisherZInv(coef(model_og_LBCI)[1])#rho reported in text


# Conservative shift (Table 2)
set.seed(1234)
model_cs_LBCI <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                       data = dat_conservatism, 
                       I2=c("I2q", "ICC"),
                       intervals.type = "LB",
                       model.name = "3 level conservatism LB")
summary(model_cs_LBCI)#results presented in Table 2, as well as
#throughout the Supplementary Information (esp. Table B.3).
FisherZInv(coef(model_cs_LBCI)[1])#rho reported in text


# Rally-around-the-flag (Table 2)
set.seed(1234)
model_rf_LBCI <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                       data = dat_rally, 
                       I2=c("I2q", "ICC"),
                       intervals.type = "LB",
                       model.name = "3 level rally LB")
summary(model_rf_LBCI)#results presented in Table 2, as well as
#throughout the Supplementary Information (esp. Table B.3).
FisherZInv(coef(model_rf_LBCI)[1])#rho reported in text


## Substantive interpretation --------------------------------------------------
## i.e., comparison to 2016 ANES Feeling Thermometers
anes2016 <- read_dta(file = "data/anes_timeseries_2016.dta")
#select variables to summarize
anes2016 <- anes2016 %>%
  select(immterm = V162313, 
         reptherm = V161096) 
#format(round(psych::describe(anes2016),3), nsmall=3): counts NA's
#recode NA's (!)
anes2016$immterm <- ifelse(anes2016$immterm %in% -9:-1, NA, 
                           identity(anes2016$immterm)) 
anes2016$reptherm <- ifelse(anes2016$reptherm %in% -99:-1, NA, 
                            identity(anes2016$reptherm))
#find correct SDs to use in calculations based on Paluck et al. (2021)
format(round(psych::describe(anes2016),3), nsmall=3) 
#the effect of terror on the feeling therm. toward immigrants. 6.9
cohens_d(r = FisherZInv(coef(model_og_LBCI)[1]))*27.285 
#the effect of terror on the feeling therm. toward Rep. party: 7.2
cohens_d(r = FisherZInv(coef(model_cs_LBCI)[1]))*27.330 


## Practical implication -------------------------------------------------------
## i.e., what sample size is recommended based on this meta-analysis?
# Outgroup hostility
pwr.r.test(r = FisherZInv(coef(model_og_LBCI)[1]), sig.level = 0.05, power = 0.8) 
# Conservative shift
pwr.r.test(r = FisherZInv(coef(model_cs_LBCI)[1]), sig.level = 0.05, power = 0.8) 
# Rally-around-the-flag
pwr.r.test(r = FisherZInv(coef(model_rf_LBCI)[1]), sig.level = 0.05, power = 0.8) 



##### 2. Test variance components #####

### Note: These results are briefly mentioned in footnote 13, but 
### full results are discussed and presented in SI B.2.3.

# Outgroup hostility

#tau^2_3=0 (reported in Table B.2)
set.seed(1234)
model_og_t3 <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                      data = dat_outgroup, 
                      RE3.constraints = 0, # constraining between-study level
                      model.name = "Tau2_3 eq 0",
                      intervals.type = "LB")
summary(model_og_t3)#this equals a conventional random-effects model
anova(model_og_LBCI, model_og_t3)

#tau^2_2=0 (reported in Table B.2)
set.seed(1234)
model_og_t2 <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                     data = dat_outgroup,
                     RE2.constraints = 0, # constraining within-study level
                     model.name = "Tau2_2 eq 0",
                     intervals.type = "LB")
summary(model_og_t2)
anova(model_og_LBCI, model_og_t2)


# Conservative shift

#tau^2_3=0 (reported in Table B.2)
set.seed(1234)
model_cs_t3 <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                     data = dat_conservatism, 
                     RE3.constraints = 0, # constraining between-study level
                     model.name = "Tau2_3 eq 0",
                     intervals.type = "LB")
summary(model_cs_t3)
anova(model_cs_LBCI, model_cs_t3)

#tau^2_2=0 (reported in Table B.2)
set.seed(1234)
model_cs_t2 <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                     data = dat_conservatism,
                     RE2.constraints = 0, # constraining within-study level
                     model.name = "Tau2_2 eq 0",
                     intervals.type = "LB")
summary(model_cs_t2)
anova(model_cs_LBCI, model_cs_t2)


# Rally-around-the-flag

#tau^2_3=0 (reported in Table B.2)
set.seed(1234)
model_rf_t3 <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                     data = dat_rally, 
                     RE3.constraints = 0, # constraining between-study level
                     model.name = "Tau2_3 eq 0",
                     intervals.type = "LB")
summary(model_rf_t3)
anova(model_rf_LBCI, model_rf_t3)

#tau^2_2=0 (reported in Table B.2)
set.seed(1234)
model_rf_t2 <- meta3(y=Fisher, v=Variance_F, cluster=ID_R,
                     data = dat_rally,
                     RE2.constraints = 0, # constraining within-study level
                     model.name = "Tau2_2 eq 0",
                     intervals.type = "LB")
summary(model_rf_t2)
anova(model_rf_LBCI, model_rf_t2)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
######## C. MODERATOR ANALYSES: for each of the hypotheses ######## 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 

### Some of the results below are highlighted in the paper. 
      ## To replicate the correlation coefficients reported in the text, call the relevant FisherZinv() function.
      ## To replicate the b estimates reported in the text, call the relevant summary(model...) function
      ## To replicate Figure 3, go to line 2026
      ## To replicate Figure 4, go to line 2475

### Full results are reported in the Supplementary Information (SI)
      ## To replicate Tables C.1, C.2 and C.3:
         # 1) go to the relevant moderator and 
         # 2) call summary() for the overall effect size (b) and 95% LBCI values
         # 3) run the line of code to calculate the R squared (sum score)
         # 4) run the anova() function to get the LRT statistic
         # 5) run the lines of code to calculate the adjusted p-values to get the reported subscripts
         # Note: the number of effect sizes per category (k) can be found 
         # via the summary() function for each of the relevant moderators.


##### 1. Outgroup Hostility Hypothesis ##### 

##### 1.1. Terrorism-Related Moderators ##### 

### Type/Ideology of Terrorism

#check data
summary(as.factor(dat_outgroup$Type))
summary(as.factor(dat_outgroup$TerrorType))
#recode data
dat_outgroup$TerrorType2[dat_outgroup$Type==0] <- "No ideology"
dat_outgroup$TerrorType2[dat_outgroup$Type==1] <- "Islamist"
dat_outgroup$TerrorType2[is.na(dat_outgroup$TerrorType2)] <- "Other"
dat_outgroup$TerrorType2 <- as.factor(dat_outgroup$TerrorType2)
table(dat_outgroup$TerrorType2)
#convert into dummy variables
dat_outgroup$Islam <- ifelse(dat_outgroup$TerrorType2=="Islamist", yes=1, no=0)
dat_outgroup$No_Ideology <- ifelse(dat_outgroup$TerrorType2=="No ideology", yes=1, no=0)
dat_outgroup$Other_Ideology <- ifelse(dat_outgroup$TerrorType2=="Other", yes=1, no=0)

#fit model with effect coding
set.seed(1234)
model_og_type <- meta3(y=Fisher, v=Variance_F,
                       intercept.constraints = 0,
                       x=cbind(Islam, Other_Ideology, No_Ideology),
                       cluster=ID_R, 
                       data = dat_outgroup,
                       intervals.type = "LB",
                       model.name = "og ideology effect")
summary(model_og_type) #(mainly reported in SI Table C.1)
FisherZInv(coef(model_og_type)[1]) #Islam (reported in text)
FisherZInv(coef(model_og_type)[2]) #No Islam (reported in text)
#explained variance (only reported in SI Table C.1)
summary(model_og_type)$R2.values[3] + summary(model_og_type)$R2.values[6]
#model fit (only reported SI Table C.1)
anova(model_og_type, model_og_LBCI)

#fit models with dummy coding
set.seed(1234)
model_og_type_d1 <- meta3(y=Fisher, v=Variance_F,
                               x=cbind(Islam, Other_Ideology),
                               cluster=ID_R,
                               data = dat_outgroup,
                               model.name = "og target sim dummy 1")
summary(model_og_type_d1)
set.seed(1234)
model_og_type_d2 <- meta3(y=Fisher, v=Variance_F,
                               x=cbind(Other_Ideology, No_Ideology),
                               cluster=ID_R,
                               data = dat_outgroup,
                               model.name = "og target sim dummy 1")
summary(model_og_type_d2)
#adjust p-values for multiple comparisons
p <- cbind(summary(model_og_type_d1)$coefficients$`Pr(>|z|)`[2],#Islam vs. no ideology
           summary(model_og_type_d1)$coefficients$`Pr(>|z|)`[3],#Other ideology vs. no ideology
           summary(model_og_type_d2)$coefficients$`Pr(>|z|)`[2])#Other ideology vs. Islam
p.adjust(p, "BH")


### Terrorism Measure/Manipulation

#check data
summary(as.factor(dat_outgroup$TerrorMeasure))
#convert characters into dummy variables
dat_outgroup$Exposure <- ifelse(dat_outgroup$TerrorMeasure=="Acts of violence", yes=1, no=0)
dat_outgroup$R_Exposure <- ifelse(dat_outgroup$TerrorMeasure=="Reported exposure", yes=1, no=0)
dat_outgroup$Threat <- ifelse(dat_outgroup$TerrorMeasure=="Threat", yes=1, no=0)
table(as.factor(dat_outgroup$Threat))
dat_outgroup$Emotions <- ifelse(dat_outgroup$TerrorMeasure=="Fear" | dat_outgroup$TerrorMeasure=="Anger", yes=1, no=0)
dat_outgroup$Fear <- ifelse(dat_outgroup$TerrorMeasure=="Fear", yes=1, no=0)
dat_outgroup$Anger <- ifelse(dat_outgroup$TerrorMeasure=="Anger", yes=1, no=0)
dat_outgroup$Other_IV <- ifelse(dat_outgroup$TerrorMeasure=="Other", yes=1, no=0)

#fit model with effect coding
set.seed(1234)
model_og_ivm <- meta3(y=Fisher, v=Variance_F,
                      intercept.constraints = 0,
                      x=cbind(Exposure, R_Exposure, Threat, Emotions, Other_IV),
                      cluster=ID_R, 
                      data = dat_outgroup,
                      intervals.type = "LB",
                      model.name = "og iv measure")
summary(model_og_ivm)#(mainly reported in SI Table C.1)
FisherZInv(coef(model_og_ivm)[1]) #Objective exposure (reported in text)
#explained variance (only reported in SI Table C.1)
summary(model_og_ivm)$R2.values[3] + summary(model_og_ivm)$R2.values[6]
#model fit (only reported in SI Table C.1)
anova(model_og_ivm, model_og_LBCI)
#test anger and fear as well (as pre-registered)
set.seed(1234)
model_og_ivm2 <- meta3(y=Fisher, v=Variance_F,
                      intercept.constraints = 0,
                      x=cbind(Exposure, R_Exposure, Threat, Fear, Anger, Other_IV),
                      cluster=ID_R, 
                      data = dat_outgroup,
                      intervals.type = "LB",
                      model.name = "og iv measure2")
summary(model_og_ivm2) #(only reported in SI)

#fit models with dummy coding
set.seed(1234)
model_og_ivm_d1 <- meta3(y=Fisher, v=Variance_F,
                      x=cbind(R_Exposure, Threat, Emotions, Other_IV),
                      cluster=ID_R, 
                      data = dat_outgroup,
                      model.name = "og iv measure d1")
summary(model_og_ivm_d1)
set.seed(1234)
model_og_ivm_d2 <- meta3(y=Fisher, v=Variance_F,
                         x=cbind(Exposure, Threat, Emotions, Other_IV),
                         cluster=ID_R, 
                         data = dat_outgroup,
                         model.name = "og iv measure d2")
summary(model_og_ivm_d2)
set.seed(1234)
model_og_ivm_d3 <- meta3(y=Fisher, v=Variance_F,
                         x=cbind(Exposure, R_Exposure, Emotions, Other_IV),
                         cluster=ID_R, 
                         data = dat_outgroup,
                         model.name = "og iv measure d3")
summary(model_og_ivm_d3)
set.seed(1234)
model_og_ivm_d4 <- meta3(y=Fisher, v=Variance_F,
                         x=cbind(Exposure, R_Exposure, Threat, Other_IV),
                         cluster=ID_R, 
                         data = dat_outgroup,
                         model.name = "og iv measure d4")
summary(model_og_ivm_d4)
#adjust p-values for multiple comparisons
p <- cbind(0.0650055, 0.0000001914305559, 0.0000203826870788, 0.0019952, #exposure vs. other measures
           0.0000127213971284, 0.0020052, 0.0656952, #self-reported exposure vs. other measures (excl. exp.)
           0.3485339, 0.0080471, #threat vs. emotions, threat vs. other 
           0.1110643) #emotions vs. other 
p.adjust(p, "BH")



###### 1.2. Outcome-Related Moderators ##### 

### Target Group of Outgroup Hostility

#check data
summary(as.factor(dat_outgroup$OA_Target_A))
#recode data
dat_outgroup$TargetType[dat_outgroup$OA_Target_A==3] <- "Religious Outgroup"
dat_outgroup$TargetType[dat_outgroup$OA_Target_A==2] <- "Immigrants"
dat_outgroup$TargetType[dat_outgroup$OA_Target_A==1] <- "Other Outgroup"
dat_outgroup$TargetType[dat_outgroup$OA_Target_A==4] <- "Other Outgroup"
dat_outgroup$TargetType[dat_outgroup$OA_Target_A==5] <- "Other Outgroup"
dat_outgroup$TargetType <- as.factor(dat_outgroup$TargetType)
table(dat_outgroup$TargetType)
#convert into dummy variables
dat_outgroup$ReligiousOG <- ifelse(dat_outgroup$TargetType=="Religious Outgroup", yes=1, no=0)
dat_outgroup$ImmigrantOG <- ifelse(dat_outgroup$TargetType=="Immigrants", yes=1, no=0)
dat_outgroup$OtherOG <- ifelse(dat_outgroup$TargetType=="Other Outgroup", yes=1, no=0)

#fit model with effect coding
set.seed(1234)
model_og_target <- meta3(y=Fisher, v=Variance_F,
                         intercept.constraints = 0,
                         x=cbind(ReligiousOG, ImmigrantOG, OtherOG),
                         cluster=ID_R,
                         data = dat_outgroup, 
                         intervals.type = "LB",
                         model.name = "og target")
summary(model_og_target) #(mainly reported in SI Table C.1)
FisherZInv(coef(model_og_target)[1]) #Religious (reported in text)
FisherZInv(coef(model_og_target)[2]) #Migrant (reported in text)
FisherZInv(coef(model_og_target)[3]) #Other (reported in text)
#explained variance (only reported in SI Table C.1)
summary(model_og_target)$R2.values[3] + summary(model_og_target)$R2.values[6]
#model fit (only reported in SI Table C.1)
anova(model_og_target, model_og_LBCI)

#fit model with dummy coding
set.seed(1234)
model_og_target_d1 <- meta3(y=Fisher, v=Variance_F,
                          x=cbind(ReligiousOG, ImmigrantOG),
                          cluster=ID_R,
                          data = dat_outgroup,
                          model.name = "og target d1")
summary(model_og_target_d1)
set.seed(1234)
model_og_target_d2 <- meta3(y=Fisher, v=Variance_F,
                           x=cbind(ImmigrantOG, OtherOG),
                           cluster=ID_R,
                           data = dat_outgroup,
                           model.name = "og target d2")
summary(model_og_target_d2)
#adjust p-values for multiple comparisons
p <- cbind(summary(model_og_target_d1)$coefficients$`Pr(>|z|)`[2], #Religious vs. Other
           summary(model_og_target_d1)$coefficients$`Pr(>|z|)`[3], #Immigrant vs. Other
           summary(model_og_target_d2)$coefficients$`Pr(>|z|)`[2]) #Religious vs. Immigrant
p.adjust(p, "BH")


### Target Group Similarity

#check data
summary(as.factor(dat_outgroup$OA_Similarity))
#recode data
dat_outgroup$TargetSim[dat_outgroup$OA_Similarity==2] <- "A lot"
dat_outgroup$TargetSim[dat_outgroup$OA_Similarity==1] <- "A bit"
dat_outgroup$TargetSim[dat_outgroup$OA_Similarity==0] <- "No"
dat_outgroup$TargetSim <- as.factor(dat_outgroup$TargetSim)
table(dat_outgroup$TargetSim)
#convert into dummy variables
dat_outgroup$Lot <- ifelse(dat_outgroup$TargetSim=="A lot", yes=1, no=0)
dat_outgroup$Bit <- ifelse(dat_outgroup$TargetSim=="A bit", yes=1, no=0)
dat_outgroup$No <- ifelse(dat_outgroup$TargetSim=="No", yes=1, no=0)

#fit model with effect coding
set.seed(1234)
model_og_targetsim <- meta3(y=Fisher, v=Variance_F,
                            intercept.constraints = 0,
                            x=cbind(Lot, Bit, No),
                            cluster=ID_R,
                            data = dat_outgroup, 
                            intervals.type = "LB",
                            model.name = "og target sim")
summary(model_og_targetsim) #(mainly reported in SI Table C.1)
FisherZInv(coef(model_og_targetsim)[1]) #Hard association (reported in text)
FisherZInv(coef(model_og_targetsim)[2]) #Soft association (reported in text)
FisherZInv(coef(model_og_targetsim)[3]) #No association (reported in text)
#explained variance (only reported in SI Table C.1)
summary(model_og_targetsim)$R2.values[3] + summary(model_og_targetsim)$R2.values[6]
#model fit (only reported in SI Table C.1)
anova(model_og_targetsim, model_og_LBCI)

#fit model with dummy coding
set.seed(1234)
model_og_targetsim_d1 <- meta3(y=Fisher, v=Variance_F,
                            x=cbind(Lot, Bit),
                            cluster=ID_R,
                            data = dat_outgroup,
                            model.name = "og target sim dummy 1")
summary(model_og_targetsim_d1)
set.seed(1234)
model_og_targetsim_d2 <- meta3(y=Fisher, v=Variance_F,
                              x=cbind(No, Bit),
                              cluster=ID_R,
                              data = dat_outgroup,
                              model.name = "og target sim dummy 1")
summary(model_og_targetsim_d2)
#adjust p-values for multiple comparisons
p <- cbind(summary(model_og_targetsim_d1)$coefficients$`Pr(>|z|)`[2], #lot vs. no
           summary(model_og_targetsim_d1)$coefficients$`Pr(>|z|)`[3], #bit vs. no
           summary(model_og_targetsim_d2)$coefficients$`Pr(>|z|)`[3]) #lot vs. bit
p.adjust(p, "BH")


##### 1.3. Methodological Moderators ##### 

### Research Design

#check data
summary(as.factor(dat_outgroup$TypeStudy))
#recode data
dat_outgroup$Design[dat_outgroup$TypeStudy==1] <- "Experiments"
dat_outgroup$Design[dat_outgroup$TypeStudy==2] <- "Quasi-Long"
dat_outgroup$Design[dat_outgroup$TypeStudy==3] <- "Correlations"
dat_outgroup$Design[dat_outgroup$TypeStudy==4] <- "Quasi-Long"
dat_outgroup$Design <- as.factor(dat_outgroup$Design)
table(dat_outgroup$Design)
#convert characters into dummy variables
dat_outgroup$Exp <- ifelse(dat_outgroup$TypeStudy==1, yes=1, no=0)
dat_outgroup$QuasiLong <- ifelse(dat_outgroup$TypeStudy==2 | 
                                   dat_outgroup$TypeStudy==4, yes=1, no=0)
dat_outgroup$Corr <- ifelse(dat_outgroup$TypeStudy==3, yes=1, no=0)

#fit model with effect coding
model_og_design <- meta3(y=Fisher, v=Variance_F,
                         intercept.constraints = 0,
                         x=cbind(Exp, QuasiLong, Corr),
                         cluster=ID_R,
                         data = dat_outgroup,
                         intervals.type = "LB",
                         model.name = "og design")
summary(model_og_design) #(mainly reported in SI Table C.1)
FisherZInv(coef(model_og_design)[3]) #Correlation (reported in text)
FisherZInv(coef(model_og_design)[1]) #Experiment (reported in text)
FisherZInv(coef(model_og_design)[2]) #Other designs (reported in text)
#explained variance (only reported in SI Table C.1)
summary(model_og_design)$R2.values[3] + summary(model_og_design)$R2.values[6]
#model fit (only reported in SI Table C.1)
anova(model_og_design, model_og_LBCI)

#fit model with effect coding
set.seed(1234)
model_og_design_d1 <- meta3(y=Fisher, v=Variance_F,
                         x=cbind(QuasiLong, Corr),
                         cluster=ID_R,
                         data = dat_outgroup,
                         model.name = "og design d1")
summary(model_og_design_d1)
set.seed(1234)
model_og_design_d2 <- meta3(y=Fisher, v=Variance_F,
                            x=cbind(Exp, Corr),
                            cluster=ID_R,
                            data = dat_outgroup,
                            model.name = "og design d2")
summary(model_og_design_d2)
#adjust p-values for multiple comparisons
p <- cbind(summary(model_og_design_d1)$coefficients$`Pr(>|z|)`[2], #exp. vs. quasi-long
           summary(model_og_design_d1)$coefficients$`Pr(>|z|)`[3], #exp. vs. corr. 
           summary(model_og_design_d2)$coefficients$`Pr(>|z|)`[3]) #quasi-long vs. corr.
p.adjust(p, "BH")


### Time Interval 

#check data
summary(dat_outgroup$DV_Time)
#convert characters into dummy variables
dat_outgroup$Direct <- ifelse(dat_outgroup$DV_Time==0, yes=1, no=0)
dat_outgroup$Interval <- ifelse(dat_outgroup$Direct==0, yes=1, no=0)
table(as.factor(dat_outgroup$Direct))
table(as.factor(dat_outgroup$Interval))

#fit model with effect coding
set.seed(1234)
model_og_interval <- meta3(y=Fisher, v=Variance_F,
                           intercept.constraints = 0, #if not, unidentified
                           x=cbind(Interval, Direct),
                           cluster=ID_R,
                           data = dat_outgroup,
                           intervals.type = "LB",
                           model.name = "og interval")
summary(model_og_interval) #(mainly reported in SI Table C.1)
FisherZInv(coef(model_og_interval)[1]) #Delay (reported in text)
FisherZInv(coef(model_og_interval)[2]) #No delay (reported in text)
#explained variance (only reported in SI Table C.1)
summary(model_og_interval)$R2.values[3] + summary(model_og_interval)$R2.values[6]
#model fit (only reported in SI Table C.1)
anova(model_og_interval, model_og_LBCI)

#fit model with dummy coding for studies with interval (=1)
set.seed(1234)
model_og_interval_d <- meta3(y=Fisher, v=Variance_F,
                           x=Interval,
                           cluster=ID_R,
                           data = dat_outgroup,
                           #intervals.type = "LB",
                           model.name = "og interval d")
summary(model_og_interval_d)


### Sampling Protocol

#recode data
dat_outgroup$Sample[dat_outgroup$GeneralPop==1] <- "General"
dat_outgroup$Sample[dat_outgroup$StudentPop==1] <- "Student"
dat_outgroup$Sample[is.na(dat_outgroup$Sample)] <- "Convenience"
dat_outgroup$Sample <- as.factor(dat_outgroup$Sample)
table(dat_outgroup$Sample)
#convert characters into dummy variables
dat_outgroup$General <- ifelse(dat_outgroup$Sample=="General", yes=1, no=0)
dat_outgroup$Student <- ifelse(dat_outgroup$Sample=="Student", yes=1, no=0)
dat_outgroup$Convenience <- ifelse(dat_outgroup$Sample=="Convenience", yes=1, no=0)

#fit model with effect coding
set.seed(1234)
model_og_sample <- meta3(y=Fisher, v=Variance_F,
                         intercept.constraints = 0,
                         x=cbind(General, Student, Convenience),
                         cluster=ID_R,
                         data = dat_outgroup,
                         intervals.type = "LB",
                         model.name = "og sample")
summary(model_og_sample) #(mainly reported in SI Table C.1)
FisherZInv(coef(model_og_sample)[1]) #General pop (reported in text)
FisherZInv(coef(model_og_sample)[2]) #Students (reported in text)
FisherZInv(coef(model_og_sample)[3]) #Convenience s. (reported in text)
#explained variance (only reported in SI Table C.1)
summary(model_og_sample)$R2.values[3] + summary(model_og_sample)$R2.values[6]
#model fit (only reported in SI Table C.1)
anova(model_og_sample, model_og_LBCI)

#fit model with dummy coding
set.seed(1234)
model_og_sample_d1 <- meta3(y=Fisher, v=Variance_F,
                         x=cbind(General, Convenience),
                         cluster=ID_R,
                         data = dat_outgroup,
                         model.name = "og sample d1")
summary(model_og_sample_d1)
set.seed(1234)
model_og_sample_d2 <- meta3(y=Fisher, v=Variance_F,
                            x=cbind(Student, Convenience),
                            cluster=ID_R,
                            data = dat_outgroup,
                            model.name = "og sample d2")
summary(model_og_sample_d2)
#adjust p-values for multiple comparisons
p <- cbind(summary(model_og_sample_d1)$coefficients$`Pr(>|z|)`[2], #general vs. student
           summary(model_og_sample_d1)$coefficients$`Pr(>|z|)`[3], #convenience vs. student
           summary(model_og_sample_d2)$coefficients$`Pr(>|z|)`[3]) #convenience vs. general
p.adjust(p, "BH")


### Location

#check data
summary(as.factor(dat_outgroup$Country))
#recode data
dat_outgroup$Country2[dat_outgroup$Country=="US"] <- "US"
dat_outgroup$Country2[dat_outgroup$Country=="Israel"] <- "Israel"
dat_outgroup$Country2[is.na(dat_outgroup$Country2)] <- "Other"
dat_outgroup$Country2 <- as.factor(dat_outgroup$Country2)
table(dat_outgroup$Country2)
#convert characters into dummy variables
dat_outgroup$US <- ifelse(dat_outgroup$Country=="US", yes=1, no=0)
dat_outgroup$Israel <- ifelse(dat_outgroup$Country=="Israel", yes=1, no=0)
dat_outgroup$Other_c <- ifelse(dat_outgroup$US==0 
                               & dat_outgroup$Israel==0, yes=1, no=0)

#fit model with effect coding
set.seed(1234)
model_og_geo <- meta3(y=Fisher, v=Variance_F,
                      intercept.constraints = 0,
                      x=cbind(US, Israel, Other_c),
                      cluster=ID_R,
                      data = dat_outgroup,
                      intervals.type = "LB",
                      model.name = "og geo")
summary(model_og_geo) #(mainly reported in SI Table C.1)
FisherZInv(coef(model_og_geo)[1]) #US (reported in text)
FisherZInv(coef(model_og_geo)[2]) #Israel (reported in text)
FisherZInv(coef(model_og_geo)[3]) #other countries
#explained variance (only reported in SI Table C.1)
summary(model_og_geo)$R2.values[3] + summary(model_og_geo)$R2.values[6]
#model fit (only reported in SI Table C.1)
anova(model_og_geo, model_og_LBCI)

#fit model with dummy coding
set.seed(1234)
model_og_geo_d1 <- meta3(y=Fisher, v=Variance_F,
                      x=cbind(US, Israel),
                      cluster=ID_R,
                      data = dat_outgroup,
                      model.name = "og geo d1")
summary(model_og_geo_d1)
set.seed(1234)
model_og_geo_d2 <- meta3(y=Fisher, v=Variance_F,
                         x=cbind(Israel, Other_c),
                         cluster=ID_R,
                         data = dat_outgroup,
                         model.name = "og geo d2")
summary(model_og_geo_d2)
#adjust p-values for multiple comparisons
p <- cbind(summary(model_og_geo_d1)$coefficients$`Pr(>|z|)`[2], #US vs. Other
           summary(model_og_geo_d1)$coefficients$`Pr(>|z|)`[3], #Israel vs. Other
           summary(model_og_geo_d2)$coefficients$`Pr(>|z|)`[2]) #Israel vs. US
p.adjust(p, "BH")


### Mean Age
#check data
summary(dat_outgroup$AgeMean) 
#fit model with age centered, but not standardized
set.seed(1234)
model_og_age <- meta3(y=Fisher, v=Variance_F,
                      x=scale(AgeMean, center = T, scale = F),
                      cluster=ID_R,
                      data = dat_outgroup, 
                      intervals.type = "LB",
                      model.name = "og age")
summary(model_og_age)
#explained variance (only reported in SI Table C.1)
summary(model_og_age)$R2.values[3] + summary(model_og_age)$R2.values[6]


### Percentage Women
#check data
summary(dat_outgroup$FemaleP) 
#fit model and check model fit improvement
set.seed(1234)
model_og_fem <- meta3(y=Fisher, v=Variance_F,
                      x=scale(FemaleP, center = T, scale = F),
                      cluster=ID_R,
                      data = dat_outgroup, 
                      intervals.type = "LB",
                      model.name = "og female")
summary(model_og_fem)#not significant 
#explained variance (only reported in SI Table C.1)
summary(model_og_fem)$R2.values[3] + summary(model_og_fem)$R2.values[6]


### Year of Data Collection
#check data
summary(dat_outgroup$StudyYear)
#fit model and check model fit improvement
set.seed(1234)
model_og_datayear <- meta3(y=Fisher, v=Variance_F,
                           x=scale(StudyYear, center = T, scale = F),
                           cluster=ID_R,
                           data = dat_outgroup, 
                           intervals.type = "LB",
                           model.name = "og data year")
summary(model_og_datayear)#not significant 
#explained variance (only reported in SI Table C.1)
summary(model_og_datayear)$R2.values[3] + summary(model_og_datayear)$R2.values[6]


### Year of Publication
#check data
summary(dat_outgroup$Year)
#fit model and check model fit improvement
set.seed(1234)
model_og_year <- meta3(y=Fisher, v=Variance_F,
                       x=scale(Year, center = T, scale = F),
                       cluster=ID_R,
                       data = dat_outgroup, 
                       intervals.type = "LB",
                       model.name = "og year")
summary(model_og_year)#not significant 
anova(model_og_year, model_og_LBCI) #(only reported in SI Table C.1)
#explained variance (only reported in SI Table C.1)
summary(model_og_geo)$R2.values[3] + summary(model_og_geo)$R2.values[6]


##### 2. Conservative Shift Hypothesis ##### 

##### 2.1. Terrorism-Related Moderators ##### 

### Type/Ideology of Terrorism

#check data
summary(as.factor(dat_conservatism$Type))
summary(as.factor(dat_conservatism$TerrorType))
#recode data
dat_conservatism$TerrorType2[dat_conservatism$Type==0] <- "No ideology"
dat_conservatism$TerrorType2[dat_conservatism$Type==1] <- "Islamist"
dat_conservatism$TerrorType2[is.na(dat_conservatism$TerrorType2)] <- "Other"
dat_conservatism$TerrorType2 <- as.factor(dat_conservatism$TerrorType2)
table(dat_conservatism$TerrorType2)
#convert into dummy variables
dat_conservatism$Islam <- ifelse(dat_conservatism$TerrorType2=="Islamist", yes=1, no=0)
dat_conservatism$No_Ideology <- ifelse(dat_conservatism$TerrorType2=="No ideology", yes=1, no=0)
dat_conservatism$Other_Ideology <- ifelse(dat_conservatism$TerrorType2=="Other", yes=1, no=0)

#fit model with effect coding
set.seed(1234)
model_cs_type <- meta3(y=Fisher, v=Variance_F,
                       intercept.constraints = 0,
                       x=cbind(Islam, Other_Ideology, No_Ideology),
                       cluster=ID_R, 
                       data = dat_conservatism,
                       intervals.type = "LB",
                       model.name = "cs ideology effect")
summary(model_cs_type)#(mainly reported in SI Table C.2)
FisherZInv(coef(model_cs_type)[1]) #Islam (reported in text)
FisherZInv(coef(model_cs_type)[2]) #No Islam (reported in text)
FisherZInv(coef(model_cs_type)[3]) #Non-spec. (reported in text)
#explained variance (only reported in SI Table C.2)
summary(model_cs_type)$R2.values[3] + summary(model_cs_type)$R2.values[6]
#model fit (only reported SI Table C.2)
anova(model_cs_type, model_cs_LBCI)

#fit models with dummy coding
set.seed(1234)
model_cs_type_d1 <- meta3(y=Fisher, v=Variance_F,
                          x=cbind(Islam, Other_Ideology),
                          cluster=ID_R,
                          data = dat_conservatism,
                          model.name = "cs target sim dummy 1")
summary(model_cs_type_d1)
set.seed(1234)
model_cs_type_d2 <- meta3(y=Fisher, v=Variance_F,
                          x=cbind(Other_Ideology, No_Ideology),
                          cluster=ID_R,
                          data = dat_conservatism,
                          model.name = "cs target sim dummy 1")
summary(model_cs_type_d2)
#adjust p-values for multiple comparisons
p <- cbind(summary(model_cs_type_d1)$coefficients$`Pr(>|z|)`[2], #Islam vs. no ideology
           summary(model_cs_type_d1)$coefficients$`Pr(>|z|)`[3], #Other ideology vs. no ideology
           summary(model_cs_type_d2)$coefficients$`Pr(>|z|)`[2]) #Other ideology vs. Islam
p.adjust(p, "BH")


### Terrorism Measure/Manipulation

#check data
summary(as.factor(dat_conservatism$TerrorMeasure))
#convert characters into dummy variables
dat_conservatism$Exposure <- ifelse(dat_conservatism$TerrorMeasure=="Acts of violence", yes=1, no=0)
dat_conservatism$R_Exposure <- ifelse(dat_conservatism$TerrorMeasure=="Reported exposure", yes=1, no=0)
dat_conservatism$Threat <- ifelse(dat_conservatism$TerrorMeasure=="Threat", yes=1, no=0)
dat_conservatism$Emotions <- ifelse(dat_conservatism$TerrorMeasure=="Fear" | 
                                      dat_conservatism$TerrorMeasure=="Anger", yes=1, no=0)
dat_conservatism$Fear <- ifelse(dat_conservatism$TerrorMeasure=="Fear", yes=1, no=0)
dat_conservatism$Anger <- ifelse(dat_conservatism$TerrorMeasure=="Anger", yes=1, no=0)
dat_conservatism$Other_IV <- ifelse(dat_conservatism$TerrorMeasure=="Other", yes=1, no=0)

#fit model with effect coding
set.seed(1234)
model_cs_ivm <- meta3(y=Fisher, v=Variance_F,
                      intercept.constraints = 0,
                      x=cbind(Exposure, R_Exposure, Threat, Emotions, Other_IV),
                      cluster=ID_R, 
                      data = dat_conservatism,
                      intervals.type = "LB",
                      model.name = "cs iv measure")
summary(model_cs_ivm) #(mainly reported in SI Table C.2)
FisherZInv(coef(model_cs_ivm)[1])#Objective (reported in text)
FisherZInv(coef(model_cs_ivm)[3])#Cognition (reported in text)
FisherZInv(coef(model_cs_ivm)[4])#Emotions (reported in text)
#explained variance (only reported in SI Table C.2)
summary(model_cs_ivm)$R2.values[3] + summary(model_cs_ivm)$R2.values[6]
#model fit (only reported SI Table C.2)
anova(model_cs_ivm, model_cs_LBCI)
#test anger and fear as well (as pre-registered)
set.seed(1234)
model_cs_ivm2 <- meta3(y=Fisher, v=Variance_F,
                       intercept.constraints = 0,
                       x=cbind(Exposure, R_Exposure, Threat, Fear, Anger, Other_IV),
                       cluster=ID_R, 
                       data = dat_conservatism,
                       intervals.type = "LB",
                       model.name = "cs iv measure2")
summary(model_cs_ivm2)

#fit models with dummy coding
set.seed(1234)
model_cs_ivm_d1 <- meta3(y=Fisher, v=Variance_F,
                         x=cbind(R_Exposure, Threat, Emotions, Other_IV),
                         cluster=ID_R, 
                         data = dat_conservatism,
                         model.name = "cs iv measure d1")
summary(model_cs_ivm_d1)
set.seed(1234)
model_cs_ivm_d2 <- meta3(y=Fisher, v=Variance_F,
                         x=cbind(Exposure, Threat, Emotions, Other_IV),
                         cluster=ID_R, 
                         data = dat_conservatism,
                         model.name = "cs iv measure d2")
summary(model_cs_ivm_d2)
set.seed(1234)
model_cs_ivm_d3 <- meta3(y=Fisher, v=Variance_F,
                         x=cbind(Exposure, R_Exposure, Emotions, Other_IV),
                         cluster=ID_R, 
                         data = dat_conservatism,
                         model.name = "cs iv measure d3")
summary(model_cs_ivm_d3)
set.seed(1234)
model_cs_ivm_d4 <- meta3(y=Fisher, v=Variance_F,
                         x=cbind(Exposure, R_Exposure, Threat, Other_IV),
                         cluster=ID_R, 
                         data = dat_conservatism,
                         model.name = "cs iv measure d4")
summary(model_cs_ivm_d4)
#adjust p-values for multiple comparisons
p <- cbind(summary(model_cs_ivm_d1)$coefficients$`Pr(>|z|)`[2], #obj. exp. vs. subj. exp.:no
           summary(model_cs_ivm_d1)$coefficients$`Pr(>|z|)`[3], #obj. exp. vs. threat:yes
           summary(model_cs_ivm_d1)$coefficients$`Pr(>|z|)`[4], #obj. exp. vs. emotions:yes
           summary(model_cs_ivm_d1)$coefficients$`Pr(>|z|)`[5], #obj. exp. vs. other:no
           summary(model_cs_ivm_d2)$coefficients$`Pr(>|z|)`[3], #subj. exp. vs. threat:yes
           summary(model_cs_ivm_d2)$coefficients$`Pr(>|z|)`[4], #subj. exp. vs. emotions:yes
           summary(model_cs_ivm_d2)$coefficients$`Pr(>|z|)`[5], #subj. exp. vs. other:no
           summary(model_cs_ivm_d3)$coefficients$`Pr(>|z|)`[4], #threat vs. emotions:no
           summary(model_cs_ivm_d3)$coefficients$`Pr(>|z|)`[5], #threat vs. other:yes
           summary(model_cs_ivm_d4)$coefficients$`Pr(>|z|)`[5]) #emotions vs. other:yes
p.adjust(p, "BH")


###### 2.2. Outcome-Related Moderators ##### 

### Measurement

#check data
summary(as.factor(dat_conservatism$PA_Category))
#recode data
dat_conservatism$DV[dat_conservatism$PA_Category==5] <- "RWA"
dat_conservatism$DV[dat_conservatism$PA_Category==6] <- "SDO"
dat_conservatism$DV[dat_conservatism$PA_Category==7] <- "Ideology"
dat_conservatism$DV[dat_conservatism$PA_Category==8] <- "Military"
dat_conservatism$DV[dat_conservatism$PA_Category==9] <- "Civil"
dat_conservatism$DV[dat_conservatism$PA_Category==11] <- "Other"
dat_conservatism$DV <- as.factor(dat_conservatism$DV)
table(dat_conservatism$DV)
#convert into dummy variables
dat_conservatism$RWA <- ifelse(dat_conservatism$DV=="RWA", yes=1, no=0)
dat_conservatism$SDO <- ifelse(dat_conservatism$DV=="SDO", yes=1, no=0)
dat_conservatism$Ideology <- ifelse(dat_conservatism$DV=="Ideology", yes=1, no=0)
dat_conservatism$Military <- ifelse(dat_conservatism$DV=="Military", yes=1, no=0)
dat_conservatism$Liberties <- ifelse(dat_conservatism$DV=="Civil", yes=1, no=0)
dat_conservatism$Other <- ifelse(dat_conservatism$DV=="Other", yes=1, no=0)

#fit model with effect coding
set.seed(1234)
model_cs_outcome <- meta3(y=Fisher, v=Variance_F,
                         intercept.constraints = 0,
                         x=cbind(RWA, SDO, Ideology, Military, Liberties, Other),
                         cluster=ID_R,
                         data = dat_conservatism, 
                         intervals.type = "LB",
                         model.name = "cs target")
summary(model_cs_outcome) #(mainly reported in SI Table C.2)
FisherZInv(coef(model_cs_outcome)[1]) #RWA (reported in text)
FisherZInv(coef(model_cs_outcome)[5]) #CL (reported in text)
FisherZInv(coef(model_cs_outcome)[2]) #SDO (reported in text)
FisherZInv(coef(model_cs_outcome)[6]) #Other (reported in text)
#explained variance (only reported in SI Table C.2)
summary(model_cs_outcome)$R2.values[3] + summary(model_cs_outcome)$R2.values[6]
#model fit (only reported SI Table C.2)
anova(model_cs_outcome, model_cs_LBCI)

#fit model with dummy coding
set.seed(1234)
model_cs_outcome_d1 <- meta3(y=Fisher, v=Variance_F,
                            x=cbind(SDO, Ideology, Military, Liberties, Other),
                            cluster=ID_R,
                            data = dat_conservatism,
                            model.name = "cs target d1")
summary(model_cs_outcome_d1)
set.seed(1234)
model_cs_outcome_d2 <- meta3(y=Fisher, v=Variance_F,
                             x=cbind(RWA, Ideology, Military, Liberties, Other),
                             cluster=ID_R,
                             data = dat_conservatism,
                             model.name = "cs target d1")
summary(model_cs_outcome_d2)
set.seed(1234)
model_cs_outcome_d3 <- meta3(y=Fisher, v=Variance_F,
                             x=cbind(RWA, SDO, Military, Liberties, Other),
                             cluster=ID_R,
                             data = dat_conservatism,
                             model.name = "cs target d1")
summary(model_cs_outcome_d3)
set.seed(1234)
model_cs_outcome_d4 <- meta3(y=Fisher, v=Variance_F,
                             x=cbind(RWA, SDO, Ideology, Liberties, Other),
                             cluster=ID_R,
                             data = dat_conservatism,
                             model.name = "cs target d1")
summary(model_cs_outcome_d4)
set.seed(1234)
model_cs_outcome_d5 <- meta3(y=Fisher, v=Variance_F,
                             x=cbind(RWA, SDO, Ideology, Military, Other),
                             cluster=ID_R,
                             data = dat_conservatism,
                             model.name = "cs target d1")
summary(model_cs_outcome_d5)
#adjust p-values for multiple comparisons
p <- cbind(summary(model_cs_outcome_d1)$coefficients$`Pr(>|z|)`[2],#yes
           summary(model_cs_outcome_d1)$coefficients$`Pr(>|z|)`[3],#yes 
           summary(model_cs_outcome_d1)$coefficients$`Pr(>|z|)`[4],#no
           summary(model_cs_outcome_d1)$coefficients$`Pr(>|z|)`[5],#no
           summary(model_cs_outcome_d1)$coefficients$`Pr(>|z|)`[6],#yes
           summary(model_cs_outcome_d2)$coefficients$`Pr(>|z|)`[3],#no 
           summary(model_cs_outcome_d2)$coefficients$`Pr(>|z|)`[4],#no
           summary(model_cs_outcome_d2)$coefficients$`Pr(>|z|)`[5],#yes
           summary(model_cs_outcome_d2)$coefficients$`Pr(>|z|)`[6],#yes
           summary(model_cs_outcome_d3)$coefficients$`Pr(>|z|)`[4],#no
           summary(model_cs_outcome_d3)$coefficients$`Pr(>|z|)`[5],#yes
           summary(model_cs_outcome_d3)$coefficients$`Pr(>|z|)`[6],#yes
           summary(model_cs_outcome_d4)$coefficients$`Pr(>|z|)`[5],#no
           summary(model_cs_outcome_d4)$coefficients$`Pr(>|z|)`[6],#yes
           summary(model_cs_outcome_d5)$coefficients$`Pr(>|z|)`[6])#yes
p.adjust(p, "BH")


##### 2.3. Methodological Moderators ##### 

### Study Design

#check data
summary(as.factor(dat_conservatism$TypeStudy))
#recode data
dat_conservatism$Design[dat_conservatism$TypeStudy==1] <- "Experiments"
dat_conservatism$Design[dat_conservatism$TypeStudy==2] <- "Quasi-Long"
dat_conservatism$Design[dat_conservatism$TypeStudy==3] <- "Correlations"
dat_conservatism$Design[dat_conservatism$TypeStudy==4] <- "Quasi-Long"
dat_conservatism$Design <- as.factor(dat_conservatism$Design)
table(dat_conservatism$Design)
#convert characters into dummy variables
dat_conservatism$Exp <- ifelse(dat_conservatism$TypeStudy==1, yes=1, no=0)
dat_conservatism$QuasiLong <- ifelse(dat_conservatism$TypeStudy==2 | 
                                       dat_conservatism$TypeStudy==4, yes=1, no=0)
dat_conservatism$Corr <- ifelse(dat_conservatism$TypeStudy==3, yes=1, no=0)

#fit model with effect coding
model_cs_design <- meta3(y=Fisher, v=Variance_F,
                         intercept.constraints = 0,
                         x=cbind(Exp, QuasiLong, Corr),
                         cluster=ID_R,
                         data = dat_conservatism,
                         intervals.type = "LB",
                         model.name = "cs design")
summary(model_cs_design) #(mainly reported in SI Table C.2)
FisherZInv(coef(model_cs_design)[3]) #Correlation (reported in text)
FisherZInv(coef(model_cs_design)[1]) #Experiment (reported in text)
FisherZInv(coef(model_cs_design)[2]) #Other designs (reported in text)
#explained variance (only reported in SI Table C.2)
summary(model_cs_design)$R2.values[3] + summary(model_cs_design)$R2.values[6]
#model fit (only reported SI Table C.2)
anova(model_cs_design, model_cs_LBCI)

#fit model with effect coding
set.seed(1234)
model_cs_design_d1 <- meta3(y=Fisher, v=Variance_F,
                            x=cbind(QuasiLong, Corr),
                            cluster=ID_R,
                            data = dat_outgroup,
                            model.name = "cs design d1")
summary(model_cs_design_d1)
set.seed(1234)
model_cs_design_d2 <- meta3(y=Fisher, v=Variance_F,
                            x=cbind(Exp, Corr),
                            cluster=ID_R,
                            data = dat_outgroup,
                            model.name = "cs design d2")
summary(model_cs_design_d2)
#adjust p-values for multiple comparisons
p <- cbind(summary(model_cs_design_d1)$coefficients$`Pr(>|z|)`[2], #exp. vs. quasi-long: no
           summary(model_cs_design_d1)$coefficients$`Pr(>|z|)`[3], #exp. vs. corr.: yes
           summary(model_cs_design_d2)$coefficients$`Pr(>|z|)`[3]) #quasi-long vs. corr.; yes
p.adjust(p, "BH")


### Time Interval 

#check data
summary(dat_conservatism$DV_Time)
#convert characters into dummy variables
dat_conservatism$Direct <- ifelse(dat_conservatism$DV_Time==0, yes=1, no=0)
dat_conservatism$Interval <- ifelse(dat_conservatism$Direct==0, yes=1, no=0)
table(as.factor(dat_conservatism$Direct))
table(as.factor(dat_conservatism$Interval))

#fit model with effect coding
set.seed(1234)
model_cs_interval <- meta3(y=Fisher, v=Variance_F,
                           intercept.constraints = 0,
                           x=cbind(Interval, Direct),
                           cluster=ID_R,
                           data = dat_conservatism,
                           intervals.type = "LB",
                           model.name = "cs interval")
summary(model_cs_interval) #(mainly reported in SI Table C.2)
FisherZInv(coef(model_cs_interval)[1]) #Delay (reported in footnote 15)
FisherZInv(coef(model_cs_interval)[2]) #No delay (reported in footnote 15)
#explained variance (only reported in SI Table C.2)
summary(model_cs_interval)$R2.values[3] + summary(model_cs_interval)$R2.values[6]
#model fit (only reported SI Table C.2)
anova(model_cs_interval, model_cs_LBCI)

#fit model with dummy coding for studies with interval (=1)
set.seed(1234)
model_cs_interval_d <- meta3(y=Fisher, v=Variance_F,
                             x=Interval,
                             cluster=ID_R,
                             data = dat_conservatism,
                             model.name = "cs interval")
summary(model_cs_interval_d)


### Sampling Protocol

#recode data
dat_conservatism$Sample[dat_conservatism$GeneralPop==1] <- "General"
dat_conservatism$Sample[dat_conservatism$StudentPop==1] <- "Student"
dat_conservatism$Sample[is.na(dat_conservatism$Sample)] <- "Convenience"
dat_conservatism$Sample <- as.factor(dat_conservatism$Sample)
table(dat_conservatism$Sample)
#convert characters into dummy variables
dat_conservatism$General <- ifelse(dat_conservatism$Sample=="General", yes=1, no=0)
dat_conservatism$Student <- ifelse(dat_conservatism$Sample=="Student", yes=1, no=0)
dat_conservatism$Convenience <- ifelse(dat_conservatism$Sample=="Convenience", yes=1, no=0)

#fit model with effect coding
set.seed(1234)
model_cs_sample <- meta3(y=Fisher, v=Variance_F,
                         intercept.constraints = 0,
                         x=cbind(General, Student, Convenience),
                         cluster=ID_R,
                         data = dat_conservatism,
                         intervals.type = "LB",
                         model.name = "cs sample")
summary(model_cs_sample) #(mainly reported in SI Table C.2)
FisherZInv(coef(model_cs_sample)[1]) #General pop (reported in text)
FisherZInv(coef(model_cs_sample)[2]) #Students (reported in text)
FisherZInv(coef(model_cs_sample)[3]) #Convenience s. (reported in text)
#explained variance (only reported in SI Table C.2)
summary(model_cs_sample)$R2.values[3] + summary(model_cs_sample)$R2.values[6]
#model fit (only reported SI Table C.2)
anova(model_cs_sample, model_cs_LBCI)

#fit model with dummy coding
set.seed(1234)
model_cs_sample_d1 <- meta3(y=Fisher, v=Variance_F,
                            x=cbind(General, Convenience),
                            cluster=ID_R,
                            data = dat_conservatism,
                            model.name = "cs sample d1")
summary(model_cs_sample_d1)
set.seed(1234)
model_cs_sample_d2 <- meta3(y=Fisher, v=Variance_F,
                            x=cbind(Student, Convenience),
                            cluster=ID_R,
                            data = dat_conservatism,
                            model.name = "cs sample d2")
summary(model_cs_sample_d2)
#adjust p-values for multiple comparisons
p <- cbind(summary(model_cs_sample_d1)$coefficients$`Pr(>|z|)`[2], #general vs. students: no
           summary(model_cs_sample_d1)$coefficients$`Pr(>|z|)`[3], #convenience vs. students: loses significance after correction
           summary(model_cs_sample_d2)$coefficients$`Pr(>|z|)`[3]) #convenience vs. general: yes
p.adjust(p, "BH")


### Location

#check data
summary(as.factor(dat_conservatism$Country))
#recode data
dat_conservatism$Country2[dat_conservatism$Country=="US"] <- "US"
dat_conservatism$Country2[dat_conservatism$Country=="Israel"] <- "Israel"
dat_conservatism$Country2[is.na(dat_conservatism$Country2)] <- "Other"
dat_conservatism$Country2 <- as.factor(dat_conservatism$Country2)
table(dat_conservatism$Country2)
#convert characters into dummy variables
dat_conservatism$US <- ifelse(dat_conservatism$Country=="US", yes=1, no=0)
dat_conservatism$Israel <- ifelse(dat_conservatism$Country=="Israel", yes=1, no=0)
dat_conservatism$Other_c <- ifelse(dat_conservatism$US==0 
                               & dat_conservatism$Israel==0, yes=1, no=0)

#fit model with effect coding
set.seed(1234)
model_cs_geo <- meta3(y=Fisher, v=Variance_F,
                      intercept.constraints = 0,
                      x=cbind(US, Israel, Other_c),
                      cluster=ID_R,
                      data = dat_conservatism,
                      intervals.type = "LB",
                      model.name = "cs geo")
summary(model_cs_geo) #(mainly reported in SI Table C.2)
FisherZInv(coef(model_cs_geo)[1]) #US (reported in text)
FisherZInv(coef(model_cs_geo)[2]) #Israel (reported in text)
FisherZInv(coef(model_cs_geo)[3]) #other countries (reported in text)
#explained variance (only reported in SI Table C.2)
summary(model_cs_geo)$R2.values[3] + summary(model_cs_geo)$R2.values[6]
#model fit (only reported SI Table C.2)
anova(model_cs_geo, model_cs_LBCI)

#fit model with dummy coding
set.seed(1234)
model_cs_geo_d1 <- meta3(y=Fisher, v=Variance_F,
                         x=cbind(US, Israel),
                         cluster=ID_R,
                         data = dat_conservatism,
                         model.name = "cs geo d1")
summary(model_cs_geo_d1)
set.seed(1234)
model_cs_geo_d2 <- meta3(y=Fisher, v=Variance_F,
                         x=cbind(Israel, Other_c),
                         cluster=ID_R,
                         data = dat_conservatism,
                         model.name = "cs geo d2")
summary(model_cs_geo_d2)
#adjust p-values for multiple comparisons
p <- cbind(summary(model_cs_geo_d1)$coefficients$`Pr(>|z|)`[2], 
           summary(model_cs_geo_d1)$coefficients$`Pr(>|z|)`[3], 
           summary(model_cs_geo_d2)$coefficients$`Pr(>|z|)`[2])
p.adjust(p, "BH")


### Mean Age
#check data
summary(dat_conservatism$AgeMean) 
#fit model with age centered, but not standardized
set.seed(1234)
model_cs_age <- meta3(y=Fisher, v=Variance_F,
                      x=scale(AgeMean, center = T, scale = F),
                      cluster=ID_R,
                      data = dat_conservatism, 
                      intervals.type = "LB",
                      model.name = "cs age")
summary(model_cs_age)#not significant
#explained variance (only reported in SI Table C.2)
summary(model_cs_age)$R2.values[3] + summary(model_cs_age)$R2.values[6]


### Percentage Women
#check data
summary(dat_conservatism$FemaleP) 
#fit model and check model fit improvement
set.seed(1234)
model_cs_fem <- meta3(y=Fisher, v=Variance_F,
                      x=scale(FemaleP, center = T, scale = F),
                      cluster=ID_R,
                      data = dat_conservatism, 
                      intervals.type = "LB",
                      model.name = "cs female")
summary(model_cs_fem)#not significant 
#summary(rerun(model_cs_fem))
#explained variance (only reported in SI Table C.2)
summary(model_cs_fem)$R2.values[3] + summary(model_cs_fem)$R2.values[6]


### Year of Data Collection
#check data
summary(dat_conservatism$StudyYear)
#fit model and check model fit improvement
set.seed(1234)
model_cs_datayear <- meta3(y=Fisher, v=Variance_F,
                           x=scale(StudyYear, center = T, scale = F),
                           cluster=ID_R,
                           data = dat_conservatism, 
                           intervals.type = "LB",
                           model.name = "cs data year")
summary(model_cs_datayear)#not significant 
#explained variance (only reported in SI Table C.2)
summary(model_cs_datayear)$R2.values[3] + summary(model_cs_datayear)$R2.values[6]


### Year of Publication
#check data
summary(dat_conservatism$Year)
#fit model and check model fit improvement
set.seed(1234)
model_cs_year <- meta3(y=Fisher, v=Variance_F,
                       x=scale(Year, center = T, scale = F),
                       cluster=ID_R,
                       data = dat_conservatism, 
                       intervals.type = "LB",
                       model.name = "cs year")
summary(model_cs_year)#not significant 
#explained variance (only reported in SI Table C.2)
summary(model_cs_year)$R2.values[3] + summary(model_cs_year)$R2.values[6]
#model fit (reported in SI Table C.2)
anova(model_cs_year, model_cs_LBCI)


##### 3. Rally-Around-The-Flag Hypothesis ##### 

##### 3.1. Terrorism-Related Moderators ##### 

### Type/Ideology of Terrorism

#check data
summary(as.factor(dat_rally$Type))
summary(as.factor(dat_rally$TerrorType))
#recode data
dat_rally$TerrorType2[dat_rally$Type==0] <- "No ideology"
dat_rally$TerrorType2[dat_rally$Type==1] <- "Islamist"
dat_rally$TerrorType2[is.na(dat_rally$TerrorType2)] <- "Other"
dat_rally$TerrorType2 <- as.factor(dat_rally$TerrorType2)
table(dat_rally$TerrorType2)
#convert into dummy variables
dat_rally$Islam <- ifelse(dat_rally$TerrorType2=="Islamist", yes=1, no=0)
dat_rally$No_Ideology <- ifelse(dat_rally$TerrorType2=="No ideology", yes=1, no=0)
dat_rally$Other_Ideology <- ifelse(dat_rally$TerrorType2=="Other", yes=1, no=0)

#fit model with effect coding
set.seed(1234)
model_rf_type <- meta3(y=Fisher, v=Variance_F,
                       intercept.constraints = 0,
                       x=cbind(Islam, Other_Ideology, No_Ideology),
                       cluster=ID_R, 
                       data = dat_rally,
                       intervals.type = "LB",
                       model.name = "rf ideology effect")
summary(model_rf_type)#(mainly reported in Table C.3)
FisherZInv(coef(model_rf_type)[1]) #Islam
FisherZInv(coef(model_rf_type)[2]) #Other ideology (reported in text)
FisherZInv(coef(model_rf_type)[3]) #Unspecified
#explained variance (only reported SI Table C.3)
summary(model_rf_type)$R2.values[3] + summary(model_rf_type)$R2.values[6]
#model fit (only reported SI Table C.3)
anova(model_rf_type, model_rf_LBCI)

#fit models with dummy coding
set.seed(1234)
model_rf_type_d1 <- meta3(y=Fisher, v=Variance_F,
                          x=cbind(Islam, Other_Ideology),
                          cluster=ID_R,
                          data = dat_rally,
                          model.name = "rf ideology dummy 1")
summary(model_rf_type_d1)
set.seed(1234)
model_rf_type_d2 <- meta3(y=Fisher, v=Variance_F,
                          x=cbind(Other_Ideology, No_Ideology),
                          cluster=ID_R,
                          data = dat_rally,
                          model.name = "rf ideology dummy 1")
summary(model_rf_type_d2)


### Terrorism Measure/Manipulation

#check data
summary(as.factor(dat_rally$TerrorMeasure))
#convert characters into dummy variables
dat_rally$Exposure <- ifelse(dat_rally$TerrorMeasure=="Acts of violence", yes=1, no=0)
dat_rally$R_Exposure <- ifelse(dat_rally$TerrorMeasure=="Reported exposure", yes=1, no=0)
dat_rally$Threat <- ifelse(dat_rally$TerrorMeasure=="Threat", yes=1, no=0)
dat_rally$Emotions <- ifelse(dat_rally$TerrorMeasure=="Fear" | 
                                      dat_rally$TerrorMeasure=="Anger", yes=1, no=0)
dat_rally$Fear <- ifelse(dat_rally$TerrorMeasure=="Fear", yes=1, no=0)
dat_rally$Anger <- ifelse(dat_rally$TerrorMeasure=="Anger", yes=1, no=0)
dat_rally$Other_IV <- ifelse(dat_rally$TerrorMeasure=="Other", yes=1, no=0)

#fit model with effect coding
set.seed(1234)
model_rf_ivm <- meta3(y=Fisher, v=Variance_F,
                      intercept.constraints = 0,
                      x=cbind(Exposure, R_Exposure, Threat, Emotions, Other_IV),
                      cluster=ID_R, 
                      data = dat_rally,
                      intervals.type = "LB",
                      model.name = "rf iv measure")
summary(model_rf_ivm) #(mainly reported in Table C.3)
#explained variance (only reported SI Table C.3)
summary(model_rf_ivm)$R2.values[3] + summary(model_rf_ivm)$R2.values[6]
#model fit (only reported SI Table C.3)
anova(model_rf_ivm, model_rf_LBCI)
#test anger and fear as well (as pre-registered; reported in SI)
set.seed(1234)
model_rf_ivm2 <- meta3(y=Fisher, v=Variance_F,
                       intercept.constraints = 0,
                       x=cbind(Exposure, R_Exposure, Threat, Fear, Anger, Other_IV),
                       cluster=ID_R,
                       data = dat_rally,
                       intervals.type = "LB",
                       model.name = "rf iv measure2")
summary(model_rf_ivm2)

#fit models with dummy coding
set.seed(1234)
model_rf_ivm_d1 <- meta3(y=Fisher, v=Variance_F,
                         x=cbind(R_Exposure, Threat, Emotions, Other_IV),
                         cluster=ID_R, 
                         data = dat_rally,
                         model.name = "rf iv measure d1")
summary(model_rf_ivm_d1)
set.seed(1234)
model_rf_ivm_d2 <- meta3(y=Fisher, v=Variance_F,
                         x=cbind(Exposure, Threat, Emotions, Other_IV),
                         cluster=ID_R, 
                         data = dat_rally,
                         model.name = "rf iv measure d2")
summary(model_rf_ivm_d2)
set.seed(1234)
model_rf_ivm_d3 <- meta3(y=Fisher, v=Variance_F,
                         x=cbind(Exposure, R_Exposure, Emotions, Other_IV),
                         cluster=ID_R, 
                         data = dat_rally,
                         model.name = "rf iv measure d3")
summary(model_rf_ivm_d3)
set.seed(1234)
model_rf_ivm_d4 <- meta3(y=Fisher, v=Variance_F,
                         x=cbind(Exposure, R_Exposure, Threat, Other_IV),
                         cluster=ID_R, 
                         data = dat_rally,
                         model.name = "rf iv measure d4")
summary(model_rf_ivm_d4)


### 9/11 Effect

#check data
dat_rally$ExactAttack[is.na(dat_rally$ExactAttack)] <- 0
summary(as.factor(dat_rally$ExactAttack))
#convert characters into dummy variables
dat_rally$Nine11 <- ifelse(dat_rally$ExactAttack=="9/11 Attacks" |
                             dat_rally$ExactAttack=="9/11 Attacks and 11/03 Madrid", yes=1, no=0)
dat_rally$NoNine11 <- ifelse(dat_rally$Nine11==0, yes=1, no=0)
table(as.factor(dat_rally$Nine11))

#fit model with 9/11 effect coding
set.seed(1234)
model_rf_nine <- meta3(y=Fisher, v=Variance_F,
                       intercept.constraints = 0, #if not, unidentified
                       x=cbind(NoNine11, Nine11),
                       cluster=ID_R,
                       data = dat_rally,
                       intervals.type = "LB",
                       model.name = "rf interval")
summary(model_rf_nine) #(mainly reported in Table C.3)
FisherZInv(coef(model_rf_nine)[2]) #Reference to 9/11 (reported in text)
FisherZInv(coef(model_rf_nine)[1]) #No reference to 9/11 (reported in text)
#explained variance (only reported SI Table C.3)
summary(model_rf_nine)$R2.values[3] + summary(model_rf_nine)$R2.values[6]
#model fit (only reported SI Table C.3)
anova(model_rf_nine, model_rf_LBCI)

#fit model with 9/11 dummy coding
set.seed(1234)
model_rf_nine_d <- meta3(y=Fisher, v=Variance_F,
                         #intercept.constraints = 0,
                         x=Nine11,
                         cluster=ID_R, 
                         data = dat_rally,
                         intervals.type = "LB",
                         model.name = "rf nine eleven d")
summary(model_rf_nine_d)


###### 3.2. Outcome-Related Moderators ##### 

### DV measurements used

#check the data
summary(as.factor(dat_rally$PA_Category))
#convert into dummy variables
dat_rally$Politicians <- ifelse(dat_rally$PA_Category==4, yes=1, no=0)
dat_rally$PolTrust <- ifelse(dat_rally$PA_Category==2 |
                                  dat_rally$PA_Category==3 , yes=1, no=0)
dat_rally$Nationalism <- ifelse(dat_rally$PA_Category==1, yes=1, no=0)
table(as.factor(dat_rally$Politicians))
table(as.factor(dat_rally$PolTrust))
table(as.factor(dat_rally$Nationalism))

#fit model with effect coding
set.seed(1234)
model_rf_dv <- meta3(y=Fisher, v=Variance_F,
                     intercept.constraints = 0,
                     x=cbind(Politicians, PolTrust, Nationalism),
                     cluster=ID_R,
                     data = dat_rally, 
                     intervals.type = "LB",
                     model.name = "rf outcome")
summary(model_rf_dv) #(mainly reported in SI Table C.3)
#explained variance (only reported in SI Table C.3)
summary(model_rf_dv)$R2.values[3] + summary(model_rf_dv)$R2.values[6]
#model fit (only reported in SI Table C.3)
anova(model_rf_dv, model_rf_LBCI)

#fit model with dummy coding
set.seed(1234)
model_rf_dv_d1 <- meta3(y=Fisher, v=Variance_F,
                        x=cbind(PolTrust, Nationalism),
                        cluster=ID_R,
                        data = dat_rally, 
                        model.name = "rf outcome d1")
summary(model_rf_dv_d1)
set.seed(1234)
model_rf_dv_d2 <- meta3(y=Fisher, v=Variance_F,
                        x=cbind(Politicians, Nationalism),
                        cluster=ID_R,
                        data = dat_rally, 
                        model.name = "rf outcome d1")
summary(model_rf_dv_d2)


### Republican Effect

#recode data
summary(as.factor(dat_rally$Republican))
dat_rally$Republican2 <- ifelse(dat_rally$Republican==1, yes=1, no=0)
dat_rally$Republican2[is.na(dat_rally$Republican)] <- 0
summary(as.factor(dat_rally$Republican2))

#fit model with dummy coding
set.seed(1234)
model_rf_rep_d <- meta3(y=Fisher, v=Variance_F,
                      x=Republican2,
                      cluster=ID_R,
                      data = dat_rally,
                      intervals.type = "LB",
                      model.name = "rf rep dummy")
summary(model_rf_rep_d)
anova(model_rf_rep_d, model_rf_LBCI)


### Incumbency Effect

#recode data
summary(as.factor(dat_rally$Incumbent))
dat_rally$Incumbent2 <- ifelse(dat_rally$Incumbent==1, yes=1, no=0)
dat_rally$Incumbent2[is.na(dat_rally$Incumbent)] <- 0
summary(as.factor(dat_rally$Incumbent2))

#fit model with dummy coding
set.seed(1234)
model_rf_inc_d <- meta3(y=Fisher, v=Variance_F,
                      x=Incumbent2,
                      cluster=ID_R,
                      data = dat_rally,
                      intervals.type = "LB",
                      model.name = "rf incumbent dummy")
summary(model_rf_inc_d) 
anova(model_rf_inc_d, model_rf_LBCI)


### Bush Effect

#recode data
dat_rally$Bush <- ifelse(dat_rally$Name=="George W. Bush", yes=1, no=0)
dat_rally$Bush[is.na(dat_rally$Bush)] <- 0
summary(as.factor(dat_rally$Bush))

#fit model with dummy coding
set.seed(1234)
model_rf_Bush <- meta3(y=Fisher, v=Variance_F,
                       x=Bush,
                       cluster=ID_R,
                       data = dat_rally,
                       intervals.type = "LB",
                       model.name = "rf Bush")
summary(model_rf_Bush) #(mainly reported in SI Table C.3)
#explained variance (only reported in SI Table C.3)
summary(model_rf_Bush)$R2.values[3] + summary(model_rf_Bush)$R2.values[6]
#model fit
anova(model_rf_Bush, model_rf_LBCI)


##### 3.3. Methodological Moderators ##### 

### Study Design

#check data
summary(as.factor(dat_rally$TypeStudy))
#recode data
dat_rally$Design[dat_rally$TypeStudy==1] <- "Experiments"
dat_rally$Design[dat_rally$TypeStudy==2] <- "Quasi-Long"
dat_rally$Design[dat_rally$TypeStudy==3] <- "Correlations"
dat_rally$Design[dat_rally$TypeStudy==4] <- "Quasi-Long"
dat_rally$Design <- as.factor(dat_rally$Design)
table(dat_rally$Design)
#convert characters into dummy variables
dat_rally$Exp <- ifelse(dat_rally$TypeStudy==1, yes=1, no=0)
dat_rally$QuasiLong <- ifelse(dat_rally$TypeStudy==2 | 
                                dat_rally$TypeStudy==4, yes=1, no=0)
dat_rally$Corr <- ifelse(dat_rally$TypeStudy==3, yes=1, no=0)

#fit model with effect coding
model_rf_design <- meta3(y=Fisher, v=Variance_F,
                         intercept.constraints = 0,
                         x=cbind(Exp, QuasiLong, Corr),
                         cluster=ID_R,
                         data = dat_rally,
                         intervals.type = "LB",
                         model.name = "rf design")
summary(model_rf_design) #(mainly reported in SI Table C.3)
FisherZInv(coef(model_rf_design)[1]) #(quasi)Experiment (reported in text)
FisherZInv(coef(model_rf_design)[2]) #Longitudinal (reported in text)
FisherZInv(coef(model_rf_design)[3]) #Correlation (reported in text)
#explained variance (only reported in SI Table C.3)
summary(model_rf_design)$R2.values[3] + summary(model_rf_design)$R2.values[6]
#model fit (only reported in SI Table C.3)
anova(model_rf_design, model_rf_LBCI)

#fit model with effect coding
set.seed(1234)
model_rf_design_d1 <- meta3(y=Fisher, v=Variance_F,
                            x=cbind(QuasiLong, Corr),
                            cluster=ID_R,
                            data = dat_rally,
                            model.name = "og design d1")
summary(model_rf_design_d1)
set.seed(1234)
model_rf_design_d2 <- meta3(y=Fisher, v=Variance_F,
                            x=cbind(Exp, Corr),
                            cluster=ID_R,
                            data = dat_rally,
                            model.name = "og design d2")
summary(model_rf_design_d2)


### Time Interval 

#check data
summary(dat_rally$DV_Time)
#convert characters into dummy variables
dat_rally$Direct <- ifelse(dat_rally$DV_Time==0, yes=1, no=0)
dat_rally$Interval <- ifelse(dat_rally$Direct==0, yes=1, no=0)
table(as.factor(dat_rally$Direct))
table(as.factor(dat_rally$Interval))

#fit model with effect coding
set.seed(1234)
model_rf_interval <- meta3(y=Fisher, v=Variance_F,
                           intercept.constraints = 0, #if not, unidentified
                           x=cbind(Interval, Direct),
                           cluster=ID_R,
                           data = dat_rally,
                           intervals.type = "LB",
                           model.name = "rf interval")
summary(model_rf_interval) #(mainly reported in SI Table C.3)
FisherZInv(coef(model_rf_interval)[1]) #Delay (reported in text)
FisherZInv(coef(model_rf_interval)[2]) #No delay (reported in text)
#explained variance (only reported in SI Table C.3)
summary(model_rf_interval)$R2.values[3] + summary(model_rf_interval)$R2.values[6]
#model fit (only reported in SI Table C.3)
anova(model_rf_interval, model_rf_LBCI)


### Sampling Protocol

#recode data
dat_rally$Sample[dat_rally$GeneralPop==1] <- "General"
dat_rally$Sample[dat_rally$StudentPop==1] <- "Student"
dat_rally$Sample[is.na(dat_rally$Sample)] <- "Convenience"
dat_rally$Sample <- as.factor(dat_rally$Sample)
table(dat_rally$Sample)
#convert characters into dummy variables
dat_rally$General <- ifelse(dat_rally$Sample=="General", yes=1, no=0)
dat_rally$Student <- ifelse(dat_rally$Sample=="Student", yes=1, no=0)
dat_rally$Convenience <- ifelse(dat_rally$Sample=="Convenience", yes=1, no=0)

#fit model with effect coding
set.seed(1234)
model_rf_sample <- meta3(y=Fisher, v=Variance_F,
                         intercept.constraints = 0,
                         x=cbind(General, Student, Convenience),
                         cluster=ID_R,
                         data = dat_rally,
                         intervals.type = "LB",
                         model.name = "rf sample")
summary(model_rf_sample) #(mainly reported in SI Table C.3)
FisherZInv(coef(model_rf_sample)[1]) #General pop
FisherZInv(coef(model_rf_sample)[2]) #Students
FisherZInv(coef(model_rf_sample)[3]) #Convenience s.
#explained variance (only reported in SI Table C.2)
summary(model_rf_sample)$R2.values[3] + summary(model_rf_sample)$R2.values[6]
#model fit (only reported SI Table C.2)
anova(model_rf_sample, model_rf_LBCI)

#fit model with dummy coding
set.seed(1234)
model_rf_sample_d1 <- meta3(y=Fisher, v=Variance_F,
                            x=cbind(General, Convenience),
                            cluster=ID_R,
                            data = dat_rally,
                            model.name = "rf sample d1")
summary(model_rf_sample_d1)
set.seed(1234)
model_rf_sample_d2 <- meta3(y=Fisher, v=Variance_F,
                            x=cbind(Student, Convenience),
                            cluster=ID_R,
                            data = dat_rally,
                            model.name = "rf sample d2")
summary(model_rf_sample_d2)


### Location

#check data
summary(as.factor(dat_rally$Country))
prop.table(table(
  as.factor(dat_rally$Country)))
length(unique(as.factor(dat_rally$Country)))
#recode data
dat_rally$Country2[dat_rally$Country=="US"] <- "US"
dat_rally$Country2[dat_rally$Country=="Israel"] <- "Israel"
dat_rally$Country2[is.na(dat_rally$Country2)] <- "Other"
dat_rally$Country2 <- as.factor(dat_rally$Country2)
table(dat_rally$Country2)
prop.table(table(dat_rally$Country2))
#convert characters into dummy variables
dat_rally$US <- ifelse(dat_rally$Country=="US", yes=1, no=0)
dat_rally$Israel <- ifelse(dat_rally$Country=="Israel", yes=1, no=0)
dat_rally$Other_c <- ifelse(dat_rally$US==0 
                                   & dat_rally$Israel==0, yes=1, no=0)

#fit model with effect coding
set.seed(1234)
model_rf_geo <- meta3(y=Fisher, v=Variance_F,
                      intercept.constraints = 0,
                      x=cbind(US, Israel, Other_c),
                      cluster=ID_R,
                      data = dat_rally,
                      intervals.type = "LB",
                      model.name = "rf geo")
summary(model_rf_geo) #(mainly reported in SI Table C.3)
FisherZInv(coef(model_rf_geo)[1]) #US (reported in text)
FisherZInv(coef(model_rf_geo)[2]) #Israel (reported in text)
FisherZInv(coef(model_rf_geo)[3]) #Other countries (reported in text)
#explained variance (only reported in SI Table C.2)
summary(model_rf_geo)$R2.values[3] + summary(model_rf_geo)$R2.values[6]
#model fit (only reported SI Table C.2)
anova(model_rf_geo, model_rf_LBCI)

#fit model with dummy coding
set.seed(1234)
model_rf_geo_d1 <- meta3(y=Fisher, v=Variance_F,
                         x=cbind(US, Israel),
                         cluster=ID_R,
                         data = dat_rally,
                         model.name = "rf geo d1")
summary(model_rf_geo_d1)
set.seed(1234)
model_rf_geo_d2 <- meta3(y=Fisher, v=Variance_F,
                         x=cbind(Israel, Other_c),
                         cluster=ID_R,
                         data = dat_rally,
                         model.name = "rf geo d2")
summary(model_rf_geo_d2)
#adjust p-values for multiple comparisons
p <- cbind(summary(model_rf_geo_d1)$coefficients$`Pr(>|z|)`[2],#US vs. other: yes
           summary(model_rf_geo_d1)$coefficients$`Pr(>|z|)`[3],#Israel vs. other: no
           summary(model_rf_geo_d2)$coefficients$`Pr(>|z|)`[2])#US vs. Israel
p.adjust(p, "BH")


### Mean Age
#check data
summary(dat_rally$AgeMean) 
#fit model with age centered, but not standardized
set.seed(1234)
model_rf_age <- meta3(y=Fisher, v=Variance_F,
                         x=scale(AgeMean, center = T, scale = F),
                         cluster=ID_R,
                         data = dat_rally,
                         intervals.type = "LB",
                         model.name = "rf age")
summary(model_rf_age)


### Percentage Women
#check data
summary(dat_rally$FemaleP) 
#fit model and check model fit improvement
set.seed(1234)
model_rf_fem <- meta3(y=Fisher, v=Variance_F,
                      x=scale(FemaleP, center = T, scale = F),
                      cluster=ID_R,
                      data = dat_rally,
                      intervals.type = "LB",
                      model.name = "rf fem")
summary(model_rf_fem)


### Year of Data Collection
#check data
summary(dat_rally$StudyYear)
#fit model and check model fit improvement
set.seed(1234)
model_rf_datayear <- meta3(y=Fisher, v=Variance_F,
                      x=scale(StudyYear, center = T, scale = F),
                      cluster=ID_R,
                      data = dat_rally,
                      intervals.type = "LB",
                      model.name = "rf data year")
summary(model_rf_datayear)
#explained variance (only reported in SI Table C.2)
summary(model_rf_datayear)$R2.values[3] + summary(model_rf_datayear)$R2.values[6]


### Year of Publication
#check data
summary(dat_rally$Year)
#fit model and check model fit improvement
set.seed(1234)
model_rf_year <- meta3(y=Fisher, v=Variance_F,
                       x=scale(Year, center = T, scale = F),
                       cluster=ID_R,
                       data = dat_rally, 
                       intervals.type = "LB",
                       model.name = "rf year")
summary(model_rf_year)
#explained variance (only reported in SI Table C.2)
summary(model_rf_year)$R2.values[3] + summary(model_rf_year)$R2.values[6]
#model fit
anova(model_rf_year, model_rf_LBCI)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
################################# FIGURE 3 ###############################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
theme_set(theme_bw() + 
            theme(panel.grid.major.y = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  axis.text = element_text(size = (14), colour = "black"),
                  axis.ticks.length = unit(0.3, "cm"),
                  axis.title.y = element_text(size = (14), colour = "black")))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
##### Fig3A. Outgroup Hostility Hypothesis ##### 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
##### Fig3A.1 Overall ##### 
set.seed(1234)
og_overall <- rma.mv(y=Fisher, V=Variance_F,
                     random = ~ 1 | ID_R/ID_ES_Unique, data=dat_outgroup,
                     method = "ML")
og_overall_df <- tidy(og_overall) # create data.frame of regression results
og_overall_df # a tidy data.frame available for dwplot
og_overall_df <- og_overall_df %>%  #re-label variable
  relabel_predictors(overall = "Overall Fisher's Z \n Correlation Coef." )           
og_overall_plot <- dwplot(og_overall_df, #plot
                          dot_args = list(size = 3),
                          whisker_args = list(size = 1),
                          vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2, size = 1))+
  ggtitle("Outgroup Hostility") +
  scale_color_grey() +
  xlim(-0.1, .25) +
  theme(axis.title.x=element_blank(),
        plot.margin=grid::unit(c(0,0,0,0), "mm"),
        axis.text.y = element_text(face = "bold",size = 14),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        panel.background = element_rect(fill = "grey90",
                                        colour = "grey90"),
        panel.grid = element_line(colour = "grey60"))


##### Fig3A.2 Terror type ##### 
#fit model and plot results
set.seed(1234)
og_terror <- rma.mv(y=Fisher, V=Variance_F, mods = cbind(Islam, Other_Ideology, No_Ideology),
                    intercept=F, 
                    random = ~ 1 | ID_R/ID_ES_Unique, data=dat_outgroup,
                    method = "ML")
summary(og_terror)
og_terror_df <- tidy(og_terror) # create data.frame of regression results
og_terror_df # a tidy data.frame available for dwplot
og_terror_df <- og_terror_df %>%  #re-label variable
  relabel_predictors(c(Islam = "Islamist",
                       Other_Ideology  = "Other ideology",
                       No_Ideology  = "No ideology"))
og_terror_plot <- dwplot(og_terror_df, #plot
                         dot_args = list(size = 2),
                         whisker_args = list(size = 0.5),
                         vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2, size = 0.7)) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        plot.margin=grid::unit(c(0,0,0,0), "mm")) +
  scale_color_grey() +
  xlim(-0.1, .25) 

##### Fig3A.3 Research Design ##### 
#fit model and plot results
set.seed(1234)
og_design <- rma.mv(y=Fisher, V=Variance_F, mods = cbind(Exp, QuasiLong, Corr),
                    intercept=F, 
                    random = ~ 1 | ID_R/ID_ES_Unique, data=dat_outgroup,
                    method = "ML")
summary(og_design)
og_design_df <- tidy(og_design) # create data.frame of regression results
og_design_df # a tidy data.frame available for dwplot
og_design_df <- og_design_df %>%  #re-label variable
  relabel_predictors(c(Exp = "Experiment",
                       Corr  = "Cross-sectional",
                       QuasiLong  = "Other designs"))
og_design_plot <- dwplot(og_design_df, #plot
                         dot_args = list(size = 2),
                         whisker_args = list(size = 0.5),
                         vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2, size = 0.7)) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        plot.margin=grid::unit(c(0,0,0,0), "mm")) +
  scale_color_grey() +
  xlim(-0.1, .25) 

##### Fig3A.4 Sample  ##### 
#fit model and plot results
set.seed(1234)
og_sample <- rma.mv(y=Fisher, V=Variance_F, mods = cbind(General, Student, Convenience),
                    intercept=F, 
                    random = ~ 1 | ID_R/ID_ES_Unique, data=dat_outgroup,
                    method = "ML")
summary(og_sample)
og_sample_df <- tidy(og_sample) # create data.frame of regression results
og_sample_df # a tidy data.frame available for dwplot
og_sample_df <- og_sample_df %>%  #re-label variable
  relabel_predictors(c(General = "General population",
                       Student  = "Student sample",
                       Convenience  = "Convenience sample"))
og_sample_plot <- dwplot(og_sample_df, #plot
                         dot_args = list(size = 2),
                         whisker_args = list(size = 0.5),
                         vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2, size = 0.7)) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        plot.margin=grid::unit(c(0,0,0,0), "mm")) +
  scale_color_grey() +
  xlim(-0.1, .25) 

##### Fig3A.5 Geography ##### 
#fit model and plot results
set.seed(1234)
og_location <- rma.mv(y=Fisher, V=Variance_F, mods = cbind(US, Israel, Other_c),
                    intercept=F, 
                    random = ~ 1 | ID_R/ID_ES_Unique, data=dat_outgroup,
                    method = "ML")
summary(og_location)
og_location_df <- tidy(og_location) # create data.frame of regression results
og_location_df # a tidy data.frame available for dwplot
og_location_df <- og_location_df %>%  #re-label variable
  relabel_predictors(c(US = "United States",
                       Israel  = "Israel",
                       Other_c  = "Other countries"))
og_location_plot <- dwplot(og_location_df, #plot
                         dot_args = list(size = 2),
                         whisker_args = list(size = 0.5),
                         vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2, size = 0.7)) +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) +
  scale_color_grey() +
  xlim(-0.1, .25) 

##### Fig3A.6 Guilt-by-Association ##### 
#fit model and plot results
set.seed(1234)
og_association <- rma.mv(y=Fisher, V=Variance_F, mods = cbind(Lot, Bit, No),
                      intercept=F, 
                      random = ~ 1 | ID_R/ID_ES_Unique, data=dat_outgroup,
                      method = "ML")
summary(og_association)
og_association_df <- tidy(og_association) # create data.frame of regression results
og_association_df # a tidy data.frame available for dwplot
og_association_df <- og_association_df %>%  #re-label variable
  relabel_predictors(c(Lot = "Strong association",
                       Bit  = "Moderate association",
                       No  = "No association"))
og_association_plot <- dwplot(og_association_df, #plot
                           dot_args = list(size = 2),
                           whisker_args = list(size = 0.5),
                           vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2, size = 0.7)) +
  ggtitle("Guilt-By-Association") +
  scale_y_discrete(position = "left") +
  scale_color_grey() +
  theme(axis.title=element_text(face = "bold"),
        plot.margin=grid::unit(c(0,0,0,0), "mm"))+
  xlim(-0.1, .25) 

##### Fig3A Outgroup Hostility Effect Sizes Plot ##### 
og_plots <- plot_grid(og_overall_plot, og_terror_plot, og_design_plot, og_sample_plot, og_location_plot,
                      og_association_plot,
                      ncol=1, nrow=6, 
                      rel_heights=c(0.8,1,1,1,1,1), align = "v")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
##### Fig3B. Conservative Shift Hypothesis ##### 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
##### Fig3B.1 Overall ##### 
set.seed(1234)
cs_overall <- rma.mv(y=Fisher, V=Variance_F,
                     random = ~ 1 | ID_R/ID_ES_Unique, data=dat_conservatism,
                     method = "ML")
cs_overall_df <- tidy(cs_overall) # create data.frame of regression results
cs_overall_df # a tidy data.frame available for dwplot
cs_overall_df <- cs_overall_df %>% 
  relabel_predictors(c(overall = "Overall Effect Size")) #re-label variable
cs_overall_plot <- dwplot(cs_overall_df, #plot
                          dot_args = list(size = 3),
                          whisker_args = list(size = 1),
                          vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2, size = 1))+
  ggtitle("Conservative Shift") +
  scale_color_grey() +
  xlim(-0.1, .25) +
  theme(axis.title.x=element_blank(),
        plot.margin=grid::unit(c(0,0,0,0), "mm"),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        panel.background = element_rect(fill = "grey90",
                                        colour = "grey90"),
        panel.grid = element_line(colour = "grey60"))

##### Fig3B.2 Terror type ##### 
#fit model and plot results
set.seed(1234)
cs_terror <- rma.mv(y=Fisher, V=Variance_F, mods = cbind(Islam, Other_Ideology, No_Ideology),
                    intercept=F, 
                    random = ~ 1 | ID_R/ID_ES_Unique, data=dat_conservatism,
                    method = "ML")
summary(cs_terror)
cs_terror_df <- tidy(cs_terror) # create data.frame of regression results
cs_terror_df # a tidy data.frame available for dwplot
cs_terror_df <- cs_terror_df %>%  #re-label variable
  relabel_predictors(c(Islam = "Islamist",
                       Other_Ideology  = "Other Ideology",
                       No_Ideology  = "No Ideology"))
cs_terror_plot <- dwplot(cs_terror_df, #plot
                         dot_args = list(size = 2),
                         whisker_args = list(size = 0.5),
                         vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2, size = 0.7)) +
  theme(plot.title = element_text(size = 14), axis.text=element_blank(), axis.ticks=element_blank(),
        plot.margin=grid::unit(c(0,0,0,0), "mm")) +
  scale_color_grey() +
  xlim(-0.1, .25) 

##### Fig3B.3 Research Design ##### 
#fit model and plot results
set.seed(1234)
cs_design <- rma.mv(y=Fisher, V=Variance_F, mods = cbind(Exp, QuasiLong, Corr),
                    intercept=F, 
                    random = ~ 1 | ID_R/ID_ES_Unique, data=dat_conservatism,
                    method = "ML")
summary(cs_design)
cs_design_df <- tidy(cs_design) # create data.frame of regression results
cs_design_df # a tidy data.frame available for dwplot
cs_design_df <- cs_design_df %>%  #re-label variable
  relabel_predictors(c(Exp = "Experiments",
                       Corr  = "Cross-Sectional",
                       QuasiLong  = "Other Designs"))
cs_design_plot <- dwplot(cs_design_df, #plot
                         dot_args = list(size = 2),
                         whisker_args = list(size = 0.5),
                         vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2, size = 0.7)) +
  theme(plot.title = element_text(size = 14), axis.text=element_blank(), 
        axis.ticks=element_blank(),
        plot.margin=grid::unit(c(0,0,0,0), "mm")) +
  scale_color_grey() +
  xlim(-0.1, .25) 

##### Fig3B.4 Sample  ##### 
#fit model and plot results
set.seed(1234)
cs_sample <- rma.mv(y=Fisher, V=Variance_F, mods = cbind(General, Student, Convenience),
                    intercept=F, 
                    random = ~ 1 | ID_R/ID_ES_Unique, data=dat_conservatism,
                    method = "ML")
summary(cs_sample)
cs_sample_df <- tidy(cs_sample) # create data.frame of regression results
cs_sample_df # a tidy data.frame available for dwplot
cs_sample_df <- cs_sample_df %>%  #re-label variable
  relabel_predictors(c(General = "General Population",
                       Student  = "Student Samples",
                       Convenience  = "Convenience Samples"))
cs_sample_plot <- dwplot(cs_sample_df, #plot
                         dot_args = list(size = 2),
                         whisker_args = list(size = 0.5),
                         vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2, size = 0.7)) +
  theme(plot.title = element_text(size = 14), 
        axis.text=element_blank(), axis.ticks=element_blank(),
        plot.margin=grid::unit(c(0,0,0,0), "mm")) +
  scale_color_grey() +
  xlim(-0.1, .25) 

##### Fig3B.5 Geography ##### 
#fit model and plot results
set.seed(1234)
cs_location <- rma.mv(y=Fisher, V=Variance_F, mods = cbind(US, Israel, Other_c),
                      intercept=F, 
                      random = ~ 1 | ID_R/ID_ES_Unique, data=dat_conservatism,
                      method = "ML")
summary(cs_location)
cs_location_df <- tidy(cs_location) # create data.frame of regression results
cs_location_df # a tidy data.frame available for dwplot
cs_location_df <- cs_location_df %>%  #re-label variable
  relabel_predictors(c(US = "United States",
                       Israel  = "Israel",
                       Other_c  = "Other Countries"))
cs_location_plot <- dwplot(cs_location_df, #plot
                           dot_args = list(size = 2),
                           whisker_args = list(size = 0.5),
                           vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2, size = 0.7)) +
  theme(plot.title = element_text(size = 14), 
        axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        plot.margin=grid::unit(c(0,0,0,0), "mm")) +
  scale_color_grey() +
  xlim(-0.1, .25) 

##### Fig3B Conservative Shift Effect Sizes Plot ##### 
cs_plots <- plot_grid(cs_overall_plot, cs_terror_plot, cs_design_plot, cs_sample_plot, cs_location_plot, 
                      ncol=1, nrow=6, 
                      rel_heights=c(0.8,1,1,1,1,1), align = "v")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
##### Fig3C. Rally Effects Hypothesis ##### 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
##### Fig3C.1 Overall ##### 
set.seed(1234)
rf_overall <- rma.mv(y=Fisher, V=Variance_F,
                     random = ~ 1 | ID_R/ID_ES_Unique, data=dat_rally,
                     method = "ML")
rf_overall_df <- tidy(rf_overall) # create data.frame of regression results
rf_overall_df # a tidy data.frame available for dwplot
rf_overall_df <- rf_overall_df %>% 
  relabel_predictors(c(overall = "Overall Effect Size")) #re-label variable
rf_overall_plot <- dwplot(rf_overall_df, #plot
                          dot_args = list(size = 3),
                          whisker_args = list(size = 1),
                          vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2, size = 1))+
  ggtitle("Rally Effects") +
  scale_color_grey() +
  xlim(-0.1, .25) +
  theme(axis.title.x=element_blank(),
        plot.margin=grid::unit(c(0,0,0,0), "mm"),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
        panel.background = element_rect(fill = "grey90",
                                        colour = "grey90"),
        panel.grid = element_line(colour = "grey60"))


##### Fig3C.2 Terror type ##### 
#fit model and plot results
set.seed(1234)
rf_terror <- rma.mv(y=Fisher, V=Variance_F, mods = cbind(Islam, Other_Ideology, No_Ideology),
                    intercept=F, 
                    random = ~ 1 | ID_R/ID_ES_Unique, data=dat_rally,
                    method = "ML")
summary(rf_terror)
rf_terror_df <- tidy(rf_terror) # create data.frame of regression results
rf_terror_df # a tidy data.frame available for dwplot
rf_terror_df <- rf_terror_df %>%  #re-label variable
  relabel_predictors(c(Islam = "Islamist",
                       Other_Ideology  = "Other Ideology",
                       No_Ideology  = "No Ideology"))
rf_terror_plot <- dwplot(rf_terror_df, #plot
                         dot_args = list(size = 2),
                         whisker_args = list(size = 0.5),
                         vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2, size = 0.7)) +
  ylab("Ideology") +
  scale_y_discrete(position = "right") +
  scale_color_grey() +
  theme(axis.text=element_blank(), axis.ticks=element_blank(), 
        axis.title.y=element_text(face = "bold", size = 14),
        plot.margin=grid::unit(c(0,0,0,0), "mm"))+
  xlim(-0.1, .25) 


##### Fig3C.3 Research Design ##### 
#fit model and plot results
set.seed(1234)
rf_design <- rma.mv(y=Fisher, V=Variance_F, mods = cbind(Exp, QuasiLong, Corr),
                    intercept=F, 
                    random = ~ 1 | ID_R/ID_ES_Unique, data=dat_rally,
                    method = "ML")
summary(rf_design)
rf_design_df <- tidy(rf_design) # create data.frame of regression results
rf_design_df # a tidy data.frame available for dwplot
rf_design_df <- rf_design_df %>%  #re-label variable
  relabel_predictors(c(Exp = "Experiments",
                       Corr  = "Cross-Sectional",
                       QuasiLong  = "Other Designs"))
rf_design_plot <- dwplot(rf_design_df, #plot
                         dot_args = list(size = 2),
                         whisker_args = list(size = 0.5),
                         vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2, size = 0.7)) +
  theme(axis.text=element_blank(), axis.ticks=element_blank())+
  ylab("Design") +
  scale_y_discrete(position = "right") +
  scale_color_grey() +
  theme(axis.text=element_blank(), axis.ticks=element_blank(), 
        axis.title.y=element_text(face = "bold", size = 14),
        plot.margin=grid::unit(c(0,0,0,0), "mm"))+
  xlim(-0.1, .25) 

##### Fig3C.4 Sample  ##### 
#fit model and plot results
set.seed(1234)
rf_sample <- rma.mv(y=Fisher, V=Variance_F, mods = cbind(General, Student, Convenience),
                    intercept=F, 
                    random = ~ 1 | ID_R/ID_ES_Unique, data=dat_rally,
                    method = "ML")
summary(rf_sample)
rf_sample_df <- tidy(rf_sample) # create data.frame of regression results
rf_sample_df # a tidy data.frame available for dwplot
rf_sample_df <- rf_sample_df %>%  #re-label variable
  relabel_predictors(c(General = "General Population",
                       Student  = "Student Samples",
                       Convenience  = "Convenience Samples"))
rf_sample_plot <- dwplot(rf_sample_df, #plot
                         dot_args = list(size = 2),
                         whisker_args = list(size = 0.5),
                         vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2, size = 0.7)) +
  theme(axis.text=element_blank(), axis.ticks=element_blank())+
  ylab("Sample") +
  scale_y_discrete(position = "right") +
  scale_color_grey() +
  theme(axis.text=element_blank(), axis.ticks=element_blank(), 
        axis.title.y=element_text(face = "bold", size = 14),
        plot.margin=grid::unit(c(0,0,0,0), "mm"))+
  xlim(-0.1, .25) 


##### Fig3C.5 Geography ##### 
#fit model and plot results
set.seed(1234)
rf_location <- rma.mv(y=Fisher, V=Variance_F, mods = cbind(US, Israel, Other_c),
                      intercept=F, 
                      random = ~ 1 | ID_R/ID_ES_Unique, data=dat_rally,
                      method = "ML")
summary(rf_location)
rf_location_df <- tidy(rf_location) # create data.frame of regression results
rf_location_df # a tidy data.frame available for dwplot
rf_location_df <- rf_location_df %>%  #re-label variable
  relabel_predictors(c(US = "United States",
                       Israel  = "Israel",
                       Other_c  = "Other Countries"))
rf_location_plot <- dwplot(rf_location_df, #plot
                           dot_args = list(size = 2),
                           whisker_args = list(size = 0.5),
                           vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2, size = 0.7)) +
  ylab("Country") +
  scale_y_discrete(position = "right") +
  scale_color_grey() +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), 
        axis.title.y=element_text(face = "bold", size = 14),
        plot.margin=grid::unit(c(0,0,0,0), "mm"))+
  xlim(-0.1, .25) 

##### Fig3C Rally Effects Effect Sizes Plot ##### 
rf_plots <- plot_grid(rf_overall_plot, rf_terror_plot, rf_design_plot, rf_sample_plot, rf_location_plot,
                      ncol=1, nrow=6, 
                      rel_heights=c(0.8,1,1,1,1,1), align = "v")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
################################# FIGURE 3 ###############################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
# 1. Open jpeg file
jpeg("Figures/Fig3.jpeg", width = 10, height = 11.3, units = 'in', res = 400)
# 2. Create the plot
ggarrange(og_plots, cs_plots, rf_plots,
          align = "hv",
          ncol = 3, nrow = 1,
          widths=c(2.2,1.3,1.4)) #arrange all plots
# 3. Close the file
dev.off()
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
################################# FIGURE 4 ###############################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
dat_rally2 <- dat_rally[ which(dat_rally$Year > 2000), ] #select data
dat_rally2$Country2 <- factor(dat_rally2$Country2, 
                             levels = c("US", "Israel", "Other")) #relevel data
dat_rally2$Nine11 <- factor(dat_rally2$Nine11,                 
                           levels = c("1", "0")) #relevel data
mytheme2 <- theme(plot.title = element_text(size = (14), face = "bold", colour = "black"),
                  axis.text = element_text(size = (12), colour = "black"),
                  axis.ticks.length = unit(0.4, "cm"),
                  axis.title = element_text(size = (14), colour = "black"),
                  legend.title = element_text(size = 11),
                  legend.text = element_text(size = 11)) #set theme for plots
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
fig4A <- ggplot(dat_rally2, aes(x=Year, y=Fisher, color=Nine11, shape=Nine11)) +
  geom_point(size=2.5) + 
  geom_smooth(method=lm, se=F, fullrange=F)+
  scale_shape_manual(values=c(1,4),
                     labels=c("Yes", "No"))+
  scale_color_grey(start = 0, end = 0.7,
                     labels=c("Yes", "No"))+
  labs(title = "Evolution over time for rally effects, by reference to 9/11", 
       x = "Year of publication\n", 
       y = "Fisher's Z correlation coefficient", 
       shape="Reference to 9/11", colour="Reference to 9/11") +
  scale_y_continuous(breaks = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1)) +
  scale_x_continuous(breaks = seq(2000, 2020, by = 2), limits = c(2004, 2020)) +
  theme_classic()
fig4A <- fig4A + mytheme2+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))
fig4A
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
fig4B <- ggplot(dat_rally2, aes(x=Year, y=Fisher, color=TerrorType2, shape=TerrorType2)) +
  geom_point(size=2.5) + 
  geom_smooth(method=lm, se=F, fullrange=F)+
  scale_shape_manual(values=c(1,4,17),
                     labels=c("Islamist", "No ideology", "Other ideology"))+
  scale_color_grey(start = 0, end = 0.8,
                   labels=c("Islamist", "No ideology", "Other ideology"))+
  labs(title = "Evolution over time for rally effects, by ideology of terrorism", 
       x = "Year of publication", 
       y = "Fisher's Z correlation coefficient", 
       shape="Ideology", colour="Ideology") +
  scale_y_continuous(breaks = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1)) +
  scale_x_continuous(breaks = seq(2000, 2020, by = 2), limits = c(2004, 2020)) +
  theme_classic()
fig4B <- fig4B + mytheme2 +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))
fig4B
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
# 1. Open jpeg file
jpeg("Figures/Fig4.jpeg", width = 7, height = 9, units = 'in', res = 300)
# 2. Create the plot
ggarrange(fig4A, fig4B,
          labels = c("(A)", "(B)"),
          #font.label = list(size = 22, color = "black"),
          ncol = 1, nrow = 2) #arrange all plots
# 3. Close the file
dev.off()
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
