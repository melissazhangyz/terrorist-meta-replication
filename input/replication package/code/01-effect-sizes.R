# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
# PROJECT:    How Terrorism Does (and Does Not) Affect Citizens' 
#             Political Attitudes: A Meta-Analysis. 
# AUTHOR:     Amelie Godefroidt
# CONTACT:    amelie.godefroidt@ntnu.no
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
# This R file contains the code necessary to compute the common effect sizes
# 
# Last successful replication: 01.11.2021
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 

# first things first: set your settings.

rm(list=ls())
setwd(".../Godefroidt_2021_ReplicationFiles/") #set your directory
getwd()

###### Install and load packages ###### 

# ipak function: install and load multiple R packages.
# check to see if packages are installed. Install them if they are not, then load them into the R session.

ipak <- function(pkg){  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
if(length(new.pkg)) install.packages(new.pkg, dependencies=TRUE)
sapply(pkg, require, character.only=TRUE)
}

packages <- c("esc", "compute.es", "metaSEM", "data.table", "pastecs", "ltm", 
              "pwr", "psych", "tidyr", "haven", "dplyr", "Hmisc", "srvyr", "readxl") 
ipak(packages)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 
###### A. Functions to calculate effect sizes ###### 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ # 

#dataframe to collect effect sizes
Effectsizes <- data.frame(matrix(ncol=6, nrow=0))

#function to calculate (most) effect sizes
effect_size <- function( eff.type, n, n.1, n.2, 
                         f,
                         r, 
                         R, q, 
                         lor, var.lor, 
                         m.1, m.2, sd.1, sd.2,
                         p1, p2, n.ab, n.cd,
                         t,
                         p, tail,
                         d,
                         chi.sq){
  
  # ANCOVA F-test
  if (eff.type == "a.fes"){
    es <- a.fes(f, n.1, n.2, R, q, dig=6)
  }
  # Correlation
  if (eff.type == "corr"){ 
    es <- res(r, var.r=NULL, n, dig=6)
  }
  # F-test
  if (eff.type == "f.test"){ 
    es <- fes(f, n.1, n.2, dig=6)
  }
  # Log Odds Ratio: Logit
  if (eff.type == "lor"){
    es <- lores(lor, var.lor, n.1, n.2, dig=6)
  }
  # Means
  if (eff.type == "means"){ 
    es <- mes(m.1, m.2, sd.1, sd.2, n.1, n.2, dig=6)
  }
  # Proportions
  if (eff.type == "prop"){ 
    es <- propes(p1, p2, n.ab, n.cd, dig=6)
  }
  # T-test
  if (eff.type == "t.test"){ 
    es <- tes(t, n.1, n.2, dig=6)
  }
  # P-value
  if (eff.type == "p.val"){ 
    es <- pes(p, n.1, n.2, tail, dig=6)
  }
  # Cohen's D
  if (eff.type == "coh.d"){ 
    es <- des(d, n.1, n.2, dig=6)
  }
  # Chi square
  if (eff.type == "chi.sq"){ 
    es <- chies(chi.sq, n, dig=6)
  }
  
  # standard error of the correlation coefficient is the square root of its variance
  st.err.r <- round(sqrt(es$var.r), digits=6)
  
  # standard error of Fisher's Z is the square root of its variance
  st.err.z <- round(sqrt(es$var.z), digits=6)
  
  # post-hoc power calculation
  pwr <- pwr.r.test(n=es$N.total, r=es$r, sig.level=.05)
  
  # print everything in data.frame "Effectsizes" created above
  results <- data.frame(c(es$r, es$var.r, st.err.r,
                          es$fisher.z, es$var.z, st.err.z, es$N.total, pwr$power))
  results <- transpose(results)
  rbind(Effectsizes, results)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#### 0n the use of Beta Coefficients in Meta-Analysis ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Note 1: 
# corrected beta coefficients from multivariate regressions are used  
# when correlation matrix is not reported or obtained via e-mail using:
# r=0.98 beta + 0.05 lambda, where lambda is an indicator variable 
# that equals 1 when beta is nonnegative and 0 when beta is negative.

beta2cor <- function(b,
                     l) {
   return((0.98 * b) + (0.05 * l))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Note 2: 
# b_effect_size function is used when b/beta coefficients are 
# reported using experimental design/two-group designs

b_effect_size <- function(eff.type, 
                          b, beta, sdy,
                          grp1n, grp2n){ 
  
  # Unstandardized Regression Coefficient
  if (eff.type == "b"){ 
    es <- esc_B(b, sdy, grp1n, grp2n, es.type="r")
  }
  
  # Unstandardized Regression Coefficient
  if (eff.type == "beta"){ 
    es <- esc_beta(beta, sdy, grp1n, grp2n, es.type="r")
  }
  
  # standard error of Fisher's Z
  z.se <- round(1/sqrt((es$totaln)-3), digits=6)
  
  # variance of Fisher's Z
  z.var <- z.se^2
  
  # power calculation
  pwr <- pwr.r.test(n=es$totaln, r=es$es, sig.level=.05)
  
  # print everything in data.frame "Effectsizes" created above
  results <- data.frame(c(es$es, es$var, es$se,
                          es$zr, z.var, z.se, es$totaln, pwr$power))
  results <- transpose(results)
  rbind(Effectsizes, results)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#### 0n the use of Logit and Probit Coefficients in Meta-Analysis ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Note 3: 
## First, convert probit into OR via logit transformation based on
## Birnbaum, A. (1968). Some latent trait models and their use in inferring an examinee's ability. In F.M. Lord and M.R. Novick (Eds.), Statistical theries of mental test scores. Addison-Wesly. 395-479.
##### OR <- exp(1.7 * probit) !
## Second, convert odd ratio to pearson correlation based on
## Bonett, D. G. (2007). Transforming odds ratios into correlations for meta-analytic research American Psychologist, 62(3), 254-255. doi:10.1037/0003-066X.62.3.254:

OR2cor <- function(OR) {
  return(pearson=cos(pi/(1 + sqrt(OR))))
}

## Hence, add OR2cor(exp(1.7 * probit)) in function for correlation for probit models
## Hence, add OR2cor(exp(logit)) in function for correlation for logit models


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#### B. Calculating the ES's and binding them into the dataframe #### 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Aber et al.(2004)
Effectsizes <- effect_size(eff.type="corr", r=-0.06, n=768)
Effectsizes <- effect_size(eff.type="corr", r=-0.03, n=768)
Effectsizes <- effect_size(eff.type="corr", r=-0.07, n=768)

# Abrams et. al. (2017)
Effectsizes <- effect_size(eff.type="f.test", f=17.53, n.1=1068, n.2=869)

# Aizapurua et al. (2017)
Effectsizes <- effect_size(eff.type="corr", r=0.16, n=418)
Effectsizes <- effect_size(eff.type="corr", r=0.08, n=418)
Effectsizes <- effect_size(eff.type="corr", r=0.22, n=411)

# Anonymous (2012)
# Ideology
Effectsizes <- effect_size(eff.type="corr", r= 0.01, n=297)
Effectsizes <- effect_size(eff.type="corr", r= 0.03, n=297)
Effectsizes <- effect_size(eff.type="corr", r=-0.05, n=297)
Effectsizes <- effect_size(eff.type="corr", r= 0.12, n=297)
Effectsizes <- effect_size(eff.type="corr", r= 0.01, n=297)
# Institutional Trust
Effectsizes <- effect_size(eff.type="corr", r=-0.03, n=297)
Effectsizes <- effect_size(eff.type="corr", r= 0.15, n=297)
Effectsizes <- effect_size(eff.type="corr", r= 0.05, n=297)
Effectsizes <- effect_size(eff.type="corr", r= 0.13, n=297)
Effectsizes <- effect_size(eff.type="corr", r= 0.02, n=297)
# Immigration policies
Effectsizes <- effect_size(eff.type="corr", r= 0.05, n=297)
Effectsizes <- effect_size(eff.type="corr", r=-0.17, n=297)
Effectsizes <- effect_size(eff.type="corr", r=-0.04, n=297)
Effectsizes <- effect_size(eff.type="corr", r= 0.07, n=297)
Effectsizes <- effect_size(eff.type="corr", r= 0.03, n=297)

# Argyrides & Downey (2012)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=20.65, m.2=17.53,
                           sd.1=7.91, sd.2=4.67,
                           n.1 =34, n.2=30)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=18.20, m.2=17.53,
                           sd.1=4.58, sd.2=4.67,
                           n.1=30, n.2 =30)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=19.61, m.2=17.53,
                           sd.1=5.79, sd.2=4.67,
                           n.1=32, n.2 =30)
Effectsizes <- effect_size(eff.type="means",
                           m.1=16.19, m.2=17.53,
                           sd.1=4.80, sd.2=4.67, 
                           n.1=32, n.2 =30)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=19.68, m.2=17.53, 
                           sd.1=5.15, sd.2=4.67, 
                           n.1 =45, n.2=30)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=19.81, m.2=17.53, 
                           sd.1=6.46, sd.2=4.67, 
                           n.1 =40, n.2=30)

# Arian (1989)
Effectsizes <- effect_size(eff.type="corr", r=0.43, n=1116)
Effectsizes <- effect_size(eff.type="corr", r=0.12, n=1116)

# Arvanitidis, Economou, & Kollias (2016)
Effectsizes <- effect_size(eff.type="corr", r=0.0293, n=1148)
Effectsizes <- effect_size(eff.type="corr", r=-0.0108, n=1303)
Effectsizes <- effect_size(eff.type="corr", r=0.0243, n=1593)
Effectsizes <- effect_size(eff.type="corr", r=-0.0673, n=250)

# Asbrock & Fritsche (2013)
# a) Study 1
Effectsizes <- effect_size(eff.type="means", 
                           m.1=-0.72, m.2=-1.37, 
                           sd.1=1.49, sd.2=1.75, 
                           n.1=48, n.2 =48)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=-1.52, m.2=-1.37, 
                           sd.1=1.83, sd.2=1.75, 
                           n.1=48, n.2 =48)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.30, m.2=2.89, 
                           sd.1=0.73, sd.2=0.91, 
                           n.1=48, n.2 =48)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=2.97, m.2=2.89, 
                           sd.1=0.86, sd.2=0.91, 
                           n.1=48, n.2 =48)
# b) Study 2
Effectsizes <- effect_size(eff.type="corr", r=0.19, n=99)
Effectsizes <- effect_size(eff.type="corr", r=0.14, n=99)
Effectsizes <- effect_size(eff.type="corr", r=-0.08, n=99)
Effectsizes <- effect_size(eff.type="corr", r=0.04, n=99)
Effectsizes <- effect_size(eff.type="corr", r=0.03, n=99)
Effectsizes <- effect_size(eff.type="corr", r=-0.02, n=99)

# Aytac & Carkoglu (2019)
Effectsizes <- effect_size(eff.type="corr", r=(OR2cor(exp(0.15))), n=1776)
Effectsizes <- effect_size(eff.type="corr", r=-(OR2cor(exp(0.001))), n=932)
Effectsizes <- effect_size(eff.type="corr", r=(OR2cor(exp(0.59))), n=1191)
Effectsizes <- effect_size(eff.type="corr", r=(OR2cor(exp(0.65))), n=648)

# Baele et al. (2019)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=5.70625, m.2=5.559524, 
                           sd.1=1.251901, sd.2=1.400181, 
                           n.1=160, n.2=168)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=5.736842, m.2=5.559524, 
                           sd.1=1.454551, sd.2=1.400181, 
                           n.1=152, n.2=168)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=5.158228, m.2=5.01227, 
                           sd.1=1.541453, sd.2=1.507134, 
                           n.1=158, n.2=163)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=5.226667, m.2=5.01227, 
                           sd.1=1.400607, sd.2=1.507134, 
                           n.1=150, n.2=163)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.551515, m.2=3.478261, 
                           sd.1=1.437333, sd.2=1.541456, 
                           n.1=165, n.2=161)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.551515, m.2=3.469799, 
                           sd.1=1.437333, sd.2=1.392825, 
                           n.1=165, n.2=149)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.797297, m.2=3.596386, 
                           sd.1=1.507169, sd.2=1.676862, 
                           n.1=148, n.2=166)
Effectsizes <- effect_size(eff.type="means",
                           m.1=3.748428, m.2=3.596386, 
                           sd.1=1.698728, sd.2=1.676862, 
                           n.1=159, n.2=166)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.26875, m.2=3.473054, 
                           sd.1=1.882464, sd.2=1.887863, 
                           n.1=160, n.2=167)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.337748, m.2=3.473054, 
                           sd.1=1.7393, sd.2=1.887863, 
                           n.1=151, n.2=167)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.41875, m.2=3.214286, 
                           sd.1=2.111579, sd.2=2.021273, 
                           n.1=151, n.2=168)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.556291, m.2=3.214286, 
                           sd.1=2.005445, sd.2=2.021273, 
                           n.1=160, n.2=168)

# Bali (2007)
# a) CIS sample
Effectsizes <- effect_size(eff.type="prop", 
                           p1=0.914, p2=0.869, 
                           n.ab=1554, n.cd=3845) #Turnout
Effectsizes <- effect_size(eff.type="prop", 
                           p1=0.448, p2=0.656, 
                           n.ab=3845, n.cd=1554) #PSOE 
Effectsizes <- effect_size(eff.type="prop", 
                           p1=0.194, p2=0.361, 
                           n.ab=1554, n.cd=3845) #PP
Effectsizes <- effect_size(eff.type="prop", 
                           p1=0.056, p2=0.047, 
                           n.ab=3845, n.cd=1554) #IU 
# b) OPINA sample
Effectsizes <- effect_size(eff.type="prop", 
                           p1=0.389, p2=0.541, 
                           n.ab=716, n.cd=284) #PSOE 
Effectsizes <- effect_size(eff.type="prop", 
                           p1=0.224, p2=0.431, 
                           n.ab=284, n.cd=716) #PP 
Effectsizes <- effect_size(eff.type="prop", 
                           p1=0.045, p2=0.040, 
                           n.ab=284, n.cd=716) #IU

# Bali (2009)
Effectsizes <- effect_size(eff.type="prop", 
                           p2=0.686, p1=0.738, 
                           n.cd=175, n.ab=189) #Overall support
Effectsizes <- effect_size(eff.type="prop", 
                           p2=0.745, p1=0.769, 
                           n.cd=175, n.ab=189) #Documentation
Effectsizes <- effect_size(eff.type="prop", 
                           p2=0.399, p1=0.420, 
                           n.cd=175, n.ab=189) #Costs of validation
Effectsizes <- effect_size(eff.type="prop", 
                           p2=0.835, p1=0.915, 
                           n.cd=175, n.ab=189) #Proof legal status
Effectsizes <- effect_size(eff.type="prop", 
                           p2=0.659, p1=0.714, 
                           n.cd=175, n.ab=189) #Bar codes
Effectsizes <- effect_size(eff.type="prop", 
                           p2=0.708, p1=0.762, 
                           n.cd=175, n.ab=189) #States share

# Barnes et al. (2012)
# a) Study 1
Effectsizes <- effect_size(eff.type="corr", r=0.23, n=150)
Effectsizes <- effect_size(eff.type="corr", r=0.18, n=150)
Effectsizes <- effect_size(eff.type="corr", r=0.29, n=150)
Effectsizes <- effect_size(eff.type="corr", r=0.31, n=150)
Effectsizes <- effect_size(eff.type="corr", r=0.33, n=150)
Effectsizes <- effect_size(eff.type="corr", r=0.36, n=150)
Effectsizes <- effect_size(eff.type="corr", r=0.26, n=150)

# Bar-Tal & Labin (2001)
# a) social-evaluative traits
### 1 day after attacks
Effectsizes <- effect_size(eff.type="means", 
                           m.2=2.71, m.1=2.89, 
                           sd.2=0.70, sd.1=0.58, 
                           n.2=82, n.1 =82) # Palestinians
Effectsizes <- effect_size(eff.type="means", 
                           m.2=2.95, m.1=3.00, 
                           sd.2=0.58, sd.1=0.57, 
                           n.2=82, n.1 =82) # Arabs
Effectsizes <- effect_size(eff.type="means", 
                           m.2=3.34, m.1=3.45, 
                           sd.2=0.53, sd.1=0.59, 
                           n.2=82, n.1 =82) # Jordans
### three months after attacks
Effectsizes <- effect_size(eff.type="means", 
                           m.2=2.82, m.1=2.89, 
                           sd.2=0.61, sd.1=0.58, 
                           n.2=82, n.1 =82) # Palestinians
Effectsizes <- effect_size(eff.type="means", 
                           m.2=2.94, m.1=3.00, 
                           sd.2=0.57, sd.1=0.57, 
                           n.2=82, n.1 =82) # Arabs
Effectsizes <- effect_size(eff.type="means", 
                           m.2=3.29, m.1=3.45, 
                           sd.2=0.49, sd.1=0.59, 
                           n.2=82, n.1 =82) # Jordans
# b) Potency traits
### 1 day after attacks
Effectsizes <- effect_size(eff.type="means", 
                           m.2=3.27, m.1=3.37, 
                           sd.2=0.78, sd.1=0.67, 
                           n.2=82, n.1 =82) # Palestinians
Effectsizes <- effect_size(eff.type="means", 
                           m.2=3.32, m.1=3.45, 
                           sd.2=0.70, sd.1=0.72, 
                           n.2=82, n.1 =82) # Arabs
Effectsizes <- effect_size(eff.type="means", 
                           m.2=3.25, m.1=3.39, 
                           sd.2=0.58, sd.1=0.59, 
                           n.2=82, n.1 =82) # Jordans
### three months after attacks
Effectsizes <- effect_size(eff.type="means", 
                           m.2=3.21, m.1=3.37,
                           sd.2=0.61, sd.1=0.67, 
                           n.2=82, n.1 =82) # Palestinians
Effectsizes <- effect_size(eff.type="means", 
                           m.2=3.17, m.1=3.45, 
                           sd.2=0.71, sd.1=0.72, 
                           n.2=82, n.1 =82) # Arabs
Effectsizes <- effect_size(eff.type="means", 
                           m.2=3.14, m.1=3.39, 
                           sd.2=0.45, sd.1=0.59, 
                           n.2=82, n.1 =82) # Jordans
# c) Intellectual traits
### 1 day after attacks
Effectsizes <- effect_size(eff.type="means", 
                           m.2=2.80, m.1=2.76, 
                           sd.2=0.80, sd.1=0.81, 
                           n.2=82, n.1 =82) # Palestinians
Effectsizes <- effect_size(eff.type="means", 
                           m.2=2.89, m.1=2.88, 
                           sd.2=0.71, sd.1=0.85, 
                           n.2=82, n.1 =82) # Arabs
Effectsizes <- effect_size(eff.type="means", 
                           m.2=3.19, m.1=3.30, 
                           sd.2=0.65, sd.1=0.79, 
                           n.2=82, n.1 =82) # Jordans
### three months after attacks
Effectsizes <- effect_size(eff.type="means", 
                           m.2=2.62, m.1=2.76, 
                           sd.2=0.76, sd.1=0.81, 
                           n.2=82, n.1 =82) # Palestinians
Effectsizes <- effect_size(eff.type="means", 
                           m.2=2.80, m.1=2.88, 
                           sd.2=0.76, sd.1=0.85, 
                           n.2=82, n.1 =82) # Arabs
Effectsizes <- effect_size(eff.type="means", 
                           m.2=3.11, m.1=3.30, 
                           sd.2=0.69, sd.1=0.79, 
                           n.2=82, n.1 =82) # Jordans
# d) Social distance
### 1 day after attacks
Effectsizes <- effect_size(eff.type="means", 
                           m.2=1.72, m.1=1.96, 
                           sd.2=1.17, sd.1=1.12, 
                           n.2=78, n.1 =78) # Palestinians
Effectsizes <- effect_size(eff.type="means", 
                           m.2=2.03, m.1=2.14, 
                           sd.2=1.07, sd.1=0.99, 
                           n.2=78, n.1 =78) # Arabs
Effectsizes <- effect_size(eff.type="means", 
                           m.2=2.28, m.1=2.36, 
                           sd.2=0.92, sd.1=0.87, 
                           n.2=78, n.1 =78) # Jordans
### three months after attacks
Effectsizes <- effect_size(eff.type="means", 
                           m.2=1.73, m.1=1.96, 
                           sd.2=1.17, sd.1=1.12, 
                           n.2=78, n.1 =78) # Palestinians
Effectsizes <- effect_size(eff.type="means", 
                           m.2=1.94, m.1=2.14, 
                           sd.2=1.07, sd.1=0.99, 
                           n.2=78, n.1 =78) # Arabs
Effectsizes <- effect_size(eff.type="means", 
                           m.2=2.22, m.1=2.36, 
                           sd.2=0.95, sd.1=0.87, 
                           n.2=78, n.1 =78) # Jordans
# e) Negative feelings
### 1 day after attacks
Effectsizes <- effect_size(eff.type="means", 
                           m.2=1.89, m.1=2.01, 
                           sd.2=0.52, sd.1=0.57, 
                           n.2=80, n.1 =80) # Palestinians
Effectsizes <- effect_size(eff.type="means", 
                           m.2=1.71, m.1=1.73, 
                           sd.2=0.51, sd.1=0.54, 
                           n.2=80, n.1 =80) # Arabs
Effectsizes <- effect_size(eff.type="means", 
                           m.2=1.26, m.1=1.29, 
                           sd.2=0.35, sd.1=0.37, 
                           n.2=80, n.1 =80) # Jordans
### three months after attacks
Effectsizes <- effect_size(eff.type="means", 
                           m.2=1.89, m.1=1.82, 
                           sd.2=0.52, sd.1=0.57, 
                           n.2=80, n.1 =80) # Palestinians
Effectsizes <- effect_size(eff.type="means", 
                           m.2=1.71, m.1=1.16, 
                           sd.2=0.51, sd.1=0.54, 
                           n.2=80, n.1 =80) # Arabs
Effectsizes <- effect_size(eff.type="means", 
                           m.2=1.26, m.1=1.34, 
                           sd.2=0.35, sd.1=0.45, 
                           n.2=80, n.1 =80) # Jordans
# f) Positive feelings
### 1 day after attacks
Effectsizes <- effect_size(eff.type="means",
                           m.2=1.41, m.1=1.62, 
                           sd.2=0.40, sd.1=0.42, 
                           n.2=81, n.1 =81) # Palestinians
Effectsizes <- effect_size(eff.type="means", 
                           m.2=1.50, m.1=1.65, 
                           sd.2=0.46, sd.1=0.43, 
                           n.2=81, n.1 =81) # Arabs
Effectsizes <- effect_size(eff.type="means", 
                           m.2=1.75, m.1=1.84, 
                           sd.2=0.51, sd.1=0.47,
                           n.2=81, n.1 =81) # Jordans
### three months after attacks
Effectsizes <- effect_size(eff.type="means", 
                           m.2=1.49, m.1=1.62, 
                           sd.2=0.44, sd.1=0.42, 
                           n.2=81, n.1 =81) # Palestinians
Effectsizes <- effect_size(eff.type="means",
                           m.2=1.55, m.1=1.65, 
                           sd.2=0.50, sd.1=0.43, 
                           n.2=81, n.1 =81) # Arabs
Effectsizes <- effect_size(eff.type="means", 
                           m.2=1.66, m.1=1.84, 
                           sd.2=0.50, sd.1=0.47, 
                           n.2=81, n.1 =81) # Jordans
# g) Attribution of intention
### 1 day after attacks
Effectsizes <- effect_size(eff.type="means", 
                           m.2=3.06, m.1=3.17, 
                           sd.2=0.89, sd.1=0.79, 
                           n.2=80, n.1 =80) # Palestinians
Effectsizes <- effect_size(eff.type="means", 
                           m.2=3.24, m.1=3.20, 
                           sd.2=0.83, sd.1=1.86, 
                           n.2=80, n.1 =80) # Arabs
Effectsizes <- effect_size(eff.type="means",
                           m.2=3.99, m.1=3.91, 
                           sd.2=0.56, sd.1=0.60, 
                           n.2=80, n.1 =80) # Jordans
### three months after attacks
Effectsizes <- effect_size(eff.type="means", 
                           m.2=3.16, m.1=3.17, 
                           sd.2=0.82, sd.1=0.79, 
                           n.2=80, n.1 =80) # Palestinians
Effectsizes <- effect_size(eff.type="means", 
                           m.2=3.44, m.1=3.20, 
                           sd.2=0.78, sd.1=1.86, 
                           n.2=80, n.1 =80) # Arabs
Effectsizes <- effect_size(eff.type="means", 
                           m.2=3.93, m.1=3.91, 
                           sd.2=0.63, sd.1=0.60, 
                           n.2=80, n.1 =80) # Jordans

# Berg (2010)
Effectsizes <- effect_size(eff.type="corr", r=0.37, n=425)

# Besser & Neria (2009)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=6.78, m.2=6.63, 
                           sd.1=1.92, sd.2=2.13, 
                           n.1=160, n.2=181)

# Best, Krueger & Pearson-Merkowitz (2012)
Effectsizes <- effect_size(eff.type="corr", r=beta2cor(0.100, 1), n=915) 
Effectsizes <- effect_size(eff.type="corr", r=beta2cor(0.129, 1), n=915) 

# Bilali (2015)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=2.91, m.2=2.79, 
                           sd.1=1.18, sd.2=1.18, 
                           n.1=49, n.2=49) # general
Effectsizes <- effect_size(eff.type="t.test", t=2.33, n.1=49, n.2 =49) # low in-group homogeneity
Effectsizes <- effect_size(eff.type="t.test", t=1.89, n.1=49, n.2 =49) # high in-group homogeneity
Effectsizes <- effect_size(eff.type="corr", r=0.30, n=147)
Effectsizes <- effect_size(eff.type="corr", r=0.19, n=147)

# Blanchar (2016)
# a) Correlations
Effectsizes <- effect_size(eff.type="corr", r=0.27, n=300)
Effectsizes <- effect_size(eff.type="corr", r=0.13, n=300)
Effectsizes <- effect_size(eff.type="corr", r=0.14, n=300)
Effectsizes <- effect_size(eff.type="corr", r=0.30, n=300)
# b) Experiment
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.64, m.2=3.73, 
                           sd.1=.95, sd.2=.93, 
                           n.1=149, n.2 =151)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.83, m.2=3.76, 
                           sd.1=0.9, sd.2=.9, 
                           n.1=149, n.2 =151)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.22, m.2=3.29, 
                           sd.1=1.02, sd.2=0.98, 
                           n.1=149, n.2 =151)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.07, m.2=4.11, 
                           sd.1=1.59, sd.2=1.46, 
                           n.1=149, n.2 =151)
# c) Quasi-experiment
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.72, m.2=3.56, 
                           sd.1=0.93, sd.2=0.97, 
                           n.1=235, n.2 =65)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.78, m.2=3.88, 
                           sd.1=0.88, sd.2=0.96, 
                           n.1=235, n.2 =65)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.30, m.2=3.09, 
                           sd.1=1.00, sd.2=0.99,  
                           n.1=235, n.2 =65)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.09, m.2=4.09, 
                           sd.1=1.54, sd.2=1.49,  
                           n.1=235, n.2 =65)

# Blauwkamp et al. (2018)
Effectsizes <- effect_size(eff.type="t.test", t=4.653, n.1=389, n.2 =411)
Effectsizes <- effect_size(eff.type="t.test", t=3.247, n.1=388, n.2 =411)
Effectsizes <- effect_size(eff.type="t.test", t=4.681, n.1=395, n.2 =411)
Effectsizes <- effect_size(eff.type="t.test", t=3.713, n.1=398, n.2 =411) 
Effectsizes <- effect_size(eff.type="t.test", t=0.879, n.1=349, n.2 =411)
Effectsizes <- effect_size(eff.type="t.test", t=1.318, n.1=369, n.2 =411)

# Bonanno & Jost (2006)
Effectsizes <- effect_size(eff.type="corr", r=0.21, n=45)
Effectsizes <- effect_size(eff.type="corr", r=0.23, n=45)
Effectsizes <- effect_size(eff.type="corr", r=0.41, n=45)
Effectsizes <- effect_size(eff.type="corr", r=0.36, n=45)

# Boomgaarden & de Vreese (2007)
# a) Main effects
Effectsizes <- effect_size(eff.type="p.val", p=0.114, n.1=148, n.2=128, tail="one") # Index
Effectsizes <- effect_size(eff.type="t.test", t=2.80, n.1=148, n.2=128) # Threat index
Effectsizes <- effect_size(eff.type="prop", p1=0.286, p2=0.129, n.ab=148, n.cd=128) # Migration importance

# Boydstun et al. (2019)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=58.46181, m.2=58.85031, 
                           sd.1=25.21559, sd.2=25.87983, 
                           n.1=838, n.2=648) # After Paris
Effectsizes <- effect_size(eff.type="means", 
                           m.1=58.46181, m.2=59.58028, 
                           sd.1=25.21559, sd.2=27.12645, 
                           n.1=838, n.2=436) # After San Bernardino

# Breugelmans et al. (2009)
Effectsizes <- effect_size(eff.type="coh.d", d=0.63, n.1=333, n.2=170) # T1
Effectsizes <- effect_size(eff.type="coh.d", d=-0.22, n.1=333, n.2=306) # T2
Effectsizes <- effect_size(eff.type="coh.d", d=-0.58, n.1=333, n.2=462) # T3

# Brinson & Stohl (2012)
Effectsizes <- effect_size(eff.type="corr", r=sqrt(0.19-0.01), n=371) # Domestic frame
Effectsizes <- effect_size(eff.type="corr", r=sqrt(0.24-0.09), n=371) # International frame
Effectsizes <- effect_size(eff.type="corr", r=sqrt(0.22-0.02), n=371) # Domestic frame (others)
Effectsizes <- effect_size(eff.type="corr", r=sqrt(0.26-0.08), n=371) # International frame (others)
Effectsizes <- effect_size(eff.type="corr", r=sqrt(0.11-0.01), n=371) # Domestic frame (self)
Effectsizes <- effect_size(eff.type="corr", r=sqrt(0.12-0.03), n=371) # International frame (self)

# Bruneau et al. (2019)
Effectsizes <- effect_size(eff.type="corr", r=0.71, n=193) 
Effectsizes <- effect_size(eff.type="corr", r=0.73, n=193) 
Effectsizes <- effect_size(eff.type="corr", r=0.68, n=193) 
Effectsizes <- effect_size(eff.type="corr", r=0.72, n=193) 
Effectsizes <- effect_size(eff.type="corr", r=0.56, n=193) 

# Canetti et al. (2009)
Effectsizes <- effect_size(eff.type="corr", r=0.34, n=190) 
Effectsizes <- effect_size(eff.type="corr", r=0.64, n=190) 

# Canetti et al. (2009)
Effectsizes <- effect_size(eff.type="corr", r=0.05, n=469) 
Effectsizes <- effect_size(eff.type="corr", r=0.09, n=469) 
Effectsizes <- effect_size(eff.type="corr", r=0.18, n=469) 

# Canetti et al. (2016)
# Study 2
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.20, m.2=4.06, 
                           sd.1=0.99, sd.2=0.98, 
                           n.1=137, n.2=93)

# Canetti et al. (2017)
# a) Support for Compromise
Effectsizes <- effect_size(eff.type="corr", r=0.04, n=781) # Israeli
Effectsizes <- effect_size(eff.type="corr", r=-0.01, n=1196) # Jews
Effectsizes <- effect_size(eff.type="corr", r=-0.01, n=781) # Israeli
Effectsizes <- effect_size(eff.type="corr", r=-0.04, n=1196) # Jews
Effectsizes <- effect_size(eff.type="corr", r=0.08, n=781) # Israeli
Effectsizes <- effect_size(eff.type="corr", r=0.06, n=1196) # Jews
# b) Political orientation
Effectsizes <- effect_size(eff.type="corr", r=0.10, n=781) # Israeli
Effectsizes <- effect_size(eff.type="corr", r=-0.03, n=1196) # Jews
Effectsizes <- effect_size(eff.type="corr", r=0.11, n=781) # Israeli
Effectsizes <- effect_size(eff.type="corr", r=-0.09, n=1196) # Jews
Effectsizes <- effect_size(eff.type="corr", r=0.06, n=781) # Israeli
Effectsizes <- effect_size(eff.type="corr", r=-0.00, n=1196) # Jews

# Canetti et al. (2018)
# a) Study 1
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.54, m.2=3.83,
                           sd.1=1.04, sd.2=0.99, 
                           n.1=19, n.2=19)
# b) Study 2
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.87, m.2=4.22, 
                           sd.1=0.93, sd.2=0.88, 
                           n.1=23, n.2=23)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.07, m.2=4.22, 
                           sd.1=1, sd.2=0.88, 
                           n.1=23, n.2=23)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.67, m.2=3.15, 
                           sd.1=2.15, sd.2=2.03, 
                           n.1=23, n.2=23)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.67, m.2=5.17, 
                           sd.1=2.15, sd.2=1.99, 
                           n.1=23, n.2=23)
# c) Study 3
Effectsizes <- effect_size(eff.type="t.test", t=1.84, n.1=105, n.2 =47) 
# d) Study 4
# Holocaust group vs. non-Holocaust group
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.99, m.2=4.08, 
                           sd.1=1.23, sd.2=1.13, 
                           n.1=323, n.2=539) #militancy
Effectsizes <- effect_size(eff.type="means", 
                           m.1=2.85, m.2=3.08, 
                           sd.1=1.76, sd.2=1.80, 
                           n.1=512, n.2=312) #compromise
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.29, m.2=3.63, 
                           sd.1=1.26, sd.2=1.41, 
                           n.1=474, n.2=296) #pol. orientation
# bivariate correlations within non-H and H group: self-reported exposure
Effectsizes <- effect_size(eff.type="corr", r=0.01, n=539) #militancy N-H
Effectsizes <- effect_size(eff.type="corr", r=0.08, n=323) #militancy H
Effectsizes <- effect_size(eff.type="corr", r=0.05, n=512) #compromise N-H
Effectsizes <- effect_size(eff.type="corr", r=0.091, n=312) #compromise H
Effectsizes <- effect_size(eff.type="corr", r=0.08, n=474) #pol. orientation N-H
Effectsizes <- effect_size(eff.type="corr", r=0.185, n=296) #pol. orientation H
# bivariate correlations within non-H and H group: perceived threat
Effectsizes <- effect_size(eff.type="corr", r=0.081, n=539) #militancy N-H
Effectsizes <- effect_size(eff.type="corr", r=0.154, n=323) #militancy H
Effectsizes <- effect_size(eff.type="corr", r=0.090, n=512) #compromise N-H
Effectsizes <- effect_size(eff.type="corr", r=0.222, n=312) #compromise H
Effectsizes <- effect_size(eff.type="corr", r=0.140, n=474) #pol. orientation N-H
Effectsizes <- effect_size(eff.type="corr", r=0.127, n=296) #pol. orientation H

# Canetti et al. (2018)
# Full sample
Effectsizes <- effect_size(eff.type="means", 
                           m.1=0.54, m.2=0.49, 
                           sd.1=0.28, sd.2=0.30, 
                           n.1=177, n.2=243)

# Canetti et al. (2019)
Effectsizes <- effect_size(eff.type="corr", r=0.23, n=83)
Effectsizes <- effect_size(eff.type="corr", r=0.24, n=41)

# Canetti et al. (2008)
Effectsizes <- effect_size(eff.type="corr", r=0.51, n=504) 
Effectsizes <- effect_size(eff.type="corr", r=0.60, n=504)
Effectsizes <- effect_size(eff.type="corr", r=0.34, n=504) 
Effectsizes <- effect_size(eff.type="corr", r=0.47, n=504) 

# Carnagey & Anderson (2007)
Effectsizes <- effect_size(eff.type="coh.d", d=0.37, n.1=780, n.2=1034) # T1: war
Effectsizes <- effect_size(eff.type="coh.d", d=0.14, n.1=780, n.2=1034) # T1: prison
Effectsizes <- effect_size(eff.type="coh.d", d=0.24, n.1=780, n.2=889) # T2: war
Effectsizes <- effect_size(eff.type="coh.d", d=0.07, n.1=780, n.2=889) # T2: prison

# Carriere, Hendricks, & Moghaddam (2019)
Effectsizes <- effect_size(eff.type="corr", r=0.080, n=2019)
Effectsizes <- effect_size(eff.type="corr", r=0.020, n=2037)
Effectsizes <- effect_size(eff.type="corr", r=0.100, n=1896)
Effectsizes <- effect_size(eff.type="corr", r=0.077, n=1925)
Effectsizes <- effect_size(eff.type="corr", r=0.100, n=4508)
Effectsizes <- effect_size(eff.type="corr", r=0.183, n=5283)
Effectsizes <- effect_size(eff.type="corr", r=0.107, n=5320)
Effectsizes <- effect_size(eff.type="corr", r=0.159, n=5072)
Effectsizes <- effect_size(eff.type="corr", r=0.168, n=5246)
Effectsizes <- effect_size(eff.type="corr", r=0.133, n=3606)
Effectsizes <- effect_size(eff.type="corr", r=0.093, n=3590)
Effectsizes <- effect_size(eff.type="corr", r=0.064, n=3550)
Effectsizes <- effect_size(eff.type="corr", r=0.121, n=3546)

# Castanho Silva (2018)
# a) Study 1
Effectsizes <- effect_size(eff.type="means", 
                           m.1=5.19, m.2=5.07, 
                           sd.1=2.14, sd.2= 2.38, 
                           n.1=267, n.2=430)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.86, m.2=4.99, 
                           sd.1=2.01, sd.2=2.10, 
                           n.1=2852, n.2=2033)
# b) Study 2
Effectsizes <- effect_size(eff.type="means", 
                           m.1=2.55, m.2=2.74, 
                           sd.1=1.00, sd.2=0.95,  
                           n.1=97, n.2=928)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=2.88, m.2=2.79, 
                           sd.1=0.99, sd.2=0.97, 
                           n.1=7537, n.2=19058)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=2.31, m.2=2.34, 
                           sd.1=0.96, sd.2=0.94, 
                           n.1=97, n.2=928)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=2.43, m.2=2.19, 
                           sd.1=1.02, sd.2=0.95, 
                           n.1=7537, n.2=19058)

# Celebi et al. (2009)
# a) Outgroup trust
Effectsizes <- effect_size(eff.type="corr", r=-0.24, n=288) # Kurds
Effectsizes <- effect_size(eff.type="corr", r=0.10, n=337) # Turks
# b) National identification
Effectsizes <- effect_size(eff.type="corr", r=0.23, n=288) # Kurds
Effectsizes <- effect_size(eff.type="corr", r=0.11, n=337) # Turks

# Cheung-Blunden (2020)
# a) US sample
Effectsizes <- effect_size(eff.type="corr", r=0.18, n=220)
Effectsizes <- effect_size(eff.type="corr", r=0.27, n=220)
Effectsizes <- effect_size(eff.type="corr", r=0.03, n=220)
Effectsizes <- effect_size(eff.type="corr", r=0.21, n=220)
Effectsizes <- effect_size(eff.type="corr", r=0.12, n=220)
Effectsizes <- effect_size(eff.type="corr", r=-0.02, n=220)
# b) German sample
Effectsizes <- effect_size(eff.type="corr", r=0.21, n=230)
Effectsizes <- effect_size(eff.type="corr", r=0.22, n=230)
Effectsizes <- effect_size(eff.type="corr", r=0.21, n=230)
Effectsizes <- effect_size(eff.type="corr", r=0.15, n=230)
Effectsizes <- effect_size(eff.type="corr", r=0.31, n=230)
Effectsizes <- effect_size(eff.type="corr", r=-0.17, n=230)
# c) British sample
Effectsizes <- effect_size(eff.type="corr", r=0.10, n=151)
Effectsizes <- effect_size(eff.type="corr", r=0.07, n=151)
Effectsizes <- effect_size(eff.type="corr", r=0.31, n=151)
Effectsizes <- effect_size(eff.type="corr", r=0.03, n=151)
Effectsizes <- effect_size(eff.type="corr", r=0.11, n=151)
Effectsizes <- effect_size(eff.type="corr", r=0.18, n=151)
Effectsizes <- effect_size(eff.type="corr", r=0.02, n=151)
# d) Spanish sample
Effectsizes <- effect_size(eff.type="corr", r=-0.01, n=183)
Effectsizes <- effect_size(eff.type="corr", r=-0.01, n=183)
Effectsizes <- effect_size(eff.type="corr", r=0.14, n=183)
Effectsizes <- effect_size(eff.type="corr", r=0.26, n=183)
Effectsizes <- effect_size(eff.type="corr", r=0.07, n=183)
Effectsizes <- effect_size(eff.type="corr", r=-0.03, n=183)
Effectsizes <- effect_size(eff.type="corr", r=0.00, n=183)

# Cheung-Blunden & Blunden (2010)
# a) Anger
Effectsizes <- effect_size(eff.type="corr", r=0.45, n=588) 
Effectsizes <- effect_size(eff.type="corr", r=0.37, n=588) 
Effectsizes <- effect_size(eff.type="corr", r=0.19, n=588) 
# b) Relevance scale
Effectsizes <- effect_size(eff.type="corr", r=0.16, n=588) 
Effectsizes <- effect_size(eff.type="corr", r=0.04, n=588) 
Effectsizes <- effect_size(eff.type="corr", r=0.13, n=588) 

# Choma et al. (2015)
Effectsizes <- effect_size(eff.type="t.test", t=1.68, n.1=58, n.2=(57+59)) 
Effectsizes <- effect_size(eff.type="t.test", t=2.81, n.1=58, n.2=(57+59)) 

# Choma et al. (2018)
# a) Study 1
Effectsizes <- effect_size(eff.type="corr", r=0.07, n=215) 
Effectsizes <- effect_size(eff.type="corr", r=0.04, n=215) 
Effectsizes <- effect_size(eff.type="corr", r=-0.05, n=215) 
Effectsizes <- effect_size(eff.type="corr", r=0.15, n=215) 
Effectsizes <- effect_size(eff.type="corr", r=0.11, n=215) 
# b) Study 2
Effectsizes <- effect_size(eff.type="corr", r=-0.01, n=468) 
Effectsizes <- effect_size(eff.type="corr", r=0.07, n=468) 
Effectsizes <- effect_size(eff.type="corr", r=0.03, n=468) 
Effectsizes <- effect_size(eff.type="corr", r=0.06, n=468) 
Effectsizes <- effect_size(eff.type="coh.d", d=0.27, n.1=301, n.2=176)
Effectsizes <- effect_size(eff.type="coh.d", d=0.28, n.1=301, n.2=176)
Effectsizes <- effect_size(eff.type="coh.d", d=0.22, n.1=301, n.2=176)
Effectsizes <- effect_size(eff.type="coh.d", d=0.27, n.1=301, n.2=176)

# Cohrs et al. (2005a)
Effectsizes <- effect_size(eff.type="corr", r=0.16, n=275) 
Effectsizes <- effect_size(eff.type="corr", r=0.15, n=275) 
Effectsizes <- effect_size(eff.type="corr", r=0.17, n=275) 
Effectsizes <- effect_size(eff.type="corr", r=0.07, n=275) 
Effectsizes <- effect_size(eff.type="corr", r=0.13, n=275) 

# Cohrs et al. (2005b)
Effectsizes <- effect_size(eff.type="corr", r=0.117, n=547) 
Effectsizes <- effect_size(eff.type="corr", r=0.149, n=548) 
Effectsizes <- effect_size(eff.type="corr", r=0.209, n=549) 

# Cohu et al. (2016)
# a) Prejudice
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.32, m.2=3.49, 
                           sd.1=1.12, sd.2=1.26, 
                           n.1=56, n.2=65) # T1
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.12, m.2=3.49, 
                           sd.1=1.37, sd.2=1.26, 
                           n.1=41, n.2=65) # T2
# b) SDO
Effectsizes <- effect_size(eff.type="means", 
                           m.1=2.54, m.2=2.10, 
                           sd.1=1.33, sd.2=1.03, 
                           n.1=56, n.2=65) # T1
Effectsizes <- effect_size(eff.type="means", 
                           m.1=2.23, m.2=2.10, 
                           sd.1=0.82, sd.2=1.03, 
                           n.1=41, n.2=65) # T2
# c) Laicite
Effectsizes <- effect_size(eff.type="means", 
                           m.1=5.25, m.2=5.44, 
                           sd.1=1.80, sd.2=1.56, 
                           n.1=56, n.2=65) # T1
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.44, m.2=5.44, 
                           sd.1=1.98, sd.2=1.56, 
                           n.1=41, n.2=65) # T2

# Conejero & Etxebarria (2007)
Effectsizes <- effect_size(eff.type="corr", r=0.34, n=1346)
Effectsizes <- effect_size(eff.type="corr", r=-0.14, n=232)
Effectsizes <- effect_size(eff.type="corr", r=sqrt(0.03), n=1196)

# Conrad (2018): sample descriptives calculated via replication data
Effectsizes <- effect_size(eff.type="t.test", t=-3.47, n.1=292, n.2=304) 


# Coryn, Baele, & Myers (2004)
Effectsizes <- effect_size(eff.type="corr", r=0.428, n=301)
Effectsizes <- effect_size(eff.type="corr", r=0.290, n=301)
Effectsizes <- effect_size(eff.type="corr", r=0.310, n=301)

# Crowson (2009)
Effectsizes <- effect_size(eff.type="corr", r=0.557, n=176)
Effectsizes <- effect_size(eff.type="corr", r=0.377, n=176)
Effectsizes <- effect_size(eff.type="corr", r=0.295, n=176)
Effectsizes <- effect_size(eff.type="corr", r=0.241, n=176)
Effectsizes <- effect_size(eff.type="corr", r=0.422, n=176)
Effectsizes <- effect_size(eff.type="corr", r=0.286, n=176)

# Crowson et al. (2005)
Effectsizes <- effect_size(eff.type="corr", r=0.21, n=136) 
Effectsizes <- effect_size(eff.type="corr", r=0.01, n=136) 
Effectsizes <- effect_size(eff.type="corr", r=0.23, n=136) 
Effectsizes <- effect_size(eff.type="corr", r=0.20, n=136) 
Effectsizes <- effect_size(eff.type="corr", r=0.20, n=136) 
Effectsizes <- effect_size(eff.type="corr", r=0.06, n=136) 

# Da Costa Silva (2019)
Effectsizes <- effect_size(eff.type="t.test", t=2.49, n.1=(196/2), n.2 =(196/2))
Effectsizes <- effect_size(eff.type="t.test", t=2.06, n.1=(300/2), n.2 =(300/2))

# Das et al. (2009)
# Study 1
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.92, m.2=5.13, 
                           sd.1=1.77, sd.2=1.46,
                           n.1=28, n.2=16) #before murder
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.92, m.2=4.65,
                           sd.1=2.23, sd.2=2.23,
                           n.1=22, n.2=34) #after murder
# Study 3
Effectsizes <- effect_size(eff.type="means", 
                           m.1=-0.12, m.2=-0.43, 
                           sd.1=0.85, sd.2=0.65,
                           n.1=28, n.2=35) # Muslim sample
Effectsizes <- effect_size(eff.type="means", 
                           m.1=0.37, m.2=0.2, 
                           sd.1=0.61, sd.2=0.68,
                           n.1=39, n.2=38) # Non-muslim sample
# Study 2
Effectsizes <- effect_size(eff.type="corr", r=0.26, n=96)

# Davies, Steele, & Markus
Effectsizes <- effect_size(eff.type="f.test", f=7.06, n.1=60, n.2=60)

# Davis (2006) and Davis & Silver (2004)
Effectsizes <- effect_size(eff.type="corr", r=0.22, n=514+653+164+48)
Effectsizes <- effect_size(eff.type="corr", r=0.04, n=1081)

# Dawson (2011)
# a) Fear
Effectsizes <- effect_size(eff.type="means", 
                           m.1=2.1146, m.2=2.0980, 
                           sd.1=1.28129, sd.2=1.31336, 
                           n.1=96, n.2=68) # Muslims
Effectsizes <- effect_size(eff.type="means", 
                           m.1=1.9282, m.2=2.0343, 
                           sd.1=1.27107, sd.2=1.42838, 
                           n.1=96, n.2=68) # Immigrants
# b) Dehumanization
Effectsizes <- effect_size(eff.type="means", 
                           m.1=0.8984, m.2=0.9449, 
                           sd.1=0.9926, sd.2=0.73409, 
                           n.1=96, n.2=68) # Muslims
Effectsizes <- effect_size(eff.type="means", 
                           m.1=1.3464, m.2=1.3078, 
                           sd.1=1.08298, sd.2=.97991, 
                           n.1=96, n.2=68) # Immigrants
# c) Moral Depravity
Effectsizes <- effect_size(eff.type="means", 
                           m.1=1.6423, m.2=1.5356, 
                           sd.1=.94057, sd.2=.92845, 
                           n.1=96, n.2=68) # Muslims
Effectsizes <- effect_size(eff.type="means", 
                           m.1=2.0843, m.2=2.1471, 
                           sd.1=.89698, sd.2=.94950, 
                           n.1=96, n.2=68) # Immigrants

# De Coninck et al. (2018)
Effectsizes <- effect_size(eff.type="corr", r=0.27, n=1500) # Immigrants
Effectsizes <- effect_size(eff.type="corr", r=0.30, n=1500) # Refugees

# Dinesen & Jaeger (2013)
# a) Post 0-7 days
Effectsizes <- effect_size(eff.type="prop", p1=0.58, p2=0.43, n.ab=72, n.cd=928)
Effectsizes <- effect_size(eff.type="prop", p1=0.59, p2=0.44, n.ab=72, n.cd=928)
Effectsizes <- effect_size(eff.type="prop", p1=0.50, p2=0.27, n.ab=72, n.cd=928)
Effectsizes <- effect_size(eff.type="prop", p1=0.62, p2=0.48, n.ab=72, n.cd=928)
Effectsizes <- effect_size(eff.type="prop", p1=0.73, p2=0.61, n.ab=72, n.cd=928)
# b) Post 7 months
Effectsizes <- effect_size(eff.type="prop", p1=0.55, p2=0.43, n.ab=1023, n.cd=928)
Effectsizes <- effect_size(eff.type="prop", p1=0.53, p2=0.44, n.ab=1023, n.cd=928)
Effectsizes <- effect_size(eff.type="prop", p1=0.29, p2=0.27, n.ab=1023, n.cd=928)
Effectsizes <- effect_size(eff.type="prop", p1=0.48, p2=0.48, n.ab=1023, n.cd=928)
Effectsizes <- effect_size(eff.type="prop", p1=0.64, p2=0.61, n.ab=1023, n.cd=928)
# c) Post 14 months
Effectsizes <- effect_size(eff.type="prop", p1=0.43, p2=0.43, n.ab=1024, n.cd=928)
Effectsizes <- effect_size(eff.type="prop", p1=0.42, p2=0.44, n.ab=1024, n.cd=928)
Effectsizes <- effect_size(eff.type="prop", p1=0.26, p2=0.27, n.ab=1024, n.cd=928)
Effectsizes <- effect_size(eff.type="prop", p1=0.50, p2=0.48, n.ab=1024, n.cd=928)
# d) Post 19 months
Effectsizes <- effect_size(eff.type="prop", p1=0.47, p2=0.43, n.ab=1015, n.cd=928)
Effectsizes <- effect_size(eff.type="prop", p1=0.47, p2=0.44, n.ab=1015, n.cd=928)
Effectsizes <- effect_size(eff.type="prop", p1=0.31, p2=0.27, n.ab=1015, n.cd=928)
Effectsizes <- effect_size(eff.type="prop", p1=0.50, p2=0.48, n.ab=1015, n.cd=928)
Effectsizes <- effect_size(eff.type="prop", p1=0.62, p2=0.61, n.ab=1015, n.cd=928)

# Doosje et al. (2009)
Effectsizes <- effect_size(eff.type="corr", r=0.13, n=690)
Effectsizes <- effect_size(eff.type="corr", r=0.22, n=690)  
Effectsizes <- effect_size(eff.type="corr", r=0.06, n=690) 
Effectsizes <- effect_size(eff.type="corr", r=0.39, n=669) 
Effectsizes <- effect_size(eff.type="corr", r=0.16, n=679) 
Effectsizes <- effect_size(eff.type="corr", r=0.20, n=673) 
Effectsizes <- effect_size(eff.type="corr", r=0.17, n=690) 

# Doosje, van der Veen, & Klaver (2018)
# a) Political trust
Effectsizes <- effect_size(eff.type="means", m.1=4.42, m.2=4.49, sd.1=2.03, sd.2=2.12, n.1=1768, n.2=1859) #T1 (Belgium)
Effectsizes <- effect_size(eff.type="means", m.1=4.62, m.2=4.49, sd.1=1.95, sd.2=2.12, n.1=1791, n.2=1859) #T2 (Belgium)
Effectsizes <- effect_size(eff.type="means", m.1=3.01, m.2=2.97, sd.1=2.08, sd.2=2.09, n.1=2000, n.2=2927) #T1 (Czech Republic)
Effectsizes <- effect_size(eff.type="means", m.1=2.95, m.2=2.97, sd.1=2.14, sd.2=2.09, n.1=2348, n.2=2927) #T2 (Czech Republic)
Effectsizes <- effect_size(eff.type="means", m.1=5.74, m.2=5.91, sd.1=1.75, sd.2=1.94, n.1=41, n.2=3050) #T1 (Denmark)
Effectsizes <- effect_size(eff.type="means", m.1=5.37, m.2=5.91, sd.1=1.97, sd.2=1.94, n.1=1553, n.2=3050) #T2 (Denmark)
Effectsizes <- effect_size(eff.type="means", m.1=5.46, m.2=5.44, sd.1=1.93, sd.2=1.99, n.1=1770, n.2=2300) #T1 (Finland)
Effectsizes <- effect_size(eff.type="means", m.1=4.88, m.2=5.44, sd.1=2.09, sd.2=1.99, n.1=1866, n.2=2300) #T2 (Finland)
Effectsizes <- effect_size(eff.type="means", m.1=4.06, m.2=3.85, sd.1=1.97, sd.2=2.09, n.1=1061, n.2=4532) #T1 (Germany)
Effectsizes <- effect_size(eff.type="means", m.1=3.81, m.2=3.85, sd.1=2.11, sd.2=2.09, n.1=2998, n.2=4532) #T2 (Germany)
Effectsizes <- effect_size(eff.type="means", m.1=2.21, m.2=2.92, sd.1=2.01, sd.2=2.09, n.1=1502, n.2=1562) #T1 (Hungary)
Effectsizes <- effect_size(eff.type="means", m.1=3.66, m.2=2.92, sd.1=2.32, sd.2=2.09, n.1=1533, n.2=1562) #T2 (Hungary)
Effectsizes <- effect_size(eff.type="means", m.1=4.43, m.2=4.26, sd.1=2.13, sd.2=2.16, n.1=690, n.2=3342) #T1 (Ireland)
Effectsizes <- effect_size(eff.type="means", m.1=3.53, m.2=4.26, sd.1=2.14, sd.2=2.16, n.1=1282, n.2=3342) #T2 (Ireland)
Effectsizes <- effect_size(eff.type="means", m.1=5.47, m.2=5.06, sd.1=1.91, sd.2=2.16, n.1=1544, n.2=1541) #T1 (Norway)
Effectsizes <- effect_size(eff.type="means", m.1=5.66, m.2=5.06, sd.1=1.90, sd.2=2.16, n.1=1610, n.2=1541) #T2 (Norway)
Effectsizes <- effect_size(eff.type="means", m.1=4.29, m.2=4.01, sd.1=2.04, sd.2=2.15, n.1=1621, n.2=1646) #T1 (Spain)
Effectsizes <- effect_size(eff.type="means", m.1=4.22, m.2=4.01, sd.1=2.09, sd.2=2.15, n.1=1848, n.2=1646) #T2 (Spain)
Effectsizes <- effect_size(eff.type="means", m.1=5.73, m.2=5.32, sd.1=1.78, sd.2=1.94, n.1=547, n.2=2748) #T1 (Sweden)
Effectsizes <- effect_size(eff.type="means", m.1=5.37, m.2=5.32, sd.1=2.02, sd.2=1.94, n.1=1829, n.2=2748) #T2 (Sweden)
Effectsizes <- effect_size(eff.type="means", m.1=5.36, m.2=5.34, sd.1=1.78, sd.2=1.81, n.1=1774, n.2=1775) #T1 (Switzerland)
Effectsizes <- effect_size(eff.type="means", m.1=5.39, m.2=5.34, sd.1=1.95, sd.2=1.81, n.1=1480, n.2=1775) #T2 (Switzerland)
Effectsizes <- effect_size(eff.type="means", m.1=4.97, m.2=5.00, sd.1=1.78, sd.2=1.84, n.1=2622, n.2=2493) #T1 (The Netherlands)
Effectsizes <- effect_size(eff.type="means", m.1=5.23, m.2=5.00, sd.1=1.69, sd.2=1.84, n.1=989, n.2=2493) #T2 (The Netherlands)
Effectsizes <- effect_size(eff.type="means", m.1=3.83, m.2=3.96, sd.1=2.09, sd.2=2.14, n.1=2936, n.2=1884) #T1 (The UK)
Effectsizes <- effect_size(eff.type="means", m.1=3.97, m.2=3.96, sd.1=2.21, sd.2=2.14, n.1=2334, n.2=1884) #T2 (The UK)
# b) Institutional trust
Effectsizes <- effect_size(eff.type="means", m.1=5.29, m.2=4.97, sd.1=2.05, sd.2=2.15, n.1=1768, n.2=1859) #T1 (Belgium)
Effectsizes <- effect_size(eff.type="means", m.1=5.37, m.2=2.04, sd.1=1.95, sd.2=2.15, n.1=1791, n.2=1859) #T2 (Belgium)
Effectsizes <- effect_size(eff.type="means", m.1=4.51, m.2=3.98, sd.1=2.24, sd.2=2.26, n.1=2000, n.2=2927) #T1 (Czech Republic)
Effectsizes <- effect_size(eff.type="means", m.1=4.52, m.2=3.98, sd.1=2.20, sd.2=2.26, n.1=2348, n.2=2927) #T2 (Czech Republic)
Effectsizes <- effect_size(eff.type="means", m.1=6.33, m.2=7.36, sd.1=2.52, sd.2=1.80, n.1=41, n.2=3050) #T1 (Denmark)
Effectsizes <- effect_size(eff.type="means", m.1=7.45, m.2=7.36, sd.1=1.74, sd.2=1.80, n.1=1553, n.2=3050) #T2 (Denmark)
Effectsizes <- effect_size(eff.type="means", m.1=7.55, m.2=7.53, sd.1=1.67, sd.2=1.63, n.1=1770, n.2=2300) #T1 (Finland)
Effectsizes <- effect_size(eff.type="means", m.1=7.44, m.2=7.53, sd.1=1.67, sd.2=1.63, n.1=1866, n.2=2300) #T2 (Finland)
Effectsizes <- effect_size(eff.type="means", m.1=6.40, m.2=6.15, sd.1=2.00, sd.2=2.09, n.1=1061, n.2=4532) #T1 (Germany)
Effectsizes <- effect_size(eff.type="means", m.1=6.28, m.2=6.15, sd.1=2.09, sd.2=2.09, n.1=2998, n.2=4532) #T2 (Germany)
Effectsizes <- effect_size(eff.type="means", m.1=3.94, m.2=4.71, sd.1=2.36, sd.2=2.49, n.1=1502, n.2=1562) #T1 (Hungary)
Effectsizes <- effect_size(eff.type="means", m.1=4.87, m.2=4.71, sd.1=2.30, sd.2=2.49, n.1=1533, n.2=1562) #T2 (Hungary)
Effectsizes <- effect_size(eff.type="means", m.1=5.61, m.2=5.78, sd.1=2.22, sd.2=2.13, n.1=690, n.2=3342) #T1 (Ireland)
Effectsizes <- effect_size(eff.type="means", m.1=5.78, m.2=5.78, sd.1=1.99, sd.2=2.13, n.1=1282, n.2=3342) #T2 (Ireland)
Effectsizes <- effect_size(eff.type="means", m.1=7.01, m.2=6.79, sd.1=1.88, sd.2=1.89, n.1=1544, n.2=1541) #T1 (Norway)
Effectsizes <- effect_size(eff.type="means", m.1=7.15, m.2=6.79, sd.1=1.75, sd.2=1.89, n.1=1610, n.2=1541) #T2 (Norway)
Effectsizes <- effect_size(eff.type="means", m.1=5.26, m.2=4.83, sd.1=2.06, sd.2=2.22, n.1=1621, n.2=1646) #T1 (Spain)
Effectsizes <- effect_size(eff.type="means", m.1=5.55, m.2=4.83, sd.1=2.03, sd.2=2.22, n.1=1848, n.2=1646) #T2 (Spain)
Effectsizes <- effect_size(eff.type="means", m.1=6.74, m.2=6.45, sd.1=1.85, sd.2=1.91, n.1=547, n.2=2748) #T1 (Sweden)
Effectsizes <- effect_size(eff.type="means", m.1=6.53, m.2=6.45, sd.1=1.93, sd.2=1.91, n.1=1829, n.2=2748) #T2 (Sweden)
Effectsizes <- effect_size(eff.type="means", m.1=6.59, m.2=6.61, sd.1=1.84, sd.2=1.89, n.1=1774, n.2=1775) #T1 (Switzerland)
Effectsizes <- effect_size(eff.type="means", m.1=6.67, m.2=6.61, sd.1=1.96, sd.2=1.89, n.1=1480, n.2=1775) #T2 (Switzerland)
Effectsizes <- effect_size(eff.type="means", m.1=5.80, m.2=5.70, sd.1=1.76, sd.2=1.84, n.1=2622, n.2=2493) #T1 (The Netherlands)
Effectsizes <- effect_size(eff.type="means", m.1=5.94, m.2=5.70, sd.1=1.68, sd.2=1.84, n.1=989, n.2=2493) #T2 (The Netherlands)
Effectsizes <- effect_size(eff.type="means", m.1=5.51, m.2=5.63, sd.1=2.10, sd.2=2.08, n.1=2936, n.2=1884) #T1 (The UK)
Effectsizes <- effect_size(eff.type="means", m.1=5.75, m.2=5.63, sd.1=2.14, sd.2=2.08, n.1=2334, n.2=1884) #T2 (The UK)

# Dvir-Gvirsman et al. (2016)
# Palestinian sample
Effectsizes <- effect_size(eff.type="corr", r=-0.158, n=569)
Effectsizes <- effect_size(eff.type="corr", r=0.133, n=569)
Effectsizes <- effect_size(eff.type="corr", r=0.194, n=572)
Effectsizes <- effect_size(eff.type="corr", r=0.141, n=572)
# Palestinian sample
Effectsizes <- effect_size(eff.type="corr", r=0.093, n=385)
Effectsizes <- effect_size(eff.type="corr", r=0.152, n=385)
Effectsizes <- effect_size(eff.type="corr", r=0.075, n=385)
Effectsizes <- effect_size(eff.type="corr", r=0.033, n=385)
# Palestinian sample
Effectsizes <- effect_size(eff.type="corr", r=0.082, n=282)
Effectsizes <- effect_size(eff.type="corr", r=0.057, n=282)
Effectsizes <- effect_size(eff.type="corr", r=0.124, n=282)
Effectsizes <- effect_size(eff.type="corr", r=0.065, n=282)

# Eadeh et al. (2020)
# a) Experiment 1A
Effectsizes <- effect_size(eff.type="f.test", f=8.04, n.1=142, n.2=142) # Republican politician
Effectsizes <- effect_size(eff.type="f.test", f=4.49, n.1=142, n.2=142) # Democratic politician
# b) Experiment 1B
Effectsizes <- effect_size(eff.type="f.test", f=3.90, n.1=110, n.2=110) 
# c) Experiment 2
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.15, m.2=3.95,
                           sd.1=1.36, sd.2=1.51, 
                           n.1=88, n.2=88)
# d) Experiment 3
# anger manipulation
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.11, m.2=4.22, 
                           sd.1=1.29, sd.2=1.27, 
                           n.1=104, n.2=104)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.78, m.2=3.94, 
                           sd.1=1.43, sd.2=1.36, 
                           n.1=104, n.2=104)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=2.94, m.2=3.01, 
                           sd.1=1.13, sd.2=1.09, 
                           n.1=104, n.2=104)
# fear manipulation
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.14, m.2=4.22, 
                           sd.1=1.32, sd.2=1.27, 
                           n.1=104, n.2=104)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.09, m.2=3.94, 
                           sd.1=1.46, sd.2=1.36, 
                           n.1=104, n.2=104)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.16, m.2=3.01, 
                           sd.1=1.15, sd.2=1.09, 
                           n.1=104, n.2=104)
# correlation with risk of terrorism
Effectsizes <- effect_size(eff.type="corr", r=0.22, n=312)
Effectsizes <- effect_size(eff.type="corr", r=0.06, n=312)
Effectsizes <- effect_size(eff.type="corr", r=0.02, n=312)
# e) Experiment 4
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.68, m.2=4.07, 
                           sd.1=1.60, sd.2=1.42, 
                           n.1=91, n.2=91)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.47, m.2=4.07, 
                           sd.1=1.56, sd.2=1.42, 
                           n.1=91, n.2=91)

# Eadeh, Godefroidt, & Trojan (2020)
# Eadehetal2020 <- read_sav("https://www.dropbox.com/s/px4scj614s6kmgb/Eadeh%2C%20Godefroidt%2C%20%26%20Trojan%20%282020%29%20Unpublished%20dataset.sav?dl=1")
# to replicate results
# describeBy(Eadehetal2020$rwa, group=Eadehetal2020$expcondition)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.542071, m.2=3.608974, 
                           sd.1=1.080305, sd.2=1.185532, 
                           n.1=104, n.2=103)
# describeBy(Eadehetal2020$sdo, group=Eadehetal2020$expcondition)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=1.567961, m.2=1.382212, 
                           sd.1=1.858735, sd.2=1.963859, 
                           n.1=103, n.2=104)
# describeBy(Eadehetal2020$hawkishness, group=Eadehetal2020$expcondition)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.532362, m.2=3.479167, 
                           sd.1=1.136888, sd.2=1.149243, 
                           n.1=103, n.2=104)
# describeBy(Eadehetal2020$hawkishness_2, group=Eadehetal2020$expcondition)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.171521, m.2=4.035256, 
                           sd.1=1.474448, sd.2=1.624328, 
                           n.1=103, n.2=104)
# describeBy(Eadehetal2020$dovish, group=Eadehetal2020$expcondition)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=5.076923, m.2=5.106796, 
                           sd.1=1.138927, sd.2=1.138927, 
                           n.1=104, n.2=103)

# Earle & Hodson (2019)
Effectsizes <- effect_size(eff.type="corr", r=0.61, n=310)
Effectsizes <- effect_size(eff.type="corr", r=0.46, n=310)
Effectsizes <- effect_size(eff.type="corr", r=0.32, n=310)
Effectsizes <- effect_size(eff.type="corr", r=0.53, n=310)
Effectsizes <- effect_size(eff.type="corr", r=0.57, n=310)
Effectsizes <- effect_size(eff.type="corr", r=0.54, n=310)
Effectsizes <- effect_size(eff.type="corr", r=0.68, n=310)
Effectsizes <- effect_size(eff.type="corr", r=0.43, n=310)
Effectsizes <- effect_size(eff.type="corr", r=0.42, n=310)

# Echebarria-Echabe & Fernandez-Guede (2006)
Effectsizes <- effect_size(eff.type="f.test", f=5.73, n.1=103, n.2=103)
Effectsizes <- effect_size(eff.type="f.test", f=4.80, n.1=103, n.2=103)
Effectsizes <- effect_size(eff.type="f.test", f=5.77, n.1=103, n.2=103)
Effectsizes <- effect_size(eff.type="f.test", f=2.92, n.1=103, n.2=103)
Effectsizes <- effect_size(eff.type="f.test", f=4.91, n.1=103, n.2=103)

# Echebarria-Echabe (2009)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=2.07, m.2=1.60, 
                           sd.1=0.86, sd.2=0.62, 
                           n.1=30, n.2=30) # Outgoup evaluation
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.80, m.2=2.53, 
                           sd.1=1.37, sd.2=1.43, 
                           n.1=30, n.2=30) # Intergroup bias

# Ellis et al. (2015)
Effectsizes <- effect_size(eff.type="corr", r=-0.243, n=77) 
Effectsizes <- effect_size(eff.type="corr", r=0.364, n=77) 
Effectsizes <- effect_size(eff.type="corr", r=0.168, n=77) 
Effectsizes <- effect_size(eff.type="corr", r=0.285, n=77) 
Effectsizes <- effect_size(eff.type="corr", r=-0.273, n=77) 
Effectsizes <- effect_size(eff.type="corr", r=0.131, n=77) 

# Elsayed & de Grip (2018)
# a) muslims
Effectsizes <- effect_size(eff.type="means", 
                           m.1=-0.004, m.2=-0.51, 
                           sd.1=0.51, sd.2=0.70, 
                           n.1=140, n.2=140) 
# b) non-muslims
Effectsizes <- effect_size(eff.type="means", 
                           m.1=0.11, m.2=-0.11, 
                           sd.1=0.57, sd.2=0.64, 
                           n.1=76, n.2=76) 

# Eshel & Kimhi (2016)
Effectsizes <- effect_size(eff.type="corr", r=-0.377, n=510) 
Effectsizes <- effect_size(eff.type="corr", r=-0.409, n=510) 
Effectsizes <- effect_size(eff.type="corr", r=-0.288, n=510) 

# Falk & Kenski (2016)
Effectsizes <- effect_size(eff.type="prop", 
                           p1=0.423, p2=0.262, 
                           n.ab=189, n.cd=488)
Effectsizes <- effect_size(eff.type="prop", 
                           p1=0.162, p2=0.127, 
                           n.ab=189, n.cd=488)

# Feinstein (2018)
Effectsizes <- effect_size(eff.type="corr", r=0.121, n=758)
Effectsizes <- effect_size(eff.type="corr", r=0.070, n=758)
Effectsizes <- effect_size(eff.type="corr", r=0.063, n=758)
Effectsizes <- effect_size(eff.type="corr", r=0.010, n=758)
Effectsizes <- effect_size(eff.type="corr", r=0.050, n=758)
Effectsizes <- effect_size(eff.type="corr", r=0.043, n=758)
Effectsizes <- effect_size(eff.type="corr", r=0.036, n=758)
Effectsizes <- effect_size(eff.type="corr", r=0.047, n=758)
Effectsizes <- effect_size(eff.type="corr", r=0.071, n=758)
Effectsizes <- effect_size(eff.type="corr", r=0.012, n=758)
Effectsizes <- effect_size(eff.type="corr", r=0.082, n=758)
Effectsizes <- effect_size(eff.type="corr", r=0.009, n=758)
Effectsizes <- effect_size(eff.type="corr", r=0.005, n=758)
Effectsizes <- effect_size(eff.type="corr", r=-0.009, n=758)
Effectsizes <- effect_size(eff.type="corr", r=0.137, n=758)
Effectsizes <- effect_size(eff.type="corr", r=-0.035, n=758)
Effectsizes <- effect_size(eff.type="corr", r=0.058, n=758)
Effectsizes <- effect_size(eff.type="corr", r=0.013, n=758)
Effectsizes <- effect_size(eff.type="corr", r=-0.035, n=758)
Effectsizes <- effect_size(eff.type="corr", r=0.016, n=758)

# Ferrin, Mancosu, & Cappiali (2018)
Effectsizes <- effect_size(eff.type="corr", r=OR2cor(exp(0.15)), n=26064)
Effectsizes <- effect_size(eff.type="corr", r=OR2cor(exp(0.06)), n=26121)

# Ferwerda, Flynn, & Horiuchi (2017)  
# main effects calculated using replication data
Effectsizes <- effect_size(eff.type="means", 
                           m.1=5.21, m.2=4.82, 
                           sd.1=3.36, sd.2=3.37,
                           n.1=761, n.2=765)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=5.55, m.2=5.04, 
                           sd.1=3.29, sd.2=3.34,
                           n.1=761, n.2=765)

# Finseraas & Listhaug (2013)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=0.79, m.2=0.80, 
                           sd.1=0.41, sd.2=0.40, 
                           n.1=2471, n.2=2945) 
Effectsizes <- effect_size(eff.type="means", 
                           m.1=0.70, m.2=0.69, 
                           sd.1=0.46, sd.2=0.46, 
                           n.1=2932, n.2=2469) 
Effectsizes <- effect_size(eff.type="means", 
                           m.1=0.62, m.2=0.61, 
                           sd.1=0.49, sd.2=0.49, 
                           n.1=2932, n.2=2465) 
Effectsizes <- effect_size(eff.type="means", 
                           m.1=5.03, m.2=5.05, 
                           sd.1=2.05, sd.2=2.05, 
                           n.1=2361, n.2=2820) 

# Finseraas, Jakobsson, & (2011)
Effectsizes <- effect_size(eff.type="t.test", t=1.87, n.1=10376, n.2=9183)
Effectsizes <- effect_size(eff.type="corr", r=beta2cor(0.043,1), 
                           n=561+457)#spain
Effectsizes <- effect_size(eff.type="corr", r=beta2cor(0.036,1),
                           n=682+638)#slovakia
Effectsizes <- effect_size(eff.type="corr", r=beta2cor(0.031,1),
                           n=364+641)#norway
Effectsizes <- effect_size(eff.type="corr", r=beta2cor(0.025,1),
                           n=399+540)#luxemb.
Effectsizes <- effect_size(eff.type="corr", r=beta2cor(0.022,1),
                           n=320+931)#UK
Effectsizes <- effect_size(eff.type="corr", r=beta2cor(0.019,1),
                           n=766+1354)#Czech
Effectsizes <- effect_size(eff.type="corr", r=beta2cor(0.013,1),
                           n=576+860)#finland
Effectsizes <- effect_size(eff.type="corr", r=beta2cor(0.013,1),
                           n=683+406)#belgium
Effectsizes <- effect_size(eff.type="corr", r=beta2cor(0.011,1),
                           n=728+927)#poland
Effectsizes <- effect_size(eff.type="corr", r=beta2cor(0.009,1),
                           n=637+784)#sweden
Effectsizes <- effect_size(eff.type="corr", r=beta2cor(0.009,1),
                           n=810+488)#slovenia
Effectsizes <- effect_size(eff.type="corr", r=beta2cor(0.002,1),
                           n=540+555)#switz.
Effectsizes <- effect_size(eff.type="corr", r=beta2cor(-0.011,0),
                           n=453+583)#germany
Effectsizes <- effect_size(eff.type="corr", r=beta2cor(-0.011,0),
                           n=252+620)#netherl.
Effectsizes <- effect_size(eff.type="corr", r=beta2cor(-0.027,0),
                           n=547+370)#estonia
Effectsizes <- effect_size(eff.type="corr", r=beta2cor(-0.030,0),
                           n=865+222)#denmark

# Fischer et al. (2010)
Effectsizes <- effect_size(eff.type="coh.d", d=0.50, n.1=40, n.2=40)
Effectsizes <- effect_size(eff.type="coh.d", d=0.74, n.1=40, n.2=40)
Effectsizes <- effect_size(eff.type="coh.d", d=0.83, n.1=11, n.2=11)

# Fischer et al. (2010)
Effectsizes <- effect_size(eff.type="corr", r=0.26, n=80) 

# Fischer, Oswald & Seiler (2013)
Effectsizes <- effect_size(eff.type="a.fes", f=5.22, 
                           n.1=110, n.2=110, R=0.31, q=1)
Effectsizes <- effect_size(eff.type="a.fes", f=6.62, 
                           n.1=110, n.2=110, R=0.30, q=1)     
Effectsizes <- effect_size(eff.type="a.fes", f=10.57, 
                           n.1=110, n.2=110, R=0.33, q=1)


# Fisk, Merolla, & Ramos (2019)
# main effects calculated using replication data
# a) France
# Fisk2019_france <- read_dta("https://www.dropbox.com/s/7wb2v0rxa91ij7o/Fisk2019_france.dta?dl=1")
# Fisk2019_france <- Fisk2019_france %>% drop_na(treatment, DRONES)
# describeBy(Fisk2019_france$DRONES, group=Fisk2019_france$treatment)
# control group: neutral
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.27, m.2=4.74, 
                           sd.1=1.84, sd.2=1.68, 
                           n.1=157, n.2=172)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.28, m.2=4.74, 
                           sd.1=1.90, sd.2=1.68, 
                           n.1=164, n.2=172)
# control group: economic threat
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.27, m.2=4.62, 
                           sd.1=1.84, sd.2=1.95, 
                           n.1=157, n.2=162)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.28, m.2=4.62, 
                           sd.1=1.90, sd.2=1.95, 
                           n.1=164, n.2=162)
# b) Turkey
# Fisk2019_turkey <- read_dta("https://www.dropbox.com/s/pbjfq04d1dskmsh/Fisk2019_turkey.dta?dl=1")
# Fisk2019_turkey <- Fisk2019_turkey %>% drop_na(treatment, DRONES1)
# describeBy(Fisk2019_turkey$DRONES1, group=Fisk2019_turkey$treatment)
# control group: neutral
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.04, m.2=2.68,
                           sd.1=1.96, sd.2=1.91, 
                           n.1=181, n.2=175)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=2.81, m.2=2.68, 
                           sd.1=1.93, sd.2=1.91, 
                           n.1=170, n.2=175)
# control group: economic threat
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.04, m.2=2.74, 
                           sd.1=1.96, sd.2=1.83, 
                           n.1=181, n.2=162)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=2.81, m.2=2.74, 
                           sd.1=1.93, sd.2=1.83, 
                           n.1=170, n.2=162)
# c) US
# Fisk2019_us <- read_dta("https://www.dropbox.com/s/hvng3gtaxwnnz3u/Fisk2019_UK.dta?dl=1")
# Fisk2019_us <- Fisk2019_us %>% drop_na(treatment, drones)
# describeBy(Fisk2019_us$drones, group=Fisk2019_us$treatment)
# control group: neutral
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.86, m.2=4.53, 
                           sd.1=1.68, sd.2=1.67, 
                           n.1=179, n.2=172)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.81, m.2=4.53, 
                           sd.1=1.57, sd.2=1.67, 
                           n.1=170, n.2=172)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=5.05, m.2=4.53, 
                           sd.1=1.61, sd.2=1.67, 
                           n.1=165, n.2=172)
# control group: economic threat
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.86, m.2=5.04, 
                           sd.1=1.68, sd.2=1.63, 
                           n.1=179, n.2=172)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.81, m.2=5.04, 
                           sd.1=1.57, sd.2=1.63, 
                           n.1=179, n.2=172)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=5.05, m.2=5.04, 
                           sd.1=1.61, sd.2=1.63, 
                           n.1=179, n.2=172)

# Frink et al. (2004)
# a) Time 1
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.75, m.2=4.72, 
                           sd.1=1.00, sd.2=0.94, 
                           n.1=228, n.2=129) 
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.64, m.2=4.61, 
                           sd.1=1.07, sd.2=1.01, 
                           n.1=228, n.2=129) 
# b) Time 2
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.75, m.2=4.73, 
                           sd.1=1.00, sd.2=0.96, 
                           n.1=228, n.2=123) 
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.64, m.2=4.68, 
                           sd.1=1.07, sd.2=1.00, 
                           n.1=228, n.2=123) 

# Frissen et al. (2018)
Effectsizes <- effect_size(eff.type="corr", r=-0.082, n=703) 
Effectsizes <- effect_size(eff.type="corr", r=0.037, n=703) 
Effectsizes <- effect_size(eff.type="corr", r=0.058, n=703) 
Effectsizes <- effect_size(eff.type="corr", r=0.117, n=703) 
Effectsizes <- effect_size(eff.type="corr", r=0.132, n=703) 
Effectsizes <- effect_size(eff.type="corr", r=0.134, n=703)
Effectsizes <- effect_size(eff.type="corr", r=0.197, n=703) 
Effectsizes <- effect_size(eff.type="corr", r=0.170, n=703) 
Effectsizes <- effect_size(eff.type="corr", r=-0.184, n=703)
Effectsizes <- effect_size(eff.type="corr", r=-0.044, n=703)

# Gadarian (2010)
# calculated via mail correspondence
# DV: militarism
Effectsizes <- effect_size(eff.type="prop", 
                           p1=0.360, p2=0.310, 
                           n.ab=424, n.cd=390) # Visuals
Effectsizes <- effect_size(eff.type="prop", 
                           p1=0.329, p2=0.310, 
                           n.ab=411, n.cd=390) # Non-visuals
# DV: Sudan
Effectsizes <- effect_size(eff.type="prop", 
                           p1=0.337, p2=0.319, 
                           n.ab=185, n.cd=176) # Visuals
Effectsizes <- effect_size(eff.type="prop", 
                           p1=0.344, p2=0.319, 
                           n.ab=194, n.cd=176) # Non-visuals
# DV: Iraq
Effectsizes <- effect_size(eff.type="means", 
                           m.1=.290, m.2=.281969, 
                           sd.1=.34192, sd.2=.33099,
                           n.1=425, n.2=391) # Visuals
Effectsizes <- effect_size(eff.type="means", 
                           m.1=.291971, m.2=.281969, 
                           sd.1=.34558, sd.2=.33099,
                           n.1=411, n.2=391) # Non-visuals
# DV: Government approval
Effectsizes <- effect_size(eff.type="means", 
                           m.1=0.024241, m.2=-0.0014651, 
                           sd.1=0.885111, sd.2=0.90385927,
                           n.1=424, n.2=390) # Visuals
Effectsizes <- effect_size(eff.type="means", 
                           m.1=-0.0236749, m.2=-0.0014651, 
                           sd.1=0.92421899, sd.2=0.90385927,
                           n.1=410, n.2=390) # Non-visuals
# DV: Liberties
Effectsizes <- effect_size(eff.type="means", 
                           m.1=0.375294, m.2=0.384860, 
                           sd.1=0.328539, sd.2=0.340608,
                           n.1=425, n.2=393) # Visuals
Effectsizes <- effect_size(eff.type="means", 
                           m.1=0.383212, m.2=0.384860, 
                           sd.1=0.33457, sd.2=0.340608,
                           n.1=411, n.2=393) # Non-visuals

# Garzon (2018)
Effectsizes <- effect_size(eff.type="corr", r=0.11, n=1200) 
Effectsizes <- effect_size(eff.type="corr", r=0.08, n=1200) 
Effectsizes <- effect_size(eff.type="corr", r=0.08, n=1200) 
Effectsizes <- effect_size(eff.type="corr", r=0.09, n=1200) 
Effectsizes <- effect_size(eff.type="corr", r=0.08, n=1200) 

# Gilboa (1990)
Effectsizes <- effect_size(eff.type="prop", 
                           p1=0.62, p2=0.57, 
                           n.ab=508, n.cd=508)

# Giner-Sorolla & Maitner (2013)
# a) Study 1
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.59, m.2=2.75, 
                           sd.1=1.19, sd.2=1.33, 
                           n.1=30, n.2=30) 
# b) Study 2: Aggregate across power conditions
Effectsizes <- effect_size(eff.type="p.val", 
                           p=0.001, n.1=352, n.2=352, tail="two") 
Effectsizes <- effect_size(eff.type="p.val", 
                           p=0.001, n.1=352, n.2=352, tail="two") 
Effectsizes <- effect_size(eff.type="p.val", 
                           p=0.001, n.1=352, n.2=352, tail="two") 
Effectsizes <- effect_size(eff.type="p.val", 
                           p=0.001, n.1=352, n.2=352, tail="two") 
Effectsizes <- effect_size(eff.type="p.val",
                           p=0.006, n.1=352, n.2=352, tail="two") 
Effectsizes <- effect_size(eff.type="p.val", 
                           p=0.001, n.1=352, n.2=352, tail="two") 

# Godefroidt & Langer (2017)
# CHRISTIAN subsample
# GodefroidtLanger2017_C <- read_sav("https://www.dropbox.com/s/2dvci8zqzcs6nxv/Godefroidt%20%26%20Langer%20%282017%29%20Unpublished%20Dataset_Christians.sav?dl=1") # to replicate results
# a) Nationalism
Effectsizes <- effect_size(eff.type="corr", r=-0.103, n=1973) #exposure
Effectsizes <- effect_size(eff.type="corr", r=-0.047, n=1970) #worry
Effectsizes <- effect_size(eff.type="corr", r=-0.021, n=1953) #anger
Effectsizes <- effect_size(eff.type="corr", r=-0.042, n=1954) #fear
# b) Cognitive stereotypes
Effectsizes <- effect_size(eff.type="corr", r=-0.060, n=1899) #exposure
Effectsizes <- effect_size(eff.type="corr", r=0.036, n=1896) #worry
Effectsizes <- effect_size(eff.type="corr", r=-0.016, n=1882) #anger
Effectsizes <- effect_size(eff.type="corr", r=0.011, n=1883) #fear
# c) Affective stereotypes 1
Effectsizes <- effect_size(eff.type="corr", r=0.109, n=1948) #exposure
Effectsizes <- effect_size(eff.type="corr", r=0.102, n=1945) #worry
Effectsizes <- effect_size(eff.type="corr", r=0.050, n=1930) #anger
Effectsizes <- effect_size(eff.type="corr", r=0.135, n=123) #fear
# d) Affective stereotypes 2
Effectsizes <- effect_size(eff.type="corr", r=0.065, n=1973) #exposure
Effectsizes <- effect_size(eff.type="corr", r=0.119, n=1970) #worry
Effectsizes <- effect_size(eff.type="corr", r=0.119, n=1954) #fear
Effectsizes <- effect_size(eff.type="corr", r=0.125, n=1955) #anger
# MUSLIM subsample
# GodefroidtLanger2017_M <- read_sav("https://www.dropbox.com/s/gc8wajj2m2i147y/Godefroidt%20%26%20Langer%20%282017%29%20Unpublished%20Dataset_Muslims.sav?dl=1") # to replicate results
# a) Nationalism
Effectsizes <- effect_size(eff.type="corr", r=0.118, n=503) #exposure
Effectsizes <- effect_size(eff.type="corr", r=0.171, n=502) #worry
Effectsizes <- effect_size(eff.type="corr", r=0.154, n=502) #anger
Effectsizes <- effect_size(eff.type="corr", r=0.134, n=502) #fear
# b) Cognitive stereotypes
Effectsizes <- effect_size(eff.type="corr", r=0.002, n=489) #exposure
Effectsizes <- effect_size(eff.type="corr", r=-0.077, n=488) #worry
Effectsizes <- effect_size(eff.type="corr", r=-0.121, n=488) #anger
Effectsizes <- effect_size(eff.type="corr", r=-0.095, n=488) #fear
# c) Affective stereotypes 1
Effectsizes <- effect_size(eff.type="corr", r=-0.100, n=496) #exposure
Effectsizes <- effect_size(eff.type="corr", r=-0.013, n=496) #worry
Effectsizes <- effect_size(eff.type="corr", r=-0.117, n=495) #anger
Effectsizes <- effect_size(eff.type="corr", r=-0.087, n=495) #fear
# d) Affective stereotypes 2
Effectsizes <- effect_size(eff.type="corr", r=0.056, n=502) #exposure
Effectsizes <- effect_size(eff.type="corr", r=0.033, n=501) #worry
Effectsizes <- effect_size(eff.type="corr", r=-0.010, n=501) #anger
Effectsizes <- effect_size(eff.type="corr", r=-0.015, n=501) #fear

# Godefroidt & Van Assche (2019)
#GodefroidtVanAssche2019 <- read_sav("https://www.dropbox.com/s/wxsgfdogcsbjtzl/Godefroidt%20%26%20Van%20Assche%20%282019%29%20Unpublished%20dataset.sav?dl=1") # to replicate results
# a) immigration policies
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.165000, m.2=4.16748, 
                           sd.1=1.575212, sd.2=1.550045, 
                           n.1=200, n.2=205)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.883838, m.2=4.16748, 
                           sd.1=1.491649, sd.2=1.550045, 
                           n.1=198, n.2=205)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.955224, m.2=4.16748, 
                           sd.1=1.450167, sd.2=1.550045, 
                           n.1=201, n.2=205)
# b) domestic counterterrorism policies
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.093333, m.2=2.839024, 
                           sd.1=1.530687, sd.2=1.325794, 
                           n.1=200, n.2=205)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=2.887205, m.2=2.839024, 
                           sd.1=1.432506, sd.2=1.325794, 
                           n.1=198, n.2=205)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=2.779436, m.2=2.839024, 
                           sd.1=1.266008, sd.2=1.325794, 
                           n.1=201, n.2=205)
# c) foreign counterterrorism policies
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.785, m.2=3.74878, 
                           sd.1=1.362187, sd.2=1.241268, 
                           n.1=200, n.2=205)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.613636, m.2=3.74878, 
                           sd.1=1.270903, sd.2=1.241268, 
                           n.1=200, n.2=205)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.636816, m.2=3.74878, 
                           sd.1=1.22574, sd.2=1.241268, 
                           n.1=200, n.2=205)
# d) environmental policies
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.884146, m.2=4.91125, 
                           sd.1=1.072863, sd.2=1.07153, 
                           n.1=205, n.2=200)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.884146, m.2=4.938131, 
                           sd.1=1.072863, sd.2=1.083575, 
                           n.1=205, n.2=198)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.884146, m.2=4.962687, 
                           sd.1=1.072863, sd.2=1.134064, 
                           n.1=205, n.2=201)
# e) national political trust
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.931250, m.2=4.014634, 
                           sd.1=2.092541, sd.2=2.124156, 
                           n.1=200, n.2=205)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.255076, m.2=4.014634, 
                           sd.1=1.950759, sd.2=2.124156, 
                           n.1=197, n.2=205)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.185323, m.2=4.014634, 
                           sd.1=2.081957, sd.2=2.124156, 
                           n.1=201, n.2=205)
# f) international political trust
Effectsizes <- effect_size(eff.type="means", 
                           m.1=5.611667, m.2=5.878049, 
                           sd.1=2.320924, sd.2=1.986002, 
                           n.1=200, n.2=205)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=5.853535, m.2=5.878049, 
                           sd.1=2.034221, sd.2=1.986002, 
                           n.1=198, n.2=205)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=5.895522, m.2=5.878049, 
                           sd.1=2.170926, sd.2=1.986002, 
                           n.1=201, n.2=205)
# g) feeling thermometer: moslims
Effectsizes <- effect_size(eff.type="means", 
                           m.1=52.985294, m.2=56.070707, 
                           sd.1=26.898730, sd.2=28.299971, 
                           n.1=204, n.2=198)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=52.985294, m.2=57.619289, 
                           sd.1=26.898730, sd.2=25.430221, 
                           n.1=204, n.2=197)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=52.985294, m.2=55.910448, 
                           sd.1=26.898730, sd.2=26.661807, 
                           n.1=204, n.2=201)

# Godefroidt, Schroyens, & Langer (2018)
# STUDENT subsample
#Godefroidtetal2018_Student <- read_sav("https://www.dropbox.com/s/1th3lxnppbv6m05/Godefroidt%2C%20Schroyens%2C%20%26%20Langer%20%282018%29%20Unpublished%20dataset_Students.sav?dl=1") # to replicate results
# a) Feeling thermometer: Muslims
Effectsizes <- effect_size(eff.type="corr", r=-0.044883, n=2312) #worry
Effectsizes <- effect_size(eff.type="corr", r=-0.134857, n=2309) #sadness
Effectsizes <- effect_size(eff.type="corr", r=0.062, n=2311) #anger
Effectsizes <- effect_size(eff.type="corr", r=-0.027, n=2309) #fear
# b) Feeling thermometer: Refugees
Effectsizes <- effect_size(eff.type="corr", r=-0.000, n=2301) #worry
Effectsizes <- effect_size(eff.type="corr", r=-0.122, n=2298) #sadness
Effectsizes <- effect_size(eff.type="corr", r=0.037, n=2300) #anger
Effectsizes <- effect_size(eff.type="corr", r=-0.012, n=2298) #fear
# c) Prejudice towards immigrants
Effectsizes <- effect_size(eff.type="corr", r=-0.054, n=2338) #worry
Effectsizes <- effect_size(eff.type="corr", r=-0.162, n=2335) #sadness
Effectsizes <- effect_size(eff.type="corr", r=0.057, n=2337) #anger
Effectsizes <- effect_size(eff.type="corr", r=-0.076, n=2335) #fear
# d) Outgroup collaboration in school
Effectsizes <- effect_size(eff.type="corr", r=-0.118, n=2338) #worry
Effectsizes <- effect_size(eff.type="corr", r=-0.214, n=2335) #sadness
Effectsizes <- effect_size(eff.type="corr", r=-0.023, n=2337) #anger
Effectsizes <- effect_size(eff.type="corr", r=-0.118, n=2335) #fear
# e) Political trust
Effectsizes <- effect_size(eff.type="corr", r=-0.152, n=2330) #worry
Effectsizes <- effect_size(eff.type="corr", r=-0.200, n=2327) #sadness
Effectsizes <- effect_size(eff.type="corr", r=-0.123, n=2329) #anger
Effectsizes <- effect_size(eff.type="corr", r=-0.200, n=2327) #fear
# f) Political orientation
Effectsizes <- effect_size(eff.type="corr", r=-0.032, n=1664) #worry
Effectsizes <- effect_size(eff.type="corr", r=-0.105, n=1661) #sadness
Effectsizes <- effect_size(eff.type="corr", r=0.118, n=1662) #anger
Effectsizes <- effect_size(eff.type="corr", r=-0.047, n=1661) #fear
# TEACHER subsample
#Godefroidtetal2018_Student <- read_sav("https://www.dropbox.com/s/c5g23f2civaksrp/Godefroidt%2C%20Schroyens%2C%20%26%20Langer%20%282018%29%20Unpublished%20dataset_Teachers.sav?dl=1") # to replicate results
# a) Feeling thermometer: Muslims
Effectsizes <- effect_size(eff.type="corr", r=0.273089, n=850) #worry
Effectsizes <- effect_size(eff.type="corr", r=0.033368, n=850) #sadness
Effectsizes <- effect_size(eff.type="corr", r=0.225837, n=850) #anger
Effectsizes <- effect_size(eff.type="corr", r=0.128372, n=850) #fear
# b) Feeling thermometer: Refugees
Effectsizes <- effect_size(eff.type="corr", r=0.263705, n=849) #worry
Effectsizes <- effect_size(eff.type="corr", r=0.004529, n=849) #sadness
Effectsizes <- effect_size(eff.type="corr", r=0.228492, n=849) #anger
Effectsizes <- effect_size(eff.type="corr", r=0.145896, n=849) #fear
# c) Prejudice towards immigrants
Effectsizes <- effect_size(eff.type="corr", r=0.282966, n=864) #worry
Effectsizes <- effect_size(eff.type="corr", r=0.061498, n=861) #sadness
Effectsizes <- effect_size(eff.type="corr", r=0.229941, n=861) #anger
Effectsizes <- effect_size(eff.type="corr", r=0.155593, n=861) #fear
# e) Political trust
Effectsizes <- effect_size(eff.type="corr", r=0.087318, n=867) #worry
Effectsizes <- effect_size(eff.type="corr", r=-0.035850, n=863) #sadness
Effectsizes <- effect_size(eff.type="corr", r=0.103546, n=863) #anger
Effectsizes <- effect_size(eff.type="corr", r=-0.000821, n=863) #fear
# f) Political orientation
Effectsizes <- effect_size(eff.type="corr", r=0.279094, n=788) #worry
Effectsizes <- effect_size(eff.type="corr", r=0.018899, n=787) #sadness
Effectsizes <- effect_size(eff.type="corr", r=0.270816, n=787) #anger
Effectsizes <- effect_size(eff.type="corr", r=0.233206, n=787) #fear

# De Zavala et al. (2010)
Effectsizes <- effect_size(eff.type="corr", r=0.30, n=187) # Out-group hostility
Effectsizes <- effect_size(eff.type="corr", r=0.21, n=187) # Political conservatism

# Goodwin et al. (2017): Regression coefficients
Effectsizes <- effect_size(eff.type="corr", r=beta2cor(0.10,1), n=3666)
Effectsizes <- effect_size(eff.type="corr", r=beta2cor(0.09,1), n=3666)
Effectsizes <- effect_size(eff.type="corr", r=beta2cor(0.08,1), n=3666)
Effectsizes <- effect_size(eff.type="corr", r=beta2cor(0.08,1), n=3666)

# Goren & Neter (2016)
Effectsizes <- effect_size(eff.type="corr", r=0.26, n=263)
Effectsizes <- effect_size(eff.type="coh.d", d=1.06, n.1=159, n.2=104)
Effectsizes <- effect_size(eff.type="coh.d", d=0.77, n.1=159, n.2=104)

# Greenaway et al. (2014)
Effectsizes <- effect_size(eff.type="corr", r=beta2cor(0.03,1), n=87)
Effectsizes <- effect_size(eff.type="corr", r=beta2cor(0.10,1), n=2394)

# Grizzard et al. (2017)
# a) study 1
Effectsizes <- effect_size(eff.type="means", 
                           m.2=4.66, m.1=4.733333333, 
                           sd.2=0.94, sd.1=0.853333333,
                           n.2=81, n.1=234)
Effectsizes <- effect_size(eff.type="means", 
                           m.2=6.54, m.1=6.166666667, 
                           sd.2=1.92, sd.1=2.093333333,
                           n.2=81, n.1=234)
Effectsizes <- effect_size(eff.type="means", 
                           m.2=6.276666667, m.1=6.21, 
                           sd.2=2.086666667,sd.1=2.06,
                           n.2=234, n.1=81)
# b) study 2
Effectsizes <- effect_size(eff.type="means", 
                           m.2=4.74, m.1=4.686666667, 
                           sd.2=0.84, sd.1=0.85,
                           n.2=67, n.1=195)
Effectsizes <- effect_size(eff.type="means", 
                           m.2=5.86, m.1=6.266666667, 
                           sd.2=2.05, sd.1=2.213333333,
                           n.2=67, n.1=195)
Effectsizes <- effect_size(eff.type="means", 
                           m.2=6.65, m.1=6.44, 
                           sd.2=1.933333333, sd.1=2.01,
                           n.2=195, n.1=67)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=0.15, m.2=0.11, 
                           sd.1=1.4, sd.2=0.37,
                           n.1=195, n.2=67)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=2.936666667, m.2=2.58, 
                           sd.1=1.4966667, sd.2=1.49,
                           n.1=195, n.2=67)

# Gross, Canetti, & Vashdi (2017)
# a) Experiment 1
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.5, m.2=4.5, 
                           sd.1=1.23, sd.2=1.2,
                           n.1=346, n.2=334)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.5, m.2=4.61, 
                           sd.1=1.23, sd.2=1.09,
                           n.1=346, n.2=346)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.55, m.2=3.56, 
                           sd.1=1.28, sd.2=1.23,
                           n.1=346, n.2=334)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.55, m.2=3.87, 
                           sd.1=1.28, sd.2=1.15,
                           n.1=346, n.2=345)
# b) Experiment 2
#regulation attitudes
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.321, m.2=4.1409, 
                           sd.1=1.34166, sd.2=1.42793,
                           n.1=162, n.2=298)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.321, m.2=4.3316, 
                           sd.1=1.34166, sd.2=1.34517,
                           n.1=162, n.2=291)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.321, m.2=4.2583, 
                           sd.1=1.34166, sd.2=1.26604,
                           n.1=162, n.2=151)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.5926, m.2=3.6057, 
                           sd.1=1.41555, sd.2=1.44997,
                           n.1=162, n.2=298)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.5926, m.2=3.7818, 
                           sd.1=1.41555, sd.2=1.42017,
                           n.1=162, n.2=291)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.5926, m.2=3.5033, 
                           sd.1=1.41555, sd.2=1.56684,
                           n.1=162, n.2=151)
#confidence attitudes
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.1734, m.2=4.1409, 
                           sd.1=0.95354, sd.2=0.97529,
                           n.1=298, n.2=162)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.1506, m.2=4.1409, 
                           sd.1=0.95137, sd.2=0.97529,
                           n.1=291, n.2=162)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.1523, m.2=4.1409, 
                           sd.1=0.92715, sd.2=0.97529,
                           n.1=151, n.2=162)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.5134, m.2=3.5802, 
                           sd.1=1.12284, sd.2=1.05707,
                           n.1=298, n.2=162)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.5601, m.2=3.5802, 
                           sd.1=1.07769, sd.2=1.05707,
                           n.1=291, n.2=162)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.6887, m.2=3.5802, 
                           sd.1=1.07681, sd.2=1.05707,
                           n.1=151, n.2=162)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.5289, m.2=4.3716, 
                           sd.1=0.95504, sd.2=0.96597,
                           n.1=298, n.2=162)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.5753, m.2=4.3716, 
                           sd.1=0.94756, sd.2=0.96597,
                           n.1=291, n.2=162)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.5298, m.2=4.3716, 
                           sd.1=1.03877, sd.2=0.96597,
                           n.1=151, n.2=162)
#hawkish attitudes
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.64, m.2=4.47, 
                           sd.1=1.447, sd.2=1.566,
                           n.1=162, n.2=296)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.64, m.2=4.59, 
                           sd.1=1.447, sd.2=1.541,
                           n.1=162, n.2=291)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.64, m.2=4.65, 
                           sd.1=1.447, sd.2=1.524,
                           n.1=162, n.2=151)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.87, m.2=4.83, 
                           sd.1=1.271, sd.2=1.333,
                           n.1=162, n.2=296)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.87, m.2=4.87, 
                           sd.1=1.271, sd.2=1.326,
                           n.1=162, n.2=291)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.87, m.2=4.95, 
                           sd.1=1.271, sd.2=1.341,
                           n.1=162, n.2=151)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.80, m.2=3.84, 
                           sd.1=1.56, sd.2=1.60,
                           n.1=162, n.2=296)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.80, m.2=3.84, 
                           sd.1=1.56, sd.2=1.67,
                           n.1=162, n.2=291)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.80, m.2=4.01, 
                           sd.1=1.56, sd.2=1.596,
                           n.1=162, n.2=151)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.99, m.2=3.98, 
                           sd.1=1.522, sd.2=1.557,
                           n.1=162, n.2=296)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.99, m.2=4.10, 
                           sd.1=1.522, sd.2=1.502,
                           n.1=162, n.2=291)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.99, m.2=4.38, 
                           sd.1=1.522, sd.2=1.436,
                           n.1=162, n.2=151)
# c) Experiment 3
#regulation attitudes
es <- esc_mean_sd(grp1m=4.26, grp1sd=1.22, grp1n=454, 
                  grp2m=4.17, grp2sd=1.27, grp2n=454, r=0.63, es.type="r")
z.se <- round(1/sqrt((es$totaln)-3), digits=6)
z.var <- z.se^2
pwr <- pwr.r.test(n=es$totaln, r=es$es, sig.level=.05)
results <- data.frame(c(es$es, es$var, es$se,
                        es$zr, z.var, z.se, es$totaln, pwr$power))
results <- transpose(results)
Effectsizes <- rbind(Effectsizes, results)
es <- esc_mean_sd(grp1m=3.07, grp1sd=1.09, grp1n=454, 
                  grp2m=3.15, grp2sd=1.14, grp2n=454, r=0.75, es.type="r")
z.se <- round(1/sqrt((es$totaln)-3), digits=6)
z.var <- z.se^2
pwr <- pwr.r.test(n=es$totaln, r=es$es, sig.level=.05)
results <- data.frame(c(es$es, es$var, es$se,
                        es$zr, z.var, z.se, es$totaln, pwr$power))
results <- transpose(results)
Effectsizes <- rbind(Effectsizes, results)
#hawkish attitudes
es <- esc_mean_sd(grp1m=4.75, grp1sd=1.18, grp1n=455, 
                  grp2m=4.71, grp2sd=1.17, grp2n=455, r=0.38, es.type="r")
z.se <- round(1/sqrt((es$totaln)-3), digits=6)
z.var <- z.se^2
pwr <- pwr.r.test(n=es$totaln, r=es$es, sig.level=.05)
results <- data.frame(c(es$es, es$var, es$se,
                        es$zr, z.var, z.se, es$totaln, pwr$power))
results <- transpose(results)
Effectsizes <- rbind(Effectsizes, results)
es <- esc_mean_sd(grp1m=4.44, grp1sd=1.36, grp1n=455, 
                  grp2m=4.12, grp2sd=1.49, grp2n=455, r=0.48, es.type="r")
z.se <- round(1/sqrt((es$totaln)-3), digits=6)
z.var <- z.se^2
pwr <- pwr.r.test(n=es$totaln, r=es$es, sig.level=.05)
results <- data.frame(c(es$es, es$var, es$se,
                        es$zr, z.var, z.se, es$totaln, pwr$power))
results <- transpose(results)
Effectsizes <- rbind(Effectsizes, results)
es <- esc_mean_sd(grp1m=3.06, grp1sd=1.52, grp1n=455, 
                  grp2m=2.95, grp2sd=1.54, grp2n=455, r=0.62, es.type="r")
z.se <- round(1/sqrt((es$totaln)-3), digits=6)
z.var <- z.se^2
pwr <- pwr.r.test(n=es$totaln, r=es$es, sig.level=.05)
results <- data.frame(c(es$es, es$var, es$se,
                        es$zr, z.var, z.se, es$totaln, pwr$power))
results <- transpose(results)
Effectsizes <- rbind(Effectsizes, results)
es <- esc_mean_sd(grp1m=2.89, grp1sd=1.58, grp1n=455, 
                  grp2m=2.46, grp2sd=1.53, grp2n=455, r=0.66, es.type="r")
z.se <- round(1/sqrt((es$totaln)-3), digits=6)
z.var <- z.se^2
pwr <- pwr.r.test(n=es$totaln, r=es$es, sig.level=.05)
results <- data.frame(c(es$es, es$var, es$se,
                        es$zr, z.var, z.se, es$totaln, pwr$power))
results <- transpose(results)
Effectsizes <- rbind(Effectsizes, results)
#change sign (entered reversed above)
Effectsizes[702:738,1] <- -Effectsizes[702:738,1]
Effectsizes[702:738,4] <- -Effectsizes[702:738,4]

# Hall et al. (2009)
Effectsizes <- effect_size(eff.type="corr", r=0.09, n=961)
Effectsizes <- effect_size(eff.type="corr", r=0.01, n=961)
Effectsizes <- effect_size(eff.type="corr", r=0.21, n=961)
Effectsizes <- effect_size(eff.type="corr", r=0.21, n=961)
Effectsizes <- effect_size(eff.type="corr", r=0.16, n=961)
Effectsizes <- effect_size(eff.type="corr", r=0.20, n=961)
Effectsizes <- effect_size(eff.type="corr", r=0.09, n=961)
Effectsizes <- effect_size(eff.type="corr", r=-0.01, n=961)
Effectsizes <- effect_size(eff.type="corr", r=0.22, n=961)
Effectsizes <- effect_size(eff.type="corr", r=0.23, n=961)
Effectsizes <- effect_size(eff.type="corr", r=0.17, n=961)
Effectsizes <- effect_size(eff.type="corr", r=0.18, n=961)

# Hall et al. (2018)
# a) Amnesty
Effectsizes <- b_effect_size(eff.type="b", b=0.393, sdy=1.657, grp1n=250, grp2n=253) #never displaced (vs. still displaced)
Effectsizes <- b_effect_size(eff.type="b", b=0.876, sdy=1.657, grp1n=250, grp2n=504) #returnee (vs. still displaced)
Effectsizes <- b_effect_size(eff.type="b", b=1.071, sdy=1.657, grp1n=29, grp2n=706)  #physical injury
Effectsizes <- b_effect_size(eff.type="b", b=0.594, sdy=1.657, grp1n=366, grp2n=371) #lost a loved one
Effectsizes <- b_effect_size(eff.type="b", b=0.214, sdy=1.657, grp1n=49, grp2n=685)  #imprisonment
Effectsizes <- b_effect_size(eff.type="b", b=-0.867, sdy=1.657, grp1n=16, grp2n=719) #torture
Effectsizes <- b_effect_size(eff.type="b", b=0.003, sdy=1.657, grp1n=593, grp2n=148) #loss of property
# b) Forgiveness
Effectsizes <- b_effect_size(eff.type="b", b=0.126, sdy=1.501, grp1n=250, grp2n=253) #never displaced (vs. still displaced)
Effectsizes <- b_effect_size(eff.type="b", b=0.353, sdy=1.501, grp1n=250, grp2n=504) #returnee (vs. still displaced)
Effectsizes <- b_effect_size(eff.type="b", b=0.557, sdy=1.501, grp1n=29, grp2n=706)  #physical injury
Effectsizes <- b_effect_size(eff.type="b", b=0.277, sdy=1.501, grp1n=366, grp2n=371) #lost a loved one
Effectsizes <- b_effect_size(eff.type="b", b=0.453, sdy=1.501, grp1n=49, grp2n=685)  #imprisonment
Effectsizes <- b_effect_size(eff.type="b", b=-0.453, sdy=1.501, grp1n=16, grp2n=719) #torture
Effectsizes <- b_effect_size(eff.type="b", b=0.304, sdy=1.501, grp1n=593, grp2n=148) #loss of property
# c) Trail and punish
Effectsizes <- b_effect_size(eff.type="b", b=-0.005, sdy=0.800, grp1n=250, grp2n=253)#never displaced (vs. still displaced)
Effectsizes <- b_effect_size(eff.type="b", b=0.264, sdy=0.800, grp1n=250, grp2n=504) #returnee (vs. still displaced)
Effectsizes <- b_effect_size(eff.type="b", b=0.317, sdy=0.800, grp1n=29, grp2n=706)  #physical injury
Effectsizes <- b_effect_size(eff.type="b", b=0.255, sdy=0.800, grp1n=366, grp2n=371) #lost a loved one
Effectsizes <- b_effect_size(eff.type="b", b=0.349, sdy=0.800, grp1n=49, grp2n=685)  #Imprisonment
Effectsizes <- b_effect_size(eff.type="b", b=0.417, sdy=0.800, grp1n=16, grp2n=719)  #torture
Effectsizes <- b_effect_size(eff.type="b", b=0.132, sdy=0.800, grp1n=593, grp2n=148) #loss of property
# d) Pay compensation; war criminals
Effectsizes <- b_effect_size(eff.type="b", b=0.116, sdy=1.169, grp1n=250, grp2n=253) #displaced
Effectsizes <- b_effect_size(eff.type="b", b=-0.201, sdy=1.169, grp1n=250, grp2n=504) #returnee
Effectsizes <- b_effect_size(eff.type="b", b=-0.201, sdy=1.169, grp1n=29, grp2n=706)  #physical injury
Effectsizes <- b_effect_size(eff.type="b", b=-0.160, sdy=1.169, grp1n=366, grp2n=371) #lost a loved one
Effectsizes <- b_effect_size(eff.type="b", b=0.113, sdy=1.169, grp1n=49, grp2n=685)  #imprisonment
Effectsizes <- b_effect_size(eff.type="b", b=-1.286, sdy=1.169, grp1n=16, grp2n=719)  #torture
Effectsizes <- b_effect_size(eff.type="b", b=0.293, sdy=1.169, grp1n=593, grp2n=148) #loss of property

# Halperin et al. (2009)
Effectsizes <- effect_size(eff.type="corr", r= 0.10, n=772)

# Halperin, Canetti-Nisim, & Hirsch-Hoefler (2009)
Effectsizes <- effect_size(eff.type="corr", r=-0.46, n=847) # Hatred survey 1
Effectsizes <- effect_size(eff.type="corr", r=-0.50, n=453) # Hatred survey 2
Effectsizes <- effect_size(eff.type="corr", r=-0.31, n=847) # Anger survey 1
Effectsizes <- effect_size(eff.type="corr", r=-0.33, n=453) # Anger survey 2
Effectsizes <- effect_size(eff.type="corr", r=-0.29, n=847) # Fear survey 1
Effectsizes <- effect_size(eff.type="corr", r=-0.33, n=453) # Fear survey 2

# Halperin et al. (2013)
Effectsizes <- effect_size(eff.type="coh.d", d=0.89, n.1=31, n.2=63)
Effectsizes <- effect_size(eff.type="coh.d", d=-0.57, n.1=27, n.2=27)

# Hamama-raz et al. (2008)
Effectsizes <- effect_size(eff.type="corr", r= 0.09, n=276) 
Effectsizes <- effect_size(eff.type="corr", r= 0.08, n=276) 
Effectsizes <- effect_size(eff.type="corr", r= 0.06, n=276) 
Effectsizes <- effect_size(eff.type="corr", r=-0.01, n=1469) 
Effectsizes <- effect_size(eff.type="corr", r=-0.02, n=1469) 
Effectsizes <- effect_size(eff.type="corr", r=-0.06, n=1469) 

# Haridakis & Rubin (2005)
Effectsizes <- effect_size(eff.type="corr", r=beta2cor(0.06,1), n=218) 
Effectsizes <- effect_size(eff.type="corr", r=beta2cor(0.13,1), n=218)

# Hawi et al. (2019)
# a) feeling thermometers
Effectsizes <- effect_size(eff.type="corr", r=0.164, n=16328)
Effectsizes <- effect_size(eff.type="corr", r=0.120, n=16328)
Effectsizes <- effect_size(eff.type="corr", r=0.077, n=16328)
Effectsizes <- effect_size(eff.type="corr", r=0.084, n=16328)
Effectsizes <- effect_size(eff.type="corr", r=0.110, n=16328)
Effectsizes <- effect_size(eff.type="corr", r=0.138, n=16328)
# b) political variables
Effectsizes <- effect_size(eff.type="corr", r=0.199, n=16328)
Effectsizes <- effect_size(eff.type="corr", r=0.161, n=16328)
Effectsizes <- effect_size(eff.type="corr", r=0.053, n=16328)
Effectsizes <- effect_size(eff.type="corr", r=0.170, n=16328)

# Henderson-King et al. (2009)
Effectsizes <- effect_size(eff.type="corr", r=0.10, n=266) 
Effectsizes <- effect_size(eff.type="corr", r=0.39, n=266) 
Effectsizes <- effect_size(eff.type="corr", r=-0.07, n=266) 
Effectsizes <- effect_size(eff.type="corr", r=0.11, n=266) 
Effectsizes <- effect_size(eff.type="corr", r=0.40, n=266) 
Effectsizes <- effect_size(eff.type="corr", r=0.36, n=266) 
Effectsizes <- effect_size(eff.type="corr", r=0.28, n=266) 
Effectsizes <- effect_size(eff.type="corr", r=0.33, n=266) 
Effectsizes <- effect_size(eff.type="corr", r=0.18, n=266) 

# Hendricks (2019)
Effectsizes <- effect_size(eff.type="corr", r=0.52, n=410)
Effectsizes <- effect_size(eff.type="corr", r=0.34, n=410)
Effectsizes <- effect_size(eff.type="corr", r=0.34, n=410)
Effectsizes <- effect_size(eff.type="corr", r=0.57, n=410)

# Hetherington & Suhay (2011)
# correlations: partisanship and authoritarianism in 2006 and 2008
Effectsizes <- effect_size(eff.type="corr", r=0.08, n=1000) 
Effectsizes <- effect_size(eff.type="corr", r=0.16, n=1000) 
Effectsizes <- effect_size(eff.type="corr", r=0.07, n=1500) 
Effectsizes <- effect_size(eff.type="corr", r=0.18, n=1500) 

# Hindman (2008)
Effectsizes <- effect_size(eff.type="prop", p1=0.837, p2=0.507, n.ab=1500, n.cd=2001)

# Hirschberger et al. (2017)
# a) study 1
Effectsizes <- effect_size(eff.type="f.test", f=1.54, n.1=57, n.2=57) #US support
Effectsizes <- effect_size(eff.type="f.test", f=5.32, n.1=57, n.2=57) #no US support
Effectsizes <- effect_size(eff.type="f.test", f=0.08, n.1=57, n.2=57) #US support
Effectsizes <- effect_size(eff.type="f.test", f=6.77, n.1=57, n.2=57) #no US support
Effectsizes <- effect_size(eff.type="f.test", f=7.68, n.1=57, n.2=57)
# b) study 2
Effectsizes <- effect_size(eff.type="p.val", p=0.012, n.1=44, n.2=44, tail="two") # vs. human
Effectsizes <- effect_size(eff.type="p.val", p=0.002, n.1=44, n.2=44, tail="two") # vs. pain
# c) study 3
Effectsizes <- effect_size(eff.type="corr", r=beta2cor(0.03,1), n=478)
Effectsizes <- effect_size(eff.type="corr", r=beta2cor(0.09,1), n=478)

# Hirsch-Hoefler et al. (2014)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=2.58, m.2=2.12, 
                           sd.1=1.63, sd.2=1.51, 
                           n.1=830, n.2=797) 
Effectsizes <- effect_size(eff.type="corr", r=beta2cor(0.10,0), n=1627)

# Hitlan et al. (2014)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.51, m.2=4.48, 
                           sd.1=1.07, sd.2=1.27,
                           n.1=140, n.2=84)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.28, m.2=3.70, 
                           sd.1=1.47, sd.2=1.56, 
                           n.1=140, n.2=84)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.06, m.2=2.86, 
                           sd.1=1.36, sd.2=1.40, 
                           n.1=140, n.2=84)

# Hobfoll et al. (2006)
Effectsizes <- effect_size(eff.type="corr", r=0.040, n=185) 
Effectsizes <- effect_size(eff.type="corr", r=0.039, n=720) 
Effectsizes <- effect_size(eff.type="corr", r=-0.072, n=185) 
Effectsizes <- effect_size(eff.type="corr", r=0.158, n=720) 
Effectsizes <- effect_size(eff.type="corr", r=0.129, n=185) 
Effectsizes <- effect_size(eff.type="corr", r=0.007, n=720) 

# Holman et al. (2019)
# a) Experiment 1
Effectsizes <- effect_size(eff.type="means", 
                           m.1=32.424581, m.2=39.534626, 
                           sd.1=30.466387, sd.2=33.420628,
                           n.1=179, n.2=362)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=28.050279, m.2=22.404432, 
                           sd.1=34.625975, sd.2=30.52071,
                           n.1=179, n.2=362)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=17.011173, m.2=18.952646, 
                           sd.1=23.132007, sd.2=26.412827,
                           n.1=179, n.2=362)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=61.642458, m.2=63.911357, 
                           sd.1=32.877635, sd.2=33.94281,
                           n.1=179, n.2=362)
# b) Experiment 2
Effectsizes <- effect_size(eff.type="means", 
                           m.1=53.67382, m.2=51.004566, 
                           sd.1=36.4353, sd.2=36.594166,
                           n.1=250, n.2=228)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=31.806034, m.2=32.990826, 
                           sd.1=33.185056, sd.2=34.674655,
                           n.1=250, n.2=228)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.3198381, m.2=4.30837, 
                           sd.1=2.2879437, sd.2=2.1763147,
                           n.1=250, n.2=228)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.1902834, m.2=3.1497797, 
                           sd.1=2.2970967, sd.2=2.1808778,
                           n.1=250, n.2=228)

# Holman, Merolla, & Zechmeister (2011)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=2.91, m.2=2.32, 
                           sd.1=1.56, sd.2=1.37, 
                           n.1=49, n.2=48) # Gender-trait stereotype: terror vs. control
Effectsizes <- effect_size(eff.type="means", 
                           m.1=2.91, m.2=2.18, 
                           sd.1=1.56, sd.2=1.35, 
                           n.1=49, n.2=52) # Gender-trait stereotype: terror vs. economic threat
Effectsizes <- effect_size(eff.type="means", 
                           m.1=63.17, m.2=53.76, 
                           sd.1=27.62, sd.2=31.42,
                           n.1=48, n.2=49) # Clinton: terror vs. control
Effectsizes <- effect_size(eff.type="means", 
                           m.1=58.28, m.2=53.76, 
                           sd.1=31.16, sd.2=31.42,
                           n.1=52, n.2=49) # Clinton: terror vs. economic threat
Effectsizes <- effect_size(eff.type="means", 
                           m.1=51.09, m.2=43.06, 
                           sd.1=26.65, sd.2=26.57,
                           n.1=49, n.2=48) # Rice: terror vs. control
Effectsizes <- effect_size(eff.type="means", 
                           m.1=51.09, m.2=46.69, 
                           sd.1=26.65, sd.2=29.56,
                           n.1=49, n.2=52) # Rice: terror vs. economic threat					
Effectsizes <- effect_size(eff.type="means", 
                           m.1=42.61, m.2=22.79, 
                           sd.1=31.17, sd.2=26.06,
                           n.1=49, n.2=51) # Bush: terror vs. control
Effectsizes <- effect_size(eff.type="means", 
                           m.1=42.61, m.2=31.7, 
                           sd.1=31.17, sd.2=29.53,
                           n.1=49, n.2=51) # Bush: terror vs. economic threat					
Effectsizes <- effect_size(eff.type="means", 
                           m.1=53.4, m.2=50.06, 
                           sd.1=23.63, sd.2=26.58,
                           n.1=48, n.2=49) # Kerry: terror vs. control
Effectsizes <- effect_size(eff.type="means", 
                           m.1=49.48, m.2=50.06, 
                           sd.1=24.2, sd.2=26.58,
                           n.1=52, n.2=49) # Kerry: terror vs. economic threat

# Holman, Merolla, & Zechmeister (2016)
# Study 3
#Holmanetal2016 <- read_sav("https://www.dropbox.com/s/u7gsti6kl8ttp7o/Holman%20et%20al.%20%282016%29%20Replication.sav?dl=1") # to replicate results (watch out: put on weights! via dplyr)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=51.59, m.2=50.44, 
                           sd.1=21.204, sd.2=18.926, 
                           n.1=514, n.2=525) #feeling thermometer democratic candidate (calculated via replication material)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=47.94, m.2=51.52, 
                           sd.1=17.624, sd.2=18.437, 
                           n.1=525, n.2=514) #feeling thermometer republican candidate (calculated via replication material)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.52, m.2=3.56, 
                           sd.1=1.217, sd.2=1.292, 
                           n.1=531, n.2=543) #leadership democratic candidate (calculated via replication material)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.73, m.2=3.59, 
                           sd.1=1.193, sd.2=1.073, 
                           n.1=543, n.2=531) #leadership republican candidate (calculated via replication material)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.60, m.2=3.70, 
                           sd.1=1.277, sd.2=1.145, 
                           n.1=531, n.2=543) #trust democratic candidate (calculated via replication material)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.85, m.2=3.69, 
                           sd.1=1.050, sd.2=1.022, 
                           n.1=543, n.2=531) #trust republican candidate (calculated via replication material)

# Holman, Merolla, & Zechmeister (2017)
# a) H. Clinton
Effectsizes <- effect_size(eff.type="p.val", p=0.21, n.1=(201+191+177), n.2=194, tail="two")
Effectsizes <- effect_size(eff.type="p.val", p=0.15, n.1=(201+191+177), n.2=186, tail="two")
# b) N. Pelosi
Effectsizes <- effect_size(eff.type="p.val", p=0.04, n.1=(201+191+177), n.2=194, tail="one")
Effectsizes <- effect_size(eff.type="p.val", p=0.04, n.1=(201+191+177), n.2=186, tail="one")
# c) B. Obama
Effectsizes <- effect_size(eff.type="p.val", p=0.18, n.1=(201+191+177), n.2=194, tail="two")
Effectsizes <- effect_size(eff.type="p.val", p=0.31, n.1=(201+191+177), n.2=186, tail="two")
# reverse sign for democratic candidates
Effectsizes[868:873,1] <- -Effectsizes[868:873,1]
Effectsizes[868:873,4] <- -Effectsizes[868:873,4]
# d) C. Rice
Effectsizes <- effect_size(eff.type="p.val", p=0.38, n.1=(201+191+177), n.2=194, tail="two")
Effectsizes <- effect_size(eff.type="p.val", p=0.51, n.1=(201+191+177), n.2=186, tail="two")
# e) S. Palin
Effectsizes <- effect_size(eff.type="p.val", p=0.52, n.1=(201+191+177), n.2=194, tail="two")
Effectsizes <- effect_size(eff.type="p.val", p=0.74, n.1=(201+191+177), n.2=186, tail="two")
# f) M. Romney
Effectsizes <- effect_size(eff.type="p.val", p=0.80, n.1=(201+191+177), n.2=194, tail="two")
Effectsizes <- effect_size(eff.type="p.val", p=0.84, n.1=(201+191+177), n.2=186, tail="two")

# Horry & Wright (2009)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=710.22, m.2=640.25, 
                           sd.1=178.16, sd.2=101.40, 
                           n.1=47, n.2=50)

# Huesmann et al. (2012)
Effectsizes <- effect_size(eff.type="corr", r=0.06, n=269)
Effectsizes <- effect_size(eff.type="corr", r=-0.34, n=269)
Effectsizes <- effect_size(eff.type="corr", r=0.28, n=269)

# Iyer et al. (2014)
# a) Help for victims
Effectsizes <- effect_size(eff.type="corr", r=-0.29, n=235)
Effectsizes <- effect_size(eff.type="corr", r=-0.09, n=235)
Effectsizes <- effect_size(eff.type="corr", r=-0.16, n=235)
Effectsizes <- effect_size(eff.type="corr", r=-0.38, n=235)
Effectsizes <- effect_size(eff.type="corr", r=-0.23, n=235)
Effectsizes <- effect_size(eff.type="corr", r=-0.16, n=235)
Effectsizes <- effect_size(eff.type="corr", r=-0.15, n=235)
# b) Negotiations
Effectsizes <- effect_size(eff.type="corr", r=-0.04, n=235)
Effectsizes <- effect_size(eff.type="corr", r=0.03, n=235)
Effectsizes <- effect_size(eff.type="corr", r=-0.12, n=235)
Effectsizes <- effect_size(eff.type="corr", r=-0.02, n=235)
Effectsizes <- effect_size(eff.type="corr", r=0.14, n=235)
Effectsizes <- effect_size(eff.type="corr", r=-0.28, n=235)
Effectsizes <- effect_size(eff.type="corr", r=-0.02, n=235)
# b) Aggression
Effectsizes <- effect_size(eff.type="corr", r=0.03, n=235)
Effectsizes <- effect_size(eff.type="corr", r=0.16, n=235)
Effectsizes <- effect_size(eff.type="corr", r=-0.08, n=235)
Effectsizes <- effect_size(eff.type="corr", r=-0.05, n=235)
Effectsizes <- effect_size(eff.type="corr", r=0.30, n=235)
Effectsizes <- effect_size(eff.type="corr", r=0.13, n=235)
Effectsizes <- effect_size(eff.type="corr", r=0.14, n=235)

# Iyer et al. (2015)
# a) study 1
Effectsizes <- effect_size(eff.type="corr", r=-0.36, n=106)
Effectsizes <- effect_size(eff.type="corr", r=-0.25, n=106)
Effectsizes <- effect_size(eff.type="corr", r=0.22, n=106)
Effectsizes <- effect_size(eff.type="corr", r=0.01, n=106)
Effectsizes <- effect_size(eff.type="corr", r=0.04, n=106)
Effectsizes <- effect_size(eff.type="corr", r=0.19, n=106)
# b) study 2
Effectsizes <- effect_size(eff.type="corr", r=-0.28, n=332)
Effectsizes <- effect_size(eff.type="corr", r=-0.33, n=332)
Effectsizes <- effect_size(eff.type="corr", r=0.33, n=332)
Effectsizes <- effect_size(eff.type="corr", r=0.12, n=332)
Effectsizes <- effect_size(eff.type="corr", r=0.08, n=332)
Effectsizes <- effect_size(eff.type="corr", r=-0.03, n=332)

# Jakobsson & Blom (2014)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=25.146, m.2=27.189, 
                           sd.1=5.36, sd.2=4.30, 
                           n.1=824, n.2=164)

# Johnson et al. (2009)
# a) Jewish sample
Effectsizes <- effect_size(eff.type="corr", r=0.002, n=560)
Effectsizes <- effect_size(eff.type="corr", r=0.039, n=560)
Effectsizes <- effect_size(eff.type="corr", r=0.013, n=560)
Effectsizes <- effect_size(eff.type="corr", r=0.038, n=560)
Effectsizes <- effect_size(eff.type="corr", r=0.023, n=560)
Effectsizes <- effect_size(eff.type="corr", r=0.040, n=560)
Effectsizes <- effect_size(eff.type="corr", r=0.039, n=560)
Effectsizes <- effect_size(eff.type="corr", r=0.140, n=560)
Effectsizes <- effect_size(eff.type="corr", r=0.096, n=560)
Effectsizes <- effect_size(eff.type="corr", r=0.021, n=560)
Effectsizes <- effect_size(eff.type="corr", r=0.009, n=560)
Effectsizes <- effect_size(eff.type="corr", r=0.004, n=560)
Effectsizes <- effect_size(eff.type="corr", r=0.116, n=560)
Effectsizes <- effect_size(eff.type="corr", r=0.135, n=560)
Effectsizes <- effect_size(eff.type="corr", r=0.137, n=560)
Effectsizes <- effect_size(eff.type="corr", r=0.111, n=560)
Effectsizes <- effect_size(eff.type="corr", r=0.058, n=560)
Effectsizes <- effect_size(eff.type="corr", r=0.089, n=560)
Effectsizes <- effect_size(eff.type="corr", r=0.194, n=560)
Effectsizes <- effect_size(eff.type="corr", r=0.172, n=560)
Effectsizes <- effect_size(eff.type="corr", r=0.193, n=560)
Effectsizes <- effect_size(eff.type="corr", r=0.150, n=560)
Effectsizes <- effect_size(eff.type="corr", r=0.047, n=560)
Effectsizes <- effect_size(eff.type="corr", r=0.044, n=560)
# a) Arab sample
Effectsizes <- effect_size(eff.type="corr", r=0.030, n=182)
Effectsizes <- effect_size(eff.type="corr", r=0.016, n=182)
Effectsizes <- effect_size(eff.type="corr", r=0.041, n=182)
Effectsizes <- effect_size(eff.type="corr", r=-0.040, n=182)
Effectsizes <- effect_size(eff.type="corr", r=-0.031, n=182)
Effectsizes <- effect_size(eff.type="corr", r=0.036, n=182)
Effectsizes <- effect_size(eff.type="corr", r=0.038, n=182)
Effectsizes <- effect_size(eff.type="corr", r=-0.005, n=182)
Effectsizes <- effect_size(eff.type="corr", r=0.039, n=182)
Effectsizes <- effect_size(eff.type="corr", r=0.159, n=182)
Effectsizes <- effect_size(eff.type="corr", r=0.007, n=182)
Effectsizes <- effect_size(eff.type="corr", r=0.119, n=182)
Effectsizes <- effect_size(eff.type="corr", r=0.183, n=182)
Effectsizes <- effect_size(eff.type="corr", r=0.095, n=182)
Effectsizes <- effect_size(eff.type="corr", r=0.205, n=182)
Effectsizes <- effect_size(eff.type="corr", r=0.123, n=182)
Effectsizes <- effect_size(eff.type="corr", r=-0.051, n=182)
Effectsizes <- effect_size(eff.type="corr", r=-0.013, n=182)
Effectsizes <- effect_size(eff.type="corr", r=0.160, n=182)
Effectsizes <- effect_size(eff.type="corr", r=0.148, n=182)
Effectsizes <- effect_size(eff.type="corr", r=0.200, n=182)
Effectsizes <- effect_size(eff.type="corr", r=0.062, n=182)
Effectsizes <- effect_size(eff.type="corr", r=0.021, n=182)
Effectsizes <- effect_size(eff.type="corr", r=0.029, n=182)

# Johnson et al. (2011)
Effectsizes <- effect_size(eff.type="chi.sq", 
                           chi.sq=0.494, n =(302+270))  #general
Effectsizes <- effect_size(eff.type="chi.sq", 
                           chi.sq=2.67, n =(302+271)) #stop and question
Effectsizes <- effect_size(eff.type="chi.sq", 
                           chi.sq=4.53, n =(300+269)) #bag searches
Effectsizes <- effect_size(eff.type="chi.sq", 
                           chi.sq=2.46, n =(301+271)) #wiretap
Effectsizes <- effect_size(eff.type="chi.sq", 
                           chi.sq=1.39, n =(301+271)) #house searches
Effectsizes[967:970,1] <- -Effectsizes[967:970,1]
Effectsizes[967:970,4] <- -Effectsizes[967:970,4]

# Jonathan (2010)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.77, m.2=3.26, 
                           sd.1=1.06, sd.2=1.03, 
                           n.1=1197, n.2=812)

# Jost at al. (2007)
Effectsizes <- effect_size(eff.type="corr", r=0.340, n=161)
Effectsizes <- effect_size(eff.type="corr", r=0.390, n=182)

# Jungkunz, Helbling, & Schwemmer (2019)
# a) Feeling thermometer
Effectsizes <- effect_size(eff.type="means", 
                           m.1=74.42857, m.2=70.47368, 
                           sd.1=16.56509, sd.2=23.05779,
                           n.1=28, n.2=38)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=71.89189, m.2=78.60417, 
                           sd.1=21.53006, sd.2=16.38076,
                           n.1=37, n.2=48)
# b) Social welfare benefits
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.607143, m.2=4.526316, 
                           sd.1=0.5669467, sd.2=0.9791569,
                           n.1=28, n.2=38)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.459459, m.2=4.680851, 
                           sd.1=0.7671953, sd.2=0.7831481,
                           n.1=37, n.2=47)
# c) Voting rights
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.789474, m.2=3.571429, 
                           sd.1=1.454869, sd.2=1.199647,
                           n.1=38, n.2=28)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.978723, m.2=3.243243, 
                           sd.1=1.326803, sd.2=1.233905,
                           n.1=47, n.2=37)
# e) Public office
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.000000, m.2=4.105263, 
                           sd.1=1.054093, sd.2=1.247473,
                           n.1=28, n.2=38)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.72973, m.2=4.212766, 
                           sd.1=1.170201, sd.2=1.082191,
                           n.1=37, n.2=47)

# Kaakinen et al. (2011)
Effectsizes <- effect_size(eff.type="corr", r=0.268, n=(555+192)) #general
Effectsizes <- effect_size(eff.type="corr", r=0.282, n=(555+192)) #ethnicity/nationality
Effectsizes <- effect_size(eff.type="corr", r=0.341, n=(555+192)) #religion
Effectsizes <- effect_size(eff.type="corr", r=0.191, n=(555+192)) #political view
Effectsizes <- effect_size(eff.type="corr", r=-0.073, n=(555+192)) #sexual orientation
Effectsizes <- effect_size(eff.type="p.val", p=0.041, n.1=555, n.2=192, tail="two") #gender
Effectsizes <- effect_size(eff.type="p.val", p=0.780, n.1=555, n.2=192, tail="two") #disability 

# Kalagy et al. (2007)
# a) National Religious sample (n=53)
Effectsizes <- effect_size(eff.type="corr", r=-0.36, n=53) #anxiety-peace
Effectsizes <- effect_size(eff.type="corr", r=-0.13, n=53) #anger-peace
Effectsizes <- effect_size(eff.type="corr", r=-0.16, n=53) #anxiety-war
Effectsizes <- effect_size(eff.type="corr", r=0.12, n=53) #anger-war
# b) Ultra-Orthodox sample (n=53)
Effectsizes <- effect_size(eff.type="corr", r=-0.03, n=54) #anxiety-peace
Effectsizes <- effect_size(eff.type="corr", r=-0.03, n=54) #anger-peace
Effectsizes <- effect_size(eff.type="corr", r=-0.02, n=54) #anxiety-war
Effectsizes <- effect_size(eff.type="corr", r=-0.02, n=54) #anger-war

# Kastenmueller (2011)
# a) Study 1
Effectsizes <- effect_size(eff.type="means", 
                           m.1=.65, m.2=1.24, 
                           sd.1=1.01, sd.2=1.21, 
                           n.1=20, n.2=19)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=1.71, m.2=1.34, 
                           sd.1=1.17, sd.2=1.25, 
                           n.1=19, n.2=19)
# b) Study 2
Effectsizes <- effect_size(eff.type="means", 
                           m.1=.96, m.2=.88, 
                           sd.1=.85, sd.2=.86, 
                           n.1=24, n.2=25)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=1.44, m.2=1.00, 
                           sd.1=.98, sd.2=.81, 
                           n.1=24, n.2=23)

# Kim (2016)
# a) Study 1
Effectsizes <- effect_size(eff.type="corr", r=0.220, n=394)
Effectsizes <- effect_size(eff.type="corr", r=0.340, n=394)
Effectsizes <- effect_size(eff.type="corr", r=0.050, n=394)
Effectsizes <- effect_size(eff.type="corr", r=0.210, n=394)
Effectsizes <- effect_size(eff.type="corr", r=0.130, n=394)
Effectsizes <- effect_size(eff.type="corr", r=0.290, n=394)
Effectsizes <- effect_size(eff.type="corr", r=0.060, n=394)
Effectsizes <- effect_size(eff.type="corr", r=0.200, n=394)
Effectsizes <- effect_size(eff.type="corr", r=0.190, n=394)
Effectsizes <- effect_size(eff.type="corr", r=0.270, n=394)
Effectsizes <- effect_size(eff.type="corr", r=0.100, n=394)
Effectsizes <- effect_size(eff.type="corr", r=0.270, n=394)
# b) Study 2
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.41, m.2=3.39, 
                           sd.1=3.36, sd.2=1.39, 
                           n.1=40, n.2=40)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=2.71, m.2=2.14, 
                           sd.1=1.45, sd.2=1.35, 
                           n.1=40, n.2=40)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=1.88, m.2=1.40, 
                           sd.1=0.91, sd.2=0.63, 
                           n.1=40, n.2=40)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.16, m.2=2.57, 
                           sd.1=1.55, sd.2=1.52, 
                           n.1=40, n.2=40)

# Kimhi et al. (2019)
# a) IV: Sense of danger
Effectsizes <- effect_size(eff.type="corr", r=-0.229, n=1022) #Overall:national resilience
Effectsizes <- effect_size(eff.type="corr", r=-0.178, n=1022) #1:Identification with the state
Effectsizes <- effect_size(eff.type="corr", r=0.203, n=1022) #3:Trust in national institutions
Effectsizes <- effect_size(eff.type="corr", r=-0.128, n=1022) #4:Trust in public justice
# b) IV: Exposure
Effectsizes <- effect_size(eff.type="corr", r=-0.064, n=1022) #Overall:national resilience
Effectsizes <- effect_size(eff.type="corr", r=-0.087, n=1022) #1:Identification with the state
Effectsizes <- effect_size(eff.type="corr", r=-0.030, n=1022) #3:Trust in national institutions
Effectsizes <- effect_size(eff.type="corr", r=-0.095, n=1022) #4:Trust in public justice

# Klar & Schori-Eyal (2015)
# a) Study 1
Effectsizes <- effect_size(eff.type="corr", r=0.05, n=80)
Effectsizes <- effect_size(eff.type="corr", r=0.03, n=80)
Effectsizes <- effect_size(eff.type="corr", r=0.13, n=80)
Effectsizes <- effect_size(eff.type="corr", r=0.22, n=80)
# b) Study 2
Effectsizes <- effect_size(eff.type="corr", r=-0.03, n=123)
Effectsizes <- effect_size(eff.type="corr", r=0.14, n=123)
Effectsizes <- effect_size(eff.type="corr", r=0.06, n=123)
Effectsizes <- effect_size(eff.type="corr", r=0.08, n=123)
Effectsizes <- effect_size(eff.type="corr", r=0.23, n=123)

# Kossowska et al. (2011)
Effectsizes <- effect_size(eff.type="corr", r=0.377, n=86) 
Effectsizes <- effect_size(eff.type="corr", r=0.573, n=86) 
Effectsizes <- effect_size(eff.type="corr", r=0.210, n=68)
Effectsizes <- effect_size(eff.type="corr", r=0.313, n=68)
Effectsizes <- effect_size(eff.type="corr", r=0.256, n=121) 
Effectsizes <- effect_size(eff.type="corr", r=0.174, n=121) 
Effectsizes <- effect_size(eff.type="corr", r=0.350, n=175)
Effectsizes <- effect_size(eff.type="corr", r=0.349, n=175)

# Kossowska, de Zavala, & Kubik (2010)
Effectsizes <- effect_size(eff.type="corr", r= 0.35, n=52) 
Effectsizes <- effect_size(eff.type="corr", r= 0.23, n=52) 
Effectsizes <- effect_size(eff.type="corr", r=-0.21, n=52)
Effectsizes <- effect_size(eff.type="p.val", p=0.99, n.1=(52/2), n.2=(52/2), tail="two")

# Kteily et al. (2015)
# Study 3A
Effectsizes <- effect_size(eff.type="means", 
                           m.1=15.58, m.2=10.77, 
                           sd.1=25.53, sd.2=23.44, 
                           n.1=345, n.2=420)#3 days after attack
Effectsizes <- effect_size(eff.type="means", 
                           m.1=10.14, m.2=10.77, 
                           sd.1=21.99, sd.2=23.44, 
                           n.1=204, n.2=420)#6 months after attack
Effectsizes <- effect_size(eff.type="corr", r=0.36, n=345)

# Ladd (2007)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=66.8, m.2=57.9, 
                           sd.1=27.3, sd.2=27.3, 
                           n.1=1071, n.2=1071)

# Lahav & Courtemanche (2012)
# a) overall sample
Effectsizes <- effect_size(eff.type="means", 
                           m.1=0.584, m.2=0.551, 
                           sd.1=0.27, sd.2=0.26,
                           n.1=99, n.2=99)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=0.358, m.2=0.300, 
                           sd.1=0.24, sd.2=0.23,
                           n.1=99, n.2=99)
# b) liberal sample
Effectsizes <- effect_size(eff.type="means", 
                           m.1=0.523, m.2=0.511, 
                           sd.1=0.25, sd.2=0.27,
                           n.1=57, n.2=57)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=.279, m.2=.238, 
                           sd.1=0.21, sd.2=0.21,
                           n.1=57, n.2=57)
# c) conservative sample
Effectsizes <- effect_size(eff.type="means", 
                           m.1=.660, m.2=.583, 
                           sd.1=0.28, sd.2=0.24,
                           n.1=20, n.2=20)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=0.449, m.2=0.442, 
                           sd.1=0.26, sd.2=0.24,
                           n.1=20, n.2=20)
# d) moderate sample
Effectsizes <- effect_size(eff.type="means", 
                           m.1=.645, m.2=.619, 
                           sd.1=0.28, sd.2=0.23,
                           n.1=22, n.2=22)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=.450, m.2=.373, 
                           sd.1=0.18, sd.2=0.24,
                           n.1=22, n.2=22)

# Lahav, Shahrabani, & Benzion (2019)
# a) Anger at Hamas
Effectsizes <- effect_size(eff.type="corr", r=0.10, n=515)
Effectsizes <- effect_size(eff.type="corr", r=0.06, n=515)
Effectsizes <- effect_size(eff.type="corr", r=0.08, n=515)
# b) Anger at Government
Effectsizes <- effect_size(eff.type="corr", r=-0.06, n=515)
Effectsizes <- effect_size(eff.type="corr", r=-0.09, n=515)
Effectsizes <- effect_size(eff.type="corr", r=-0.29, n=515)
# c) Perceived risk
Effectsizes <- effect_size(eff.type="corr", r=0.01, n=515)
Effectsizes <- effect_size(eff.type="corr", r=0.02, n=515)
Effectsizes <- effect_size(eff.type="corr", r=0.03, n=515)

# Lambert et al. (2010)
# a) Study 1
Effectsizes <- effect_size(eff.type="f.test", f=11.76, n.1=66, n.2=66) 
Effectsizes <- effect_size(eff.type="f.test", f=6.07, n.1=66, n.2=66)    
Effectsizes <- effect_size(eff.type="f.test", f=4.60, n.1=66, n.2=66)    
Effectsizes <- effect_size(eff.type="f.test", f=0.01, n.1=66, n.2=66)
# b) Study 2
Effectsizes <- effect_size(eff.type="f.test", f=4.96, n.1=23, n.2=23) 
Effectsizes <- effect_size(eff.type="f.test", f=0.01, n.1=23, n.2=23)
Effectsizes <- effect_size(eff.type="f.test", f=5.36, n.1=23, n.2=23) 
# c) Study 3
Effectsizes <- effect_size(eff.type="f.test", f=3.91, n.1=25, n.2=25) 
Effectsizes <- effect_size(eff.type="f.test", f=4.42, n.1=25, n.2=25) 
Effectsizes <- effect_size(eff.type="f.test", f=0.01, n.1=25, n.2=25)

# Landau et al. (2004)
Effectsizes <- effect_size(eff.type="t.test", t=9.31, n.1=25, n.2=49) 
Effectsizes <- effect_size(eff.type="f.test", f=0.01, n.1=25, n.2=49)

# Larsen, Cuttes, & Goodwin (2019)
# Larsenetal2019_ch <- read.csv("https://www.dropbox.com/s/8mrme3a1nxx90a3/Larsenetal2019_ch.csv?dl=1")
# Larsenetal2019_ess <- read.csv("https://www.dropbox.com/s/4a6olkpoodz9nzi/Larsenetal2019_ess.csv?dl=1")
# article focuses on germany:
# ch_germany <- Larsenetal2019_ch %>%
#   filter(germany == 1)
# describe(ch_germany$age)
# describe(ch_germany$male)
# describe(ch_germany$lrscale)
# ess_germany <- Larsenetal2019_ess %>%
#   filter(germany == 1)
# describe(ess_germany$age)
# describe(ess_germany$male)
# describe(ess_germany$lrscale)
# a) The effect of the Berlin attack on attitudes towards the EU in Germany, bivariate OLS regressions
# cor.test(ch_germany$tr, ch_germany$eu1)
Effectsizes <- effect_size(eff.type="corr", r=-0.09332056, n=817)
# cor.test(ch_germany$tr, ch_germany$eu2)
Effectsizes <- effect_size(eff.type="corr", r=-0.07990301, n=817)
# cor.test(ch_germany$tr, ch_germany$eu3)
Effectsizes <- effect_size(eff.type="corr", r=-0.09810114, n=816)
# cor.test(ch_germany$tr, ch_germany$eu4)
Effectsizes <- effect_size(eff.type="corr", r=-0.12176, n=816)
# cor.test(ess_germany$tr, ess_germany$eu5)
Effectsizes <- effect_size(eff.type="corr", r=-0.1256395, n=462) #ESS
# b) The effect of the Berlin attack on attitudes towards immigration in Germany, bivariate OLS regressions
# cor.test(ch_germany$tr, ch_germany$imm1)
Effectsizes <- effect_size(eff.type="corr", r=0.01404887, n=817)
# cor.test(ch_germany$tr, ch_germany$imm2)
Effectsizes <- effect_size(eff.type="corr", r=0.04762521, n=817)
# cor.test(ch_germany$tr, ch_germany$imm3)
Effectsizes <- effect_size(eff.type="corr", r=0.006415198, n=817)
# cor.test(ch_germany$tr, ch_germany$imm4)
Effectsizes <- effect_size(eff.type="corr", r=0.05272005, n=817)
# cor.test(ch_germany$tr, ch_germany$imm5)
Effectsizes <- effect_size(eff.type="corr", r=0.07036447, n=817)
# cor.test(ch_germany$tr, ch_germany$imm6)
Effectsizes <- effect_size(eff.type="corr", r=0.004653309, n=816)
# cor.test(ess_germany$tr, ess_germany$imm7)
Effectsizes <- effect_size(eff.type="corr", r=0.05834226, n=468) #ESS
# c) The effect of the Berlin attack on attitudes towards refugees, bivariate OLS regressions in Germany, bivariate OLS regressions
# cor.test(ch_germany$tr, ch_germany$ref1)
Effectsizes <- effect_size(eff.type="corr", r=-0.003984094, n=817)
# cor.test(ess_germany$tr, ess_germany$ref2)
Effectsizes <- effect_size(eff.type="corr", r=0.0319909, n=471) #ESS
# cor.test(ess_germany$tr, ess_germany$ref3)
Effectsizes <- effect_size(eff.type="corr", r=0.1045562, n=456) #ESS
# cor.test(ess_germany$tr, ess_germany$ref4)
Effectsizes <- effect_size(eff.type="corr", r=0.02608492, n=470) #ESS

# Laufer & Solomon (2010)
Effectsizes <- effect_size(eff.type="corr", r=0.35, n=2525) 
Effectsizes <- effect_size(eff.type="corr", r=0.08, n=2525) 
Effectsizes <- effect_size(eff.type="corr", r=0.34, n=1772)
Effectsizes <- effect_size(eff.type="corr", r=0.08, n=1772)
Effectsizes <- effect_size(eff.type="corr", r=0.25, n=569) 
Effectsizes <- effect_size(eff.type="corr", r=0.21, n=569) 
Effectsizes <- effect_size(eff.type="corr", r=0.11, n=184)
Effectsizes <- effect_size(eff.type="corr", r=-0.03, n=184)

# Lavi & Solomon (2005)
# a) Israeli-Palestinian sample
Effectsizes <- effect_size(eff.type="corr", r=0.051, n=300)
Effectsizes <- effect_size(eff.type="corr", r=0.070, n=300)
Effectsizes <- effect_size(eff.type="corr", r=0.100, n=300)
Effectsizes <- effect_size(eff.type="corr", r=0.195, n=300)
Effectsizes <- effect_size(eff.type="corr", r=0.223, n=300)
Effectsizes <- effect_size(eff.type="corr", r=0.073, n=300)
Effectsizes <- effect_size(eff.type="corr", r=0.099, n=300)
# b) PA sample
Effectsizes <- effect_size(eff.type="corr", r=0.330, n=245)
Effectsizes <- effect_size(eff.type="corr", r=0.107, n=245)
Effectsizes <- effect_size(eff.type="corr", r=0.143, n=245)
Effectsizes <- effect_size(eff.type="corr", r=0.041, n=245)
Effectsizes <- effect_size(eff.type="corr", r=0.150, n=245)
Effectsizes <- effect_size(eff.type="corr", r=0.144, n=245)
Effectsizes <- effect_size(eff.type="corr", r=0.158, n=245)

# Lavi et al. (2014)
# a)  Israeli sample
##political attitudes
Effectsizes <- effect_size(eff.type="corr", r=0.08, n=781) #General exposure
Effectsizes <- effect_size(eff.type="corr", r=0.09, n=781) #Financial loss
Effectsizes <- effect_size(eff.type="corr", r=0.01, n=781) #Property damage
Effectsizes <- effect_size(eff.type="corr", r=0.11, n=781) #PSS terror
Effectsizes <- effect_size(eff.type="corr", r=0.06, n=781) #National threat
Effectsizes <- effect_size(eff.type="corr", r=0.07, n=781) #Personal threat
##fear Palestinians
Effectsizes <- effect_size(eff.type="corr", r=0.06, n=781) #General exposure
Effectsizes <- effect_size(eff.type="corr", r=0.06, n=781) #Financial loss
Effectsizes <- effect_size(eff.type="corr", r=0.03, n=781) #Property damage
Effectsizes <- effect_size(eff.type="corr", r=0.12, n=781) #PSS terror
Effectsizes <- effect_size(eff.type="corr", r=-0.03, n=781) #National threat
Effectsizes <- effect_size(eff.type="corr", r=0.08, n=781) #Personal threat
##hate Palestinians
Effectsizes <- effect_size(eff.type="corr", r=0.08, n=781) #General exposure
Effectsizes <- effect_size(eff.type="corr", r=0.05, n=781) #Financial loss
Effectsizes <- effect_size(eff.type="corr", r=0.03, n=781) #Property damage
Effectsizes <- effect_size(eff.type="corr", r=0.03, n=781) #PSS terror
Effectsizes <- effect_size(eff.type="corr", r=0.01, n=781) #National threat
Effectsizes <- effect_size(eff.type="corr", r=0.03, n=781) #Personal threat
# b)  Arab sample
##political attitudes
Effectsizes <- effect_size(eff.type="corr", r=-0.02, n=1196) #Death of friend
Effectsizes <- effect_size(eff.type="corr", r=-0.01, n=1196) #Injury to self
Effectsizes <- effect_size(eff.type="corr", r=0.04, n=1196) #Injury to friend
Effectsizes <- effect_size(eff.type="corr", r=0.08, n=1196) #Witnessing
Effectsizes <- effect_size(eff.type="corr", r=-0.04, n=1196) #Home demolished
Effectsizes <- effect_size(eff.type="corr", r=-0.08, n=1196) #Financial loss
Effectsizes <- effect_size(eff.type="corr", r=-0.08, n=1196) #PSS terror
Effectsizes <- effect_size(eff.type="corr", r=0.00, n=1196) #National threat
Effectsizes <- effect_size(eff.type="corr", r=-0.05, n=1196) #Personal threat
##fear Jews
Effectsizes <- effect_size(eff.type="corr", r=0.02, n=1196) #Death of friend
Effectsizes <- effect_size(eff.type="corr", r=-0.02, n=1196) #Injury to self
Effectsizes <- effect_size(eff.type="corr", r=0.00, n=1196) #Injury to friend
Effectsizes <- effect_size(eff.type="corr", r=-0.06, n=1196) #Witnessing
Effectsizes <- effect_size(eff.type="corr", r=0.09, n=1196) #Home demolished
Effectsizes <- effect_size(eff.type="corr", r=-0.01, n=1196) #Financial loss
Effectsizes <- effect_size(eff.type="corr", r=0.19, n=1196) #PSS terror
Effectsizes <- effect_size(eff.type="corr", r=0.20, n=1196) #National threat
Effectsizes <- effect_size(eff.type="corr", r=0.15, n=1196) #Personal threat
##hate Jews
Effectsizes <- effect_size(eff.type="corr", r=0.03, n=1196) #Death of friend
Effectsizes <- effect_size(eff.type="corr", r=0.01, n=1196) #Injury to self
Effectsizes <- effect_size(eff.type="corr", r=0.05, n=1196) #Injury to friend
Effectsizes <- effect_size(eff.type="corr", r=0.08, n=1196) #Witnessing
Effectsizes <- effect_size(eff.type="corr", r=0.04, n=1196) #Home demolished
Effectsizes <- effect_size(eff.type="corr", r=0.02, n=1196) #Financial loss
Effectsizes <- effect_size(eff.type="corr", r=0.20, n=1196) #PSS terror
Effectsizes <- effect_size(eff.type="corr", r=0.22, n=1196) #National threat
Effectsizes <- effect_size(eff.type="corr", r=0.18, n=1196) #Personal threat

# Legewie (2013): information at country level via mail correspondence
Effectsizes <- effect_size(eff.type="means", 
                           m.1=0.319, m.2=0.250, 
                           sd.1=0.892, sd.2=0.814,
                           n.1=205, n.2=277) # Belgium
Effectsizes <- effect_size(eff.type="means", 
                           m.1=-0.0754, m.2=-0.125, 
                           sd.1=0.798, sd.2=0.785,
                           n.1=57, n.2=312) # Switserland
Effectsizes <- effect_size(eff.type="means", 
                           m.1=-0.0458, m.2=-0.201, 
                           sd.1=0.851, sd.2=0.831,
                           n.1=136, n.2=711) # Finland	
Effectsizes <- effect_size(eff.type="means", 
                           m.1=0.312, m.2=0.435, 
                           sd.1=0.925, sd.2=1.08,
                           n.1=304, n.2=333) # UK
Effectsizes <- effect_size(eff.type="means", 
                           m.1=0.0297, m.2=0.0855, 
                           sd.1=0.780, sd.2=0.806,
                           n.1=133, n.2=514) # Netherlands
Effectsizes <- effect_size(eff.type="means", 
                           m.1=-0.114, m.2=-0.0411, 
                           sd.1=0.753, sd.2=0.838,
                           n.1=163, n.2=796) # Norway					
Effectsizes <- effect_size(eff.type="means", 
                           m.1=-0.148, m.2=0.0787, 
                           sd.1=0.912, sd.2=1.01,
                           n.1=316, n.2=67) # Poland						
Effectsizes <- effect_size(eff.type="means", 
                           m.1=0.572, m.2=0.113, 
                           sd.1=0.884, sd.2=0.932,
                           n.1=69, n.2=137) # Portugal					
Effectsizes <- effect_size(eff.type="means", 
                           m.1=-0.352, m.2=-0.367, 
                           sd.1=0.962, sd.2=0.865,
                           n.1=187, n.2=120) # Sweden		

# Lerner et al. (2003)
Effectsizes <- effect_size(eff.type="corr", r=-0.06, n=973)
Effectsizes <- effect_size(eff.type="corr", r=0.28, n=973)
Effectsizes <- effect_size(eff.type="corr", r=0.02, n=973)
Effectsizes <- effect_size(eff.type="corr", r=0.13, n=973)

# Lett, DiPietro, & Johnson (2004)
Effectsizes <- effect_size(eff.type="corr", r=beta2cor(0.21,1), n=234)

# Liberman & Skitka (2019)
Effectsizes <- effect_size(eff.type="corr", r=0.29, n=595)
Effectsizes <- effect_size(eff.type="corr", r=0.06, n=595)
Effectsizes <- effect_size(eff.type="corr", r=0.22, n=595)
Effectsizes <- effect_size(eff.type="corr", r=0.25, n=595)
Effectsizes <- effect_size(eff.type="corr", r=0.23, n=595)
Effectsizes <- effect_size(eff.type="corr", r=0.19, n=595)
Effectsizes <- effect_size(eff.type="corr", r=0.27, n=595)
Effectsizes <- effect_size(eff.type="corr", r=0.03, n=595)
Effectsizes <- effect_size(eff.type="corr", r=0.22, n=595)
Effectsizes <- effect_size(eff.type="corr", r=.2418, n=595)
Effectsizes <- effect_size(eff.type="corr", r=.0265, n=595)
Effectsizes <- effect_size(eff.type="corr", r=.0966, n=595)

# Liendo & Braitwaite (2018)
Effectsizes <- effect_size(eff.type="corr", r= OR2cor(exp(-0.018)), n=1136) 
Effectsizes <- effect_size(eff.type="corr", r= OR2cor(exp(0.014)), n=1136) 
Effectsizes <- effect_size(eff.type="corr", r= OR2cor(exp(-0.056)), n=1136) 

# Linden, Bjorklund, & Backstrom (2018)
Effectsizes <- effect_size(eff.type="coh.d", d=0.31, n.1=152, n.2=140)
Effectsizes <- effect_size(eff.type="coh.d", d=0.13, n.1=152, n.2=140)
Effectsizes <- effect_size(eff.type="coh.d", d=0.06, n.1=152, n.2=140)
Effectsizes <- effect_size(eff.type="coh.d", d=0.42, n.1=152, n.2=140)

# Lizotte (2017)
# Load data and subset female/male dataframe
# Lizotte2017 <- read_dta("https://www.dropbox.com/s/djdetf66y6n0bkx/Lizotte2017.dta?dl=1")
# Lizotte2017_female <- subset(Lizotte2017, female==1, 
#                              select=c(threat, protorturedich, age, white, pid))
# Lizotte2017_male <- subset(Lizotte2017, female==0, 
#                            select=c(threat, protorturedich, age, white, pid))
# # Descriptives (covariates)
# summary(factor(Lizotte2017$protorturedich))
# summary(Lizotte2017$female) 
# summary(Lizotte2017$age) 
# sqrt(var(Lizotte2017$age))
# summary(Lizotte2017$pid) 
# summary(Lizotte2017_female$age) 
# sqrt(var(Lizotte2017_female$age))
# summary(Lizotte2017_female$pid) 
# summary(Lizotte2017_female$white) 
# summary(Lizotte2017_male$age)
# sqrt(var(Lizotte2017_male$age))
# summary(Lizotte2017_male$pid) 
# summary(Lizotte2017_male$white) 
# # Biserial correlation between threat perceptions and pro-torture attitudes
# biserial.cor(Lizotte2017$threat, (factor(Lizotte2017$protorturedich)), 
#              use="complete.obs", level=1) #all respondents
# biserial.cor(Lizotte2017_female$threat, Lizotte2017_female$protorturedich, 
#              use="complete.obs", level=1) #female respondents
# biserial.cor(Lizotte2017_male$threat, Lizotte2017_male$protorturedich, 
#              use="complete.obs", level=1) #male respondents
# Effect size calculation
Effectsizes <- effect_size(eff.type="corr", r=0.0729115, n=2009) #all respondents
Effectsizes <- effect_size(eff.type="corr", r=0.0237114, n=1136) #female respondents
Effectsizes <- effect_size(eff.type="corr", r=0.1283340, n=873) #male respondents

# Lopes & Jaspal (2015)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=61.37, m.2=56.35, 
                           sd.1=40.5, sd.2=36.8, 
                           n.1=30, n.2=23) 
Effectsizes <- effect_size(eff.type="means", 
                           m.1=28.37, m.2=29.35, 
                           sd.1=12.3, sd.2=12.9,
                           n.1=30, n.2=23) 

# Lupu & Peisakhin (2017)
Effectsizes <- effect_size(eff.type="corr", r=sqrt(0.010), n=974) # Turnout
Effectsizes <- effect_size(eff.type="corr", r=sqrt(0.015), n=727) # Participation

# Maeseele et al. (2008)
Effectsizes <- effect_size(eff.type="corr", r=-0.139, n=1040)
Effectsizes <- effect_size(eff.type="corr", r=-0.039, n=1040)
Effectsizes <- effect_size(eff.type="corr", r=-0.143, n=1040)
Effectsizes <- effect_size(eff.type="corr", r=-0.095, n=1040)
Effectsizes <- effect_size(eff.type="corr", r=-0.127, n=1040)
Effectsizes <- effect_size(eff.type="corr", r=-0.149, n=1040)

# Maitner, Mackie, & Smith (2006)
Effectsizes <- effect_size(eff.type="corr", r=beta2cor(.504,1), n=318) #anger
Effectsizes <- effect_size(eff.type="corr", r=beta2cor(-.216,0), n=318) #fear

# Malhotra & Popp (2012)
Effectsizes <- effect_size(eff.type="corr", r=0.23, n=227) #corr in control group
Effectsizes <- effect_size(eff.type="corr", r=beta2cor(0.08,1), n=646) #exp in overall sample

# Maoz & Eidelson (2007)
Effectsizes <- effect_size(eff.type="corr", r=0.30, n=504) # Transfer vs. Compromise
Effectsizes <- effect_size(eff.type="corr", r=0.21, n=504) # Hawk 

# Maoz & McCauley (2008)
# a) study 1: threat from Palestinian terrorism
Effectsizes <- effect_size(eff.type="corr", r=0.48, n=504) # Transfer
Effectsizes <- effect_size(eff.type="corr", r=0.37, n=504) # Dehumanization 
Effectsizes <- effect_size(eff.type="corr", r=0.51, n=504) # Hawk 
# b) study 1: threat from Palestinian terrorism
Effectsizes <- effect_size(eff.type="corr", r=0.42, n=501) # CCA
Effectsizes <- effect_size(eff.type="corr", r=0.34, n=501) # Dehumanization 
Effectsizes <- effect_size(eff.type="corr", r=0.54, n=501) # Hawk 

# Maoz & McCauley (2009)
Effectsizes <- effect_size(eff.type="corr", r=0.49, n=504)
Effectsizes <- effect_size(eff.type="corr", r=0.10, n=504)
Effectsizes <- effect_size(eff.type="corr", r=0.49, n=504)
Effectsizes <- effect_size(eff.type="corr", r=0.27, n=504)
Effectsizes <- effect_size(eff.type="corr", r=0.54, n=504)
Effectsizes <- effect_size(eff.type="corr", r=0.04, n=504)

# Matthes, Schmunck & von Sikorski (2019)
Effectsizes <- effect_size(eff.type="corr", r=0.48, n=501)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.12233009708738, m.2=4.28979591836735, 
                           sd.1=1.62341518730108, sd.2=1.55108806478724, 
                           n.1=103, n.2=98)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.42526315789474, m.2=4.28979591836735, 
                           sd.1=1.55529079315985, sd.2=1.55108806478724, 
                           n.1=95, n.2=98)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.94038461538461, m.2=4.28979591836735, 
                           sd.1=1.74308727512223, sd.2=1.55108806478724, 
                           n.1=104, n.2=98)
Effectsizes <- effect_size(eff.type="means",
                           m.1=4.16435643564356, m.2=4.28979591836735, 
                           sd.1=1.46243524016729, sd.2=1.55108806478724, 
                           n.1=101, n.2=98)

#  Merolla, Ramos, & Zechmeister (2007)
Effectsizes <- b_effect_size(eff.type="beta", beta=0.184, sdy=1.12, grp1n=102, grp2n=94)

# Moaddel & Latif
# a) Egypt
Effectsizes <- effect_size(eff.type="corr", r=-0.225, n=2580)
Effectsizes <- effect_size(eff.type="corr", r=0.060, n=2580)
Effectsizes <- effect_size(eff.type="corr", r=-0.088, n=2580)
# b) Morocco
Effectsizes <- effect_size(eff.type="corr", r=-0.094, n=528)
Effectsizes <- effect_size(eff.type="corr", r=-0.090, n=528)
Effectsizes <- effect_size(eff.type="corr", r=-0.311, n=528)

# Mondak & Hurwitz
Effectsizes <- effect_size(eff.type="means", 
                           m.1=2.93, m.2=2.47, 
                           sd.1=1.17, sd.2=1.05, 
                           n.1=522, n.2=564) 

# Morales-Marente et al. (2009)
Effectsizes <- effect_size(eff.type="corr", r=0.03, n=326)
Effectsizes <- effect_size(eff.type="corr", r=0.04, n=326)
Effectsizes <- effect_size(eff.type="corr", r=0.03, n=326)
Effectsizes <- effect_size(eff.type="corr", r=-0.04, n=326)

# Mosher (2008)
# exp. 1
Effectsizes <- effect_size(eff.type="means", 
                           m.1=(6.97+5.58)/2, m.2=(6.66+5.86)/2, 
                           sd.1=1, sd.2=1, 
                           n.1=65, n.2=65) #use of force
Effectsizes <- effect_size(eff.type="f.test", f=4.57, n.1=45, n.2=45) #negotations
Effectsizes[1244,1] <- -Effectsizes[1244,1]
Effectsizes[1244,4] <- -Effectsizes[1244,4]
# exp. 2
Effectsizes <- effect_size(eff.type="f.test", f=2.72, n.1=66, n.2=66) #use of force
Effectsizes[1245,1] <- -Effectsizes[1245,1]
Effectsizes[1245,4] <- -Effectsizes[1245,4]
Effectsizes <- effect_size(eff.type="f.test", f=3.61, n.1=65, n.2=65) #negotations
Effectsizes[1246,1] <- -Effectsizes[1246,1]
Effectsizes[1246,4] <- -Effectsizes[1246,4]

# Moskalenko et al. (2006)
# a) 9-item Country Identification Scale
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.72, m.2=4.37, 
                           sd.1=0.93, sd.2=1.04, 
                           n.1=159, n.2=99)#T1
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.26, m.2=4.37, 
                           sd.1=1.13, sd.2=1.04, 
                           n.1=351, n.2=99)#T2
# b) Importance country
Effectsizes <- effect_size(eff.type="means", 
                           m.1=5.6, m.2=4.8, 
                           sd.1=1.2, sd.2=1.7, 
                           n.1=159, n.2=99)#T1
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.9, m.2=4.8, 
                           sd.1=1.6, sd.2=1.7, 
                           n.1=351, n.2=99)#T2

# Nagoshi et al. (2007)
# a) RWA
Effectsizes <- effect_size(eff.type="f.test", f=6.27, n.1=86, n.2=91) #full sample
Effectsizes <- effect_size(eff.type="means", 
                           m.1=0.80, m.2=0.41, 
                           sd.1=0.72, sd.2=0.64, 
                           n.1=20, n.2=34)#male
Effectsizes <- effect_size(eff.type="means", 
                           m.1=0.55, m.2=0.33, 
                           sd.1=0.65, sd.2=0.76, 
                           n.1=66, n.2=57)#female
# b) SDO
Effectsizes <- effect_size(eff.type="means", 
                           m.1=(4.11+3.30)/2, m.2=(3.86+3.20)/2, 
                           sd.1=(0.63+0.65)/2, sd.2=(0.78+0.75)/2, 
                           n.1=(20+66), n.2=(34+57))#full sample
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.11, m.2=3.86, 
                           sd.1=0.63, sd.2=0.78, 
                           n.1=20, n.2=34)#male
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.30, m.2=3.20, 
                           sd.1=0.65, sd.2=0.75, 
                           n.1=66, n.2=57)#female

# Nail & McGregor (2009)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=1.75, m.2=0.81, 
                           sd.1=1.75, sd.2=1.81, 
                           n.1=70, n.2=69)
Effectsizes <- effect_size(eff.type="f.test", f=1.74, n.1=70, n.2=69)
Effectsizes <- effect_size(eff.type="f.test", f=79.95, n.1=70, n.2=69)
Effectsizes <- effect_size(eff.type="f.test", f=3.63, n.1=70, n.2=69)
Effectsizes <- effect_size(eff.type="f.test", f=32.51, n.1=70, n.2=69)
Effectsizes <- effect_size(eff.type="f.test", f=29.75, n.1=70, n.2=69)
Effectsizes <- effect_size(eff.type="f.test", f=0.20, n.1=70, n.2=69)
Effectsizes <- effect_size(eff.type="f.test", f=1.93, n.1=70, n.2=69)
Effectsizes <- effect_size(eff.type="f.test", f=3.13, n.1=70, n.2=69)
Effectsizes <- effect_size(eff.type="f.test", f=1.46, n.1=70, n.2=69)
Effectsizes <- effect_size(eff.type="f.test", f=1.16, n.1=70, n.2=69)


# Nilsen et al. (2019)
# a) Survivors
Effectsizes <- effect_size(eff.type="corr", r=-0.30, n=325)
Effectsizes <- effect_size(eff.type="corr", r=-0.20, n=325)
Effectsizes <- effect_size(eff.type="corr", r=-0.33, n=325)
Effectsizes <- effect_size(eff.type="corr", r=-0.27, n=325)
Effectsizes <- effect_size(eff.type="corr", r=-0.33, n=285)
Effectsizes <- effect_size(eff.type="corr", r=-0.29, n=285)
Effectsizes <- effect_size(eff.type="corr", r=-0.32, n=285)
Effectsizes <- effect_size(eff.type="corr", r=-0.29, n=285)
# b) Parents of survivors
Effectsizes <- effect_size(eff.type="corr", r=-0.09, n=453)
Effectsizes <- effect_size(eff.type="corr", r=-0.11, n=453)
Effectsizes <- effect_size(eff.type="corr", r=-0.09, n=453)
Effectsizes <- effect_size(eff.type="corr", r=-0.10, n=453)
Effectsizes <- effect_size(eff.type="corr", r=-0.24, n=435)
Effectsizes <- effect_size(eff.type="corr", r=-0.15, n=435)
Effectsizes <- effect_size(eff.type="corr", r=-0.21, n=435)
Effectsizes <- effect_size(eff.type="corr", r=-0.14, n=435)

# Niwa et al. (2016)
Effectsizes <- effect_size(eff.type="corr", r=0.21, n=600+901)
Effectsizes <- effect_size(eff.type="corr", r=0.15, n=600+901)
Effectsizes <- effect_size(eff.type="corr", r=0.25, n=600+901)

# Norris (2017)
Effectsizes <- effect_size(eff.type="t.test", t=-0.466, n.1=98, n.2 =129)

# Nugier et al. (2016)
# a) Overall threat
Effectsizes <- effect_size(eff.type="means", 
                           m.1=2.96, m.2=2.45, 
                           sd.1=0.61, sd.2=0.39, 
                           n.1=24, n.2=6)#Control
Effectsizes <- effect_size(eff.type="means", 
                           m.1=2.95, m.2=2.62, 
                           sd.1=0.86, sd.2=0.82, 
                           n.1=19, n.2=13)#Laicite
Effectsizes <- effect_size(eff.type="means", 
                           m.1=2.57, m.2=3.25, 
                           sd.1=0.69, sd.2=0.62, 
                           n.1=18, n.2=12)#Equality
Effectsizes <- effect_size(eff.type="means",
                           m.1=2.59, m.2=2.95, 
                           sd.1=0.81, sd.2=0.76, 
                           n.1=17, n.2=12)#Assimilation
Effectsizes <- effect_size(eff.type="means", 
                           m.1=2.57, m.2=3.12, 
                           sd.1=0.90, sd.2=0.78, 
                           n.1=17, n.2=11)#Mc
# a) Symbolic threat
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.07, m.2=2.45, 
                           sd.1=0.64, sd.2=0.34, 
                           n.1=24, n.2=6)#Control
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.06, m.2=2.66,
                           sd.1=0.80, sd.2=0.86, 
                           n.1=19, n.2=13)#Laicite
Effectsizes <- effect_size(eff.type="means", 
                           m.1=2.72, m.2=3.36, 
                           sd.1=0.69, sd.2=0.69,
                           n.1=18, n.2=12)#Equality
Effectsizes <- effect_size(eff.type="means", 
                           m.1=2.65, m.2=3.11, 
                           sd.1=0.81, sd.2=0.71, 
                           n.1=17, n.2=12)#Assimilation
Effectsizes <- effect_size(eff.type="means", 
                           m.1=2.65, m.2=3.27, 
                           sd.1=1.07, sd.2=1.00, 
                           n.1=17, n.2=11)#Mc
# a) Realistic threat
Effectsizes <- effect_size(eff.type="means",  
                           m.1=2.86, m.2=2.44, 
                           sd.1=0.78, sd.2=0.69, 
                           n.1=24, n.2=6)#Control
Effectsizes <- effect_size(eff.type="means", 
                           m.1=2.70, m.2=2.56, 
                           sd.1=0.89, sd.2=0.96, 
                           n.1=19, n.2=13)#Laicite
Effectsizes <- effect_size(eff.type="means", 
                           m.1=2.41, m.2=3.14, 
                           sd.1=0.87, sd.2=0.81,
                           n.1=18, n.2=12)#Equality
Effectsizes <- effect_size(eff.type="means", 
                           m.1=2.52, m.2=2.78, 
                           sd.1=0.85, sd.2=0.94, 
                           n.1=17, n.2=12)#Assimilation
Effectsizes <- effect_size(eff.type="means", 
                           m.1=2.46, m.2=2.97, 
                           sd.1=0.94, sd.2=0.69, 
                           n.1=17, n.2=11)#Mc
# a) Ingroup bias
Effectsizes <- effect_size(eff.type="means",
                           m.1=1.71, m.2=-0.83, 
                           sd.1=2.14, sd.2=1.33, 
                           n.1=24, n.2=6)#Control
Effectsizes <- effect_size(eff.type="means", 
                           m.1=1.70, m.2=0.07, 
                           sd.1=3.06, sd.2=1.97, 
                           n.1=19, n.2=13)#Laicite
Effectsizes <- effect_size(eff.type="means", 
                           m.1=0.82, m.2=3.17, 
                           sd.1=2.35, sd.2=2.59, 
                           n.1=18, n.2=12)#Equality
Effectsizes <- effect_size(eff.type="means", 
                           m.1=1.32, m.2=1.92, 
                           sd.1=1.79, sd.2=2.61,
                           n.1=17, n.2=12)#Assimilation
Effectsizes <- effect_size(eff.type="means", 
                           m.1=1.72, m.2=3.18, 
                           sd.1=3.05, sd.2=3.12, 
                           n.1=17, n.2=11)#Mc

# Nussio, Bove, & Steele (2019)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=0.59, m.2=0.52, 
                           sd.1=0.31, sd.2=0.33,
                           n.1=3082, n.2=2703)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=0.41, m.2=0.35, 
                           sd.1=0.32, sd.2=0.32,
                           n.1=3001, n.2=2703)

# Obaidi, Kunst, & Kteily (2006)
Effectsizes <- effect_size(eff.type="corr", r=0.36, n=204)
Effectsizes <- effect_size(eff.type="corr", r=0.20, n=204)
Effectsizes <- effect_size(eff.type="corr", r=0.19, n=204)
Effectsizes <- effect_size(eff.type="corr", r=0.44, n=205)
Effectsizes <- effect_size(eff.type="corr", r=0.44, n=205)
Effectsizes <- effect_size(eff.type="corr", r=0.39, n=205)

# Onraet & Van Hiel (2013)
Effectsizes <- effect_size(eff.type="corr", r=0.22, n=588)
Effectsizes <- effect_size(eff.type="corr", r=0.29, n=588)
Effectsizes <- effect_size(eff.type="corr", r=0.35, n=588)
Effectsizes <- effect_size(eff.type="corr", r=0.49, n=588)
Effectsizes <- effect_size(eff.type="corr", r=0.35, n=588)
Effectsizes <- effect_size(eff.type="corr", r=0.30, n=588)
Effectsizes <- effect_size(eff.type="corr", r=0.44, n=588)
Effectsizes <- effect_size(eff.type="corr", r=0.57, n=588)

# Onraet et al. (2013)
# a) Study 1
Effectsizes <- effect_size(eff.type="corr", r=0.14, n=300)
Effectsizes <- effect_size(eff.type="corr", r=0.48, n=300)
Effectsizes <- effect_size(eff.type="corr", r=0.49, n=300)
# b) Study 3
Effectsizes <- effect_size(eff.type="corr", r=0.24, n=800)
Effectsizes <- effect_size(eff.type="corr", r=0.48, n=800)
Effectsizes <- effect_size(eff.type="corr", r=0.52, n=800)
Effectsizes <- effect_size(eff.type="corr", r=0.23, n=800)
Effectsizes <- effect_size(eff.type="corr", r=0.24, n=800)
Effectsizes <- effect_size(eff.type="corr", r=0.22, n=800)

# Oswald (2015)
# a) threat item 1
Effectsizes <- effect_size(eff.type="corr", r=0.13, n=201)
Effectsizes <- effect_size(eff.type="corr", r=-0.01, n=201)
Effectsizes <- effect_size(eff.type="corr", r=0.04, n=201)
Effectsizes <- effect_size(eff.type="corr", r=-0.10, n=201)
# b) threat item 2
Effectsizes <- effect_size(eff.type="corr", r=0.15, n=201)
Effectsizes <- effect_size(eff.type="corr", r=0.09, n=201)
Effectsizes <- effect_size(eff.type="corr", r=0.09, n=201)
Effectsizes <- effect_size(eff.type="corr", r=0.10, n=201)
# c) threat item 3
Effectsizes <- effect_size(eff.type="corr", r=0.07, n=201)
Effectsizes <- effect_size(eff.type="corr", r=0.03, n=201)
Effectsizes <- effect_size(eff.type="corr", r=0.05, n=201)
Effectsizes <- effect_size(eff.type="corr", r=0.03, n=201)
# d) threat item 4
Effectsizes <- effect_size(eff.type="corr", r=0.26, n=201)
Effectsizes <- effect_size(eff.type="corr", r=0.21, n=201)
Effectsizes <- effect_size(eff.type="corr", r=0.26, n=201)
Effectsizes <- effect_size(eff.type="corr", r=0.09, n=201)

# Pearce et al. (2019)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.01, m.2=3.77, 
                           sd.1=0.93, sd.2=0.97, 
                           n.1=484, n.2=484)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.09, m.2=3.95, 
                           sd.1=0.87, sd.2=0.89,
                           n.1=489, n.2=492)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.08, m.2=3.77, 
                           sd.1=0.90, sd.2=0.97, 
                           n.1=490, n.2=484)
Effectsizes <- effect_size(eff.type="means",
                           m.1=4.12, m.2=3.95, 
                           sd.1=0.90, sd.2=0.89, 
                           n.1=488, n.2=492)

# Quirin et al. (2014)
# a) Study 1
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.15, m.2=(5.5+5.32)/2,
                           sd.1=0.71, sd.2=(0.51+0.48)/2, 
                           n.1=33, n.2=32+31)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.74, m.2=(5.7+5.7)/2,
                           sd.1=0.87, sd.2=(0.54+0.78)/2, 
                           n.1=33, n.2=32+31)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=.46, m.2=(.62+.62)/2,
                           sd.1=0.25, sd.2=(0.21+0.22)/2, 
                           n.1=33, n.2=32+31)
# b) Study 2
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.43, m.2=(3.7+3.72)/2,
                           sd.1=0.60, sd.2=(0.36+0.39)/2, 
                           n.1=29, n.2=30+30)

# Rehman & Vanin (2017)
Effectsizes <- effect_size(eff.type="corr", r=(sqrt(0.036)), n=5626)

# Roccas, Klar, & Liviatan (2006)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=5.25, m.2=5.42, 
                           sd.1=1.07, sd.2=1.07,
                           n.1=165, n.2=216)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.39, m.2=3.45, 
                           sd.1=1, sd.2=0.91,
                           n.1=165, n.2=216)

# Routledge et al. (2010)
Effectsizes <- effect_size(eff.type="p.val", p=0.92, n.1=25, n.2=25, tail="two")

# Rovenpor et al. (2019)
# a) Study 1 (B)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=6.56, m.2=6.31, 
                           sd.1=1.91, sd.2=1.83, 
                           n.1=127, n.2=127)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.91, m.2=4.62, 
                           sd.1=1.69, sd.2=1.56, 
                           n.1=127, n.2=127)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=5.42, m.2=4.79, 
                           sd.1=2.15, sd.2=1.99, 
                           n.1=127, n.2=127)
# b) Study 2
Effectsizes <- effect_size(eff.type="means", 
                           m.1=6.07, m.2=5.73, 
                           sd.1=1.65, sd.2=1.83, 
                           n.1=33, n.2=33)#longitudinal
Effectsizes <- effect_size(eff.type="means", 
                           m.1=7.41, m.2=6.95, 
                           sd.1=1.48, sd.2=1.55, 
                           n.1=134, n.2=134)#crosssect.
# c) Study 3
Effectsizes <- effect_size(eff.type="corr", r=0.291, n=203)
Effectsizes <- effect_size(eff.type="corr", r=0.232, n=203)
Effectsizes <- effect_size(eff.type="corr", r=0.168, n=203)

# Sadler et al. (2005)
Effectsizes <- effect_size(eff.type="corr", r=0.43, n=120)
Effectsizes <- effect_size(eff.type="corr", r=-0.29, n=120)
Effectsizes <- effect_size(eff.type="corr", r=-0.29, n=120)
Effectsizes <- effect_size(eff.type="corr", r=0.15, n=120)
Effectsizes <- effect_size(eff.type="corr", r=-0.09, n=120)
Effectsizes <- effect_size(eff.type="corr", r=-0.11, n=120)
Effectsizes <- effect_size(eff.type="corr", r=0.30, n=120)
Effectsizes <- effect_size(eff.type="corr", r=-0.31, n=120)
Effectsizes <- effect_size(eff.type="corr", r=-0.11, n=120)
Effectsizes <- effect_size(eff.type="corr", r=0.18, n=120)
Effectsizes <- effect_size(eff.type="corr", r=-0.17, n=120)
Effectsizes <- effect_size(eff.type="corr", r=-0.08, n=120)
Effectsizes <- effect_size(eff.type="corr", r=0.16, n=120)
Effectsizes <- effect_size(eff.type="corr", r=-0.23, n=120)
Effectsizes <- effect_size(eff.type="corr", r=0.00, n=120)
Effectsizes <- effect_size(eff.type="corr", r=0.12, n=110)
Effectsizes <- effect_size(eff.type="corr", r=-0.13, n=110)
Effectsizes <- effect_size(eff.type="corr", r=-0.04, n=110)
Effectsizes <- effect_size(eff.type="corr", r=0.21, n=120)
Effectsizes <- effect_size(eff.type="corr", r=-0.18, n=120)
Effectsizes <- effect_size(eff.type="corr", r=-0.10, n=120)

# Sahar (2008)
# a) 2001 data
Effectsizes <- effect_size(eff.type="corr", r=0.22, n=113)
Effectsizes <- effect_size(eff.type="corr", r=0.22, n=113)
Effectsizes <- effect_size(eff.type="corr", r=0.18, n=113)
# b) 2005 data
Effectsizes <- effect_size(eff.type="corr", r=0.18, n=166)
Effectsizes <- effect_size(eff.type="corr", r=0.16, n=166)
Effectsizes <- effect_size(eff.type="corr", r=0.07, n=166)
Effectsizes <- effect_size(eff.type="corr", r=0.08, n=166)

# Saleem & Anderson (2013)
# a) Study 1
Effectsizes <- effect_size(eff.type="means", 
                           m.1=0.58, m.2=0.41, 
                           sd.1=0.31, sd.2=0.36, 
                           n.1=59, n.2=59)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=0.50, m.2=0.41, 
                           sd.1=0.27, sd.2=0.36,
                           n.1=59, n.2=59)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=0.15, m.2=-0.22, 
                           sd.1=0.88, sd.2=0.95, 
                           n.1=62, n.2=62)
Effectsizes <- effect_size(eff.type="means",
                           m.1=0.10, m.2=-0.22, 
                           sd.1=0.76, sd.2=0.95, 
                           n.1=62, n.2=62)
Effectsizes <- effect_size(eff.type="means",
                           m.1=1.23, m.2=0.77,
                           sd.1=0.59, sd.2=0.72,
                           n.1=51, n.2=51)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=0.91, m.2=0.77,
                           sd.1=0.69, sd.2=0.72, 
                           n.1=51, n.2=51)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=0.80, m.2=0.16, 
                           sd.1=0.46, sd.2=0.83, 
                           n.1=55, n.2=55)
Effectsizes <- effect_size(eff.type="means",
                           m.1=0.80, m.2=0.52, 
                           sd.1=0.46, sd.2=0.65, 
                           n.1=55, n.2=55)
# B) Study 2
Effectsizes <- effect_size(eff.type="coh.d", d=0.44, n.1=100, n.2=100)
Effectsizes <- effect_size(eff.type="coh.d", d=0.27, n.1=50, n.2=100)
Effectsizes <- effect_size(eff.type="coh.d", d=0.49, n.1=50, n.2=50)
Effectsizes <- effect_size(eff.type="coh.d", d=0.49, n.1=100, n.2=100)
Effectsizes <- effect_size(eff.type="coh.d", d=0.28, n.1=50, n.2=50)
Effectsizes <- effect_size(eff.type="coh.d", d=0.55, n.1=100, n.2=100)
Effectsizes <- effect_size(eff.type="coh.d", d=0.31, n.1=100, n.2=100)
Effectsizes <- effect_size(eff.type="coh.d", d=0.46, n.1=50, n.2=50)

# Saleem et al. (2017)
# a) Study 1
Effectsizes <- effect_size(eff.type="corr", r=0.15, n=715)
Effectsizes <- effect_size(eff.type="corr", r=0.05, n=715)
Effectsizes <- effect_size(eff.type="corr", r=-0.03, n=715)
# b) Study 2
Effectsizes <- effect_size(eff.type="corr", r=0.24, n=200)
Effectsizes <- effect_size(eff.type="corr", r=0.35, n=200)
Effectsizes <- effect_size(eff.type="corr", r=0.32, n=200)
Effectsizes <- effect_size(eff.type="corr", r=0.21, n=200)
# c) Study 3
Effectsizes <- effect_size(eff.type="f.test", f=8.50, n.1=99, n.2=99)
Effectsizes <- effect_size(eff.type="f.test", f=5.25, n.1=99, n.2=99)
Effectsizes <- effect_size(eff.type="f.test", f=4.41, n.1=99, n.2=99)

# Schildkraut (2009)
Effectsizes <- effect_size(eff.type="prop", p1=0.66, p2=0.23, n.ab=1400, n.cd=1400)

# Schuller (2015)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=0.262, m.2=0.281, 
                           sd.1=0.440, sd.2=0.450, 
                           n.1=36982, n.2=688)
Effectsizes <- effect_size(eff.type="means",
                           m.2=0.198, m.1=0.221 , 
                           sd.2=0.399, sd.1=0.415, 
                           n.2=36982, n.1=688)

# Sekerdej & Kossowska (2011)
# Study 1
Effectsizes <- effect_size(eff.type="corr", r=0.33, n=80)
Effectsizes <- effect_size(eff.type="corr", r=0.27, n=80)
# Study 2
Effectsizes <- effect_size(eff.type="corr", r=0.07, n=139)
Effectsizes <- effect_size(eff.type="corr", r=0.34, n=139)

# Shamai & Kimhi
Effectsizes <- effect_size(eff.type="t.test", t=3.89, n.1=353, n.2=65)

# Sharvit (2014)
# a) Study 1
Effectsizes <- effect_size(eff.type="f.test", f=3.43, n.1=41, n.2=41)
Effectsizes <- effect_size(eff.type="f.test", f=4.43, n.1=41, n.2=41)
Effectsizes <- effect_size(eff.type="f.test", f=2.86, n.1=41, n.2=41)

# Sharvit (2014)
Effectsizes <- effect_size(eff.type="f.test", f=20.73, 
                           n.1=(387+389+378+396+392+410), n.2=(393+400+396+403+405+390))
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.22, m.2=3.27, 
                           sd.1=1.16, sd.2=1.14, 
                           n.1=(369+355), n.2=(323+347))

# Shoshani & Slone (2008)
Effectsizes <- effect_size(eff.type="f.test", f=34.60, n.1=148, n.2=148)
Effectsizes <- effect_size(eff.type="f.test", f=37.55, n.1=148, n.2=148)
Effectsizes <- effect_size(eff.type="f.test", f=52.51, n.1=148, n.2=148)
Effectsizes <- effect_size(eff.type="f.test", f=37.42, n.1=148, n.2=148)
Effectsizes <- effect_size(eff.type="f.test", f=28.54, n.1=148, n.2=148)

# Shoshani & Slone (2016)
# a) Stereotype attribution: Jewish sample
es <- esc_mean_gain(pre1mean=91.24, pre1sd=8.47, post1mean=95.82, post1sd=7.15, grp1n=59,
                    pre2mean=93.81, pre2sd=7.41, post2mean=92.85, post2sd=6.88, grp2n=59, 
                    es.type="r")
z.se <- round(1/sqrt((es$totaln)-3), digits=6)
z.var <- z.se^2
pwr <- pwr.r.test(n=es$totaln, r=es$es, sig.level=.05)
results <- data.frame(c(es$es, es$var, es$se,
                        es$zr, z.var, z.se, es$totaln, pwr$power))
results <- transpose(results)
Effectsizes <- rbind(Effectsizes, results)
# b) Stereotype attribution: Arab sample
es <- esc_mean_gain(pre1mean=87.39, pre1sd=6.51, post1mean=90.58, post1sd=7.25, grp1n=55,
                    pre2mean=85.92, pre2sd=6.34, post2mean=86.60, post2sd=6.90, grp2n=55, 
                    es.type="r")
z.se <- round(1/sqrt((es$totaln)-3), digits=6)
z.var <- z.se^2
pwr <- pwr.r.test(n=es$totaln, r=es$es, sig.level=.05)
results <- data.frame(c(es$es, es$var, es$se,
                        es$zr, z.var, z.se, es$totaln, pwr$power))
results <- transpose(results)
Effectsizes <- rbind(Effectsizes, results)
# c) Trust in the adversary: Jewish sample
es <- esc_mean_gain(pre1mean=26.25, pre1sd=5.28, post1mean=21.76, post1sd=4.94, grp1n=59,
                    pre2mean=25.75, pre2sd=5.19, post2mean=26.41, post2sd=4.82, grp2n=59, 
                    es.type="r")
z.se <- round(1/sqrt((es$totaln)-3), digits=6)
z.var <- z.se^2
pwr <- pwr.r.test(n=es$totaln, r=es$es, sig.level=.05)
results <- data.frame(c(-es$es, es$var, es$se,
                        -es$zr, z.var, z.se, es$totaln, pwr$power))
results <- transpose(results)
Effectsizes <- rbind(Effectsizes, results)
# d) Trust in the adversary: Arab sample
es <- esc_mean_gain(pre1mean=24.20, pre1sd=5.30, post1mean=21.30, post1sd=5.58, grp1n=55,
                    pre2mean=25.37, pre2sd=4.62, post2mean=25.63, post2sd=3.69, grp2n=55, 
                    es.type="r")
z.se <- round(1/sqrt((es$totaln)-3), digits=6)
z.var <- z.se^2
pwr <- pwr.r.test(n=es$totaln, r=es$es, sig.level=.05)
results <- data.frame(c(-es$es, es$var, es$se,
                        -es$zr, z.var, z.se, es$totaln, pwr$power))
results <- transpose(results)
Effectsizes <- rbind(Effectsizes, results)
# e) Empathy toward the adversary: Jewish sample
es <- esc_mean_gain(pre1mean=53.62, pre1sd=6.64, post1mean=49.80, post1sd=5.75, grp1n=59,
                    pre2mean=52.15, pre2sd=5.90, post2mean=51.67, post2sd=7.01, grp2n=59, 
                    es.type="r")
z.se <- round(1/sqrt((es$totaln)-3), digits=6)
z.var <- z.se^2
pwr <- pwr.r.test(n=es$totaln, r=es$es, sig.level=.05)
results <- data.frame(c(-es$es, es$var, es$se,
                        -es$zr, z.var, z.se, es$totaln, pwr$power))
results <- transpose(results)
Effectsizes <- rbind(Effectsizes, results)
# f) Empathy toward the adversary: Arab sample
es <- esc_mean_gain(pre1mean=55.25, pre1sd=6.76, post1mean=52.27, post1sd=5.35, grp1n=55,
                    pre2mean=56.80, pre2sd=7.15, post2mean=57.41, post2sd=6.90, grp2n=55, 
                    es.type="r")
z.se <- round(1/sqrt((es$totaln)-3), digits=6)
z.var <- z.se^2
pwr <- pwr.r.test(n=es$totaln, r=es$es, sig.level=.05)
results <- data.frame(c(-es$es, es$var, es$se,
                        -es$zr, z.var, z.se, es$totaln, pwr$power))
results <- transpose(results)
Effectsizes <- rbind(Effectsizes, results)
# g) Willingness to negotiate: Jewish sample
es <- esc_mean_gain(pre1mean=70.41, pre1sd=8.75, post1mean=71.67, post1sd=7.80, grp1n=59,
                    pre2mean=72.15, pre2sd=7.93, post2mean=73.08, post2sd=7.62, grp2n=59, 
                    es.type="r")
z.se <- round(1/sqrt((es$totaln)-3), digits=6)
z.var <- z.se^2
pwr <- pwr.r.test(n=es$totaln, r=es$es, sig.level=.05)
results <- data.frame(c(-es$es, es$var, es$se,
                        -es$zr, z.var, z.se, es$totaln, pwr$power))
results <- transpose(results)
Effectsizes <- rbind(Effectsizes, results)
# h) Willingness to negotiate: Arab sample
es <- esc_mean_gain(pre1mean=74.20, pre1sd=7.90, post1mean=79.80, post1sd=7.15, grp1n=55,
                    pre2mean=73.67, pre2sd=8.10, post2mean=72.90, post2sd=8.40, grp2n=55, 
                    es.type="r")
z.se <- round(1/sqrt((es$totaln)-3), digits=6)
z.var <- z.se^2
pwr <- pwr.r.test(n=es$totaln, r=es$es, sig.level=.05)
results <- data.frame(c(-es$es, es$var, es$se,
                        -es$zr, z.var, z.se, es$totaln, pwr$power))
results <- transpose(results)
Effectsizes <- rbind(Effectsizes, results)
# i) Adversary hostility: Jewish sample
es <- esc_mean_gain(pre1mean=21.80, pre1sd=4.67, post1mean=22.10, post1sd=4.53, grp1n=59, 
                    pre2mean=23.14, pre2sd=4.10, post2mean=25.65, post2sd=4.48, grp2n=59,
                    es.type="r")
z.se <- round(1/sqrt((es$totaln)-3), digits=6)
z.var <- z.se^2
pwr <- pwr.r.test(n=es$totaln, r=es$es, sig.level=.05)
results <- data.frame(c(-es$es, es$var, es$se,
                        -es$zr, z.var, z.se, es$totaln, pwr$power))
results <- transpose(results)
Effectsizes <- rbind(Effectsizes, results)
# j) Adversary hostility: Arab sample
es <- esc_mean_gain(pre1mean=20.60, pre1sd=4.95, post1mean=20.35, post1sd=5.10, grp1n=55, 
                    pre2mean=19.87, pre2sd=5.01, post2mean=23.10, post2sd=5.36, grp2n=55,
                    es.type="r")
z.se <- round(1/sqrt((es$totaln)-3), digits=6)
z.var <- z.se^2
pwr <- pwr.r.test(n=es$totaln, r=es$es, sig.level=.05)
results <- data.frame(c(-es$es, es$var, es$se,
                        -es$zr, z.var, z.se, es$totaln, pwr$power))
results <- transpose(results)
Effectsizes <- rbind(Effectsizes, results)
#reverse all effect sizes to obtain correct direction
Effectsizes[1441:1450,1] <- -Effectsizes[1441:1450,1]
Effectsizes[1441:1450,4] <- -Effectsizes[1441:1450,4]

# Sinclair & LoCicero
Effectsizes <- effect_size(eff.type="corr", r=sqrt(0.021), n=146)
Effectsizes <- effect_size(eff.type="corr", r=sqrt(0.037), n=146)

# Sirin & Geva (2013)
Effectsizes <- effect_size(eff.type="corr", r=(OR2cor(exp(1.436))), n=83)

# Skitka (2005)
Effectsizes <- effect_size(eff.type="corr", r=0.22, n=605)
Effectsizes <- effect_size(eff.type="corr", r=0.15, n=605)
Effectsizes <- effect_size(eff.type="corr", r=0.14, n=605)
Effectsizes <- effect_size(eff.type="corr", r=0.16, n=605)
Effectsizes <- effect_size(eff.type="corr", r=0.20, n=605)
Effectsizes <- effect_size(eff.type="corr", r=0.13, n=605)
Effectsizes <- effect_size(eff.type="corr", r=0.17, n=605)
Effectsizes <- effect_size(eff.type="corr", r=0.22, n=605)
Effectsizes <- effect_size(eff.type="corr", r=-0.12, n=605)

# Skitka et al. (2006)
Effectsizes <- effect_size(eff.type="corr", r=0.08, n=550)
Effectsizes <- effect_size(eff.type="corr", r=0.19, n=550)
Effectsizes <- effect_size(eff.type="corr", r=0.19, n=550)
Effectsizes <- effect_size(eff.type="corr", r=-0.16, n=550)
Effectsizes <- effect_size(eff.type="corr", r=0.34, n=550)
Effectsizes <- effect_size(eff.type="corr", r=0.02, n=550)
Effectsizes <- effect_size(eff.type="corr", r=0.10, n=550)
Effectsizes <- effect_size(eff.type="corr", r=0.19, n=550)

# SSkitka, Bauman, & Mullen (2004)
Effectsizes <- effect_size(eff.type="corr", r=0.36, n=550)
Effectsizes <- effect_size(eff.type="corr", r=0.30, n=550)
Effectsizes <- effect_size(eff.type="corr", r=0.38, n=550)
Effectsizes <- effect_size(eff.type="corr", r=-0.07, n=550)
Effectsizes <- effect_size(eff.type="corr", r=-0.11, n=550)
Effectsizes <- effect_size(eff.type="corr", r=-0.09, n=550)
Effectsizes <- effect_size(eff.type="corr", r=0.21, n=550)
Effectsizes <- effect_size(eff.type="corr", r=0.20, n=550)
Effectsizes <- effect_size(eff.type="corr", r=0.17, n=550)
Effectsizes <- effect_size(eff.type="corr", r=0.12, n=550)
Effectsizes <- effect_size(eff.type="corr", r=0.26, n=550)
Effectsizes <- effect_size(eff.type="corr", r=0.20, n=550)
Effectsizes <- effect_size(eff.type="corr", r=0.14, n=550)
Effectsizes <- effect_size(eff.type="corr", r=0.13, n=550)
Effectsizes <- effect_size(eff.type="corr", r=0.16, n=550)

# Sniderman et al. (2019)
# a) Danish study
Effectsizes <- effect_size(eff.type="p.val", p=0.397, n.1=151, n.2=977, tail="two")
Effectsizes <- effect_size(eff.type="p.val", p=0.512, n.1=137, n.2=913, tail="two")
Effectsizes <- effect_size(eff.type="p.val", p=0.153, n.1=150, n.2=968, tail="two")
Effectsizes <- effect_size(eff.type="p.val", p=0.380, n.1=147, n.2=962, tail="two")
Effectsizes <- effect_size(eff.type="p.val", p=0.193, n.1=128, n.2=854, tail="two")
# b) UK Study
Effectsizes <- effect_size(eff.type="p.val", p=0.9, n.1=1707, n.2=1621, tail="two")
Effectsizes <- effect_size(eff.type="p.val", p=0.9, n.1=1707, n.2=1621, tail="two")

# Solheim (2018)
Effectsizes <- effect_size(eff.type="corr", r=beta2cor(0.047,0), n= 1873)
Effectsizes <- effect_size(eff.type="corr", r=-beta2cor(0.014,0), n= 1873)

# Solomon & Lavi (2018)
Effectsizes <- effect_size(eff.type="corr", r=0.261, n= 307)
Effectsizes <- effect_size(eff.type="corr", r=0.155, n= 269)
Effectsizes <- effect_size(eff.type="corr", r=-0.063, n= 164)
Effectsizes <- effect_size(eff.type="corr", r=-0.167, n= 307)
Effectsizes <- effect_size(eff.type="corr", r=-0.018, n= 269)
Effectsizes <- effect_size(eff.type="corr", r=-0.173, n= 164)

# Steele et al. (2019)
Effectsizes <- effect_size(eff.type="f.test", f=2.60, n.1=101, n.2=83)
Effectsizes <- effect_size(eff.type="f.test", f=0, n.1=101, n.2=83)

# Steele, Parker, & Lickel (2015)
Effectsizes <- effect_size(eff.type="coh.d", d=0.46, n.1=32, n.2=32)
Effectsizes <- effect_size(eff.type="coh.d", d=1.10, n.1=32, n.2=32)

# Steen-Johnson & (2020)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.85, m.2=4.09, 
                           sd.1=1.46, sd.2=1.50,
                           n.1=522, n.2=1019)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.00, m.2=4.34, 
                           sd.1=1.32, sd.2=1.32,
                           n.1=522, n.2=1019)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.26, m.2=4.47, 
                           sd.1=1.21, sd.2=1.24,
                           n.1=522, n.2=1019)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.41, m.2=4.58, 
                           sd.1=1.18, sd.2=1.15,
                           n.1=522, n.2=1019)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.87, m.2=4.69, 
                           sd.1=1, sd.2=1,
                           n.1=1019, n.2=522)

# Stollberg, Fritsche, & Jonas (2017)
Effectsizes <- effect_size(eff.type="corr", r=-0.184, n=74)
Effectsizes <- effect_size(eff.type="corr", r=-0.051, n=74)

# Stroessner et al. (2015)
# a) Study 3
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.1455, m.2=3.2126, 
                           sd.1=1.3715, sd.2=1.4649,
                           n.1=134, n.2=107)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.3333, m.2=3.3112, 
                           sd.1=1.3664, sd.2=1.4018,
                           n.1=75, n.2=49)#Arab
Effectsizes <- effect_size(eff.type="means", 
                           m.1=2.9068, m.2=3.1293, 
                           sd.1=1.3518, sd.2=1.5234,
                           n.1=59, n.2=58)#White
# a) Study 4
Effectsizes <- effect_size(eff.type="means", 
                           m.1=2.6358, m.2=2.8239, 
                           sd.1=1.2619, sd.2=1.1865,
                           n.1=81, n.2=71)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=2.8041, m.2=2.9189, 
                           sd.1=1.2939, sd.2=1.3084,
                           n.1=37, n.2=37)#Arab
Effectsizes <- effect_size(eff.type="means", 
                           m.1=2.4943, m.2=2.7206, 
                           sd.1=1.2312, sd.2=1.0477,
                           n.1=44, n.2=34)#Latino

# Sun, Wu, & Poteyeva (2011) 
Effectsizes <- effect_size(eff.type="corr", r=(OR2cor(1.21)), n=810) 
Effectsizes <- effect_size(eff.type="corr", r=(OR2cor(0.94)), n=810) 
Effectsizes <- effect_size(eff.type="corr", r=(OR2cor(1.07)), n=810) 
Effectsizes <- effect_size(eff.type="corr", r=(OR2cor(0.95)), n=810) 
Effectsizes <- effect_size(eff.type="corr", r=(OR2cor(1.17)), n=810) 
Effectsizes <- effect_size(eff.type="corr", r=(OR2cor(0.78)), n=810) 

# Tabri, Wohl, & Caouette (2018)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.0056, m.2=3.1034, 
                           sd.1=1.5986, sd.2=1.3525, 
                           n.1=60, n.2=58)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.0923, m.2=2.9931, 
                           sd.1=1.4307, sd.2=1.39145, 
                           n.1=48, n.2=65)

# Tamborini et al. (2017)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.10, m.2=3.82, 
                           sd.1=0.61, sd.2=0.75, 
                           n.1=77+84, n.2=77)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=5.66, m.2=3.05, 
                           sd.1=4.03, sd.2=4.40,
                           n.1=85, n.2=35)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=5.29, m.2=5.81, 
                           sd.1=3.51, sd.2=3.24, 
                           n.1=44, n.2=90)

# Tartaglia, Rollero, & Bergagna (2019)
# a) National threat
Effectsizes <- effect_size(eff.type="corr", r=0.19, n= 366)
Effectsizes <- effect_size(eff.type="corr", r=0.43, n= 366)
Effectsizes <- effect_size(eff.type="corr", r=0.21, n= 366)
Effectsizes <- effect_size(eff.type="corr", r=0.47, n= 366)
Effectsizes <- effect_size(eff.type="corr", r=0.17, n= 366)
# b) Personal threat
Effectsizes <- effect_size(eff.type="corr", r=0.23, n= 366)
Effectsizes <- effect_size(eff.type="corr", r=0.38, n= 366)
Effectsizes <- effect_size(eff.type="corr", r=0.25, n= 366)
Effectsizes <- effect_size(eff.type="corr", r=0.35, n= 366)
Effectsizes <- effect_size(eff.type="corr", r=0.04, n= 366)

# Thorisdottir & Jost (2011)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=-0.75, m.2=-1.84, 
                           sd.1=2.37, sd.2=2.01, 
                           n.1=35, n.2=35)

# Tropp et al. (2017)
# a) Study 1: Northern Ireland
Effectsizes <- effect_size(eff.type="corr", r=0.203, n=133+152)
Effectsizes <- effect_size(eff.type="corr", r=0.204, n=133+152)
Effectsizes <- effect_size(eff.type="corr", r=0.219, n=133+152)
Effectsizes <- effect_size(eff.type="corr", r=-0.023, n=133+152)
# a) Study 2: South-Africa
Effectsizes <- effect_size(eff.type="corr", r=0.034, n=205)
Effectsizes <- effect_size(eff.type="corr", r=-0.004, n=205)
Effectsizes <- effect_size(eff.type="corr", r=-0.075, n=205)
Effectsizes <- effect_size(eff.type="corr", r=-0.053, n=205)
Effectsizes <- effect_size(eff.type="corr", r=0.05, n=205)
Effectsizes <- effect_size(eff.type="corr", r=0.003, n=205)
Effectsizes <- effect_size(eff.type="corr", r=0.015, n=205)
Effectsizes <- effect_size(eff.type="corr", r=0.016, n=205)
Effectsizes <- effect_size(eff.type="corr", r=-0.137, n=205)
Effectsizes <- effect_size(eff.type="corr", r=-0.169, n=205)

# Ullrich & Cohrs (2007)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.17, m.2=2.70, 
                           sd.1=0.89, sd.2=0.75, 
                           n.1=38, n.2=40)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.33, m.2=3.06, 
                           sd.1=0.78, sd.2=0.71, 
                           n.1=68, n.2=83)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.17, m.2=3.00, 
                           sd.1=0.86, sd.2=0.59,
                           n.1=13, n.2=(13+14))
Effectsizes <- effect_size(eff.type="means",
                           m.1=3.62, m.2=3.00, 
                           sd.1=0.68, sd.2=0.59,
                           n.1=11, n.2=(13+14))
Effectsizes <- effect_size(eff.type="means", 
                           m.1=2.82, m.2=3.00, 
                           sd.1=0.64, sd.2=0.59, 
                           n.1=11, n.2=(13+14))
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.61, m.2=3.28, 
                           sd.1=0.82, sd.2=0.97, 
                           n.1=24, n.2=26)
# Van Assche (2015)
# VanAssche2015 <- read_sav("https://www.dropbox.com/s/27b0oidnoa09a5m/Van%20Assche%20%282015%29%20Unpublished%20dataset.sav?dl=1") # to replicate results
Effectsizes <- effect_size(eff.type="corr", r=0.342335, n=155)
Effectsizes <- effect_size(eff.type="corr", r=0.388186, n=155)
Effectsizes <- effect_size(eff.type="corr", r=0.178786, n=155)
Effectsizes <- effect_size(eff.type="corr", r=0.213883, n=155)
Effectsizes <- effect_size(eff.type="corr", r=0.251846, n=154)
Effectsizes <- effect_size(eff.type="corr", r=0.330597, n=154)
Effectsizes <- effect_size(eff.type="corr", r=0.308942, n=154)
Effectsizes <- effect_size(eff.type="corr", r=0.369054, n=154)
Effectsizes <- effect_size(eff.type="corr", r=0.249811, n=154)
Effectsizes <- effect_size(eff.type="corr", r=0.255891, n=154)
Effectsizes <- effect_size(eff.type="corr", r=0.172062, n=154)
Effectsizes <- effect_size(eff.type="corr", r=0.234243, n=154)
Effectsizes <- effect_size(eff.type="corr", r=0.253778, n=155)
Effectsizes <- effect_size(eff.type="corr", r=0.340585, n=155)
Effectsizes <- effect_size(eff.type="corr", r=0.122746, n=154)
Effectsizes <- effect_size(eff.type="corr", r=0.238028, n=154)

# Van Assche & Dierckx (2019)
# note: correlation for within-subject design obtained via email correspondence
# a) Study 1a
## time 1 vs. time 2: short-term changes
## positive feelings towards terrorists
es <- esc_mean_sd(grp1m=1.24, grp1sd=0.95, grp1n=86, 
                  grp2m=1.44, grp2sd=1.39, grp2n=86, r=0.161, es.type="r")
z.se <- round(1/sqrt((es$totaln)-3), digits=6)
z.var <- z.se^2
pwr <- pwr.r.test(n=es$totaln, r=es$es, sig.level=.05)
results <- data.frame(c(es$es, es$var, es$se,
                        es$zr, z.var, z.se, es$totaln, pwr$power))
results <- transpose(results)
Effectsizes <- rbind(Effectsizes, results)
## positive feelings towards Muslims
es <- esc_mean_sd(grp1m=5.66, grp1sd=3.32, grp1n=86, 
                  grp2m=5.26, grp2sd=2.44, grp2n=86, r=0.761, es.type="r")
z.se <- round(1/sqrt((es$totaln)-3), digits=6)
z.var <- z.se^2
pwr <- pwr.r.test(n=es$totaln, r=es$es, sig.level=.05)
results <- data.frame(c(es$es, es$var, es$se,
                        es$zr, z.var, z.se, es$totaln, pwr$power))
results <- transpose(results)
Effectsizes <- rbind(Effectsizes, results)
## positive feelings towards refugees
es <- esc_mean_sd(grp1m=5.19, grp1sd=2.37, grp1n=86, 
                  grp2m=5.30, grp2sd=2.29, grp2n=86, r=0.759, es.type="r")
z.se <- round(1/sqrt((es$totaln)-3), digits=6)
z.var <- z.se^2
pwr <- pwr.r.test(n=es$totaln, r=es$es, sig.level=.05)
results <- data.frame(c(es$es, es$var, es$se,
                        es$zr, z.var, z.se, es$totaln, pwr$power))
results <- transpose(results)
Effectsizes <- rbind(Effectsizes, results)
## trust towards immigrants
es <- esc_mean_sd(grp1m=3.73, grp1sd=1.46, grp1n=86, 
                  grp2m=3.72, grp2sd=1.39, grp2n=86, r=0.635, es.type="r")
z.se <- round(1/sqrt((es$totaln)-3), digits=6)
z.var <- z.se^2
pwr <- pwr.r.test(n=es$totaln, r=es$es, sig.level=.05)
results <- data.frame(c(es$es, es$var, es$se,
                        es$zr, z.var, z.se, es$totaln, pwr$power))
results <- transpose(results)
Effectsizes <- rbind(Effectsizes, results)
## threat towards immigrants
es <- esc_mean_sd(grp1m=3.62, grp1sd=1.53, grp1n=86, 
                  grp2m=3.69, grp2sd=1.32, grp2n=86, r=0.686, es.type="r")
z.se <- round(1/sqrt((es$totaln)-3), digits=6)
z.var <- z.se^2
pwr <- pwr.r.test(n=es$totaln, r=es$es, sig.level=.05)
results <- data.frame(c(-es$es, es$var, es$se,
                        -es$zr, z.var, z.se, es$totaln, pwr$power))
results <- transpose(results)
Effectsizes <- rbind(Effectsizes, results)
## prejudice towards immigrants
es <- esc_mean_sd(grp1m=2.11, grp1sd=1.18, grp1n=86, 
                  grp2m=2.03, grp2sd=1.17, grp2n=86, r=0.761, es.type="r")
z.se <- round(1/sqrt((es$totaln)-3), digits=6)
z.var <- z.se^2
pwr <- pwr.r.test(n=es$totaln, r=es$es, sig.level=.05)
results <- data.frame(c(-es$es, es$var, es$se,
                        -es$zr, z.var, z.se, es$totaln, pwr$power))
results <- transpose(results)
Effectsizes <- rbind(Effectsizes, results)

## time 1 vs. time 3: long-term changes
## positive feelings towards terrorists
es <- esc_mean_sd(grp1m=1.04, grp1sd=0.19, grp1n=55, 
                  grp2m=1.44, grp2sd=1.39, grp2n=86, r=-0.055, es.type="r")
z.se <- round(1/sqrt((es$totaln)-3), digits=6)
z.var <- z.se^2
pwr <- pwr.r.test(n=es$totaln, r=es$es, sig.level=.05)
results <- data.frame(c(es$es, es$var, es$se,
                        es$zr, z.var, z.se, es$totaln, pwr$power))
results <- transpose(results)
Effectsizes <- rbind(Effectsizes, results)
## positive feelings towards Muslims
es <- esc_mean_sd(grp1m=6.08, grp1sd=2.13, grp1n=55, 
                  grp2m=5.26, grp2sd=2.44, grp2n=86, r=0.513, es.type="r")
z.se <- round(1/sqrt((es$totaln)-3), digits=6)
z.var <- z.se^2
pwr <- pwr.r.test(n=es$totaln, r=es$es, sig.level=.05)
results <- data.frame(c(es$es, es$var, es$se,
                        es$zr, z.var, z.se, es$totaln, pwr$power))
results <- transpose(results)
Effectsizes <- rbind(Effectsizes, results)
## positive feelings towards refugees
es <- esc_mean_sd(grp1m=5.49, grp1sd=2.34, grp1n=55, 
                  grp2m=5.30, grp2sd=2.29, grp2n=86, r=0.587, es.type="r")
z.se <- round(1/sqrt((es$totaln)-3), digits=6)
z.var <- z.se^2
pwr <- pwr.r.test(n=es$totaln, r=es$es, sig.level=.05)
results <- data.frame(c(es$es, es$var, es$se,
                        es$zr, z.var, z.se, es$totaln, pwr$power))
results <- transpose(results)
Effectsizes <- rbind(Effectsizes, results)
## trust towards immigrants
es <- esc_mean_sd(grp1m=3.98, grp1sd=1.49, grp1n=55, 
                  grp2m=3.72, grp2sd=1.39, grp2n=86, r=0.574, es.type="r")
z.se <- round(1/sqrt((es$totaln)-3), digits=6)
z.var <- z.se^2
pwr <- pwr.r.test(n=es$totaln, r=es$es, sig.level=.05)
results <- data.frame(c(es$es, es$var, es$se,
                        es$zr, z.var, z.se, es$totaln, pwr$power))
results <- transpose(results)
Effectsizes <- rbind(Effectsizes, results)
## threat towards immigrants
es <- esc_mean_sd(grp1m=3.32, grp1sd=1.37, grp1n=55, 
                  grp2m=3.69, grp2sd=1.32, grp2n=86, r=0.577, es.type="r")
z.se <- round(1/sqrt((es$totaln)-3), digits=6)
z.var <- z.se^2
pwr <- pwr.r.test(n=es$totaln, r=es$es, sig.level=.05)
results <- data.frame(c(-es$es, es$var, es$se,
                        -es$zr, z.var, z.se, es$totaln, pwr$power))
results <- transpose(results)
Effectsizes <- rbind(Effectsizes, results)
## prejudice towards immigrants
es <- esc_mean_sd(grp1m=2.13, grp1sd=1.20, grp1n=55, 
                  grp2m=2.03, grp2sd=1.17, grp2n=86, r=0.696, es.type="r")
z.se <- round(1/sqrt((es$totaln)-3), digits=6)
z.var <- z.se^2
pwr <- pwr.r.test(n=es$totaln, r=es$es, sig.level=.05)
results <- data.frame(c(-es$es, es$var, es$se,
                        -es$zr, z.var, z.se, es$totaln, pwr$power))
results <- transpose(results)
Effectsizes <- rbind(Effectsizes, results)

# b) Study 1b
## time 1 vs. time 2: short-term changes
## positive feelings towards refugees
es <- esc_mean_sd(grp1m=5.13, grp1sd=1.28, grp1n=38, 
                  grp2m=5.34, grp2sd=1.24, grp2n=38, r=0.570, es.type="r")
z.se <- round(1/sqrt((es$totaln)-3), digits=6)
z.var <- z.se^2
pwr <- pwr.r.test(n=es$totaln, r=es$es, sig.level=.05)
results <- data.frame(c(es$es, es$var, es$se,
                        es$zr, z.var, z.se, es$totaln, pwr$power))
results <- transpose(results)
Effectsizes <- rbind(Effectsizes, results)
## trust towards refugees
es <- esc_mean_sd(grp1m=4.92, grp1sd=1.36, grp1n=38, 
                  grp2m=5.00, grp2sd=1.34, grp2n=38, r=0.683, es.type="r")
z.se <- round(1/sqrt((es$totaln)-3), digits=6)
z.var <- z.se^2
pwr <- pwr.r.test(n=es$totaln, r=es$es, sig.level=.05)
results <- data.frame(c(es$es, es$var, es$se,
                        es$zr, z.var, z.se, es$totaln, pwr$power))
results <- transpose(results)
Effectsizes <- rbind(Effectsizes, results)
## Threat towards refugees
es <- esc_mean_sd(grp1m=3.17, grp1sd=1.42, grp1n=38, 
                  grp2m=3.26, grp2sd=1.33, grp2n=38, r=0.756, es.type="r")
z.se <- round(1/sqrt((es$totaln)-3), digits=6)
z.var <- z.se^2
pwr <- pwr.r.test(n=es$totaln, r=es$es, sig.level=.05)
results <- data.frame(c(-es$es, es$var, es$se,
                        -es$zr, z.var, z.se, es$totaln, pwr$power))
results <- transpose(results)
Effectsizes <- rbind(Effectsizes, results)
## Avoidance of contact with refugees
es <- esc_mean_sd(grp1m=2.26, grp1sd=1.35, grp1n=38, 
                  grp2m=2.26, grp2sd=1.33, grp2n=38, r=0.503, es.type="r")
z.se <- round(1/sqrt((es$totaln)-3), digits=6)
z.var <- z.se^2
pwr <- pwr.r.test(n=es$totaln, r=es$es, sig.level=.05)
results <- data.frame(c(-es$es, es$var, es$se,
                        -es$zr, z.var, z.se, es$totaln, pwr$power))
results <- transpose(results)
Effectsizes <- rbind(Effectsizes, results)

## time 1 vs. time 3: long-term changes
## positive feelings towards refugees
es <- esc_mean_sd(grp1m=4.87, grp1sd=1.33, grp1n=24, 
                  grp2m=5.34, grp2sd=1.24, grp2n=38, r=0.517, es.type="r")
z.se <- round(1/sqrt((es$totaln)-3), digits=6)
z.var <- z.se^2
pwr <- pwr.r.test(n=es$totaln, r=es$es, sig.level=.05)
results <- data.frame(c(es$es, es$var, es$se,
                        es$zr, z.var, z.se, es$totaln, pwr$power))
results <- transpose(results)
Effectsizes <- rbind(Effectsizes, results)
## trust towards refugees
es <- esc_mean_sd(grp1m=4.83, grp1sd=1.27, grp1n=24, 
                  grp2m=5.00, grp2sd=1.34, grp2n=38, r=0.592, es.type="r")
z.se <- round(1/sqrt((es$totaln)-3), digits=6)
z.var <- z.se^2
pwr <- pwr.r.test(n=es$totaln, r=es$es, sig.level=.05)
results <- data.frame(c(es$es, es$var, es$se,
                        es$zr, z.var, z.se, es$totaln, pwr$power))
results <- transpose(results)
Effectsizes <- rbind(Effectsizes, results)
## Threat towards refugees
es <- esc_mean_sd(grp1m=2.90, grp1sd=1.30, grp1n=24, 
                  grp2m=3.26, grp2sd=1.33, grp2n=38, r=0.399, es.type="r")
z.se <- round(1/sqrt((es$totaln)-3), digits=6)
z.var <- z.se^2
pwr <- pwr.r.test(n=es$totaln, r=es$es, sig.level=.05)
results <- data.frame(c(-es$es, es$var, es$se,
                        -es$zr, z.var, z.se, es$totaln, pwr$power))
results <- transpose(results)
Effectsizes <- rbind(Effectsizes, results)
## Avoidance of contact with refugees
es <- esc_mean_sd(grp1m=2.67, grp1sd=1.66, grp1n=38, 
                  grp2m=2.26, grp2sd=1.33, grp2n=38, r=0.755, es.type="r")
z.se <- round(1/sqrt((es$totaln)-3), digits=6)
z.var <- z.se^2
pwr <- pwr.r.test(n=es$totaln, r=es$es, sig.level=.05)
results <- data.frame(c(-es$es, es$var, es$se,
                        -es$zr, z.var, z.se, es$totaln, pwr$power))
results <- transpose(results)
Effectsizes <- rbind(Effectsizes, results)
Effectsizes[1576:1595,1] <- -Effectsizes[1576:1595,1]
Effectsizes[1576:1595,4] <- -Effectsizes[1576:1595,4]


# Van de Vyver et al. (2016)
# a) Full sample
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.495455, m.2=3.419441, 
                           sd.1=0.913441, sd.2=0.916667, 
                           n.1=1100, n.2=931)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.572273, m.2=3.405478, 
                           sd.1=0.919164, sd.2=0.976984, 
                           n.1=1100, n.2=931)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.087273, m.2=3.949517, 
                           sd.1=0.858306, sd.2=0.963097, 
                           n.1=1100, n.2=931)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.116364, m.2=4.137487, 
                           sd.1=1.476262, sd.2=1.486282, 
                           n.1=1100, n.2=931)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=5.106337, m.2=5.075455, 
                           sd.1=0.994865, sd.2=0.934982, 
                           n.1=931, n.2=1100)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.056928, m.2=3.954545, 
                           sd.1=0.882313, sd.2=0.842311, 
                           n.1=931, n.2=1100)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.362645, m.2=3.443524, 
                           sd.1=1.067454, sd.2=1.133173, 
                           n.1=1100, n.2=931)
# b) Liberals
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.438799, m.2=3.300000, 
                           sd.1=0.957582, sd.2=0.960257, 
                           n.1=433, n.2=345)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.532333, m.2=3.173913, 
                           sd.1=0.962307, sd.2=1.002997, 
                           n.1=433, n.2=345)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.133949, m.2=3.881159, 
                           sd.1=0.825282, sd.2=0.958632, 
                           n.1=433, n.2=345)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.0000, m.2=3.994203, 
                           sd.1=1.465656, sd.2=1.50386, 
                           n.1=433, n.2=345)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=5.194203, m.2=5.083141, 
                           sd.1=0.982386, sd.2=0.921712, 
                           n.1=345, n.2=433)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.150725, m.2=3.95843, 
                           sd.1=0.835551, sd.2=0.859656, 
                           n.1=345, n.2=433)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=2.3903, m.2=2.37971, 
                           sd.1=0.728137, sd.2=0.733853, 
                           n.1=433, n.2=345)
# c) Conservatives
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.595368, m.2=3.646417, 
                           sd.1=0.866281, sd.2=0.895507, 
                           n.1=367, n.2=321)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.667575, m.2=3.632399, 
                           sd.1=0.857216, sd.2=0.947238, 
                           n.1=367, n.2=321)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.1281, m.2=4.077882, 
                           sd.1=0.818136, sd.2=0.930479, 
                           n.1=367, n.2=321)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.286104, m.2=4.23053, 
                           sd.1=1.470065, sd.2=1.480014, 
                           n.1=367, n.2=321)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=5.043614, m.2=5.114441, 
                           sd.1=1.023764, sd.2=0.87652, 
                           n.1=321, n.2=367)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.915888, m.2=3.950954, 
                           sd.1=0.988574, sd.2=0.828319, 
                           n.1=321, n.2=367)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.479564, m.2=4.623053, 
                           sd.1=0.652018, sd.2=0.731841, 
                           n.1=367, n.2=931)

# Van Hauwaert & Huber (2020)
Effectsizes <- effect_size(eff.type="p.val", p=0.0001, 
                           n.1=101, n.2=193, tail="two")
Effectsizes <- effect_size(eff.type="p.val", p=0.677, 
                           n.1=191, n.2=207, tail="two")
Effectsizes <- effect_size(eff.type="p.val", p=0.501, 
                           n.1=191, n.2=207, tail="two")
Effectsizes <- effect_size(eff.type="p.val", p=0.772, 
                           n.1=191, n.2=207, tail="two")

# Vasilopoulos (2018)
# a) Study 1
es <- esc_mean_sd(grp1m=0.60, grp1sd=0.3, grp1n=1513, 
                  grp2m=0.57, grp2sd=0.3, grp2n=1845, r=0.53, es.type="r")
z.se <- round(1/sqrt((es$totaln)-3), digits=6)
z.var <- z.se^2
pwr <- pwr.r.test(n=es$totaln, r=es$es, sig.level=.05)
results <- data.frame(c(es$es, es$var, es$se,
                        es$zr, z.var, z.se, es$totaln, pwr$power))
results <- transpose(results)
Effectsizes <- rbind(Effectsizes, results)
Effectsizes <- effect_size(eff.type="corr", r=-0.13, n=1483)
Effectsizes <- effect_size(eff.type="corr", r=0.16, n=1494)
# b) Study 2
Effectsizes <- effect_size(eff.type="corr", r=-0.02, n=24322)
Effectsizes <- effect_size(eff.type="corr", r=0.01, n=24032)
Effectsizes <- effect_size(eff.type="corr", r=0.05, n=24322)
Effectsizes <- effect_size(eff.type="corr", r=0.04, n=24032)

# Vasilopoulos & Brouard (2019); more correct n from Jost (2019)
# a) Fear
Effectsizes <- effect_size(eff.type="corr", r=-0.045, n=24369)
Effectsizes <- effect_size(eff.type="corr", r=0.1748, n=24325)
Effectsizes <- effect_size(eff.type="corr", r=0.1038, n=22777)
Effectsizes <- effect_size(eff.type="corr", r=0.1230, n=24369)
# b) Anger
Effectsizes <- effect_size(eff.type="corr", r=-0.130, n=24369)
Effectsizes <- effect_size(eff.type="corr", r=0.2377, n=24325)
Effectsizes <- effect_size(eff.type="corr", r=0.2552, n=22777)
Effectsizes <- effect_size(eff.type="corr", r=0.1440, n=24369)

# Vasilopoulos et al. (2019)
Effectsizes <- effect_size(eff.type="corr", r=0.07, n=21311)
Effectsizes <- effect_size(eff.type="corr", r=0.23, n=21311)

# Vasilopoulos, Marcus, & Foucault (2018)
# Attack - Authoritarianism
es <- esc_mean_sd(grp1m=0.46, grp1sd=0.23, grp1n=1431, 
                  grp2m=0.45, grp2sd=0.23, grp2n=1431, r=0.81, es.type="r")
z.se <- round(1/sqrt((es$totaln)-3), digits=6)
z.var <- z.se^2
pwr <- pwr.r.test(n=es$totaln, r=es$es, sig.level=.05)
results <- data.frame(c(es$es, es$var, es$se,
                        es$zr, z.var, z.se, es$totaln, pwr$power))
results <- transpose(results)
Effectsizes <- rbind(Effectsizes, results)
# Fear - Authoritarianism
Effectsizes <- effect_size(eff.type="means", 
                           m.1=0.53, m.2=0.39, 
                           sd.1=0.21, sd.2=0.24,
                           n.1=743, n.2=683)
# Anger - Authoritarianism
Effectsizes <- effect_size(eff.type="means", 
                           m.1=0.48, m.2=0.39, 
                           sd.1=0.23, sd.2=0.23,
                           n.1=1207, n.2=229)
# Attack - left/right self-placement
Effectsizes <- effect_size(eff.type="p.val", p=0.9, n.1=1431, n.2=1431, tail="two")

# Vazquez & Hervas (2010)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.09, m.2=3.08, 
                           sd.1=1.04, sd.2=1.04, 
                           n.1=250, n.2=247)

# Vergani & Tacchi (2016)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.57, m.2=4.03, 
                           sd.1=3.081, sd.2=3.058,
                           n.1=70, n.2=66)

# Vergani et al. (2019)
# a) Study 1
Effectsizes <- effect_size(eff.type="corr", r=0.02, n=26540)
Effectsizes <- effect_size(eff.type="corr", r=0.03, n=26214)
# b) Study 2
# UK
Effectsizes <- effect_size(eff.type="p.val", p=0.664, n.1=67, n.2=63, tail="two")
Effectsizes <- effect_size(eff.type="p.val", p=0.341, n.1=67, n.2=63, tail="two")
Effectsizes <- effect_size(eff.type="p.val", p=0.342, n.1=67, n.2=63, tail="two")
Effectsizes <- effect_size(eff.type="p.val", p=0.125, n.1=67, n.2=63, tail="two")
# France
Effectsizes <- effect_size(eff.type="p.val", p=0.99, n.1=56, n.2=53, tail="two")
Effectsizes <- effect_size(eff.type="p.val", p=0.99, n.1=56, n.2=53, tail="two")
Effectsizes <- effect_size(eff.type="p.val", p=0.99, n.1=54, n.2=53, tail="two")
Effectsizes <- effect_size(eff.type="p.val", p=0.99, n.1=54, n.2=53, tail="two")
# Italy
Effectsizes <- effect_size(eff.type="p.val", p=0.529, n.1=167, n.2=67, tail="two")
Effectsizes <- effect_size(eff.type="p.val", p=0.99, n.1=78, n.2=67, tail="two")
Effectsizes <- effect_size(eff.type="p.val", p=0.441, n.1=167, n.2=67, tail="two")
Effectsizes <- effect_size(eff.type="p.val", p=0.479, n.1=78, n.2=67, tail="two")
# Romania
Effectsizes <- effect_size(eff.type="p.val", p=0.99, n.1=52, n.2=44, tail="two")
Effectsizes <- effect_size(eff.type="p.val", p=0.088, n.1=48, n.2=44, tail="two")
Effectsizes <- effect_size(eff.type="p.val", p=0.385, n.1=52, n.2=44, tail="two")
Effectsizes <- effect_size(eff.type="p.val", p=0.256, n.1=48, n.2=44, tail="two")

# Victoroff et al. (2010)
Effectsizes <- effect_size(eff.type="corr", r=0.31, n=52)
Effectsizes <- effect_size(eff.type="corr", r=0.22, n=52)

# von Sikorski et al. (2017)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=1.99, m.2=1.61, 
                           sd.1=1.26, sd.2=0.89, 
                           n.1=35, n.2=34)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=1.97, m.2=1.61, 
                           sd.1=1.00, sd.2=0.89, 
                           n.1=32, n.2=34)

# Vorsina et al. (2019)
Effectsizes <- effect_size(eff.type="corr", r=0.19, n=1200)
Effectsizes <- effect_size(eff.type="corr", r=0.26, n=1200)
Effectsizes <- effect_size(eff.type="corr", r=0.40, n=1200)
Effectsizes <- effect_size(eff.type="corr", r=0.50, n=1200)
Effectsizes <- effect_size(eff.type="corr", r=0.27, n=1200)
Effectsizes <- effect_size(eff.type="corr", r=0.41, n=1200)

# Washburn & Skitka (2017)
Effectsizes <- effect_size(eff.type="means",
                           m.1=-0.66, m.2=-1.03, 
                           sd.1=2.63, sd.2=2.70, 
                           n.1=62, n.2=63)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=0.01, m.2=-1.03, 
                           sd.1=2.97, sd.2=2.70,
                           n.1=62, n.2=63)

# Wayne (2020)
# note: dataset obtained via email correspondence
# a) study 1
# Effect of terror manipulation on 4 outcome variables.
describeBy(Wayne2020_study1$concessions, group=Wayne2020_study1$terror)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=0.46, m.2=0.45, 
                           sd.1=0.29, sd.2=0.28,
                           n.1=215, n.2=792)
describeBy(Wayne2020_study1$statusquo, group=Wayne2020_study1$terror)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=0.50, m.2=0.50, 
                           sd.1=0.29, sd.2=0.28,
                           n.1=215, n.2=797)
describeBy(Wayne2020_study1$eyeforeye, group=Wayne2020_study1$terror)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=0.76, m.2=0.75, 
                           sd.1=0.26, sd.2=0.26,
                           n.1=801, n.2=216)
describeBy(Wayne2020_study1$onlystrength, group=Wayne2020_study1$terror)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=0.64, m.2=0.63, 
                           sd.1=0.28, sd.2=0.29,
                           n.1=803, n.2=217)
# Correlations between emotions and outcome variables within the terrorism treatments.
# anger
Effectsizes <- effect_size(eff.type="corr", r=0.15113558, n=781)
Effectsizes <- effect_size(eff.type="corr", r=0.07140478, n=786)
Effectsizes <- effect_size(eff.type="corr", r=0.26765966, n=790)
Effectsizes <- effect_size(eff.type="corr", r=0.36397795, n=792)
# fear
Effectsizes <- effect_size(eff.type="corr", r=-0.05161689, n=780)
Effectsizes <- effect_size(eff.type="corr", r=-0.08141443, n=785)
Effectsizes <- effect_size(eff.type="corr", r=0.18842940, n=789)
Effectsizes <- effect_size(eff.type="corr", r=0.27655833, n=790)

# b) study 2
# Effect of terror manipulation on 3 outcome variables.
describeBy(Wayne2020_study2$rankdiplomacy, group=Wayne2020_study2$terror)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=0.44, m.2=0.46, 
                           sd.1=0.41, sd.2=0.42,
                           n.1=306, n.2=598)
describeBy(Wayne2020_study2$vengeance, group=Wayne2020_study2$terror)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=0.70, m.2=0.64, 
                           sd.1=0.29, sd.2=0.31,
                           n.1=616, n.2=320)
describeBy(Wayne2020_study2$deterrence, group=Wayne2020_study2$terror)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=0.70, m.2=0.65, 
                           sd.1=0.28, sd.2=0.31,
                           n.1=616, n.2=320)

# Willer & Adams (2008)
Effectsizes <- effect_size(eff.type="t.test", t=0.05, n.1=645, n.2=637)
Effectsizes[1689,1] <- -Effectsizes[1689,1]
Effectsizes[1689,4] <- -Effectsizes[1689,4]
Effectsizes <- effect_size(eff.type="prop", 
                           p1=0.315, p2=0.310, n.ab=645, n.cd=637)
Effectsizes <- effect_size(eff.type="t.test", t=1.14, n.1=645, n.2=637)
Effectsizes[1691,1] <- -Effectsizes[1691,1]
Effectsizes[1691,4] <- -Effectsizes[1691,4]

# Williamson (2019)
Effectsizes <- effect_size(eff.type="corr", r=0.20, n=1199)
Effectsizes <- effect_size(eff.type="corr", r=0.18, n=1199)
Effectsizes <- effect_size(eff.type="corr", r=0.41, n=1199)
Effectsizes <- effect_size(eff.type="corr", r=0.46, n=1199)

# Wohl & Branscombe (2008)
Effectsizes <- effect_size(eff.type="coh.d", d=0.62, n.1=(54/2), n.2=(54/2))

# Wohl & Branscombe (2009)
# a) Experiment 1
# American sample
Effectsizes <- effect_size(eff.type="means", 
                           m.1=5.5, m.2=4.55, 
                           sd.1=1.18, sd.2=1.73,
                           n.1=25, n.2=25)
Effectsizes <- effect_size(eff.type="f.test", f=0.28, n.1=25, n.2=25)
# Canadian sample
Effectsizes <- effect_size(eff.type="means", 
                           m.1=3.44, m.2=3.88, 
                           sd.1=1.71, sd.2=1.58,
                           n.1=20, n.2=20)
Effectsizes <- effect_size(eff.type="f.test", f=0.28, n.1=20, n.2=20)
# b) Experiment 2
Effectsizes <- effect_size(eff.type="means", 
                           m.1=6.14, m.2=5.14, 
                           sd.1=1.62, sd.2=1.59,
                           n.1=52, n.2=52)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=5.00, m.2=4.32, 
                           sd.1=1.58, sd.2=1.54,
                           n.1=52, n.2=52)

# Woods & Marciniak (2016)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=39.0, m.2=32.6, 
                           sd.1=16.8, sd.2=16.8,
                           n.1=225, n.2=185)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=39.3, m.2=32.9, 
                           sd.1=16.3, sd.2=17.1,
                           n.1=199, n.2=231)

# Zanbar, Kaniasty & Ben-Tzur (2018)
Effectsizes <- effect_size(eff.type="corr", r=0.20, n=764)
Effectsizes <- effect_size(eff.type="corr", r=0.04, n=764)

# Zeitzoff (2019)
Effectsizes <- effect_size(eff.type="p.val", p=0.440, 
                           n.1=(297/2), n.2=(297/2), tail="two")

# Zipris et al. (2019)
Effectsizes <- effect_size(eff.type="corr", r=0.13, n=121)
Effectsizes <- effect_size(eff.type="corr", r=0.10, n=440)
Effectsizes <- effect_size(eff.type="corr", r=0.13, n=440)

# Munoz, Falc?-Gimeno, & Hern?ndez (2020)
Effectsizes <- b_effect_size(eff.type="b", b=0.88, sdy=2, grp1n=(1867/2), grp2n=(1867/2))

# Nussio (2020)
# 3-Day bandwidth
Effectsizes <- effect_size(eff.type="p.val", p=0.815, 
                           n.1=1148, n.2=245, tail="two")
Effectsizes <- effect_size(eff.type="p.val", p=0.801, 
                           n.1=1156, n.2=245, tail="two")
Effectsizes <- effect_size(eff.type="p.val", p=0.026, 
                           n.1=1145, n.2=243, tail="two")
Effectsizes <- effect_size(eff.type="p.val", p=0.298, 
                           n.1=1150, n.2=245, tail="two")
Effectsizes <- effect_size(eff.type="p.val", p=0.550, 
                           n.1=1070, n.2=226, tail="two")
Effectsizes <- effect_size(eff.type="p.val", p=0.565, 
                           n.1=1145, n.2=240, tail="two")
Effectsizes <- effect_size(eff.type="p.val", p=0.675, 
                           n.1=1146, n.2=239, tail="two")
# 30-Day bandwidth
Effectsizes <- effect_size(eff.type="p.val", p=0.376, 
                           n.1=1148, n.2=769, tail="two")
Effectsizes <- effect_size(eff.type="p.val", p=0.491, 
                           n.1=1156, n.2=772, tail="two")
Effectsizes <- effect_size(eff.type="p.val", p=0.138, 
                           n.1=1145, n.2=765, tail="two")
Effectsizes <- effect_size(eff.type="p.val", p=0.925, 
                           n.1=1150, n.2=770, tail="two")
Effectsizes <- effect_size(eff.type="p.val", p=0.852, 
                           n.1=1070, n.2=681, tail="two")
Effectsizes <- effect_size(eff.type="p.val", p=0.308, 
                           n.1=1145, n.2=766, tail="two")
Effectsizes <- effect_size(eff.type="p.val", p=0.139, 
                           n.1=1146, n.2=756, tail="two")
Effectsizes[1716,1] <- -Effectsizes[1716,1]
Effectsizes[1716,4] <- -Effectsizes[1716,4]
Effectsizes[1719:1720,1] <- -Effectsizes[1719:1720,1]
Effectsizes[1719:1720,4] <- -Effectsizes[1719:1720,4]
Effectsizes[1723,1] <- -Effectsizes[1723,1]
Effectsizes[1723,4] <- -Effectsizes[1723,4]

# Balcells and Torrats-Espinosa: one-day estimates without controls 
#(main model in paper, controls account for too much of the r-squared
# to include information from models with controls.
# insignificance used of incumbent effect.)
Effectsizes <- effect_size(eff.type="corr", r=sqrt(0.009), n=3810)
Effectsizes <- effect_size(eff.type="p.val", p=0.9, 
                           n.1=(3810/2), n.2=(3810/2), tail="two") 

# Hatton & Nielsen (2016)
Effectsizes <- effect_size(eff.type="f.test", f=0.005, n.1=93, n.2=187)
Effectsizes <- effect_size(eff.type="f.test", f=0.074, n.1=152, n.2=128)
Effectsizes <- effect_size(eff.type="f.test", f=0.483, n.1=93, n.2=187)
Effectsizes <- effect_size(eff.type="f.test", f=0.763, n.1=152, n.2=128)

# Nagel & Lutter (2020)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=2.923, m.2=2.882, 
                           sd.1=0.647, sd.2=0.654,
                           n.1=271, n.2=532)
Effectsizes <- effect_size(eff.type="means", 
                           m.1=4.300, m.2=4.600, 
                           sd.1=1.879, sd.2=1.850,
                           n.1=287, n.2=571)

					
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#### C. Saving the effect sizes into a csv #### 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
colnames(Effectsizes) <- c("corr", "v_r", "se_r", 
                          "fisher_z", "v_z","se_z", 
                          "ess", "power")
write.csv(Effectsizes, file="data/00-effect-sizes-2.csv")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~FINISHED~~~~~~~~~~~~~~~~~~~ #