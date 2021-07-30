### Structural equation models (SEM)----

# This script requires that the Matlab script mixingtest_mainscript.m has been run
# and the output saved in a table (.csv format). To reproduce results from Ahrends et al (2021),
# the Matlab main script has to be run for all combinations of variables shown in the last part of the script (DO NOT RUN).  
#
# For the real dataset, these are: 
# parcellation = groupICA50, groupICA100, PROFUMO50, Yeo100, DK80
# nsubs = 50, 100, 200
# nts = 200, 500, all
# sr = 1, 2, 3
# nregions = 10, 25, 50, all
#
# For the simulated dataset, these are:
# these_regions = 1:10, 1:50
# n_subj = 100, 20
#
# Christine Ahrends
# (Aarhus University 2020)


# load necessary packages

packages = c("piecewiseSEM", "sem", "nlme", "lme4", "ggplot2")

package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)


### SEM for real data----

# load table with results for real data
data <- read.csv("realdata_mixing.csv", header = TRUE, sep=",")

# make sure that parcellation is ordered factor
data$parcellation = factor(data$parcellation, levels = c("groupICA50", "groupICA100", "PROFUMO50", "DK80", "Yeo100"))
# calculate number of free parameters (and inverse) and number of observations
K <- 12
data$DF <- K*(K-1)+(K-1)+(K*data$nregions*(data$nregions+1)/2)
data$inv_DF = data$DF^(-1)
data$observations = data$nsubs*data$nts/data$sr

# SEM: full (incl. random effects for parcellations) and reduced (without random effects)
SEM_real_full <- psem(
  lme(staticFC_similarity ~ observations, random = ~1|parcellation, data = data, method = "REML"),
  lme(mean_maxFO ~ staticFC_similarity + observations*inv_DF, random = ~staticFC_similarity|parcellation, data = data, method = "REML"),
  data = data, standardize = "scale"
)
summary(SEM_real_full)
rsquared(SEM_real_full)
plot(SEM_real_full) # note that the plot does not include random effects

SEM_real_reduced <- psem(
  lm(staticFC_similarity ~ observations, data = data),
  lm(mean_maxFO ~ staticFC_similarity + observations*inv_DF, data = data),
  data = data, standardize = "scale"
)
summary(SEM_real_reduced)
rsquared(SEM_real_reduced)

# get mean and standard deviation of mean maxFO by parcellation
tapply(data$mean_maxFO, data$parcellation, mean)
tapply(data$mean_maxFO, data$parcellation, sd)


### SEM for simulated data----

# load table containing results for simulated data
data_simu <- read.csv("simudata_mixing.csv", header = TRUE, sep=",")

# calculate number of free parameters (and inverse) and number of observations
K <- 6
data_simu$DF <- K*(K-1)+(K-1)+(K*data_simu$nregions*(data_simu$nregions+1)/2)
data_simu$inv_DF = data_simu$DF^(-1)
data_simu$observations = data_simu$nsubs*data_simu$nts/data_simu$sr

# SEM
simu_SEM <- psem(
  lm(staticFC_similarity ~ subject_inconsistency + observations, data = data_simu),
  lm(mean_maxFO ~ state_inconsistency + staticFC_similarity + inv_DF + observations, data = data_simu),
  data = data_simu, standardize = "scale"
)

# Overview of SEM results and plot
summary(simu_SEM)
plot(simu_SEM)
