### Figures 2B and 3B----
# Paper "Data and model considerations for estimating time-varying functional connectivity in fMRI"
# (Ahrends et al., 2021)
#
# This script requires that the Matlab script mixing_mainscript.m has been run
# and the output saved in a table (.csv format). To reproduce results from Ahrends et al. (2021),
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

packages = c("ggplot2", "wesanderson")

package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

# load table with results for real data
data <- read.csv("realdata_mixing.csv", header = TRUE, sep=",")

# make sure that parcellation is ordered factor
data$parcellation = factor(data$parcellation, levels = c("groupICA50", "groupICA100", "PROFUMO50", "DK80", "Yeo100"))
# calculate number of free parameters (and inverse) and number of observations
K <- 12
data$DF <- K*(K-1)+(K-1)+(K*data$nregions*(data$nregions+1)/2)
data$inv_DF = data$DF^(-1)
data$observations = data$nsubs*data$nts/data$sr

### Figure 2B----

# top panel
fig2b_top <- ggplot(aes(y=mean_maxFO, x=staticFC_similarity, colour=parcellation), data = data)
fig2b_top + geom_point(size=1) + geom_smooth(method = "lm") + 
  scale_color_manual(name = "Parcellation", values = wes_palettes$Darjeeling1[c(4,3,5,1,2)]) + 
  theme_bw() + ylim(c(0.1,NA)) + xlab("FC similarity") + ylab("Mean maxFO")


# bottom panel
fig2b_bottom <- ggplot(aes(y=mean_maxFO, x = parcellation, fill = parcellation), data = data)
fig2b_bottom + theme_bw() + 
  geom_jitter(width = 0.3, size=1, aes(colour = parcellation), alpha = 0.7) + 
  scale_colour_manual(name = "Parcellation", values = wes_palettes$Darjeeling1[c(4,3,5,1,2)]) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5) + 
  scale_fill_manual(name = "Parcellation", values = wes_palettes$Darjeeling1[c(4,3,5,1,2)]) +
  xlab("Parcellation") + ylab("Mean maxFO")

### Figure 3B----

# left panel
fig3B_left <- ggplot(aes(y=mean_maxFO, x=inv_DF, colour=parcellation), data = data)
fig3B_left + geom_jitter(width = 0.00005, size=1) + 
  scale_color_manual(name = "Parcellation", values = wes_palettes$Darjeeling1[c(4,3,5,1,2)]) + 
  theme_bw() + ylim(c(0.1,NA)) + geom_smooth(method = "lm") + 
  xlab("1/Free parameters") + ylab("Mean maxFO")

# middle panel
fig3B_middle <- ggplot(aes(y=mean_maxFO, x=observations, colour=parcellation), data = data)
fig3B_middle + geom_point(size=1) + 
  scale_color_manual(name = "Parcellation", values = wes_palettes$Darjeeling1[c(4,3,5,1,2)]) + 
  theme_bw() + ylim(c(0.1,NA)) + geom_smooth(method = "lm") + 
  xlab("Number of Observations") + ylab("Mean maxFO")

# right panel
fig3B_right <- ggplot(aes(y=mean_maxFO, x=observations*inv_DF, colour=parcellation), data = data)
fig3B_right + geom_point(size=1) + 
  scale_color_manual(name = "Parcellation", values = wes_palettes$Darjeeling1[c(4,3,5,1,2)]) + 
  theme_bw() + ylim(c(0.1,NA)) + geom_smooth(formula = y~log(x)^2) + 
  xlab("Observations/free parameters") + ylab("Mean maxFO")

