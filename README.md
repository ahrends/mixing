# Mixing
## Code for paper "Data and model considerations for estimating time-varying functional connectivity in fMRI" (Ahrends et al. NeuroImage 2022)

Contains all code to run and evaluate HMMs on the HCP data with varying secondary parameters, to simulate fMRI timeseries, run and evaluate HMMs on these simulated timeseries, to quantify the influence of different parameters on model stasis using SEM, and to reproduce figures. Dependencies: [HMM-MAR toolbox](https://github.com/OHBA-analysis/HMM-MAR) and [FSLnets](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSLNets) in Matlab, and [piecewiseSEM](https://cran.r-project.org/web/packages/piecewiseSEM/) in R. 

The main script to run and evaluate HMMs on real and simulated data is mixing_mainscript.m and the main script for SEM is SEM.R. Functions to run and evaluate HMMs on real data are in the folder /real and functions for simulations are in the folder /simulate. Code to reproduce figures 2A and 3A (surface plots simulations) is in the script Fig_2A_3A.m and code to reproduce figures 2B and 3B (scatter and violin plots real data) is in the script Fig_2B_3B.R.
