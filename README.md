This repository contains the code to fit BayesDTE on survival data. BayesDTE is a Bayesian method for the analysis of clinical trials with a delayed treatment that takes into account the delay through the structure of the hazard ratio.

The Stan and R code for BayesDTE is found in the folder `Stan-R`. Note that dte.stan and dte_funs.stan are coded to work with cmdstanr; it is recommended that you use these files for fitting BayesDTE. 

There is an additional folder called `Stan-R/RStan` within  containing the file `dte_funs_for_rstan.stan`, which includes functions that work with Rstan. They are slightly different because Rstan handles arrays differently than cmdstandr. The purpose of this file is solely to extract the functions into R to visualize the model fit.

An example model fit is included in the folder `Example`. In the file `dte_example.R`, data are first generated using BayesDTE. BayesDTE is then fit on the generated data, and the model fit is visualized by overlaying the posterior means of the survival curves on top of the Kaplan-Meier curves.

Finally, the folder `Data` contains the coordinates used to extract the individual patient data from the published Kaplan-Meier curves of the CASPIAN trial (Paz-Ares et al., 2019, 2022), as well as the R code to turn those coordinates into individual patient data. The CASPIAN trial was a phase III small cell lung cancer clinical trial that exhibited a delayed treatment effect.