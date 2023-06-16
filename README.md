# Calibrated-Imputation-for-Multivariate-Categorical-Data

These files contain R code for mass imputation of categorial data given known totals. The code has been developed for the paper "Calibrated Imputation for Multivariate Categorical Data".

The file "ImputationCode parallel final.R" contains the code for (mass) imputing categorial data given known totals and user-specified edit rules. This code has been parallelized in order to speed up the simulation study.

The file "ImputationCode parallel fast pseudo pop bootstrap final.R" contains the code for estimating the variance of estimates for categories of the variable Occupation by means of a pseudo-population bootstrap approach. This code has been parallelized in order to speed up the simulation study.

The file "Calculate quality measures final.R" contains R code for calculating quality measures. This code has been used for a simulation study.
