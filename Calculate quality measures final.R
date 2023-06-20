#############################################################################################################################
# CODE FOR CALCULATING QUALITY MEASURES AFTER MASS IMPUTATION OF CATEGORICAL DATA GIVEN KNOWN TOTALS 
# PROJECT: Multivariate Mass Imputation for the Population Census Given Known Totals (2018-NL-ESS.VIP.BUS.ADMIN)  
# CALL: Multipurpose statistics for efficiency gains in production
#############################################################################################################################

################################
# MAIN PATHS TO FOLDERS WITH DATASETS WITH RESULTS
################################
path <- "PATH TO WORK DiRECTORY"
pathdata <- "PATH TO FOLDER WITH DATA"
pathresults <- "PATH TO FOLDER WITH RESULTS"
setwd(path)

################################
# NUMBERS OF DATASETS PER PART OF SIMULATION STUDY
################################
numdatasets_low_1_iter <- 250
numdatasets_low_5_iter <- 250
numdatasets_low_10_iter <- 250
numdatasets_low_20_iter <- 250
numdatasets_low_MAR_edu <- 250
numdatasets_low_MAR_occ <- 250
numdatasets_low_NMAR_edu <- 250
numdatasets_low_NMAR_occ <- 250
numdatasets_low_NMAR_edu_no_tot <- 250
numdatasets_low_NMAR_occ_no_tot <- 250
numdatasets_low_more_totals_known <- 250
numdatasets_high_1_iter <- 250
numdatasets_high_5_iter <- 250
numdatasets_high_10_iter <- 250
numdatasets_high_20_iter <- 250
numdatasets_high_MAR_edu <- 250
numdatasets_high_MAR_occ <- 250
numdatasets_high_NMAR_edu <- 250
numdatasets_high_NMAR_occ <- 250
numdatasets_high_NMAR_edu_no_tot <- 250
numdatasets_high_NMAR_occ_no_tot <- 250
numdatasets_high_more_totals_known <- 0

################################
# SPECIFIC PATHS TO FOLDERS WITH DATASETS WITH RESULTS / SUFFIX 2 REFERS TO VARIABLE "AGE", 8 TO "EDUCATIONAL LEVEL", 9 TO "ECONOMIC STATUS", AND 10 TO "OCCUPATION"
################################
Low_path <- paste0(pathresults, "Low missing percentages/")
High_path <- paste0(pathresults, "High missing percentages/")
text_1_iter <- "1 iteration/"
text_5_iter <- "5 iterations/"
text_10_iter <- "10 iterations/"
text_20_iter <- "20 iterations/"
text_MAR_edu <- "Selective MAR edu 10 iterations/"
text_MAR_occ <- "Selective MAR occ 10 iterations/"
text_NMAR_edu <- "Selective NMAR edu 10 iterations/"
text_NMAR_occ <- "Selective NMAR occ 10 iterations/"
text_NMAR_edu_no_tot <- "Selective NMAR edu 10 iterations no totals/"
text_NMAR_occ_no_tot <- "Selective NMAR occ 10 iterations no totals/"
text_more_totals_known <- "more totals known 5 iterations/"

text_correctly <- "CorrectlyImputed/CorrectlyImputed_"

text_bias <- "bias/bias_"
text_true <- "bias/truetotals"

text_bias_8_10 <- "bias_8_10/bias_8_10_"
text_true_8_10 <- "bias_8_10/truecrosstotals_8_10"
text_bias_8_2 <- "bias_8_2/bias_8_2_"
text_true_8_2 <- "bias_8_2/truecrosstotals_8_2"
text_bias_8_9 <- "bias_8_9/bias_8_9_"
text_true_8_9 <- "bias_8_9/truecrosstotals_8_9"
text_bias_10_2 <- "bias_10_2/bias_10_2_"
text_true_10_2 <- "bias_10_2/truecrosstotals_10_2"
text_bias_10_9 <- "bias_10_9/bias_10_9_"
text_true_10_9 <- "bias_10_9/truecrosstotals_10_9"

correctly_imputed_low_1_iter <- array(0, dim=c(3, 12))
correctly_imputed_low_5_iter <- array(0, dim=c(3, 12))
correctly_imputed_low_10_iter <- array(0, dim=c(3, 12))
correctly_imputed_low_20_iter <- array(0, dim=c(3, 12))
correctly_imputed_low_MAR_edu <- array(0, dim=c(3, 12))
correctly_imputed_low_MAR_occ <- array(0, dim=c(3, 12))
correctly_imputed_low_NMAR_edu <- array(0, dim=c(3, 12))
correctly_imputed_low_NMAR_occ <- array(0, dim=c(3, 12))
correctly_imputed_low_NMAR_edu_no_tot <- array(0, dim=c(3, 12))
correctly_imputed_low_NMAR_occ_no_tot <- array(0, dim=c(3, 12))
correctly_imputed_low_more_totals_known <- array(0, dim=c(3, 12))
correctly_imputed_high_1_iter <- array(0, dim=c(3, 12))
correctly_imputed_high_5_iter <- array(0, dim=c(3, 12))
correctly_imputed_high_10_iter <- array(0, dim=c(3, 12))
correctly_imputed_high_20_iter <- array(0, dim=c(3, 12))
correctly_imputed_high_MAR_edu <- array(0, dim=c(3, 12))
correctly_imputed_high_MAR_occ <- array(0, dim=c(3, 12))
correctly_imputed_high_NMAR_edu <- array(0, dim=c(3, 12))
correctly_imputed_high_NMAR_occ <- array(0, dim=c(3, 12))
correctly_imputed_high_NMAR_edu_no_tot <- array(0, dim=c(3, 12))
correctly_imputed_high_NMAR_occ_no_tot <- array(0, dim=c(3, 12))
correctly_imputed_high_more_totals_known <- array(0, dim=c(3, 12))

################################
# ARRAYS WITH QUALITY MEASURES
################################
total_bias_low_1_iter <- array(0, dim=c(3, 12, 17))
total_var_low_1_iter <- array(0, dim=c(3, 12, 17))
total_bias_low_5_iter <- array(0, dim=c(3, 12, 17))
total_var_low_5_iter <- array(0, dim=c(3, 12, 17))
total_bias_low_10_iter <- array(0, dim=c(3, 12, 17))
total_var_low_10_iter <- array(0, dim=c(3, 12, 17))
total_bias_low_20_iter <- array(0, dim=c(3, 12, 17))
total_var_low_20_iter <- array(0, dim=c(3, 12, 17))
total_bias_low_MAR_edu <- array(0, dim=c(3, 12, 17))
total_var_low_MAR_edu <- array(0, dim=c(3, 12, 17))
total_bias_low_MAR_occ <- array(0, dim=c(3, 12, 17))
total_var_low_MAR_occ <- array(0, dim=c(3, 12, 17))
total_bias_low_NMAR_edu <- array(0, dim=c(3, 12, 17))
total_var_low_NMAR_edu <- array(0, dim=c(3, 12, 17))
total_bias_low_NMAR_occ <- array(0, dim=c(3, 12, 17))
total_var_low_NMAR_occ <- array(0, dim=c(3, 12, 17))
total_bias_low_NMAR_edu_no_tot <- array(0, dim=c(3, 12, 17))
total_var_low_NMAR_edu_no_tot <- array(0, dim=c(3, 12, 17))
total_bias_low_NMAR_occ_no_tot <- array(0, dim=c(3, 12, 17))
total_var_low_NMAR_occ_no_tot <- array(0, dim=c(3, 12, 17))
total_var_low_more_totals_known <- array(0, dim=c(3, 12, 17))
total_bias_low_more_totals_known <- array(0, dim=c(3, 12, 17))
total_bias_high_1_iter <- array(0, dim=c(3, 12, 17))
total_var_high_1_iter <- array(0, dim=c(3, 12, 17))
total_bias_high_5_iter <- array(0, dim=c(3, 12, 17))
total_var_high_5_iter <- array(0, dim=c(3, 12, 17))
total_bias_high_10_iter <- array(0, dim=c(3, 12, 17))
total_var_high_10_iter <- array(0, dim=c(3, 12, 17))
total_bias_high_20_iter <- array(0, dim=c(3, 12, 17))
total_var_high_20_iter <- array(0, dim=c(3, 12, 17))
total_bias_high_MAR_edu <- array(0, dim=c(3, 12, 17))
total_var_high_MAR_edu <- array(0, dim=c(3, 12, 17))
total_bias_high_MAR_occ <- array(0, dim=c(3, 12, 17))
total_var_high_MAR_occ <- array(0, dim=c(3, 12, 17))
total_bias_high_NMAR_edu <- array(0, dim=c(3, 12, 17))
total_var_high_NMAR_edu <- array(0, dim=c(3, 12, 17))
total_bias_high_NMAR_occ <- array(0, dim=c(3, 12, 17))
total_var_high_NMAR_occ <- array(0, dim=c(3, 12, 17))
total_bias_high_NMAR_edu_no_tot <- array(0, dim=c(3, 12, 17))
total_var_high_NMAR_edu_no_tot <- array(0, dim=c(3, 12, 17))
total_bias_high_NMAR_occ_no_tot <- array(0, dim=c(3, 12, 17))
total_var_high_NMAR_occ_no_tot <- array(0, dim=c(3, 12, 17))
total_bias_high_more_totals_known <- array(0, dim=c(3, 12, 17))
total_var_high_more_totals_known <- array(0, dim=c(3, 12, 17))
truetotals <- array(0, dim=c(12, 17))

total_bias_8_10_low_1_iter <- array(0, dim=c(3, 7, 10))
total_var_8_10_low_1_iter <- array(0, dim=c(3, 7, 10))
total_bias_8_10_low_5_iter <- array(0, dim=c(3, 7, 10))
total_var_8_10_low_5_iter <- array(0, dim=c(3, 7, 10))
total_bias_8_10_low_10_iter <- array(0, dim=c(3, 7, 10))
total_var_8_10_low_10_iter <- array(0, dim=c(3, 7, 10))
total_bias_8_10_low_20_iter <- array(0, dim=c(3, 7, 10))
total_var_8_10_low_20_iter <- array(0, dim=c(3, 7, 10))
total_bias_8_10_low_MAR_edu <- array(0, dim=c(3, 7, 10))
total_var_8_10_low_MAR_edu <- array(0, dim=c(3, 7, 10))
total_bias_8_10_low_MAR_occ <- array(0, dim=c(3, 7, 10))
total_var_8_10_low_MAR_occ <- array(0, dim=c(3, 7, 10))
total_bias_8_10_low_NMAR_edu <- array(0, dim=c(3, 7, 10))
total_var_8_10_low_NMAR_edu <- array(0, dim=c(3, 7, 10))
total_bias_8_10_low_NMAR_occ <- array(0, dim=c(3, 7, 10))
total_var_8_10_low_NMAR_occ <- array(0, dim=c(3, 7, 10))
total_bias_8_10_low_NMAR_edu_no_tot <- array(0, dim=c(3, 7, 10))
total_var_8_10_low_NMAR_edu_no_tot <- array(0, dim=c(3, 7, 10))
total_bias_8_10_low_NMAR_occ_no_tot <- array(0, dim=c(3, 7, 10))
total_var_8_10_low_NMAR_occ_no_tot <- array(0, dim=c(3, 7, 10))
total_bias_8_10_low_more_totals_known <- array(0, dim=c(3, 7, 10))
total_var_8_10_low_more_totals_known <- array(0, dim=c(3, 7, 10))
total_bias_8_10_high_1_iter <- array(0, dim=c(3, 7, 10))
total_var_8_10_high_1_iter <- array(0, dim=c(3, 7, 10))
total_bias_8_10_high_5_iter <- array(0, dim=c(3, 7, 10))
total_var_8_10_high_5_iter <- array(0, dim=c(3, 7, 10))
total_bias_8_10_high_10_iter <- array(0, dim=c(3, 7, 10))
total_var_8_10_high_10_iter <- array(0, dim=c(3, 7, 10))
total_bias_8_10_high_20_iter <- array(0, dim=c(3, 7, 10))
total_var_8_10_high_20_iter <- array(0, dim=c(3, 7, 10))
total_bias_8_10_high_MAR_edu <- array(0, dim=c(3, 7, 10))
total_var_8_10_high_MAR_edu <- array(0, dim=c(3, 7, 10))
total_bias_8_10_high_MAR_occ <- array(0, dim=c(3, 7, 10))
total_var_8_10_high_MAR_occ <- array(0, dim=c(3, 7, 10))
total_bias_8_10_high_NMAR_edu <- array(0, dim=c(3, 7, 10))
total_var_8_10_high_NMAR_edu <- array(0, dim=c(3, 7, 10))
total_bias_8_10_high_NMAR_occ <- array(0, dim=c(3, 7, 10))
total_var_8_10_high_NMAR_occ <- array(0, dim=c(3, 7, 10))
total_bias_8_10_high_NMAR_edu_no_tot <- array(0, dim=c(3, 7, 10))
total_var_8_10_high_NMAR_edu_no_tot <- array(0, dim=c(3, 7, 10))
total_bias_8_10_high_NMAR_occ_no_tot <- array(0, dim=c(3, 7, 10))
total_var_8_10_high_NMAR_occ_no_tot <- array(0, dim=c(3, 7, 10))
total_bias_8_10_high_more_totals_known <- array(0, dim=c(3, 7, 10))
total_var_8_10_high_more_totals_known <- array(0, dim=c(3, 7, 10))
truetotals_8_10 <- array(0, dim=c(7, 10))

total_bias_8_2_low_1_iter <- array(0, dim=c(3, 7, 17))
total_var_8_2_low_1_iter <- array(0, dim=c(3, 7, 17))
total_bias_8_2_low_5_iter <- array(0, dim=c(3, 7, 17))
total_var_8_2_low_5_iter <- array(0, dim=c(3, 7, 17))
total_bias_8_2_low_10_iter <- array(0, dim=c(3, 7, 17))
total_var_8_2_low_10_iter <- array(0, dim=c(3, 7, 17))
total_bias_8_2_low_20_iter <- array(0, dim=c(3, 7, 17))
total_var_8_2_low_20_iter <- array(0, dim=c(3, 7, 17))
total_bias_8_2_low_MAR_edu <- array(0, dim=c(3, 7, 17))
total_var_8_2_low_MAR_edu <- array(0, dim=c(3, 7, 17))
total_bias_8_2_low_MAR_occ <- array(0, dim=c(3, 7, 17))
total_var_8_2_low_MAR_occ <- array(0, dim=c(3, 7, 17))
total_bias_8_2_low_NMAR_edu <- array(0, dim=c(3, 7, 17))
total_var_8_2_low_NMAR_edu <- array(0, dim=c(3, 7, 17))
total_bias_8_2_low_NMAR_occ <- array(0, dim=c(3, 7, 17))
total_var_8_2_low_NMAR_occ <- array(0, dim=c(3, 7, 17))
total_bias_8_2_low_NMAR_edu_no_tot <- array(0, dim=c(3, 7, 17))
total_var_8_2_low_NMAR_edu_no_tot <- array(0, dim=c(3, 7, 17))
total_bias_8_2_low_NMAR_occ_no_tot <- array(0, dim=c(3, 7, 17))
total_var_8_2_low_NMAR_occ_no_tot <- array(0, dim=c(3, 7, 17))
total_bias_8_2_low_more_totals_known <- array(0, dim=c(3, 7, 17))
total_var_8_2_low_more_totals_known <- array(0, dim=c(3, 7, 17))
total_bias_8_2_high_1_iter <- array(0, dim=c(3, 7, 17))
total_var_8_2_high_1_iter <- array(0, dim=c(3, 7, 17))
total_bias_8_2_high_5_iter <- array(0, dim=c(3, 7, 17))
total_var_8_2_high_5_iter <- array(0, dim=c(3, 7, 17))
total_bias_8_2_high_10_iter <- array(0, dim=c(3, 7, 17))
total_var_8_2_high_10_iter <- array(0, dim=c(3, 7, 17))
total_bias_8_2_high_20_iter <- array(0, dim=c(3, 7, 17))
total_var_8_2_high_20_iter <- array(0, dim=c(3, 7, 17))
total_bias_8_2_high_MAR_edu <- array(0, dim=c(3, 7, 17))
total_var_8_2_high_MAR_edu <- array(0, dim=c(3, 7, 17))
total_bias_8_2_high_MAR_occ <- array(0, dim=c(3, 7, 17))
total_var_8_2_high_MAR_occ <- array(0, dim=c(3, 7, 17))
total_bias_8_2_high_NMAR_edu <- array(0, dim=c(3, 7, 17))
total_var_8_2_high_NMAR_edu <- array(0, dim=c(3, 7, 17))
total_bias_8_2_high_NMAR_occ <- array(0, dim=c(3, 7, 17))
total_var_8_2_high_NMAR_occ <- array(0, dim=c(3, 7, 17))
total_bias_8_2_high_NMAR_edu_no_tot <- array(0, dim=c(3, 7, 17))
total_var_8_2_high_NMAR_edu_no_tot <- array(0, dim=c(3, 7, 17))
total_bias_8_2_high_NMAR_occ_no_tot <- array(0, dim=c(3, 7, 17))
total_var_8_2_high_NMAR_occ_no_tot <- array(0, dim=c(3, 7, 17))
total_bias_8_2_high_more_totals_known <- array(0, dim=c(3, 7, 17))
total_var_8_2_high_more_totals_known <- array(0, dim=c(3, 7, 17))
truetotals_8_2 <- array(0, dim=c(7, 17))

total_bias_8_9_low_1_iter <- array(0, dim=c(3, 7, 8))
total_var_8_9_low_1_iter <- array(0, dim=c(3, 7, 8))
total_bias_8_9_low_5_iter <- array(0, dim=c(3, 7, 8))
total_var_8_9_low_5_iter <- array(0, dim=c(3, 7, 8))
total_bias_8_9_low_10_iter <- array(0, dim=c(3, 7, 8))
total_var_8_9_low_10_iter <- array(0, dim=c(3, 7, 8))
total_bias_8_9_low_20_iter <- array(0, dim=c(3, 7, 8))
total_var_8_9_low_20_iter <- array(0, dim=c(3, 7, 8))
total_bias_8_9_low_MAR_edu <- array(0, dim=c(3, 7, 8))
total_var_8_9_low_MAR_edu <- array(0, dim=c(3, 7, 8))
total_bias_8_9_low_MAR_occ <- array(0, dim=c(3, 7, 8))
total_var_8_9_low_MAR_occ <- array(0, dim=c(3, 7, 8))
total_bias_8_9_low_NMAR_edu <- array(0, dim=c(3, 7, 8))
total_var_8_9_low_NMAR_edu <- array(0, dim=c(3, 7, 8))
total_bias_8_9_low_NMAR_occ <- array(0, dim=c(3, 7, 8))
total_var_8_9_low_NMAR_occ <- array(0, dim=c(3, 7, 8))
total_bias_8_9_low_NMAR_edu_no_tot <- array(0, dim=c(3, 7, 8))
total_var_8_9_low_NMAR_edu_no_tot <- array(0, dim=c(3, 7, 8))
total_bias_8_9_low_NMAR_occ_no_tot <- array(0, dim=c(3, 7, 8))
total_var_8_9_low_NMAR_occ_no_tot <- array(0, dim=c(3, 7, 8))
total_bias_8_9_low_more_totals_known <- array(0, dim=c(3, 7, 8))
total_var_8_9_low_more_totals_known <- array(0, dim=c(3, 7, 8))
total_bias_8_9_high_1_iter <- array(0, dim=c(3, 7, 8))
total_var_8_9_high_1_iter <- array(0, dim=c(3, 7, 8))
total_bias_8_9_high_5_iter <- array(0, dim=c(3, 7, 8))
total_var_8_9_high_5_iter <- array(0, dim=c(3, 7, 8))
total_bias_8_9_high_10_iter <- array(0, dim=c(3, 7, 8))
total_var_8_9_high_10_iter <- array(0, dim=c(3, 7, 8))
total_bias_8_9_high_20_iter <- array(0, dim=c(3, 7, 8))
total_var_8_9_high_20_iter <- array(0, dim=c(3, 7, 8))
total_bias_8_9_high_MAR_edu <- array(0, dim=c(3, 7, 8))
total_var_8_9_high_MAR_edu <- array(0, dim=c(3, 7, 8))
total_bias_8_9_high_MAR_occ <- array(0, dim=c(3, 7, 8))
total_var_8_9_high_MAR_occ <- array(0, dim=c(3, 7, 8))
total_bias_8_9_high_NMAR_edu <- array(0, dim=c(3, 7, 8))
total_var_8_9_high_NMAR_edu <- array(0, dim=c(3, 7, 8))
total_bias_8_9_high_NMAR_occ <- array(0, dim=c(3, 7, 8))
total_var_8_9_high_NMAR_occ <- array(0, dim=c(3, 7, 8))
total_bias_8_9_high_NMAR_edu_no_tot <- array(0, dim=c(3, 7, 8))
total_var_8_9_high_NMAR_edu_no_tot <- array(0, dim=c(3, 7, 8))
total_bias_8_9_high_NMAR_occ_no_tot <- array(0, dim=c(3, 7, 8))
total_var_8_9_high_NMAR_occ_no_tot <- array(0, dim=c(3, 7, 8))
total_bias_8_9_high_more_totals_known <- array(0, dim=c(3, 7, 8))
total_var_8_9_high_more_totals_known <- array(0, dim=c(3, 7, 8))
truetotals_8_9 <- array(0, dim=c(7, 8))

total_bias_10_2_low_1_iter <- array(0, dim=c(3, 10, 17))
total_var_10_2_low_1_iter <- array(0, dim=c(3, 10, 17))
total_bias_10_2_low_5_iter <- array(0, dim=c(3, 10, 17))
total_var_10_2_low_5_iter <- array(0, dim=c(3, 10, 17))
total_bias_10_2_low_10_iter <- array(0, dim=c(3, 10, 17))
total_var_10_2_low_10_iter <- array(0, dim=c(3, 10, 17))
total_bias_10_2_low_20_iter <- array(0, dim=c(3, 10, 17))
total_var_10_2_low_20_iter <- array(0, dim=c(3, 10, 17))
total_bias_10_2_low_MAR_edu <- array(0, dim=c(3, 10, 17))
total_var_10_2_low_MAR_edu <- array(0, dim=c(3, 10, 17))
total_bias_10_2_low_MAR_occ <- array(0, dim=c(3, 10, 17))
total_var_10_2_low_MAR_occ <- array(0, dim=c(3, 10, 17))
total_bias_10_2_low_NMAR_edu <- array(0, dim=c(3, 10, 17))
total_var_10_2_low_NMAR_edu <- array(0, dim=c(3, 10, 17))
total_bias_10_2_low_NMAR_occ <- array(0, dim=c(3, 10, 17))
total_var_10_2_low_NMAR_occ <- array(0, dim=c(3, 10, 17))
total_bias_10_2_low_NMAR_edu_no_tot <- array(0, dim=c(3, 10, 17))
total_var_10_2_low_NMAR_edu_no_tot <- array(0, dim=c(3, 10, 17))
total_bias_10_2_low_NMAR_occ_no_tot <- array(0, dim=c(3, 10, 17))
total_var_10_2_low_NMAR_occ_no_tot <- array(0, dim=c(3, 10, 17))
total_bias_10_2_low_more_totals_known <- array(0, dim=c(3, 10, 17))
total_var_10_2_low_more_totals_known <- array(0, dim=c(3, 10, 17))
total_bias_10_2_high_1_iter <- array(0, dim=c(3, 10, 17))
total_var_10_2_high_1_iter <- array(0, dim=c(3, 10, 17))
total_bias_10_2_high_5_iter <- array(0, dim=c(3, 10, 17))
total_var_10_2_high_5_iter <- array(0, dim=c(3, 10, 17))
total_bias_10_2_high_10_iter <- array(0, dim=c(3, 10, 17))
total_var_10_2_high_10_iter <- array(0, dim=c(3, 10, 17))
total_bias_10_2_high_20_iter <- array(0, dim=c(3, 10, 17))
total_var_10_2_high_20_iter <- array(0, dim=c(3, 10, 17))
total_bias_10_2_high_MAR_edu <- array(0, dim=c(3, 10, 17))
total_var_10_2_high_MAR_edu <- array(0, dim=c(3, 10, 17))
total_bias_10_2_high_MAR_occ <- array(0, dim=c(3, 10, 17))
total_var_10_2_high_MAR_occ <- array(0, dim=c(3, 10, 17))
total_bias_10_2_high_NMAR_edu <- array(0, dim=c(3, 10, 17))
total_var_10_2_high_NMAR_edu <- array(0, dim=c(3, 10, 17))
total_bias_10_2_high_NMAR_occ <- array(0, dim=c(3, 10, 17))
total_var_10_2_high_NMAR_occ <- array(0, dim=c(3, 10, 17))
total_bias_10_2_high_NMAR_edu_no_tot <- array(0, dim=c(3, 10, 17))
total_var_10_2_high_NMAR_edu_no_tot <- array(0, dim=c(3, 10, 17))
total_bias_10_2_high_NMAR_occ_no_tot <- array(0, dim=c(3, 10, 17))
total_var_10_2_high_NMAR_occ_no_tot <- array(0, dim=c(3, 10, 17))
total_bias_10_2_high_more_totals_known <- array(0, dim=c(3, 10, 17))
total_var_10_2_high_more_totals_known <- array(0, dim=c(3, 10, 17))
truetotals_10_2 <- array(0, dim=c(10, 17))

total_bias_10_9_low_1_iter <- array(0, dim=c(3, 10, 8))
total_var_10_9_low_1_iter <- array(0, dim=c(3, 10, 8))
total_bias_10_9_low_5_iter <- array(0, dim=c(3, 10, 8))
total_var_10_9_low_5_iter <- array(0, dim=c(3, 10, 8))
total_bias_10_9_low_10_iter <- array(0, dim=c(3, 10, 8))
total_var_10_9_low_10_iter <- array(0, dim=c(3, 10, 8))
total_bias_10_9_low_20_iter <- array(0, dim=c(3, 10, 8))
total_var_10_9_low_20_iter <- array(0, dim=c(3, 10, 8))
total_bias_10_9_low_MAR_edu <- array(0, dim=c(3, 10, 8))
total_var_10_9_low_MAR_edu <- array(0, dim=c(3, 10, 8))
total_bias_10_9_low_MAR_occ <- array(0, dim=c(3, 10, 8))
total_var_10_9_low_MAR_occ <- array(0, dim=c(3, 10, 8))
total_bias_10_9_low_NMAR_edu <- array(0, dim=c(3, 10, 8))
total_var_10_9_low_NMAR_edu <- array(0, dim=c(3, 10, 8))
total_bias_10_9_low_NMAR_occ <- array(0, dim=c(3, 10, 8))
total_var_10_9_low_NMAR_occ <- array(0, dim=c(3, 10, 8))
total_bias_10_9_low_NMAR_edu_no_tot <- array(0, dim=c(3, 10, 8))
total_var_10_9_low_NMAR_edu_no_tot <- array(0, dim=c(3, 10, 8))
total_bias_10_9_low_NMAR_occ_no_tot <- array(0, dim=c(3, 10, 8))
total_var_10_9_low_NMAR_occ_no_tot <- array(0, dim=c(3, 10, 8))
total_bias_10_9_low_more_totals_known <- array(0, dim=c(3, 10, 8))
total_var_10_9_low_more_totals_known <- array(0, dim=c(3, 10, 8))
total_bias_10_9_high_1_iter <- array(0, dim=c(3, 10, 8))
total_var_10_9_high_1_iter <- array(0, dim=c(3, 10, 8))
total_bias_10_9_high_5_iter <- array(0, dim=c(3, 10, 8))
total_var_10_9_high_5_iter <- array(0, dim=c(3, 10, 8))
total_bias_10_9_high_10_iter <- array(0, dim=c(3, 10, 8))
total_var_10_9_high_10_iter <- array(0, dim=c(3, 10, 8))
total_bias_10_9_high_20_iter <- array(0, dim=c(3, 10, 8))
total_var_10_9_high_20_iter <- array(0, dim=c(3, 10, 8))
total_bias_10_9_high_MAR_edu <- array(0, dim=c(3, 10, 8))
total_var_10_9_high_MAR_edu <- array(0, dim=c(3, 10, 8))
total_bias_10_9_high_MAR_occ <- array(0, dim=c(3, 10, 8))
total_var_10_9_high_MAR_occ <- array(0, dim=c(3, 10, 8))
total_bias_10_9_high_NMAR_edu <- array(0, dim=c(3, 10, 8))
total_var_10_9_high_NMAR_edu <- array(0, dim=c(3, 10, 8))
total_bias_10_9_high_NMAR_occ <- array(0, dim=c(3, 10, 8))
total_var_10_9_high_NMAR_occ <- array(0, dim=c(3, 10, 8))
total_bias_10_9_high_NMAR_edu_no_tot <- array(0, dim=c(3, 10, 8))
total_var_10_9_high_NMAR_edu_no_tot <- array(0, dim=c(3, 10, 8))
total_bias_10_9_high_NMAR_occ_no_tot <- array(0, dim=c(3, 10, 8))
total_var_10_9_high_NMAR_occ_no_tot <- array(0, dim=c(3, 10, 8))
total_bias_10_9_high_more_totals_known <- array(0, dim=c(3, 10, 8))
total_var_10_9_high_more_totals_known <- array(0, dim=c(3, 10, 8))
truetotals_10_9 <- array(0, dim=c(10, 8))

################################
# FUNCTIONS FOR CALCULATING AND WRITING QUALITY MEASURES
################################
Calculate_bias_var_singleset <- function(intdatasetpath, intnrows, intncolumns) {
  intbias_singlesettext <- array(NA, dim=c(3, intnrows, intncolumns))
  intbias_singleset <- array(0, dim=c(3, intnrows, intncolumns))
  intbiastext <- read.csv2(intdatasetpath, sep=";", header=FALSE)
  skiplines <- 1
  for (j in (1:intnrows)) {
    for (k in (1:intncolumns)) {
      intbias_singlesettext[1, j, k] <- as.character(intbiastext[j + skiplines, k])
    }
  }
  skiplines <- 2 + intnrows
  for (j in (1:intnrows)) {
    for (k in (1:intncolumns)) {
      intbias_singlesettext[2, j, k] <- as.character(intbiastext[j + skiplines, k])
    }
  }
  skiplines <- 3 + 2*intnrows
  for (j in (1:intnrows)) {
    for (k in (1:intncolumns)) {
      intbias_singlesettext[3, j, k] <- as.character(intbiastext[j + skiplines, k])
    }
  }
  intbias_singleset[,,] <- as.integer(intbias_singlesettext[,,])  
  rm(intbias_singlesettext)
  return(intbias_singleset[,,])
}

Calculate_bias_var_allsets <- function(intdatasetpath, intnumdatasets, intnrows, intncolumns) {
  inttotal_bias <- array(0, dim=c(3, intnrows, intncolumns))
  inttotal_var <- array(0, dim=c(3, intnrows, intncolumns))
  for (s in (1:intnumdatasets)) {
    intbias_singleset <- array(0, dim=c(3, intnrows, intncolumns))
    intdatasetpath2 <- paste0(intdatasetpath, s, ".csv")
    intbias_singleset <- Calculate_bias_var_singleset(intdatasetpath2, intnrows, intncolumns)
    inttotal_bias[,,] <- inttotal_bias[,,] + intbias_singleset[,,]
    inttotal_var[,,] <- inttotal_var[,,] + intbias_singleset[,,]*intbias_singleset[,,]
  }
  return(list("total_bias"=inttotal_bias[,,], "total_var"=inttotal_var))
} 

Find_totals <- function(intdatasetpath, intnrows, intncolumns) {
  inttruetext <- array(NA, dim=c(intnrows, intncolumns))
  inttruetotals <- array(0, dim=c(intnrows, intncolumns))
  inttext <- read.csv2(intdatasetpath, sep=";", header=FALSE)
  for (j in (1:intnrows)) {
    for (k in (1:intncolumns)) {
      inttruetext[j, k] <- as.character(inttext[j, k])
    }
  }
  inttruetotals[,] <- as.integer(inttruetext[,])  
  rm(inttruetext)
  return(inttruetotals[,])
}

Calculate_correctly_imputed_singleset <- function(intdatasetpath, intnvars) {
  intcorrectlytext <- array(NA, dim=c(3, intnvars))
  intcorrectlyvalues <- array(0, dim=c(3, intnvars))
  inttext <- read.csv2(intdatasetpath, sep=";", header=FALSE)
  skiplines <- 1
  for (j in (1:intnvars)) {
    intcorrectlytext[1, j] <- as.character(inttext[j + skiplines, 1])
  }
  skiplines <- 2 + intnvars
  for (j in (1:intnvars)) {
    intcorrectlytext[2, j] <- as.character(inttext[j + skiplines, 1])
  }
  skiplines <- 3 + 2*intnvars
  for (j in (1:intnvars)) {
    intcorrectlytext[3, j] <- as.character(inttext[j + skiplines, 1])
  }
  intcorrectlyvalues[,] <- as.integer(intcorrectlytext[,])  
  rm(intcorrectlytext)
  return(intcorrectlyvalues[,])
}

Calculate_correctly_imputed_all <- function(intdatasetpath, intnumdatasets, intnvars) {
  inttotal_correctly <- array(0, dim=c(3, intnvars))
  for (s in (1:intnumdatasets)) {
    intcorrectly_singleset <- array(0, dim=c(3, intnvars))
    intdatasetpath2 <- paste0(intdatasetpath, s, ".csv")
    intcorrectly_singleset <- Calculate_correctly_imputed_singleset(intdatasetpath2, intnvars)
    inttotal_correctly[,] <- inttotal_correctly[,] + intcorrectly_singleset[,]
  }
  return(inttotal_correctly[,])
}

Write_results <- function(intresultsfile, intresults, intnrows, intncolumns, intNumImputedDataSets) {
  for (approxmethod in (1:3)) {
    write.table(paste0("approx method ", approxmethod), file=intresultsfile, append=TRUE, col.names = FALSE)    
    write.table(intresults[approxmethod,,], file=intresultsfile, append=TRUE, sep=";", col.names = FALSE)
    intsum <- 0
    for (i in (1:intnrows)) {
      for (j in (1:intncolumns)) {
        intsum <- intsum + abs(intresults[approxmethod, i, j])
      }
    }
    inttotresults <- rep(0,3)
    inttotresults[1] <- intsum
    inttotresults[2] <- intsum/intNumImputedDataSets
    inttotresults[3] <- intsum/(intNumImputedDataSets*intnrows*intncolumns)
    write.table(inttotresults, file=intresultsfile, append=TRUE, sep=";", col.names = FALSE)
#    write.table(paste0(intsum, ";", intsum/NumImputedDataSets, ";", intsum/(NumImputedDataSets*intnrows*intncolumns)), file=intresultsfile, append=TRUE, sep=";", col.names = FALSE)
  }  
}

Write_correctly_imp <- function(intresultsfile, intcorrectlyimputed, intNumImputedDataSets) {
  for (approxmethod in (1:3)) {
    write.table(paste0("approx method ", approxmethod), file=intresultsfile, append=TRUE, col.names = FALSE)    
    write.table(intcorrectlyimputed[approxmethod,], file=intresultsfile, append=TRUE, sep=";", col.names = FALSE)
    intsum <- 0
    for (i in (1:12)) {
      intsum <- intsum + abs(intcorrectlyimputed[approxmethod, i])
    }
    write.table("Total", file=intresultsfile, append=TRUE, sep=";", col.names = FALSE)
#    write.table(paste0(intsum, " ", intsum/NumImputedDataSets), file=intresultsfile, append=TRUE, sep=";", col.names = FALSE)
    inttotresults <- rep(0,2)
    inttotresults[1] <- intsum
    inttotresults[2] <- intsum/intNumImputedDataSets
    write.table(inttotresults, file=intresultsfile, append=TRUE, sep=";", col.names = FALSE)
  }  
}

Write_univariate_results <- function(intresultsfile, intresults) {
  for (approxmethod in (1:3)) {
    write.table(paste0("approx method ", approxmethod), file=intresultsfile, append=TRUE, col.names = FALSE)    
    write.table(intresults[approxmethod,,], file=intresultsfile, append=TRUE, sep=";", col.names = FALSE)
    intsum <- 0
    for (i in (1:12)) {
      for (j in (1:17)) {
        intsum <- intsum + abs(intresults[approxmethod,i,j])    
      }
    }
    write.table(intsum, file=intresultsfile, append=TRUE, sep=";", col.names = FALSE)
  }  
}

################################
# ACTUALLY COMPUTING AND WRITING QUALITY MEASURES FOR VARIOUS PARTS OF THE SIMULATION STUDY
################################

##############################
# CORRECTLYIMPUTED: LOW / 1 ITER
##############################
numdatasets <- numdatasets_low_1_iter
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_1_iter, text_correctly)
  correctly_imputed_low_1_iter[,] <- Calculate_correctly_imputed_all(datasetpath, numdatasets, 12) 
  resultsfile <- paste0(datasetpath, "results_correctly_imp.csv")
  write.table("CORRECTLY IMPUTED", file=resultsfile, append=FALSE, col.names = FALSE)
  write.table("LOW; 1 ITER", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_correctly_imp(resultsfile, correctly_imputed_low_1_iter, numdatasets_low_1_iter)
}
##############################
# CORRECTLYIMPUTED: LOW / 5 ITER
##############################
numdatasets <- numdatasets_low_5_iter
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_5_iter, text_correctly)
  correctly_imputed_low_5_iter[,] <- Calculate_correctly_imputed_all(datasetpath, numdatasets, 12)  
  resultsfile <- paste0(datasetpath, "results_correctly_imp.csv")
  write.table("CORRECTLY IMPUTED", file=resultsfile, append=FALSE, col.names = FALSE)
  write.table("LOW; 5 ITER", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_correctly_imp(resultsfile, correctly_imputed_low_5_iter, numdatasets_low_5_iter)
}
##############################
# CORRECTLYIMPUTED: LOW / 10 ITER
##############################
numdatasets <- numdatasets_low_10_iter
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_10_iter, text_correctly) 
  correctly_imputed_low_10_iter[,] <- Calculate_correctly_imputed_all(datasetpath, numdatasets, 12)
  resultsfile <- paste0(datasetpath, "results_correctly_imp.csv")
  write.table("CORRECTLY IMPUTED", file=resultsfile, append=FALSE, col.names = FALSE)
  write.table("LOW; 10 ITER", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_correctly_imp(resultsfile, correctly_imputed_low_10_iter, numdatasets_low_10_iter)
}
##############################
# CORRECTLYIMPUTED: LOW / 20 ITER
##############################
numdatasets <- numdatasets_low_20_iter
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_20_iter, text_correctly) 
  correctly_imputed_low_20_iter[,] <- Calculate_correctly_imputed_all(datasetpath, numdatasets, 12)
  resultsfile <- paste0(datasetpath, "results_correctly_imp.csv")
  write.table("CORRECTLY IMPUTED", file=resultsfile, append=FALSE, col.names = FALSE)
  write.table("LOW; 20 ITER", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_correctly_imp(resultsfile, correctly_imputed_low_20_iter, numdatasets_low_20_iter)
}
##############################
# CORRECTLYIMPUTED: LOW / "Selective MAR edu 10 iterations/"
##############################
numdatasets <- numdatasets_low_MAR_edu
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_MAR_edu, text_correctly) 
  correctly_imputed_low_MAR_edu[,] <- Calculate_correctly_imputed_all(datasetpath, numdatasets, 12)
  resultsfile <- paste0(datasetpath, "results_correctly_imp.csv")
  write.table("CORRECTLY IMPUTED", file=resultsfile, append=FALSE, col.names = FALSE)
  write.table("LOW; 10 ITER", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_correctly_imp(resultsfile, correctly_imputed_low_MAR_edu, numdatasets_low_MAR_edu)
}
##############################
# CORRECTLYIMPUTED: LOW / "Selective MAR occ 10 iterations/"
##############################
numdatasets <- numdatasets_low_MAR_occ
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_MAR_occ, text_correctly) 
  correctly_imputed_low_MAR_occ[,] <- Calculate_correctly_imputed_all(datasetpath, numdatasets, 12)
  resultsfile <- paste0(datasetpath, "results_correctly_imp.csv")
  write.table("CORRECTLY IMPUTED", file=resultsfile, append=FALSE, col.names = FALSE)
  write.table("LOW; 10 ITER", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_correctly_imp(resultsfile, correctly_imputed_low_MAR_occ, numdatasets_low_MAR_occ)
}
##############################
# CORRECTLYIMPUTED: LOW / "Selective NMAR edu 10 iterations/"
##############################
numdatasets <- numdatasets_low_NMAR_edu
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_NMAR_edu, text_correctly) 
  correctly_imputed_low_NMAR_edu[,] <- Calculate_correctly_imputed_all(datasetpath, numdatasets, 12)
  resultsfile <- paste0(datasetpath, "results_correctly_imp.csv")
  write.table("CORRECTLY IMPUTED", file=resultsfile, append=FALSE, col.names = FALSE)
  write.table("LOW; 10 ITER", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_correctly_imp(resultsfile, correctly_imputed_low_NMAR_edu, numdatasets_low_NMAR_edu)
}
##############################
# CORRECTLYIMPUTED: LOW / "Selective NMAR occ 10 iterations/"
##############################
numdatasets <- numdatasets_low_NMAR_occ
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_NMAR_occ, text_correctly) 
  correctly_imputed_low_NMAR_occ[,] <- Calculate_correctly_imputed_all(datasetpath, numdatasets, 12)
  resultsfile <- paste0(datasetpath, "results_correctly_imp.csv")
  write.table("CORRECTLY IMPUTED", file=resultsfile, append=FALSE, col.names = FALSE)
  write.table("LOW; 10 ITER", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_correctly_imp(resultsfile, correctly_imputed_low_NMAR_occ, numdatasets_low_NMAR_occ)
}
##############################
# CORRECTLYIMPUTED: LOW / "Selective NMAR edu 10 iterations no totals/"
##############################
numdatasets <- numdatasets_low_NMAR_edu_no_tot
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_NMAR_edu_no_tot, text_correctly) 
  correctly_imputed_low_NMAR_edu_no_tot[,] <- Calculate_correctly_imputed_all(datasetpath, numdatasets, 12)
  resultsfile <- paste0(datasetpath, "results_correctly_imp.csv")
  write.table("CORRECTLY IMPUTED", file=resultsfile, append=FALSE, col.names = FALSE)
  write.table("LOW; 10 ITER", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_correctly_imp(resultsfile, correctly_imputed_low_NMAR_edu_no_tot, numdatasets_low_NMAR_edu_no_tot)
}
##############################
# CORRECTLYIMPUTED: LOW / "Selective NMAR occ 10 iterations no totals/"
##############################
numdatasets <- numdatasets_low_NMAR_occ_no_tot
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_NMAR_occ_no_tot, text_correctly) 
  correctly_imputed_low_NMAR_occ_no_tot[,] <- Calculate_correctly_imputed_all(datasetpath, numdatasets, 12)
  resultsfile <- paste0(datasetpath, "results_correctly_imp.csv")
  write.table("CORRECTLY IMPUTED", file=resultsfile, append=FALSE, col.names = FALSE)
  write.table("LOW; 10 ITER", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_correctly_imp(resultsfile, correctly_imputed_low_NMAR_occ_no_tot, numdatasets_low_NMAR_occ_no_tot)
}
##############################
# CORRECTLYIMPUTED: LOW / "more totals known 10 iterations/"
##############################
numdatasets <- numdatasets_low_more_totals_known
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_more_totals_known, text_correctly) 
  correctly_imputed_low_more_totals_known[,] <- Calculate_correctly_imputed_all(datasetpath, numdatasets, 12)
  resultsfile <- paste0(datasetpath, "results_correctly_imp.csv")
  write.table("CORRECTLY IMPUTED", file=resultsfile, append=FALSE, col.names = FALSE)
  write.table("LOW; 10 ITER", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_correctly_imp(resultsfile, correctly_imputed_low_more_totals_known, numdatasets_low_more_totals_known)
}
##############################
# CORRECTLYIMPUTED: HIGH / 1 ITER
##############################
numdatasets <- numdatasets_high_1_iter
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_1_iter, text_correctly)
  correctly_imputed_high_1_iter[,] <- Calculate_correctly_imputed_all(datasetpath, numdatasets, 12) 
  resultsfile <- paste0(datasetpath, "results_correctly_imp.csv")
  write.table("CORRECTLY IMPUTED", file=resultsfile, append=FALSE, col.names = FALSE)
  write.table("HIGH; 1 ITER", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_correctly_imp(resultsfile, correctly_imputed_high_1_iter, numdatasets_high_1_iter)
}
##############################
# CORRECTLYIMPUTED: HIGH / 5 ITER
##############################
numdatasets <- numdatasets_high_5_iter
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_5_iter, text_correctly)
  correctly_imputed_high_5_iter[,] <- Calculate_correctly_imputed_all(datasetpath, numdatasets, 12)  
  resultsfile <- paste0(datasetpath, "results_correctly_imp.csv")
  write.table("CORRECTLY IMPUTED", file=resultsfile, append=FALSE, col.names = FALSE)
  write.table("HIGH; 5 ITER", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_correctly_imp(resultsfile, correctly_imputed_high_5_iter, numdatasets_high_5_iter)
}
##############################
# CORRECTLYIMPUTED: HIGH / 10 ITER
##############################
numdatasets <- numdatasets_high_10_iter
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_10_iter, text_correctly) 
  correctly_imputed_high_10_iter[,] <- Calculate_correctly_imputed_all(datasetpath, numdatasets, 12)
  resultsfile <- paste0(datasetpath, "results_correctly_imp.csv")
  write.table("CORRECTLY IMPUTED", file=resultsfile, append=FALSE, col.names = FALSE)
  write.table("HIGH; 10 ITER", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_correctly_imp(resultsfile, correctly_imputed_high_10_iter, numdatasets_high_10_iter)
}
##############################
# CORRECTLYIMPUTED: HIGH / 20 ITER
##############################
numdatasets <- numdatasets_high_20_iter
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_20_iter, text_correctly) 
  correctly_imputed_high_20_iter[,] <- Calculate_correctly_imputed_all(datasetpath, numdatasets, 12)
  resultsfile <- paste0(datasetpath, "results_correctly_imp.csv")
  write.table("CORRECTLY IMPUTED", file=resultsfile, append=FALSE, col.names = FALSE)
  write.table("HIGH; 20 ITER", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_correctly_imp(resultsfile, correctly_imputed_high_20_iter, numdatasets_high_20_iter)
}
##############################
# CORRECTLYIMPUTED: HIGH / "Selective MAR edu 10 iterations/"
##############################
numdatasets <- numdatasets_high_MAR_edu
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_MAR_edu, text_correctly) 
  correctly_imputed_high_MAR_edu[,] <- Calculate_correctly_imputed_all(datasetpath, numdatasets, 12)
  resultsfile <- paste0(datasetpath, "results_correctly_imp.csv")
  write.table("CORRECTLY IMPUTED", file=resultsfile, append=FALSE, col.names = FALSE)
  write.table("HIGH; 10 ITER", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_correctly_imp(resultsfile, correctly_imputed_high_MAR_edu, numdatasets_high_MAR_edu)
}
##############################
# CORRECTLYIMPUTED: HIGH / "Selective MAR occ 10 iterations/"
##############################
numdatasets <- numdatasets_high_MAR_occ
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_MAR_occ, text_correctly) 
  correctly_imputed_high_MAR_occ[,] <- Calculate_correctly_imputed_all(datasetpath, numdatasets, 12)
  resultsfile <- paste0(datasetpath, "results_correctly_imp.csv")
  write.table("CORRECTLY IMPUTED", file=resultsfile, append=FALSE, col.names = FALSE)
  write.table("HIGH; 10 ITER", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_correctly_imp(resultsfile, correctly_imputed_high_MAR_occ, numdatasets_high_MAR_occ)
}
##############################
# CORRECTLYIMPUTED: HIGH / "Selective NMAR edu 10 iterations/"
##############################
numdatasets <- numdatasets_high_NMAR_edu
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_NMAR_edu, text_correctly) 
  correctly_imputed_high_NMAR_edu[,] <- Calculate_correctly_imputed_all(datasetpath, numdatasets, 12)
  resultsfile <- paste0(datasetpath, "results_correctly_imp.csv")
  write.table("CORRECTLY IMPUTED", file=resultsfile, append=FALSE, col.names = FALSE)
  write.table("HIGH; 10 ITER", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_correctly_imp(resultsfile, correctly_imputed_high_NMAR_edu, numdatasets_high_NMAR_edu)
}
##############################
# CORRECTLYIMPUTED: HIGH / "Selective NMAR occ 10 iterations/"
##############################
numdatasets <- numdatasets_high_NMAR_occ
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_NMAR_occ, text_correctly) 
  correctly_imputed_high_NMAR_occ[,] <- Calculate_correctly_imputed_all(datasetpath, numdatasets, 12)
  resultsfile <- paste0(datasetpath, "results_correctly_imp.csv")
  write.table("CORRECTLY IMPUTED", file=resultsfile, append=FALSE, col.names = FALSE)
  write.table("HIGH; 10 ITER", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_correctly_imp(resultsfile, correctly_imputed_high_NMAR_occ, numdatasets_high_NMAR_occ)
}
##############################
# CORRECTLYIMPUTED: HIGH / "Selective NMAR edu 10 iterations no totals/"
##############################
numdatasets <- numdatasets_high_NMAR_edu_no_tot
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_NMAR_edu_no_tot, text_correctly) 
  correctly_imputed_high_NMAR_edu_no_tot[,] <- Calculate_correctly_imputed_all(datasetpath, numdatasets, 12)
  resultsfile <- paste0(datasetpath, "results_correctly_imp.csv")
  write.table("CORRECTLY IMPUTED", file=resultsfile, append=FALSE, col.names = FALSE)
  write.table("HIGH; 10 ITER", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_correctly_imp(resultsfile, correctly_imputed_high_NMAR_edu_no_tot, numdatasets_high_NMAR_edu_no_tot)
}
##############################
# CORRECTLYIMPUTED: HIGH / "Selective NMAR occ 10 iterations no totals/"
##############################
numdatasets <- numdatasets_high_NMAR_occ_no_tot
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_NMAR_occ_no_tot, text_correctly) 
  correctly_imputed_high_NMAR_occ_no_tot[,] <- Calculate_correctly_imputed_all(datasetpath, numdatasets, 12)
  resultsfile <- paste0(datasetpath, "results_correctly_imp.csv")
  write.table("CORRECTLY IMPUTED", file=resultsfile, append=FALSE, col.names = FALSE)
  write.table("HIGH; 10 ITER", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_correctly_imp(resultsfile, correctly_imputed_high_NMAR_occ_no_tot, numdatasets_high_NMAR_occ_no_tot)
}
##############################
# CORRECTLYIMPUTED: HIGH / "more totals known 10 iterations/"
##############################
numdatasets <- numdatasets_high_more_totals_known
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_more_totals_known, text_correctly) 
  correctly_imputed_high_more_totals_known[,] <- Calculate_correctly_imputed_all(datasetpath, numdatasets, 12)
  resultsfile <- paste0(datasetpath, "results_correctly_imp.csv")
  write.table("CORRECTLY IMPUTED", file=resultsfile, append=FALSE, col.names = FALSE)
  write.table("HIGH; 10 ITER", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_correctly_imp(resultsfile, correctly_imputed_high_more_totals_known, numdatasets_high_more_totals_known)
}
##############################

##############################
# UNIVARIATE:
##############################
nrows <- 12
ncolumns <- 17
##############################
# BIAS: TRUETOTALS
##############################
datasetpath <- paste0(Low_path, text_1_iter, text_true, '.csv')
truetotals[,] <- Find_totals(datasetpath, nrows, ncolumns)
##############################
# BIAS: LOW / 1 ITER
##############################
numdatasets <- numdatasets_low_1_iter
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_1_iter, text_bias)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_low_1_iter[,,] <- bias_var_list$total_bias
  total_var_low_1_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath, "quality measures.csv")
  write.table("BIAS", file=resultsfile, append=FALSE, col.names = FALSE)
  write.table("LOW MISSINGNESS; 1 ITER", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_univariate_results(resultsfile, total_bias_low_1_iter)
  write.table("MSE", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_univariate_results(resultsfile, total_var_low_1_iter)
}
##############################
# BIAS: LOW / 5 ITER
##############################
numdatasets <- numdatasets_low_5_iter
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_5_iter, text_bias)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_low_5_iter[,,] <- bias_var_list$total_bias
  total_var_low_5_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath, "quality measures.csv")
  write.table("BIAS", file=resultsfile, append=FALSE, col.names = FALSE)
  write.table("LOW MISSINGNESS; 5 ITER", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_univariate_results(resultsfile, total_bias_low_5_iter)
  write.table("MSE", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_univariate_results(resultsfile, total_var_low_5_iter)
}
##############################
# BIAS: LOW / 10 ITER
##############################
numdatasets <- numdatasets_low_10_iter
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_10_iter, text_bias)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_low_10_iter[,,] <- bias_var_list$total_bias
  total_var_low_10_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath, "quality measures.csv")
  write.table("BIAS", file=resultsfile, append=FALSE, col.names = FALSE)
  write.table("LOW MISSINGNESS; 10 ITER", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_univariate_results(resultsfile, total_bias_low_10_iter)
  write.table("MSE", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_univariate_results(resultsfile, total_var_low_10_iter)
}
##############################
# BIAS: LOW / 20 ITER
##############################
numdatasets <- numdatasets_low_20_iter
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_20_iter, text_bias)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_low_20_iter[,,] <- bias_var_list$total_bias
  total_var_low_20_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath, "quality measures.csv")
  write.table("BIAS", file=resultsfile, append=FALSE, col.names = FALSE)
  write.table("LOW MISSINGNESS; 20 ITER", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_univariate_results(resultsfile, total_bias_low_20_iter)
  write.table("MSE", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_univariate_results(resultsfile, total_var_low_20_iter)
}
##############################
# BIAS: LOW / "Selective MAR edu 10 iterations/"
##############################
numdatasets <- numdatasets_low_MAR_edu
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_MAR_edu, text_bias)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_low_MAR_edu[,,] <- bias_var_list$total_bias
  total_var_low_MAR_edu[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath, "quality measures.csv")
  write.table("BIAS", file=resultsfile, append=FALSE, col.names = FALSE)
  write.table("LOW MISSINGNESS; 10 ITER", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_univariate_results(resultsfile, total_bias_low_MAR_edu)
  write.table("MSE", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_univariate_results(resultsfile, total_var_low_MAR_edu)
}
##############################
# BIAS: LOW / "Selective MAR occ 10 iterations/"
##############################
numdatasets <- numdatasets_low_MAR_occ
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_MAR_occ, text_bias)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_low_MAR_occ[,,] <- bias_var_list$total_bias
  total_var_low_MAR_occ[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath, "quality measures.csv")
  write.table("BIAS", file=resultsfile, append=FALSE, col.names = FALSE)
  write.table("LOW MISSINGNESS; 10 ITER", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_univariate_results(resultsfile, total_bias_low_MAR_occ)
  write.table("MSE", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_univariate_results(resultsfile, total_var_low_MAR_occ)
}
##############################
# BIAS: LOW / "Selective NMAR edu 10 iterations/"
##############################
numdatasets <- numdatasets_low_NMAR_edu
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_NMAR_edu, text_bias)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_low_NMAR_edu[,,] <- bias_var_list$total_bias
  total_var_low_NMAR_edu[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath, "quality measures.csv")
  write.table("BIAS", file=resultsfile, append=FALSE, col.names = FALSE)
  write.table("LOW MISSINGNESS; 10 ITER", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_univariate_results(resultsfile, total_bias_low_NMAR_edu)
  write.table("MSE", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_univariate_results(resultsfile, total_var_low_NMAR_edu)
}
##############################
# BIAS: LOW / "Selective NMAR occ 10 iterations/"
##############################
numdatasets <- numdatasets_low_NMAR_occ
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_NMAR_occ, text_bias)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_low_NMAR_occ[,,] <- bias_var_list$total_bias
  total_var_low_NMAR_occ[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath, "quality measures.csv")
  write.table("BIAS", file=resultsfile, append=FALSE, col.names = FALSE)
  write.table("LOW MISSINGNESS; 10 ITER", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_univariate_results(resultsfile, total_bias_low_NMAR_occ)
  write.table("MSE", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_univariate_results(resultsfile, total_var_low_NMAR_occ)
}
##############################
# BIAS: LOW / "Selective NMAR edu 10 iterations no totals/"
##############################
numdatasets <- numdatasets_low_NMAR_edu_no_tot
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_NMAR_edu_no_tot, text_bias)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_low_NMAR_edu_no_tot[,,] <- bias_var_list$total_bias
  total_var_low_NMAR_edu_no_tot[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath, "quality measures.csv")
  write.table("BIAS", file=resultsfile, append=FALSE, col.names = FALSE)
  write.table("LOW MISSINGNESS; 10 ITER", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_univariate_results(resultsfile, total_bias_low_NMAR_edu_no_tot)
  write.table("MSE", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_univariate_results(resultsfile, total_var_low_NMAR_edu_no_tot)
}
##############################
# BIAS: LOW / "Selective NMAR occ 10 iterations no totals/"
##############################
numdatasets <- numdatasets_low_NMAR_occ_no_tot
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_NMAR_occ_no_tot, text_bias)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_low_NMAR_occ_no_tot[,,] <- bias_var_list$total_bias
  total_var_low_NMAR_occ_no_tot[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath, "quality measures.csv")
  write.table("BIAS", file=resultsfile, append=FALSE, col.names = FALSE)
  write.table("LOW MISSINGNESS; 10 ITER", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_univariate_results(resultsfile, total_bias_low_NMAR_occ_no_tot)
  write.table("MSE", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_univariate_results(resultsfile, total_var_low_NMAR_occ_no_tot)
}
##############################
# BIAS: LOW / "more totals known 10 iterations/"
##############################
numdatasets <- numdatasets_low_more_totals_known
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_more_totals_known, text_bias)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_low_more_totals_known[,,] <- bias_var_list$total_bias
  total_var_low_more_totals_known[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath, "quality measures.csv")
  write.table("BIAS", file=resultsfile, append=FALSE, col.names = FALSE)
  write.table("LOW MISSINGNESS; 10 ITER", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_univariate_results(resultsfile, total_bias_low_more_totals_known)
  write.table("MSE", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_univariate_results(resultsfile, total_var_low_more_totals_known)
}
##############################
# BIAS: HIGH / 1 ITER
##############################
numdatasets <- numdatasets_high_1_iter
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_1_iter, text_bias)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_high_1_iter[,,] <- bias_var_list$total_bias
  total_var_high_1_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath, "quality measures.csv")
  write.table("BIAS", file=resultsfile, append=FALSE, col.names = FALSE)
  write.table("HIGH MISSINGNESS; 1 ITER", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_univariate_results(resultsfile, total_bias_high_1_iter)
  write.table("MSE", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_univariate_results(resultsfile, total_var_high_1_iter)
}
##############################
# BIAS: HIGH / 5 ITER
##############################
numdatasets <- numdatasets_high_5_iter
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_5_iter, text_bias)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_high_5_iter[,,] <- bias_var_list$total_bias
  total_var_high_5_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath, "quality measures.csv")
  write.table("BIAS", file=resultsfile, append=FALSE, col.names = FALSE)
  write.table("HIGH MISSINGNESS; 5 ITER", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_univariate_results(resultsfile, total_bias_high_5_iter)
  write.table("MSE", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_univariate_results(resultsfile, total_var_high_5_iter)
}
##############################
# BIAS: HIGH / 10 ITER
##############################
numdatasets <- numdatasets_high_10_iter
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_10_iter, text_bias)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_high_10_iter[,,] <- bias_var_list$total_bias
  total_var_high_10_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath, "quality measures.csv")
  write.table("BIAS", file=resultsfile, append=FALSE, col.names = FALSE)
  write.table("HIGH MISSINGNESS; 10 ITER", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_univariate_results(resultsfile, total_bias_high_10_iter)
  write.table("MSE", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_univariate_results(resultsfile, total_var_high_10_iter)
}
##############################
# BIAS: HIGH / 20 ITER
##############################
numdatasets <- numdatasets_high_20_iter
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_20_iter, text_bias)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_high_20_iter[,,] <- bias_var_list$total_bias
  total_var_high_20_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath, "quality measures.csv")
  write.table("BIAS", file=resultsfile, append=FALSE, col.names = FALSE)
  write.table("HIGH MISSINGNESS; 20 ITER", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_univariate_results(resultsfile, total_bias_high_20_iter)
  write.table("MSE", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_univariate_results(resultsfile, total_var_high_20_iter)
}
##############################
# BIAS: HIGH / "Selective MAR edu 10 iterations/"
##############################
numdatasets <- numdatasets_high_MAR_edu
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_MAR_edu, text_bias)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_high_MAR_edu[,,] <- bias_var_list$total_bias
  total_var_high_MAR_edu[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath, "quality measures.csv")
  write.table("BIAS", file=resultsfile, append=FALSE, col.names = FALSE)
  write.table("HIGH MISSINGNESS; 10 ITER", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_univariate_results(resultsfile, total_bias_high_MAR_edu)
  write.table("MSE", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_univariate_results(resultsfile, total_var_high_MAR_edu)
}
##############################
# BIAS: HIGH / "Selective MAR occ 10 iterations/"
##############################
numdatasets <- numdatasets_high_MAR_occ
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_MAR_occ, text_bias)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_high_MAR_occ[,,] <- bias_var_list$total_bias
  total_var_high_MAR_occ[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath, "quality measures.csv")
  write.table("BIAS", file=resultsfile, append=FALSE, col.names = FALSE)
  write.table("HIGH MISSINGNESS; 10 ITER", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_univariate_results(resultsfile, total_bias_high_MAR_occ)
  write.table("MSE", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_univariate_results(resultsfile, total_var_high_MAR_occ)
}
##############################
# BIAS: HIGH / "Selective NMAR edu 10 iterations/"
##############################
numdatasets <- numdatasets_high_NMAR_edu
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_NMAR_edu, text_bias)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_high_NMAR_edu[,,] <- bias_var_list$total_bias
  total_var_high_NMAR_edu[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath, "quality measures.csv")
  write.table("BIAS", file=resultsfile, append=FALSE, col.names = FALSE)
  write.table("HIGH MISSINGNESS; 10 ITER", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_univariate_results(resultsfile, total_bias_high_NMAR_edu)
  write.table("MSE", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_univariate_results(resultsfile, total_var_high_NMAR_edu)
}
##############################
# BIAS: HIGH / "Selective NMAR occ 10 iterations/"
##############################
numdatasets <- numdatasets_high_NMAR_occ
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_NMAR_occ, text_bias)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_high_NMAR_occ[,,] <- bias_var_list$total_bias
  total_var_high_NMAR_occ[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath, "quality measures.csv")
  write.table("BIAS", file=resultsfile, append=FALSE, col.names = FALSE)
  write.table("HIGH MISSINGNESS; 10 ITER", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_univariate_results(resultsfile, total_bias_high_NMAR_occ)
  write.table("MSE", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_univariate_results(resultsfile, total_var_high_NMAR_occ)
}
##############################
# BIAS: HIGH / "Selective NMAR edu 10 iterations no totals/"
##############################
numdatasets <- numdatasets_high_NMAR_edu_no_tot
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_NMAR_edu_no_tot, text_bias)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_high_NMAR_edu_no_tot[,,] <- bias_var_list$total_bias
  total_var_high_NMAR_edu_no_tot[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath, "quality measures.csv")
  write.table("BIAS", file=resultsfile, append=FALSE, col.names = FALSE)
  write.table("HIGH MISSINGNESS; 10 ITER", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_univariate_results(resultsfile, total_bias_high_NMAR_edu_no_tot)
  write.table("MSE", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_univariate_results(resultsfile, total_var_high_NMAR_edu_no_tot)
}
##############################
# BIAS: HIGH / "Selective NMAR occ 10 iteration no totals/"
##############################
numdatasets <- numdatasets_high_NMAR_occ_no_tot
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_NMAR_occ_no_tot, text_bias)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_high_NMAR_occ_no_tot[,,] <- bias_var_list$total_bias
  total_var_high_NMAR_occ_no_tot[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath, "quality measures.csv")
  write.table("BIAS", file=resultsfile, append=FALSE, col.names = FALSE)
  write.table("HIGH MISSINGNESS; 10 ITER", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_univariate_results(resultsfile, total_bias_high_NMAR_occ_no_tot)
  write.table("MSE", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_univariate_results(resultsfile, total_var_high_NMAR_occ_no_tot)
}
##############################
# BIAS: HIGH / "more totals known 10 iteration/"
##############################
numdatasets <- numdatasets_high_more_totals_known
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_more_totals_known, text_bias)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_high_more_totals_known[,,] <- bias_var_list$total_bias
  total_var_high_more_totals_known[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath, "quality measures.csv")
  write.table("BIAS", file=resultsfile, append=FALSE, col.names = FALSE)
  write.table("HIGH MISSINGNESS; 10 ITER", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_univariate_results(resultsfile, total_bias_high_more_totals_known)
  write.table("MSE", file=resultsfile, append=TRUE, col.names = FALSE)
  Write_univariate_results(resultsfile, total_var_high_more_totals_known)
}
##############################

##############################
# CROSS-TABLE 8_10:
##############################
nrows <- 7
ncolumns <- 10
##############################
# BIAS_8_10: TRUETOTALS
##############################
datasetpath <- paste0(Low_path, text_1_iter, text_true_8_10, '.csv')
truetotals_8_10[,] <- Find_totals(datasetpath, nrows, ncolumns)
##############################
# BIAS_8_10: LOW / 1 ITER
##############################
numdatasets <- numdatasets_low_1_iter
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_1_iter, text_bias_8_10)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_10_low_1_iter[,,] <- bias_var_list$total_bias
  total_var_8_10_low_1_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_10; LOW MISSINGNESS; 1 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_10_low_1_iter, nrows, ncolumns, numdatasets_low_1_iter)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_10_low_1_iter, nrows, ncolumns, numdatasets_low_1_iter)
}
##############################
# BIAS_8_10: LOW / 5 ITER
##############################
numdatasets <- numdatasets_low_5_iter
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_5_iter, text_bias_8_10)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_10_low_5_iter[,,] <- bias_var_list$total_bias
  total_var_8_10_low_5_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_10; LOW MISSINGNESS; 5 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_10_low_5_iter, nrows, ncolumns, numdatasets_low_5_iter)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_10_low_5_iter, nrows, ncolumns, numdatasets_low_5_iter)
}
##############################
# BIAS_8_10: LOW / 10 ITER
##############################
numdatasets <- numdatasets_low_10_iter
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_10_iter, text_bias_8_10)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_10_low_10_iter[,,] <- bias_var_list$total_bias
  total_var_8_10_low_10_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_10; LOW MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_10_low_10_iter, nrows, ncolumns, numdatasets_low_10_iter)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_10_low_10_iter, nrows, ncolumns, numdatasets_low_10_iter)
}
##############################
# BIAS_8_10: LOW / 20 ITER
##############################
numdatasets <- numdatasets_low_20_iter
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_20_iter, text_bias_8_10)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_10_low_20_iter[,,] <- bias_var_list$total_bias
  total_var_8_10_low_20_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_10; LOW MISSINGNESS; 20 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_10_low_20_iter, nrows, ncolumns, numdatasets_low_20_iter)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_10_low_20_iter, nrows, ncolumns, numdatasets_low_20_iter)
}
##############################
# BIAS_8_10: LOW / "Selective MAR edu 10 iterations/"
##############################
numdatasets <- numdatasets_low_MAR_edu
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_MAR_edu, text_bias_8_10)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_10_low_MAR_edu[,,] <- bias_var_list$total_bias
  total_var_8_10_low_MAR_edu[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_10; LOW MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_10_low_MAR_edu, nrows, ncolumns, numdatasets_low_MAR_edu)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_10_low_MAR_edu, nrows, ncolumns, numdatasets_low_MAR_edu)
}
##############################
# BIAS_8_10: LOW / "Selective MAR occ 10 iterations/"
##############################
numdatasets <- numdatasets_low_MAR_occ
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_MAR_occ, text_bias_8_10)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_10_low_MAR_occ[,,] <- bias_var_list$total_bias
  total_var_8_10_low_MAR_occ[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_10; LOW MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_10_low_MAR_occ, nrows, ncolumns, numdatasets_low_MAR_occ)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_10_low_MAR_occ, nrows, ncolumns, numdatasets_low_MAR_occ)
}
##############################
# BIAS_8_10: LOW / "Selective NMAR edu 10 iterations/"
##############################
numdatasets <- numdatasets_low_NMAR_edu
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_NMAR_edu, text_bias_8_10)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_10_low_NMAR_edu[,,] <- bias_var_list$total_bias
  total_var_8_10_low_NMAR_edu[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_10; LOW MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_10_low_NMAR_edu, nrows, ncolumns, numdatasets_low_NMAR_edu)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_10_low_NMAR_edu, nrows, ncolumns, numdatasets_low_NMAR_edu)
}
##############################
# BIAS_8_10: LOW / "Selective NMAR occ 10 iterations/"
##############################
numdatasets <- numdatasets_low_NMAR_occ
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_NMAR_occ, text_bias_8_10)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_10_low_NMAR_occ[,,] <- bias_var_list$total_bias
  total_var_8_10_low_NMAR_occ[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_10; LOW MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_10_low_NMAR_occ, nrows, ncolumns, numdatasets_low_NMAR_occ)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_10_low_NMAR_occ, nrows, ncolumns, numdatasets_low_NMAR_occ)
}
##############################
# BIAS_8_10: LOW / "Selective NMAR edu 10 iterations no totals/"
##############################
numdatasets <- numdatasets_low_NMAR_edu_no_tot
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_NMAR_edu_no_tot, text_bias_8_10)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_10_low_NMAR_edu_no_tot[,,] <- bias_var_list$total_bias
  total_var_8_10_low_NMAR_edu_no_tot[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_10; LOW MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_10_low_NMAR_edu_no_tot, nrows, ncolumns, numdatasets_low_NMAR_edu_no_tot)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_10_low_NMAR_edu_no_tot, nrows, ncolumns, numdatasets_low_NMAR_edu_no_tot)
}
##############################
# BIAS_8_10: LOW / "Selective NMAR occ 10 iterations no totals/"
##############################
numdatasets <- numdatasets_low_NMAR_occ_no_tot
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_NMAR_occ_no_tot, text_bias_8_10)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_10_low_NMAR_occ_no_tot[,,] <- bias_var_list$total_bias
  total_var_8_10_low_NMAR_occ_no_tot[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_10; LOW MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_10_low_NMAR_occ_no_tot, nrows, ncolumns, numdatasets_low_NMAR_occ_no_tot)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_10_low_NMAR_occ_no_tot, nrows, ncolumns, numdatasets_low_NMAR_occ_no_tot)
}
##############################
# BIAS_8_10: LOW / "more totals known 10 iterations/"
##############################
numdatasets <- numdatasets_low_more_totals_known
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_more_totals_known, text_bias_8_10)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_10_low_more_totals_known[,,] <- bias_var_list$total_bias
  total_var_8_10_low_more_totals_known[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_10; LOW MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_10_low_more_totals_known, nrows, ncolumns, numdatasets_low_more_totals_known)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_10_low_more_totals_known, nrows, ncolumns, numdatasets_low_more_totals_known)
}
##############################
# BIAS_8_10: HIGH / 1 ITER
##############################
numdatasets <- numdatasets_high_1_iter
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_1_iter, text_bias_8_10)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_10_high_1_iter[,,] <- bias_var_list$total_bias
  total_var_8_10_high_1_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_10; HIGH MISSINGNESS; 1 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_10_high_1_iter, nrows, ncolumns, numdatasets_high_1_iter)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_10_high_1_iter, nrows, ncolumns, numdatasets_high_1_iter)
}
##############################
# BIAS_8_10: HIGH / 5 ITER
##############################
numdatasets <- numdatasets_high_5_iter
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_5_iter, text_bias_8_10)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_10_high_5_iter[,,] <- bias_var_list$total_bias
  total_var_8_10_high_5_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_10; HIGH MISSINGNESS; 5 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_10_high_5_iter, nrows, ncolumns, numdatasets_high_5_iter)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_10_high_5_iter, nrows, ncolumns, numdatasets_high_5_iter)
}
##############################
# BIAS_8_10: HIGH / 10 ITER
##############################
numdatasets <- numdatasets_high_10_iter
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_10_iter, text_bias_8_10)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_10_high_10_iter[,,] <- bias_var_list$total_bias
  total_var_8_10_high_10_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_10; HIGH MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_10_high_10_iter, nrows, ncolumns, numdatasets_high_10_iter)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_10_high_10_iter, nrows, ncolumns, numdatasets_high_10_iter)
}
##############################
# BIAS_8_10: HIGH / 20 ITER
##############################
numdatasets <- numdatasets_high_20_iter
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_20_iter, text_bias_8_10)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_10_high_20_iter[,,] <- bias_var_list$total_bias
  total_var_8_10_high_20_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_10; HIGH MISSINGNESS; 20 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_10_high_20_iter, nrows, ncolumns, numdatasets_high_20_iter)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_10_high_20_iter, nrows, ncolumns, numdatasets_high_20_iter)
}
##############################
# BIAS_8_10: HIGH / "Selective MAR edu 10 iterations/"
##############################
numdatasets <- numdatasets_high_MAR_edu
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_MAR_edu, text_bias_8_10)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_10_high_MAR_edu[,,] <- bias_var_list$total_bias
  total_var_8_10_high_MAR_edu[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_10; HIGH MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_10_high_MAR_edu, nrows, ncolumns, numdatasets_high_MAR_edu)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_10_high_MAR_edu, nrows, ncolumns, numdatasets_high_MAR_edu)
}
##############################
# BIAS_8_10: HIGH / "Selective MAR occ 10 iterations/"
##############################
numdatasets <- numdatasets_high_MAR_occ
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_MAR_occ, text_bias_8_10)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_10_high_MAR_occ[,,] <- bias_var_list$total_bias
  total_var_8_10_high_MAR_occ[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_10; HIGH MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_10_high_MAR_occ, nrows, ncolumns, numdatasets_high_MAR_occ)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_10_high_MAR_occ, nrows, ncolumns, numdatasets_high_MAR_occ)
}
##############################
# BIAS_8_10: HIGH / "Selective NMAR edu 10 iterations/"
##############################
numdatasets <- numdatasets_high_NMAR_edu
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_NMAR_edu, text_bias_8_10)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_10_high_NMAR_edu[,,] <- bias_var_list$total_bias
  total_var_8_10_high_NMAR_edu[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_10; HIGH MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_10_high_NMAR_edu, nrows, ncolumns, numdatasets_high_NMAR_edu)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_10_high_NMAR_edu, nrows, ncolumns, numdatasets_high_NMAR_edu)
}
##############################
# BIAS_8_10: HIGH / "Selective NMAR occ 10 iterations/"
##############################
numdatasets <- numdatasets_high_NMAR_occ
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_NMAR_occ, text_bias_8_10)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_10_high_NMAR_occ[,,] <- bias_var_list$total_bias
  total_var_8_10_high_NMAR_occ[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_10; HIGH MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_10_high_NMAR_occ, nrows, ncolumns, numdatasets_high_NMAR_occ)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_10_high_NMAR_occ, nrows, ncolumns, numdatasets_high_NMAR_occ)
}
##############################
# BIAS_8_10: HIGH / "Selective NMAR edu 10 iterations no totals/"
##############################
numdatasets <- numdatasets_high_NMAR_edu_no_tot
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_NMAR_edu_no_tot, text_bias_8_10)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_10_high_NMAR_edu_no_tot[,,] <- bias_var_list$total_bias
  total_var_8_10_high_NMAR_edu_no_tot[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_10; HIGH MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_10_high_NMAR_edu_no_tot, nrows, ncolumns, numdatasets_high_NMAR_edu_no_tot)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_10_high_NMAR_edu_no_tot, nrows, ncolumns, numdatasets_high_NMAR_edu_no_tot)
}
##############################
# BIAS_8_10: HIGH / "Selective NMAR occ 10 iterations no totals/"
##############################
numdatasets <- numdatasets_high_NMAR_occ_no_tot
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_NMAR_occ_no_tot, text_bias_8_10)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_10_high_NMAR_occ_no_tot[,,] <- bias_var_list$total_bias
  total_var_8_10_high_NMAR_occ_no_tot[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_10; HIGH MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_10_high_NMAR_occ_no_tot, nrows, ncolumns, numdatasets_high_NMAR_occ_no_tot)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_10_high_NMAR_occ_no_tot, nrows, ncolumns, numdatasets_high_NMAR_occ_no_tot)
}
##############################
# BIAS_8_10: HIGH / "more totals known 10 iterations/"
##############################
numdatasets <- numdatasets_high_more_totals_known
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_more_totals_known, text_bias_8_10)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_10_high_more_totals_known[,,] <- bias_var_list$total_bias
  total_var_8_10_high_more_totals_known[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_10; HIGH MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_10_high_more_totals_known, nrows, ncolumns, numdatasets_high_more_totals_known)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_10_high_more_totals_known, nrows, ncolumns, numdatasets_high_more_totals_known)
}
##############################

##############################
# CROSS-TABLE 8_2:
##############################
nrows <- 7
ncolumns <- 17
##############################
# BIAS_8_2: TRUETOTALS
##############################
datasetpath <- paste0(Low_path, text_1_iter, text_true_8_2, '.csv')
truetotals_8_2[,] <- Find_totals(datasetpath, nrows, ncolumns)
##############################
# BIAS_8_2: LOW / 1 ITER
##############################
numdatasets <- numdatasets_low_1_iter
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_1_iter, text_bias_8_2)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_2_low_1_iter[,,] <- bias_var_list$total_bias
  total_var_8_2_low_1_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_2; LOW MISSINGNESS; 1 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_2_low_1_iter, nrows, ncolumns, numdatasets_low_1_iter)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE, numdatasets_low_1_iter)
  Write_results(resultsfile, total_var_8_2_low_1_iter, nrows, ncolumns, numdatasets_low_1_iter)
}
##############################
# BIAS_8_2: LOW / 5 ITER
##############################
numdatasets <- numdatasets_low_5_iter
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_5_iter, text_bias_8_2)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_2_low_5_iter[,,] <- bias_var_list$total_bias
  total_var_8_2_low_5_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_2; LOW MISSINGNESS; 5 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_2_low_5_iter, nrows, ncolumns, numdatasets_low_5_iter)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_2_low_5_iter, nrows, ncolumns, numdatasets_low_5_iter)
}
##############################
# BIAS_8_2: LOW / 10 ITER
##############################
numdatasets <- numdatasets_low_10_iter
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_10_iter, text_bias_8_2)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_2_low_10_iter[,,] <- bias_var_list$total_bias
  total_var_8_2_low_10_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_2; LOW MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_2_low_10_iter, nrows, ncolumns, numdatasets_low_10_iter)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_2_low_10_iter, nrows, ncolumns, numdatasets_low_10_iter)
}
##############################
# BIAS_8_2: LOW / 20 ITER
##############################
numdatasets <- numdatasets_low_20_iter
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_20_iter, text_bias_8_2)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_2_low_20_iter[,,] <- bias_var_list$total_bias
  total_var_8_2_low_20_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_2; LOW MISSINGNESS; 20 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_2_low_20_iter, nrows, ncolumns, numdatasets_low_20_iter)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_2_low_20_iter, nrows, ncolumns, numdatasets_low_20_iter)
}
##############################
# BIAS_8_2: LOW / "Selective MAR edu 10 iterations/"
##############################
numdatasets <- numdatasets_low_MAR_edu
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_MAR_edu, text_bias_8_2)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_2_low_MAR_edu[,,] <- bias_var_list$total_bias
  total_var_8_2_low_MAR_edu[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_2; LOW MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_2_low_MAR_edu, nrows, ncolumns, numdatasets_low_MAR_edu)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_2_low_MAR_edu, nrows, ncolumns, numdatasets_low_MAR_edu)
}
##############################
# BIAS_8_2: LOW / "Selective MAR occ 10 iterations/"
##############################
numdatasets <- numdatasets_low_MAR_occ
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_MAR_occ, text_bias_8_2)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_2_low_MAR_occ[,,] <- bias_var_list$total_bias
  total_var_8_2_low_MAR_occ[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_2; LOW MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_2_low_MAR_occ, nrows, ncolumns, numdatasets_low_MAR_occ)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_2_low_MAR_occ, nrows, ncolumns, numdatasets_low_MAR_occ)
}
##############################
# BIAS_8_2: LOW / "Selective NMAR edu 10 iterations/"
##############################
numdatasets <- numdatasets_low_NMAR_edu
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_NMAR_edu, text_bias_8_2)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_2_low_NMAR_edu[,,] <- bias_var_list$total_bias
  total_var_8_2_low_NMAR_edu[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_2; LOW MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_2_low_NMAR_edu, nrows, ncolumns, numdatasets_low_NMAR_edu)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_2_low_NMAR_edu, nrows, ncolumns, numdatasets_low_NMAR_edu)
}
##############################
# BIAS_8_2: LOW / "Selective NMAR occ 10 iterations/"
##############################
numdatasets <- numdatasets_low_NMAR_occ
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_NMAR_occ, text_bias_8_2)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_2_low_NMAR_occ[,,] <- bias_var_list$total_bias
  total_var_8_2_low_NMAR_occ[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_2; LOW MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_2_low_NMAR_occ, nrows, ncolumns, numdatasets_low_NMAR_occ)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_2_low_NMAR_occ, nrows, ncolumns, numdatasets_low_NMAR_occ)
}
##############################
# BIAS_8_2: LOW / "Selective NMAR edu 10 iterations no totals/"
##############################
numdatasets <- numdatasets_low_NMAR_edu_no_tot
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_NMAR_edu_no_tot, text_bias_8_2)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_2_low_NMAR_edu_no_tot[,,] <- bias_var_list$total_bias
  total_var_8_2_low_NMAR_edu_no_tot[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_2; LOW MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_2_low_NMAR_edu_no_tot, nrows, ncolumns, numdatasets_low_NMAR_edu_no_tot)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_2_low_NMAR_edu_no_tot, nrows, ncolumns, numdatasets_low_NMAR_edu_no_tot)
}
##############################
# BIAS_8_2: LOW / "Selective NMAR occ 10 iterations no totals/"
##############################
numdatasets <- numdatasets_low_NMAR_occ_no_tot
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_NMAR_occ_no_tot, text_bias_8_2)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_2_low_NMAR_occ_no_tot[,,] <- bias_var_list$total_bias
  total_var_8_2_low_NMAR_occ_no_tot[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_2; LOW MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_2_low_NMAR_occ_no_tot, nrows, ncolumns, numdatasets_low_NMAR_occ_no_tot)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_2_low_NMAR_occ_no_tot, nrows, ncolumns, numdatasets_low_NMAR_occ_no_tot)
}
##############################
# BIAS_8_2: LOW / "more totals known 10 iterations/"
##############################
numdatasets <- numdatasets_low_more_totals_known
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_more_totals_known, text_bias_8_2)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_2_low_more_totals_known[,,] <- bias_var_list$total_bias
  total_var_8_2_low_more_totals_known[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_2; LOW MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_2_low_more_totals_known, nrows, ncolumns, numdatasets_low_more_totals_known)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_2_low_more_totals_known, nrows, ncolumns, numdatasets_low_more_totals_known)
}
##############################
# BIAS_8_2: HIGH / 1 ITER
##############################
numdatasets <- numdatasets_high_1_iter
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_1_iter, text_bias_8_2)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_2_high_1_iter[,,] <- bias_var_list$total_bias
  total_var_8_2_high_1_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_2; HIGH MISSINGNESS; 1 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_2_high_1_iter, nrows, ncolumns, numdatasets_high_1_iter)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_2_high_1_iter, nrows, ncolumns, numdatasets_high_1_iter)
}
##############################
# BIAS_8_2: HIGH / 5 ITER
##############################
numdatasets <- numdatasets_high_5_iter
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_5_iter, text_bias_8_2)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_2_high_5_iter[,,] <- bias_var_list$total_bias
  total_var_8_2_high_5_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_2; HIGH MISSINGNESS; 5 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_2_high_5_iter, nrows, ncolumns, numdatasets_high_5_iter)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_2_high_5_iter, nrows, ncolumns, numdatasets_high_5_iter)
}
##############################
# BIAS_8_2: HIGH / 10 ITER
##############################
numdatasets <- numdatasets_high_10_iter
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_10_iter, text_bias_8_2)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_2_high_10_iter[,,] <- bias_var_list$total_bias
  total_var_8_2_high_10_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_2; HIGH MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_2_high_10_iter, nrows, ncolumns, numdatasets_high_10_iter)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_2_high_10_iter, nrows, ncolumns, numdatasets_high_10_iter)
}
##############################
# BIAS_8_2: HIGH / 20 ITER
##############################
numdatasets <- numdatasets_high_MAR_edu
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_20_iter, text_bias_8_2)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_2_high_20_iter[,,] <- bias_var_list$total_bias
  total_var_8_2_high_20_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_2; HIGH MISSINGNESS; 20 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_2_high_20_iter, nrows, ncolumns, numdatasets_high_20_iter)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_2_high_20_iter, nrows, ncolumns, numdatasets_high_20_iter)
}
##############################
# BIAS_8_2: HIGH / "Selective MAR edu 10 iterations/"
##############################
numdatasets <- numdatasets_high_MAR_edu
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_MAR_edu, text_bias_8_2)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_2_high_MAR_edu[,,] <- bias_var_list$total_bias
  total_var_8_2_high_MAR_edu[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_2; HIGH MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_2_high_MAR_edu, nrows, ncolumns, numdatasets_high_MAR_edu)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_2_high_MAR_edu, nrows, ncolumns, numdatasets_high_MAR_edu)
}
##############################
# BIAS_8_2: HIGH / "Selective MAR occ 10 iterations/"
##############################
numdatasets <- numdatasets_high_MAR_occ
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_MAR_occ, text_bias_8_2)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_2_high_MAR_occ[,,] <- bias_var_list$total_bias
  total_var_8_2_high_MAR_occ[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_2; HIGH MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_2_high_MAR_occ, nrows, ncolumns, numdatasets_high_MAR_occ)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_2_high_MAR_occ, nrows, ncolumns, numdatasets_high_MAR_occ)
}
##############################
# BIAS_8_2: HIGH / "Selective NMAR edu 10 iterations/"
##############################
numdatasets <- numdatasets_high_NMAR_edu
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_NMAR_edu, text_bias_8_2)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_2_high_NMAR_edu[,,] <- bias_var_list$total_bias
  total_var_8_2_high_NMAR_edu[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_2; HIGH MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_2_high_NMAR_edu, nrows, ncolumns, numdatasets_high_NMAR_edu)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_2_high_NMAR_edu, nrows, ncolumns, numdatasets_high_NMAR_edu)
}
##############################
# BIAS_8_2: HIGH / "Selective NMAR occ 10 iterations/"
##############################
numdatasets <- numdatasets_high_NMAR_occ
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_NMAR_occ, text_bias_8_2)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_2_high_NMAR_occ[,,] <- bias_var_list$total_bias
  total_var_8_2_high_NMAR_occ[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_2; HIGH MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_2_high_NMAR_occ, nrows, ncolumns, numdatasets_high_NMAR_occ)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_2_high_NMAR_occ, nrows, ncolumns, numdatasets_high_NMAR_occ)
}
##############################
# BIAS_8_2: HIGH / "Selective NMAR edu 10 iterations no totals/"
##############################
numdatasets <- numdatasets_high_NMAR_edu_no_tot
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_NMAR_edu_no_tot, text_bias_8_2)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_2_high_NMAR_edu_no_tot[,,] <- bias_var_list$total_bias
  total_var_8_2_high_NMAR_edu_no_tot[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_2; HIGH MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_2_high_NMAR_edu_no_tot, nrows, ncolumns, numdatasets_high_NMAR_edu_no_tot)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_2_high_NMAR_edu_no_tot, nrows, ncolumns, numdatasets_high_NMAR_edu_no_tot)
}
##############################
# BIAS_8_2: HIGH / "Selective NMAR occ 10 iterations no totals/"
##############################
numdatasets <- numdatasets_high_NMAR_occ_no_tot
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_NMAR_occ_no_tot, text_bias_8_2)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_2_high_NMAR_occ_no_tot[,,] <- bias_var_list$total_bias
  total_var_8_2_high_NMAR_occ_no_tot[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_2; HIGH MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_2_high_NMAR_occ_no_tot, nrows, ncolumns, numdatasets_high_NMAR_occ_no_tot)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_2_high_NMAR_occ_no_tot, nrows, ncolumns, numdatasets_high_NMAR_occ_no_tot)
}
##############################
# BIAS_8_2: HIGH / "more totals known 10 iterations/"
##############################
numdatasets <- numdatasets_high_more_totals_known
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_more_totals_known, text_bias_8_2)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_2_high_more_totals_known[,,] <- bias_var_list$total_bias
  total_var_8_2_high_more_totals_known[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_2; HIGH MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_2_high_more_totals_known, nrows, ncolumns, numdatasets_high_more_totals_known)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_2_high_more_totals_known, nrows, ncolumns, numdatasets_high_more_totals_known)
}
##############################

##############################
# CROSS-TABLE 8_9:
##############################
nrows <- 7
ncolumns <- 8
##############################
# BIAS__8_9: TRUETOTALS
##############################
datasetpath <- paste0(Low_path, text_1_iter, text_true_8_9, '.csv')
truetotals_8_9[,] <- Find_totals(datasetpath, nrows, ncolumns)
##############################
# BIAS_8_9: LOW / 1 ITER
##############################
numdatasets <- numdatasets_low_1_iter
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_1_iter, text_bias_8_9)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_9_low_1_iter[,,] <- bias_var_list$total_bias
  total_var_8_9_low_1_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_9; LOW MISSINGNESS; 1 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_9_low_1_iter, nrows, ncolumns, numdatasets_low_1_iter)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_9_low_1_iter, nrows, ncolumns, numdatasets_low_1_iter)
}
##############################
# BIAS__8_9: LOW / 5 ITER
##############################
numdatasets <- numdatasets_low_5_iter
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_5_iter, text_bias_8_9)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_9_low_5_iter[,,] <- bias_var_list$total_bias
  total_var_8_9_low_5_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_9; LOW MISSINGNESS; 5 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_9_low_5_iter, nrows, ncolumns, numdatasets_low_5_iter)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_9_low_5_iter, nrows, ncolumns, numdatasets_low_5_iter)
}
##############################
# BIAS_8_9: LOW / 10 ITER
##############################
numdatasets <- numdatasets_low_10_iter
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_10_iter, text_bias_8_9)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_9_low_10_iter[,,] <- bias_var_list$total_bias
  total_var_8_9_low_10_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_9; LOW MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_9_low_10_iter, nrows, ncolumns, numdatasets_low_10_iter)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_9_low_10_iter, nrows, ncolumns, numdatasets_low_10_iter)
}
##############################
# BIAS_8_9: LOW / 20 ITER
##############################
numdatasets <- numdatasets_low_20_iter
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_20_iter, text_bias_8_9)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_9_low_20_iter[,,] <- bias_var_list$total_bias
  total_var_8_9_low_20_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_9; LOW MISSINGNESS; 20 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_9_low_20_iter, nrows, ncolumns, numdatasets_low_20_iter)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_9_low_20_iter, nrows, ncolumns, numdatasets_low_20_iter)
}
##############################
# BIAS_8_9: LOW / "Selective MAR edu 10 iterations/"
##############################
numdatasets <- numdatasets_low_MAR_edu
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_MAR_edu, text_bias_8_9)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_9_low_MAR_edu[,,] <- bias_var_list$total_bias
  total_var_8_9_low_MAR_edu[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_9; LOW MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_9_low_MAR_edu, nrows, ncolumns, numdatasets_low_MAR_edu)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_9_low_MAR_edu, nrows, ncolumns, numdatasets_low_MAR_edu)
}
##############################
# BIAS_8_9: LOW / "Selective MAR occ 10 iterations/"
##############################
numdatasets <- numdatasets_low_MAR_occ
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_MAR_occ, text_bias_8_9)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_9_low_MAR_occ[,,] <- bias_var_list$total_bias
  total_var_8_9_low_MAR_occ[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_9; LOW MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_9_low_MAR_occ, nrows, ncolumns, numdatasets_low_MAR_occ)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_9_low_MAR_occ, nrows, ncolumns, numdatasets_low_MAR_occ)
}
##############################
# BIAS_8_9: LOW / "Selective NMAR edu 10 iterations/"
##############################
numdatasets <- numdatasets_low_NMAR_edu
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_NMAR_edu, text_bias_8_9)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_9_low_NMAR_edu[,,] <- bias_var_list$total_bias
  total_var_8_9_low_NMAR_edu[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_9; LOW MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_9_low_NMAR_edu, nrows, ncolumns, numdatasets_low_NMAR_edu)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_9_low_NMAR_edu, nrows, ncolumns, numdatasets_low_NMAR_edu)
}
##############################
# BIAS_8_9: LOW / "Selective NMAR occ 10 iterations/"
##############################
numdatasets <- numdatasets_low_NMAR_occ
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_NMAR_occ, text_bias_8_9)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_9_low_NMAR_occ[,,] <- bias_var_list$total_bias
  total_var_8_9_low_NMAR_occ[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_9; LOW MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_9_low_NMAR_occ, nrows, ncolumns, numdatasets_low_NMAR_occ)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_9_low_NMAR_occ, nrows, ncolumns, numdatasets_low_NMAR_occ)
}
##############################
# BIAS_8_9: low / "Selective NMAR edu 10 iterations no totals/"
##############################
numdatasets <- numdatasets_low_NMAR_edu_no_tot
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_NMAR_edu_no_tot, text_bias_8_9)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_9_low_NMAR_edu_no_tot[,,] <- bias_var_list$total_bias
  total_var_8_9_low_NMAR_edu_no_tot[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_9; low MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_9_low_NMAR_edu_no_tot, nrows, ncolumns, numdatasets_low_NMAR_edu_no_tot)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_9_low_NMAR_edu_no_tot, nrows, ncolumns, numdatasets_low_NMAR_edu_no_tot)
}
##############################
# BIAS_8_9: LOW / "Selective NMAR occ 10 iterations no totals/"
##############################
numdatasets <- numdatasets_low_NMAR_occ_no_tot
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_NMAR_occ_no_tot, text_bias_8_9)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_9_low_NMAR_occ_no_tot[,,] <- bias_var_list$total_bias
  total_var_8_9_low_NMAR_occ_no_tot[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_9; LOW MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_9_low_NMAR_occ_no_tot, nrows, ncolumns, numdatasets_low_NMAR_occ_no_tot)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_9_low_NMAR_occ_no_tot, nrows, ncolumns, numdatasets_low_NMAR_occ_no_tot)
}
##############################
# BIAS_8_9: LOW / "more totals known 10 iterations/"
##############################
numdatasets <- numdatasets_low_more_totals_known
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_more_totals_known, text_bias_8_9)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_9_low_more_totals_known[,,] <- bias_var_list$total_bias
  total_var_8_9_low_more_totals_known[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_9; LOW MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_9_low_more_totals_known, nrows, ncolumns, numdatasets_low_more_totals_known)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_9_low_more_totals_known, nrows, ncolumns, numdatasets_low_more_totals_known)
}
##############################
# BIAS_8_9: HIGH / 1 ITER
##############################
numdatasets <- numdatasets_high_1_iter
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_1_iter, text_bias_8_9)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_9_high_1_iter[,,] <- bias_var_list$total_bias
  total_var_8_9_high_1_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_9; HIGH MISSINGNESS; 1 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_9_high_1_iter, nrows, ncolumns, numdatasets_high_1_iter)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_9_high_1_iter, nrows, ncolumns, numdatasets_high_1_iter)
}
##############################
# BIAS__8_9: HIGH / 5 ITER
##############################
numdatasets <- numdatasets_high_5_iter
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_5_iter, text_bias_8_9)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_9_high_5_iter[,,] <- bias_var_list$total_bias
  total_var_8_9_high_5_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_9; HIGH MISSINGNESS; 5 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_9_high_5_iter, nrows, ncolumns, numdatasets_high_5_iter)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_9_high_5_iter, nrows, ncolumns, numdatasets_high_5_iter)
}
##############################
# BIAS_8_9: HIGH / 10 ITER
##############################
numdatasets <- numdatasets_high_10_iter
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_10_iter, text_bias_8_9)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_9_high_10_iter[,,] <- bias_var_list$total_bias
  total_var_8_9_high_10_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_9; HIGH MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_9_high_10_iter, nrows, ncolumns, numdatasets_high_10_iter)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_9_high_10_iter, nrows, ncolumns, numdatasets_high_10_iter)
}
##############################
# BIAS_8_9: HIGH / 20 ITER
##############################
numdatasets <- numdatasets_high_20_iter
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_20_iter, text_bias_8_9)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_9_high_20_iter[,,] <- bias_var_list$total_bias
  total_var_8_9_high_20_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_9; HIGH MISSINGNESS; 20 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_9_high_20_iter, nrows, ncolumns, numdatasets_high_20_iter)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_9_high_20_iter, nrows, ncolumns, numdatasets_high_20_iter)
}
##############################
# BIAS_8_9: HIGH / "Selective MAR edu 10 iterations/"
##############################
numdatasets <- numdatasets_high_MAR_edu
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_MAR_edu, text_bias_8_9)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_9_high_MAR_edu[,,] <- bias_var_list$total_bias
  total_var_8_9_high_MAR_edu[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_9; HIGH MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_9_high_MAR_edu, nrows, ncolumns, numdatasets_high_MAR_edu)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_9_high_MAR_edu, nrows, ncolumns, numdatasets_high_MAR_edu)
}
##############################
# BIAS_8_9: HIGH / "Selective MAR occ 10 iterations/"
##############################
numdatasets <- numdatasets_high_MAR_occ
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_MAR_occ, text_bias_8_9)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_9_high_MAR_occ[,,] <- bias_var_list$total_bias
  total_var_8_9_high_MAR_occ[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_9; HIGH MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_9_high_MAR_occ, nrows, ncolumns, numdatasets_high_MAR_occ)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_9_high_MAR_occ, nrows, ncolumns, numdatasets_high_MAR_occ)
}
##############################
# BIAS_8_9: HIGH / "Selective NMAR edu 10 iterations/"
##############################
numdatasets <- numdatasets_high_NMAR_edu
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_NMAR_edu, text_bias_8_9)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_9_high_NMAR_edu[,,] <- bias_var_list$total_bias
  total_var_8_9_high_NMAR_edu[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_9; HIGH MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_9_high_NMAR_edu, nrows, ncolumns, numdatasets_high_NMAR_edu)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_9_high_NMAR_edu, nrows, ncolumns, numdatasets_high_NMAR_edu)
}
##############################
# BIAS_8_9: HIGH / "Selective NMAR occ 10 iterations/"
##############################
numdatasets <- numdatasets_high_NMAR_occ
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_NMAR_occ, text_bias_8_9)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_9_high_NMAR_occ[,,] <- bias_var_list$total_bias
  total_var_8_9_high_NMAR_occ[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_9; HIGH MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_9_high_NMAR_occ, nrows, ncolumns, numdatasets_high_NMAR_occ)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_9_high_NMAR_occ, nrows, ncolumns, numdatasets_high_NMAR_occ)
}
##############################
# BIAS_8_9: HIGH / "Selective NMAR edu 10 iterations no totals/"
##############################
numdatasets <- numdatasets_high_NMAR_edu_no_tot
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_NMAR_edu_no_tot, text_bias_8_9)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_9_high_NMAR_edu_no_tot[,,] <- bias_var_list$total_bias
  total_var_8_9_high_NMAR_edu_no_tot[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_9; HIGH MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_9_high_NMAR_edu_no_tot, nrows, ncolumns, numdatasets_high_NMAR_edu_no_tot)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_9_high_NMAR_edu_no_tot, nrows, ncolumns, numdatasets_high_NMAR_edu_no_tot)
}
##############################
# BIAS_8_9: HIGH / "Selective NMAR occ 10 iterations no totals/"
##############################
numdatasets <- numdatasets_high_NMAR_occ_no_tot
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_NMAR_occ_no_tot, text_bias_8_9)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_9_high_NMAR_occ_no_tot[,,] <- bias_var_list$total_bias
  total_var_8_9_high_NMAR_occ_no_tot[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_9; HIGH MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_9_high_NMAR_occ_no_tot, nrows, ncolumns, numdatasets_high_NMAR_occ_no_tot)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_9_high_NMAR_occ_no_tot, nrows, ncolumns, numdatasets_high_NMAR_occ_no_tot)
}
##############################
# BIAS_8_9: HIGH / "more totals known 10 iterations/"
##############################
numdatasets <- numdatasets_high_more_totals_known
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_more_totals_known, text_bias_8_9)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_8_9_high_more_totals_known[,,] <- bias_var_list$total_bias
  total_var_8_9_high_more_totals_known[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 8_9; HIGH MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_8_9_high_more_totals_known, nrows, ncolumns, numdatasets_high_more_totals_known)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_8_9_high_more_totals_known, nrows, ncolumns, numdatasets_high_more_totals_known)
}
##############################

##############################
# CROSS-TABLE 10_2:
##############################
nrows <- 10
ncolumns <- 17
##############################
# BIAS_10_2: TRUETOTALS
##############################
datasetpath <- paste0(Low_path, text_1_iter, text_true_10_2, '.csv')
truetotals_10_2[,] <- Find_totals(datasetpath, nrows, ncolumns)
##############################
# BIAS_10_2: LOW / 1 ITER
##############################
numdatasets <- numdatasets_low_1_iter
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_1_iter, text_bias_10_2)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_10_2_low_1_iter[,,] <- bias_var_list$total_bias
  total_var_10_2_low_1_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 10_2; LOW MISSINGNESS; 1 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_10_2_low_1_iter, nrows, ncolumns, numdatasets_low_1_iter)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_10_2_low_1_iter, nrows, ncolumns, numdatasets_low_1_iter)
}
##############################
# BIAS_10_2: LOW / 5 ITER
##############################
numdatasets <- numdatasets_low_5_iter
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_5_iter, text_bias_10_2)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_10_2_low_5_iter[,,] <- bias_var_list$total_bias
  total_var_10_2_low_5_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 10_2; LOW MISSINGNESS; 5 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_10_2_low_5_iter, nrows, ncolumns, numdatasets_low_5_iter)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_10_2_low_5_iter, nrows, ncolumns, numdatasets_low_5_iter)
}
##############################
# BIAS_10_2: LOW / 10 ITER
##############################
numdatasets <- numdatasets_low_10_iter
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_10_iter, text_bias_10_2)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_10_2_low_10_iter[,,] <- bias_var_list$total_bias
  total_var_10_2_low_10_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 10_2; LOW MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_10_2_low_10_iter, nrows, ncolumns, numdatasets_low_10_iter)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_10_2_low_10_iter, nrows, ncolumns, numdatasets_low_10_iter)
}
##############################
# BIAS_10_2: LOW / 20 ITER
##############################
numdatasets <- numdatasets_low_20_iter
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_20_iter, text_bias_10_2)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_10_2_low_20_iter[,,] <- bias_var_list$total_bias
  total_var_10_2_low_20_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 10_2; LOW MISSINGNESS; 20 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_10_2_low_20_iter, nrows, ncolumns, numdatasets_low_20_iter)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_10_2_low_20_iter, nrows, ncolumns, numdatasets_low_20_iter)
}
##############################
# BIAS_10_2: LOW / "Selective MAR edu 10 iterations/"
##############################
numdatasets <- numdatasets_low_MAR_edu
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_MAR_edu, text_bias_10_2)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_10_2_low_MAR_edu[,,] <- bias_var_list$total_bias
  total_var_10_2_low_MAR_edu[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 10_2; LOW MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_10_2_low_MAR_edu, nrows, ncolumns, numdatasets_low_MAR_edu)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_10_2_low_MAR_edu, nrows, ncolumns, numdatasets_low_MAR_edu)
}
##############################
# BIAS_10_2: LOW / "Selective MAR occ 10 iterations/"
##############################
numdatasets <- numdatasets_low_MAR_occ
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_MAR_occ, text_bias_10_2)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_10_2_low_MAR_occ[,,] <- bias_var_list$total_bias
  total_var_10_2_low_MAR_occ[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 10_2; LOW MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_10_2_low_MAR_occ, nrows, ncolumns, numdatasets_low_MAR_occ)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_10_2_low_MAR_occ, nrows, ncolumns, numdatasets_low_MAR_occ)
}
##############################
# BIAS_10_2: LOW / "Selective NMAR edu 10 iterations/"
##############################
numdatasets <- numdatasets_low_NMAR_edu
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_NMAR_edu, text_bias_10_2)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_10_2_low_NMAR_edu[,,] <- bias_var_list$total_bias
  total_var_10_2_low_NMAR_edu[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 10_2; LOW MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_10_2_low_NMAR_edu, nrows, ncolumns, numdatasets_low_NMAR_edu)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_10_2_low_NMAR_edu, nrows, ncolumns, numdatasets_low_NMAR_edu)
}
##############################
# BIAS_10_2: LOW / "Selective NMAR occ 10 iterations/"
##############################
numdatasets <- numdatasets_low_NMAR_occ
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_NMAR_occ, text_bias_10_2)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_10_2_low_NMAR_occ[,,] <- bias_var_list$total_bias
  total_var_10_2_low_NMAR_occ[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 10_2; LOW MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_10_2_low_NMAR_occ, nrows, ncolumns, numdatasets_low_NMAR_occ)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_10_2_low_NMAR_occ, nrows, ncolumns, numdatasets_low_NMAR_occ)
}
##############################
# BIAS_10_2: LOW / "Selective NMAR edu 10 iterations no totals/"
##############################
numdatasets <- numdatasets_low_NMAR_edu_no_tot
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_NMAR_edu_no_tot, text_bias_10_2)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_10_2_low_NMAR_edu_no_tot[,,] <- bias_var_list$total_bias
  total_var_10_2_low_NMAR_edu_no_tot[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 10_2; LOW MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_10_2_low_NMAR_edu_no_tot, nrows, ncolumns, numdatasets_low_NMAR_edu_no_tot)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_10_2_low_NMAR_edu_no_tot, nrows, ncolumns, numdatasets_low_NMAR_edu_no_tot)
}
##############################
# BIAS_10_2: LOW / "Selective NMAR occ 10 iterations no totals/"
##############################
numdatasets <- numdatasets_low_NMAR_occ_no_tot
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_NMAR_occ_no_tot, text_bias_10_2)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_10_2_low_NMAR_occ_no_tot[,,] <- bias_var_list$total_bias
  total_var_10_2_low_NMAR_occ_no_tot[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 10_2; LOW MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_10_2_low_NMAR_occ_no_tot, nrows, ncolumns, numdatasets_low_NMAR_occ_no_tot)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_10_2_low_NMAR_occ_no_tot, nrows, ncolumns, numdatasets_low_NMAR_occ_no_tot)
}
##############################
# BIAS_10_2: LOW / "more totals known 10 iterations/"
##############################
numdatasets <- numdatasets_low_more_totals_known
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_more_totals_known, text_bias_10_2)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_10_2_low_more_totals_known[,,] <- bias_var_list$total_bias
  total_var_10_2_low_more_totals_known[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 10_2; LOW MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_10_2_low_more_totals_known, nrows, ncolumns, numdatasets_low_more_totals_known)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_10_2_low_more_totals_known, nrows, ncolumns, numdatasets_low_more_totals_known)
}
##############################
# BIAS_10_2: HIGH / 1 ITER
##############################
numdatasets <- numdatasets_high_1_iter
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_1_iter, text_bias_10_2)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_10_2_high_1_iter[,,] <- bias_var_list$total_bias
  total_var_10_2_high_1_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 10_2; HIGH MISSINGNESS; 1 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_10_2_high_1_iter, nrows, ncolumns, numdatasets_high_1_iter)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_10_2_high_1_iter, nrows, ncolumns, numdatasets_high_1_iter)
}
##############################
# BIAS_10_2: HIGH / 5 ITER
##############################
numdatasets <- numdatasets_high_5_iter
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_5_iter, text_bias_10_2)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_10_2_high_5_iter[,,] <- bias_var_list$total_bias
  total_var_10_2_high_5_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 10_2; HIGH MISSINGNESS; 5 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_10_2_high_5_iter, nrows, ncolumns, numdatasets_high_5_iter)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_10_2_high_5_iter, nrows, ncolumns, numdatasets_high_5_iter)
}
##############################
# BIAS_10_2: HIGH / 10 ITER
##############################
numdatasets <- numdatasets_high_10_iter
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_10_iter, text_bias_10_2)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_10_2_high_10_iter[,,] <- bias_var_list$total_bias
  total_var_10_2_high_10_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 10_2; HIGH MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_10_2_high_10_iter, nrows, ncolumns, numdatasets_high_10_iter)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_10_2_high_10_iter, nrows, ncolumns, numdatasets_high_10_iter)
}
##############################
# BIAS_10_2: HIGH / 20 ITER
##############################
numdatasets <- numdatasets_high_20_iter
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_20_iter, text_bias_10_2)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_10_2_high_20_iter[,,] <- bias_var_list$total_bias
  total_var_10_2_high_20_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 10_2; HIGH MISSINGNESS; 20 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_10_2_high_20_iter, nrows, ncolumns, numdatasets_high_20_iter)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_10_2_high_20_iter, nrows, ncolumns, numdatasets_high_20_iter)
}
##############################
# BIAS_10_2: HIGH / "Selective MAR edu 10 iterations/"
##############################
numdatasets <- numdatasets_high_MAR_edu
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_MAR_edu, text_bias_10_2)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_10_2_high_MAR_edu[,,] <- bias_var_list$total_bias
  total_var_10_2_high_MAR_edu[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 10_2; HIGH MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_10_2_high_MAR_edu, nrows, ncolumns, numdatasets_high_MAR_edu)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_10_2_high_MAR_edu, nrows, ncolumns, numdatasets_high_MAR_edu)
}
##############################
# BIAS_10_2: HIGH / "Selective MAR occ 10 iterations/"
##############################
numdatasets <- numdatasets_high_MAR_occ
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_MAR_occ, text_bias_10_2)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_10_2_high_MAR_occ[,,] <- bias_var_list$total_bias
  total_var_10_2_high_MAR_occ[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 10_2; HIGH MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_10_2_high_MAR_occ, nrows, ncolumns, numdatasets_high_MAR_occ)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_10_2_high_MAR_occ, nrows, ncolumns, numdatasets_high_MAR_occ)
}
##############################
# BIAS_10_2: HIGH / "Selective NMAR edu 10 iterations/"
##############################
numdatasets <- numdatasets_high_NMAR_edu
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_NMAR_edu, text_bias_10_2)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_10_2_high_NMAR_edu[,,] <- bias_var_list$total_bias
  total_var_10_2_high_NMAR_edu[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 10_2; HIGH MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_10_2_high_NMAR_edu, nrows, ncolumns, numdatasets_high_NMAR_edu)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_10_2_high_NMAR_edu, nrows, ncolumns, numdatasets_high_NMAR_edu)
}
##############################
# BIAS_10_2: HIGH / "Selective NMAR occ 10 iterations/"
##############################
numdatasets <- numdatasets_high_NMAR_occ
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_NMAR_occ, text_bias_10_2)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_10_2_high_NMAR_occ[,,] <- bias_var_list$total_bias
  total_var_10_2_high_NMAR_occ[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 10_2; HIGH MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_10_2_high_NMAR_occ, nrows, ncolumns, numdatasets_high_NMAR_occ)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_10_2_high_NMAR_occ, nrows, ncolumns, numdatasets_high_NMAR_occ)
}
##############################
# BIAS_10_2: HIGH / "Selective NMAR edu 10 iterations no totals/"
##############################
numdatasets <- numdatasets_high_NMAR_edu_no_tot
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_NMAR_edu_no_tot, text_bias_10_2)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_10_2_high_NMAR_edu_no_tot[,,] <- bias_var_list$total_bias
  total_var_10_2_high_NMAR_edu_no_tot[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 10_2; HIGH MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_10_2_high_NMAR_edu_no_tot, nrows, ncolumns, numdatasets_high_NMAR_edu_no_tot)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_10_2_high_NMAR_edu_no_tot, nrows, ncolumns, numdatasets_high_NMAR_edu_no_tot)
}
##############################
# BIAS_10_2: HIGH / "Selective NMAR occ 10 iterations no totals/"
##############################
numdatasets <- numdatasets_high_NMAR_occ_no_tot
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_NMAR_occ_no_tot, text_bias_10_2)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_10_2_high_NMAR_occ_no_tot[,,] <- bias_var_list$total_bias
  total_var_10_2_high_NMAR_occ_no_tot[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 10_2; HIGH MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_10_2_high_NMAR_occ_no_tot, nrows, ncolumns, numdatasets_high_NMAR_occ_no_tot)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_10_2_high_NMAR_occ_no_tot, nrows, ncolumns, numdatasets_high_NMAR_occ_no_tot)
}
##############################
# BIAS_10_2: HIGH / "more totals known 10 iterations/"
##############################
numdatasets <- numdatasets_high_more_totals_known
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_more_totals_known, text_bias_10_2)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_10_2_high_more_totals_known[,,] <- bias_var_list$total_bias
  total_var_10_2_high_more_totals_known[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 10_2; HIGH MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_10_2_high_more_totals_known, nrows, ncolumns, numdatasets_high_more_totals_known)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_10_2_high_more_totals_known, nrows, ncolumns, numdatasets_high_more_totals_known)
}
##############################

##############################
# CROSS-TABLE 10_9:
##############################
nrows <- 10
ncolumns <- 8
##############################
# BIAS_10_9: TRUETOTALS
##############################
datasetpath <- paste0(Low_path, text_1_iter, text_true_10_9, '.csv')
truetotals_10_9[,] <- Find_totals(datasetpath, nrows, ncolumns)
##############################
# BIAS_10_9: LOW / 1 ITER
##############################
numdatasets <- numdatasets_low_1_iter
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_1_iter, text_bias_10_9)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_10_9_low_1_iter[,,] <- bias_var_list$total_bias
  total_var_10_9_low_1_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 10_9; LOW MISSINGNESS; 1 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_10_9_low_1_iter, nrows, ncolumns, numdatasets_low_1_iter)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_10_9_low_1_iter, nrows, ncolumns, numdatasets_low_1_iter)
}
##############################
# BIAS_10_9: LOW / 5 ITER
##############################
numdatasets <- numdatasets_low_5_iter
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_5_iter, text_bias_10_9)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_10_9_low_5_iter[,,] <- bias_var_list$total_bias
  total_var_10_9_low_5_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 10_9; LOW MISSINGNESS; 5 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_10_9_low_5_iter, nrows, ncolumns, numdatasets_low_5_iter)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_10_9_low_5_iter, nrows, ncolumns, numdatasets_low_5_iter)
}
##############################
# BIAS_10_9: LOW / 10 ITER
##############################
numdatasets <- numdatasets_low_10_iter
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_10_iter, text_bias_10_9)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_10_9_low_10_iter[,,] <- bias_var_list$total_bias
  total_var_10_9_low_10_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 10_9; LOW MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_10_9_low_10_iter, nrows, ncolumns, numdatasets_low_10_iter)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_10_9_low_10_iter, nrows, ncolumns, numdatasets_low_10_iter)
}
##############################
# BIAS_10_9: LOW / 20 ITER
##############################
numdatasets <- numdatasets_low_20_iter
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_20_iter, text_bias_10_9)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_10_9_low_20_iter[,,] <- bias_var_list$total_bias
  total_var_10_9_low_20_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 10_9; LOW MISSINGNESS; 20 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_10_9_low_20_iter, nrows, ncolumns, numdatasets_low_20_iter)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_10_9_low_20_iter, nrows, ncolumns, numdatasets_low_20_iter)
}
##############################
# BIAS_10_9: LOW / "Selective MAR edu 10 iterations/"
##############################
numdatasets <- numdatasets_low_MAR_edu
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_MAR_edu, text_bias_10_9)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_10_9_low_MAR_edu[,,] <- bias_var_list$total_bias
  total_var_10_9_low_MAR_edu[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 10_9; LOW MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_10_9_low_MAR_edu, nrows, ncolumns, numdatasets_low_MAR_edu)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_10_9_low_MAR_edu, nrows, ncolumns, numdatasets_low_MAR_edu)
}
##############################
# BIAS_10_9: LOW / "Selective MAR occ 10 iterations/"
##############################
numdatasets <- numdatasets_low_MAR_occ
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_MAR_occ, text_bias_10_9)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_10_9_low_MAR_occ[,,] <- bias_var_list$total_bias
  total_var_10_9_low_MAR_occ[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 10_9; LOW MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_10_9_low_MAR_occ, nrows, ncolumns, numdatasets_low_MAR_occ)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_10_9_low_MAR_occ, nrows, ncolumns, numdatasets_low_MAR_occ)
}
##############################
# BIAS_10_9: LOW / "Selective NMAR edu 10 iterations/"
##############################
numdatasets <- numdatasets_low_NMAR_edu
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_NMAR_edu, text_bias_10_9)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_10_9_low_NMAR_edu[,,] <- bias_var_list$total_bias
  total_var_10_9_low_NMAR_edu[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 10_9; LOW MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_10_9_low_NMAR_edu, nrows, ncolumns, numdatasets_low_NMAR_edu)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_10_9_low_NMAR_edu, nrows, ncolumns, numdatasets_low_NMAR_edu)
}
##############################
# BIAS_10_9: LOW / "Selective NMAR occ 10 iterations/"
##############################
numdatasets <- numdatasets_low_NMAR_occ
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_NMAR_occ, text_bias_10_9)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_10_9_low_NMAR_occ[,,] <- bias_var_list$total_bias
  total_var_10_9_low_NMAR_occ[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 10_9; LOW MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_10_9_low_NMAR_occ, nrows, ncolumns, numdatasets_low_NMAR_occ)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_10_9_low_NMAR_occ, nrows, ncolumns, numdatasets_low_NMAR_occ)
}
##############################
# BIAS_10_9: LOW / "Selective NMAR edu 10 iterations no totals/"
##############################
numdatasets <- numdatasets_low_NMAR_edu_no_tot
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_NMAR_edu_no_tot, text_bias_10_9)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_10_9_low_NMAR_edu_no_tot[,,] <- bias_var_list$total_bias
  total_var_10_9_low_NMAR_edu_no_tot[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 10_9; LOW MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_10_9_low_NMAR_edu_no_tot, nrows, ncolumns, numdatasets_low_NMAR_edu_no_tot)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_10_9_low_NMAR_edu_no_tot, nrows, ncolumns, numdatasets_low_NMAR_edu_no_tot)
}
##############################
# BIAS_10_9: LOW / "Selective NMAR occ 10 iterations no totals/"
##############################
numdatasets <- numdatasets_low_NMAR_occ_no_tot
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_NMAR_occ_no_tot, text_bias_10_9)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_10_9_low_NMAR_occ_no_tot[,,] <- bias_var_list$total_bias
  total_var_10_9_low_NMAR_occ_no_tot[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 10_9; LOW MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_10_9_low_NMAR_occ_no_tot, nrows, ncolumns, numdatasets_low_NMAR_occ_no_tot)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_10_9_low_NMAR_occ_no_tot, nrows, ncolumns, numdatasets_low_NMAR_occ_no_tot)
}
##############################
# BIAS_10_9: LOW / "more totals known 10 iterations/"
##############################
numdatasets <- numdatasets_low_more_totals_known
if (numdatasets > 0) {
  datasetpath <- paste0(Low_path, text_more_totals_known, text_bias_10_9)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_10_9_low_more_totals_known[,,] <- bias_var_list$total_bias
  total_var_10_9_low_more_totals_known[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 10_9; LOW MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_10_9_low_more_totals_known, nrows, ncolumns, numdatasets_low_more_totals_known)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_10_9_low_more_totals_known, nrows, ncolumns, numdatasets_low_more_totals_known)
}
##############################
# BIAS_10_9: HIGH / 1 ITER
##############################
numdatasets <- numdatasets_high_1_iter
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_1_iter, text_bias_10_9)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_10_9_high_1_iter[,,] <- bias_var_list$total_bias
  total_var_10_9_high_1_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 10_9; HIGH MISSINGNESS; 1 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_10_9_high_1_iter, nrows, ncolumns, numdatasets_high_1_iter)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_10_9_high_1_iter, nrows, ncolumns, numdatasets_high_1_iter)
}
##############################
# BIAS_10_9: HIGH / 5 ITER
##############################
numdatasets <- numdatasets_high_5_iter
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_5_iter, text_bias_10_9)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_10_9_high_5_iter[,,] <- bias_var_list$total_bias
  total_var_10_9_high_5_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 10_9; HIGH MISSINGNESS; 5 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_10_9_high_5_iter, nrows, ncolumns, numdatasets_high_5_iter)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_10_9_high_5_iter, nrows, ncolumns, numdatasets_high_5_iter)
}
##############################
# BIAS_10_9: HIGH / 10 ITER
##############################
numdatasets <- numdatasets_high_10_iter
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_10_iter, text_bias_10_9)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_10_9_high_10_iter[,,] <- bias_var_list$total_bias
  total_var_10_9_high_10_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 10_9; HIGH MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_10_9_high_10_iter, nrows, ncolumns, numdatasets_high_10_iter)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_10_9_high_10_iter, nrows, ncolumns, numdatasets_high_10_iter)
}
##############################
# BIAS_10_9: HIGH / 20 ITER
##############################
numdatasets <- numdatasets_high_20_iter
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_20_iter, text_bias_10_9)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_10_9_high_20_iter[,,] <- bias_var_list$total_bias
  total_var_10_9_high_20_iter[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 10_9; HIGH MISSINGNESS; 20 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_10_9_high_20_iter, nrows, ncolumns, numdatasets_high_20_iter)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_10_9_high_20_iter, nrows, ncolumns, numdatasets_high_20_iter)
}
##############################
# BIAS_10_9: HIGH / "Selective MAR edu 10 iterations/"
##############################
numdatasets <- numdatasets_high_MAR_edu
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_MAR_edu, text_bias_10_9)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_10_9_high_MAR_edu[,,] <- bias_var_list$total_bias
  total_var_10_9_high_MAR_edu[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 10_9; HIGH MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_10_9_high_MAR_edu, nrows, ncolumns, numdatasets_high_MAR_edu)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_10_9_high_MAR_edu, nrows, ncolumns, numdatasets_high_MAR_edu)
}
##############################
# BIAS_10_9: HIGH / "Selective MAR occ 10 iterations/"
##############################
numdatasets <- numdatasets_high_MAR_occ
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_MAR_occ, text_bias_10_9)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_10_9_high_MAR_occ[,,] <- bias_var_list$total_bias
  total_var_10_9_high_MAR_occ[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 10_9; HIGH MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_10_9_high_MAR_occ, nrows, ncolumns, numdatasets_high_MAR_occ)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_10_9_high_MAR_occ, nrows, ncolumns, numdatasets_high_MAR_occ)
}
##############################
# BIAS_10_9: HIGH / "Selective NMAR edu 10 iterations/"
##############################
numdatasets <- numdatasets_high_NMAR_edu
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_NMAR_edu, text_bias_10_9)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_10_9_high_NMAR_edu[,,] <- bias_var_list$total_bias
  total_var_10_9_high_NMAR_edu[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 10_9; HIGH MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_10_9_high_NMAR_edu, nrows, ncolumns, numdatasets_high_NMAR_edu)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_10_9_high_NMAR_edu, nrows, ncolumns, numdatasets_high_NMAR_edu)
}
##############################
# BIAS_10_9: HIGH / "Selective NMAR occ 10 iterations/"
##############################
numdatasets <- numdatasets_high_NMAR_occ
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_NMAR_occ, text_bias_10_9)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_10_9_high_NMAR_occ[,,] <- bias_var_list$total_bias
  total_var_10_9_high_NMAR_occ[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 10_9; HIGH MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_10_9_high_NMAR_occ, nrows, ncolumns, numdatasets_high_NMAR_occ)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_10_9_high_NMAR_occ, nrows, ncolumns, numdatasets_high_NMAR_occ)
}
##############################
# BIAS_10_9: HIGH / "Selective NMAR edu 10 iterations no totals/"
##############################
numdatasets <- numdatasets_high_NMAR_edu_no_tot
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_NMAR_edu_no_tot, text_bias_10_9)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_10_9_high_NMAR_edu_no_tot[,,] <- bias_var_list$total_bias
  total_var_10_9_high_NMAR_edu_no_tot[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 10_9; HIGH MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_10_9_high_NMAR_edu_no_tot, nrows, ncolumns, numdatasets_high_NMAR_edu_no_tot)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_10_9_high_NMAR_edu_no_tot, nrows, ncolumns, numdatasets_high_NMAR_edu_no_tot)
}
##############################
# BIAS_10_9: HIGH / "Selective NMAR occ 10 iterations no totals/"
##############################
numdatasets <- numdatasets_high_NMAR_occ_no_tot
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_NMAR_occ_no_tot, text_bias_10_9)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_10_9_high_NMAR_occ_no_tot[,,] <- bias_var_list$total_bias
  total_var_10_9_high_NMAR_occ_no_tot[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 10_9; HIGH MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_10_9_high_NMAR_occ_no_tot, nrows, ncolumns, numdatasets_high_NMAR_occ_no_tot)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_10_9_high_NMAR_occ_no_tot, nrows, ncolumns, numdatasets_high_NMAR_occ_no_tot)
}
##############################
# BIAS_10_9: HIGH / "more totals known 10 iterations/"
##############################
numdatasets <- numdatasets_high_more_totals_known
if (numdatasets > 0) {
  datasetpath <- paste0(High_path, text_more_totals_known, text_bias_10_9)
  bias_var_list <- Calculate_bias_var_allsets(datasetpath, numdatasets, nrows, ncolumns)
  total_bias_10_9_high_more_totals_known[,,] <- bias_var_list$total_bias
  total_var_10_9_high_more_totals_known[,,] <- bias_var_list$total_var
  resultsfile <- paste0(datasetpath,"quality measures.csv")
  write.table(file=resultsfile, "CROSS-TABLE 10_9; HIGH MISSINGNESS; 10 ITER", append = FALSE, col.names = FALSE)
  write.table(file=resultsfile, "BIAS", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_bias_10_9_high_more_totals_known, nrows, ncolumns, numdatasets_high_more_totals_known)
  write.table(file=resultsfile, "MSE", append = TRUE, col.names = FALSE)
  Write_results(resultsfile, total_var_10_9_high_more_totals_known, nrows, ncolumns, numdatasets_high_more_totals_known)
}
##############################
