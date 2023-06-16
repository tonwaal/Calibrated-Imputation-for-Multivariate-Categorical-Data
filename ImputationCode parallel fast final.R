#############################################################################################################################
# CODE FOR MASS IMPUTATION OF CATEGORICAL DATA GIVEN KNOWN TOTALS 
# PROJECT: Multivariate Mass Imputation for the Population Census Given Known Totals (2018-NL-ESS.VIP.BUS.ADMIN)  
# CALL: Multipurpose statistics for efficiency gains in production
#
# THIS CODE IS PARALLELIZED IN ORDER TO SPEED UP THE SIMULATION STUDY
#############################################################################################################################

# Onder RDS-Werk:
Sys.setenv(PATH=paste0("Y:/3.5/R-3.5.2/bin;Y:/3.5/Rtools/bin;", Sys.getenv("PATH")))
Sys.setenv(BINPREF="Y:/3.5/Rtools/mingw_$(WIN)/bin/" )

library(foreach)
library(doParallel)

init <- Sys.time()

path <- "//CBSP.NL/Productie/Secundair/MPOnderzoek/Werk/Combineren/Medewerkers/TWAL/Eigen projecten/Calibrated Imputation/Categorical data/R scripts/"
setwd(path)

myCluster <- makeCluster(detectCores() - 1, type="PSOCK")
registerDoParallel(myCluster)

########################
# NUMBER OF MISSING DATASETS FOR THIS PART OF SIMULATION STUDY
########################
NumMissingSets <- 250

########################
# LOOP OVER MISSING DATASETS
########################
foreach(s = (1:NumMissingSets)) %dopar% {  
  
  library(editrules)
  library(stats)
  library(nnet)
  library(foreign)

  path <- "//CBSP.NL/Productie/Secundair/MPOnderzoek/Werk/Combineren/Medewerkers/TWAL/Eigen projecten/Calibrated Imputation/Categorical data/R scripts/"
  pathdata <- "//CBSP.NL/Productie/Secundair/MPOnderzoek/Werk/Combineren/Medewerkers/TWAL/Eigen projecten/Calibrated Imputation/Categorical data/Test data/Corrected data/"
  pathresults <- "//CBSP.NL/Productie/Secundair/MPOnderzoek/Werk/Combineren/Medewerkers/TWAL/Eigen projecten/Calibrated Imputation/Categorical data/Simulation results/"
  
  logfilename <- paste0(path,"logfile.txt")
  cat(paste0(" Starting missing set ", s), file=logfilename, append = TRUE)
  
  #############################################################################################################################
  # READING TRUE DATA
  #############################################################################################################################
  truedata <- read.csv2(paste(pathdata,"IPUMS klein corrected.csv", sep=""), sep=";", header=TRUE)
  truedata$Geslacht <- as.factor(truedata$Geslacht)
  truedata$Leeftijd <- as.factor(truedata$Leeftijd)
  truedata$HH_Pos <- as.factor(truedata$HH_Pos)
  truedata$HH_grootte <- as.factor(truedata$HH_grootte)
  truedata$Woonregio.vorig.jaar <- as.factor(truedata$Woonregio.vorig.jaar)
  truedata$Nationaliteit <- as.factor(truedata$Nationaliteit)
  truedata$Geboorteland <- as.factor(truedata$Geboorteland)
  truedata$Onderwijsniveau <- as.factor(truedata$Onderwijsniveau)
  truedata$Econ..status <- as.factor(truedata$Econ..status)
  truedata$Beroep <- as.factor(truedata$Beroep)
  truedata$SBI <- as.factor(truedata$SBI)
  truedata$Burg..Staat <- as.factor(truedata$Burg..Staat)
  # ADJUST BASE LEVEL OF VARIABLES WITH ONLY TWO CATEGORIES
  truedata$Geslacht <- relevel(truedata$Geslacht, ref="1")
  
  #######################################################################################
  # STATEMENTS AND FUNCTIONS FOR INITIALISING AND READING DATA (INCLUDING TOTALEN) AND EDITS, AND CREATING MISSINGS
  #######################################################################################  
  
  nvars <- 12
  nrows <- 3784
  
  ncats <- rep(0,nvars)
  ncats[1] <- 2
  ncats[2] <- 17
  ncats[3] <- 8
  ncats[4] <- 6
  ncats[5] <- 3
  ncats[6] <- 3
  ncats[7] <- 3  
  ncats[8] <- 7
  ncats[9] <- 8
  ncats[10] <- 10
  ncats[11] <- 13  
  ncats[12] <- 4
  
  max_ncats <- 17
  
  codes <- matrix(NA, nvars, max_ncats)
  codes[1,] <- c('1','2', NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
  codes[2,] <- c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17')
  codes[3,] <- c('1110','1121','1122','1131','1132','1140','1210','1220', NA, NA, NA, NA, NA, NA, NA, NA, NA)
  codes[4,] <- c('111','112','113','114','125','126', NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
  codes[5,] <- c('1','2','9', NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
  codes[6,] <- c('1','2','3', NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
  codes[7,] <- c('1','2','3', NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
  codes[8,] <- c('0','1','2','3','4','5','9',  NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
  codes[9,] <- c('111','112','120','210','221','222','223','224',  NA, NA, NA, NA, NA, NA, NA, NA, NA)
  codes[10,] <- c('1','2','3','4','5','6','7','8','9','999', NA, NA, NA, NA, NA, NA, NA)
  codes[11,] <- c('111','122','124','131','132','133','134','135','136','137','138','139','200', NA, NA, NA, NA)
  codes[12,] <- c('1','2','3','4', NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
  
  max_num_missings <- 2000
  
  ########################
  # MARGINAL FREQUENCIES PER CATEGORY / WHEN NO MARGINAL FREQUENCIES ARE ASSUMED TO BE KNOWN FOR VARIABLE var, totals[var, cat] WILL BE SET TO "NA" FOR ALL CATEGORIES cat
  ########################
  totals <- matrix(NA, nvars, max_ncats)
  
  missings <- matrix(FALSE, nrows, nvars)
  num_missings <- rep(0, nvars)
  nummissingsrec <- rep(0,nrows)
  missing_indices <- matrix(0, nvars, max_num_missings)
  nummissingsrecinedits <- rep(0,nrows)
  
  imptotals <- matrix(NA,nvars,max_ncats)
  fractions <- matrix(0,nvars,max_ncats)
  
  Create_Missings <- function(intmissfactor) {
    intoriginaldata <- truedata
    MissingFractions <- rep(0,nvars)
    MissingFractions[1] <- 0.005*intmissfactor
    MissingFractions[2] <- 0.05*intmissfactor
    MissingFractions[3] <- 0.05*intmissfactor
    MissingFractions[4] <- 0.05*intmissfactor
    MissingFractions[5] <- 0.005*intmissfactor
    MissingFractions[6] <- 0.005*intmissfactor
    MissingFractions[7] <- 0.005*intmissfactor
    MissingFractions[8] <- 0.25*intmissfactor
    MissingFractions[9] <- 0.05*intmissfactor
    MissingFractions[10] <- 0.25*intmissfactor
    MissingFractions[11] <- 0.05*intmissfactor
    MissingFractions[12] <- 0.05*intmissfactor
    for (j in (1:nvars)) {
      nummissings_j <- round(nrows*MissingFractions[j])
      testsample <- sample((1:nrows), nummissings_j, replace = FALSE)
      for (i in (1:nummissings_j)) {
        intoriginaldata[testsample[i],j + 1] <- NA 
      }
    }
    return(intoriginaldata)
  }
  
  Create_MAR_Selective_Missings <- function(intmissfactor, intvarname) {
    intnummissing <- rep(0,12)
    if (intvarname == "edu") {
      intyoungdata <- truedata[truedata$Leeftijd %in% c("1","2","3","4","5","6"),]
      intnyoungrows <- dim(intyoungdata)[1]
      if (intmissfactor == 1) {
        intnummissing[1] <- 6
        intnummissing[2] <- 64
        intnummissing[3] <- 64
        intnummissing[4] <- 64
        intnummissing[5] <- 6
        intnummissing[6] <- 6
        intnummissing[7] <- 6
        intnummissing[8] <- 161
        intnummissing[9] <- 64
        intnummissing[10] <- 321
        intnummissing[11] <- 64
        intnummissing[12] <- 64
        for (j in (1:nvars)) {
          testsample <- sample((1:intnyoungrows), intnummissing[j], replace = FALSE)
          for (i in (1:intnummissing[j])) {
            intyoungdata[testsample[i],j + 1] <- NA 
          }
        } 
      }
      if (intmissfactor == 2) {
        intnummissing[1] <- 13
        intnummissing[2] <- 128
        intnummissing[3] <- 128
        intnummissing[4] <- 128
        intnummissing[5] <- 13
        intnummissing[6] <- 13
        intnummissing[7] <- 13
        intnummissing[8] <- 321
        intnummissing[9] <- 128
        intnummissing[10] <- 642
        intnummissing[11] <- 128
        intnummissing[12] <- 128
        for (j in (1:nvars)) {
          testsample <- sample((1:intnyoungrows), intnummissing[j], replace = FALSE)
          for (i in (1:intnummissing[j])) {
            intyoungdata[testsample[i],j + 1] <- NA 
          }
        } 
      }
      intmiddledata <- truedata[truedata$Leeftijd %in% c("7","8","9","10","11"),]
      intnmiddlerows <- dim(intmiddledata)[1]
      if (intmissfactor == 1) {
        intnummissing[1] <- 9
        intnummissing[2] <- 86
        intnummissing[3] <- 86
        intnummissing[4] <- 86
        intnummissing[5] <- 9
        intnummissing[6] <- 9
        intnummissing[7] <- 9
        intnummissing[8] <- 471
        intnummissing[9] <- 86
        intnummissing[10] <- 429
        intnummissing[11] <- 86
        intnummissing[12] <- 86
        for (j in (1:nvars)) {
          testsample <- sample((1:intnmiddlerows), intnummissing[j], replace = FALSE)
          for (i in (1:intnummissing[j])) {
            intmiddledata[testsample[i],j + 1] <- NA 
          }
        } 
      }
      if (intmissfactor == 2) {
        intnummissing[1] <- 17
        intnummissing[2] <- 171
        intnummissing[3] <- 171
        intnummissing[4] <- 171
        intnummissing[5] <- 17
        intnummissing[6] <- 17
        intnummissing[7] <- 17
        intnummissing[8] <- 943
        intnummissing[9] <- 171
        intnummissing[10] <- 857
        intnummissing[11] <- 171
        intnummissing[12] <- 171
        for (j in (1:nvars)) {
          testsample <- sample((1:intnmiddlerows), intnummissing[j], replace = FALSE)
          for (i in (1:intnummissing[j])) {
            intmiddledata[testsample[i],j + 1] <- NA 
          }
        } 
      }
      intolddata <- truedata[truedata$Leeftijd %in% c("12","13","14","15","16","17"),]
      intnoldrows <- dim(intolddata)[1]
      if (intmissfactor == 1) {
        intnummissing[1] <- 4
        intnummissing[2] <- 39
        intnummissing[3] <- 39
        intnummissing[4] <- 39
        intnummissing[5] <- 4
        intnummissing[6] <- 4
        intnummissing[7] <- 4
        intnummissing[8] <- 314
        intnummissing[9] <- 39
        intnummissing[10] <- 196
        intnummissing[11] <- 39
        intnummissing[12] <- 39
        for (j in (1:nvars)) {
          testsample <- sample((1:intnoldrows), intnummissing[j], replace = FALSE)
          for (i in (1:intnummissing[j])) {
            intolddata[testsample[i],j + 1] <- NA 
          }
        } 
      }
      if (intmissfactor == 2) {
        intnummissing[1] <- 8
        intnummissing[2] <- 79
        intnummissing[3] <- 79
        intnummissing[4] <- 79
        intnummissing[5] <- 8
        intnummissing[6] <- 8
        intnummissing[7] <- 8
        intnummissing[8] <- 628
        intnummissing[9] <- 79
        intnummissing[10] <- 393
        intnummissing[11] <- 79
        intnummissing[12] <- 79
        for (j in (1:nvars)) {
          testsample <- sample((1:intnoldrows), intnummissing[j], replace = FALSE)
          for (i in (1:intnummissing[j])) {
            intolddata[testsample[i],j + 1] <- NA 
          }
        } 
      }
    }
    if (intvarname == "occ") {
      intyoungdata <- truedata[truedata$Leeftijd %in% c("1","2","3","4","5","6"),]
      intnyoungrows <- dim(intyoungdata)[1]
      if (intmissfactor == 1) {
        intnummissing[1] <- 6
        intnummissing[2] <- 64
        intnummissing[3] <- 64
        intnummissing[4] <- 64
        intnummissing[5] <- 6
        intnummissing[6] <- 6
        intnummissing[7] <- 6
        intnummissing[8] <- 321
        intnummissing[9] <- 64
        intnummissing[10] <- 161
        intnummissing[11] <- 64
        intnummissing[12] <- 64
        for (j in (1:nvars)) {
          testsample <- sample((1:intnyoungrows), intnummissing[j], replace = FALSE)
          for (i in (1:intnummissing[j])) {
            intyoungdata[testsample[i],j + 1] <- NA 
          }
        } 
      }
      if (intmissfactor == 2) {
        intnummissing[1] <- 13
        intnummissing[2] <- 128
        intnummissing[3] <- 128
        intnummissing[4] <- 128
        intnummissing[5] <- 13
        intnummissing[6] <- 13
        intnummissing[7] <- 13
        intnummissing[8] <- 642
        intnummissing[9] <- 128
        intnummissing[10] <- 321
        intnummissing[11] <- 128
        intnummissing[12] <- 128
        for (j in (1:nvars)) {
          testsample <- sample((1:intnyoungrows), intnummissing[j], replace = FALSE)
          for (i in (1:intnummissing[j])) {
            intyoungdata[testsample[i],j + 1] <- NA 
          }
        } 
      }
      intmiddledata <- truedata[truedata$Leeftijd %in% c("7","8","9","10","11"),]
      intnmiddlerows <- dim(intmiddledata)[1]
      if (intmissfactor == 1) {
        intnummissing[1] <- 9
        intnummissing[2] <- 86
        intnummissing[3] <- 86
        intnummissing[4] <- 86
        intnummissing[5] <- 9
        intnummissing[6] <- 9
        intnummissing[7] <- 9
        intnummissing[8] <- 429
        intnummissing[9] <- 86
        intnummissing[10] <- 471
        intnummissing[11] <- 86
        intnummissing[12] <- 86
        for (j in (1:nvars)) {
          testsample <- sample((1:intnmiddlerows), intnummissing[j], replace = FALSE)
          for (i in (1:intnummissing[j])) {
            intmiddledata[testsample[i],j + 1] <- NA 
          }
        } 
      }
      if (intmissfactor == 2) {
        intnummissing[1] <- 17
        intnummissing[2] <- 171
        intnummissing[3] <- 171
        intnummissing[4] <- 171
        intnummissing[5] <- 17
        intnummissing[6] <- 17
        intnummissing[7] <- 17
        intnummissing[8] <- 857
        intnummissing[9] <- 171
        intnummissing[10] <- 943
        intnummissing[11] <- 171
        intnummissing[12] <- 171
        for (j in (1:nvars)) {
          testsample <- sample((1:intnmiddlerows), intnummissing[j], replace = FALSE)
          for (i in (1:intnummissing[j])) {
            intmiddledata[testsample[i],j + 1] <- NA 
          }
        } 
      }
      intolddata <- truedata[truedata$Leeftijd %in% c("12","13","14","15","16","17"),]
      intnoldrows <- dim(intolddata)[1]
      if (intmissfactor == 1) {
        intnummissing[1] <- 4
        intnummissing[2] <- 39
        intnummissing[3] <- 39
        intnummissing[4] <- 39
        intnummissing[5] <- 4
        intnummissing[6] <- 4
        intnummissing[7] <- 4
        intnummissing[8] <- 196
        intnummissing[9] <- 39
        intnummissing[10] <- 314
        intnummissing[11] <- 39
        intnummissing[12] <- 39
        for (j in (1:nvars)) {
          testsample <- sample((1:intnoldrows), intnummissing[j], replace = FALSE)
          for (i in (1:intnummissing[j])) {
            intolddata[testsample[i],j + 1] <- NA 
          }
        } 
      }
      if (intmissfactor == 2) {
        intnummissing[1] <- 8
        intnummissing[2] <- 79
        intnummissing[3] <- 79
        intnummissing[4] <- 79
        intnummissing[5] <- 8
        intnummissing[6] <- 8
        intnummissing[7] <- 8
        intnummissing[8] <- 393
        intnummissing[9] <- 79
        intnummissing[10] <- 628
        intnummissing[11] <- 79
        intnummissing[12] <- 79
        for (j in (1:nvars)) {
          testsample <- sample((1:intnoldrows), intnummissing[j], replace = FALSE)
          for (i in (1:intnummissing[j])) {
            intolddata[testsample[i],j + 1] <- NA 
          }
        } 
      }
    }
    intoriginaldata <- rbind(intyoungdata, intmiddledata, intolddata)
    return(intoriginaldata)
  }
  
  Create_NMAR_Selective_Missings <- function(intmissfactor, intvarname) {
    intnummissing <- rep(0,12)
    if (intvarname == "edu") {
      intedu1data <- truedata[truedata$Onderwijsniveau %in% c("0","1","9"),]
      intnedu1rows <- dim(intedu1data)[1]
      if (intmissfactor == 1) {
        intnummissing[1] <- 6
        intnummissing[2] <- 55
        intnummissing[3] <- 55
        intnummissing[4] <- 55
        intnummissing[5] <- 6
        intnummissing[6] <- 6
        intnummissing[7] <- 6
        intnummissing[8] <- 109
        intnummissing[9] <- 55
        intnummissing[10] <- 272
        intnummissing[11] <- 55
        intnummissing[12] <- 55
        for (j in (1:nvars)) {
          testsample <- sample((1:intnedu1rows), intnummissing[j], replace = FALSE)
          for (i in (1:intnummissing[j])) {
            intedu1data[testsample[i],j + 1] <- NA 
          }
        } 
      }
      if (intmissfactor == 2) {
        intnummissing[1] <- 12
        intnummissing[2] <- 110
        intnummissing[3] <- 110
        intnummissing[4] <- 110
        intnummissing[5] <- 12
        intnummissing[6] <- 12
        intnummissing[7] <- 12
        intnummissing[8] <- 218
        intnummissing[9] <- 110
        intnummissing[10] <- 544
        intnummissing[11] <- 110
        intnummissing[12] <- 110
        for (j in (1:nvars)) {
          testsample <- sample((1:intnedu1rows), intnummissing[j], replace = FALSE)
          for (i in (1:intnummissing[j])) {
            intedu1data[testsample[i],j + 1] <- NA 
          }
        } 
      }
      intedu2data <- truedata[truedata$Onderwijsniveau %in% c("4","5"),]
      intnedu2rows <- dim(intedu2data)[1]      
      if (intmissfactor == 1) {
        intnummissing[1] <- 3
        intnummissing[2] <- 33
        intnummissing[3] <- 33
        intnummissing[4] <- 33
        intnummissing[5] <- 33
        intnummissing[6] <- 3
        intnummissing[7] <- 3
        intnummissing[8] <- 26
        intnummissing[9] <- 33
        intnummissing[10] <- 167
        intnummissing[11] <- 33
        intnummissing[12] <- 33
        for (j in (1:nvars)) {
          testsample <- sample((1:intnedu2rows), intnummissing[j], replace = FALSE)
          for (i in (1:intnummissing[j])) {
            intedu2data[testsample[i],j + 1] <- NA 
          }
        } 
      }
      if (intmissfactor == 2) {
        intnummissing[1] <- 6
        intnummissing[2] <- 66
        intnummissing[3] <- 66
        intnummissing[4] <- 66
        intnummissing[5] <- 6
        intnummissing[6] <- 6
        intnummissing[7] <- 6
        intnummissing[8] <- 52
        intnummissing[9] <- 66
        intnummissing[10] <- 334
        intnummissing[11] <- 66
        intnummissing[12] <- 66
        for (j in (1:nvars)) {
          testsample <- sample((1:intnedu2rows), intnummissing[j], replace = FALSE)
          for (i in (1:intnummissing[j])) {
            intedu2data[testsample[i],j + 1] <- NA 
          }
        } 
      }
      intedu3data <- truedata[truedata$Onderwijsniveau %in% c("2","3"),]
      intnedu3rows <- dim(intedu3data)[1]
      if (intmissfactor == 1) {
        intnummissing[1] <- 10
        intnummissing[2] <- 101
        intnummissing[3] <- 101
        intnummissing[4] <- 101
        intnummissing[5] <- 10
        intnummissing[6] <- 10
        intnummissing[7] <- 10
        intnummissing[8] <- 811
        intnummissing[9] <- 101
        intnummissing[10] <- 507
        intnummissing[11] <- 101
        intnummissing[12] <- 101
        for (j in (1:nvars)) {
          testsample <- sample((1:intnedu3rows), intnummissing[j], replace = FALSE)
          for (i in (1:intnummissing[j])) {
            intedu3data[testsample[i],j + 1] <- NA 
          }
        } 
      }
      if (intmissfactor == 2) {
        intnummissing[1] <- 20
        intnummissing[2] <- 202
        intnummissing[3] <- 202
        intnummissing[4] <- 202
        intnummissing[5] <- 20
        intnummissing[6] <- 20
        intnummissing[7] <- 20
        intnummissing[8] <- 1622
        intnummissing[9] <- 202
        intnummissing[10] <- 1014
        intnummissing[11] <- 202
        intnummissing[12] <- 202
        for (j in (1:nvars)) {
          testsample <- sample((1:intnedu3rows), intnummissing[j], replace = FALSE)
          for (i in (1:intnummissing[j])) {
            intedu3data[testsample[i],j + 1] <- NA 
          }
        } 
      }
    }
    if (intvarname == "occ") {
      intocc1data <- truedata[truedata$Beroep %in% c("1","2","3","4","7"),]
      intnocc1rows <- dim(intocc1data)[1]
      if (intmissfactor == 1) {
        intnummissing[1] <- 7
        intnummissing[2] <- 71
        intnummissing[3] <- 71
        intnummissing[4] <- 71
        intnummissing[5] <- 7
        intnummissing[6] <- 7
        intnummissing[7] <- 7
        intnummissing[8] <- 355
        intnummissing[9] <- 71
        intnummissing[10] <- 355
        intnummissing[11] <- 71
        intnummissing[12] <- 71
        for (j in (1:nvars)) {
          testsample <- sample((1:intnocc1rows), intnummissing[j], replace = FALSE)
          for (i in (1:intnummissing[j])) {
            intocc1data[testsample[i],j + 1] <- NA 
          }
        } 
      }
      if (intmissfactor == 2) {
        intnummissing[1] <- 14
        intnummissing[2] <- 142
        intnummissing[3] <- 142
        intnummissing[4] <- 142
        intnummissing[5] <- 14
        intnummissing[6] <- 14
        intnummissing[7] <- 14
        intnummissing[8] <- 710
        intnummissing[9] <- 142
        intnummissing[10] <- 710
        intnummissing[11] <- 142
        intnummissing[12] <- 142
        for (j in (1:nvars)) {
          testsample <- sample((1:intnocc1rows), intnummissing[j], replace = FALSE)
          for (i in (1:intnummissing[j])) {
            intocc1data[testsample[i],j + 1] <- NA 
          }
        } 
      }
      intocc2data <- truedata[truedata$Beroep %in% c("5","6","8","9"),]
      intnocc2rows <- dim(intocc2data)[1]
      if (intmissfactor == 1) {
        intnummissing[1] <- 3
        intnummissing[2] <- 29
        intnummissing[3] <- 29
        intnummissing[4] <- 29
        intnummissing[5] <- 3
        intnummissing[6] <- 3
        intnummissing[7] <- 3
        intnummissing[8] <- 147
        intnummissing[9] <- 29
        intnummissing[10] <- 236
        intnummissing[11] <- 29
        intnummissing[12] <- 29
        for (j in (1:nvars)) {
          testsample <- sample((1:intnocc2rows), intnummissing[j], replace = FALSE)
          for (i in (1:intnummissing[j])) {
            intocc2data[testsample[i],j + 1] <- NA 
          }
        } 
      }
      if (intmissfactor == 2) {
        intnummissing[1] <- 6
        intnummissing[2] <- 58
        intnummissing[3] <- 58
        intnummissing[4] <- 58
        intnummissing[5] <- 6
        intnummissing[6] <- 6
        intnummissing[7] <- 6
        intnummissing[8] <- 294
        intnummissing[9] <- 58
        intnummissing[10] <- 472
        intnummissing[11] <- 58
        intnummissing[12] <- 58
        for (j in (1:nvars)) {
          testsample <- sample((1:intnocc2rows), intnummissing[j], replace = FALSE)
          for (i in (1:intnummissing[j])) {
            intocc2data[testsample[i],j + 1] <- NA 
          }
        } 
      }
      intocc3data <- truedata[truedata$Beroep %in% c("999"),]
      intnocc3rows <- dim(intocc3data)[1]
      if (intmissfactor == 1) {
        intnummissing[1] <- 9
        intnummissing[2] <- 89
        intnummissing[3] <- 89
        intnummissing[4] <- 89
        intnummissing[5] <- 9
        intnummissing[6] <- 9
        intnummissing[7] <- 9
        intnummissing[8] <- 444
        intnummissing[9] <- 89
        intnummissing[10] <- 355
        intnummissing[11] <- 89
        intnummissing[12] <- 89
        for (j in (1:nvars)) {
          testsample <- sample((1:intnocc3rows), intnummissing[j], replace = FALSE)
          for (i in (1:intnummissing[j])) {
            intocc3data[testsample[i],j + 1] <- NA 
          }
        } 
      }
      if (intmissfactor == 2) {
        intnummissing[1] <- 18
        intnummissing[2] <- 178
        intnummissing[3] <- 178
        intnummissing[4] <- 178
        intnummissing[5] <- 18
        intnummissing[6] <- 18
        intnummissing[7] <- 18
        intnummissing[8] <- 888
        intnummissing[9] <- 178
        intnummissing[10] <- 710
        intnummissing[11] <- 178
        intnummissing[12] <- 178
        for (j in (1:nvars)) {
          testsample <- sample((1:intnocc3rows), intnummissing[j], replace = FALSE)
          for (i in (1:intnummissing[j])) {
            intocc3data[testsample[i],j + 1] <- NA 
          }
        } 
      }
    }
    if (intvarname == "edu") {
      intoriginaldata <- rbind(intedu1data, intedu2data, intedu3data)
    }
    if (intvarname == "occ") {
      intoriginaldata <- rbind(intocc1data, intocc2data, intocc3data)
    }
    return(intoriginaldata)
  }
  
  ########################
  # SPECIFY EXPLICIT EDITS
  ########################
  ExpEdits <- editarray(expression(
    Geslacht %in% c('1','2'),
    Leeftijd %in% c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17'),
    HH_Pos %in% c('1110','1121','1122','1131','1132','1140','1210','1220'),
    HH_grootte %in% c('111','112','113','114','125','126'),
    Woonregio.vorig.jaar %in% c('1','2','9'),
    Nationaliteit %in% c('1','2','3'),
    Geboorteland %in% c('1','2','3'),
    Onderwijsniveau %in% c('0','1','2','3','4','5','9'),
    Econ..status %in% c('111','112','120','210','221','222','223','224'),
    Beroep %in% c('1','2','3','4','5','6','7','8','9','999'),
    SBI %in% c('111','122','124','131','132','133','134','135','136','137','138','139','200'),
    Burg..Staat %in% c('1','2','3','4'),
    if (Leeftijd %in% c('1','2','3')) (Burg..Staat == '1'),
    if (HH_Pos %in% c('1121','1122')) (Burg..Staat == '2'),
    if (Leeftijd %in% c('1','2')) (Onderwijsniveau %in% c('0','9')),
    if (Leeftijd %in% c('1','2','3')) (Onderwijsniveau %in% c('0','1','9')),
    if (Leeftijd %in% c('1','2','3')) (HH_Pos %in% c('1110','1220')),
    if (Leeftijd %in% c('1','2','3','4')) (Onderwijsniveau %in% c('0','1','2','3','4','9')),
    if (Leeftijd %in% c('1','2','3')) (Econ..status %in% c('112','221','224')),
    if (Leeftijd %in% c('1','2','3')) (Beroep %in% c('999')),
    if (Leeftijd %in% c('1','2','3')) (SBI %in% c('200')),
    if (HH_grootte == '111') (HH_Pos == '1210'),
    if (HH_grootte == '112') (HH_Pos %in% c('1110','1121','1131','1140','1220')),
    if (HH_grootte %in% c('113','114','125','126')) (HH_Pos != '1210')
  ))
  
  Determine_totals <- function(intvarnr) {
    inttotalsintvarnr <- rep(NA, max_ncats)
    mytrueTable <- table(truedata[,intvarnr + 1])
    for (k in (1:ncats[intvarnr])) {
      inttotalsintvarnr[k] <- margin.table(mytrueTable,1)[k]
    }
    return(inttotalsintvarnr)
  }
  
  ########################
  # MARGINAL TOTALS THAT ARE ASSUMED TO BE KNOWN / IN OUR SIMULATION STUDY USUALLY ONLY TOTALS FOR VARIABLE 8 ("EDUCATIONAL LEVEL") ARE ASSUMED TO BE KNOWN
  ########################
  # totals[1,] <- Determine_totals(1)
  # totals[2,] <- Determine_totals(2)
  # totals[3,] <- Determine_totals(3)
  # totals[4,] <- Determine_totals(4)
  # totals[5,] <- Determine_totals(5)
  # totals[6,] <- Determine_totals(6)
  # totals[7,] <- Determine_totals(7)
  totals[8,] <- Determine_totals(8)
  # totals[9,] <- Determine_totals(9)
  # totals[10,] <- Determine_totals(10)
  # totals[11,] <- Determine_totals(11)
  # totals[12,] <- Determine_totals(12)


  #############################################################################################################################
  # FUNCTIONS FOR START-UP PHASE
  #############################################################################################################################
  
  ########################
  # SIMPLE IMPUTATION METHOD / IMPUTATION PROBABILITIES BASED ON OBSERVED PROPORTIONS
  ########################
  SimpleImpute <- function(intvarnr, intfractions) {
    draw <- rmultinom(1, 1, intfractions[])
    for (kk in (1:length(draw))) {
      if (draw[kk] == 1) {
        drawvalue <- kk
      }
    }
    return(drawvalue)
  }
  
  ########################
  # ENSURING THAT IMPUTATIONS SATISFY EDITS
  ########################
  ImputeUntilSatisfied <- function(intvarnr, intEditsInRec) {
    intsatisfied <- FALSE
    intfractions <- rep(0, ncats[intvarnr])
    intsumfractions <- 0
    for (jj in (1:ncats[intvarnr])) {
      intfractions[jj] <- fractions[intvarnr, jj]
      intsumfractions <- intsumfractions + intfractions[jj]
    }
    adjustedfracs <- 0
    while (!intsatisfied && (adjustedfracs < 2)) {
      if (intsumfractions < 1e-8) {
        for (jj in (1:ncats[intvarnr])) {
          intfractions[jj] <- 1/ncats[intvarnr]
        }
        adjustedfracs <- adjustedfracs + 1
      }
      intimpvalue <- SimpleImpute(intvarnr, intfractions)
      codedintimpvalue <- codes[intvarnr, intimpvalue]
      if (intvarnr == 1) {
        testEdits <- substValue(intEditsInRec,"Geslacht", codedintimpvalue, reduce = FALSE)  
      }
      if (intvarnr == 2) {
        testEdits <- substValue(intEditsInRec,"Leeftijd", codedintimpvalue, reduce = FALSE)  
      }
      if (intvarnr == 3) {
        testEdits <- substValue(intEditsInRec,"HH_Pos", codedintimpvalue, reduce = FALSE)  
      }
      if (intvarnr == 4) {
        testEdits <- substValue(intEditsInRec,"HH_grootte", codedintimpvalue, reduce = FALSE) 
      }
      if (intvarnr == 5) {
        testEdits <- substValue(intEditsInRec,"Woonregio.vorig.jaar", codedintimpvalue, reduce = FALSE)  
      }
      if (intvarnr == 6) {
        testEdits <- substValue(intEditsInRec,"Nationaliteit", codedintimpvalue, reduce = FALSE)  
      }
      if (intvarnr == 7) {
        testEdits <- substValue(intEditsInRec,"Geboorteland", codedintimpvalue, reduce = FALSE)  
      }
      if (intvarnr == 8) {
        testEdits <- substValue(intEditsInRec,"Onderwijsniveau", codedintimpvalue, reduce = FALSE)  
      }
      if (intvarnr == 9) {
        testEdits <- substValue(intEditsInRec,"Econ..status", codedintimpvalue, reduce = FALSE)  
      }
      if (intvarnr == 10) {
        testEdits <- substValue(intEditsInRec,"Beroep", codedintimpvalue, reduce = FALSE)  
      }
      if (intvarnr == 11) {
        testEdits <- substValue(intEditsInRec,"SBI", codedintimpvalue, reduce = FALSE)  
      }
      if (intvarnr == 12) {
        testEdits <- substValue(intEditsInRec,"Burg..Staat", codedintimpvalue, reduce = FALSE)  
      }
      intsatisfied <- isFeasible(testEdits)
      if (!intsatisfied) {
        intsumfractions <- 1 - intfractions[intimpvalue]
        if (intsumfractions < 1e-12) {
          for (k in (1:ncats[intvarnr])) {
            intfractions[k] <- 1/ncats[intvarnr]
          }
          adjustedfracs <- adjustedfracs + 1
        }
        else {
          intfractions[intimpvalue] <- 0
          for (k in (1:ncats[intvarnr])) {
            if (k != intimpvalue) {
              intfractions[k] <- intfractions[k]/intsumfractions  
            }
          }
        }
      }
    }
    rm(intfractions)
    gc()
    return(codedintimpvalue)
  }
  
  ReduceEdits <- function(intrownr, intvarnr, intEditsInRec, intnummissingsrecinedits) {
    intEditsInRecReduced <- intEditsInRec
    if ((intvarnr %in% (2:4)) | (intvarnr %in% (8:12))) {
      if (missings[intrownr, intvarnr]) {
        if (intnummissingsrecinedits > 1) {
          if (intvarnr == 2) {
            intEditsInRecReduced <- eliminate(intEditsInRec, "Leeftijd", reduce = TRUE)
          }
          if (intvarnr == 3) {
            intEditsInRecReduced <- eliminate(intEditsInRec, "HH_Pos", reduce = TRUE)
          }        
          if (intvarnr == 4) {
            intEditsInRecReduced <- eliminate(intEditsInRec, "HH_grootte", reduce = TRUE)
          }        
          if (intvarnr == 8) {
            intEditsInRecReduced <- eliminate(intEditsInRec, "Onderwijsniveau", reduce = TRUE)
          }        
          if (intvarnr == 9) {
            intEditsInRecReduced <- eliminate(intEditsInRec, "Econ..status", reduce = TRUE)
          }        
          if (intvarnr == 10) {
            intEditsInRecReduced <- eliminate(intEditsInRec, "Beroep", reduce = TRUE)
          }        
          if (intvarnr == 11) {
            intEditsInRecReduced <- eliminate(intEditsInRec, "SBI", reduce = TRUE)
          }        
          if (intvarnr == 12) {
            intEditsInRecReduced <- eliminate(intEditsInRec, "Burg..Staat", reduce = TRUE)
          }
          intnummissingsrecinedits <- intnummissingsrecinedits - 1  
        }
      }
      else {
        if (!missings[intrownr, intvarnr]) {
          if (intvarnr == 2) {
            intEditsInRecReduced <- substValue(intEditsInRec,"Leeftijd",originaldata$Leeftijd[intrownr], reduce = TRUE)  
          }
          if (intvarnr == 3) {
            intEditsInRecReduced <- substValue(intEditsInRec, "HH_Pos",originaldata$HH_Pos[intrownr], reduce = TRUE)
          }        
          if (intvarnr == 4) {
            intEditsInRecReduced <- substValue(intEditsInRec, "HH_grootte", originaldata$HH_grootte[intrownr], reduce = TRUE)
          }        
          if (intvarnr == 8) {
            intEditsInRecReduced <- substValue(intEditsInRec, "Onderwijsniveau",originaldata$Onderwijsniveau[intrownr], reduce = TRUE)
          }        
          if (intvarnr == 9) {
            intEditsInRecReduced <- substValue(intEditsInRec, "Econ..status",originaldata$Econ..status[intrownr], reduce = TRUE)
          }        
          if (intvarnr == 10) {
            intEditsInRecReduced <- substValue(intEditsInRec, "Beroep",originaldata$Beroep[intrownr], reduce = TRUE)
          }        
          if (intvarnr == 11) {
            intEditsInRecReduced <- substValue(intEditsInRec, "SBI",originaldata$SBI[intrownr], reduce = TRUE)
          }        
          if (intvarnr == 12) {
            intEditsInRecReduced <- substValue(intEditsInRec, "Burg..Staat",originaldata$Burg..Staat[intrownr], reduce = TRUE)
          }
        }
      }
    }
    else{
      if (!missings[intrownr, intvarnr]) {
        if (intvarnr == 1) {
          intEditsInRecReduced <- substValue(intEditsInRec,"Geslacht",originaldata$Geslacht[intrownr], reduce = TRUE)  
        }
        if (intvarnr == 5) {
          intEditsInRecReduced <- substValue(intEditsInRec, "Woonregio.vorig.jaar",originaldata$Woonregio.vorig.jaar[intrownr], reduce = TRUE)
        }        
        if (intvarnr == 6) {
          intEditsInRecReduced <- substValue(intEditsInRec, "Nationaliteit", originaldata$Nationaliteit[intrownr], reduce = TRUE)
        }        
        if (intvarnr == 7) {
          intEditsInRecReduced <- substValue(intEditsInRec, "Geboorteland",originaldata$Geboorteland[intrownr], reduce = TRUE)
        }        
      }
    }
    return(intEditsInRecReduced)
  }
  
  #############################################################################################################################
  # FUNCTIONS FOR IMPUTATION PHASE
  #############################################################################################################################
  
  ########################
  # ESTIMATING MULTINOMIAL MODEL / DETERMINING IMPUTATION PROBABILITIES ACCORDING TO MULTINOMIAL MODEL
  ########################
  EstimateProbabilities <- function(intvarnr) {
    maxiterations <- 1000
    # VARIABLE Geslacht  
    # FOR VARIABLES WITH ONLY TWO CATEGORIES WE NEED TO MAKE VECTOR WITH TWO PROBABITIES (ONE FOR EACH CATEGORY)
    if (intvarnr == 1) {
      logresults <- multinom(Geslacht ~ Leeftijd + HH_Pos + HH_grootte + Woonregio.vorig.jaar + Nationaliteit + Geboorteland + Onderwijsniveau + Econ..status + Beroep + SBI + Burg..Staat, maxit = maxiterations, data = currentdata)
      pp2_Vrouw <- fitted(logresults)  # MODEL PRODUCES PROBABILITY  FOR "female" (CATEGORY 2 ; "Vrouw" IN DUTCH), SINCE BASIC LEVEL IS male (CATEGORY 1; "man" IN DUTCH)
      pp2_Man <- 1 - pp2_Vrouw
      pp <- cbind(pp2_Man, pp2_Vrouw)    
    }
    # VARIABLE Leeftijd
    if (intvarnr == 2) {
      logresults <- multinom(Leeftijd ~ Geslacht + HH_Pos + HH_grootte + Woonregio.vorig.jaar + Nationaliteit + Geboorteland + Onderwijsniveau + Econ..status + Beroep + SBI + Burg..Staat, maxit = maxiterations, data = currentdata)
      pp <- fitted(logresults)
    }
    # VARIABLE HH_Pos
    if (intvarnr == 3) {
      logresults_HH_Pos <- multinom(HH_Pos ~ Geslacht + Leeftijd + HH_grootte + Woonregio.vorig.jaar + Nationaliteit + Geboorteland + Onderwijsniveau + Econ..status + Beroep + SBI + Burg..Staat, maxit = maxiterations, data = currentdata)
      pp <- fitted(logresults_HH_Pos)
    }
    # VARIABLE HH_grootte
    if (intvarnr == 4) {
      logresults <- multinom(HH_grootte ~ Geslacht + Leeftijd + HH_Pos + Woonregio.vorig.jaar + Nationaliteit + Geboorteland + Onderwijsniveau + Econ..status + Beroep + SBI + Burg..Staat, maxit = maxiterations, data = currentdata)
      pp <- fitted(logresults)  
    }
    # VARIABLE Woonregio vorig jaar
    if (intvarnr == 5) {
      logresults <- multinom(Woonregio.vorig.jaar ~ Geslacht + Leeftijd + HH_Pos + HH_grootte + Nationaliteit + Geboorteland + Onderwijsniveau + Econ..status + Beroep + SBI + Burg..Staat, maxit = maxiterations, data = currentdata)
      pp <- fitted(logresults)
    }
    # VARIABLE Nationaliteit
    if (intvarnr == 6) {
      logresults <- multinom(Nationaliteit ~ Geslacht + Leeftijd + HH_Pos + HH_grootte + Woonregio.vorig.jaar + Geboorteland + Onderwijsniveau + Econ..status + Beroep + SBI + Burg..Staat, maxit = maxiterations, data = currentdata)
      pp <- fitted(logresults)
    }
    # VARIABLE Geboorteland
    if (intvarnr == 7) {
      logresults <- multinom(Geboorteland ~ Geslacht + Leeftijd + HH_Pos + HH_grootte + Woonregio.vorig.jaar + Nationaliteit + Onderwijsniveau + Econ..status + Beroep + SBI + Burg..Staat, maxit = maxiterations, data = currentdata)
      pp <- fitted(logresults)
    }
    # VARIABLE Onderwijsniveau
    if (intvarnr == 8) {
      logresults <- multinom(Onderwijsniveau ~ Geslacht + Leeftijd + HH_Pos + HH_grootte + Woonregio.vorig.jaar + Nationaliteit + Geboorteland + Econ..status + Beroep + SBI + Burg..Staat, maxit = maxiterations, data = currentdata)
      pp <- fitted(logresults)
    }
    # VARIABLE Econ..Status
    if (intvarnr == 9) {
      logresults <- multinom(Econ..status ~ Geslacht + Leeftijd + HH_Pos + HH_grootte + Woonregio.vorig.jaar + Nationaliteit + Geboorteland + Onderwijsniveau + Beroep + SBI + Burg..Staat, maxit = maxiterations, data = currentdata)
      pp <- fitted(logresults)
    }
    # VARIABLE Beroep
    if (intvarnr == 10) {
      logresults <- multinom(Beroep ~ Geslacht + Leeftijd + HH_Pos + HH_grootte + Woonregio.vorig.jaar + Nationaliteit + Geboorteland + Onderwijsniveau + Econ..status + SBI + Burg..Staat, maxit = maxiterations, data = currentdata)
      pp <- fitted(logresults)
    }
    # VARIABLE SBI
    if (intvarnr == 11) {
      logresults <- multinom(SBI ~ Geslacht + Leeftijd + HH_Pos + HH_grootte + Woonregio.vorig.jaar + Nationaliteit + Geboorteland + Onderwijsniveau + Econ..status + Beroep + Burg..Staat, maxit = maxiterations, data = currentdata)
      pp <- fitted(logresults)
    }
    # VARIABLE Burg..Staat
    if (intvarnr == 12) {
      logresults <- multinom(Burg..Staat ~ Geslacht + Leeftijd + HH_Pos + HH_grootte + Woonregio.vorig.jaar + Nationaliteit + Geboorteland + Onderwijsniveau + Econ..status + Beroep + SBI, maxit = maxiterations, data = currentdata)
      pp <- fitted(logresults)
    }
    return(pp)
  }
  
  ########################
  # ADJUSTING IMPUTATION PROBABILITIES FOR STRUCTURAL ZEROES
  ########################
    UseTestCatsToAdjustProbs <- function(intvarnr, intEditsInRec, intprobs) {
    for (k in (1:ncats[intvarnr])) {
      intsumfractions <- 1
      intsatisfied <- FALSE
      if (intvarnr == 1) {
        testEdits <- substValue(intEditsInRec,"Geslacht", codes[intvarnr, k], reduce = FALSE)  
      }
      if (intvarnr == 2) {
        testEdits <- substValue(intEditsInRec,"Leeftijd", codes[intvarnr, k], reduce = FALSE)  
      }
      if (intvarnr == 3) {
        testEdits <- substValue(intEditsInRec,"HH_Pos", codes[intvarnr, k], reduce = FALSE)  
      }
      if (intvarnr == 4) {
        testEdits <- substValue(intEditsInRec,"HH_grootte", codes[intvarnr, k], reduce = FALSE) 
      }
      if (intvarnr == 5) {
        testEdits <- substValue(intEditsInRec,"Woonregio.vorig.jaar", codes[intvarnr, k], reduce = FALSE)  
      }
      if (intvarnr == 6) {
        testEdits <- substValue(intEditsInRec,"Nationaliteit", codes[intvarnr, k], reduce = FALSE)  
      }
      if (intvarnr == 7) {
        testEdits <- substValue(intEditsInRec,"Geboorteland", codes[intvarnr, k], reduce = FALSE)  
      }
      if (intvarnr == 8) {
        testEdits <- substValue(intEditsInRec,"Onderwijsniveau", codes[intvarnr, k], reduce = FALSE)  
      }
      if (intvarnr == 9) {
        testEdits <- substValue(intEditsInRec,"Econ..status", codes[intvarnr, k], reduce = FALSE)  
      }
      if (intvarnr == 10) {
        testEdits <- substValue(intEditsInRec,"Beroep", codes[intvarnr, k], reduce = FALSE)  
      }
      if (intvarnr == 11) {
        testEdits <- substValue(intEditsInRec,"SBI", codes[intvarnr, k], reduce = FALSE)  
      }
      if (intvarnr == 12) {
        testEdits <- substValue(intEditsInRec,"Burg..Staat", codes[intvarnr, k], reduce = FALSE)  
      }
      intsatisfied <- isFeasible(testEdits)
      if (!intsatisfied) {
        intsumfractions <- intsumfractions - intprobs[k]
        intprobs[k] <- 0
        if (intsumfractions > 1e-10) {
          for (kk in (1:ncats[intvarnr])) {
            if (kk != k) {
              intprobs[kk] <- intprobs[kk]/intsumfractions  
            }
          }        
        }
        else {
          for (kk in (1:ncats[intvarnr])) {
            if (kk != k) {
              intprobs[kk] <- intprobs[kk]/(1 - ncats[intvarnr])  
            }
          }  
        }
      }
    }
    return(intprobs)
  }
  
  AdjustProbs_StructZeroes_AllRecs <- function(intpp, intvarnr) {
    intpp_small <- matrix(0,num_missings[intvarnr],ncats[intvarnr])
    # NO EDITS, SO NO STRUCTURAL ZEROES
    if ((intvarnr == 1) | (intvarnr %in% (5:7))) {
      for (ii in (1:num_missings[intvarnr])) {
        i <- missing_indices[intvarnr, ii]
        intpp_small[ii,] <- intpp[i,]
      }    
    }
    # EDITS, SO PERHAPS STRUCTURAL ZEROES
    if ((intvarnr %in% (2:4)) | (intvarnr %in% (8:12))) {
      for (ii in (1:num_missings[intvarnr])) {
        i <- missing_indices[intvarnr, ii]
        EditsInRec <- ExpEdits
        if (intvarnr != 2) {
          EditsInRec <- substValue(EditsInRec,"Leeftijd",currentdata$Leeftijd[i], reduce = TRUE)
        }  
        if (intvarnr != 3) {
          EditsInRec <- substValue(EditsInRec,"HH_Pos",currentdata$HH_Pos[i], reduce = TRUE)
        }  
        if (intvarnr != 4) {
          EditsInRec <- substValue(EditsInRec,"HH_grootte",currentdata$HH_grootte[i], reduce = TRUE)
        }
        if (intvarnr != 8) {
          EditsInRec <- substValue(EditsInRec,"Onderwijsniveau",currentdata$Onderwijsniveau[i], reduce = TRUE)
        }
        if (intvarnr != 9) {
          EditsInRec <- substValue(EditsInRec,"Econ..status",currentdata$Econ..status[i], reduce = TRUE)
        }
        if (intvarnr != 10) {
          EditsInRec <- substValue(EditsInRec,"Beroep",currentdata$Beroep[i], reduce = TRUE)
        }
        if (intvarnr != 11) {
          EditsInRec <- substValue(EditsInRec,"SBI",currentdata$SBI[i], reduce = TRUE)
        }
        if (intvarnr != 12) {
          EditsInRec <- substValue(EditsInRec,"Burg..Staat",currentdata$Burg..Staat[i], reduce = TRUE)
        }
        intpp_small[ii,] <- UseTestCatsToAdjustProbs(intvarnr, EditsInRec, intpp[missing_indices[intvarnr,ii],])
      }  
    }
    return(intpp_small)
  }
  
  AdjustProbs_StructZeroes_Rec <- function(intpp_rec, intvarnr) {
    intpp_small_rec <- rep(0,ncats[intvarnr])
    # NO EDITS, SO NO STRUCTURAL ZEROES
    if ((intvarnr == 1) | (intvarnr %in% (5:7))) {
      intpp_small_rec[ii,] <- intpp_rec[]
    }    
    # EDITS, SO PERHAPS STRUCTURAL ZEROES
    if ((intvarnr %in% (2:4)) | (intvarnr %in% (8:12))) {
      EditsInRec <- ExpEdits
      if (intvarnr != 2) {
        EditsInRec <- substValue(EditsInRec,"Leeftijd",currentdata$Leeftijd[i], reduce = TRUE)
      }  
      if (intvarnr != 3) {
        EditsInRec <- substValue(EditsInRec,"HH_Pos",currentdata$HH_Pos[i], reduce = TRUE)
      }  
      if (intvarnr != 4) {
        EditsInRec <- substValue(EditsInRec,"HH_grootte",currentdata$HH_grootte[i], reduce = TRUE)
      }
      if (intvarnr != 8) {
        EditsInRec <- substValue(EditsInRec,"Onderwijsniveau",currentdata$Onderwijsniveau[i], reduce = TRUE)
      }
      if (intvarnr != 9) {
        EditsInRec <- substValue(EditsInRec,"Econ..status",currentdata$Econ..status[i], reduce = TRUE)
      }
      if (intvarnr != 10) {
        EditsInRec <- substValue(EditsInRec,"Beroep",currentdata$Beroep[i], reduce = TRUE)
      }
      if (intvarnr != 11) {
        EditsInRec <- substValue(EditsInRec,"SBI",currentdata$SBI[i], reduce = TRUE)
      }
      if (intvarnr != 12) {
        EditsInRec <- substValue(EditsInRec,"Burg..Staat",currentdata$Burg..Staat[i], reduce = TRUE)
      }
      intpp_small_rec[] <- UseTestCatsToAdjustProbs(intvarnr, EditsInRec, intpp_rec[missing_indices[intvarnr,ii],])
    }
    return(intpp_small_rec)
  }
  
  ########################
  # FUNCTION: BINOMIAL APPROXIMATION FOR IMPUTATION PROBABILITIES
  ########################
  Bin_Approximation <- function(int_pp, int_nummissings, intvarnr, int_ncats) {
    int_adjust_pp <- int_pp
    if (int_nummissings > 1) {
      for (jj in (1:int_ncats)) {
        Sum_column <- 0
        for (ii in (1:int_nummissings)) {
          Sum_column <- Sum_column + int_pp[ii,jj]
        }
        for (ii in (1:int_nummissings)) {
          if (imptotals[intvarnr,jj] > 0.5) {
            BinomialProb <- (Sum_column - int_pp[ii,jj])/(int_nummissings - 1)
            if (BinomialProb < 1) {
              BinomialCorrection <- (int_nummissings - imptotals[intvarnr,jj])/imptotals[intvarnr,jj]
              BinomialCorrection <- BinomialCorrection*BinomialProb/(1 - BinomialProb)
              int_adjust_pp[ii,jj] <- int_pp[ii,jj]/(int_pp[ii,jj] + (1 - int_pp[ii,jj])*BinomialCorrection)
            }
            else {
              int_adjust_pp[ii,jj] <- 1
            }
          }
          else {
            int_adjust_pp[ii,jj] <- 0
          }
        }    
      }
    }
    return(int_adjust_pp)
  }
  
  ########################
  # FUNCTION: POISSON APPROXIMATION FOR IMPUTATION PROBABILITIES
  ########################  
  Pois_Approximation <- function(int_pp, int_nummissings, intvarnr, int_ncats) {
    int_adjust_pp <- int_pp
    for (jj in (1:int_ncats)) {
      Sum_column <- 0
      for (ii in (1:int_nummissings)) {
        Sum_column <- Sum_column + int_pp[ii,jj]
      }
      for (ii in (1:int_nummissings))  {
        if (imptotals[intvarnr, jj] > 0.5) {
          PoissonCorrection <- (Sum_column - int_pp[ii,jj])/imptotals[intvarnr, jj]
          int_adjust_pp[ii,jj] <- int_pp[ii,jj]/(int_pp[ii,jj] + (1 - int_pp[ii,jj])*PoissonCorrection)          
        }
        else {
          int_adjust_pp[ii,jj] <- 0
        }
      }    
    }
    return(int_adjust_pp)
  }
  
  ########################
  # FUNCTION IPD
  ########################
  IPF <- function(int_pp, int_nummissings, intvarnr, int_ncats) {
    MaxDiff <- 1
    while (MaxDiff > IPFMaxCriterion) {
      MaxDiff <- 0
      for (ii in (1:int_nummissings)) {
        Sum_row <- 0
        for (jj in (1:int_ncats)) {
          Sum_row <- Sum_row + int_pp[ii,jj]
        }
        if (abs(Sum_row - 1) > MaxDiff) {
          MaxDiff <- abs(Sum_row - 1)
        }
        if (Sum_row != 0) {
          AdjustmentFactor <- 1/Sum_row
          for (jj in (1:int_ncats)) {
            int_pp[ii,jj] <- int_pp[ii,jj]*AdjustmentFactor
          }
        }
      }
      for (jj in (1:int_ncats)) {
        Sum_column <- 0
        for (ii in (1:int_nummissings)) {
          Sum_column <- Sum_column + int_pp[ii,jj]
        }
        if (abs(Sum_column - imptotals[intvarnr, jj]) > MaxDiff) {
          MaxDiff <- abs(Sum_column - imptotals[intvarnr, jj])
        }
        if (Sum_column != 0) {
          AdjustmentFactor <- imptotals[intvarnr,jj]/Sum_column
          for (ii in (1:int_nummissings)) {
            int_pp[ii,jj] <- int_pp[ii,jj]*AdjustmentFactor
          }
        }
      }
    }
    return(int_pp)
  }
  
  ########################
  # CALCULATING ADJUSTED IMPUTATION PROBABILITIES USING ONE OF THE THREE APPROXIMATION METHODS
  #######################
  ApproximateProbs <- function(intpp_small, intvarnr, intmethod) {
    if (intmethod == 1) {
      intpp_small <- Bin_Approximation(intpp_small, num_missings[intvarnr], intvarnr, ncats[intvarnr])  
    }
    if (intmethod == 2) {
      intpp_small <- Pois_Approximation(intpp_small, num_missings[intvarnr], intvarnr, ncats[intvarnr])
    }
    intpp_small <- IPF(intpp_small, num_missings[intvarnr], intvarnr, ncats[intvarnr])  
    return(intpp_small)
  }
  
  ##########################################################################
  ##  Cox' ALGORITHM FOR CONTROLLED RANDOM ROUNDING (TO BASE 1) (developed by Sander Scholtus, 2018; update September 2019)
  ##########################################################################
  
  
  ##### R code for controlled rounding of a two-dimensional tabel using the method of Cox (1987) #####
  ## Sander Scholtus, 2018-2019
  
  library(Rcpp)
  
  ##########
  #### Functions
  
  
  sourceCpp(code = '
            #include "Rcpp.h"
            using namespace Rcpp;
            // [[Rcpp::plugins("cpp11")]]
            
            /* Function that finds the first two elements of a row-column path;
            always look for an element in the first possible row and column. */
            // [[Rcpp::export]]
            IntegerVector find_first(IntegerMatrix Frac, int nr, int nk) {
            IntegerVector res = {-1, -1};
            IntegerVector v(nk);
            int rij = 0;
            while ( rij < nr ) {
            v = Frac( rij, _ );
            if (max(v) == 1) {
            res[0] = rij;
            res[1] = which_max(v);
            return res;
            } else {
            rij++;
            }
            }
            return res;
            }
            
            // Function that constructs a row-column path.
            // [[Rcpp::export]]
            IntegerMatrix construct_cycle(IntegerMatrix Frac, int maxpos) {
            int nr = Frac.nrow(), nk = Frac.ncol();
            IntegerVector Lrij = rep(-1, maxpos);
            IntegerVector Lkol = rep(-1, maxpos);
            
            IntegerVector knoop = find_first(Frac, nr, nk);
            if (knoop[0] == -1) return -1;
            Lrij[0] = knoop[0];
            Lkol[0] = knoop[1];
            Frac(Lrij[0],Lkol[0]) = 0;
            
            // for the second element, stay in row and take first possible new column
            Lrij[1] = Lrij[0];
            Lkol[1] = which_max(Frac(Lrij[0], _ ));
            Frac(Lrij[1],Lkol[1]) = 0;
            
            // now add new elements to row-column path until cycle is formed
            int s = 1;
            int current;
            int finished = 0;
            LogicalVector samevec(maxpos);
            while ((finished == 0) && (s < maxpos)) {
            
            // within current column, choose first available row
            s++;
            Lkol[s] = Lkol[s-1];
            current = which_max(Frac( _ , Lkol[s-1] ));
            Frac(current, Lkol[s]) = 0;
            
            // check whether cycle has been obtained
            samevec = (Lrij == current);
            Lrij[s] = current;
            if (is_true(any(samevec))) {
            finished = 1;
            break;
            }
            // otherwise, continue while loop
            
            // within current row, choose first available column
            s++;
            Lrij[s] = Lrij[s-1];
            current = which_max(Frac(Lrij[s-1], _ ));
            Frac(Lrij[s], current) = 0;
            
            // check whether cycle has been obtained
            samevec = (Lkol == current);
            Lkol[s] = current;
            if (is_true(any(samevec))) finished = 2;
            // otherwise, continue while loop
            }
            
            if (finished == 0) {
            // stop("No cycle was obtained.");
            return (-10);
            }
            
            // reduce L to just the cycle
            IntegerVector telvec = seq_len(maxpos);
            IntegerVector sametelvec = telvec[samevec];
            int same = max(sametelvec) - 1;
            IntegerMatrix L (2, s - same + 1);
            if (finished == 1) {
            // re-order L so the cycle is along a row-column path
            IntegerVector idx = seq(same + 1, s);
            IntegerVector auxrij = Lrij[idx];
            auxrij.push_back(Lrij[same]);
            L( 0 , _ ) = auxrij;
            IntegerVector auxkol = Lkol[idx];
            auxkol.push_back(Lkol[same]);
            L( 1 , _ ) = auxkol;
            return L;
            }
            if (finished == 2) {
            IntegerVector idx = seq(same, s);
            IntegerVector aux = Lrij[idx];
            L( 0 , _ ) = aux;
            aux = Lkol[idx];
            L( 1 , _ ) = aux;
            return L;
            }
            
            }
            
            // Function that constructs a consistent rounding.
            // [[Rcpp::export]]
            NumericMatrix apply_cycle(NumericMatrix C, IntegerMatrix L, double round_base) {
            
            // determine dmin and dpls
            int s = L.ncol();
            LogicalVector tf = {1, 0};
            LogicalVector odd = rep(tf, s/2);
            NumericVector dminvec(s);
            NumericVector dplsvec(s);
            bool oddi;
            for (int i = 0; i < s; ++i) {
            oddi = odd[i];
            if (oddi) {
            dminvec[i] = C(L(0,i), L(1,i));
            dplsvec[i] = round_base - C(L(0,i), L(1,i));
            } else {
            dminvec[i] = round_base - C(L(0,i), L(1,i));
            dplsvec[i] = C(L(0,i), L(1,i));
            }
            }
            double dmin = min(dminvec);
            double dpls = min(dplsvec);
            
            // determine probability of selecting dmin rather than dpls
            double probmin = dpls/(dmin+dpls);
            
            RNGScope scope;
            double u = runif(1)[0];
            if (u < probmin) { // select dmin
            for (int j = 0; j < s; ++j) {
            oddi = odd[j];
            if (oddi) {
            C(L(0,j), L(1,j)) = C(L(0,j), L(1,j)) - dmin;
            } else {
            C(L(0,j), L(1,j)) = C(L(0,j), L(1,j)) + dmin;
            }
            }
            } else { // select dpls
            for (int j = 0; j < s; ++j) {
            oddi = odd[j];
            if (oddi) {
            C(L(0,j), L(1,j)) = C(L(0,j), L(1,j)) + dpls;
            } else {
            C(L(0,j), L(1,j)) = C(L(0,j), L(1,j)) - dpls;
            }
            }
            }
            
            return C;
            }
            
            // Function that checks which elements have been rounded.
            // [[Rcpp::export]]
            int check_rounding(NumericMatrix C, LogicalMatrix frac, IntegerMatrix L, int remaining, double round_base, double tol) {
            
            int s = L.ncol();
            double val;
            bool check;
            for (int i = 0; i < s; ++i) {
            val = C(L(0,i), L(1,i));
            check = ((val < tol) || (val > round_base - tol));
            if (check == 1) {
            frac(L(0,i), L(1,i)) = 0;
            remaining--;
            }
            }
            
            return remaining;
            }
            ')
  
  
  
  ## function that performs the method of Cox (1987)
  ## for unbiased controlled rounding (to any rounding base)
  
  roundCox_cpp <- function(A, round.base = 1L, return.margins = FALSE, tol = 1e-8,
                           max.cycle = 50L, max.iter = 100000L, seed = NULL, verbose = FALSE) {
    
    if (!is.null(seed)) set.seed(seed)
    
    A.rs <- rowSums(A)
    A.cs <- colSums(A)
    A.gt <- sum(A)
    r <- nrow(A)
    k <- ncol(A)
    
    C <- rbind(cbind(A %% round.base, round.base - (A.rs %% round.base)),
               c(round.base - (A.cs %% round.base), (A.gt %% round.base)))
    
    fractions <- ((C %% round.base > tol) & (C %% round.base < round.base - tol))
    remaining <- sum(fractions)
    max.cycle <- min(max.cycle, remaining)
    iter <- 1L
    
    while (remaining > 0 & iter <= max.iter) {
      
      if (((iter == 1) | (iter %% 10 == 0)) & verbose) {
        cat(sprintf('Iteration %d - Number of cells left to round off: %d\n', iter, remaining))
      }
      
      # construct a cycle L of fractions in C
      L <- construct_cycle(fractions, max.cycle)
      #print(L)
      
      C <- apply_cycle(C, L, round.base)
      
      # stop if there are no fractions left in C (or the maximum number of
      # iterations has been reached), otherwise continue to next iteration
      # (use the fact that the C++ function updates fractions internally)
      remaining <- check_rounding(C, fractions, L, remaining, round.base, tol)
      if (remaining > 0) iter <- iter + 1
      
    }
    
    if (iter >= max.iter) stop('Maximum number of iterations has been reached!')
    if (verbose) cat(sprintf('Iteration %d - A controlled rounding has been obtained\n', iter))
    
    C <- round.base * round(C / round.base)
    C0 <- C[1:r,1:k]
    
    # inner part and grand total of A: round down if C == 0, round up if C == 1
    A <- round.base * floor(A / round.base)
    A[C0 == round.base] <- A[C0 == round.base] + round.base
    A.gt <- round.base * floor(A.gt / round.base)
    if (C[r+1,k+1] == round.base) A.gt <- A.gt + round.base
    
    # row and column sums of A: round down if C == 1, round up if C == 0
    A.rs <- round.base * floor(A.rs / round.base)
    A.rs[C[-(r+1),k+1] == 0] <- A.rs[C[-(r+1),k+1] == 0] + round.base
    A.cs <- round.base * floor(A.cs / round.base)
    A.cs[C[r+1,-(k+1)] == 0] <- A.cs[C[r+1,-(k+1)] == 0] + round.base
    
    if (return.margins) {
      res <- rbind(cbind(A, A.rs), c(A.cs, A.gt))
      row.names(res) <- c(1:r,'colSums')
      colnames(res) <- c(1:k,'rowSums')
    } else {
      res <- A
    }
    
    return(res)
    
  }
  
  #############################################################
  # CHECKING WHETHER DATA SATISFY EDITS
  #############################################################
  EditsCheck <- function(intcurrentdata) {
    intsatisfied <- TRUE
    for (i in (1:nrows)) {
      testEdits <- substValue(ExpEdits,"Geslacht", intcurrentdata$Geslacht[i], reduce = FALSE)  
      testEdits <- substValue(testEdits,"Leeftijd", intcurrentdata$Leeftijd[i], reduce = FALSE)  
      testEdits <- substValue(testEdits,"HH_Pos", intcurrentdata$HH_Pos[i], reduce = FALSE)  
      testEdits <- substValue(testEdits,"HH_grootte", intcurrentdata$HH_grootte[i], reduce = FALSE) 
      testEdits <- substValue(testEdits,"Woonregio.vorig.jaar", intcurrentdata$Woonregio.vorig.jaar[i], reduce = FALSE)  
      testEdits <- substValue(testEdits,"Nationaliteit", intcurrentdata$Nationaliteit[i], reduce = FALSE)  
      testEdits <- substValue(testEdits,"Geboorteland", intcurrentdata$Geboorteland[i], reduce = FALSE)  
      testEdits <- substValue(testEdits,"Onderwijsniveau", intcurrentdata$Onderwijsniveau[i], reduce = FALSE)  
      testEdits <- substValue(testEdits,"Econ..status", intcurrentdata$Econ..status[i], reduce = FALSE)  
      testEdits <- substValue(testEdits,"Beroep", intcurrentdata$Beroep[i], reduce = FALSE)  
      testEdits <- substValue(testEdits,"SBI", intcurrentdata$SBI[i], reduce = FALSE)  
      testEdits <- substValue(testEdits,"Burg..Staat", intcurrentdata$Burg..Staat[i], reduce = FALSE)  
      testsatisfied <- isFeasible(testEdits)
      if (!testsatisfied) {
        intsatisfied <- FALSE
      }
    }
    return(intsatisfied)
  }
  
  #############################################################
  # CHECKING WHETHER DATA SATISFY TOTALS
  #############################################################
  TotalsCheck <- function(intcurrentdata) {
    intsatisfied <- TRUE
    finaltotals <- matrix(NA,nvars,max_ncats)
    for (j in (1:nvars)) {
      if (!is.na(totals[j,1])) {
        myfinalTable <- table(intcurrentdata[,j + 1])
        for (k in (1:ncats[j])) {
          finaltotals[j,k] <- margin.table(myfinalTable,1)[k]
          if (finaltotals[j,k] != totals[j,k]) {
            intsatisfied <- FALSE
          }
        }
      }
    }
    rm(finaltotals)
    gc()
    return(intsatisfied)
  }

  #############################################################
  # CHECKING WHETHER ALL MISSING VALUES ARE REALLY IMPUTED
  #############################################################
  NAsCheck <- function(intcurrentdata) {
    intnoNAs <- TRUE
    for (i in (1:nrows)) {
      for (j in (1:nvars)) {
        if (is.na(intcurrentdata[i,j])) {
          intnoNAs <- FALSE
        }
      }
    }
    return(intnoNAs)
  }
  
  ###############################
  # FUNCTIONS FOR CALCULATING QUALITY MEASURES
  ###############################
  Determine_Quality <- function(intapproxmethod, intbias, intbias_8_10, intbias_8_2, intbias_8_9, intbias_10_2, intbias_10_9, intCorrectlyImputed) {
    for (j in (1:nvars)) {
      mytrueTable <- table(truedata[,j + 1])
      myimpTable <- table(currentdata[,j + 1])
      for (k in (1:ncats[j])) {
        intdiff <- margin.table(myimpTable,1)[k] - margin.table(mytrueTable,1)[k]
        intbias[intapproxmethod, j, k] <- intdiff + intbias[intapproxmethod, j, k]
        rm(intdiff)
        gc()
      }
      rm(mytrueTable)
      rm(myimpTable)
      gc()
    }
    mytruecrossTable <- table(truedata$Onderwijsniveau, truedata$Beroep)
    myimpcrossTable <- table(currentdata$Onderwijsniveau, truedata$Beroep)
    diffcrosstables <- myimpcrossTable - mytruecrossTable
    for (ii in (1:ncats[8])) {
      for (jj in (1:ncats[10])) {
        intbias_8_10[intapproxmethod, ii,jj] <- diffcrosstables[ii,jj] + intbias_8_10[intapproxmethod, ii,jj]
      }
    }
    mytruecrossTable <- table(truedata$Onderwijsniveau, truedata$Leeftijd)
    myimpcrossTable <- table(currentdata$Onderwijsniveau, truedata$Leeftijd)
    diffcrosstables <- myimpcrossTable - mytruecrossTable
    for (ii in (1:ncats[8])) {
      for (jj in (1:ncats[2])) {
        intbias_8_2[intapproxmethod, ii,jj] <- diffcrosstables[ii,jj] + intbias_8_2[intapproxmethod, ii,jj]
      }
    }
    mytruecrossTable <- table(truedata$Onderwijsniveau, truedata$Econ..status)
    myimpcrossTable <- table(currentdata$Onderwijsniveau, truedata$Econ..status)
    diffcrosstables <- myimpcrossTable - mytruecrossTable
    for (ii in (1:ncats[8])) {
      for (jj in (1:ncats[9])) {
        intbias_8_9[intapproxmethod, ii,jj] <- diffcrosstables[ii,jj] + intbias_8_9[intapproxmethod, ii,jj]
      }
    }
    mytruecrossTable <- table(truedata$Beroep, truedata$Leeftijd)
    myimpcrossTable <- table(currentdata$Beroep, truedata$Leeftijd)
    diffcrosstables <- myimpcrossTable - mytruecrossTable
    for (ii in (1:ncats[10])) {
      for (jj in (1:ncats[2])) {
        intbias_10_2[intapproxmethod, ii,jj] <- diffcrosstables[ii,jj] + intbias_10_2[intapproxmethod, ii,jj]
      }
    }
    mytruecrossTable <- table(truedata$Beroep, truedata$Econ..status)
    myimpcrossTable <- table(currentdata$Beroep, truedata$Econ..status)
    diffcrosstables <- myimpcrossTable - mytruecrossTable
    for (ii in (1:ncats[10])) {
      for (jj in (1:ncats[9])) {
        intbias_10_9[intapproxmethod, ii,jj] <- diffcrosstables[ii,jj] + intbias_10_9[intapproxmethod, ii,jj]
      }
    }
    for (j in (1:nvars)) {
      for (i in (1:num_missings[j])) {
        if (currentdata[missing_indices[j,i],j + 1] == truedata[missing_indices[j,i],j + 1]) {
          intCorrectlyImputed[intapproxmethod, j] <- intCorrectlyImputed[intapproxmethod, j] + 1  
        }
      }
    }
    list("bias"=intbias, "bias_8_10"=intbias_8_10, "bias_8_2"=intbias_8_2, "bias_8_9"=intbias_8_9, 
         "bias_10_2"=intbias_10_2, "bias_10_9"=intbias_10_9, "CorrectlyImputed"=intCorrectlyImputed)
  }
  
  #############################################################
  # DETERMINE MISSING VALUES IN THE DATA
  #############################################################
    Determine_missings <- function(intoriginaldata) {
    intmissings <- matrix(FALSE, nrows, nvars)
    intnum_missings <- rep(0, nvars)
    intnummissingsrec <- rep(0,nrows)
    intmissing_indices <- matrix(0, nvars, max_num_missings)
    intnummissingsrecinedits <- rep(0,nrows)
    for (i in (1:nrows)) {
      for (j in (1:nvars)) {
        if (is.na(intoriginaldata[i,j + 1])) {
          intmissings[i,j] <- TRUE
          intnum_missings[j] <- intnum_missings[j] + 1
          intmissing_indices[j,intnum_missings[j]] <- i
          intnummissingsrec[i] <- intnummissingsrec[i] + 1
          if ((j %in% (2:4)) | (j %in% (8:12))) {
            intnummissingsrecinedits[i] <- intnummissingsrecinedits[i] + 1
          }
        } 
      }
    }
    return(list("missings"=intmissings, "num_missings"=intnum_missings, "missing_indices"=intmissing_indices, "nummissingsrec"=intnummissingsrec, "nummissingsrecinedits"=intnummissingsrecinedits))
  }
  
  
  #############################################################
  # DETERMING FRACTIONS (PROPORTIONS) IN THE DATA AND NUMBER OF UNITS TO BE IMPUTED
  #############################################################
    Determine_fractions_imptotals <- function(intoriginaldata) {
    intimptotals <- matrix(NA,nvars,max_ncats)
    intfractions <- matrix(0,nvars,max_ncats)
    for (j in (1:nvars)) {
      myTable <- table(intoriginaldata[,j + 1])
      intsum <- 0
      for (k in (1:ncats[j])) {
        intsum <- intsum + margin.table(myTable,1)[k]
      }
      for (k in (1:ncats[j])) {
        intfractions[j,k] <- margin.table(myTable,1)[k]/intsum
      }
      if (!is.na(totals[j,1])) {
        for (k in (1:ncats[j])) {
          intimptotals[j,k] <- totals[j,k] - margin.table(myTable,1)[k]
        }
      }
    }
    return(list("fractions"=intfractions,"imptotals"=intimptotals))
  }
  
  
  ##############################################################################
  # MAIN PROGRAM
  ##############################################################################
  
  ###########################
  # SETTING PARAMETERS
  ###########################
  IPFMaxCriterion <- 1e-10
  IterDiff <- 1
  
  ###########################
  # MaxIteR: NUMBER OF ITERATIONS; MissFactor: 1 LOW MISSINGNESS, 2 HIGH MISSINGNESS, s_done: NUMBER OF DATASETS ALREADY IMPUTED IN A PREVIOUS RUN (USEFUL IN CASE SIMULATION WAS SOMEHOW INTERRUPTED)
  ###########################
  MaxIter <- 5
  Missfactor <- 1
  s_done <- 325
  
  ###########################
  # FOR MAR AND NMAR MECHANISMS: SELECT VARIABLE IN WHICH MAR OR NMAR MISSINGNESS IS INTRODUCED
  ###########################
  # varname <- "edu"
  # varname <- "occ"


  ###########################
  # SELECTING KIND OF MISSINGNESS MECHANISM: MCAR, MAR OR NMAR
  ###########################  
  originaldata <- Create_Missings(Missfactor)
  # originaldata <- Create_MAR_Selective_Missings(Missfactor, varname)  
  # originaldata <- Create_NMAR_Selective_Missings(Missfactor, varname)   
 
  ###########################
  # INITIALIZING QUALITY MEASURES
  ###########################
  bias <- array(0, dim=c(3, nvars, max_ncats))
  variance <- array(0,dim=c(3, nvars, max_ncats))
  bias_8_10 <- array(0, dim=c(3, ncats[8], ncats[10]))
  var_8_10  <- array(0, dim=c(3, ncats[8], ncats[10]))
  bias_8_2  <- array(0, dim=c(3, ncats[8], ncats[2]))
  var_8_2  <- array(0, dim=c(3, ncats[8], ncats[2]))
  bias_8_9  <- array(0, dim=c(3, ncats[8], ncats[9]))
  var_8_9  <- array(0, dim=c(3, ncats[8], ncats[9]))
  bias_10_2  <- array(0, dim=c(3, ncats[10], ncats[2]))
  var_10_2  <- array(0, dim=c(3, ncats[10], ncats[2]))
  bias_10_9  <- array(0, dim=c(3, ncats[10], ncats[9]))
  var_10_9  <- array(0, dim=c(3, ncats[10], ncats[9]))
  CorrectlyImputed <- matrix(0, 3, nvars)
  FracCorrectlyImputed <- matrix(0, 3, nvars)
 
  ########################
  # DETERMINING MISSINGS
  ########################
  missinglist <- Determine_missings(originaldata)
  missings <- missinglist$missings
  num_missings <- missinglist$num_missings
  nummissingsrec <- missinglist$nummissingsrec
  missing_indices <- missinglist$missing_indices
  nummissingsrecinedits <- missinglist$nummissingsrecinedits
  
  fractions_imptotalslist <- Determine_fractions_imptotals(originaldata)
  fractions <- fractions_imptotalslist$fractions
  imptotals <- fractions_imptotalslist$imptotals

##############################################################
# START-UP PHASE
##############################################################
  cat("Start of start-up phase ", file=logfilename, append = TRUE)
  cat(paste0("Dataset with missings ", s), file=logfilename, append = TRUE)
  cat("Start of start-up phase ")
  cat(paste0("Dataset with missings ", s))
  
  ########################
  # INITIAL DATA: DATASET WIH MISSING DATA
  ########################
  initialdata <- originaldata
  ###############################
  # FOR EACH RECORD: SUBSTITURE OBSERVED VALUES INTO EXPLICIT EDITS AND ELIMINIATE VARIBALES WITH MISSING VALUES FROM SET OF EXPLICIT EDITS
  ###############################
  for (i in (1:nrows)) {
    if (nummissingsrec[i] > 0) {
      E0 <- ExpEdits
      E1 <- ReduceEdits(i, 1, E0, nummissingsrec[i])
      E2 <- ReduceEdits(i, 2, E1, nummissingsrec[i]) 
      E3 <- ReduceEdits(i, 3, E2, nummissingsrec[i]) 
      E4 <- ReduceEdits(i, 4, E3, nummissingsrec[i]) 
      E5 <- ReduceEdits(i, 5, E4, nummissingsrec[i]) 
      E6 <- ReduceEdits(i, 6, E5, nummissingsrec[i]) 
      E7 <- ReduceEdits(i, 7, E6, nummissingsrec[i]) 
      E8 <- ReduceEdits(i, 8, E7, nummissingsrec[i]) 
      E9 <- ReduceEdits(i, 9, E8, nummissingsrec[i]) 
      E10 <- ReduceEdits(i, 10, E9, nummissingsrec[i]) 
      E11 <- ReduceEdits(i, 11, E10, nummissingsrec[i]) 
      if (missings[i,12]) {
        initialdata$Burg..Staat[i] <- ImputeUntilSatisfied(12, E11)
      }
      E10 <- substValue(E10, "Burg..Staat", initialdata$Burg..Staat[i], reduce = TRUE)
      E9 <- substValue(E9, "Burg..Staat", initialdata$Burg..Staat[i], reduce = TRUE)
      E8 <- substValue(E8, "Burg..Staat", initialdata$Burg..Staat[i], reduce = TRUE)
      E7 <- substValue(E7, "Burg..Staat", initialdata$Burg..Staat[i], reduce = TRUE)
      E6 <- substValue(E6, "Burg..Staat", initialdata$Burg..Staat[i], reduce = TRUE)
      E5 <- substValue(E5, "Burg..Staat", initialdata$Burg..Staat[i], reduce = TRUE)
      E4 <- substValue(E4, "Burg..Staat", initialdata$Burg..Staat[i], reduce = TRUE)
      E3 <- substValue(E3, "Burg..Staat", initialdata$Burg..Staat[i], reduce = TRUE)
      E2 <- substValue(E2, "Burg..Staat", initialdata$Burg..Staat[i], reduce = TRUE)
      E1 <- substValue(E1, "Burg..Staat", initialdata$Burg..Staat[i], reduce = TRUE)
      E0 <- substValue(E0, "Burg..Staat", initialdata$Burg..Staat[i], reduce = TRUE)
      if (missings[i,11]) {
        initialdata$SBI[i] <- ImputeUntilSatisfied(11, E10)
      }
      E9 <- substValue(E9, "SBI", initialdata$SBI[i], reduce = TRUE)
      E8 <- substValue(E8, "SBI", initialdata$SBI[i], reduce = TRUE) 
      E7 <- substValue(E7, "SBI", initialdata$SBI[i], reduce = TRUE) 
      E6 <- substValue(E6, "SBI", initialdata$SBI[i], reduce = TRUE) 
      E5 <- substValue(E5, "SBI", initialdata$SBI[i], reduce = TRUE) 
      E4 <- substValue(E4, "SBI", initialdata$SBI[i], reduce = TRUE) 
      E3 <- substValue(E3, "SBI", initialdata$SBI[i], reduce = TRUE) 
      E2 <- substValue(E2, "SBI", initialdata$SBI[i], reduce = TRUE) 
      E1 <- substValue(E1, "SBI", initialdata$SBI[i], reduce = TRUE) 
      E0 <- substValue(E0, "SBI", initialdata$SBI[i], reduce = TRUE) 
      if (missings[i,10]) {
        initialdata$Beroep[i] <- ImputeUntilSatisfied(10, E8)
      }
      E8 <- substValue(E8, "Beroep", initialdata$Beroep[i], reduce = TRUE)
      E7 <- substValue(E7, "Beroep", initialdata$Beroep[i], reduce = TRUE)
      E6 <- substValue(E6, "Beroep", initialdata$Beroep[i], reduce = TRUE)
      E5 <- substValue(E5, "Beroep", initialdata$Beroep[i], reduce = TRUE)
      E4 <- substValue(E4, "Beroep", initialdata$Beroep[i], reduce = TRUE)
      E3 <- substValue(E3, "Beroep", initialdata$Beroep[i], reduce = TRUE)
      E2 <- substValue(E2, "Beroep", initialdata$Beroep[i], reduce = TRUE)
      E1 <- substValue(E1, "Beroep", initialdata$Beroep[i], reduce = TRUE)
      E0 <- substValue(E0, "Beroep", initialdata$Beroep[i], reduce = TRUE)
      if (missings[i,9]) {
        initialdata$Econ..status[i] <- ImputeUntilSatisfied(9, E8)
      }
      E7 <- substValue(E7, "Econ..status", initialdata$Econ..status[i], reduce = TRUE)
      E6 <- substValue(E6, "Econ..status", initialdata$Econ..status[i], reduce = TRUE)
      E5 <- substValue(E5, "Econ..status", initialdata$Econ..status[i], reduce = TRUE)
      E4 <- substValue(E4, "Econ..status", initialdata$Econ..status[i], reduce = TRUE)
      E3 <- substValue(E3, "Econ..status", initialdata$Econ..status[i], reduce = TRUE)
      E2 <- substValue(E2, "Econ..status", initialdata$Econ..status[i], reduce = TRUE)
      E1 <- substValue(E1, "Econ..status", initialdata$Econ..status[i], reduce = TRUE)
      E0 <- substValue(E0, "Econ..status", initialdata$Econ..status[i], reduce = TRUE)
      if (missings[i,8]) {
        initialdata$Onderwijsniveau[i] <- ImputeUntilSatisfied(8, E7)
      }
      E6 <- substValue(E6, "Onderwijsniveau", initialdata$Onderwijsniveau[i], reduce = TRUE)
      E5 <- substValue(E5, "Onderwijsniveau", initialdata$Onderwijsniveau[i], reduce = TRUE)
      E4 <- substValue(E4, "Onderwijsniveau", initialdata$Onderwijsniveau[i], reduce = TRUE)
      E3 <- substValue(E3, "Onderwijsniveau", initialdata$Onderwijsniveau[i], reduce = TRUE)
      E2 <- substValue(E2, "Onderwijsniveau", initialdata$Onderwijsniveau[i], reduce = TRUE)
      E1 <- substValue(E1, "Onderwijsniveau", initialdata$Onderwijsniveau[i], reduce = TRUE)
      E0 <- substValue(E0, "Onderwijsniveau", initialdata$Onderwijsniveau[i], reduce = TRUE)
      if (missings[i,7]) {
        initialdata$Geboorteland[i] <- ImputeUntilSatisfied(7, E6)
      }
      E5 <- substValue(E5, "Geboorteland", initialdata$Geboorteland[i], reduce = TRUE)
      E4 <- substValue(E4, "Geboorteland", initialdata$Geboorteland[i], reduce = TRUE)
      E3 <- substValue(E3, "Geboorteland", initialdata$Geboorteland[i], reduce = TRUE)
      E2 <- substValue(E2, "Geboorteland", initialdata$Geboorteland[i], reduce = TRUE)
      E1 <- substValue(E1, "Geboorteland", initialdata$Geboorteland[i], reduce = TRUE)
      E0 <- substValue(E0, "Geboorteland", initialdata$Geboorteland[i], reduce = TRUE)
      if (missings[i,6]) {
        initialdata$Nationaliteit[i] <- ImputeUntilSatisfied(6, E5)
      }
      E4 <- substValue(E4, "Nationaliteit", initialdata$Nationaliteit[i], reduce = TRUE)
      E3 <- substValue(E3, "Nationaliteit", initialdata$Nationaliteit[i], reduce = TRUE)
      E2 <- substValue(E2, "Nationaliteit", initialdata$Nationaliteit[i], reduce = TRUE)
      E1 <- substValue(E1, "Nationaliteit", initialdata$Nationaliteit[i], reduce = TRUE)
      E0 <- substValue(E0, "Nationaliteit", initialdata$Nationaliteit[i], reduce = TRUE)
      if (missings[i,5]) {
        initialdata$Woonregio.vorig.jaar[i] <- ImputeUntilSatisfied(5, E4)
      }
      E3 <- substValue(E3, "Woonregio.vorig.jaar", initialdata$Woonregio.vorig.jaar[i], reduce = TRUE)
      E2 <- substValue(E2, "Woonregio.vorig.jaar", initialdata$Woonregio.vorig.jaar[i], reduce = TRUE)
      E1 <- substValue(E1, "Woonregio.vorig.jaar", initialdata$Woonregio.vorig.jaar[i], reduce = TRUE)
      E0 <- substValue(E0, "Woonregio.vorig.jaar", initialdata$Woonregio.vorig.jaar[i], reduce = TRUE)
      if (missings[i,4]) {
        initialdata$HH_grootte[i] <- ImputeUntilSatisfied(4, E3)
      }
      E2 <- substValue(E2, "HH_grootte", initialdata$HH_grootte[i], reduce = TRUE)
      E1 <- substValue(E1, "HH_grootte", initialdata$HH_grootte[i], reduce = TRUE)
      E0 <- substValue(E0, "HH_grootte", initialdata$HH_grootte[i], reduce = TRUE)
      if (missings[i,3]) {
        initialdata$HH_Pos[i] <- ImputeUntilSatisfied(3, E2)
      }
      E1 <- substValue(E1, "HH_Pos", initialdata$HH_Pos[i], reduce = TRUE)
      E0 <- substValue(E0, "HH_Pos", initialdata$HH_Pos[i], reduce = TRUE)
      if (missings[i,2]) {
        initialdata$Leeftijd[i] <- ImputeUntilSatisfied(2, E1)
      }
      E0 <- substValue(E1,"Leeftijd", initialdata$Leeftijd[i], reduce = TRUE) 
      if (missings[i,1]) {
        initialdata$Geslacht[i] <- ImputeUntilSatisfied(1, E0)
      }
    }
  }

################################################
# IMPUTATION PHASE
################################################
  ########################
  # LOOP OVER APPROXIMATION METHODS
  ########################
    for (approxmethod in (1:3)) {
    currentdata <- initialdata
    num_Iter <- 0

    #################
    # LOOP OVER ITERATIONS OF THE IMPUTATION APPPROACH
    #################
    while (num_Iter < MaxIter) {
      num_Iter <- num_Iter + 1

    ########################
    # LOOP OVER VARIABLES
    ########################
      for (j in (1:nvars)) {
        cat(paste0("Dataset with missings ", s), file=logfilename, append = TRUE)
        cat(paste0("Approximation method ", approxmethod), file=logfilename, append = TRUE)
        cat(paste0("Start iteration ", num_Iter), file=logfilename, append = TRUE)
        cat(paste0("Variable ", j), file=logfilename, append = TRUE)
        
        probs <- matrix(0, nrows, ncats[j])
        probs_small <- matrix(0, num_missings[j], ncats[j])
        ImpMatrix <- matrix(0, num_missings[j], ncats[j])
        probs <- EstimateProbabilities(j)
        probs_small <- AdjustProbs_StructZeroes_AllRecs(probs, j)
        for (i in (1:num_missings[j])) {
          sum_probs <- 0
          for (k in (1:ncats[j])) {
            if (probs_small[i,k] < 1e-10) {
              probs_small[i,k] <- 0
            }
            if (is.na(probs_small[i,k])) {
              probs_small[i,k] <- 0
            }
            sum_probs <- sum_probs + probs_small[i,k]
          }
          if ((sum_probs < 1e-10) | (sum_probs > (1 + 1e-10))) {
            for (k in (1:ncats[j])) {
              probs_small[i,k] <- 1/ncats[j]
            }
          }  
        }
        if (!is.na(totals[j,1])) {
          probs_approx <- matrix(0, num_missings[j], ncats[j])
          probs_approx <- ApproximateProbs(probs_small, j, approxmethod)
          for (i in (1:num_missings[j])) {
            largestprob <- 0
            largestprobindex <- 0
            sum_probs <- 0
            for (k in (1:ncats[j])) {
              probs_approx[i,k] <- round(probs_approx[i,k], digits = 5)
              if (probs_approx[i,k] > largestprob) {
                largestprobindex <- k
              }
            }
            if (largestprobindex > 1) {
              for (kk in (1:(largestprobindex - 1))) {
                sum_probs <- sum_probs + probs_approx[i,kk]
              }  
            }
            if (largestprobindex < ncats[j]) {
              for (kk in ((largestprobindex + 1):ncats[j])) {
                sum_probs <- sum_probs + probs_approx[i,kk]
              }
            }
            probs_approx[i, largestprobindex] <- 1 - sum_probs
            probs_approx[i, largestprobindex] <- round(probs_approx[i, largestprobindex], digits = 5)
          } 
          RoundedMatrix <- FALSE
          while (!RoundedMatrix) {
            RoundedMatrix <- TRUE
            ImpMatrix <- roundCox_cpp(A = probs_approx, verbose = TRUE)
#            ImpMatrix <- roundCox_cpp(A = probs_approx, cycle = 'first', verbose = TRUE)
            for (ii in (1:num_missings[j]))
            {
              testrows <- 0
              for (kk in (1:ncats[j])) {
                testrows <- testrows + ImpMatrix[ii, kk]
              }
              if (abs(testrows - 1) > 0.5) {
                RoundedMatrix <- FALSE
              }
            }
            for (kk in (1:ncats[j]))
            {
              testcolumns <- 0
              for (ii in (1:num_missings[j])) {
                testcolumns <- testcolumns + ImpMatrix[ii, kk]
              }
              if (abs(testcolumns - imptotals[j, kk]) > 0.5) {
                RoundedMatrix <- FALSE
              }
            }
          }
          for (ii in (1:num_missings[j])) {
            for (k in (1:ncats[j])) {
              if (ImpMatrix[ii, k] == 1) {
                currentdata[missing_indices[j, ii], j + 1] <- codes[j, k]
              }
            }
          }
          rm(probs_approx)
          gc()
        }
        else {
          for (ii in (1:num_missings[j])) {
            sum_probs <- 0
            for (kk in (1:ncats[j])) {
              sum_probs <- sum_probs + probs_small[ii,kk]
            }
            if (is.na (sum_probs)) {
              sum_probs <- 0
            }
            if (sum_probs < 1e-10) {
              for (kk in (1:ncats[j])) {
                probs_small[ii,kk] <- 1/ncats[j]
              }
              probs_small[ii,] <- AdjustProbs_StructZeroes_Rec(probs, j)
            }
            currentdata[missing_indices[j, ii], j + 1] <- codes[j, SimpleImpute(j, probs_small[ii,])]
          }  
        }
        rm(ImpMatrix)
        rm(probs)
        rm(probs_small)
        gc()
      }
      EditsSatisfied <- EditsCheck(currentdata)
      if (!EditsSatisfied) {
        cat("Edits failed", file=logfilename, append = TRUE)
      }
      TotalsSatisfied <- TotalsCheck(currentdata)
      if (!TotalsSatisfied) {
        cat("Totals failed", file=logfilename, append = TRUE)
      }
      NAsSatisfied <- NAsCheck(currentdata)
      if (!NAsSatisfied) {
        cat("NAs in data", file=logfilename, append = TRUE)
      }
    }
    Qualitylist <- Determine_Quality(approxmethod, bias, bias_8_10, bias_8_2, bias_8_9, bias_10_2, bias_10_9, CorrectlyImputed)
    bias[approxmethod,,] <- Qualitylist$bias[approxmethod,,]
    bias_8_10[approxmethod,,] <- Qualitylist$bias_8_10[approxmethod,,]
    bias_8_2[approxmethod,,] <- Qualitylist$bias_8_2[approxmethod,,]
    bias_8_9[approxmethod,,] <- Qualitylist$bias_8_9[approxmethod,,]
    bias_10_2[approxmethod,,] <- Qualitylist$bias_10_2[approxmethod,,]
    bias_10_9[approxmethod,,] <- Qualitylist$bias_10_9[approxmethod,,]
    CorrectlyImputed[approxmethod,] <- Qualitylist$CorrectlyImputed[approxmethod,]
  }
  filename <- paste0(pathresults,"bias/bias_",s+s_done,".csv")
  for (approxmethod in (1:3)) {
    if (approxmethod == 1) {
      write(paste0("approxmethod: ", approxmethod), filename)        
    }
    else {
      write(paste0("approxmethod: ", approxmethod), filename, append = TRUE)
    }
    write.table(bias[approxmethod,,], filename, append = TRUE, row.names = FALSE, col.names = FALSE, sep = ";")
    if (approxmethod != 3) {
      write.table("", filename, append = TRUE, row.names = FALSE, col.names = FALSE, sep=";")    
    }
  }
  filename <- paste0(pathresults,"bias_8_10/bias_8_10_",s+s_done,".csv")
  for (approxmethod in (1:3)) {
    if (approxmethod == 1) {
      write(paste0("approxmethod: ", approxmethod), filename)        
    }
    else {
      write(paste0("approxmethod: ", approxmethod), filename, append = TRUE)
    }
    write.table(bias_8_10[approxmethod,,], filename, append = TRUE, row.names = FALSE, col.names = FALSE, sep = ";")
    if (approxmethod != 3) {
      write.table("",filename, append = TRUE, row.names = FALSE, col.names = FALSE, sep=";")    
    }
  }
  filename <- paste0(pathresults,"bias_8_2/bias_8_2_",s+s_done,".csv")
  for (approxmethod in (1:3)) {
    if (approxmethod == 1) {
      write(paste0("approxmethod: ", approxmethod), filename)        
    }
    else {
      write(paste0("approxmethod: ", approxmethod), filename, append = TRUE)
    }
    write.table(bias_8_2[approxmethod,,], filename, append = TRUE, row.names = FALSE, col.names = FALSE, sep = ";")
    if (approxmethod != 3) {
      write.table("",filename, append = TRUE, row.names = FALSE, col.names = FALSE, sep=";")    
    }
  }
  filename <- paste0(pathresults,"bias_8_9/bias_8_9_",s+s_done,".csv")
  for (approxmethod in (1:3)) {
    if (approxmethod == 1) {
      write(paste0("approxmethod: ", approxmethod), filename)        
    }
    else {
      write(paste0("approxmethod: ", approxmethod), filename, append = TRUE)
    }
    write.table(bias_8_9[approxmethod,,], filename, append = TRUE, row.names = FALSE, col.names = FALSE, sep = ";")
    if (approxmethod != 3) {
      write.table("",filename, append = TRUE, row.names = FALSE, col.names = FALSE, sep=";")    
    }
  }
  filename <- paste0(pathresults,"bias_10_2/bias_10_2_",s+s_done,".csv")
  for (approxmethod in (1:3)) {
    if (approxmethod == 1) {
      write(paste0("approxmethod: ", approxmethod), filename)        
    }
    else {
      write(paste0("approxmethod: ", approxmethod), filename, append = TRUE)
    }
    write.table(bias_10_2[approxmethod,,], filename, append = TRUE, row.names = FALSE, col.names = FALSE, sep = ";")
    if (approxmethod != 3) {
      write.table("",filename, append = TRUE, row.names = FALSE, col.names = FALSE, sep=";")    
    }
  }
  filename <- paste0(pathresults,"bias_10_9/bias_10_9_",s+s_done,".csv")
  for (approxmethod in (1:3)) {
    if (approxmethod == 1) {
      write(paste0("approxmethod: ", approxmethod), filename)        
    }
    else {
      write(paste0("approxmethod: ", approxmethod), filename, append = TRUE)
    }
    write.table(bias_10_9[approxmethod,,], filename, append = TRUE, row.names = FALSE, col.names = FALSE, sep = ";")
    if (approxmethod != 3) {
      write.table("",filename, append = TRUE, row.names = FALSE, col.names = FALSE, sep=";")    
    }
  }    
  for (j in (1:nvars)) {
    for (approxmethod in (1:3)) {
      FracCorrectlyImputed[approxmethod,j] <- round(CorrectlyImputed[approxmethod,j]/num_missings[j], digits = 4)    
    }
  }
  filename <- paste0(pathresults,"CorrectlyImputed/CorrectlyImputed_",s+s_done,".csv")
  for (approxmethod in (1:3)) {
    if (approxmethod == 1) {
      write(paste0("approxmethod: ", approxmethod), filename)        
    }
    else {
      write(paste0("approxmethod: ", approxmethod), filename, append = TRUE)
    }
    write.table(CorrectlyImputed[approxmethod,], filename, append = TRUE, row.names = FALSE, col.names = FALSE, sep = ";")
    write.table("",filename, append = TRUE, row.names = FALSE, col.names = FALSE, sep=";")
  }
  for (approxmethod in (1:3)) {
    write(paste0("approxmethod: ", approxmethod), filename, append = TRUE)
    write.table(FracCorrectlyImputed[approxmethod,], filename, append = TRUE, row.names = FALSE, col.names = FALSE, sep = ";")
    if (approxmethod != 3) {
      write.table("",filename, append = TRUE, row.names = FALSE, col.names = FALSE, sep=";")    
    }
  }    
}

stopCluster(myCluster)

################################################
path <- "//CBSP.NL/Productie/Secundair/MPOnderzoek/Werk/Combineren/Medewerkers/TWAL/Eigen projecten/Calibrated Imputation/Categorical data/R scripts/"
pathdata <- "//CBSP.NL/Productie/Secundair/MPOnderzoek/Werk/Combineren/Medewerkers/TWAL/Eigen projecten/Calibrated Imputation/Categorical data/Test data/Corrected data/"
pathresults <- "//CBSP.NL/Productie/Secundair/MPOnderzoek/Werk/Combineren/Medewerkers/TWAL/Eigen projecten/Calibrated Imputation/Categorical data/Simulation results/"

#############################################################################################################################
## READING TRUE DATA FOR THE NON-PARALLEL PART OF THE CODE, NEEDED TO DETERMINE THE TRUE TOTALS
#############################################################################################################################
nvars <- 12
max_ncats <- 17

ncats <- rep(0,nvars)
ncats[1] <- 2
ncats[2] <- 17
ncats[3] <- 8
ncats[4] <- 6
ncats[5] <- 3
ncats[6] <- 3
ncats[7] <- 3  
ncats[8] <- 7
ncats[9] <- 8
ncats[10] <- 10
ncats[11] <- 13  
ncats[12] <- 4

totals <- matrix(NA, nvars, max_ncats)

truedata <- read.csv2(paste(pathdata,"IPUMS klein corrected.csv", sep=""), sep=";", header=TRUE)
truedata$Geslacht <- as.factor(truedata$Geslacht)
truedata$Leeftijd <- as.factor(truedata$Leeftijd)
truedata$HH_Pos <- as.factor(truedata$HH_Pos)
truedata$HH_grootte <- as.factor(truedata$HH_grootte)
truedata$Woonregio.vorig.jaar <- as.factor(truedata$Woonregio.vorig.jaar)
truedata$Nationaliteit <- as.factor(truedata$Nationaliteit)
truedata$Geboorteland <- as.factor(truedata$Geboorteland)
truedata$Onderwijsniveau <- as.factor(truedata$Onderwijsniveau)
truedata$Econ..status <- as.factor(truedata$Econ..status)
truedata$Beroep <- as.factor(truedata$Beroep)
truedata$SBI <- as.factor(truedata$SBI)
truedata$Burg..Staat <- as.factor(truedata$Burg..Staat)
# VARIABLEs WITH ONLY TWO CATEGORIES / SPECIFYING BASIC LEVEL
truedata$Geslacht <- relevel(truedata$Geslacht, ref="1")

Determine_totals <- function(intvarnr) {
  inttotalsintvarnr <- rep(NA, max_ncats)
  mytrueTable <- table(truedata[,intvarnr + 1])
  for (k in (1:ncats[intvarnr])) {
    inttotalsintvarnr[k] <- margin.table(mytrueTable,1)[k]
  }
  return(inttotalsintvarnr)
}

filename <- paste0(pathresults,"bias/truetotals.csv")
for (j in (1:nvars)) {
  totals[j,] <- Determine_totals(j)
}
write.table(totals[,], filename, row.names = FALSE, col.names = FALSE, sep = ";")
filename <- paste0(pathresults,"bias_8_10/truecrosstotals_8_10.csv") 
mytruecrossTable <- table(truedata$Onderwijsniveau, truedata$Beroep)
write.table(mytruecrossTable, filename, row.names = FALSE, col.names = FALSE, sep = ";")
rm(mytruecrossTable)
filename <- paste0(pathresults,"bias_8_2/truecrosstotals_8_2.csv") 
mytruecrossTable <- table(truedata$Onderwijsniveau, truedata$Leeftijd)
write.table(mytruecrossTable, filename, row.names = FALSE, col.names = FALSE, sep = ";")
rm(mytruecrossTable)
filename <- paste0(pathresults,"bias_8_9/truecrosstotals_8_9.csv") 
mytruecrossTable <- table(truedata$Onderwijsniveau, truedata$Econ..status)
write.table(mytruecrossTable, filename, row.names = FALSE, col.names = FALSE, sep = ";")
rm(mytruecrossTable)
filename <- paste0(pathresults,"bias_10_2/truecrosstotals_10_2.csv") 
mytruecrossTable <- table(truedata$Beroep, truedata$Leeftijd)
write.table(mytruecrossTable, filename, row.names = FALSE, col.names = FALSE, sep = ";")
rm(mytruecrossTable)
filename <- paste0(pathresults,"bias_10_9/truecrosstotals_10_9.csv") 
mytruecrossTable <- table(truedata$Beroep, truedata$Econ..status)
write.table(mytruecrossTable, filename, row.names = FALSE, col.names = FALSE, sep = ";")
rm(mytruecrossTable)
gc()

Sys.time() - init
