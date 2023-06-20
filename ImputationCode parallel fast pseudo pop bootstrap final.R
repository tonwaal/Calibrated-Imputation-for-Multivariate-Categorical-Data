#############################################################################################################################
# CODE FOR MASS IMPUTATION OF CATEGORICAL DATA GIVEN KNOWN TOTALS 
# PROJECT: Multivariate Mass Imputation for the Population Census Given Known Totals (2018-NL-ESS.VIP.BUS.ADMIN)  
# CALL: Multipurpose statistics for efficiency gains in production
#
# THIS CODE IS PARALLELIZED IN ORDER TO SPEED UP THE SIMULATION STUDY
#############################################################################################################################

# library(foreach)
# library(doParallel) # alleen voor doParallel

init <- Sys.time()
path <- "PATH TO WORK DIRECTORY"
pathdata <- "PATH TO FOLDER WITH DATA"
pathresults <- "PATH TO FOLDER WITH RESULTS"setwd(path)

########################
# LOOP OVER MISSING DATASETS
########################
# foreach(s = (1:NumMissingSets)) %dopar% {  # alleen voor doParallel
# install.packages("editrules")
library(parallel)
library(editrules)
library(stats)
library(nnet)
library(foreign)
library(Rcpp)

myCluster <- parallel::makeCluster(detectCores() - 1, type="PSOCK")

##############################################################################
# MAIN PROGRAM
##############################################################################

###########################
# SETTING PARAMETERS
###########################
IPFMaxCriterion <- 1e-10
IterDiff <- 1

########################
# NUMBER OF MISSING DATASETS FOR THIS PART OF SIMULATION STUDY
########################
NumMissingSets <- 6
# numPseudoPops <- 1
numBootstrapSamples <- 200 # FOR MASS IMPUTATION: BOOTSTRAP SAMPLE = PSEUDOPOPULATION

###########################
# MaxIteR: NUMBER OF ITERATIONS; MissFactor: 1 LOW MISSINGNESS, 2 HIGH MISSINGNESS, s_done: NUMBER OF DATASETS ALREADY IMPUTED IN A PREVIOUS RUN (USEFUL IN CASE SIMULATION WAS SOMEHOW INTERRUPTED)
###########################
MaxIter <- 10
Missfactor <- 1
s_done <- 247

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


###########################
# FOR MAR AND NMAR MECHANISMS: SELECT VARIABLE IN WHICH MAR OR NMAR MISSINGNESS IS INTRODUCED
###########################
# varname <- "edu"
# varname <- "occ"

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


Create_missings <- function(intmissfactor, intdata, intnrows, intnvars) {
  intoriginaldata <- intdata
  MissingFractions <- rep(0,intnvars)
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
  for (j in (1:intnvars)) {
    nummissings_j <- round(intnrows*MissingFractions[j])
    testsample <- sample((1:intnrows), nummissings_j, replace = FALSE)
    for (i in (1:nummissings_j)) {
      intoriginaldata[testsample[i],j + 1] <- NA 
    }
  }
  return(intoriginaldata)
}

Create_missings_OCC <- function(intmissfactor, intdata, intnrows, intnvars) {
  intoriginaldata <- intdata
  MissingFractions <- rep(0,intnvars)
  MissingFractions[1] <- 0
  MissingFractions[2] <- 0
  MissingFractions[3] <- 0
  MissingFractions[4] <- 0
  MissingFractions[5] <- 0
  MissingFractions[6] <- 0
  MissingFractions[7] <- 0
  MissingFractions[8] <- 0
  MissingFractions[9] <- 0
  MissingFractions[10] <- 0.25*intmissfactor
  MissingFractions[11] <- 0
  MissingFractions[12] <- 0
  for (j in (1:intnvars)) {
    nummissings_j <- round(intnrows*MissingFractions[j])
    testsample <- sample((1:intnrows), nummissings_j, replace = FALSE)
    for (i in (1:nummissings_j)) {
      intoriginaldata[testsample[i],j + 1] <- NA 
    }
  }
  return(intoriginaldata)
}

#############################################################
# DETERMINE MISSING VALUES IN THE DATA
#############################################################
Determine_missings <- function(intdata, intnrows, intnvars, intmax_num_missings) {
  intmissings <- matrix(FALSE, intnrows, intnvars)
  intnum_missings <- rep(0, intnvars)
  intnummissingsrec <- rep(0,intnrows)
  intmissing_indices <- matrix(0, intnvars, intmax_num_missings)
  intnummissingsrecinedits <- rep(0,intnrows)
  for (i in (1:intnrows)) {
    for (j in (1:intnvars)) {
      if (is.na(intdata[i,j + 1])) {
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
Determine_fractions_imptotals <- function(intdata, inttotals, intnrows, intnvars, intmax_ncats) {
  intimptotals <- matrix(NA,intnvars,intmax_ncats)
  intfractions <- matrix(0,intnvars,intmax_ncats)
  for (j in (1:intnvars)) {
    myTable <- table(intdata[,j + 1])
    intsum <- 0
    for (k in (1:ncats[j])) {
      intsum <- intsum + margin.table(myTable,1)[k]
    }
    for (k in (1:ncats[j])) {
      intfractions[j,k] <- margin.table(myTable,1)[k]/intsum
    }
    if (!is.na(inttotals[j,1])) {
      for (k in (1:ncats[j])) {
        intimptotals[j,k] <- inttotals[j,k] - margin.table(myTable,1)[k]
      }
    }
  }
  return(list("fractions"=intfractions,"imptotals"=intimptotals))
}

Create_PseudoPopulation <- function(int_nummissingsperrec, int_originaldata, intnrows) {
  out_originaldata <- int_originaldata
  int_complete_records <- rep(0,intnrows)
  int_numcompleterecs <- 0
  for (i in (1:intnrows)) {
    if (int_nummissingsperrec[i] == 0) {
      int_numcompleterecs <- int_numcompleterecs + 1
      int_complete_records[int_numcompleterecs] <- i
    }
  }
  int_numcopies <- (intnrows %/% int_numcompleterecs)
  int_numremainder <- (intnrows %% int_numcompleterecs)
  for (j  in (1:int_numcopies)) {
    for (k in (1:int_numcompleterecs)) {
      ii <- (j-1)*int_numcompleterecs + k
      out_originaldata[ii, ] <- int_originaldata[int_complete_records[k],]
    }
  }
  int_extrarecs <- sample((1:int_numcompleterecs), int_numremainder, replace = FALSE)
  for (jj in (1:int_numremainder)) {
    out_originaldata[int_numcopies*int_numcompleterecs + jj, ] <- int_originaldata[int_complete_records[jj], ]
  }
  return(out_originaldata)
}

#######################################################################################
# STATEMENTS AND FUNCTIONS FOR INITIALISING AND READING DATA (INCLUDING TOTALEN) AND EDITS, AND CREATING MISSINGS
#######################################################################################  

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

Determine_totals <- function(intvarnr, intdata, intmax_ncats) {
  inttotalsintvarnr <- rep(NA, intmax_ncats)
  myintTable <- table(intdata[,intvarnr + 1])
  for (k in (1:ncats[intvarnr])) {
    inttotalsintvarnr[k] <- margin.table(myintTable,1)[k]
  }
  return(inttotalsintvarnr)
}

#############################################################################################################################
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

ReduceEdits <- function(intrownr, intvarnr, intEditsInRec, intnummissingsrecinedits, intmissings, intdata) {
  intEditsInRecReduced <- intEditsInRec
  if ((intvarnr %in% (2:4)) | (intvarnr %in% (8:12))) {
    if (intmissings[intrownr, intvarnr]) {
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
      if (!intmissings[intrownr, intvarnr]) {
        if (intvarnr == 2) {
          intEditsInRecReduced <- substValue(intEditsInRec,"Leeftijd",intdata$Leeftijd[intrownr], reduce = TRUE)  
        }
        if (intvarnr == 3) {
          intEditsInRecReduced <- substValue(intEditsInRec, "HH_Pos",intdata$HH_Pos[intrownr], reduce = TRUE)
        }        
        if (intvarnr == 4) {
          intEditsInRecReduced <- substValue(intEditsInRec, "HH_grootte", intdata$HH_grootte[intrownr], reduce = TRUE)
        }        
        if (intvarnr == 8) {
          intEditsInRecReduced <- substValue(intEditsInRec, "Onderwijsniveau",intdata$Onderwijsniveau[intrownr], reduce = TRUE)
        }        
        if (intvarnr == 9) {
          intEditsInRecReduced <- substValue(intEditsInRec, "Econ..status",intdata$Econ..status[intrownr], reduce = TRUE)
        }        
        if (intvarnr == 10) {
          intEditsInRecReduced <- substValue(intEditsInRec, "Beroep",intdata$Beroep[intrownr], reduce = TRUE)
        }        
        if (intvarnr == 11) {
          intEditsInRecReduced <- substValue(intEditsInRec, "SBI",intdata$SBI[intrownr], reduce = TRUE)
        }        
        if (intvarnr == 12) {
          intEditsInRecReduced <- substValue(intEditsInRec, "Burg..Staat",intdata$Burg..Staat[intrownr], reduce = TRUE)
        }
      }
    }
  }
  else{
    if (!intmissings[intrownr, intvarnr]) {
      if (intvarnr == 1) {
        intEditsInRecReduced <- substValue(intEditsInRec,"Geslacht",intdata$Geslacht[intrownr], reduce = TRUE)  
      }
      if (intvarnr == 5) {
        intEditsInRecReduced <- substValue(intEditsInRec, "Woonregio.vorig.jaar",intdata$Woonregio.vorig.jaar[intrownr], reduce = TRUE)
      }        
      if (intvarnr == 6) {
        intEditsInRecReduced <- substValue(intEditsInRec, "Nationaliteit", intdata$Nationaliteit[intrownr], reduce = TRUE)
      }        
      if (intvarnr == 7) {
        intEditsInRecReduced <- substValue(intEditsInRec, "Geboorteland",intdata$Geboorteland[intrownr], reduce = TRUE)
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
EstimateProbabilities <- function(intvarnr, intcurrentdata) {
  maxiterations <- 1000
  # VARIABLE Geslacht  
  # FOR VARIABLES WITH ONLY TWO CATEGORIES WE NEED TO MAKE VECTOR WITH TWO PROBABITIES (ONE FOR EACH CATEGORY)
  if (intvarnr == 1) {
    logresults <- multinom(Geslacht ~ Leeftijd + HH_Pos + HH_grootte + Woonregio.vorig.jaar + Nationaliteit + Geboorteland + Onderwijsniveau + Econ..status + Beroep + SBI + Burg..Staat, maxit = maxiterations, data = intcurrentdata)
    pp2_Vrouw <- fitted(logresults)  # MODEL PRODUCES PROBABILITY  FOR "female" (CATEGORY 2 ; "Vrouw" IN DUTCH), SINCE BASIC LEVEL IS male (CATEGORY 1; "man" IN DUTCH)
    pp2_Man <- 1 - pp2_Vrouw
    pp <- cbind(pp2_Man, pp2_Vrouw)    
  }
  # VARIABLE Leeftijd
  if (intvarnr == 2) {
    logresults <- multinom(Leeftijd ~ Geslacht + HH_Pos + HH_grootte + Woonregio.vorig.jaar + Nationaliteit + Geboorteland + Onderwijsniveau + Econ..status + Beroep + SBI + Burg..Staat, maxit = maxiterations, data = intcurrentdata)
    pp <- fitted(logresults)
  }
  # VARIABLE HH_Pos
  if (intvarnr == 3) {
    logresults_HH_Pos <- multinom(HH_Pos ~ Geslacht + Leeftijd + HH_grootte + Woonregio.vorig.jaar + Nationaliteit + Geboorteland + Onderwijsniveau + Econ..status + Beroep + SBI + Burg..Staat, maxit = maxiterations, data = intcurrentdata)
    pp <- fitted(logresults_HH_Pos)
  }
  # VARIABLE HH_grootte
  if (intvarnr == 4) {
    logresults <- multinom(HH_grootte ~ Geslacht + Leeftijd + HH_Pos + Woonregio.vorig.jaar + Nationaliteit + Geboorteland + Onderwijsniveau + Econ..status + Beroep + SBI + Burg..Staat, maxit = maxiterations, data = intcurrentdata)
    pp <- fitted(logresults)  
  }
  # VARIABLE Woonregio vorig jaar
  if (intvarnr == 5) {
    logresults <- multinom(Woonregio.vorig.jaar ~ Geslacht + Leeftijd + HH_Pos + HH_grootte + Nationaliteit + Geboorteland + Onderwijsniveau + Econ..status + Beroep + SBI + Burg..Staat, maxit = maxiterations, data = intcurrentdata)
    pp <- fitted(logresults)
  }
  # VARIABLE Nationaliteit
  if (intvarnr == 6) {
    logresults <- multinom(Nationaliteit ~ Geslacht + Leeftijd + HH_Pos + HH_grootte + Woonregio.vorig.jaar + Geboorteland + Onderwijsniveau + Econ..status + Beroep + SBI + Burg..Staat, maxit = maxiterations, data = intcurrentdata)
    pp <- fitted(logresults)
  }
  # VARIABLE Geboorteland
  if (intvarnr == 7) {
    logresults <- multinom(Geboorteland ~ Geslacht + Leeftijd + HH_Pos + HH_grootte + Woonregio.vorig.jaar + Nationaliteit + Onderwijsniveau + Econ..status + Beroep + SBI + Burg..Staat, maxit = maxiterations, data = intcurrentdata)
    pp <- fitted(logresults)
  }
  # VARIABLE Onderwijsniveau
  if (intvarnr == 8) {
    logresults <- multinom(Onderwijsniveau ~ Geslacht + Leeftijd + HH_Pos + HH_grootte + Woonregio.vorig.jaar + Nationaliteit + Geboorteland + Econ..status + Beroep + SBI + Burg..Staat, maxit = maxiterations, data = intcurrentdata)
    pp <- fitted(logresults)
  }
  # VARIABLE Econ..Status
  if (intvarnr == 9) {
    logresults <- multinom(Econ..status ~ Geslacht + Leeftijd + HH_Pos + HH_grootte + Woonregio.vorig.jaar + Nationaliteit + Geboorteland + Onderwijsniveau + Beroep + SBI + Burg..Staat, maxit = maxiterations, data = intcurrentdata)
    pp <- fitted(logresults)
  }
  # VARIABLE Beroep
  if (intvarnr == 10) {
    logresults <- multinom(Beroep ~ Geslacht + Leeftijd + HH_Pos + HH_grootte + Woonregio.vorig.jaar + Nationaliteit + Geboorteland + Onderwijsniveau + Econ..status + SBI + Burg..Staat, maxit = maxiterations, data = intcurrentdata)
    pp <- fitted(logresults)
  }
  # VARIABLE SBI
  if (intvarnr == 11) {
    logresults <- multinom(SBI ~ Geslacht + Leeftijd + HH_Pos + HH_grootte + Woonregio.vorig.jaar + Nationaliteit + Geboorteland + Onderwijsniveau + Econ..status + Beroep + Burg..Staat, maxit = maxiterations, data = intcurrentdata)
    pp <- fitted(logresults)
  }
  # VARIABLE Burg..Staat
  if (intvarnr == 12) {
    logresults <- multinom(Burg..Staat ~ Geslacht + Leeftijd + HH_Pos + HH_grootte + Woonregio.vorig.jaar + Nationaliteit + Geboorteland + Onderwijsniveau + Econ..status + Beroep + SBI, maxit = maxiterations, data = intcurrentdata)
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

AdjustProbs_StructZeroes_AllRecs <- function(intpp, intvarnr, intcurrentdata, intnum_missings, intmissing_indices) {
  intpp_small <- matrix(0,intnum_missings[intvarnr],ncats[intvarnr])
  # NO EDITS, SO NO STRUCTURAL ZEROES
  if ((intvarnr == 1) | (intvarnr %in% (5:7))) {
    for (ii in (1:intnum_missings[intvarnr])) {
      i <- intmissing_indices[intvarnr, ii]
      intpp_small[ii,] <- intpp[i,]
    }    
  }
  # EDITS, SO PERHAPS STRUCTURAL ZEROES
  if ((intvarnr %in% (2:4)) | (intvarnr %in% (8:12))) {
    for (ii in (1:intnum_missings[intvarnr])) {
      i <- intmissing_indices[intvarnr, ii]
      EditsInRec <- ExpEdits
      if (intvarnr != 2) {
        EditsInRec <- substValue(EditsInRec,"Leeftijd",intcurrentdata$Leeftijd[i], reduce = TRUE)
      }  
      if (intvarnr != 3) {
        EditsInRec <- substValue(EditsInRec,"HH_Pos",intcurrentdata$HH_Pos[i], reduce = TRUE)
      }  
      if (intvarnr != 4) {
        EditsInRec <- substValue(EditsInRec,"HH_grootte",intcurrentdata$HH_grootte[i], reduce = TRUE)
      }
      if (intvarnr != 8) {
        EditsInRec <- substValue(EditsInRec,"Onderwijsniveau",intcurrentdata$Onderwijsniveau[i], reduce = TRUE)
      }
      if (intvarnr != 9) {
        EditsInRec <- substValue(EditsInRec,"Econ..status",intcurrentdata$Econ..status[i], reduce = TRUE)
      }
      if (intvarnr != 10) {
        EditsInRec <- substValue(EditsInRec,"Beroep",intcurrentdata$Beroep[i], reduce = TRUE)
      }
      if (intvarnr != 11) {
        EditsInRec <- substValue(EditsInRec,"SBI",intcurrentdata$SBI[i], reduce = TRUE)
      }
      if (intvarnr != 12) {
        EditsInRec <- substValue(EditsInRec,"Burg..Staat",intcurrentdata$Burg..Staat[i], reduce = TRUE)
      }
      intpp_small[ii,] <- UseTestCatsToAdjustProbs(intvarnr, EditsInRec, intpp[intmissing_indices[intvarnr,ii],])
    }  
  }
  return(intpp_small)
}

AdjustProbs_StructZeroes_Rec <- function(intpp_rec, intvarnr, intcurrentdata, intmissing_indices) {
  intpp_small_rec <- rep(0,ncats[intvarnr])
  # NO EDITS, SO NO STRUCTURAL ZEROES
  if ((intvarnr == 1) | (intvarnr %in% (5:7))) {
    intpp_small_rec[ii,] <- intpp_rec[]
  }    
  # EDITS, SO PERHAPS STRUCTURAL ZEROES
  if ((intvarnr %in% (2:4)) | (intvarnr %in% (8:12))) {
    EditsInRec <- ExpEdits
    if (intvarnr != 2) {
      EditsInRec <- substValue(EditsInRec,"Leeftijd",intcurrentdata$Leeftijd[i], reduce = TRUE)
    }  
    if (intvarnr != 3) {
      EditsInRec <- substValue(EditsInRec,"HH_Pos",intcurrentdata$HH_Pos[i], reduce = TRUE)
    }  
    if (intvarnr != 4) {
      EditsInRec <- substValue(EditsInRec,"HH_grootte",intcurrentdata$HH_grootte[i], reduce = TRUE)
    }
    if (intvarnr != 8) {
      EditsInRec <- substValue(EditsInRec,"Onderwijsniveau",intcurrentdata$Onderwijsniveau[i], reduce = TRUE)
    }
    if (intvarnr != 9) {
      EditsInRec <- substValue(EditsInRec,"Econ..status",intcurrentdata$Econ..status[i], reduce = TRUE)
    }
    if (intvarnr != 10) {
      EditsInRec <- substValue(EditsInRec,"Beroep",intcurrentdata$Beroep[i], reduce = TRUE)
    }
    if (intvarnr != 11) {
      EditsInRec <- substValue(EditsInRec,"SBI",intcurrentdata$SBI[i], reduce = TRUE)
    }
    if (intvarnr != 12) {
      EditsInRec <- substValue(EditsInRec,"Burg..Staat",intcurrentdata$Burg..Staat[i], reduce = TRUE)
    }
    intpp_small_rec[] <- UseTestCatsToAdjustProbs(intvarnr, EditsInRec, intpp_rec[intmissing_indices[intvarnr,ii],])
  }
  return(intpp_small_rec)
}

########################
# FUNCTION: BINOMIAL APPROXIMATION FOR IMPUTATION PROBABILITIES
########################
Bin_Approximation <- function(int_pp, int_nummissings, intvarnr, int_ncats, intimptotals) {
  int_adjust_pp <- int_pp
  if (int_nummissings > 1) {
    for (jj in (1:int_ncats)) {
      Sum_column <- 0
      for (ii in (1:int_nummissings)) {
        Sum_column <- Sum_column + int_pp[ii,jj]
      }
      for (ii in (1:int_nummissings)) {
        if (intimptotals[intvarnr,jj] > 0.5) {
          BinomialProb <- (Sum_column - int_pp[ii,jj])/(int_nummissings - 1)
          if (BinomialProb < 1) {
            BinomialCorrection <- (int_nummissings - intimptotals[intvarnr,jj])/intimptotals[intvarnr,jj]
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
Pois_Approximation <- function(int_pp, int_nummissings, intvarnr, int_ncats, intimptotals) {
  int_adjust_pp <- int_pp
  for (jj in (1:int_ncats)) {
    Sum_column <- 0
    for (ii in (1:int_nummissings)) {
      Sum_column <- Sum_column + int_pp[ii,jj]
    }
    for (ii in (1:int_nummissings))  {
      if (intimptotals[intvarnr, jj] > 0.5) {
        PoissonCorrection <- (Sum_column - int_pp[ii,jj])/intimptotals[intvarnr, jj]
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
# FUNCTION IPF
########################
IPF <- function(int_pp, int_nummissings, intvarnr, int_ncats, intimptotals) {
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
    MaxDiff <- 0 # "TEST"
    for (jj in (1:int_ncats)) {
      Sum_column <- 0
      for (ii in (1:int_nummissings)) {
        Sum_column <- Sum_column + int_pp[ii,jj]
      }
      if (abs(Sum_column - intimptotals[intvarnr, jj]) > MaxDiff) {
        MaxDiff <- abs(Sum_column - intimptotals[intvarnr, jj])
      }
      if (Sum_column != 0) {
        AdjustmentFactor <- intimptotals[intvarnr,jj]/Sum_column
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
ApproximateProbs <- function(intpp_small, intvarnr, intmethod, intnum_missings, intimptotals) {
  if (intmethod == 1) {
    intpp_small <- Bin_Approximation(intpp_small, intnum_missings[intvarnr], intvarnr, ncats[intvarnr], intimptotals)  
  }
  if (intmethod == 2) {
    intpp_small <- Pois_Approximation(intpp_small, intnum_missings[intvarnr], intvarnr, ncats[intvarnr], intimptotals)
  }
  intpp_small <- IPF(intpp_small, intnum_missings[intvarnr], intvarnr, ncats[intvarnr], intimptotals)  
  return(intpp_small)
}


#############################################################
# CHECKING WHETHER DATA SATISFY EDITS
#############################################################
EditsCheck <- function(intcurrentdata, intnrows) {
  intsatisfied <- TRUE
  for (i in (1:intnrows)) {
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
TotalsCheck <- function(intcurrentdata, intnvars, intmax_ncats) {
  intsatisfied <- TRUE
  finaltotals <- matrix(NA,intnvars,intmax_ncats)
  for (j in (1:intnvars)) {
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
NAsCheck <- function(intcurrentdata, intnrows, intnvars) {
  intnoNAs <- TRUE
  for (i in (1:intnrows)) {
    for (j in (1:intnvars)) {
      if (is.na(intcurrentdata[i,j])) {
        intnoNAs <- FALSE
      }
    }
  }
  return(intnoNAs)
}


for (s in (1:NumMissingSets)) {  
  ###########################
  # SELECTING KIND OF MISSINGNESS MECHANISM: MCAR, MAR OR NMAR
  ###########################  
    originaldata <- Create_missings(Missfactor, truedata, nrows, nvars)
  #  filename_missingdata <- paste0(pathresults,"dataset_with_missings_", s + s_done, ".csv")
  #  write.csv2(originaldata, filename_missingdata)
  #  originaldata <- Create_MAR_Selective_Missings(Missfactor, varname)  
  #  originaldata <- Create_NMAR_Selective_Missings(Missfactor, varname)   
  
  # Maak B keer missings in pseudopopulation ("bootstrap samples")
  # imputeer bootstrap samples
  
  ###########################
  # INITIALIZING QUALITY MEASURES
  ###########################
  #  bias <- array(0, dim=c(3, nvars, max_ncats))
  #  variance <- array(0,dim=c(3, nvars, max_ncats))
  #  bias_8_10 <- array(0, dim=c(3, ncats[8], ncats[10]))
  #  var_8_10  <- array(0, dim=c(3, ncats[8], ncats[10]))
  #  bias_8_2  <- array(0, dim=c(3, ncats[8], ncats[2]))
  #  var_8_2  <- array(0, dim=c(3, ncats[8], ncats[2]))
  #  bias_8_9  <- array(0, dim=c(3, ncats[8], ncats[9]))
  #  var_8_9  <- array(0, dim=c(3, ncats[8], ncats[9]))
  #  bias_10_2  <- array(0, dim=c(3, ncats[10], ncats[2]))
  #  var_10_2  <- array(0, dim=c(3, ncats[10], ncats[2]))
  #  bias_10_9  <- array(0, dim=c(3, ncats[10], ncats[9]))
  #  var_10_9  <- array(0, dim=c(3, ncats[10], ncats[9]))
  #  CorrectlyImputed <- matrix(0, 3, nvars)
  #  FracCorrectlyImputed <- matrix(0, 3, nvars)
  
  ########################
  # DETERMINING MISSINGS PER RECORD (NECESSARY TO CREATE PSEUDO-POPULATIONS)
  ########################
  nummissingsperrec <- rep(0, nrows)
  for (ii in (1:nrows)) {
    for (jj in (1:nvars)) {
      if (is.na(originaldata[ii, jj + 1])) {
        nummissingsperrec[ii] <- nummissingsperrec[ii] + 1
      }
    }  
  }

  outfilename <- paste0(pathresults,"totals occupation pseudopop_", s + s_done, ".csv")
  logfilename <- paste0(path,"logfile.txt")
  
  pseudopopulation_created <- Create_PseudoPopulation(nummissingsperrec, originaldata, nrows)
  col_order <- c("nr", "Geslacht", "Leeftijd", "HH_Pos", "HH_grootte", "Woonregio.vorig.jaar", "Nationaliteit", "Geboorteland", "Onderwijsniveau", "Econ..status", "Beroep", "SBI", "Burg..Staat")
  pseudopopulation <- pseudopopulation_created[, col_order]
  rm(pseudopopulation_created)
  gc()
  
  filename_pseudopop <- paste0(pathresults,"pseudopopulation_", s + s_done, ".csv")
  write.csv2(pseudopopulation, filename_pseudopop)
  
  ########################
  # MARGINAL TOTALS THAT ARE ASSUMED TO BE KNOWN / IN OUR SIMULATION STUDY USUALLY ONLY TOTALS FOR VARIABLE 8 ("EDUCATIONAL LEVEL") ARE ASSUMED TO BE KNOWN
  ########################
  totals[8,] <- Determine_totals(8, pseudopopulation, max_ncats)
  
  clusterExport(myCluster, "s")
  clusterExport(myCluster, "s_done")
  clusterExport(myCluster, "Missfactor")
  clusterExport(myCluster, "MaxIter")
  clusterExport(myCluster, "IPFMaxCriterion")
  clusterExport(myCluster, "IterDiff")
  clusterExport(myCluster, "outfilename")
  clusterExport(myCluster, "logfilename")
  clusterExport(myCluster, "originaldata")
  clusterExport(myCluster, "nrows")
  clusterExport(myCluster, "nvars")
  clusterExport(myCluster, "ncats")
  clusterExport(myCluster, "ExpEdits")
  clusterExport(myCluster, "nummissingsperrec")
  clusterExport(myCluster, "max_ncats")
  clusterExport(myCluster, "totals")
  clusterExport(myCluster, "max_num_missings")
  clusterExport(myCluster, "fractions")
  clusterExport(myCluster, "codes")
  clusterExport(myCluster, "pseudopopulation")
    
  # clusterExport(myCluster, "Create_PseudoPopulation")
  clusterExport(myCluster, "Create_missings")
  #  clusterExport(myCluster, "Create_missings_OCC")
  clusterExport(myCluster, "Determine_totals")
  clusterExport(myCluster, "SimpleImpute")
  clusterExport(myCluster, "ImputeUntilSatisfied")
  clusterExport(myCluster, "ReduceEdits")
  clusterExport(myCluster, "Determine_missings")
  clusterExport(myCluster, "Determine_fractions_imptotals")  
  clusterExport(myCluster, "EstimateProbabilities") 
  clusterExport(myCluster, "AdjustProbs_StructZeroes_AllRecs") 
  clusterExport(myCluster, "AdjustProbs_StructZeroes_Rec") 
  clusterExport(myCluster, "UseTestCatsToAdjustProbs")
  clusterExport(myCluster, "ApproximateProbs")
#  clusterExport(myCluster, "Bin_Approximation")
#  clusterExport(myCluster, "Pois_Approximation")
  clusterExport(myCluster, "IPF")
  clusterExport(myCluster, "EditsCheck")
  clusterExport(myCluster, "TotalsCheck")
  clusterExport(myCluster, "NAsCheck")
  
  parLapply(myCluster, (1:numBootstrapSamples), function(x) {
    library(editrules)
    library(stats)
    library(nnet)
    library(foreign)
    library(Rcpp)
    B <- x
    
    ##########################################################################
    ##  Cox' ALGORITHM FOR CONTROLLED RANDOM ROUNDING (TO BASE 1) (developed by Sander Scholtus, 2018; update September 2019)
    ##########################################################################
    ##### R code for controlled rounding of a two-dimensional tabel using the method of Cox (1987) #####
    ## Sander Scholtus, 2018-2019
    
    # library(Rcpp)
    
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
        if (((iter == 1) | (iter %% 10 == 0)) & verbose) cat(sprintf('Iteration %d - Number of cells left to round off: %d\n', iter, remaining))
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

        
#    pseudopopulation_created <- Create_PseudoPopulation(nummissingsperrec, originaldata, nrows)
#    col_order <- c("nr", "Geslacht", "Leeftijd", "HH_Pos", "HH_grootte", "Woonregio.vorig.jaar", "Nationaliteit", "Geboorteland", "Onderwijsniveau", "Econ..status", "Beroep", "SBI", "Burg..Staat")
#    pseudopopulation <- pseudopopulation_created[, col_order]
#    rm(pseudopopulation_created)
#    gc()

    cat(paste0(" Starting missing set ", s + s_done, " ", B), file=logfilename, append = TRUE)

    ########################
    # MARGINAL TOTALS THAT ARE ASSUMED TO BE KNOWN / IN OUR SIMULATION STUDY USUALLY ONLY TOTALS FOR VARIABLE 8 ("EDUCATIONAL LEVEL") ARE ASSUMED TO BE KNOWN
    ########################
#    totals[8,] <- Determine_totals(8, pseudopopulation, max_ncats)

    ########################
    # INITIAL DATA: DATASET WIH MISSING DATA
    ########################
    initialdata <- Create_missings(Missfactor, pseudopopulation, nrows, nvars)

    missinglist <- Determine_missings(initialdata, nrows, nvars, max_num_missings)
    missings <- missinglist$missings
    num_missings <- missinglist$num_missings
    nummissingsrec <- missinglist$nummissingsrec
    missing_indices <- missinglist$missing_indices
    nummissingsrecinedits <- missinglist$nummissingsrecinedits
    fractions_imptotalslist <- Determine_fractions_imptotals(initialdata, totals, nrows, nvars, max_ncats)
    fractions <- fractions_imptotalslist$fractions
    imptotals <- fractions_imptotalslist$imptotals

  ##############################################################
  # START-UP PHASE
  ##############################################################
    cat("Start of start-up phase ", file=logfilename, append = TRUE)
    cat(paste0("Dataset with missings ", s + s_done, " ", imptotals[8,1]), file=logfilename, append = TRUE)

  ###############################
  # FOR EACH RECORD: SUBSTITURE OBSERVED VALUES INTO EXPLICIT EDITS AND ELIMINIATE VARIBALES WITH MISSING VALUES FROM SET OF EXPLICIT EDITS
  ###############################
    for (i in (1:nrows)) {
      if (nummissingsrec[i] > 0) {
        E0 <- ExpEdits
        E1 <- ReduceEdits(i, 1, E0, nummissingsrec[i], missings, initialdata)
        E2 <- ReduceEdits(i, 2, E1, nummissingsrec[i], missings, initialdata) 
        E3 <- ReduceEdits(i, 3, E2, nummissingsrec[i], missings, initialdata) 
        E4 <- ReduceEdits(i, 4, E3, nummissingsrec[i], missings, initialdata) 
        E5 <- ReduceEdits(i, 5, E4, nummissingsrec[i], missings, initialdata) 
        E6 <- ReduceEdits(i, 6, E5, nummissingsrec[i], missings, initialdata) 
        E7 <- ReduceEdits(i, 7, E6, nummissingsrec[i], missings, initialdata) 
        E8 <- ReduceEdits(i, 8, E7, nummissingsrec[i], missings, initialdata) 
        E9 <- ReduceEdits(i, 9, E8, nummissingsrec[i], missings, initialdata) 
        E10 <- ReduceEdits(i, 10, E9, nummissingsrec[i], missings, initialdata) 
        E11 <- ReduceEdits(i, 11, E10, nummissingsrec[i], missings, initialdata) 
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
          initialdata$Beroep[i] <- ImputeUntilSatisfied(10, E9)
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
        E0 <- substValue(E0,"Leeftijd", initialdata$Leeftijd[i], reduce = TRUE) 
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
    # ONLY FOR IPF
    for (approxmethod in (3:3)) {
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
          cat(paste0("Dataset with missings ", s + s_done, " ", B), file=logfilename, append = TRUE)
#          cat(paste0("Approximation method ", approxmethod), file=logfilename, append = TRUE)
          cat(paste0("Start iteration ", num_Iter), file=logfilename, append = TRUE)
          cat(paste0("Variable ", j), file=logfilename, append = TRUE)
          probs <- matrix(0, nrows, ncats[j])
          probs_small <- matrix(0, num_missings[j], ncats[j])
          ImpMatrix <- matrix(0, num_missings[j], ncats[j])
          probs <- EstimateProbabilities(j, currentdata)
          probs_small <- AdjustProbs_StructZeroes_AllRecs(probs, j, currentdata, num_missings, missing_indices)
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
            probs_approx <- ApproximateProbs(probs_small, j, approxmethod, num_missings, imptotals)
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
              ImpMatrix <- roundCox_cpp(A = probs_approx, verbose = FALSE)
#              ImpMatrix <- roundCox_cpp(A = probs_approx, cycle = 'first', verbose = TRUE)
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
                probs_small[ii,] <- AdjustProbs_StructZeroes_Rec(probs, j, currentdata, missing_indices)
              }
              currentdata[missing_indices[j, ii], j + 1] <- codes[j, SimpleImpute(j, probs_small[ii,])]
            }  
          }
          rm(ImpMatrix)
          rm(probs)
          rm(probs_small)
          gc()
        }
        EditsSatisfied <- EditsCheck(currentdata, nrows)
        if (!EditsSatisfied) {
          cat("Edits failed", file=logfilename, append = TRUE)
        }
        TotalsSatisfied <- TotalsCheck(currentdata, nvars, max_ncats)
        if (!TotalsSatisfied) {
          cat("Totals failed", file=logfilename, append = TRUE)
        }
        NAsSatisfied <- NAsCheck(currentdata, nrows, nvars)
        if (!NAsSatisfied) {
          cat("NAs in data", file=logfilename, append = TRUE)
        }
      }

      # CALCULATE TOTALS OF OCCUPATION 
      totals_occ <- rep(0,10)
      totals_occ <- table(currentdata[,11])
      # WRITE ESTIMATE FOR TOTAL FOR OCCUPATION TO FILE
      write.table(t(totals_occ[]), file = outfilename, append = TRUE, row.names = FALSE, col.names = FALSE, sep = ";")
    }  
  })
}  
stopCluster(myCluster)
Sys.time() - init
