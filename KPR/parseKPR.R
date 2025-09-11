# Parse KPR
library(haven)
library(data.table)

# File locations ----------------------------------------------------------
# Where output files should be written
dirWorking <- "/ess/p697/cluster/users/parekh/2025-05-22_KPR"
setwd(dirWorking)

# Latest SV_INFO file
dirSVINFO   <- "/tsd/p697/data/durable/phenotypes/sv_infofiles"
fileCurrent <- list.files(path = dirSVINFO, pattern = glob2rx("PDB2445_SV_INFO*.sav"))
if (length(fileCurrent) > 1)
{
  stop('More than one current version of INFO file found; script is not designed for this')
}
SVINFOFile  <- fileCurrent[1]

# Merged and linked KPR
dirKPR      <- "/ess/p697/cluster/users/parekh/2025-05-22_KPR"
KPRFile     <- paste(as.character(Sys.Date()), "-MoBa-LinkedKPR-AllSubjects.csv", sep = "")

# Withdrawn consent
dirConsent  <- "/ess/p697/cluster/users/parekh/2025-05-22_KPR"
consentFile <- paste(as.character(Sys.Date()), "-MoBa-WithdrawnList.csv", sep = "")

# Invalid NFR
fileInvalid <- list.files(path = dirSVINFO, pattern = glob2rx("*invalid*.csv"))

# Read files --------------------------------------------------------------
# dataKPR     <- read.csv(file.path(dirKPR,      KPRFile),     sep = "\t")
dataKPR     <- fread(file = file.path(dirKPR,   KPRFile),   sep = "\t")
withdrawn   <- read.csv(file.path(dirConsent,  consentFile), sep = "\t")
dataInfo    <- read_sav(file.path(dirSVINFO,   SVINFOFile))
dataInvalid <- read.csv(file.path(dirSVINFO,  fileInvalid), sep = "\t") 

# ID2445 for children -----------------------------------------------------
dataKPR$ID_2445   <- paste(dataKPR$PREG_ID_2445,  dataKPR$ChildNumber, sep = "_")

# Replace all NA_NA with blanks
dataKPR$ID_2445[dataKPR$ID_2445 == "NA_NA"] <- ""


# Eliminate withdrawn consent ---------------------------------------------
# Make a list of parents which are present in INFO file and withdrawn list:
# these are borderline cases and for the moment we are keeping them
borderline_M_ID <- withdrawn$M_ID_2445[withdrawn$M_ID_2445 %in% dataInfo$M_ID_2445 & withdrawn$M_ID_2445 != ""]
borderline_F_ID <- withdrawn$F_ID_2445[withdrawn$F_ID_2445 %in% dataInfo$F_ID_2445 & withdrawn$F_ID_2445 != ""]

# Any pregnancy ID which is present in the withdrawn list should be removed from KPR
locEliminate <- (dataKPR$PREG_ID_2445 %in% withdrawn$PREG_ID_2445) & 
                 dataKPR$PREG_ID_2445 != "" & !is.na(dataKPR$PREG_ID_2445)
dataKPR      <- dataKPR[!locEliminate, ]

# Any F_ID or M_ID which is present in the withdrawn list should be removed from
# KPR, except for the boderline cases
locEliminate <- dataKPR$F_ID_2445 %in% withdrawn$F_ID_2445 & dataKPR$F_ID_2445 != "" & !(dataKPR$F_ID_2445 %in% borderline_F_ID)
dataKPR      <- dataKPR[!locEliminate, ]

locEliminate <- dataKPR$M_ID_2445 %in% withdrawn$M_ID_2445 & dataKPR$M_ID_2445 != "" & !(dataKPR$M_ID_2445 %in% borderline_M_ID)
dataKPR      <- dataKPR[!locEliminate, ]


# Invalid FNR should be removed
# Only parents because we know that invalid FNR only contains fathers at the moment
# Might need updating in the future? Hopefully not
locEliminate <- dataKPR$F_ID_2445 %in% dataInvalid$ID_2445 |
                dataKPR$M_ID_2445 %in% dataInvalid$ID_2445
dataKPR      <- dataKPR[!locEliminate, ]


# All codes seem to be of 3 character length
# Remove any periods (.)
dataKPR$DiagnosisCode  <- sub("\\.", "", dataKPR$DiagnosisCode)

# Convert all codes to upper case
dataKPR$DiagnosisCode <- toupper(dataKPR$DiagnosisCode)


# Split by diagnostic manual ----------------------------------------------
# Separate out ICPC codes
loc_ICPC     <- dataKPR$DiagnosisType == "ICPC-2" | dataKPR$DiagnosisType == "ICPC-2B"
dataKPR_ICPC <- dataKPR[loc_ICPC,]

# Leftover is ICD-10
dataKPR_ICD <- dataKPR[!loc_ICPC,]

stopifnot(nrow(dataKPR_ICPC) + nrow(dataKPR_ICD) - nrow(dataKPR) == 0)


## Save versions of data at this point - before counting
outputName <- paste(as.character(Sys.Date()), "-MoBa-ParsedKPR-All.csv", sep = "")
fwrite(dataKPR, file = file.path(dirWorking, outputName), sep = "\t", row.names = F)

outputName <- paste(as.character(Sys.Date()), "-MoBa-ParsedKPR-Manual-ICPC.csv", sep = "")
fwrite(dataKPR_ICPC, file = file.path(dirWorking, outputName), sep = "\t", row.names = F)

outputName <- paste(as.character(Sys.Date()), "-MoBa-ParsedKPR-Manual-ICD10.csv", sep = "")
fwrite(dataKPR_ICD, file = file.path(dirWorking, outputName), sep = "\t", row.names = F)


# Variables to ignore -----------------------------------------------------
varsToIgnore <- c("SerialNumber",              "Key", 
                  "DiagnosisType",             "isMainDiagnosis",
                  "DiffDays_Admission",        "MunicipalityNumber",
                  "TreatmentMunicipality",     "ServiceType",
                  "CountyNumber",              "TreatmentCounty")

locsIgnore <- colnames(dataKPR) %in% varsToIgnore
locsIgnore <- !locsIgnore
dataKPR    <- dataKPR[, ..locsIgnore] # data.table nonsense

locsIgnore   <- colnames(dataKPR_ICPC) %in% varsToIgnore
locsIgnore   <- !locsIgnore
dataKPR_ICPC <- dataKPR_ICPC[, ..locsIgnore] # data.table nonsense

locsIgnore   <- colnames(dataKPR_ICD) %in% varsToIgnore
locsIgnore   <- !locsIgnore
dataKPR_ICD  <- dataKPR_ICD[, ..locsIgnore] # data.table nonsense


# Get rid of KPR full table
rm(dataKPR)

# Flatten table -----------------------------------------------------------
dupRows      <- duplicated(dataKPR_ICPC)
dataKPR_ICPC <- dataKPR_ICPC[!dupRows,]

dupRows      <- duplicated(dataKPR_ICD)
dataKPR_ICD  <- dataKPR_ICD[!dupRows,]


# Make one overall ID column ----------------------------------------------
dataKPR_ICPC$ID  <- paste(dataKPR_ICPC$ID_2445,  dataKPR_ICPC$M_ID_2445,  dataKPR_ICPC$F_ID_2445,  sep = "")
dataKPR_ICD$ID   <- paste(dataKPR_ICD$ID_2445,   dataKPR_ICD$M_ID_2445,   dataKPR_ICD$F_ID_2445,   sep = "")


# Add Role column ---------------------------------------------------------
locChildren               <- dataKPR_ICPC$ID_2445    != ""
locMothers                <- dataKPR_ICPC$M_ID_2445  != ""
locFathers                <- dataKPR_ICPC$F_ID_2445  != ""
dataKPR_ICPC$Role[locChildren] <- "Child"
dataKPR_ICPC$Role[locMothers]  <- "Mother"
dataKPR_ICPC$Role[locFathers]  <- "Father"

# Add Role column ---------------------------------------------------------
locChildren               <- dataKPR_ICD$ID_2445    != ""
locMothers                <- dataKPR_ICD$M_ID_2445  != ""
locFathers                <- dataKPR_ICD$F_ID_2445  != ""
dataKPR_ICD$Role[locChildren] <- "Child"
dataKPR_ICD$Role[locMothers]  <- "Mother"
dataKPR_ICD$Role[locFathers]  <- "Father"


# Eliminate other ID columns ----------------------------------------------
toEliminate  <- !colnames(dataKPR_ICPC) %in% c("ID_2445", "M_ID_2445", "F_ID_2445")
dataKPR_ICPC <- dataKPR_ICPC[, ..toEliminate]

toEliminate <- !colnames(dataKPR_ICD) %in% c("ID_2445", "M_ID_2445", "F_ID_2445")
dataKPR_ICD <- dataKPR_ICD[, ..toEliminate]


# Rename ID to ID_2445 ----------------------------------------------------
colnames(dataKPR_ICPC)[colnames(dataKPR_ICPC) == "ID"] <- "ID_2445"
colnames(dataKPR_ICD)[colnames(dataKPR_ICD) == "ID"]   <- "ID_2445"


# Unique ICD10 codes and time points -------------------------------------
uqICDCodes    <- sort(unique(dataKPR_ICD$DiagnosisCode))
uqTimePoints  <- sort(unique(dataKPR_ICD$Year))

# Last year for which we have data? ---------------------------------------
lastYear <- as.numeric(uqTimePoints[length(uqTimePoints)])

# Combine ICD10 and time points -------------------------------------------
allUqComb_ICD   <- ""
for (time in uqTimePoints)
{
  allUqComb_ICD <- c(allUqComb_ICD, paste(uqICDCodes, time, sep = "_"))
}
allUqComb_ICD   <- allUqComb_ICD[allUqComb_ICD != ""]
allUqComb_ICD   <- sort(allUqComb_ICD)


# Reorganize ICD ----------------------------------------------------------
dataKPR_ICD     <- dataKPR_ICD[, c("ID_2445", "Role","PREG_ID_2445", "ChildNumber", 
                           "WithdrawnConsent18years", "Year", "DiagnosisCode")]

# Generate a temporary table having only unique information from KPR ICD ------
temp        <- dataKPR_ICD[, c("ID_2445", "Role","PREG_ID_2445", "ChildNumber", "WithdrawnConsent18years")]
temp        <- unique(temp)
temp        <- temp[order(temp$ID_2445), ]


# Compile information -----------------------------------------------------
# Initialize
allColNames_ICD       <- c("ID_2445", "Role", "PREG_ID_2445", "ChildNumber", 
                           "WithdrawnConsent18years", paste(uqICDCodes, "_lifetime", sep = ""), 
                           allUqComb_ICD, gsub("_", "_x", allUqComb_ICD))
DiagInfo_ICD           <- data.frame(matrix(nrow = nrow(temp), ncol = length(allColNames_ICD), data = 0))
colnames(DiagInfo_ICD) <- allColNames_ICD

# Assign information from temp
DiagInfo_ICD[, 1:5]    <- temp

for (diag in uqICDCodes)
{
  for (time in uqTimePoints)
  {
    # Find IDs which have this diagnosis at this time point (cross sectional)
    whichIDs  <- dataKPR_ICD$ID_2445[dataKPR_ICD$DiagnosisCode == diag & dataKPR_ICD$Year == time]
    
    # Columns to mark
    whichCols <- paste(diag, time, sep = "_x")
    
    # Where are these IDs - logical
    locIDs    <- DiagInfo_ICD$ID_2445 %in% whichIDs
    
    # Mark
    DiagInfo_ICD[locIDs, whichCols] <- 1
    
    # Find IDs which have this diagnosis up to this time point (cumulative)
    whichIDs  <- dataKPR_ICD$ID_2445[dataKPR_ICD$DiagnosisCode == diag & dataKPR_ICD$Year == time]
    
    # Columns to mark
    whichCols <- paste(diag, time:lastYear, sep = "_")
    
    # Mark
    DiagInfo_ICD[locIDs, whichCols] <- 1
  }
  # Info for lifetime - sum across all time points
  # Which columns to sum across=
  whichCols       <- c(paste(diag, uqTimePoints, sep = "_x"), paste(diag, uqTimePoints, sep = "_"))
  
  # Where to mark?
  toMark          <- paste(diag, "_lifetime", sep = "")
  
  # Mark
  DiagInfo_ICD[, toMark] <- as.integer(as.logical(rowSums(DiagInfo_ICD[, whichCols])))
}

# Save results ------------------------------------------------------------
outputName <- paste(as.character(Sys.Date()), "-MoBa-ParsedKPR_ICD10Codes.csv", sep = "")
fwrite(DiagInfo_ICD, file = file.path(dirWorking, outputName), sep = "\t", row.names = F)


# Unique ICPC codes and time points -------------------------------------
uqICPCCodes    <- sort(unique(dataKPR_ICPC$DiagnosisCode))
uqTimePoints  <- sort(unique(dataKPR_ICPC$Year))

# Last year for which we have data? ---------------------------------------
lastYear <- as.numeric(uqTimePoints[length(uqTimePoints)])

# Combine ICPC10 and time points -------------------------------------------
allUqComb_ICPC   <- ""
for (time in uqTimePoints)
{
  allUqComb_ICPC <- c(allUqComb_ICPC, paste(uqICPCCodes, time, sep = "_"))
}
allUqComb_ICPC   <- allUqComb_ICPC[allUqComb_ICPC != ""]
allUqComb_ICPC   <- sort(allUqComb_ICPC)


# Reorganize ICPC ----------------------------------------------------------
dataKPR_ICPC     <- dataKPR_ICPC[, c("ID_2445", "Role","PREG_ID_2445", "ChildNumber", 
                                     "WithdrawnConsent18years", "Year", "DiagnosisCode")]

# Generate a temporary table having only unique information from KPR ICPC ------
temp        <- dataKPR_ICPC[, c("ID_2445", "Role","PREG_ID_2445", "ChildNumber", "WithdrawnConsent18years")]
temp        <- unique(temp)
temp        <- temp[order(temp$ID_2445), ]


# Compile information -----------------------------------------------------
# Initialize
allColNames_ICPC       <- c("ID_2445", "Role", "PREG_ID_2445", "ChildNumber", 
                            "WithdrawnConsent18years", paste(uqICPCCodes, "_lifetime", sep = ""), 
                            allUqComb_ICPC, gsub("_", "_x", allUqComb_ICPC))
DiagInfo_ICPC           <- data.frame(matrix(nrow = nrow(temp), ncol = length(allColNames_ICPC), data = 0))
colnames(DiagInfo_ICPC) <- allColNames_ICPC

# Assign information from temp
DiagInfo_ICPC[, 1:5]    <- temp

for (diag in uqICPCCodes)
{
  for (time in uqTimePoints)
  {
    # Find IDs which have this diagnosis at this time point (cross sectional)
    whichIDs  <- dataKPR_ICPC$ID_2445[dataKPR_ICPC$DiagnosisCode == diag & dataKPR_ICPC$Year == time]
    
    # Columns to mark
    whichCols <- paste(diag, time, sep = "_x")
    
    # Where are these IDs - logical
    locIDs    <- DiagInfo_ICPC$ID_2445 %in% whichIDs
    
    # Mark
    DiagInfo_ICPC[locIDs, whichCols] <- 1
    
    # Find IDs which have this diagnosis up to this time point (cumulative)
    whichIDs  <- dataKPR_ICPC$ID_2445[dataKPR_ICPC$DiagnosisCode == diag & dataKPR_ICPC$Year == time]
    
    # Columns to mark
    whichCols <- paste(diag, time:lastYear, sep = "_")
    
    # Mark
    DiagInfo_ICPC[locIDs, whichCols] <- 1
  }
  # Info for lifetime - sum across all time points
  # Which columns to sum across=
  whichCols       <- c(paste(diag, uqTimePoints, sep = "_x"), paste(diag, uqTimePoints, sep = "_"))
  
  # Where to mark?
  toMark          <- paste(diag, "_lifetime", sep = "")
  
  # Mark
  DiagInfo_ICPC[, toMark] <- as.integer(as.logical(rowSums(DiagInfo_ICPC[, whichCols])))
}

# Save results ------------------------------------------------------------
outputName <- paste(as.character(Sys.Date()), "-MoBa-ParsedKPR_ICPCCodes.csv", sep = "")
fwrite(DiagInfo_ICPC, file = file.path(dirWorking, outputName), sep = "\t", row.names = F)