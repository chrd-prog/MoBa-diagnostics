# Parse NPR data to generate some usable file format
# Start timer
library(haven)
library(data.table)
startTime <- proc.time()

# File locations ----------------------------------------------------------
# Where output files should be written
dirOutput   <- "/ess/p697/cluster/users/parekh/2025-05-21_NPR"
setwd(dirOutput)

# Latest SV_INFO file
dirSVINFO   <- "/tsd/p697/data/durable/phenotypes/sv_infofiles"
fileCurrent <- list.files(path = dirSVINFO, pattern = glob2rx("PDB2445_SV_INFO*.sav"))
if (length(fileCurrent) > 1)
{
  stop('More than one current version of INFO file found; script is not designed for this')
}
SVINFOFile  <- fileCurrent[1]

# Merged and linked NPR
dirNPR      <- "/ess/p697/cluster/users/parekh/2025-05-21_NPR"
NPRFile     <- paste(as.character(Sys.Date()), "-MoBa-LinkedNPR-AllSubjects.csv", sep = "")

# Withdrawn consent
dirConsent  <- "/ess/p697/cluster/users/parekh/2025-05-21_NPR"
consentFile <- paste(as.character(Sys.Date()), "-MoBa-WithdrawnList.csv", sep = "")

# Invalid NFR
fileInvalid <- list.files(path = dirSVINFO, pattern = glob2rx("*invalid*.csv"))

# Read files --------------------------------------------------------------
# dataNPR     <- read.csv(file.path(dirNPR,      NPRFile),     sep = "\t")
dataNPR     <- fread(file = file.path(dirNPR,   NPRFile),   sep = "\t")
withdrawn   <- read.csv(file.path(dirConsent,  consentFile), sep = "\t")
dataInfo    <- read_sav(file.path(dirSVINFO,   SVINFOFile))
dataInvalid <- read.csv(file.path(dirSVINFO,  fileInvalid), sep = "\t") 

# ID2445 for children -----------------------------------------------------
dataNPR$ID_2445   <- paste(dataNPR$PREG_ID_2445,  dataNPR$ChildNumber, sep = "_")

# Replace all NA_NA with blanks
dataNPR$ID_2445[dataNPR$ID_2445 == "NA_NA"] <- ""


# Eliminate withdrawn consent ---------------------------------------------
# Make a list of parents which are present in INFO file and withdrawn list:
# these are borderline cases and for the moment we are keeping them
borderline_M_ID <- withdrawn$M_ID_2445[withdrawn$M_ID_2445 %in% dataInfo$M_ID_2445 & withdrawn$M_ID_2445 != ""]
borderline_F_ID <- withdrawn$F_ID_2445[withdrawn$F_ID_2445 %in% dataInfo$F_ID_2445 & withdrawn$F_ID_2445 != ""]

# Any pregnancy ID which is present in the withdrawn list should be removed from NPR
locEliminate <- (dataNPR$PREG_ID_2445 %in% withdrawn$PREG_ID_2445) & 
                 dataNPR$PREG_ID_2445 != "" & !is.na(dataNPR$PREG_ID_2445)
dataNPR      <- dataNPR[!locEliminate, ]

# Any F_ID or M_ID which is present in the withdrawn list should be removed from
# NPR, except for the boderline cases
locEliminate <- dataNPR$F_ID_2445 %in% withdrawn$F_ID_2445 & dataNPR$F_ID_2445 != "" & !(dataNPR$F_ID_2445 %in% borderline_F_ID)
dataNPR      <- dataNPR[!locEliminate, ]

locEliminate <- dataNPR$M_ID_2445 %in% withdrawn$M_ID_2445 & dataNPR$M_ID_2445 != "" & !(dataNPR$M_ID_2445 %in% borderline_M_ID)
dataNPR      <- dataNPR[!locEliminate, ]


# Invalid FNR should be removed - they don't seem to exist in latest NPR (2025-05)?
# Only parents because we know that invalid FNR only contains fathers at the moment
# Might need updating in the future? Hopefully not
locEliminate <- dataNPR$F_ID_2445 %in% dataInvalid$ID_2445 |
                dataNPR$M_ID_2445 %in% dataInvalid$ID_2445
dataNPR      <- dataNPR[!locEliminate, ]


# Clean up ICD10 codes ----------------------------------------------------
# Remove any periods (.) in ICD10 codes
dataNPR$Diagnosis_ICDCode  <- sub("\\.", "", dataNPR$Diagnosis_ICDCode)

# Remove any question marks (?) in ICD10 codes
dataNPR$Diagnosis_ICDCode  <- sub("\\?", "", dataNPR$Diagnosis_ICDCode)

# Convert all codes to upper case
dataNPR$Diagnosis_ICDCode <- toupper(dataNPR$Diagnosis_ICDCode)

# Convert all NCMP NCSP codes to upper case - not used
# Currently, generates an error because the entry at 8793620 is "S\x8fR" 
# (actually S<8f>R); manually correcting it
dataNPR$Procedures_NCMPNCSPCode[dataNPR$Procedures_NCMPNCSPCode == "S\x8fR"] <- "Unknown"
dataNPR$Procedures_NCMPNCSPCode <- toupper(dataNPR$Procedures_NCMPNCSPCode)

# Save a version of data at this point - before counting
outputName <- paste(as.character(Sys.Date()), "-MoBa-ParsedNPR.csv", sep = "")
fwrite(dataNPR, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)
# write.table(dataNPR, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)


# Variables to ignore -----------------------------------------------------
varsToIgnore <- c("SerialNumber",                    "Key", 
                  "SourceTableName",                 "DiffDays_Admission",
                  "DiffDays_Discharge",              "Duration",
                  "Somatic_LevelofCare",             "Somatic_LevelofCare_2008_2009",
                  "Somatic_LevelofCare_from2010",    "Somatic_TreatmentClinic_Derived",
                  "Somatic_TreatmentClinic_Reported","ReferralType",                    
                  "ReferralDecision",                "InstitutionID",
                  "InstitutionID_Psuedoanonymized",  "InstitutionID_MentalHealth",
                  "ManualName",                      "ManualType",
                  "Diagnosis_MainSecondary",         "Diagnosis_Axis", 
                  "ChildCodeSequence",               "Procedures_CodeNumber",           
                  "Procedures_NCMPNCSPCode",         "REH_partOf_SOM", 
                  "Unfinished_WardStay",             "Original_kodenr", 
                  "Original_kodeverdi")

locsIgnore <- colnames(dataNPR) %in% varsToIgnore
locsIgnore <- !locsIgnore
dataNPR    <- dataNPR[, ..locsIgnore] # data.table nonsense
# dataNPR     <- dataNPR[, !locsIgnore]


# Remove any row without ICD10 code ---------------------------------------
locEmpty    <- is.na(dataNPR$Diagnosis_ICDCode) | dataNPR$Diagnosis_ICDCode == ""
dataNPR     <- dataNPR[!locEmpty,]


# Flatten table -----------------------------------------------------------
dupRows     <- duplicated(dataNPR)
dataNPR     <- dataNPR[!dupRows,]
# dataNPR   <- unique(dataNPR) # alternate solution


# Make one overall ID column ----------------------------------------------
dataNPR$ID  <- paste(dataNPR$ID_2445,  dataNPR$M_ID_2445,  dataNPR$F_ID_2445,  sep = "")


# Add Role column ---------------------------------------------------------
locChildren               <- dataNPR$ID_2445    != ""
locMothers                <- dataNPR$M_ID_2445  != ""
locFathers                <- dataNPR$F_ID_2445  != ""
dataNPR$Role[locChildren] <- "Child"
dataNPR$Role[locMothers]  <- "Mother"
dataNPR$Role[locFathers]  <- "Father"


# Eliminate other ID columns ----------------------------------------------
toEliminate <- !colnames(dataNPR) %in% c("ID_2445", "M_ID_2445", "F_ID_2445")
dataNPR     <- dataNPR[, ..toEliminate]
# dataNPR <- dataNPR[, !(colnames(dataNPR) %in% c("ID_2445", "M_ID_2445", "F_ID_2445"))]


# Rename ID to ID_2445 ----------------------------------------------------
colnames(dataNPR)[colnames(dataNPR) == "ID"] <- "ID_2445"


# Unique ICD10 codes and time points -------------------------------------
uqICDCodes    <- sort(unique(dataNPR$Diagnosis_ICDCode))
uqTimePoints  <- sort(unique(dataNPR$Year))


# Last year for which we have data? ---------------------------------------
lastYear <- as.numeric(uqTimePoints[length(uqTimePoints)])

# Combine ICD10 and time points -------------------------------------------
allUqComb   <- ""
for (time in uqTimePoints)
{
  allUqComb <- c(allUqComb, paste(uqICDCodes, time, sep = "_"))
}
allUqComb   <- allUqComb[allUqComb != ""]
allUqComb   <- sort(allUqComb)


# Reorganize NPR ----------------------------------------------------------
dataNPR     <- dataNPR[, c("ID_2445", "Role","PREG_ID_2445", "ChildNumber", 
                           "WithdrawnConsent18years", "Year", "Diagnosis_ICDCode")]


# Generate a temporary table having only unique information from NPR ------
temp        <- dataNPR[, c("ID_2445", "Role","PREG_ID_2445", "ChildNumber", "WithdrawnConsent18years")]
temp        <- unique(temp)
temp        <- temp[order(temp$ID_2445), ]


# Compile information -----------------------------------------------------
# Initialize
allColNames             <- c("ID_2445", "Role", "PREG_ID_2445", "ChildNumber", 
                             "WithdrawnConsent18years", paste(uqICDCodes, "_lifetime", sep = ""), 
                             allUqComb, gsub("_", "_x", allUqComb))
MinorDiagInfo           <- data.frame(matrix(nrow = nrow(temp), ncol = length(allColNames), data = 0))
colnames(MinorDiagInfo) <- allColNames

# Assign information from temp
MinorDiagInfo[, 1:5]    <- temp

startTime <- proc.time()
for (diag in uqICDCodes)
{
  for (time in uqTimePoints)
  {
    # Find IDs which have this diagnosis at this time point (cross sectional)
    whichIDs  <- dataNPR$ID_2445[dataNPR$Diagnosis_ICDCode == diag & dataNPR$Year == time]
    
    # Columns to mark
    whichCols <- paste(diag, time, sep = "_x")
    
    # Where are these IDs - logical
    locIDs    <- MinorDiagInfo$ID_2445 %in% whichIDs
    
    # Mark
    MinorDiagInfo[locIDs, whichCols] <- 1
    
    # Find IDs which have this diagnosis up to this time point (cumulative)
    whichIDs  <- dataNPR$ID_2445[dataNPR$Diagnosis_ICDCode == diag & dataNPR$Year == time]
    
    # Columns to mark
    whichCols <- paste(diag, time:lastYear, sep = "_")
    
    # Mark
    MinorDiagInfo[locIDs, whichCols] <- 1
  }
  # Info for lifetime - sum across all time points
  # Which columns to sum across=
  whichCols       <- c(paste(diag, uqTimePoints, sep = "_x"), paste(diag, uqTimePoints, sep = "_"))
  
  # Where to mark?
  toMark          <- paste(diag, "_lifetime", sep = "")
  
  # Mark
  MinorDiagInfo[, toMark] <- as.integer(as.logical(rowSums(MinorDiagInfo[, whichCols])))
}


# Convert minor codes to major codes --------------------------------------
allMajorCodes     <- unique(substring(uqICDCodes, 1, 3))
allMajorUqComb    <- ""
for (time in uqTimePoints)
{
  allMajorUqComb  <- c(allMajorUqComb, paste(allMajorCodes, time, sep = "_"))
}
allMajorUqComb    <- allMajorUqComb[allMajorUqComb != ""]
allMajorUqComb    <- sort(allMajorUqComb)


# Compile major codes -----------------------------------------------------
# Initialize results table
allMajColNames                <- c("ID_2445", "Role","PREG_ID_2445", "ChildNumber", 
                                   "WithdrawnConsent18years", paste(allMajorCodes, "_lifetime", sep = ""),
                                   allMajorUqComb, gsub("_", "_x", allMajorUqComb))
MajorDiagInfo                 <- data.frame(matrix(nrow = nrow(temp), ncol = length(allMajColNames), data = 0))
colnames(MajorDiagInfo)       <- allMajColNames

# Copy information over from minor codes table
MajorDiagInfo[, 1:5]          <- MinorDiagInfo[, 1:5]

# Counting here is simpler - regular expression for every diagnostic category; pool
tmpColNames                   <- colnames(MinorDiagInfo)
count                         <- 1
matchedColumns                <- matrix("", nrow = length(allMajorUqComb), ncol = 1)
for (cols in allMajorUqComb)
{
  matchCols                   <- grep(glob2rx(paste(substr(cols, 1, 3), "*", substr(cols, nchar(cols)-4, nchar(cols)), sep = "")), tmpColNames)
  if (length(matchCols) > 1)
  {
    MajorDiagInfo[,cols]      <- as.integer(as.logical(rowSums(MinorDiagInfo[,matchCols])))
  }
  else
  {
    MajorDiagInfo[,cols]      <- as.integer(as.logical(MinorDiagInfo[,matchCols]))
  }
  
  if (length(matchCols) > 1)
  {
    matchedColumns[count]     <- toString(paste(tmpColNames[matchCols], "", sep = ""))
  }
  else
  {
    matchedColumns[count]     <- toString(tmpColNames[matchCols])
  }
  count <- count + 1
}


# Repeat for cross-sectional ----------------------------------------------
namesCross                    <- gsub("_", "_x", allMajorUqComb)
count                         <- 1
matchedxColumns               <- ""
for (cols in namesCross)
{
  matchCols                   <- grep(glob2rx(paste(substr(cols, 1, 3), "*", substr(cols, nchar(cols)-5, nchar(cols)), sep = "")), tmpColNames)
  if (length(matchCols) > 1)
  {
    MajorDiagInfo[,cols]      <- as.integer(as.logical(rowSums(MinorDiagInfo[,matchCols])))
  }
  else
  {
    MajorDiagInfo[,cols]      <- as.integer(as.logical(MinorDiagInfo[,matchCols]))
  }
  if (length(matchCols) > 1)
  {
    matchedxColumns[count]    <- toString(paste(tmpColNames[matchCols], "", sep = ""))
  }
  else
  {
    matchedxColumns[count]    <- toString(tmpColNames[matchCols])
  }
  count <- count + 1
}


# Finally for lifetime ---------------------------------------------------
namesLifetime                 <- paste(allMajorCodes, "_lifetime", sep = "")
count                         <- 1
matchedLColumns               <- ""
for (cols in namesLifetime)
{
  matchCols                   <- grep(glob2rx(paste(substr(cols, 1, 3), "*_lifetime", sep = "")), tmpColNames)
  if (length(matchCols) > 1)
  {
    MajorDiagInfo[,cols]      <- as.integer(as.logical(rowSums(MinorDiagInfo[,matchCols])))
  }
  else
  {
    MajorDiagInfo[,cols]      <- as.integer(as.logical(MinorDiagInfo[,matchCols]))
  }
  if (length(matchCols) > 1)
  {
    matchedLColumns[count]    <- toString(paste(tmpColNames[matchCols], "", sep = ""))
  }
  else
  {
    matchedLColumns[count]    <- toString(tmpColNames[matchCols])
  }
  count <- count + 1
}


# Save results ------------------------------------------------------------
# Minor codes
outputName <- paste(as.character(Sys.Date()), "-MoBa-ParsedNPR_allCodes.csv", sep = "")
fwrite(MinorDiagInfo, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)
# write.table(MinorDiagInfo, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)

# Major codes
outputName <- paste(as.character(Sys.Date()), "-MoBa-ParsedNPR_majorCodes.csv", sep = "")
fwrite(MajorDiagInfo, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)
# write.table(MajorDiagInfo, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)


# Split info by year ------------------------------------------------------
dirOutput     <- file.path(dirOutput, "splitInfo")
if (!dir.exists(dirOutput))
{
  dir.create(dirOutput)
}

# All (minor) codes
for (times in uqTimePoints)
{
  toKeepCols  <- c("ID_2445", "Role","PREG_ID_2445", "ChildNumber", "WithdrawnConsent18years", 
                   colnames(MinorDiagInfo[grep(glob2rx(paste("*_*", times, sep = "")), 
                                               colnames(MinorDiagInfo))]))
  temp        <- MinorDiagInfo[, toKeepCols]
  
  outputName  <- paste(as.character(Sys.Date()), "-MoBa-ParsedNPR_allCodes-", times, ".csv", sep = "")
  fwrite(temp, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)
  # write.table(temp, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)
  
}

# Minor codes - lifetime
toKeepCols  <- c("ID_2445", "Role","PREG_ID_2445", "ChildNumber", "WithdrawnConsent18years", 
                 colnames(MinorDiagInfo[grep(glob2rx("*_*lifetime"), 
                                             colnames(MinorDiagInfo))]))
temp        <- MinorDiagInfo[, toKeepCols]

outputName  <- paste(as.character(Sys.Date()), "-MoBa-ParsedNPR_allCodes-lifetime.csv", sep = "")
fwrite(temp, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)
# write.table(temp, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)

# All major codes
for (times in uqTimePoints)
{
  toKeepCols  <- c("ID_2445", "Role","PREG_ID_2445", "ChildNumber", "WithdrawnConsent18years", 
                   colnames(MajorDiagInfo[grep(glob2rx(paste("*_*", times, sep = "")), 
                                               colnames(MajorDiagInfo))]))
  temp        <- MajorDiagInfo[, toKeepCols]
  
  outputName  <- paste(as.character(Sys.Date()), "-MoBa-ParsedNPR_majorCodes-", times, ".csv", sep = "")
  fwrite(temp, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)
  # write.table(temp, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)
  
}

# Minor codes - lifetime
toKeepCols  <- c("ID_2445", "Role","PREG_ID_2445", "ChildNumber", "WithdrawnConsent18years", 
                 colnames(MajorDiagInfo[grep(glob2rx("*_*lifetime"), 
                                             colnames(MajorDiagInfo))]))
temp        <- MajorDiagInfo[, toKeepCols]

outputName  <- paste(as.character(Sys.Date()), "-MoBa-ParsedNPR_majorCodes-lifetime.csv", sep = "")
fwrite(temp, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)
# write.table(temp, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)