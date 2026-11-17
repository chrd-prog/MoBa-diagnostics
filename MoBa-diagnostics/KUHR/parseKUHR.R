# Parse KUHR data to generate some usable file format
# Start timer
library(haven)
library(data.table)

# Note to reader: the function uses a strange ..Variable to select columns, 
# instead of something standard like just using a logical variable - this is 
# because of the data.table library that this script uses. Pretty annoying - PP

# File locations ----------------------------------------------------------
# Where output files should be written
dirOutput   <- "/ess/p697/cluster/users/parekh/2023-10-24_KUHR"
setwd(dirOutput)

# Latest SV_INFO file
dirSVINFO   <- "/tsd/p697/data/durable/phenotypes/sv_infofiles"
fileCurrent <- list.files(path = dirSVINFO, pattern = glob2rx("PDB2445_SV_INFO*.sav"))
if (length(fileCurrent) > 1)
{
  stop('More than one current version of INFO file found; script is not designed for this')
}
SVINFOFile  <- fileCurrent[1]

# Merged and linked KUHR
dirKUHR      <- "/ess/p697/cluster/users/parekh/2023-10-24_KUHR"
KUHRFile     <- paste(as.character(Sys.Date()), "-MoBa-LinkedKUHR-AllSubjects.csv", sep = "")

# Withdrawn consent
dirConsent  <- "/ess/p697/cluster/users/parekh/2023-10-24_KUHR"
consentFile <- paste(as.character(Sys.Date()), "-MoBa-WithdrawnList.csv", sep = "")

# Read files --------------------------------------------------------------
# dataKUHR    <- read.csv(file.path(dirKUHR,     KUHRFile),    sep = ",")
dataKUHR    <- fread(file = file.path(dirKUHR, KUHRFile), sep = ",")
withdrawn   <- read.csv(file.path(dirConsent,  consentFile), sep = ",")
dataInfo    <- read_sav(file.path(dirSVINFO,   SVINFOFile))


# ID2445 for children -----------------------------------------------------
dataKUHR$ID_2445   <- paste(dataKUHR$PREG_ID_2445,  dataKUHR$ChildNumber, sep = "_")

# Replace all NA_NA with blanks
dataKUHR$ID_2445[dataKUHR$ID_2445 == "NA_NA"] <- ""


# Eliminate withdrawn consent ---------------------------------------------
# Make a list of parents which are present in INFO file and withdrawn list:
# these are borderline cases and for the moment we are keeping them
borderline_M_ID <- withdrawn$M_ID_2445[withdrawn$M_ID_2445 %in% dataInfo$M_ID_2445 & withdrawn$M_ID_2445 != ""]
borderline_F_ID <- withdrawn$F_ID_2445[withdrawn$F_ID_2445 %in% dataInfo$F_ID_2445 & withdrawn$F_ID_2445 != ""]

# Any pregnancy ID which is present in the withdrawn list should be removed from KUHR
locEliminate <- (dataKUHR$PREG_ID_2445 %in% withdrawn$PREG_ID_2445) & 
                 dataKUHR$PREG_ID_2445 != "" & !is.na(dataKUHR$PREG_ID_2445)
dataKUHR     <- dataKUHR[!locEliminate, ]

# Any F_ID or M_ID which is present in the withdrawn list should be removed from
# KUHR, except for the borderline cases
locEliminate <- dataKUHR$F_ID_2445 %in% withdrawn$F_ID_2445 & dataKUHR$F_ID_2445 != "" & !(dataKUHR$F_ID_2445 %in% borderline_F_ID)
dataKUHR     <- dataKUHR[!locEliminate, ]

locEliminate <- dataKUHR$M_ID_2445 %in% withdrawn$M_ID_2445 & dataKUHR$M_ID_2445 != "" & !(dataKUHR$M_ID_2445 %in% borderline_M_ID)
dataKUHR     <- dataKUHR[!locEliminate, ]


# Some clean up ------------------------------------------------------------
# Remove rows with no diagnostic information 
dataKUHR <- dataKUHR[dataKUHR$Diagnosis != "", ]

# Remove rows with just a comma in the diagnostic information
dataKUHR <- dataKUHR[dataKUHR$Diagnosis != ",", ]

# Remove rows with just zero in the diagnostic information
dataKUHR <- dataKUHR[dataKUHR$Diagnosis != "0", ]

# Remove duplicated rows
rowsDuplicated  <- duplicated(dataKUHR)
dataKUHR        <- dataKUHR[!rowsDuplicated,]

# There is at least one entry where the code ends in a comma - will create
# problems down the line - find and clean
lastChars <- substring(dataKUHR$Diagnosis, first = nchar(dataKUHR$Diagnosis), last = nchar(dataKUHR$Diagnosis))
locComma  <- lastChars == ","
dataKUHR  <- dataKUHR[!locComma, ]

# There is at least one entry with forward slash - convert to comma
dataKUHR$Diagnosis <- sub("/", ",", dataKUHR$Diagnosis)

# There is at least one entry with a hyphen or dash - remove the hyphen
dataKUHR$Diagnosis <- sub("-", "", dataKUHR$Diagnosis)

# All ICD-9 codes are actually ICD-10 codes
dataKUHR$DiagnosticManual[dataKUHR$DiagnosticManual == "ICD-9"] <- "ICD-10"

# Some ICD-10 codes have periods in them - replace commas
dataKUHR$Diagnosis <- gsub("\\.", "", dataKUHR$Diagnosis)


# Need to deal with possibly incorrect diagnostic manual:
# 2000, BUP_3, F840, 0050, 4000
# F840, BUP_3, 0050, 2000, 4000
# There are three entries where the diagnostic information mentions BUP_3 while
# the diagnostic manual is ICD10 - change the manual to Akser i BUP
tmpMisLabelManual                            <- grep(glob2rx("*BUP*"), dataKUHR$Diagnosis)
dataKUHR$DiagnosticManual[tmpMisLabelManual] <- "Akser i BUP-klassifikasjon"


# There are some additional information with some diagnostic codes - clean them
# These include the words "MARGIANL", "MARGINAL", "PLANTARFASCIALF", "MYALGI",
# "MYALGIISKULDERR", "STERKT", "LUMBAGO", "BEHANDLINGSOMOMF", "RETINERTE"
dataKUHR$Diagnosis <- gsub("MARGIANL",         "", dataKUHR$Diagnosis)
dataKUHR$Diagnosis <- gsub("MARGINAL",         "", dataKUHR$Diagnosis)
dataKUHR$Diagnosis <- gsub("PLANTARFASCIALF",  "", dataKUHR$Diagnosis)
dataKUHR$Diagnosis <- gsub("MYALGIISKULDERR",  "", dataKUHR$Diagnosis)
dataKUHR$Diagnosis <- gsub("MYALGI",           "", dataKUHR$Diagnosis)
dataKUHR$Diagnosis <- gsub("STERKT",           "", dataKUHR$Diagnosis)
dataKUHR$Diagnosis <- gsub("LUMBAGO",          "", dataKUHR$Diagnosis)
dataKUHR$Diagnosis <- gsub("BEHANDLINGSOMOMF", "", dataKUHR$Diagnosis)
dataKUHR$Diagnosis <- gsub("RETINERTE",        "", dataKUHR$Diagnosis)

# Remove rows with no diagnostic information 
dataKUHR <- dataKUHR[dataKUHR$Diagnosis != "", ]


# Sometimes multiple diagnostic information is present in a single row - these
# should be separated out as separate rows - probably very inefficient
# First, find the commas in the diagnosis column
# locsMultiple <- lengths(gregexec(pattern = ",", dataKUHR$Diagnosis))
locsMultiple <- grepl(pattern = ",", dataKUHR$Diagnosis)

# Separate out dataset that has only one diagnostic entry
dataKUHR_one <- dataKUHR[!locsMultiple, ]

# Remaining data
dataRemaining <- dataKUHR[locsMultiple, ]

# Make another column in remaining data that contains the true string
dataRemaining$ParsedDiagnosis <- dataRemaining$Diagnosis

# Same for data that has one diagnois
dataKUHR_one$ParsedDiagnosis <- dataKUHR_one$Diagnosis

# How many commas exist in the remaining dataset?
NumCommas    <- gregexec(pattern = ",", dataRemaining$Diagnosis)
LenNumCommas <- lengths(NumCommas)
tmpNumCommas <- sort(unique(LenNumCommas))

# Initialize a data frame where expanded information will be saved
dataExpanded <- data.frame(matrix(nrow = 0, ncol = length(colnames(dataRemaining))))
colnames(dataExpanded) <- colnames(dataRemaining)

# For every number of commas, go into the remaining dataset and expand
for(commas in tmpNumCommas)
{
  currDataset <- dataRemaining[LenNumCommas == commas, ]
  toRepMat    <- strsplit(currDataset$Diagnosis, split = ",")
  
  for(diag in seq(1, commas+1))
  {
    txt <- sapply(toRepMat, "[", diag)
    
    if(diag == 1)
    {
      newDataset            <- currDataset
      newDataset$Diagnosis  <- txt
    }else
    {
      tmp           <- currDataset
      tmp$Diagnosis <- txt
      newDataset    <- rbind(newDataset, tmp)
    }
  }
  
  # Size of newDataset should be nrow(currDataset)*commas
  if(nrow(currDataset) * (commas+1) != nrow(newDataset))
  {
    stop("Something went wrong when expanding diagnoses")
  }
  
  dataExpanded <- rbind(dataExpanded, newDataset)
}

# Merge dataExpanded and dataKUHR_one
dataKUHR_All <- rbind(dataKUHR_one, dataExpanded)

# Sort by serial number, year, and diffDays
dataKUHR_All <- dataKUHR_All[order(dataKUHR_All$SerialNumber, dataKUHR_All$Year, dataKUHR_All$DiffDays), ]

# Remove spaces in diagnosis
dataKUHR_All$Diagnosis <- sub(" ", "", dataKUHR_All$Diagnosis)

# It is possible that there are fresh duplicates - prune them
rowsDuplicated  <- duplicated(dataKUHR_All)
dataKUHR_All    <- dataKUHR_All[!rowsDuplicated,]


# Split by diagnostic manual ----------------------------------------------
# First assign missing diagnostic manual
dataKUHR_All$DiagnosticManual[dataKUHR_All$DiagnosticManual == ""] <- "Unspecified"
 
# For the diagnostic manual 2.16.578.1.12.4.1.1.71.7170, there are 14
# observations (11 subjects) with 11 unique diagnostic codes: "T01"  "L87"
# "R75"  "R05"  "A03"  "W12"  "W03"  "H82"  "W99"  "W78"  "W781". Unfortunately,
# these codes are both valid ICPC2 codes as well as valid ICD10 codes, thereby
# making it hard to be sure what this manual corresponds to. To be on the safer
# side, this information is saved separately as another diagnostic manual and
# nothing more is done about it
locStrange        <- dataKUHR_All$DiagnosticManual == "2.16.578.1.12.4.1.1.71.7170"
dataKUHR_Strange  <- dataKUHR_All[locStrange,]

# Separate out ICPC codes
loc_ICPC      <- dataKUHR_All$DiagnosticManual == "ICPC-1" |
                 dataKUHR_All$DiagnosticManual == "ICPC-2" |
                 dataKUHR_All$DiagnosticManual == "ICPC-2B"
dataKUHR_ICPC <- dataKUHR_All[loc_ICPC,]

# Separate out ICD codes
loc_ICD      <- dataKUHR_All$DiagnosticManual == "ICD-10"
dataKUHR_ICD <- dataKUHR_All[loc_ICD,]

# Separate out unspecified
loc_unspecified      <- dataKUHR_All$DiagnosticManual == "Unspecified"
dataKUHR_unspecified <- dataKUHR_All[loc_unspecified,]

# Separate out BUP
loc_BUP      <- dataKUHR_All$DiagnosticManual == "Akser i BUP-klassifikasjon"
dataKUHR_BUP <- dataKUHR_All[loc_BUP,]


# More trash cleaning -----------------------------------------------------
# For ICPC, the following codes are invalid: 01, 2AD, 63, A, J, P
# There is one subject where the parsed diagnosis is "W29, 99" - change this 99 to W99
is99 <- dataKUHR_ICPC$Diagnosis == "99"
if (sum(is99) == 1)
{
  dataKUHR_ICPC$Diagnosis[dataKUHR_ICPC$Diagnosis == "99"] <- "W99"  
}

# Find valid ICPC codes
isValid_ICPC <- grepl("^[A-Z][0-9]{1,3}", dataKUHR_ICPC$Diagnosis)

# Create a separate table for invalid ICPC codes
dataKUHR_ICPC_invalid <- dataKUHR_ICPC[!isValid_ICPC,]

# Delete invalid codes
dataKUHR_ICPC <- dataKUHR_ICPC[isValid_ICPC,]

# Find valid ICD codes
isValid_ICD <- grepl("^[A-Z][0-9]{1,6}", dataKUHR_ICD$Diagnosis)

# Create a separate table for invalid ICD10 codes
dataKUHR_ICD_invalid <- dataKUHR_ICD[!isValid_ICD,]

# Delete invalid codes
dataKUHR_ICD <- dataKUHR_ICD[isValid_ICD,]

# Find valid BUP codes - these should also be similar to ICD10 codes
isValid_BUP <- grepl("^[A-Z][0-9]{1,6}", dataKUHR_BUP$Diagnosis)

# Create a separate table for invalid BUP ICD10 codes
dataKUHR_BUP_invalid <- dataKUHR_BUP[!isValid_BUP,]

# Delete invalid codes
dataKUHR_BUP <- dataKUHR_BUP[isValid_BUP,]


## Save versions of data at this point - before counting
outputName <- paste(as.character(Sys.Date()), "-MoBa-ParsedKUHR-All.csv", sep = "")
fwrite(dataKUHR_All, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)
# write.table(dataKUHR_All, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)

outputName <- paste(as.character(Sys.Date()), "-MoBa-ParsedKUHR-Manual-2.16.578.1.12.4.1.1.71.7170.csv", sep = "")
fwrite(dataKUHR_Strange, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)
# write.table(dataKUHR_Strange, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)

outputName <- paste(as.character(Sys.Date()), "-MoBa-ParsedKUHR-Manual-ICPC.csv", sep = "")
fwrite(dataKUHR_ICPC, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)
# write.table(dataKUHR_ICPC, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)

outputName <- paste(as.character(Sys.Date()), "-MoBa-ParsedKUHR-Manual-ICPC-invalidCodes.csv", sep = "")
fwrite(dataKUHR_ICPC_invalid, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)
# write.table(dataKUHR_ICPC_invalid, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)

outputName <- paste(as.character(Sys.Date()), "-MoBa-ParsedKUHR-Manual-ICD10.csv", sep = "")
fwrite(dataKUHR_ICD, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)
# write.table(dataKUHR_ICD, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)

outputName <- paste(as.character(Sys.Date()), "-MoBa-ParsedKUHR-Manual-ICD10-invalidCodes.csv", sep = "")
fwrite(dataKUHR_ICD_invalid, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)
# write.table(dataKUHR_ICD_invalid, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)

outputName <- paste(as.character(Sys.Date()), "-MoBa-ParsedKUHR-Manual-Akser_i_BUP_klassifikasjon.csv", sep = "")
fwrite(dataKUHR_BUP, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)
# write.table(dataKUHR_BUP, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)

outputName <- paste(as.character(Sys.Date()), "-MoBa-ParsedKUHR-Manual-Akser_i_BUP_klassifikasjon-invalidCodes.csv", sep = "")
fwrite(dataKUHR_BUP_invalid, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)
# write.table(dataKUHR_BUP_invalid, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)

outputName <- paste(as.character(Sys.Date()), "-MoBa-ParsedKUHR-Manual-Unspecified.csv", sep = "")
fwrite(dataKUHR_unspecified, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)
# write.table(dataKUHR_unspecified, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)


# Counting module -------------------------------------------------------
# First, clear the workspace - too much floating around
toDelete <- c("currDataset", "dataExpanded", "dataInfo", "dataKUHR", "dataKUHR_All", 
              "dataKUHR_BUP_invalid", "dataKUHR_ICD_invalid", "dataKUHR_ICPC_invalid",
              "dataKUHR_one", "dataKUHR_Strange", "dataKUHR_unspecified", "dataRemaining", 
              "newDataset", "NumCommas", "tmp", "toRepMat", "withdrawn", "borderline_F_ID",
              "borderline_M_ID", "is99", "isValid_BUP", "isValid_ICD", "isValid_ICPC", 
              "lastChars", "LenNumCommas", "loc_BUP", "loc_ICD", "loc_ICPC", "loc_unspecified", 
              "locComma", "locEliminate", "locsMultiple", "locStrange", "rowsDuplicated")
rm(list = toDelete)

# Temporary section - to delete -----------------------------------------
# dataKUHR_ICD  <- fread(file = "2024-01-11-MoBa-ParsedKUHR-Manual-ICD10.csv", sep = "\t")
# dataKUHR_BUP  <- fread(file = "2024-01-11-MoBa-ParsedKUHR-Manual-Akser_i_BUP_klassifikasjon.csv", sep = "\t")
# dataKUHR_ICPC <- fread(file = "2024-01-11-MoBa-ParsedKUHR-Manual-ICPC.csv", sep = "\t")

# dataKUHR_BUP  <- read.csv(file = "2024-01-11-MoBa-ParsedKUHR-Manual-Akser_i_BUP_klassifikasjon.csv", sep = "\t")


# For counting, merge BUP to ICD
dataKUHR_ICD_BUP <- rbind(dataKUHR_ICD, dataKUHR_BUP)
rm(dataKUHR_BUP)
rm(dataKUHR_ICD)

# Variables to ignore 
varsToIgnore <- c("SerialNumber",                    "SubjectAreaCode", 
                   "PracticeType",                   "DiffDays",
                   "Municipality",                   "County",
                   "TreatmentMunicipality",          "TreatmentCounty",
                   "ContactType",                    "DiagnosticManual",
                   "ParsedDiagnosis")
locsIgnore       <- colnames(dataKUHR_ICPC) %in% varsToIgnore
locsIgnore       <- !locsIgnore
dataKUHR_ICPC    <- dataKUHR_ICPC[, ..locsIgnore] # data.table nonsense

locsIgnore       <- colnames(dataKUHR_ICD_BUP) %in% varsToIgnore
locsIgnore       <- !locsIgnore
dataKUHR_ICD_BUP <- dataKUHR_ICD_BUP[, ..locsIgnore] # data.table nonsense


# Flatten table -----------------------------------------------------------
dupRows_ICD_BUP   <- duplicated(dataKUHR_ICD_BUP)
dataKUHR_ICD_BUP  <- dataKUHR_ICD_BUP[!dupRows_ICD_BUP,]

dupRows_ICPC      <- duplicated(dataKUHR_ICPC)
dataKUHR_ICPC     <- dataKUHR_ICPC[!dupRows_ICPC,]


# Make one overall ID column ----------------------------------------------
dataKUHR_ICD_BUP$ID  <- paste(dataKUHR_ICD_BUP$ID_2445,  dataKUHR_ICD_BUP$M_ID_2445,  dataKUHR_ICD_BUP$F_ID_2445,  sep = "")
dataKUHR_ICPC$ID     <- paste(dataKUHR_ICPC$ID_2445,     dataKUHR_ICPC$M_ID_2445,     dataKUHR_ICPC$F_ID_2445,  sep = "")


# Add Role column ---------------------------------------------------------
locChildren                        <- dataKUHR_ICD_BUP$ID_2445    != ""
locMothers                         <- dataKUHR_ICD_BUP$M_ID_2445  != ""
locFathers                         <- dataKUHR_ICD_BUP$F_ID_2445  != ""
dataKUHR_ICD_BUP$Role[locChildren] <- "Child"
dataKUHR_ICD_BUP$Role[locMothers]  <- "Mother"
dataKUHR_ICD_BUP$Role[locFathers]  <- "Father"

locChildren                        <- dataKUHR_ICPC$ID_2445    != ""
locMothers                         <- dataKUHR_ICPC$M_ID_2445  != ""
locFathers                         <- dataKUHR_ICPC$F_ID_2445  != ""
dataKUHR_ICPC$Role[locChildren]    <- "Child"
dataKUHR_ICPC$Role[locMothers]     <- "Mother"
dataKUHR_ICPC$Role[locFathers]     <- "Father"


# Eliminate other ID columns ----------------------------------------------
# Deal with data.table nonsense
toEliminate      <- !colnames(dataKUHR_ICD_BUP) %in% c("ID_2445", "M_ID_2445", "F_ID_2445")
dataKUHR_ICD_BUP <- dataKUHR_ICD_BUP[, ..toEliminate]

toEliminate      <- !colnames(dataKUHR_ICPC) %in% c("ID_2445", "M_ID_2445", "F_ID_2445")
dataKUHR_ICPC    <- dataKUHR_ICPC[, ..toEliminate]

# dataKUHR_ICD_BUP <- dataKUHR_ICD_BUP[, !(colnames(dataKUHR_ICD_BUP) %in% c("ID_2445", "M_ID_2445", "F_ID_2445"))]
# dataKUHR_ICPC    <- dataKUHR_ICPC[, !(colnames(dataKUHR_ICPC) %in% c("ID_2445", "M_ID_2445", "F_ID_2445"))]


# Rename ID to ID_2445 ----------------------------------------------------
colnames(dataKUHR_ICD_BUP)[colnames(dataKUHR_ICD_BUP) == "ID"] <- "ID_2445"
colnames(dataKUHR_ICPC)[colnames(dataKUHR_ICPC) == "ID"] <- "ID_2445"


# Reorganize tables ----------------------------------------------------------
dataKUHR_ICD_BUP  <- dataKUHR_ICD_BUP[, c("ID_2445", "Role","PREG_ID_2445", "ChildNumber", 
                                          "WithdrawnConsent18years", "Year", "Diagnosis")]

dataKUHR_ICPC  <- dataKUHR_ICPC[, c("ID_2445", "Role","PREG_ID_2445", "ChildNumber", 
                                    "WithdrawnConsent18years", "Year", "Diagnosis")]


# Generate a temporary table having only unique information from both ------
temp        <- rbind(dataKUHR_ICD_BUP[, c("ID_2445", "Role","PREG_ID_2445", "ChildNumber", "WithdrawnConsent18years")],
                     dataKUHR_ICPC[, c("ID_2445", "Role","PREG_ID_2445", "ChildNumber", "WithdrawnConsent18years")])
temp        <- unique(temp)
temp        <- temp[order(temp$ID_2445), ]


# Unable to parse all minor codes - the data frame size becomes too huge; 
# for example, 420547 columns for just ICD+BUP codes
# Unique codes and time points -------------------------------------
uqICDCodes       <- sort(unique(dataKUHR_ICD_BUP$Diagnosis))
uqICDTimePoints  <- sort(unique(dataKUHR_ICD_BUP$Year))


# Convert minor codes to major codes --------------------------------------
allMajorCodesICD     <- unique(substring(uqICDCodes, 1, 3))
allMajorUqCombICD    <- ""
for (time in uqICDTimePoints)
{
  allMajorUqCombICD  <- c(allMajorUqCombICD, paste(allMajorCodesICD, time, sep = "_"))
}
allMajorUqCombICD    <- allMajorUqCombICD[allMajorUqCombICD != ""]
allMajorUqCombICD    <- sort(allMajorUqCombICD)


# Convert ICPC codes to major codes --------------------------------------
uqICPCCodes           <- sort(unique(dataKUHR_ICPC$Diagnosis))
uqICPCTimePoints      <- sort(unique(dataKUHR_ICPC$Year))

allMajorCodesICPC     <- unique(substring(uqICPCCodes, 1, 3))
allMajorUqCombICPC    <- ""
for (time in uqICPCTimePoints)
{
  allMajorUqCombICPC  <- c(allMajorUqCombICPC, paste(allMajorCodesICPC, time, sep = "_"))
}
allMajorUqCombICPC    <- allMajorUqCombICPC[allMajorUqCombICPC != ""]
allMajorUqCombICPC    <- sort(allMajorUqCombICPC)


# Compile information ICD_BUP -----------------------------------------------------
# Initialize results table
allMajColNamesICD             <- c("ID_2445", "Role","PREG_ID_2445", "ChildNumber", 
                                   "WithdrawnConsent18years", paste(allMajorCodesICD, "_lifetime", sep = ""),
                                   allMajorUqCombICD)
allMajColNamesXICD            <- c("ID_2445", "Role","PREG_ID_2445", "ChildNumber", 
                                   "WithdrawnConsent18years", gsub("_", "_x", allMajorUqCombICD))
                                   
MajorDiagInfoICD              <- data.frame(matrix(nrow = nrow(temp), ncol = length(allMajColNamesICD),  data = FALSE))
MajorDiagInfoXICD             <- data.frame(matrix(nrow = nrow(temp), ncol = length(allMajColNamesXICD), data = FALSE))

colnames(MajorDiagInfoICD)    <- allMajColNamesICD
colnames(MajorDiagInfoXICD)   <- allMajColNamesXICD

# Copy information over from temp
MajorDiagInfoICD[, 1:5]       <- temp
MajorDiagInfoXICD[, 1:5]      <- temp

# Count
for (diag in allMajorCodesICD)
{
  for (time in uqICDTimePoints)
  {
    # Diagnoses to search for
    pattern = paste("^", diag, sep = "")
    
    # Find IDs which have this diagnosis up to this time point (cumulative)
    whichIDs  <- dataKUHR_ICD_BUP$ID_2445[grepl(pattern, dataKUHR_ICD_BUP$Diagnosis) & 
                                          dataKUHR_ICD_BUP$Year == time]
    
    # Where are these IDs - logical
    locIDs    <- MajorDiagInfoICD$ID_2445  %in% whichIDs
    locIDsX   <- MajorDiagInfoXICD$ID_2445 %in% whichIDs
    
    # Columns to mark
    whichCols  <- paste(diag, time:2023, sep = "_")
    whichColsX <- paste(diag, time, sep = "_x")
    
    # Mark
    MajorDiagInfoICD[locIDs,  whichCols]  <- TRUE
    MajorDiagInfoXICD[locIDs, whichColsX] <- TRUE
  }
  # Info for lifetime - sum across all time points
  # Which columns to sum across
  whichCols       <- paste(diag, uqICDTimePoints, sep = "_")
  
  # Where to mark?
  toMark          <- paste(diag, "_lifetime", sep = "")
  
  # Mark
  MajorDiagInfoICD[, toMark] <- as.logical(rowSums(MajorDiagInfoICD[, whichCols]))
  
  # Update user
  print(paste("Finished ", diag, sep = ""))
}

MajorDiagInfoICD[, 6:ncol(MajorDiagInfoICD)][]  <- lapply(MajorDiagInfoICD[, 6:ncol(MajorDiagInfoICD)],  as.integer)
MajorDiagInfoXICD[,6:ncol(MajorDiagInfoXICD)][] <- lapply(MajorDiagInfoXICD[,6:ncol(MajorDiagInfoXICD)], as.integer)

outputName <- paste(as.character(Sys.Date()), "-MoBa-ParsedKUHR-ICD10-majorCodes.csv", sep = "")
fwrite(MajorDiagInfoICD, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)
# write.table(MajorDiagInfoICD, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)

outputName <- paste(as.character(Sys.Date()), "-MoBa-ParsedKUHR-ICD10-majorCodesX.csv", sep = "")
fwrite(MajorDiagInfoXICD, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)
# write.table(MajorDiagInfoXICD, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)

rm(MajorDiagInfoICD)
rm(MajorDiagInfoXICD)
gc()


# Compile information ICPC -----------------------------------------------------
# Initialize results table
allMajColNamesICPC            <- c("ID_2445", "Role","PREG_ID_2445", "ChildNumber", 
                                   "WithdrawnConsent18years", paste(allMajorCodesICPC, "_lifetime", sep = ""),
                                   allMajorUqCombICPC)
allMajColNamesXICPC           <- c("ID_2445", "Role","PREG_ID_2445", "ChildNumber", 
                                   "WithdrawnConsent18years", gsub("_", "_x", allMajorUqCombICPC))

MajorDiagInfoICPC             <- data.frame(matrix(nrow = nrow(temp), ncol = length(allMajColNamesICPC),  data = FALSE))
MajorDiagInfoXICPC            <- data.frame(matrix(nrow = nrow(temp), ncol = length(allMajColNamesXICPC), data = FALSE))

colnames(MajorDiagInfoICPC)   <- allMajColNamesICPC
colnames(MajorDiagInfoXICPC)  <- allMajColNamesXICPC

# Copy information over from temp
MajorDiagInfoICPC[, 1:5]      <- temp
MajorDiagInfoXICPC[, 1:5]     <- temp

# Count
startTime <- proc.time()
for (diag in allMajorCodesICPC)
{
  for (time in uqICPCTimePoints)
  {
    # Diagnoses to search for
    pattern = paste("^", diag, sep = "")
    
    # Find IDs which have this diagnosis up to this time point (cumulative)
    whichIDs  <- dataKUHR_ICPC$ID_2445[grepl(pattern, dataKUHR_ICPC$Diagnosis) & 
                                            dataKUHR_ICPC$Year == time]
    
    # Where are these IDs - logical
    locIDs    <- MajorDiagInfoICPC$ID_2445  %in% whichIDs
    locXIDs   <- MajorDiagInfoXICPC$ID_2445 %in% whichIDs
    
    # Columns to mark
    whichCols  <- paste(diag, time:2023, sep = "_")
    whichColsX <- paste(diag, time, sep = "_x")
    
    # Mark
    MajorDiagInfoICPC[locIDs,  whichCols]  <- TRUE
    MajorDiagInfoXICPC[locIDs, whichColsX] <- TRUE
  }
  
  # Info for lifetime - sum across all time points
  # Which columns to sum across
  whichCols       <- paste(diag, uqICPCTimePoints, sep = "_")
  
  # Where to mark?
  toMark          <- paste(diag, "_lifetime", sep = "")
  
  # Mark
  MajorDiagInfoICPC[, toMark] <- as.logical(rowSums(MajorDiagInfoICPC[, whichCols]))
  
  # Update user
  print(paste("Finished ", diag, sep = ""))
}

MajorDiagInfoICPC[, 6:ncol(MajorDiagInfoICPC)][]  <- lapply(MajorDiagInfoICPC[, 6:ncol(MajorDiagInfoICPC)],  as.integer)
MajorDiagInfoXICPC[,6:ncol(MajorDiagInfoXICPC)][] <- lapply(MajorDiagInfoXICPC[,6:ncol(MajorDiagInfoXICPC)], as.integer)

outputName <- paste(as.character(Sys.Date()), "-MoBa-ParsedKUHR-ICPC-majorCodes.csv", sep = "")
fwrite(MajorDiagInfoICPC, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)
# write.table(MajorDiagInfoICPC, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)

outputName <- paste(as.character(Sys.Date()), "-MoBa-ParsedKUHR-ICPC-majorCodesX.csv", sep = "")
fwrite(MajorDiagInfoXICPC, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)
# write.table(MajorDiagInfoXICPC, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)
