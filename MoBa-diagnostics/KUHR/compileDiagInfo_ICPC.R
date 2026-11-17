# Compile info on first contact for each diagnoses
library(data.table)

# Read parsed NPR
dirKUHR  <- "/ess/p697/cluster/users/parekh/2023-10-24_KUHR"
dataKUHR <- fread(file.path(dirKUHR, "2025-08-13-MoBa-ParsedKUHR-Manual-ICPC.csv"), sep = "\t")

# Variables to ignore -----------------------------------------------------
varsToIgnore <- c("SerialNumber", "SubjectAreaCode", "PracticeType",
                  "Municipality", "County", "TreatmentMunicipality",
                  "TreatmentCounty", "ContactType", "ParsedDiagnosis", "DiagnosticManual")

locsIgnore <- colnames(dataKUHR) %in% varsToIgnore
locsIgnore <- !locsIgnore
dataKUHR   <- dataKUHR[, ..locsIgnore] # data.table nonsense

# Remove any row without ICPC code ---------------------------------------
locEmpty   <- is.na(dataKUHR$Diagnosis) | dataKUHR$Diagnosis == ""
dataKUHR   <- dataKUHR[!locEmpty,]

# Flatten table -----------------------------------------------------------
dupRows     <- duplicated(dataKUHR)
dataKUHR    <- dataKUHR[!dupRows,]

# Make one overall ID column ----------------------------------------------
dataKUHR$ID  <- paste(dataKUHR$ID_2445,  dataKUHR$M_ID_2445,  dataKUHR$F_ID_2445,  sep = "")

# Add Role column ---------------------------------------------------------
locChildren               <- dataKUHR$ID_2445    != ""
locMothers                <- dataKUHR$M_ID_2445  != ""
locFathers                <- dataKUHR$F_ID_2445  != ""
dataKUHR$Role[locChildren] <- "Child"
dataKUHR$Role[locMothers]  <- "Mother"
dataKUHR$Role[locFathers]  <- "Father"

# Eliminate other ID columns ----------------------------------------------
toEliminate <- !colnames(dataKUHR) %in% c("ID_2445", "M_ID_2445", "F_ID_2445")
dataKUHR     <- dataKUHR[, ..toEliminate]

# Rename ID to ID_2445 ----------------------------------------------------
colnames(dataKUHR)[colnames(dataKUHR) == "ID"] <- "ID_2445"

# Reorganize KUHR ----------------------------------------------------------
dataKUHR     <- dataKUHR[, c("ID_2445", "Role","PREG_ID_2445", "ChildNumber", 
                             "WithdrawnConsent18years", "Year", "Diagnosis",
                             "DiffDays")]

# Unique ICPC codes and time points -------------------------------------
uqICPCCodes   <- sort(unique(dataKUHR$Diagnosis))
uqTimePoints  <- sort(unique(dataKUHR$Year))

# List of unique subjects
uqSubjs <- unique(dataKUHR$ID_2445)

# Initialize results table
diagInformation  <- as.data.frame(matrix(data = NA, nrow = length(uqSubjs), ncol = length(uqICPCCodes)+5))
tempSanityCounts <- as.data.frame(matrix(data = 0,  nrow = length(uqSubjs), ncol = length(uqICPCCodes)+5))
colnames(diagInformation) <- c("ID_2445", "Role","PREG_ID_2445", "ChildNumber", 
                               "WithdrawnConsent18years", uqICPCCodes)

colnames(tempSanityCounts) <- c("ID_2445", "Role","PREG_ID_2445", "ChildNumber", 
                                "WithdrawnConsent18years", uqICPCCodes)

# Minimal dataKUHR table
tmp_minData <- dataKUHR[, 1:5]
dupRows     <- duplicated(tmp_minData)
tmp_minData <- tmp_minData[!dupRows,]

# Assign to diagInformation
diagInformation[, 1:5]  <- tmp_minData
tempSanityCounts[, 1:5] <- tmp_minData

# Compile
for (codes in uqICPCCodes)
{
  # All rows in dataKUHR having this diagnosis
  tmpRows <- dataKUHR$Diagnosis == codes
  
  # For every found subject, work
  tmpData <- dataKUHR[tmpRows,]
  tmpUq   <- unique(tmpData$ID_2445)
  
  for (ids in tmpUq)
  {
    toMark <- diagInformation$ID_2445 == ids
    loc    <- tmpData$ID_2445 == ids
    
    # Save results
    diagInformation[toMark,  eval(codes)] <- min(tmpData$DiffDays[loc])
    tempSanityCounts[toMark, eval(codes)] <- 1
  }
  print(codes)
}

# Save
outputName  <- paste(as.character(Sys.Date()), "-MoBa-SmallestDiffDate-ICPCallCodes.csv", sep = "")
fwrite(diagInformation, file = file.path(dirKUHR, outputName), sep = "\t", row.names = F)