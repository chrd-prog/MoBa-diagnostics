# Merge KPR data
library(data.table)

# File locations ----------------------------------------------------------
# Where output files should be written
dirWorking <- "/ess/p697/cluster/users/parekh/2025-05-22_KPR"
setwd(dirWorking)

# Current location of KPR files
dirKPR <- "/tsd/p697/data/durable/phenotypes/kpr/2025-05-15/rawdata"

# The takst file contains billing information - do we care? Likely not, ignoring
data_diagnoses <- fread(file.path(dirKPR, "diagnose.csv"), sep = ";")
data_billing   <- fread(file.path(dirKPR, "regning.csv"),  sep = ";")

# Merge
data_KPR <- merge(data_diagnoses, data_billing, by = "nokkel", all = TRUE)

# Remove duplicate lopenr column and rename
wch <- !colnames(data_KPR) %in% "lopenr.y"
data_KPR  <- data_KPR[,  ..wch]
colnames(data_KPR)[colnames(data_KPR)   == "lopenr.x"] <- "lopenr"


# Renaming column names ----------------------------------------------------
# This is based mostly on the folgebrev
colnames(data_KPR)[colnames(data_KPR) == "kommunenr"]         <- "MunicipalityNumber"
colnames(data_KPR)[colnames(data_KPR) == "praksiskommune"]    <- "TreatmentMunicipality"
colnames(data_KPR)[colnames(data_KPR) == "fylkenr"]           <- "CountyNumber"
colnames(data_KPR)[colnames(data_KPR) == "praksisfylkenr"]    <- "TreatmentCounty"
colnames(data_KPR)[colnames(data_KPR) == "diagnosetabell"]    <- "DiagnosisType"
colnames(data_KPR)[colnames(data_KPR) == "diagnosekode"]      <- "DiagnosisCode"
colnames(data_KPR)[colnames(data_KPR) == "erhoveddiagnose"]   <- "isMainDiagnosis"
colnames(data_KPR)[colnames(data_KPR) == "aar"]               <- "Year"
colnames(data_KPR)[colnames(data_KPR) == "diffdager_inndato"] <- "DiffDays_Admission"
colnames(data_KPR)[colnames(data_KPR) == "tjenestetype"]      <- "ServiceType"
colnames(data_KPR)[colnames(data_KPR) == "lopenr"]            <- "SerialNumber"
colnames(data_KPR)[colnames(data_KPR) == "nokkel"]            <- "Key"

# Save merged KPR --------------------------------------------------------
outputName <- paste(as.character(Sys.Date()), "-MoBa-MergedKPR.csv", sep = "")
fwrite(data_KPR, file = file.path(dirWorking, outputName), sep = "\t", row.names = F)
