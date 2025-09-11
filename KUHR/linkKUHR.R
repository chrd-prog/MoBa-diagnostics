# Script to link KUHR data to ID2445
# Library for reading SPSS file (link files are SPSS files)
library(haven)

# File locations ----------------------------------------------------------
# Where output files should be written
dirWorking <- "/ess/p697/cluster/users/parekh/2023-10-24_KUHR"
setwd(dirWorking)

# Current location of link files
dirKUHR <- "/tsd/p697/data/durable/phenotypes/kuhr/Kobling/"

# File name for the merged KUHR
mergedFile <- paste(as.character(Sys.Date()), "-MoBa-MergedKUHR.csv", sep = "")


# Read files --------------------------------------------------------------
dataKUHR    <- read.csv(file.path(dirWorking, mergedFile), sep = ",")
link_mother <- read_sav(file.path(dirKUHR, list.files(path = dirKUHR, pattern = glob2rx("Mor*ID*2445*.sav"))))
link_father <- read_sav(file.path(dirKUHR, list.files(path = dirKUHR, pattern = glob2rx("Far*ID*2445*.sav"))))
link_child  <- read_sav(file.path(dirKUHR, list.files(path = dirKUHR, pattern = glob2rx("Barn*ID*2445*.sav"))))


# Join files --------------------------------------------------------------
# Here, we perform an inner join on the tables - entries which are not found are
# dropped; this naturally translates into three separate tables - one for
# mothers, one for fathers, and one for children
KUHR_mothers  <- merge(x = dataKUHR, y = link_mother, by.x = "SerialNumber", by.y = "PID_2445")
KUHR_fathers  <- merge(x = dataKUHR, y = link_father, by.x = "SerialNumber", by.y = "PID_2445")
KUHR_children <- merge(x = dataKUHR, y = link_child,  by.x = "SerialNumber", by.y = "PID_2445")


# Housekeeping prior to merging -------------------------------------------
# Since the children KUHR have two more columns than parents (BARN_NR and
# OppfStopp18aar), we rename these columns
colnames(KUHR_children)[colnames(KUHR_children) == "OppfStopp18aar"] <- "WithdrawnConsent18years"
colnames(KUHR_children)[colnames(KUHR_children) == "BARN_NR"]        <- "ChildNumber"

# Save files - independent dataset ----------------------------------------
outputName <- paste(as.character(Sys.Date()), "-MoBa-LinkedKUHR-Mothers.csv", sep = "")
write.table(KUHR_mothers, file = file.path(dirWorking, outputName), sep = ",", row.names = F)

outputName <- paste(as.character(Sys.Date()), "-MoBa-LinkedKUHR-Fathers.csv", sep = "")
write.table(KUHR_fathers, file = file.path(dirWorking, outputName), sep = ",", row.names = F)

outputName <- paste(as.character(Sys.Date()), "-MoBa-LinkedKUHR-Children.csv", sep = "")
write.table(KUHR_children, file = file.path(dirWorking, outputName), sep = ",", row.names = F)

# Merge tables ------------------------------------------------------------
# We add WithdrawnConsent18years and ChildNumber to the parents' KUHR as
# prior to merging; additionally, prior to merging, additional columns like 
# F_ID_2445 and M_ID_2445 need to be added to the children's KUHR (and partners' 
# ID to other partners' table)

KUHR_mothers$WithdrawnConsent18years <- ""
KUHR_mothers$ChildNumber             <- ""
KUHR_mothers$F_ID_2445               <- ""
KUHR_mothers$PREG_ID_2445            <- ""

KUHR_fathers$WithdrawnConsent18years <- ""
KUHR_fathers$ChildNumber             <- ""
KUHR_fathers$M_ID_2445               <- ""
KUHR_fathers$PREG_ID_2445            <- ""

KUHR_children$M_ID_2445              <- ""
KUHR_children$F_ID_2445              <- ""

# Merge everything --------------------------------------------------------
# Side effect: file size is huge and therefore takes a while to read...
KUHR_All <- rbind(rbind(KUHR_children, KUHR_mothers), KUHR_fathers)

# Save files --------------------------------------------------------------
outputName <- paste(as.character(Sys.Date()), "-MoBa-LinkedKUHR-AllSubjects.csv", sep = "")
write.table(KUHR_All, file = file.path(dirWorking, outputName), sep = ",", row.names = F)
