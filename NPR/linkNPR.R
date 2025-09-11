# Script to link merged NPR file with ID2445

# Library for reading SPSS file (link files are SPSS files)
library(haven)
library(data.table)

# File locations ----------------------------------------------------------
# Where output files should be written
dirWorking <- "/ess/p697/cluster/users/parekh/2025-05-21_NPR"
setwd(dirWorking)

# Current location of link files
dirNPR <- "/tsd/p697/data/durable/phenotypes/sv_infofiles/koblingsbro/"

# File name for the merged NPR
mergedFile <- paste(as.character(Sys.Date()), "-MoBa-MergedNPR.csv", sep = "")

# Read files --------------------------------------------------------------
dataNPR     <- fread(file.path(dirWorking, mergedFile), sep = "\t")
link_mother <- read_sav(file.path(dirNPR, list.files(path = dirNPR, pattern = glob2rx("Mor*ID*2445*.sav"))))
link_father <- read_sav(file.path(dirNPR, list.files(path = dirNPR, pattern = glob2rx("Far*ID*2445*.sav"))))
link_child  <- read_sav(file.path(dirNPR, list.files(path = dirNPR, pattern = glob2rx("Barn*ID*2445*.sav"))))

# Join files --------------------------------------------------------------
# Here, we perform an inner join on the tables - entries which are not found are
# dropped; this naturally translates into three separate tables - one for
# mothers, one for fathers, and one for children
NPR_mothers  <- merge(x = dataNPR, y = link_mother, by.x = "SerialNumber", by.y = "PID_1037_20250425")
NPR_fathers  <- merge(x = dataNPR, y = link_father, by.x = "SerialNumber", by.y = "PID_1037_20250425")
NPR_children <- merge(x = dataNPR, y = link_child,  by.x = "SerialNumber", by.y = "PID_1037_20250425")

stopifnot(nrow(NPR_children) + nrow(NPR_mothers) + nrow(NPR_fathers) - nrow(dataNPR) == 0)

# Housekeeping prior to merging -------------------------------------------
# Since the children NPR have two more columns than parents (BARN_NR and
# OppfStopp18aar), we rename these columns
colnames(NPR_children)[colnames(NPR_children) == "OppfStopp18aar"] <- "WithdrawnConsent18years"
colnames(NPR_children)[colnames(NPR_children) == "BARN_NR"]        <- "ChildNumber"

# Save files - independent dataset ----------------------------------------
outputName <- paste(as.character(Sys.Date()), "-MoBa-LinkedNPR-Mothers.csv", sep = "")
fwrite(NPR_mothers, file = file.path(dirWorking, outputName), sep = "\t", row.names = F)
# write.table(NPR_mothers, file = file.path(dirWorking, outputName), sep = "\t", row.names = F)

outputName <- paste(as.character(Sys.Date()), "-MoBa-LinkedNPR-Fathers.csv", sep = "")
fwrite(NPR_fathers, file = file.path(dirWorking, outputName), sep = "\t", row.names = F)
# write.table(NPR_fathers, file = file.path(dirWorking, outputName), sep = "\t", row.names = F)

outputName <- paste(as.character(Sys.Date()), "-MoBa-LinkedNPR-Children.csv", sep = "")
fwrite(NPR_children, file = file.path(dirWorking, outputName), sep = "\t", row.names = F)
# write.table(NPR_children, file = file.path(dirWorking, outputName), sep = "\t", row.names = F)

# Merge tables ------------------------------------------------------------
# We add WithdrawnConsent18years and ChildNumber to the parents' NPR as
# prior to merging; additionally, prior to merging, additional columns like 
# F_ID_2445 and M_ID_2445 need to be added to the children's NPR (and partners' 
# ID to other partners' table)

NPR_mothers$WithdrawnConsent18years <- ""
NPR_mothers$ChildNumber             <- ""
NPR_mothers$F_ID_2445               <- ""
NPR_mothers$PREG_ID_2445            <- ""

NPR_fathers$WithdrawnConsent18years <- ""
NPR_fathers$ChildNumber             <- ""
NPR_fathers$M_ID_2445               <- ""
NPR_fathers$PREG_ID_2445            <- ""

NPR_children$M_ID_2445              <- ""
NPR_children$F_ID_2445              <- ""

# Merge everything --------------------------------------------------------
# Side effect: file size is huge and therefore takes a while to read...
NPR_All <- rbind(rbind(NPR_children, NPR_mothers), NPR_fathers)

# Save files --------------------------------------------------------------
outputName <- paste(as.character(Sys.Date()), "-MoBa-LinkedNPR-AllSubjects.csv", sep = "")
fwrite(NPR_All, file = file.path(dirWorking, outputName), sep = "\t", row.names = F)
# write.table(NPR_All, file = file.path(dirWorking, outputName), sep = "\t", row.names = F)