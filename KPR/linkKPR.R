# Link KPR
# Library for reading SPSS file (link files are SPSS files)
library(haven)
library(data.table)

# File locations ----------------------------------------------------------
# Where output files should be written
dirWorking <- "/ess/p697/cluster/users/parekh/2025-05-22_KPR"
setwd(dirWorking)

# Current location of link files
dirLink <- "/tsd/p697/data/durable/phenotypes/sv_infofiles/koblingsbro/"

# File name for the merged KPR
mergedFile <- paste(as.character(Sys.Date()), "-MoBa-MergedKPR.csv", sep = "")

# Read files --------------------------------------------------------------
dataKPR     <- fread(file.path(dirWorking, mergedFile), sep = "\t")
link_mother <- read_sav(file.path(dirLink, list.files(path = dirLink, pattern = glob2rx("Mor*ID*2445*.sav"))))
link_father <- read_sav(file.path(dirLink, list.files(path = dirLink, pattern = glob2rx("Far*ID*2445*.sav"))))
link_child  <- read_sav(file.path(dirLink, list.files(path = dirLink, pattern = glob2rx("Barn*ID*2445*.sav"))))

# Join files --------------------------------------------------------------
# Here, we perform an inner join on the tables - entries which are not found are
# dropped; this naturally translates into three separate tables - one for
# mothers, one for fathers, and one for children
KPR_mothers  <- merge(x = dataKPR, y = link_mother, by.x = "SerialNumber", by.y = "PID_1037_20250425")
KPR_fathers  <- merge(x = dataKPR, y = link_father, by.x = "SerialNumber", by.y = "PID_1037_20250425")
KPR_children <- merge(x = dataKPR, y = link_child,  by.x = "SerialNumber", by.y = "PID_1037_20250425")

stopifnot(nrow(KPR_children) + nrow(KPR_mothers) + nrow(KPR_fathers) - nrow(dataKPR) == 0)


# Housekeeping prior to merging -------------------------------------------
# Since the children KPR have two more columns than parents (BARN_NR and
# OppfStopp18aar), we rename these columns
colnames(KPR_children)[colnames(KPR_children) == "OppfStopp18aar"] <- "WithdrawnConsent18years"
colnames(KPR_children)[colnames(KPR_children) == "BARN_NR"]        <- "ChildNumber"

# Save files - independent dataset ----------------------------------------
outputName <- paste(as.character(Sys.Date()), "-MoBa-LinkedKPR-Mothers.csv", sep = "")
fwrite(KPR_mothers, file = file.path(dirWorking, outputName), sep = "\t", row.names = F)

outputName <- paste(as.character(Sys.Date()), "-MoBa-LinkedKPR-Fathers.csv", sep = "")
fwrite(KPR_fathers, file = file.path(dirWorking, outputName), sep = "\t", row.names = F)

outputName <- paste(as.character(Sys.Date()), "-MoBa-LinkedKPR-Children.csv", sep = "")
fwrite(KPR_children, file = file.path(dirWorking, outputName), sep = "\t", row.names = F)


# Merge tables ------------------------------------------------------------
# We add WithdrawnConsent18years and ChildNumber to the parents' KPR as
# prior to merging; additionally, prior to merging, additional columns like 
# F_ID_2445 and M_ID_2445 need to be added to the children's KPR (and partners' 
# ID to other partners' table)

KPR_mothers$WithdrawnConsent18years <- ""
KPR_mothers$ChildNumber             <- ""
KPR_mothers$F_ID_2445               <- ""
KPR_mothers$PREG_ID_2445            <- ""

KPR_fathers$WithdrawnConsent18years <- ""
KPR_fathers$ChildNumber             <- ""
KPR_fathers$M_ID_2445               <- ""
KPR_fathers$PREG_ID_2445            <- ""

KPR_children$M_ID_2445              <- ""
KPR_children$F_ID_2445              <- ""

# Merge everything --------------------------------------------------------
# Side effect: file size is huge and therefore takes a while to read...
KPR_All <- rbind(rbind(KPR_children, KPR_mothers), KPR_fathers)

# Save files --------------------------------------------------------------
outputName <- paste(as.character(Sys.Date()), "-MoBa-LinkedKPR-AllSubjects.csv", sep = "")
fwrite(KPR_All, file = file.path(dirWorking, outputName), sep = "\t", row.names = F)