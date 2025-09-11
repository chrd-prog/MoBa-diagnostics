# Do counts ---------------------------------------------------------------
library(data.table)

# Read data
dirKPR                <- "/ess/p697/cluster/users/parekh/2025-05-22_KPR"
dirGenetics           <- "/tsd/p697/data/durable/genotype/MoBaPsychGen_v1"
dirOutput             <- file.path(dirKPR, "counts")
fnameGenCov           <- "MoBaPsychGen_v1-ec-eur-batch-basic-qc-cov.txt"
listICD               <- list.files(path = dirKPR, pattern = glob2rx("*-MoBa-ParsedKPR_ICD10Codes.csv"))
ICDCodes              <- fread(file.path(dirKPR, listICD[length(listICD)]), sep = "\t")
infoGenetics          <- fread(file.path(dirGenetics, fnameGenCov), sep = "\t")

# Get all codes and time points -------------------------------------------
# We know that the first five column names are not diagnostic categories
allNames              <- colnames(ICDCodes)[6:length(colnames(ICDCodes))]
temp                  <- strsplit(allNames, "_")
allCodes              <- unique(sapply(temp, "[", 1))
allTimes              <- unique(sapply(temp, "[", 2))

# Initialize --------------------------------------------------------------
newNames              <- c("ICD10Code", allTimes)

allCounts_children    <- as.data.frame(matrix(data = 0, nrow = length(allCodes),  ncol = length(newNames)))
allCounts_mothers     <- as.data.frame(matrix(data = 0, nrow = length(allCodes),  ncol = length(newNames)))
allCounts_fathers     <- as.data.frame(matrix(data = 0, nrow = length(allCodes),  ncol = length(newNames)))
allCounts_allSubjs    <- as.data.frame(matrix(data = 0, nrow = length(allCodes),  ncol = length(newNames)))

colnames(allCounts_children)    <- newNames
colnames(allCounts_mothers)     <- newNames
colnames(allCounts_fathers)     <- newNames
colnames(allCounts_allSubjs)    <- newNames


# Count minor/all codes ---------------------------------------------------
# Get a column-wise count across all columns
tempSums_children     <- colSums(ICDCodes[ICDCodes$Role == "Child",  6:ncol(ICDCodes)])
tempSums_mothers      <- colSums(ICDCodes[ICDCodes$Role == "Mother", 6:ncol(ICDCodes)])
tempSums_fathers      <- colSums(ICDCodes[ICDCodes$Role == "Father", 6:ncol(ICDCodes)])
tempSums_allSubjs     <- colSums(ICDCodes[,                          6:ncol(ICDCodes)])

# Get diagnoses names and times
tempTimes             <- attributes(tempSums_children)$names
tempDiag              <- sapply(strsplit(tempTimes, "_"), "[", 1)
tempTimes             <- sapply(strsplit(tempTimes, "_"), "[", 2)

# Go over time points and summarize counts
for (times in allTimes)
{
  # Find location of these times
  locTimes            <- tempTimes == times
  
  # Ensure diagnostic ordering is the same
  if (sum(tempDiag[locTimes] == allCodes) != length(allCodes))
  {
    stop("Something went wrong in the ordering of codes")
  }
  
  # Add to data frame
  allCounts_children[,times]  <- tempSums_children[locTimes]
  allCounts_mothers[,times]   <- tempSums_mothers[locTimes]
  allCounts_fathers[,times]   <- tempSums_fathers[locTimes]
  allCounts_allSubjs[,times]  <- tempSums_allSubjs[locTimes]
}

# Assign ICPC codes
allCounts_children$ICD10Code  <- allCodes
allCounts_mothers$ICD10Code   <- allCodes
allCounts_fathers$ICD10Code   <- allCodes
allCounts_allSubjs$ICD10Code  <- allCodes

# Sanity checks: allSubjs = children + mothers + fathers
if (sum(rowSums(
  as.matrix(allCounts_allSubjs[, 2:ncol(allCounts_allSubjs)]) - 
  (as.matrix(allCounts_children[,2:ncol(allCounts_children)])  + 
   as.matrix(allCounts_mothers[, 2:ncol(allCounts_mothers)])   +
   as.matrix(allCounts_fathers[, 2:ncol(allCounts_fathers)])))) != 0)
{
  stop("Mismatch between all subjects and sum of mothers, fathers, and children")
}


# Merge with genetics -----------------------------------------------------
subsetICDCodes    <- ICDCodes[ICDCodes$ID_2445 %in% infoGenetics$ID_2445, ]


# Initialize --------------------------------------------------------------
genCounts_children       <- as.data.frame(matrix(data = 0, nrow = length(allCodes),      ncol = length(newNames)))
genCounts_mothers        <- as.data.frame(matrix(data = 0, nrow = length(allCodes),      ncol = length(newNames)))
genCounts_fathers        <- as.data.frame(matrix(data = 0, nrow = length(allCodes),      ncol = length(newNames)))
genCounts_allSubjs       <- as.data.frame(matrix(data = 0, nrow = length(allCodes),      ncol = length(newNames)))

colnames(genCounts_children)       <- newNames
colnames(genCounts_mothers)        <- newNames
colnames(genCounts_fathers)        <- newNames
colnames(genCounts_allSubjs)       <- newNames

# Count minor/all codes for subjects with genetics -------------------------
# Get a column-wise count across all columns
tempSums_children     <- colSums(subsetICDCodes[subsetICDCodes$Role == "Child",  6:ncol(subsetICDCodes)])
tempSums_mothers      <- colSums(subsetICDCodes[subsetICDCodes$Role == "Mother", 6:ncol(subsetICDCodes)])
tempSums_fathers      <- colSums(subsetICDCodes[subsetICDCodes$Role == "Father", 6:ncol(subsetICDCodes)])
tempSums_allSubjs     <- colSums(subsetICDCodes[,                                6:ncol(subsetICDCodes)])

# Get diagnoses names and times
tempTimes             <- attributes(tempSums_children)$names
tempDiag              <- sapply(strsplit(tempTimes, "_"), "[", 1)
tempTimes             <- sapply(strsplit(tempTimes, "_"), "[", 2)

# Go over time points and summarize counts
for (times in allTimes)
{
  # Find location of these times
  locTimes            <- tempTimes == times
  
  # Ensure diagnostic ordering is the same
  if (sum(tempDiag[locTimes] == allCodes) != length(allCodes))
  {
    stop("Something went wrong in the ordering of codes")
  }
  
  # Add to data frame
  genCounts_children[,times]  <- tempSums_children[locTimes]
  genCounts_mothers[,times]   <- tempSums_mothers[locTimes]
  genCounts_fathers[,times]   <- tempSums_fathers[locTimes]
  genCounts_allSubjs[,times]  <- tempSums_allSubjs[locTimes]
}

# Assign ICPC codes
genCounts_children$ICD10Code  <- allCodes
genCounts_mothers$ICD10Code   <- allCodes
genCounts_fathers$ICD10Code   <- allCodes
genCounts_allSubjs$ICD10Code  <- allCodes

# Sanity checks: allSubjs = children + mothers + fathers
if (sum(rowSums(
  as.matrix(genCounts_allSubjs[, 2:ncol(genCounts_allSubjs)]) - 
  (as.matrix(genCounts_children[,2:ncol(genCounts_children)])  + 
   as.matrix(genCounts_mothers[, 2:ncol(genCounts_mothers)])   +
   as.matrix(genCounts_fathers[, 2:ncol(genCounts_fathers)])))) != 0)
{
  stop("Mismatch between all subjects and sum of mothers, fathers, and children")
}


# Save results ------------------------------------------------------------
if (!dir.exists(dirOutput))
{
  dir.create(dirOutput)
}

# Write out counts of all codes - all subjects
outputName  <- paste(as.character(Sys.Date()), "-MoBa-KPRCounts_allICDCodes-allSubjects.csv", sep = "")
write.table(allCounts_allSubjs, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)

# Write out counts of all codes - children
outputName  <- paste(as.character(Sys.Date()), "-MoBa-KPRCounts_allICDCodes-Children.csv", sep = "")
write.table(allCounts_children, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)

# Write out counts of all codes - mothers
outputName  <- paste(as.character(Sys.Date()), "-MoBa-KPRCounts_allICDCodes-Mothers.csv", sep = "")
write.table(allCounts_mothers, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)

# Write out counts of all codes - fathers
outputName  <- paste(as.character(Sys.Date()), "-MoBa-KPRCounts_allICDCodes-Fathers.csv", sep = "")
write.table(allCounts_fathers, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)


# Repeat saving for counts of subjects with genetics data -----------------
# Write out counts of all codes - all subjects
outputName  <- paste(as.character(Sys.Date()), "-MoBa-KPRGeneticsCounts_allICDCodes-allSubjects.csv", sep = "")
write.table(genCounts_allSubjs, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)

# Write out counts of all codes - children
outputName  <- paste(as.character(Sys.Date()), "-MoBa-KPRGeneticsCounts_allICDCodes-Children.csv", sep = "")
write.table(genCounts_children, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)

# Write out counts of all codes - mothers
outputName  <- paste(as.character(Sys.Date()), "-MoBa-KPRGeneticsCounts_allICDCodes-Mothers.csv", sep = "")
write.table(genCounts_mothers, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)

# Write out counts of all codes - fathers
outputName  <- paste(as.character(Sys.Date()), "-MoBa-KPRGeneticsCounts_allICDCodes-Fathers.csv", sep = "")
write.table(genCounts_fathers, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)
