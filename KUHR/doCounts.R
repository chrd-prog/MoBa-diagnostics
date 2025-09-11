# Do counts ---------------------------------------------------------------
library(data.table)

# Preliminaries
dirKUHR               <- "/ess/p697/cluster/users/parekh/2023-10-24_KUHR"
dirGenetics           <- "/tsd/p697/data/durable/genotype/MoBaPsychGen_v1"
dirOutput             <- file.path(dirKUHR, "counts")
fnameGenCov           <- "MoBaPsychGen_v1-ec-eur-batch-basic-qc-cov.txt"
listICD               <- list.files(path = dirKUHR, pattern = glob2rx("*-MoBa-ParsedKUHR-ICD10-majorCodes.csv"))
listICPC              <- list.files(path = dirKUHR, pattern = glob2rx("*-MoBa-ParsedKUHR-ICPC-majorCodes.csv"))
delimiter             <- "\t"

# Read data
dataICD               <- fread(file = listICD[length(listICD)],  sep = delimiter, header = TRUE)
dataICPC              <- fread(file = listICPC[length(listICD)], sep = delimiter, header = TRUE)

# Genetic information
infoGenetics          <- read.csv(file.path(dirGenetics, fnameGenCov), sep = "\t")

# Get all codes and time points -------------------------------------------
# We know that the first five column names are not diagnostic categories
# allNamesICD         <- colnames(dataICD)[6:length(colnames(dataICD))]
# tempICD             <- strsplit(allNamesICD, "_")
# allCodesICD         <- unique(sapply(tempICD, "[", 1))
# allTimesICD         <- unique(sapply(tempICD, "[", 2))

# All this can be done in a single line
allMajorCodesICD      <- unique(sapply(strsplit(colnames(dataICD)[6:length(colnames(dataICD))], "_"), "[", 1))
allMajorCodesICPC     <- unique(sapply(strsplit(colnames(dataICPC)[6:length(colnames(dataICPC))], "_"), "[", 1))
allTimesICD           <- unique(sapply(strsplit(colnames(dataICD)[6:length(colnames(dataICD))], "_"), "[", 2))
allTimesICPC          <- unique(sapply(strsplit(colnames(dataICPC)[6:length(colnames(dataICPC))], "_"), "[", 2))

# Initialize --------------------------------------------------------------
newNamesICD             <- c("ICD10Code", allTimesICD)
newNamesICPC            <- c("ICPCCode",  allTimesICPC)

majorCounts_ICD_children  <- as.data.frame(matrix(data = 0, nrow = length(allMajorCodesICD),  ncol = length(newNamesICD)))
majorCounts_ICD_mothers   <- as.data.frame(matrix(data = 0, nrow = length(allMajorCodesICD),  ncol = length(newNamesICD)))
majorCounts_ICD_fathers   <- as.data.frame(matrix(data = 0, nrow = length(allMajorCodesICD),  ncol = length(newNamesICD)))
majorCounts_ICD_allSubjs  <- as.data.frame(matrix(data = 0, nrow = length(allMajorCodesICD),  ncol = length(newNamesICD)))

majorCounts_ICPC_children <- as.data.frame(matrix(data = 0, nrow = length(allMajorCodesICPC), ncol = length(newNamesICPC)))
majorCounts_ICPC_mothers  <- as.data.frame(matrix(data = 0, nrow = length(allMajorCodesICPC), ncol = length(newNamesICPC)))
majorCounts_ICPC_fathers  <- as.data.frame(matrix(data = 0, nrow = length(allMajorCodesICPC), ncol = length(newNamesICPC)))
majorCounts_ICPC_allSubjs <- as.data.frame(matrix(data = 0, nrow = length(allMajorCodesICPC), ncol = length(newNamesICPC)))

colnames(majorCounts_ICD_children)  <- newNamesICD
colnames(majorCounts_ICD_mothers)   <- newNamesICD
colnames(majorCounts_ICD_fathers)   <- newNamesICD
colnames(majorCounts_ICD_allSubjs)  <- newNamesICD

colnames(majorCounts_ICPC_children) <- newNamesICPC
colnames(majorCounts_ICPC_mothers)  <- newNamesICPC
colnames(majorCounts_ICPC_fathers)  <- newNamesICPC
colnames(majorCounts_ICPC_allSubjs) <- newNamesICPC


# Count major codes ---------------------------------------------------
# Get a column-wise count across all columns
tempSums_ICD_children     <- colSums(dataICD[dataICD$Role == "Child",    6:ncol(dataICD)])
tempSums_ICD_mothers      <- colSums(dataICD[dataICD$Role == "Mother",   6:ncol(dataICD)])
tempSums_ICD_fathers      <- colSums(dataICD[dataICD$Role == "Father",   6:ncol(dataICD)])
tempSums_ICD_allSubjs     <- colSums(dataICD[,                           6:ncol(dataICD)])

tempSums_ICPC_children    <- colSums(dataICPC[dataICPC$Role == "Child",  6:ncol(dataICPC)])
tempSums_ICPC_mothers     <- colSums(dataICPC[dataICPC$Role == "Mother", 6:ncol(dataICPC)])
tempSums_ICPC_fathers     <- colSums(dataICPC[dataICPC$Role == "Father", 6:ncol(dataICPC)])
tempSums_ICPC_allSubjs    <- colSums(dataICPC[,                          6:ncol(dataICPC)])


# Get diagnoses names and times
tempTimes_ICD         <- attributes(tempSums_ICD_children)$names
tempDiag_ICD          <- sapply(strsplit(tempTimes_ICD, "_"), "[", 1)
tempTimes_ICD         <- sapply(strsplit(tempTimes_ICD, "_"), "[", 2)

tempTimes_ICPC        <- attributes(tempSums_ICPC_children)$names
tempDiag_ICPC         <- sapply(strsplit(tempTimes_ICPC, "_"), "[", 1)
tempTimes_ICPC        <- sapply(strsplit(tempTimes_ICPC, "_"), "[", 2)

# Go over time points and summarize counts - ICD version
for (times in allTimesICD)
{
  # Find location of these times
  locTimes            <- tempTimes_ICD == times
  
  # Ensure diagnostic ordering is the same
  if (sum(tempDiag_ICD[locTimes] == allMajorCodesICD) != length(allMajorCodesICD))
  {
    stop("Something went wrong in the ordering of ICD codes")
  }
  
  # Add to data frame
  majorCounts_ICD_children[,times]  <- tempSums_ICD_children[locTimes]
  majorCounts_ICD_mothers[,times]   <- tempSums_ICD_mothers[locTimes]
  majorCounts_ICD_fathers[,times]   <- tempSums_ICD_fathers[locTimes]
  majorCounts_ICD_allSubjs[,times]  <- tempSums_ICD_allSubjs[locTimes]
}

# Assign ICD10 codes
majorCounts_ICD_children$ICD10Code  <- allMajorCodesICD
majorCounts_ICD_mothers$ICD10Code   <- allMajorCodesICD
majorCounts_ICD_fathers$ICD10Code   <- allMajorCodesICD
majorCounts_ICD_allSubjs$ICD10Code  <- allMajorCodesICD

# Sanity checks: allSubjs = children + mothers + fathers
if (sum(rowSums(
  as.matrix(majorCounts_ICD_allSubjs[, 2:ncol(majorCounts_ICD_allSubjs)]) - 
  (as.matrix(majorCounts_ICD_children[,2:ncol(majorCounts_ICD_children)])  + 
   as.matrix(majorCounts_ICD_mothers[, 2:ncol(majorCounts_ICD_mothers)])   +
   as.matrix(majorCounts_ICD_fathers[, 2:ncol(majorCounts_ICD_fathers)])))) != 0)
{
  stop("Mismatch between ICD all subjects and sum of mothers, fathers, and children")
}


# Go over time points and summarize counts - ICPC version ------------------
for (times in allTimesICPC)
{
  # Find location of these times
  locTimes            <- tempTimes_ICPC == times
  
  # Ensure diagnostic ordering is the same
  if (sum(tempDiag_ICPC[locTimes] == allMajorCodesICPC) != length(allMajorCodesICPC))
  {
    stop("Something went wrong in the ordering of ICPC codes")
  }
  
  # Add to data frame
  majorCounts_ICPC_children[,times]  <- tempSums_ICPC_children[locTimes]
  majorCounts_ICPC_mothers[,times]   <- tempSums_ICPC_mothers[locTimes]
  majorCounts_ICPC_fathers[,times]   <- tempSums_ICPC_fathers[locTimes]
  majorCounts_ICPC_allSubjs[,times]  <- tempSums_ICPC_allSubjs[locTimes]
}

# Assign ICPC codes
majorCounts_ICPC_children$ICPCCode  <- allMajorCodesICPC
majorCounts_ICPC_mothers$ICPCCode   <- allMajorCodesICPC
majorCounts_ICPC_fathers$ICPCCode   <- allMajorCodesICPC
majorCounts_ICPC_allSubjs$ICPCCode  <- allMajorCodesICPC

# Sanity checks: allSubjs = children + mothers + fathers
if (sum(rowSums(
  as.matrix(majorCounts_ICPC_allSubjs[, 2:ncol(majorCounts_ICPC_allSubjs)]) - 
  (as.matrix(majorCounts_ICPC_children[,2:ncol(majorCounts_ICPC_children)])  + 
   as.matrix(majorCounts_ICPC_mothers[, 2:ncol(majorCounts_ICPC_mothers)])   +
   as.matrix(majorCounts_ICPC_fathers[, 2:ncol(majorCounts_ICPC_fathers)])))) != 0)
{
  stop("Mismatch between ICPC all subjects and sum of mothers, fathers, and children")
}


# Merge with genetics -----------------------------------------------------
subsetMajorCodes_ICD    <- dataICD[dataICD$ID_2445   %in% infoGenetics$ID_2445, ]
subsetMajorCodes_ICPC   <- dataICPC[dataICPC$ID_2445 %in% infoGenetics$ID_2445, ]

# Alternate way of merging tables
# This will, of course, add columns from genetics to the merged table
# altSubsetMinorCode  <- merge(x = minorCodes, y = infoGenetics, by = "ID_2445")
# altSubsetMajorCode  <- merge(x = majorCodes, y = infoGenetics, by = "ID_2445")


# Initialize --------------------------------------------------------------
majorCountsGen_ICD_children   <- as.data.frame(matrix(data = 0, nrow = length(allMajorCodesICD),  ncol = length(newNamesICD)))
majorCountsGen_ICD_mothers    <- as.data.frame(matrix(data = 0, nrow = length(allMajorCodesICD),  ncol = length(newNamesICD)))
majorCountsGen_ICD_fathers    <- as.data.frame(matrix(data = 0, nrow = length(allMajorCodesICD),  ncol = length(newNamesICD)))
majorCountsGen_ICD_allSubjs   <- as.data.frame(matrix(data = 0, nrow = length(allMajorCodesICD),  ncol = length(newNamesICD)))

majorCountsGen_ICPC_children  <- as.data.frame(matrix(data = 0, nrow = length(allMajorCodesICPC), ncol = length(newNamesICPC)))
majorCountsGen_ICPC_mothers   <- as.data.frame(matrix(data = 0, nrow = length(allMajorCodesICPC), ncol = length(newNamesICPC)))
majorCountsGen_ICPC_fathers   <- as.data.frame(matrix(data = 0, nrow = length(allMajorCodesICPC), ncol = length(newNamesICPC)))
majorCountsGen_ICPC_allSubjs  <- as.data.frame(matrix(data = 0, nrow = length(allMajorCodesICPC), ncol = length(newNamesICPC)))

colnames(majorCountsGen_ICD_children)  <- newNamesICD
colnames(majorCountsGen_ICD_mothers)   <- newNamesICD
colnames(majorCountsGen_ICD_fathers)   <- newNamesICD
colnames(majorCountsGen_ICD_allSubjs)  <- newNamesICD

colnames(majorCountsGen_ICPC_children)  <- newNamesICPC
colnames(majorCountsGen_ICPC_mothers)   <- newNamesICPC
colnames(majorCountsGen_ICPC_fathers)   <- newNamesICPC
colnames(majorCountsGen_ICPC_allSubjs)  <- newNamesICPC


# Count ICD10 major codes for subjects with genetics --------------------------
# Get a column-wise count across all columns
tempSums_ICD_children     <- colSums(subsetMajorCodes_ICD[subsetMajorCodes_ICD$Role == "Child",  6:ncol(subsetMajorCodes_ICD)])
tempSums_ICD_mothers      <- colSums(subsetMajorCodes_ICD[subsetMajorCodes_ICD$Role == "Mother", 6:ncol(subsetMajorCodes_ICD)])
tempSums_ICD_fathers      <- colSums(subsetMajorCodes_ICD[subsetMajorCodes_ICD$Role == "Father", 6:ncol(subsetMajorCodes_ICD)])
tempSums_ICD_allSubjs     <- colSums(subsetMajorCodes_ICD[,                                      6:ncol(subsetMajorCodes_ICD)])

# Get diagnoses names and times
tempTimes_ICD             <- attributes(tempSums_ICD_children)$names
tempDiag_ICD              <- sapply(strsplit(tempTimes_ICD, "_"), "[", 1)
tempTimes_ICD             <- sapply(strsplit(tempTimes_ICD, "_"), "[", 2)

# Go over time points and summarize counts
for (times in allTimesICD)
{
  # Find location of these times
  locTimes            <- tempTimes_ICD == times
  
  # Ensure diagnostic ordering is the same
  if (sum(tempDiag_ICD[locTimes] == allMajorCodesICD) != length(allMajorCodesICD))
  {
    stop("Something went wrong in the ordering of codes")
  }
  
  # Add to data frame
  majorCountsGen_ICD_children[,times]  <- tempSums_ICD_children[locTimes]
  majorCountsGen_ICD_mothers[,times]   <- tempSums_ICD_mothers[locTimes]
  majorCountsGen_ICD_fathers[,times]   <- tempSums_ICD_fathers[locTimes]
  majorCountsGen_ICD_allSubjs[,times]  <- tempSums_ICD_allSubjs[locTimes]
}

# Assign ICD10 codes
majorCountsGen_ICD_children$ICD10Code  <- allMajorCodesICD
majorCountsGen_ICD_mothers$ICD10Code   <- allMajorCodesICD
majorCountsGen_ICD_fathers$ICD10Code   <- allMajorCodesICD
majorCountsGen_ICD_allSubjs$ICD10Code  <- allMajorCodesICD

# Sanity checks: allSubjs = children + mothers + fathers
if (sum(rowSums(
  as.matrix(majorCountsGen_ICD_allSubjs[, 2:ncol(majorCountsGen_ICD_allSubjs)]) - 
  (as.matrix(majorCountsGen_ICD_children[,2:ncol(majorCountsGen_ICD_children)])  + 
   as.matrix(majorCountsGen_ICD_mothers[, 2:ncol(majorCountsGen_ICD_mothers)])   +
   as.matrix(majorCountsGen_ICD_fathers[, 2:ncol(majorCountsGen_ICD_fathers)])))) != 0)
{
  stop("Mismatch between all subjects and sum of mothers, fathers, and children")
}


# Count ICPC major codes for subjects with genetics --------------------------
# Get a column-wise count across all columns
tempSums_ICPC_children     <- colSums(subsetMajorCodes_ICPC[subsetMajorCodes_ICPC$Role == "Child",  6:ncol(subsetMajorCodes_ICPC)])
tempSums_ICPC_mothers      <- colSums(subsetMajorCodes_ICPC[subsetMajorCodes_ICPC$Role == "Mother", 6:ncol(subsetMajorCodes_ICPC)])
tempSums_ICPC_fathers      <- colSums(subsetMajorCodes_ICPC[subsetMajorCodes_ICPC$Role == "Father", 6:ncol(subsetMajorCodes_ICPC)])
tempSums_ICPC_allSubjs     <- colSums(subsetMajorCodes_ICPC[,                                       6:ncol(subsetMajorCodes_ICPC)])

# Get diagnoses names and times
tempTimes_ICPC             <- attributes(tempSums_ICPC_children)$names
tempDiag_ICPC              <- sapply(strsplit(tempTimes_ICPC, "_"), "[", 1)
tempTimes_ICPC             <- sapply(strsplit(tempTimes_ICPC, "_"), "[", 2)

# Go over time points and summarize counts
for (times in allTimesICPC)
{
  # Find location of these times
  locTimes            <- tempTimes_ICPC == times
  
  # Ensure diagnostic ordering is the same
  if (sum(tempDiag_ICPC[locTimes] == allMajorCodesICPC) != length(allMajorCodesICPC))
  {
    stop("Something went wrong in the ordering of codes")
  }
  
  # Add to data frame
  majorCountsGen_ICPC_children[,times]  <- tempSums_ICPC_children[locTimes]
  majorCountsGen_ICPC_mothers[,times]   <- tempSums_ICPC_mothers[locTimes]
  majorCountsGen_ICPC_fathers[,times]   <- tempSums_ICPC_fathers[locTimes]
  majorCountsGen_ICPC_allSubjs[,times]  <- tempSums_ICPC_allSubjs[locTimes]
}

# Assign ICPC10 codes
majorCountsGen_ICPC_children$ICPCCode  <- allMajorCodesICPC
majorCountsGen_ICPC_mothers$ICPCCode   <- allMajorCodesICPC
majorCountsGen_ICPC_fathers$ICPCCode   <- allMajorCodesICPC
majorCountsGen_ICPC_allSubjs$ICPCCode  <- allMajorCodesICPC

# Sanity checks: allSubjs = children + mothers + fathers
if (sum(rowSums(
  as.matrix(majorCountsGen_ICPC_allSubjs[, 2:ncol(majorCountsGen_ICPC_allSubjs)]) - 
  (as.matrix(majorCountsGen_ICPC_children[,2:ncol(majorCountsGen_ICPC_children)])  + 
   as.matrix(majorCountsGen_ICPC_mothers[, 2:ncol(majorCountsGen_ICPC_mothers)])   +
   as.matrix(majorCountsGen_ICPC_fathers[, 2:ncol(majorCountsGen_ICPC_fathers)])))) != 0)
{
  stop("Mismatch between all subjects and sum of mothers, fathers, and children")
}


# Save results ------------------------------------------------------------
if (!dir.exists(dirOutput))
{
  dir.create(dirOutput)
}

# Save ICD10 results ------------------------------------------------------
# Write out counts of ICD10 major codes - all subjects
outputName  <- paste(as.character(Sys.Date()), "-MoBa-KUHRCounts_majorCodes-ICD10-allSubjects.csv", sep = "")
fwrite(majorCounts_ICD_allSubjs, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)

# Write out counts of ICD10 major codes - children
outputName  <- paste(as.character(Sys.Date()), "-MoBa-KUHRCounts_majorCodes-ICD10-Children.csv", sep = "")
fwrite(majorCounts_ICD_children, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)

# Write out counts of ICD10 major codes - mothers
outputName  <- paste(as.character(Sys.Date()), "-MoBa-KUHRCounts_majorCodes-ICD10-Mothers.csv", sep = "")
fwrite(majorCounts_ICD_mothers, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)

# Write out counts of ICD10 major codes - fathers
outputName  <- paste(as.character(Sys.Date()), "-MoBa-KUHRCounts_majorCodes-ICD10-Fathers.csv", sep = "")
fwrite(majorCounts_ICD_fathers, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)


# Repeat saving for counts of subjects with genetics data -----------------
# Write out counts of major codes - all subjects
outputName  <- paste(as.character(Sys.Date()), "-MoBa-KUHRGeneticsCounts_majorCodes-ICD10-allSubjects.csv", sep = "")
fwrite(majorCountsGen_ICD_allSubjs, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)

# Write out counts of major codes - children
outputName  <- paste(as.character(Sys.Date()), "-MoBa-KUHRGeneticsCounts_majorCodes-ICD10-Children.csv", sep = "")
fwrite(majorCountsGen_ICD_children, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)

# Write out counts of major codes - mothers
outputName  <- paste(as.character(Sys.Date()), "-MoBa-KUHRGeneticsCounts_majorCodes-ICD10-Mothers.csv", sep = "")
fwrite(majorCountsGen_ICD_mothers, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)

# Write out counts of major codes - fathers
outputName  <- paste(as.character(Sys.Date()), "-MoBa-KUHRGeneticsCounts_majorCodes-ICD10-Fathers.csv", sep = "")
fwrite(majorCountsGen_ICD_fathers, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)


# Save ICPC results ------------------------------------------------------
# Write out counts of ICPC major codes - all subjects
outputName  <- paste(as.character(Sys.Date()), "-MoBa-KUHRCounts_majorCodes-ICPC-allSubjects.csv", sep = "")
fwrite(majorCounts_ICPC_allSubjs, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)

# Write out counts of ICPC major codes - children
outputName  <- paste(as.character(Sys.Date()), "-MoBa-KUHRCounts_majorCodes-ICPC-Children.csv", sep = "")
fwrite(majorCounts_ICPC_children, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)

# Write out counts of ICPC major codes - mothers
outputName  <- paste(as.character(Sys.Date()), "-MoBa-KUHRCounts_majorCodes-ICPC-Mothers.csv", sep = "")
fwrite(majorCounts_ICPC_mothers, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)

# Write out counts of ICPC major codes - fathers
outputName  <- paste(as.character(Sys.Date()), "-MoBa-KUHRCounts_majorCodes-ICPC-Fathers.csv", sep = "")
fwrite(majorCounts_ICPC_fathers, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)


# Repeat saving for counts of subjects with genetics data -----------------
# Write out counts of major codes - all subjects
outputName  <- paste(as.character(Sys.Date()), "-MoBa-KUHRGeneticsCounts_majorCodes-ICPC-allSubjects.csv", sep = "")
fwrite(majorCountsGen_ICPC_allSubjs, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)

# Write out counts of major codes - children
outputName  <- paste(as.character(Sys.Date()), "-MoBa-KUHRGeneticsCounts_majorCodes-ICPC-Children.csv", sep = "")
fwrite(majorCountsGen_ICPC_children, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)

# Write out counts of major codes - mothers
outputName  <- paste(as.character(Sys.Date()), "-MoBa-KUHRGeneticsCounts_majorCodes-ICPC-Mothers.csv", sep = "")
fwrite(majorCountsGen_ICPC_mothers, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)

# Write out counts of major codes - fathers
outputName  <- paste(as.character(Sys.Date()), "-MoBa-KUHRGeneticsCounts_majorCodes-ICPC-Fathers.csv", sep = "")
fwrite(majorCountsGen_ICPC_fathers, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)
