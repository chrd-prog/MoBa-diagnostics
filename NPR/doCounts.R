# Do counts ---------------------------------------------------------------
library(data.table)

# Read data
dirNPR                <- "/ess/p697/cluster/users/parekh/2025-05-21_NPR"
dirGenetics           <- "/tsd/p697/data/durable/genotype/MoBaPsychGen_v1"
dirOutput             <- file.path(dirNPR, "counts")
fnameGenCov           <- "MoBaPsychGen_v1-ec-eur-batch-basic-qc-cov.txt"
listMinor             <- list.files(path = dirNPR, pattern = glob2rx("*-MoBa-ParsedNPR_allCodes.csv"))
listMajor             <- list.files(path = dirNPR, pattern = glob2rx("*-MoBa-ParsedNPR_majorCodes.csv"))
minorCodes            <- fread(file.path(dirNPR, listMinor[length(listMinor)]), sep = "\t")
majorCodes            <- fread(file.path(dirNPR, listMajor[length(listMajor)]), sep = "\t")
infoGenetics          <- fread(file.path(dirGenetics, fnameGenCov), sep = "\t")
# minorCodes            <- read.csv(file.path(dirNPR, listMinor[length(listMinor)]), sep = "\t")
# majorCodes            <- read.csv(file.path(dirNPR, listMajor[length(listMajor)]), sep = "\t")
# infoGenetics          <- read.csv(file.path(dirGenetics, fnameGenCov), sep = "\t")

# Get all codes and time points -------------------------------------------
# We know that the first five column names are not diagnostic categories
allNames              <- colnames(minorCodes)[6:length(colnames(minorCodes))]
temp                  <- strsplit(allNames, "_")
allCodes              <- unique(sapply(temp, "[", 1))
allTimes              <- unique(sapply(temp, "[", 2))

# All this can be done in a single line
allMajorCodes         <- unique(sapply(strsplit(colnames(majorCodes)[6:length(colnames(majorCodes))], "_"), "[", 1))


# Initialize --------------------------------------------------------------
newNames              <- c("ICD10Code", allTimes)

allCounts_children    <- as.data.frame(matrix(data = 0, nrow = length(allCodes),      ncol = length(newNames)))
allCounts_mothers     <- as.data.frame(matrix(data = 0, nrow = length(allCodes),      ncol = length(newNames)))
allCounts_fathers     <- as.data.frame(matrix(data = 0, nrow = length(allCodes),      ncol = length(newNames)))
allCounts_allSubjs    <- as.data.frame(matrix(data = 0, nrow = length(allCodes),      ncol = length(newNames)))

majorCounts_children  <- as.data.frame(matrix(data = 0, nrow = length(allMajorCodes), ncol = length(newNames)))
majorCounts_mothers   <- as.data.frame(matrix(data = 0, nrow = length(allMajorCodes), ncol = length(newNames)))
majorCounts_fathers   <- as.data.frame(matrix(data = 0, nrow = length(allMajorCodes), ncol = length(newNames)))
majorCounts_allSubjs  <- as.data.frame(matrix(data = 0, nrow = length(allMajorCodes), ncol = length(newNames)))

colnames(allCounts_children)    <- newNames
colnames(allCounts_mothers)     <- newNames
colnames(allCounts_fathers)     <- newNames
colnames(allCounts_allSubjs)    <- newNames

colnames(majorCounts_children)  <- newNames
colnames(majorCounts_mothers)   <- newNames
colnames(majorCounts_fathers)   <- newNames
colnames(majorCounts_allSubjs)  <- newNames


# Count minor/all codes ---------------------------------------------------
# Get a column-wise count across all columns
tempSums_children     <- colSums(minorCodes[minorCodes$Role == "Child",  6:ncol(minorCodes)])
tempSums_mothers      <- colSums(minorCodes[minorCodes$Role == "Mother", 6:ncol(minorCodes)])
tempSums_fathers      <- colSums(minorCodes[minorCodes$Role == "Father", 6:ncol(minorCodes)])
tempSums_allSubjs     <- colSums(minorCodes[,                            6:ncol(minorCodes)])

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

# Assign ICD10 codes
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


# Count major codes ---------------------------------------------------
# Get a column-wise count across all columns
tempSums_children     <- colSums(majorCodes[majorCodes$Role == "Child",  6:ncol(majorCodes)])
tempSums_mothers      <- colSums(majorCodes[majorCodes$Role == "Mother", 6:ncol(majorCodes)])
tempSums_fathers      <- colSums(majorCodes[majorCodes$Role == "Father", 6:ncol(majorCodes)])
tempSums_allSubjs     <- colSums(majorCodes[,                            6:ncol(majorCodes)])

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
  if (sum(tempDiag[locTimes] == allMajorCodes) != length(allMajorCodes))
  {
    stop("Something went wrong in the ordering of codes")
  }
  
  # Add to data frame
  majorCounts_children[,times]  <- tempSums_children[locTimes]
  majorCounts_mothers[,times]   <- tempSums_mothers[locTimes]
  majorCounts_fathers[,times]   <- tempSums_fathers[locTimes]
  majorCounts_allSubjs[,times]  <- tempSums_allSubjs[locTimes]
}

# Assign ICD10 codes
majorCounts_children$ICD10Code  <- allMajorCodes
majorCounts_mothers$ICD10Code   <- allMajorCodes
majorCounts_fathers$ICD10Code   <- allMajorCodes
majorCounts_allSubjs$ICD10Code  <- allMajorCodes

# Sanity checks: allSubjs = children + mothers + fathers
if (sum(rowSums(
  as.matrix(majorCounts_allSubjs[, 2:ncol(majorCounts_allSubjs)]) - 
  (as.matrix(majorCounts_children[,2:ncol(majorCounts_children)])  + 
   as.matrix(majorCounts_mothers[, 2:ncol(majorCounts_mothers)])   +
   as.matrix(majorCounts_fathers[, 2:ncol(majorCounts_fathers)])))) != 0)
{
  stop("Mismatch between all subjects and sum of mothers, fathers, and children")
}


# Merge with genetics -----------------------------------------------------
subsetMinorCodes    <- minorCodes[minorCodes$ID_2445 %in% infoGenetics$ID_2445, ]
subsetMajorCodes    <- majorCodes[majorCodes$ID_2445 %in% infoGenetics$ID_2445, ]

# Alternate way of merging tables
# This will, of course, add columns from genetics to the merged table
# altSubsetMinorCode  <- merge(x = minorCodes, y = infoGenetics, by = "ID_2445")
# altSubsetMajorCode  <- merge(x = majorCodes, y = infoGenetics, by = "ID_2445")


# Initialize --------------------------------------------------------------
genCounts_children       <- as.data.frame(matrix(data = 0, nrow = length(allCodes),      ncol = length(newNames)))
genCounts_mothers        <- as.data.frame(matrix(data = 0, nrow = length(allCodes),      ncol = length(newNames)))
genCounts_fathers        <- as.data.frame(matrix(data = 0, nrow = length(allCodes),      ncol = length(newNames)))
genCounts_allSubjs       <- as.data.frame(matrix(data = 0, nrow = length(allCodes),      ncol = length(newNames)))

majorCountsGen_children  <- as.data.frame(matrix(data = 0, nrow = length(allMajorCodes), ncol = length(newNames)))
majorCountsGen_mothers   <- as.data.frame(matrix(data = 0, nrow = length(allMajorCodes), ncol = length(newNames)))
majorCountsGen_fathers   <- as.data.frame(matrix(data = 0, nrow = length(allMajorCodes), ncol = length(newNames)))
majorCountsGen_allSubjs  <- as.data.frame(matrix(data = 0, nrow = length(allMajorCodes), ncol = length(newNames)))

colnames(genCounts_children)       <- newNames
colnames(genCounts_mothers)        <- newNames
colnames(genCounts_fathers)        <- newNames
colnames(genCounts_allSubjs)       <- newNames

colnames(majorCountsGen_children)  <- newNames
colnames(majorCountsGen_mothers)   <- newNames
colnames(majorCountsGen_fathers)   <- newNames
colnames(majorCountsGen_allSubjs)  <- newNames


# Count minor/all codes for subjects with genetics -------------------------
# Get a column-wise count across all columns
tempSums_children     <- colSums(subsetMinorCodes[subsetMinorCodes$Role == "Child",  6:ncol(subsetMinorCodes)])
tempSums_mothers      <- colSums(subsetMinorCodes[subsetMinorCodes$Role == "Mother", 6:ncol(subsetMinorCodes)])
tempSums_fathers      <- colSums(subsetMinorCodes[subsetMinorCodes$Role == "Father", 6:ncol(subsetMinorCodes)])
tempSums_allSubjs     <- colSums(subsetMinorCodes[,                                  6:ncol(subsetMinorCodes)])

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

# Assign ICD10 codes
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


# Recount major codes for subjects with genetics --------------------------
# Get a column-wise count across all columns
tempSums_children     <- colSums(subsetMajorCodes[subsetMajorCodes$Role == "Child",  6:ncol(subsetMajorCodes)])
tempSums_mothers      <- colSums(subsetMajorCodes[subsetMajorCodes$Role == "Mother", 6:ncol(subsetMajorCodes)])
tempSums_fathers      <- colSums(subsetMajorCodes[subsetMajorCodes$Role == "Father", 6:ncol(subsetMajorCodes)])
tempSums_allSubjs     <- colSums(subsetMajorCodes[,                                  6:ncol(subsetMajorCodes)])

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
  if (sum(tempDiag[locTimes] == allMajorCodes) != length(allMajorCodes))
  {
    stop("Something went wrong in the ordering of codes")
  }
  
  # Add to data frame
  majorCountsGen_children[,times]  <- tempSums_children[locTimes]
  majorCountsGen_mothers[,times]   <- tempSums_mothers[locTimes]
  majorCountsGen_fathers[,times]   <- tempSums_fathers[locTimes]
  majorCountsGen_allSubjs[,times]  <- tempSums_allSubjs[locTimes]
}

# Assign ICD10 codes
majorCountsGen_children$ICD10Code  <- allMajorCodes
majorCountsGen_mothers$ICD10Code   <- allMajorCodes
majorCountsGen_fathers$ICD10Code   <- allMajorCodes
majorCountsGen_allSubjs$ICD10Code  <- allMajorCodes

# Sanity checks: allSubjs = children + mothers + fathers
if (sum(rowSums(
  as.matrix(majorCountsGen_allSubjs[, 2:ncol(majorCountsGen_allSubjs)]) - 
  (as.matrix(majorCountsGen_children[,2:ncol(majorCountsGen_children)])  + 
   as.matrix(majorCountsGen_mothers[, 2:ncol(majorCountsGen_mothers)])   +
   as.matrix(majorCountsGen_fathers[, 2:ncol(majorCountsGen_fathers)])))) != 0)
{
  stop("Mismatch between all subjects and sum of mothers, fathers, and children")
}


# Save results ------------------------------------------------------------
if (!dir.exists(dirOutput))
{
  dir.create(dirOutput)
}

# Write out counts of all codes - all subjects
outputName  <- paste(as.character(Sys.Date()), "-MoBa-NPRCounts_allCodes-allSubjects.csv", sep = "")
write.table(allCounts_allSubjs, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)

# Write out counts of all codes - children
outputName  <- paste(as.character(Sys.Date()), "-MoBa-NPRCounts_allCodes-Children.csv", sep = "")
write.table(allCounts_children, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)

# Write out counts of all codes - mothers
outputName  <- paste(as.character(Sys.Date()), "-MoBa-NPRCounts_allCodes-Mothers.csv", sep = "")
write.table(allCounts_mothers, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)

# Write out counts of all codes - fathers
outputName  <- paste(as.character(Sys.Date()), "-MoBa-NPRCounts_allCodes-Fathers.csv", sep = "")
write.table(allCounts_fathers, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)

# Write out counts of major codes - all subjects
outputName  <- paste(as.character(Sys.Date()), "-MoBa-NPRCounts_majorCodes-allSubjects.csv", sep = "")
write.table(majorCounts_allSubjs, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)

# Write out counts of major codes - children
outputName  <- paste(as.character(Sys.Date()), "-MoBa-NPRCounts_majorCodes-Children.csv", sep = "")
write.table(majorCounts_children, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)

# Write out counts of major codes - mothers
outputName  <- paste(as.character(Sys.Date()), "-MoBa-NPRCounts_majorCodes-Mothers.csv", sep = "")
write.table(majorCounts_mothers, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)

# Write out counts of major codes - fathers
outputName  <- paste(as.character(Sys.Date()), "-MoBa-NPRCounts_majorCodes-Fathers.csv", sep = "")
write.table(majorCounts_fathers, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)


# Repeat saving for counts of subjects with genetics data -----------------

# Write out counts of all codes - all subjects
outputName  <- paste(as.character(Sys.Date()), "-MoBa-NPRGeneticsCounts_allCodes-allSubjects.csv", sep = "")
write.table(genCounts_allSubjs, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)

# Write out counts of all codes - children
outputName  <- paste(as.character(Sys.Date()), "-MoBa-NPRGeneticsCounts_allCodes-Children.csv", sep = "")
write.table(genCounts_children, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)

# Write out counts of all codes - mothers
outputName  <- paste(as.character(Sys.Date()), "-MoBa-NPRGeneticsCounts_allCodes-Mothers.csv", sep = "")
write.table(genCounts_mothers, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)

# Write out counts of all codes - fathers
outputName  <- paste(as.character(Sys.Date()), "-MoBa-NPRGeneticsCounts_allCodes-Fathers.csv", sep = "")
write.table(genCounts_fathers, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)

# Write out counts of major codes - all subjects
outputName  <- paste(as.character(Sys.Date()), "-MoBa-NPRGeneticsCounts_majorCodes-allSubjects.csv", sep = "")
write.table(majorCountsGen_allSubjs, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)

# Write out counts of major codes - children
outputName  <- paste(as.character(Sys.Date()), "-MoBa-NPRGeneticsCounts_majorCodes-Children.csv", sep = "")
write.table(majorCountsGen_children, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)

# Write out counts of major codes - mothers
outputName  <- paste(as.character(Sys.Date()), "-MoBa-NPRGeneticsCounts_majorCodes-Mothers.csv", sep = "")
write.table(majorCountsGen_mothers, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)

# Write out counts of major codes - fathers
outputName  <- paste(as.character(Sys.Date()), "-MoBa-NPRGeneticsCounts_majorCodes-Fathers.csv", sep = "")
write.table(majorCountsGen_fathers, file = file.path(dirOutput, outputName), sep = "\t", row.names = F)
