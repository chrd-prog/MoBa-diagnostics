# Parse KUHR data and convert to a usable format

# File locations ----------------------------------------------------------
# Where output files should be written
dirWorking <- "/ess/p697/cluster/users/parekh/2023-10-24_KUHR"
setwd(dirWorking)

# Current location of KUHR files
dirKUHR <- "/tsd/p697/data/durable/phenotypes/kuhr/2006-2023 Data fra KUHR 21-50907/"


# Read individual files - these are .dsv files -----------------------------
data2006 <- read.csv(file.path(dirKUHR, list.files(path = dirKUHR, pattern = glob2rx("2006*.dsv"))),  sep = ";")
data2007 <- read.csv(file.path(dirKUHR, list.files(path = dirKUHR, pattern = glob2rx("2007*.dsv"))),  sep = ";")
data2008 <- read.csv(file.path(dirKUHR, list.files(path = dirKUHR, pattern = glob2rx("2008*.dsv"))),  sep = ";")
data2009 <- read.csv(file.path(dirKUHR, list.files(path = dirKUHR, pattern = glob2rx("2009*.dsv"))),  sep = ";")
data2010 <- read.csv(file.path(dirKUHR, list.files(path = dirKUHR, pattern = glob2rx("2010*.dsv"))),  sep = ";")
data2011 <- read.csv(file.path(dirKUHR, list.files(path = dirKUHR, pattern = glob2rx("2011*.dsv"))),  sep = ";")
data2012 <- read.csv(file.path(dirKUHR, list.files(path = dirKUHR, pattern = glob2rx("2012*.dsv"))),  sep = ";")
data2013 <- read.csv(file.path(dirKUHR, list.files(path = dirKUHR, pattern = glob2rx("2013*.dsv"))),  sep = ";")
data2014 <- read.csv(file.path(dirKUHR, list.files(path = dirKUHR, pattern = glob2rx("2014*.dsv"))),  sep = ";")
data2015 <- read.csv(file.path(dirKUHR, list.files(path = dirKUHR, pattern = glob2rx("2015*.dsv"))),  sep = ";")
data2016 <- read.csv(file.path(dirKUHR, list.files(path = dirKUHR, pattern = glob2rx("2016*.dsv"))),  sep = ";")
data2017 <- read.csv(file.path(dirKUHR, list.files(path = dirKUHR, pattern = glob2rx("2017*.dsv"))),  sep = ";")
data2018 <- read.csv(file.path(dirKUHR, list.files(path = dirKUHR, pattern = glob2rx("2018*.dsv"))),  sep = ";")
data2019 <- read.csv(file.path(dirKUHR, list.files(path = dirKUHR, pattern = glob2rx("2019*.dsv"))),  sep = ";")
data2020 <- read.csv(file.path(dirKUHR, list.files(path = dirKUHR, pattern = glob2rx("2020*.dsv"))),  sep = ";")
data2021 <- read.csv(file.path(dirKUHR, list.files(path = dirKUHR, pattern = glob2rx("2021*.dsv"))),  sep = ";")
data2022 <- read.csv(file.path(dirKUHR, list.files(path = dirKUHR, pattern = glob2rx("2022*.dsv"))),  sep = ";")
data2023 <- read.csv(file.path(dirKUHR, list.files(path = dirKUHR, pattern = glob2rx("2023*.dsv"))),  sep = ";")


# Add a column for year of record -----------------------------------------
data2006$Year = 2006
data2007$Year = 2007
data2008$Year = 2008
data2009$Year = 2009
data2010$Year = 2010
data2011$Year = 2011
data2012$Year = 2012
data2013$Year = 2013
data2014$Year = 2014
data2015$Year = 2015
data2016$Year = 2016
data2017$Year = 2017
data2018$Year = 2018
data2019$Year = 2019
data2020$Year = 2020
data2021$Year = 2021
data2022$Year = 2022
data2023$Year = 2023


# Vertically concatenate tables into a HUGE table --------------------------
dataAll <- rbind(data2006, data2007, data2008, data2009, data2010, data2011, 
                 data2012, data2013, data2014, data2015, data2016, data2017,
                 data2018, data2019, data2020, data2021, data2022, data2023)


# Clear up existing tables
rm(list = ls(pattern="^data2"))


# Renaming column names ----------------------------------------------------
colnames(dataAll)[colnames(dataAll) == "PASIENTLOPENUMMER"]     <- "SerialNumber"
colnames(dataAll)[colnames(dataAll) == "FAGOMRAADE_KODE"]       <- "SubjectAreaCode"
colnames(dataAll)[colnames(dataAll) == "PRAKSISTYPE"]           <- "PracticeType"
colnames(dataAll)[colnames(dataAll) == "DIFFERANSEDAGER"]       <- "DiffDays"

colnames(dataAll)[colnames(dataAll) == "PASIENTKOMMUNE"]        <- "Municipality"
colnames(dataAll)[colnames(dataAll) == "PASIENTFYLKE"]          <- "County"
colnames(dataAll)[colnames(dataAll) == "BEHANDLER_KOMMUNENR"]   <- "TreatmentMunicipality"
colnames(dataAll)[colnames(dataAll) == "BEHANDLERFYLKE"]        <- "TreatmentCounty"

colnames(dataAll)[colnames(dataAll) == "KONTAKTTYPE"]           <- "ContactType"
colnames(dataAll)[colnames(dataAll) == "DIAGNOSER"]             <- "Diagnosis"
colnames(dataAll)[colnames(dataAll) == "DIAGNOSEKODEVERK"]      <- "DiagnosticManual"
colnames(dataAll)[colnames(dataAll) == "Year"]                  <- "Year"


# Save merged KUHR --------------------------------------------------------
outputName <- paste(as.character(Sys.Date()), "-MoBa-MergedKUHR.csv", sep = "")
write.table(dataAll, file = file.path(dirWorking, outputName), sep = ",", row.names = F)