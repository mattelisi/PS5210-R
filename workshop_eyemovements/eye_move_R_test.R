rm(list=ls())

# # relies on SR Research’s EDF API being installed (comes with the EyeLink Developer’s Kit / EDF API
# # install.packages("eyelinkReader")
# 
# library(eyelinkReader)
# 
# rec <- read_edf("S1.edf", import_samples = TRUE)  # will error if EDF API isn't installed
# str(rec)


######### alternative without API

edf <- "S1.edf"

# Basic: creates subj01.edf.asc (or similar depending on platform/version)
system2("edf2asc", args = c("-y", edf))

# Then import the ASC
# install.packages("eyelinker")
library(eyelinker)
dat <- read.asc("S1.asc", parse_all=TRUE)
str(dat)

grepl("TrialData", dat$msg$text)

dat$msg$text[grepl("TrialData", dat$msg$text)]
