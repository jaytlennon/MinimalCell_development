################################################################################
#                                                                              #
#	Lennon Lab Growth Curve Analysis Example                                     #
#   Parameter Estimate Code                                                    #
#                                                                              #
################################################################################
#                                                                              #
#	Written by: M. Muscarella                                                    #
#	  Last update: 11/24/2015 by M. Muscarella & V. Kuo                          #
#                                                                              #
################################################################################

# Setup Work Environment
rm(list=ls())
setwd("~/Github/MinimalCell/scripts/6.GrowthRates/test")

# Load Dependencies
#source("../bin/modified_Gomp.R")
#source("~/Github/MinimalCell/scripts/6.GrowthRates/bin/modified_Gomp.R")
source("../bin/modified_Gomp.R")

# Examples #

# 1. Synergy = F
##  The input here should be a csv or other delimited file or an R object
## The data must have the following format: Header followed by whole or decimal numbers 
## Time  A1  A2  A3  ...
## 30    0.2 0.1 0.4 ...
GC <- read.csv("../data/MinimalCell_try2.csv")
head(GC, header = T)

growth.modGomp(input = GC, output.name = "GC_test", 
               synergy = F, temp = F, smooth = T, trim = T)


## Let the function import the CSV file
growth.modGomp("../data/MinimalCell_try2.csv", "GC_test", synergy=F, temp=F)


# 2. Synergy = T



# Run Example (Synergy)
growth.modGomp("../data/MinimalCell_try2.csv", "GC_test", skip = 38)


# These are very bad fits
growth.modGomp(input = "../data/GrowthCurve_Example2.txt", output.name = "test", skip=48,
               synergy=T, temp=F, smooth=F, trim = F)

growth.modGomp("../data/GrowthCurve_Example2.txt", "test", skip=48)
growth.modGomp("../data/RPF.txt", "RPF", skip = 39)
growth.modGomp("../test/125_1_6_GrowthCurves.txt", "BehriTest", skip = 38)

# Run Example (CSV file)
growth.modGomp("../data/MinimalCell_try2.csv", "GC_test", synergy=F, temp=F)
growth.modGomp("../data/Pseudo.csv", "csv_test", synergy=F, temp=F, smooth=F)
