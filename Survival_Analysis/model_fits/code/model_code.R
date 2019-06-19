
# This document contains the model code to create the model fits for the event history analysis book

# note: model results are different from CSSS 544 final project, fixed an error due to improper implementation
# of left-truncation correction 

#### Libs ####
library(tidyverse)
library(survival)
library(maxLik)

#### Load data ####
dat <- readRDS(here("model_fits", "survdat.rds"))

# set data up for maxLik package
x <- dat


