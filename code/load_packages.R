# ============================================= #
# script: load_packages.R
# Project: Fiber Intervention Study
# Author(s): L. Greathouse et al.
# ============================================= #
# Date Created: 2019-12-10
# Date Modified: 2019-12-10
# By: R. Noah Padgett
# ============================================= #
# ============================================= #
# Purpose:
# This R script is for loading all necessary
#   R packages
#
# No output - just loading packages into the
#   environment
# ============================================= #
# Set up directory and libraries
rm(list=ls())

# BiocManager "bioconductor manager"
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
# check if phyloseq is installed:
if(!requireNamespace("phyloseq", quietly=T)){
  BiocManager::install("phyloseq")
}


# list of packages
packages <- c("phyloseq","vegan", "lme4", "lmerTest",
              "tidyverse", "readr", "readxl", "forcats",
              "data.table", "plyr", "ggplot2",
              "kableExtra", "xtable", "gridExtra", "viridis",
              "patchwork", "gvlma", "car", "dplyr")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, quiet = T, dependencies = T)
# Load packages
lapply(packages, library, character.only = TRUE)

w.d <- getwd()

