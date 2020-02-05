#!/usr/bin/env Rscript
#args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
#if (length(args)==0) {
#  stop("At least one argument must be supplied (input phenotype file number).n", call.=FALSE)
#} #else if (length(args)==1) {
# default output file
#args[2] = "out.txt"
#}

require(tidyverse)
require(rlang)
# source("/home/pubjuenger/Alice/Github/Functions_ggplot-theme-adjustments_2018-01-03.R")
# theme_set(theme_oeco)
# source("/home/alice/CDBN/bin/R_Functions/Functions_VCF_to_GAPIT.R")
#library(multtest)
#library(gplots)
#library(LDheatmap)
library(genetics)
library(ape)
library(EMMREML)
library(compiler) #this library is already installed in R
library("scatterplot3d")
source("/home/alice/CDBN/bin/R_Functions/GAPIT_functions.R")
source("/home/alice/CDBN/bin/R_Functions/EMMA_functions.R")

## Use GAPIT to run GWAS on the BLUPs for each phenotype

#Step 1: Import and format phenotypic data; create phenotype directory
projectdir <- "/home/alice/Github/CDBNgenomics/data-raw/ASReml/GAPIT-7-PCs/"
setwd(projectdir)

gapit_phe_asreml <- read.table(file = "/home/alice/Github/CDBNgenomics/data-raw/ASReml/ASReml_coefficients_13_GAPIT_2020-01-12.txt", head = TRUE)
gapit_phe_asreml <- gapit_phe_asreml %>%
  dplyr::select(Taxa, SW)
# Run GAPIT on this LbY
myGAPIT <- GAPIT(
  Y = gapit_phe_asreml,
  PCA.total = 7,
  model = "CMLM",
  file.GD="Numerical_format_GD_CDBN_001_359_pedigree_fillin_chr",
  file.Ext.GD="txt",
  file.GM="Numerical_format_GM_CDBN_001_359_pedigree_fillin_chr",
  file.Ext.GM="txt",
  file.from = 1,
  file.to = 11,
  file.path = "/home/alice/Github/CDBNgenomics/data-raw/GAPIT_Numerical_format_files/"
)

# other options for GAPIT:
# Model.selection = TRUE
# SNP.fraction = 0.4
# args(GAPIT) # see function arguments: default for model is "MLM"
