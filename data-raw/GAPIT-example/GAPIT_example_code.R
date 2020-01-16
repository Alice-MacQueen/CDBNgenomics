require(tidyverse)
require(rlang)
#library(gplots)
#library(LDheatmap)
library(genetics)
library(ape)
library(EMMREML)
library(compiler) #this library is already installed in R
library("scatterplot3d")
source("/home/alice/Github/CDBNgenomics/data-raw/GAPIT-example/GAPIT_functions.R")
source("/home/alice/Github/CDBNgenomics/data-raw/GAPIT-example/EMMA_functions.R")

## Use GAPIT to run GWAS on the BLUPs for each phenotype

#Step 1: Import and format phenotypic data; create phenotype directory
projectdir <- "/home/alice/Github/CDBNgenomics/data-raw/GAPIT-example/"
setwd(projectdir)

gapit_phe_asreml <- read.table(file = "/home/alice/Github/CDBNgenomics/data-raw/GAPIT-example/kinBLUP_phenotypedOnly_GAPIT_PC0_2019-02-27.txt", head = TRUE)
# Run GAPIT on this LbY
myGAPIT <- GAPIT(
  Y = gapit_phe_asreml,
  PCA.total = 0,
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
