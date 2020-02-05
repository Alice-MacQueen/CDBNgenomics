library(tidyverse)
library(CDBNgenomics)


mash_output <- mash_standard_run(path = "/home/alice/Github/CDBNgenomics/data-raw/ASReml/GAPIT-for-mash", numSNPs = 2200, saveoutput = TRUE)
