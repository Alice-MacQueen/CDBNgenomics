require(tidyverse)
require(mashr)
source(file.path("~", "Alice", "2016-2017_GBS", "8_Tassel", "for_paper_2019",
                 "Functions_mash_standard_run.R"))

mash_standard_run(path = file.path("~", "Alice", "2016-2017_GBS", "8_Tassel",
                                   "for_paper_2019"), numSNPs = 4000,
                  saveoutput = TRUE, U_ed = NA)
