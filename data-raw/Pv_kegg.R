## code to prepare `Pv_kegg` dataset goes here
library(tidyverse)
Pv_kegg <- read.table(file = file.path("I:", "OneDrive", "NSF Common Bean",
                                       "CDBN Genomics", "Pvulgaris-genome",
                                       "annotation", "Pvulgaris", "v2.1",
                                       "annotation", "commonbean_kegg.txt"),
                      sep = "\t", header = FALSE)

names(Pv_kegg) <- c("transcriptName", "KEGG", "KEGGInfo")

Pv_kegg <- Pv_kegg %>%
  separate(transcriptName, into = c("GENEID", "variant"), sep = -3)

usethis::use_data(Pv_kegg)
