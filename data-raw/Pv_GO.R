## code to prepare `Pv_GO` dataset goes here
library(tidyverse)

Pv_GO <- read.table(file = file.path("I:", "OneDrive", "NSF Common Bean",
                                     "CDBN Genomics", "Pvulgaris-genome",
                                     "annotation", "Pvulgaris", "v2.1",
                                     "annotation", "commonbean_go.txt"),
                    sep = "\t", header = FALSE)


names(Pv_GO) <- c("transcriptName", "GOCategories", "GOInfo")

Pv_GO <- Pv_GO %>%
  separate(transcriptName, into = c("GENEID", "variant"), sep = -3) %>%
  mutate(GENEID = as.character(GENEID))

usethis::use_data(Pv_GO, overwrite = TRUE)
