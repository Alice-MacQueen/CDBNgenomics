## code to prepare `txdb` dataset goes here
library(AnnotationDbi)
txdb <- loadDb(file = file.path("I:", "OneDrive", "NSF Common Bean",
                                "CDBN Genomics", "Pvulgaris-genome",
                                "annotation", "Pvulgaris", "v2.1",
                                "annotation", "Pvulgaris_442_v2.1.gene.sqlite"))

usethis::use_data(txdb)
