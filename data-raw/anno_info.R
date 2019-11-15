## code to prepare `anno_info` dataset goes here
anno_info <- read.table(file = file.path("I:", "OneDrive", "NSF Common Bean",
                                         "CDBN Genomics", "Pvulgaris-genome",
                                         "annotation", "Pvulgaris", "v2.1",
                                         "annotation",
                                         "Pvulgaris_442_v2.1.annotation_info.txt"),
                        sep = "\t", fill = TRUE, header = TRUE)
usethis::use_data(anno_info)
