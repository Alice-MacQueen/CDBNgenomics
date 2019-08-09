## code to prepare `metadata` dataset goes here

metadata <- readRDS(file.path("C:", "Users", "ahm543", "OneDrive",
                              "NSF Common Bean", "R",
                              "CDBN_GxE_Analysis-1_Pheno-Ave",
                              "Metadata_for_CDBN_327_2019-02-27.rds"))

usethis::use_data(metadata)
