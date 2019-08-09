## code to prepare `BLUPs` dataset goes here
library(tidyverse)
library(rrBLUP)


METG <- readRDS(file.path("data-raw", "ignore", "CDBN_Phenotypes_V2.0.rds"))
METG <- METG %>%
  dplyr::select(Genotype:Longitude, SY, SW, DM, DF, SF, LG, HI, PH, BM, GH,
                SA, CB, RU, EV, WM, VN, CT, HB, BR, CM, RR, ZN) %>%
  mutate(VN = ifelse(VN == "a",
                     0,
                     VN),
         VN = ifelse(VN == "b",
                     1,
                     VN),
         VN = ifelse(VN >= 1,
                     1,
                     VN),
         VN = as.numeric(VN),
         RR = ifelse(RR > 9,
                     9,
                     RR)
  )

phebyloc <- readRDS(file = file.path("C:", "Users", "alice", "OneDrive",
                                     "NSF Common Bean", "R", "CDBN_mashr",
                                     paste0("Phenotype_by_location_dataframes_",
                                            "for-mashr-GWAS_312-and-327.rds")))
Seq_small327 <- phebyloc[[11]]

METG_Seq327 <- METG %>%
  dplyr::select(-Seq_ID, -Gene_pool, -Market_class_ahm, -Race) %>%
  left_join(Seq_small327) %>%
  dplyr::select(CDBN_ID, Seq_ID, Gene_pool, Race, Market_class_ahm, Det_scr,
                GHmode, everything())

V_g <- read.table(file.path("..", "..", "CDBN_GxE_Analysis-1_Pheno-Ave",
                            paste0("Filter_120Ct_MAF1per_CDBN_37yr_pedigree",
                              "_imputed_names_kinship.txt")),
             skip = 3, sep = "\t", row.names = 1)
colnames(V_g) <- rownames(V_g)
Vg <- as.matrix(V_g)
PHE <- list(quote(BM), quote(BR), quote(CB), quote(CM), quote(CT), quote(DF),
            quote(DM), quote(EV), quote(GH), quote(HB), quote(HI), quote(LG),
            quote(PH), quote(RR), quote(RU), quote(SA), quote(SF), quote(SW),
            quote(SY), quote(WM), quote(ZN))

rrkinblupout <- list()
kinseqblup <- list()

for(i in seq_along(PHE)){

  # Run rrBLUP to generate BLUPs
  testdf <- METG_Seq327 %>%
    filter(!is.na(eval(PHE[[i]])) & !is.na(Seq_ID)) %>%
    mutate(LbY = paste(Location_code, Year, sep = "_"))
  testrrk <- kin.blup(data = testdf, geno = "Taxa", pheno = deparse(PHE[[i]]),
                      K = Vg, fixed = c("Location_code", "LbY"))
  saveRDS(testrrk, paste0("data-raw/ignore/", PHE[[i]], "_BLUPs_kinship.rds"))

  # Convert BLUPs (the $g component from kin.blup) to data frames
  rrkinblupout[[i]] <- enframe(testrrk$g) %>%
    rename(Taxa = name) %>%
    left_join(Seq_small327)
  names(rrkinblupout[[i]])[2] <- deparse(PHE[[i]])

  # Which Taxa were phenotyped for this phenotype? Just keep those for GWAS.
  taxa_phenotyped <- testdf %>%
    dplyr::select(Taxa, !! PHE[[i]]) %>%
    as_tibble() %>%
    group_by(Taxa) %>%
    summarise(count = n())
  kinseqblup[[i]] <- taxa_phenotyped %>%
    left_join(rrkinblupout[[i]])
}
# Save BLUPs list

#saveRDS(rrkinblupout,
#        file = "CDBN327_21phenotype_BLUPs_with_kinship_rrBLUP_2019-08-06.rds")
#saveRDS(kinseqblup,
#        file = paste0("CDBN327_21phenotype_BLUPs_with_kinship_phenotyped_",
#                      "only_rrBLUP_2019-08-06.rds"))
#
#kinseqblup <- readRDS(file = file.path("analysis-raw",
#                                       paste0("CDBN327_21phenotype_BLUPs_with_",
#                                              "kinship_phenotyped_only_rrBLUP_",
#                                              "2019-08-06.rds")))

BLUPs <- Seq_small327 %>%
  left_join(dplyr::select(kinseqblup[[1]], Taxa, 3))

for(i in seq_along(kinseqblup)[-1]){
  BLUPs <- BLUPs %>%
    left_join(dplyr::select(kinseqblup[[i]], Taxa, 3))
}
colSums(!is.na(BLUPs))

usethis::use_data(BLUPs, overwrite = TRUE)
