library(tidyverse)
library(asreml)

METG <- readRDS(file.path("data-raw", "ignore", "CDBN_Phenotypes_V2.1.rds"))
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

phebyloc <- readRDS(file = file.path("data-raw", "ignore",
                                     paste0("Phenotype_by_location_dataframes_",
                                            "for-mashr-GWAS_312-and-327.rds")))
Seq_small327 <- phebyloc[[11]]

METG_Seq327 <- METG %>%
  dplyr::select(-Seq_ID, -Gene_pool, -Market_class_ahm, -Race) %>%
  left_join(Seq_small327) %>%
  dplyr::select(CDBN_ID, Seq_ID, Gene_pool, Race, Market_class_ahm, Det_scr,
                GHmode, everything())

V_g <- read.table(file.path("data-raw", "ignore",
                            paste0("Filter_120Ct_MAF1per_CDBN_37yr_pedigree",
                                   "_imputed_names_kinship.txt")),
                  skip = 3, sep = "\t", row.names = 1)
colnames(V_g) <- rownames(V_g)
Vg <- as.matrix(V_g)
# Seq_Vg <- which(!(colnames(Vg) %in% levels(METG_Seq327$Taxa)))
# Vg327 <- Vg[Seq_Vg, Seq_Vg]
Vg2 <- Vg[-347,-347] # remove an entry in the matrix which gives a negative eigenvalue
# VgEig <- eigen(Vg327) # look at eigenvalues of this matrix - which are negative?
# which(VgEig$values < 0 )

PHEGXE <- list(quote(BM), quote(DF), quote(DM), quote(HI), quote(LG),
               quote(PH),quote(SF), quote(SW), quote(SY))


METG_Seq327 <- METG_Seq327 %>%
  mutate(Taxa = as_factor(Taxa),
         Location_code = as_factor(Location_code),
         Year = as_factor(Year))

for(i in seq_along(PHEGXE)){
  testdf <- METG_Seq327 %>%
    filter(!is.na(eval(PHEGXE[[i]])) & !is.na(Seq_ID)) %>%
    mutate(LbY = paste(Location_code, Year, sep = "_"),
           LbY = as_factor(LbY))

  asr_out <- asreml(eval(PHEGXE[[i]]) ~ Taxa, random = ~vm(Taxa, Vg2) +
                      ~idv(Location_code) +
                      ~idv(Location_code):idv(Year) +
                      ~idv(Taxa):idv(Location_code) +
                      ~idv(Taxa):idv(Year),
                    residual = ~idv(units), data = testdf, workspace = "6gb")

  print(PHEGXE[[i]])
  asr_out$loglik
  summary(asr_out)$varcomp
  saveRDS(asr_out, file = file.path("data-raw", "ignore",
                                    paste0(PHEGXE[[i]], "_Reviewer_Model_ASReml",
                                           str_replace_all(Sys.time(), ":",
                                                           "."), ".rds")))
  asr_out_pv <- predict(asr_out, classify = "Taxa")
  saveRDS(asr_out_pv, file = file.path("data-raw", "ignore",
                                       paste0(PHEGXE[[i]],
                                              "_Reviewer_Model_Predictions",
                                              str_replace_all(Sys.time(), ":",
                                                              "."), ".rds")))
}
