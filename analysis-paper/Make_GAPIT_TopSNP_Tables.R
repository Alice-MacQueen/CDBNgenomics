#!/usr/bin/env Rscript

source(file.path("~", "Alice", "bin", "R_Functions",
                 "Functions_JLS_plots_2019-03.R"))

# ----------- PC0 Setup -------------------------------

setwd(file.path("~", "Alice", "2016-2017_GBS", "8_Tassel", "GAPIT_BLUP_PC0"))

gapitpath = paste0(file.path("~", "Alice", "2016-2017_GBS", "8_Tassel",
                             "GAPIT_BLUP_PC0"),"/")
PC0 <- c("BM", "CB", "Det_scr", "Earliest_Year_CDBN", "GHmode", "CT", "DF",
         "DM", "HB", "PH", "RU", "SA", "SF")
phenos_PC0 <- c("Biomass (kg)", "CBB score", "Determinacy",
                "First Year in CDBN", "Growth habit", "CTV presence/absence",
                "DF (days)", "DM (days)", "HB score",  "PH (cm)", "RU score",
                "SA score", "SF (days)")


## ----------- Make Top SNP DataFrames ---------------
phe_cutoff_tables <- list()

for(i in seq_along(PC0)){
  GWAS_obj <- load_GAPIT_GWAS_stats(path = gapitpath,
                                    phenotype = PC0[i])
  phe_cutoff_tables[[i]] <- gapit_table_topsnps(GWAS_obj = GWAS_obj)
}

names(phe_cutoff_tables) <- PC0
allphe_top12SNPs_within0bp <- gapit_arrange_snp_tables(phe_cutoff_tables,
                                                       "top12SNPs_within0bp")
allphe_top12SNPs_within100kbp <- gapit_arrange_snp_tables(phe_cutoff_tables,
                                                          "top12SNPs_within1e+05bp")
allphe_FDR10per_within100kbp <- gapit_arrange_snp_tables(phe_cutoff_tables,
                                                         "FDR0.1_within1e+05bp")
allphe_FDR10per_within0bp <- gapit_arrange_snp_tables(phe_cutoff_tables,
                                                      "FDR0.1_within0bp")

gapit_write_snp_tables_from_list(list = allphe_top12SNPs_within0bp,
                                 file = paste0("Top12SNPS_withAnnotatedSNP",
                                               "_GAPIT_BLUP_PC0_phenotypedonly",
                                               "_", Sys.Date(), ".xlsx"))
saveRDS(allphe_FDR10per_within0bp,
        file = paste0("10perFDR_withAnnotatedSNP",
                      "_GAPIT_BLUP_PC0_phenotypedonly",
                      "_", Sys.Date(), ".xlsx"))
saveRDS(allphe_top12SNPs_within100kbp,
        file = paste0("Top12SNPS_annosWithin100kbp",
                      "_GAPIT_BLUP_PC0_phenotypedonly",
                      "_", Sys.Date(), ".rds"))
saveRDS(allphe_FDR10per_within100kbp,
        file = paste0("10perFDR_annosWithin100kbp",
                      "_GAPIT_BLUP_PC0_phenotypedonly",
                      "_", Sys.Date(), ".rds"))

# ----------- PC1 Setup -------------------------------

setwd(file.path("~", "Alice", "2016-2017_GBS", "8_Tassel", "GAPIT_BLUP_PC1"))

gapitpath = paste0(file.path("~", "Alice", "2016-2017_GBS", "8_Tassel",
                             "GAPIT_BLUP_PC1"),"/")
PC1 <- c("RR")
phenos_PC1 <- c("Root rot score")


## ----------- Make Top SNP DataFrames ---------------
phe_cutoff_tables <- list()

for(i in seq_along(PC1)){
  GWAS_obj <- load_GAPIT_GWAS_stats(path = gapitpath,
                                    phenotype = PC1[i])
  phe_cutoff_tables[[i]] <- gapit_table_topsnps(GWAS_obj = GWAS_obj)
}

names(phe_cutoff_tables) <- PC1
allphe_top12SNPs_within0bp <- gapit_arrange_snp_tables(phe_cutoff_tables,
                                                       "top12SNPs_within0bp")
allphe_top12SNPs_within100kbp <- gapit_arrange_snp_tables(phe_cutoff_tables,
                                                          "top12SNPs_within1e+05bp")
allphe_FDR10per_within100kbp <- gapit_arrange_snp_tables(phe_cutoff_tables,
                                                         "FDR0.1_within1e+05bp")
allphe_FDR10per_within0bp <- gapit_arrange_snp_tables(phe_cutoff_tables,
                                                      "FDR0.1_within0bp")

gapit_write_snp_tables_from_list(list = allphe_top12SNPs_within0bp,
                                 file = paste0("Top12SNPS_withAnnotatedSNP",
                                               "_GAPIT_BLUP_PC1_phenotypedonly",
                                               "_", Sys.Date(), ".xlsx"))
saveRDS(allphe_FDR10per_within0bp,
        file = paste0("10perFDR_withAnnotatedSNP",
                      "_GAPIT_BLUP_PC1_phenotypedonly",
                      "_", Sys.Date(), ".rds"))
saveRDS(allphe_top12SNPs_within100kbp,
        file = paste0("Top12SNPS_annosWithin100kbp",
                      "_GAPIT_BLUP_PC1_phenotypedonly",
                      "_", Sys.Date(), ".rds"))
saveRDS(allphe_FDR10per_within100kbp,
        file = paste0("10perFDR_annosWithin100kbp",
                      "_GAPIT_BLUP_PC1_phenotypedonly",
                      "_", Sys.Date(), ".rds"))

# ----------- PC2 Setup -------------------------------

setwd(file.path("~", "Alice", "2016-2017_GBS", "8_Tassel", "GAPIT_BLUP_PC2"))

gapitpath = paste0(file.path("~", "Alice", "2016-2017_GBS", "8_Tassel",
                             "GAPIT_BLUP_PC2"),"/")
PC2 <- c("EV", "GH", "HI", "LG", "SY", "WM", "ZN", "SW", "CB", "HB")
phenos_PC2 <- c("EV score", "GH type", "HI (%)", "LG score", "SY (kg/ha)",
                "WM score", "ZN score", "Seed weight (mg)", "CBB score",
                "HB score")



## ----------- Make Top SNP DataFrames ---------------
phe_cutoff_tables <- list()

for(i in seq_along(PC2)){
  GWAS_obj <- load_GAPIT_GWAS_stats(path = gapitpath,
                                    phenotype = PC2[i])
  phe_cutoff_tables[[i]] <- gapit_table_topsnps(GWAS_obj = GWAS_obj)
}

names(phe_cutoff_tables) <- PC2
allphe_top12SNPs_within0bp <- gapit_arrange_snp_tables(phe_cutoff_tables,
                                                       "top12SNPs_within0bp")
allphe_top12SNPs_within100kbp <- gapit_arrange_snp_tables(phe_cutoff_tables,
                                                          "top12SNPs_within1e+05bp")
allphe_FDR10per_within100kbp <- gapit_arrange_snp_tables(phe_cutoff_tables,
                                                         "FDR0.1_within1e+05bp")
allphe_FDR10per_within0bp <- gapit_arrange_snp_tables(phe_cutoff_tables,
                                                      "FDR0.1_within0bp")

gapit_write_snp_tables_from_list(list = allphe_top12SNPs_within0bp,
                                 file = paste0("Top12SNPS_withAnnotatedSNP",
                                               "_GAPIT_BLUP_PC2_phenotypedonly",
                                               "_", Sys.Date(), ".xlsx"))
saveRDS(allphe_FDR10per_within0bp,
        file = paste0("10perFDR_withAnnotatedSNP",
                      "_GAPIT_BLUP_PC2_phenotypedonly",
                      "_", Sys.Date(), ".rds"))
saveRDS(allphe_top12SNPs_within100kbp,
        file = paste0("Top12SNPS_annosWithin100kbp",
                      "_GAPIT_BLUP_PC2_phenotypedonly",
                      "_", Sys.Date(), ".rds"))
saveRDS(allphe_FDR10per_within100kbp,
        file = paste0("10perFDR_annosWithin100kbp",
                      "_GAPIT_BLUP_PC2_phenotypedonly",
                      "_", Sys.Date(), ".rds"))

# ----------- PC3 Setup -------------------------------

setwd(file.path("~", "Alice", "2016-2017_GBS", "8_Tassel", "GAPIT_BLUP_PC3"))

gapitpath = paste0(file.path("~", "Alice", "2016-2017_GBS", "8_Tassel",
                             "GAPIT_BLUP_PC3"),"/")
PC3 <- c("BR", "PH", "Earliest_Year_CDBN")
phenos_PC3 <- c("BR (BCMV response)", "Plant height (cm)",
                "First Year in the CDBN")


## ----------- Make Top SNP DataFrames ---------------
phe_cutoff_tables <- list()

for(i in seq_along(PC3)){
  GWAS_obj <- load_GAPIT_GWAS_stats(path = gapitpath,
                                    phenotype = PC3[i])
  phe_cutoff_tables[[i]] <- gapit_table_topsnps(GWAS_obj = GWAS_obj)
}

names(phe_cutoff_tables) <- PC3
allphe_top12SNPs_within0bp <-
  gapit_arrange_snp_tables(phe_cutoff_tables, "top12SNPs_within0bp")
allphe_top12SNPs_within100kbp <-
  gapit_arrange_snp_tables(phe_cutoff_tables, "top12SNPs_within1e+05bp")
allphe_FDR10per_within100kbp <-
  gapit_arrange_snp_tables(phe_cutoff_tables, "FDR0.1_within1e+05bp")
allphe_FDR10per_within0bp <-
  gapit_arrange_snp_tables(phe_cutoff_tables, "FDR0.1_within0bp")

gapit_write_snp_tables_from_list(list = allphe_top12SNPs_within0bp,
                                 file = paste0("Top12SNPS_withAnnotatedSNP",
                                               "_GAPIT_BLUP_PC3_phenotypedonly",
                                               "_", Sys.Date(), ".xlsx"))
saveRDS(allphe_FDR10per_within0bp,
        file = paste0("10perFDR_withAnnotatedSNP",
                      "_GAPIT_BLUP_PC3_phenotypedonly",
                      "_", Sys.Date(), ".rds"))
saveRDS(allphe_top12SNPs_within100kbp,
        file = paste0("Top12SNPS_annosWithin100kbp",
                      "_GAPIT_BLUP_PC3_phenotypedonly",
                      "_", Sys.Date(), ".rds"))
saveRDS(allphe_FDR10per_within100kbp,
        file = paste0("10perFDR_annosWithin100kbp",
                      "_GAPIT_BLUP_PC3_phenotypedonly",
                      "_", Sys.Date(), ".rds"))

# ----------- PC10 Setup -------------------------------

setwd(file.path("~", "Alice", "2016-2017_GBS", "8_Tassel", "GAPIT_BLUP_PC10"))

gapitpath = paste0(file.path("~", "Alice", "2016-2017_GBS", "8_Tassel",
                             "GAPIT_BLUP_PC10"),"/")
PC10 <- c("RU", "SA", "GH", "SW", "SY", "Earliest_Year_CDBN", "GH_mes",
          "GH_dur", "GH_MA", "SY_mes", "SY_dur", "SY_MA", "SW_mes", "SW_dur",
          "SW_MA", "GH_det", "GH_ind", "SY_det", "SY_ind", "SW_det", "SW_ind",
          "GH_MA_ind", "SY_MA_ind", "SW_MA_ind")
phenos_PC10 <- c("Rust score", "Seed appearance score", "Growth habit",
                 "Seed weight (mg)", "Seed yield (kg/ha)",
                 "First Year in the CDBN", "Growth habit (Mesoamerican race)",
                 "Growth habit (Durango race)",
                 "Growth habit (Mesoamerican gene pool)",
                 "Seed yield (Mesoamerican race)",
                 "Seed yield (Durango race)",
                 "Seed yield (Mesoamerican gene pool)",
                 "Seed weight (Mesoamerican race)",
                 "Seed weight (Durango race)",
                 "Seed weight (Mesoamerican gene pool)",
                 "Growth habit (determinates)",
                 "Growth habit (indeterminates)",
                 "Seed yield (determinates)",
                 "Seed yield (indeterminates)",
                 "Seed weight (determinates)",
                 "Seed weight (indeterminates)",
                 "Growth habit (MA indeterminates)",
                 "Growth habit (MA indeterminates)",
                 "Growth habit (MA indeterminates)"
)


## ----------- Make Top SNP DataFrames ---------------
phe_cutoff_tables <- list()

for(i in seq_along(PC10)){
  GWAS_obj <- load_GAPIT_GWAS_stats(path = gapitpath,
                                    phenotype = PC10[i])
  phe_cutoff_tables[[i]] <- gapit_table_topsnps(GWAS_obj = GWAS_obj)
}

names(phe_cutoff_tables) <- PC10
allphe_top12SNPs_within0bp <- gapit_arrange_snp_tables(phe_cutoff_tables,
                                                       "top12SNPs_within0bp")
allphe_top12SNPs_within100kbp <- gapit_arrange_snp_tables(phe_cutoff_tables,
                                                          "top12SNPs_within1e+05bp")
allphe_FDR10per_within100kbp <- gapit_arrange_snp_tables(phe_cutoff_tables,
                                                         "FDR0.1_within1e+05bp")
allphe_FDR10per_within0bp <- gapit_arrange_snp_tables(phe_cutoff_tables,
                                                      "FDR0.1_within0bp")

gapit_write_snp_tables_from_list(list = allphe_top12SNPs_within0bp,
                                 file = paste0("Top12SNPS_withAnnotatedSNP",
                                               "_GAPIT_BLUP_PC10_phenotypedonly",
                                               "_", Sys.Date(), ".xlsx"))
saveRDS(allphe_FDR10per_within0bp,
        file = paste0("10perFDR_withAnnotatedSNP",
                      "_GAPIT_BLUP_PC10_phenotypedonly",
                      "_", Sys.Date(), ".rds"))
saveRDS(allphe_top12SNPs_within100kbp,
        file = paste0("Top12SNPS_annosWithin100kbp",
                      "_GAPIT_BLUP_PC10_phenotypedonly",
                      "_", Sys.Date(), ".rds"))
saveRDS(allphe_FDR10per_within100kbp,
        file = paste0("10perFDR_annosWithin100kbp",
                      "_GAPIT_BLUP_PC10_phenotypedonly",
                      "_", Sys.Date(), ".rds"))

## ----------- For PC10, Make plots with TopSNP Effects ---------------
for(i in seq_along(PC10)){
  GWAS_obj <- load_GAPIT_GWAS_stats(path = gapitpath,
                                    phenotype = PC10[i])

  panel_a <- gapit_plot_ggman_labeled(GWAS_obj = GWAS_obj)
  panel_b <- gapit_plot_qq(GWAS_obj = GWAS_obj)
  panel_c <- gapit_plot_resid_density(plotpheno = phenos_PC10[i],
                                      GWAS_obj = GWAS_obj)
  row1 <- plot_grid(panel_a, panel_b, panel_c, align = "hv", nrow = 1,
                    rel_widths = c(1/2, 1/4, 1/4), labels = "AUTO")
  ggsave(paste0("PC10_kinpheoBLUP", PC10[i], "Manhattan.png"), width = 9,
         height = 2.75*.9, units = "in", dpi = 400)

  panels_row2to3 <- gapit_plot_topsnp_effects(plotpheno = phenos_PC10[i],
                                              GWAS_obj = GWAS_obj)
  row2to3 <- plot_grid(panels_row2to3[[1]], panels_row2to3[[2]],
                       panels_row2to3[[3]], panels_row2to3[[4]], nrow = 2,
                       axis = "lb", rel_widths = c(2/3, 1/3),
                       labels = c("D", "E", "F", "G", "H", "I"))
  ggsave(row2to3, file = paste0("PC10_kinpheoBLUP", PC10[i],
                                "TopSNPEffects.png"),
         width = 9, height = 2.75*2*.9, units = "in", dpi = 400)

  plot_grid(row1, row2to3, axis = "lr", nrow = 2, rel_heights = c(1/3, 2/3))
  ggsave(paste0("PC10_kinpheoBLUP", PC10[i], "Manhattan_and_TopSNPEffects.png"),
         width = 10*.9, height = 2.75*3*.9, units = "in", dpi = 400)
}
