#!/usr/bin/env Rscript

source(file.path("~", "CDBN", "bin", "R_Functions",
                 "Functions_JLS_plots_2020-01.R"))
library(CDBNgenomics)

# -----------  Setup -------------------------------

setwd(file.path("~", "Github", "CDBNgenomics", "data-raw", "ASReml",
                "GAPIT-for-mash"))

gapitpath = paste0(file.path("~", "Github", "CDBNgenomics", "data-raw", "ASReml",
                             "GAPIT-for-mash"),"/")
PC0 <- gapit_results_in_filepath(path = ".") %>%
  str_sub(start = 6)

phenos_PC0 <- c("Biomass (kg)","Blackroot response", "CBB score",
                "BCMV damage score", "CTV presence/absence",
                "Days to Flowering (days)", "Days to Maturity (days)",
                "First Year in CDBN", "Early vigor score", "Growth habit",
                "Halo blight score", "Harvest index (%)", "Lodging score",
                "Plant height (cm)","Root rot damage score",
                "Rust damage score",
                "Seed appearance score", "Seedfill duration (days)",
                "Seed weight (mg)", "Seed yield (kg/ha)",
                "White mold damage score", "Zinc deficiency damage score")
#1] "CMLM.BM"                 "CMLM.BR"
#3] "CMLM.CB"                 "CMLM.CM"
#5] "CMLM.CT"                 "CMLM.DF"
#7] "CMLM.DM"                 "CMLM.Earliest_Year_CDBN"
#9] "CMLM.EV"                 "CMLM.GH"
#11] "CMLM.HB"                 "CMLM.HI"
#13] "CMLM.LG"                 "CMLM.PH"
#15] "CMLM.RR"                 "CMLM.RU"
#17] "CMLM.SA"                 "CMLM.SF"
#19] "CMLM.SW"                 "CMLM.SY"
#21] "CMLM.WM"                 "CMLM.ZN"


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
