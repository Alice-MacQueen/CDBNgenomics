#!/usr/bin/Rscript
source("/home/pubjuenger/Alice/bin/R_Functions/Functions_JLS_plots.R")

# Now kinship
tasselpath = "METG_BLUP_Output/"
PHE1 <- c("BM", "BR", "CB", "CM", "CT", "DF", "DM",  "EV", "GH", "HB",
          "HI", "LG", "PH", "RR", "RU", "SA", "SF", "SW", "SY", "WM", "ZN") #
phenopath = "METG_BLUP_Phenos/"
PHEN <- paste0("METG_kinshipBLUPs_phenotypedOnly_CDBN327_", PHE1)
pheno_name_end <- "_for-Tassel_2019-02-25_v01.txt"
phenofile <- paste0(phenopath, PHEN, pheno_name_end)
# PHE <- list(quote(BM), quote(BR), quote(CB), quote(CM), quote(CT), quote(DF), quote(DM), quote(EV), quote(GH), quote(HB), quote(HI), quote(LG), quote(PH), quote(RR), quote(RU), quote(SA), quote(SF), quote(SW), quote(SY), quote(WM), quote(ZN))
gwasname = paste0("METG_kinshipBLUPs_phenotypedOnly_CDBN327_", PHE1)
METGPHEN <- c("Biomass (kg)", "BR (BCMV response)", "CBB score", "BCMV presence/absence",
              "CTV presence/absence", "DF (days)", "DM (days)", "EV score",
              "GH type", "HB score", "HI (%)", "LG score", "PH (cm)",
              "RR score", "RU score", "SA score", "SF (days)", "SW (mg)",
              "SY (kg/ha)", "WM score", "ZN score") #


for(i in seq_along(PHE1)){
  mlm_stats <- load_Tassel_MLM_stats(path = tasselpath,
                                     gwasname = gwasname[i])

  panel_a <- tassel_plot_ggman_labeled(Model_pheno = mlm_stats)
  panel_b <- tassel_plot_qq(Model_pheno = mlm_stats)
  panel_c <- tassel_plot_resid_density(tasselplotpheno = METGPHEN[i],
                                       tasselpheno = PHE1[i],
                                       Model_pheno = mlm_stats)
  row1 <- plot_grid(panel_a, panel_b, panel_c, align = "hv", nrow = 1,
                    rel_widths = c(1/2, 1/4, 1/4), labels = "AUTO")
  ggsave(row1, file = paste0(gwasname[i], "_Manhattan_", Sys.Date(),".png"),
         width = 9, height = 2.75*.9, units = "in", dpi = 400)
}


for(i in seq_along(PHE1)){
  mlm_stats <- load_Tassel_MLM_stats(path = tasselpath,
                                     gwasname = gwasname[i])

  panel_a <- tassel_plot_ggman_labeled(Model_pheno = mlm_stats)
  panel_b <- tassel_plot_qq(Model_pheno = mlm_stats)
  panel_c <- tassel_plot_resid_density(tasselplotpheno = METGPHEN[i],
                                       tasselpheno = PHE1[i],
                                       Model_pheno = mlm_stats)
  row1 <- plot_grid(panel_a, panel_b, panel_c, align = "hv", nrow = 1,
                    rel_widths = c(1/2, 1/4, 1/4), labels = "AUTO")
  panels_row2to4 <- tassel_plot_topsnp_effects(
    phenofile = phenofile[i],
    tasselpheno = eval(PHE1[i]),
    tasselplotpheno = METGPHEN[i],
    metpheno = eval(PHE1[i]),
    metplotpheno = METGPHEN[i],
    Model_pheno = mlm_stats)


  row1 <- plot_grid(panel_a, panel_b, panel_c, align = "hv", nrow = 1,
                    rel_widths = c(1/2, 1/4, 1/4), labels = "AUTO")
  ggsave(row1, file = paste0(gwasname[i], "_Manhattan_", Sys.Date(),".png"),
         width = 9, height = 2.75*.9, units = "in", dpi = 400)

  row2to4 <- plot_grid(panels_row2to4[[1]], panels_row2to4[[2]], panels_row2to4[[3]],
                       panels_row2to4[[4]], panels_row2to4[[5]], panels_row2to4[[6]],
                       nrow = 3, axis = "lb", rel_widths = c(2/3, 1/3),
                       labels = c("D", "E", "F", "G", "H", "I"))
  ggsave(row2to4, file = paste0(gwasname, "_TopSNPEffects_", Sys.Date(),".png"),
         width = 9, height = 2.75*3*.9, units = "in", dpi = 400)

  plot_grid(row1, row2to4, axis = "lr", nrow = 2, rel_heights = c(1/4, 3/4))

  #row2 <- plot_grid(panel_d, panel_e, nrow = 1, axis = "lb", rel_widths = c(2/3, 1/3), labels = c("D", "E"))
  #plot_grid(row1, row2, axis = "lr", nrow = 2)
  ggsave(file = paste0(gwasname, "_Manhattan_TopSNPEffects_", Sys.Date(),".png"),
         width = 10*.9, height = 2.75*4*.9, units = "in", dpi = 400)

}
