#! /usr/bin/Rscript
# mashr JLS run on 22 METG BLUP phenotypes with best number of PCs for kinship
# make sure to export your LD_LIBRARY_PATH before running this script:
# LD_LIBRARY_PATH=/home/pubjuenger/Alice/bin/gsl/lib:/home/pubjuenger/Alice/bin/ed/
# export LD_LIBRARY_PATH
# Make sure you run this script using R-3.4.2, which has mashr installed, by:
# nohup ~/Alice/R_mash/bin/Rscript METG_GWAS_mash_2019-03-21.R &


# Read in the random and the strong datasets - random dataset is ~3-4x as big
# Make sure you've removed doubles and NA's from your dataset - or set Shat NA's
# to a very large  number so that that condition will be discounted by mashr.
load_mash_df <- function(path, numSNPs){
  snp_df <- readRDS(file.path(path, paste0("effects_",
                                           numSNPs, "SNPs_PartOneOutput.rds")))
  Bhat_random <- readRDS(file.path(path, paste0("B_hat_random_df_",
                                                numSNPs, "topSNPs.rds")))
  Shat_random <- readRDS(file.path(path, paste0("S_hat_random_df_",
                                                numSNPs, "topSNPs.rds")))
  Bhat_strong <- readRDS(file.path(path, paste0("B_hat_strong_df_",
                                                numSNPs, "topSNPs.rds")))
  Shat_strong <- readRDS(file.path(path, paste0("S_hat_strong_df_",
                                                numSNPs, "topSNPs.rds")))
  return(list(snp_df = snp_df, B_hat_strong = B_hat_strong,
              S_hat_strong = S_hat_strong, B_hat_random = B_hat_random,
              S_hat_random = S_hat_random))
}

mash_standard_run <- function(path, gapit2mash_obj = NA, numSNPs = NA,
                              saveoutput = FALSE, U_ed = NA){
  if(!is.na(gapit2mash_obj)){
    g2mobj <- g2mobj
  } else if (is.na(gapit2mash_obj) & !is.na(numSNPs)){
    g2mobj <- load_mash_df(path = path, numSNPs = numSNPs)
  } else stop("One of 'gapit2mash_obj' or 'numSNPs' needs to be specified to load the data frames needed for mash.")

  Vhat <- make_Vhat(g2mobj = g2mobj)

  data_strong = mash_set_data(g2mobj$Bhat_strong, g2mobj$Shat_strong, V=Vhat)
  data_random = mash_set_data(g2mobj$Bhat_random, g2mobj$Shat_random, V=Vhat)
  U_c = cov_canonical(g2mobj$data_random)
  if(!is.na(U_ed)){
  # estimate data-driven covariance using strong dataset
    U_ed <- make_U_ed(data_strong = data_strong, saveoutput = saveoutput)
  } else if(is.character(U_ed)){
    U_ed <- readRDS(U_ed)
  }
  # Run mash on the random dataset using the random data w/ correlation structure
  m = mash(data_random, Ulist = c(U_ed, U_c), outputlevel = 1)
  if(saveoutput == TRUE){
  saveRDS(m, file.path(path, "Model_of_random_tests_", numSNPs, "SNPs.rds"))
  }
  # Run mash on the strong dataset (or all data) using
  # the previous results from the random data
  m2 = mash(data_strong, g = get_fitted_g(m), fixg = TRUE)
  if(saveoutput == TRUE){
  saveRDS(m2, file.path(path, "Strong_Effects", numSNPs, "SNPs.rds"))
  }
  print("How many significant markers?")
  print(length(get_significant_results(m2)))
  if(saveoutput == TRUE){
  saveRDS(get_pairwise_sharing(m2),
          file = file.path(path, "Pairwise_sharing_Strong_Effects_", numSNPs,
                           "SNPs.rds"))
  saveRDS(get_pairwise_sharing(m2, factor = 0),
          file = file.path(path,
                           "Pairwise_sharing_sign_Strong_Effects_", numSNPs,
                           "SNPs.rds"))
  }
  return(m2)
}




make_Vhat <- function(g2mobj){
  # estimate data correlation structure using random dataset
  data_r <- mash_set_data(g2mobj$Bhat_random, g2mobj$Shat_random)
  z <- g2mobj$Bhat_random/g2mobj$Shat_random
  max_absz <- apply(abs(z), 1, max)
  nullish <-  which(max_absz < 1.5)
  nullish_z <-  z[nullish,]
  Vhat <- cor(nullish_z)
  return(Vhat)
}

make_U_ed <- function(data_strong, saveoutput = FALSE){
  U_pca = cov_pca(data_strong, 5)
  U_ed = cov_ed(data_strong, U_pca)
  if(saveoutput = TRUE){
  saveRDS(U_ed, file.path(path, paste0(
                          "Mash_data_driven_covariances_", numSNPs, ".rds")))
  # extreme decomposition takes a long time to run, so save the result.
  }
  return(U_ed)
}









