#' Read in the random and the strong datasets
#'
#' @param path File path to the rds files saved from bigsnpr, a character
#'     string. Defaults to the working directory.
#' @param numSNPs The number of most significant SNPs selected from each GWAS.
#'     Ideally this will give 1 million or fewer total cells in the resultant
#'     mash dataframes. Defaults to 1000.
#' @param suffix Character. Optional. Should be the unique suffix used to
#'     save \code{cdbn_bigsnp2mashr} output as RDS files, if it was used.
#'
#' @return A list containing five data frames: the SNPs selected, the B_hat
#'    and S_hat matrices for the strong SNP set and for a random SNP set that
#'    is twice the size.
#'
#' @note Make sure you've removed doubles and NA's from your dataset - or set
#'     Shat NA's to a very large  number so that that condition will be
#'     discounted by mashr.
load_g2m_df <- function(path, numSNPs, suffix){
  if(!(str_sub(suffix, end = 1) %in% c("", "_")))
  { suffix <- paste0("_", suffix) }
  top_set <- readRDS(file.path(path, paste0("Part-One-Output_",
                                            "Top-Effects-", numSNPs,
                                            "-SNPs", suffix, ".rds")))
  B_hat_random <- readRDS(file.path(path, paste0("B_hat_random_df_",
                                                 numSNPs, "topSNPs",
                                                 suffix, ".rds")))
  S_hat_random <- readRDS(file.path(path, paste0("S_hat_random_df_",
                                                 numSNPs, "topSNPs",
                                                 suffix, ".rds")))
  B_hat_strong <- readRDS(file.path(path, paste0("B_hat_strong_df_",
                                                 numSNPs, "topSNPs",
                                                 suffix, ".rds")))
  S_hat_strong <- readRDS(file.path(path, paste0("S_hat_strong_df_",
                                                 numSNPs, "topSNPs",
                                                 suffix, ".rds")))
  return(list(top_set = top_set, B_hat_strong = B_hat_strong,
              S_hat_strong = S_hat_strong, B_hat_random = B_hat_random,
              S_hat_random = S_hat_random))
}

make_U_ed <- function(path, data_strong, numSNPs, saveoutput = FALSE, suffix){
  requireNamespace("mashr")
  if(ncol(data_strong$Bhat < 6)){
    U_pca = mashr::cov_pca(data_strong, (ncol(data_strong$Bhat)-1))
  } else {
    U_pca = mashr::cov_pca(data_strong, 5)
  }
  U_ed = mashr::cov_ed(data_strong, U_pca)
  if(saveoutput == TRUE){
    saveRDS(U_ed, file.path(path, paste0(
      "Mash_data_driven_covariances_", numSNPs, suffix, ".rds")))
    # extreme decomposition takes a long time to run, so save the result.
  }
  return(U_ed)
}

get_loglik=function(a){a$loglik}

#' A standard run of mashr
#'
#' @description If you have prepared mash input using the bigsnp2mashr function
#'     or the gapit2mashr R package, use this function on the output to run mash
#'     as recommended by vignettes in the mashr package.
#'
#' @param path File path to the rds files saved from gapit2mashr, or
#'     bigsnpr2mashr, as a character string. Defaults to the working directory.
#' @param list_input A list containing five data frames: the SNPs selected, the
#'     B_hat and S_hat matrices for the strong SNP set and for a random SNP set
#'     chosen in gapit2mashr.
#' @param numSNPs The number of most significant SNPs selected from each GWAS.
#'     Ideally this will give 1 million or fewer total cells in the resultant
#'     mash dataframes. Defaults to 1000.
#' @param suffix Character. Optional. Should be the unique suffix used to
#'     save \code{cdbn_bigsnp2mashr} output as RDS files, if it was used.
#' @param saveoutput Logical. Should the function's output also be saved to RDS
#' files? Default is FALSE.
#' @param U_ed An optional list of data-driven covariance matrices, or a
#'     character vector containing the complete path to a .rds containing these
#'     matrices.
#' @param U_hyp An optional list of covariance matrices specified by the user.
#'
#' @return A mash result, manipulable using functions in mashr and by mash_plot
#'     functions in the CDBNgenomics package.
#'
#' @note This is a convenience function for users who have prepared their data
#'     using gapit2mashr or bigsnp2mashr. If you have not used these functions
#'     to make your mash input, you should not use this function - instead,
#'     follow the recommendations of the vignettes in the mashr package itself.
#'
#' @importFrom ashr get_fitted_g
#'
#' @export
mash_standard_run <- function(path = ".", list_input = NA, numSNPs = NA,
                              suffix = "", saveoutput = FALSE, U_ed = NA,
                              U_hyp = NA){
  requireNamespace("mashr")
  if(!(str_sub(suffix, end = 1) %in% c("", "_")))
  { suffix <- paste0("_", suffix) }
  if(!(typeof(list_input) %in% c("list", "logical"))){
    stop("list_input needs to be a list object with at least four matrices:
    B_hat_strong, S_hat_strong, B_hat_random, S_hat_random. Make this using gapit2mashr or bigsnpr2mashr.")
    } else if (!is.na(list_input[1])){
    list_input <- list_input
    numSNPs <- nrow(list_input$B_hat_strong)
    } else if (is.na(list_input[1]) & !is.na(numSNPs)){
    list_input <- load_g2m_df(path = path, numSNPs = numSNPs, suffix = suffix)
  } else stop("A bigsnpr2mashr//gapit2mashr output as 'list_input', or 'numSNPs' and 'suffix', the number of SNPs used in gapit2mashr//bigsnp2mashr and an optional suffix used in bigsnp2mashr, need to be specified to load the data frames needed for mash.")

  Bhat_strong <- as.matrix(list_input$B_hat_strong)
  Shat_strong <- as.matrix(list_input$S_hat_strong)
  Bhat_random <- as.matrix(list_input$B_hat_random)
  Shat_random <- as.matrix(list_input$S_hat_random)

  if(typeof(U_hyp) == "list"){
    # check that it contains matrices that are the same ncol & nrow as B_hat_strong so that mash will run.
  }

  data_r <- mashr::mash_set_data(Bhat_random, Shat_random)
  message(paste0("Estimating the correlation structure in the null tests from ",
                 "the random data.
                 (not the strong data because it will not necessarily contain
                 any null tests)."))
  Vhat <- mashr::estimate_null_correlation_simple(data = data_r)

  message(paste0("Setting up the main data objects with this correlation ",
                 "structure in place."))
  data_strong <- mashr::mash_set_data(Bhat_strong, Shat_strong, V=Vhat)
  data_random <- mashr::mash_set_data(Bhat_random, Shat_random, V=Vhat)
  U_c <- mashr::cov_canonical(data_random)
  if(typeof(U_ed) == "logical" & is.na(U_ed[1])){
    # estimate data-driven covariance using strong dataset
    message(paste0("Now estimating data-driven covariances using the strong",
                   " tests.
                   NB: This step may take some time to complete."))
    U_ed <- make_U_ed(path = path, data_strong = data_strong,
                      numSNPs = numSNPs, saveoutput = saveoutput,
                      suffix = suffix)
  } else if(is.character(U_ed)){
    U_ed <- readRDS(U_ed)
  } else if(typeof(U_ed) == "list"){
    U_ed <- U_ed
  }
  # Run mash on the random dataset using the random data w/ correlation structure
  message(paste0("Fit mash to the random tests using both data-driven and ",
                 "canonical covariances."))
  if(typeof(U_hyp) == "list"){
    m = mashr::mash(data_random, Ulist = c(U_ed, U_c, U_hyp), outputlevel = 1)
  } else {
    m = mashr::mash(data_random, Ulist = c(U_ed, U_c), outputlevel = 1)
  }
  if(saveoutput == TRUE){
    saveRDS(m, file.path(path, paste0("Model_of_random_tests_", numSNPs,
                                      "SNPs", suffix, ".rds")))
  }
  # Run mash on the strong dataset (or all data) using
  # the previous results from the random data
  message(paste0("Compute posterior matrices for the strong effects",
                 " using the mash fit from the
                 random tests."))
  m2 = mashr::mash(data_strong, g = get_fitted_g(m), fixg = TRUE)
  if(saveoutput == TRUE){
    saveRDS(m2, file.path(path, paste0("Strong_Effects", numSNPs, "SNPs",
                                       suffix, ".rds")))
  }
  print("Log likelihood with specified covariance matrices: ")
  print(get_loglik(m2), digits = 10)
  print("How many significant markers?")
  print(length(get_significant_results(m2)))
  if(saveoutput == TRUE){
    saveRDS(mashr::get_pairwise_sharing(m2),
            file = file.path(path, paste0("Pairwise_sharing_Strong_Effects_",
                                          numSNPs, "SNPs", suffix, ".rds")))
    saveRDS(mashr::get_pairwise_sharing(m2, factor = 0),
            file = file.path(path,
                             paste0("Pairwise_sharing_sign_Strong_Effects_",
                                    numSNPs, "SNPs", suffix, ".rds")))
  }
  return(m2)
}


