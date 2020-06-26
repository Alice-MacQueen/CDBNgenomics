#' Identify phenotype names from bigsnpr results in a folder.
#'
#' @description Creates a vector of phenotype names from bigsnpr results.
#'
#' @param path File path to the files from bigsnpr, a character string.
#'     Defaults to the current working directory.
#' @param pattern Pattern within the filename to match. Default is "*.rds".
#'
#' @return A vector of phenotype names.
#'
#' @examples
#' \dontrun{get_results_in_folder(path = system.file("inst/extdata",
#' package = "switchgrassGWAS"))}
#' \dontrun{get_results_in_folder(path = "path/to/gwas/results")}
#'
#' @export
get_results_in_folder <- function(path = ".", pattern = "*.rds"){
  result_files <- list.files(path = path, pattern = pattern)
  return(result_files)
}


#' Step One of bigsnp2mashr
#'
#' @param path Path
#' @param gwas_rds RDS file with gwas results
#' @param phenotype Character vector. Single phenotype name
#' @param numSNPs Integer. Number of top SNPs to choose.
#' @param markers Marker CHR & POS for the GWAS you ran
#'
#' @import bigsnpr
#' @importFrom stats predict
#' @importFrom dplyr left_join mutate top_n select
#' @importFrom tibble as_tibble
#' @importFrom magrittr %>%
get_top_effects_log10p <- function(path, gwas_rds, phenotype, numSNPs,
                                     markers){
  gwas_obj <- readRDS(file = file.path(path, gwas_rds))
  if(attr(gwas_obj, "class")[1] == "mhtest" & is.null(markers)){
    stop(paste0("For bigsnpr output of class mhtest, you must provide a ",
                "data frame called 'markers' containing columns 'CHR' and ",
                "'POS', the chromosome and position for each SNP tested."))
  }
  if(attr(gwas_obj, "class")[1] == "mhtest"){
    message("Finding top SNPs for a mhtest object - bigsnpr output.")
    log10p <- predict(gwas_obj)
    # Use bigsnpr predict function to generate p-values
    pre_mash_1 <- as_tibble(cbind(gwas_obj, markers, log10p))
    # Use SNP object and gwas_obj generated using SNP object to make a full
    # pre_mash object for one phenotype.
    top_set <- pre_mash_1 %>%
      top_n(-numSNPs, .data$log10p) %>%
      mutate(score_log10p = -log10p) %>%
      dplyr::select(.data$CHR, .data$POS, .data$score_log10p)
    names(top_set)[3] <- paste0(phenotype, "_score_log10p")
    # name top_set columns appropriately so that the joined file will have
    # informative column names.
    rm(log10p)

  } else if(attr(gwas_obj, "class")[1] == "tbl_df" &
            "log10p" %in% names(gwas_obj) &
            "CHR" %in% names(gwas_obj) &
            "POS" %in% names(gwas_obj) &
            "bigsnpscore" %in% names(gwas_obj) &
            "estim" %in% names(gwas_obj) &
            "std_err" %in% names(gwas_obj)){

    #message(paste0("Finding top SNPs for bigsnpr output saved by ",
    #        "'cdbn_standard_gwas()'."))
    pre_mash_1 <- tibble(CHR = gwas_obj$CHR, POS = gwas_obj$POS,
                         bigsnpscore = gwas_obj$bigsnpscore,
                         log10p = gwas_obj$log10p)
    if(gwas_obj$log10p[1] > 0){
      top_set <- pre_mash_1 %>%
        top_n(numSNPs, .data$log10p) %>%
        mutate(score_log10p = abs(.data$log10p)) %>%
        dplyr::select(.data$CHR, .data$POS, .data$score_log10p)
      names(top_set)[3] <- paste0(phenotype, "_score_log10p")
      # name top_set columns appropriately so that the joined file will have
      # informative column names.
    } else {
      top_set <- pre_mash_1 %>%
        top_n(numSNPs, .data$log10p) %>%
        mutate(score_log10p = abs(.data$log10p)) %>%
        dplyr::select(.data$CHR, .data$POS, .data$score_log10p)
      names(top_set)[3] <- paste0(phenotype, "_score_log10p")
    }
  } else {
    stop(paste0("RDS objects need to be of class 'mhtest' or of class ",
                "'tbl_df' with columns 'CHR', 'POS', 'log10p', 'estim', ",
                "'std_err', and 'bigsnpscore'."))
    # Don't need estim and std_err until s_hat_bigsnp() step but better to fail
    # loudly early.
  }
  rm(pre_mash_1)
  rm(gwas_obj)
  return(top_set)
}

#' Step Two of bigsnp2mashr
#'
#' @param path Path
#' @param gwas_rds RDS file with gwas results
#' @param phenotype Character vector. Single phenotype name
#' @param top_set Top markers chosen
#' @param random_sample Numeric vector. Random sample of SNPs
#' @param markers Marker names for the GWAS you ran
#' @param markers2 Markers to include if SNPs are LD clumped.
#' @param model One of linear or logistic. Type of GWAS model.
#' @param clump Logical. Clump SNPs?
#'
#' @import bigsnpr
#' @importFrom stats predict
#' @importFrom dplyr left_join mutate
s_hat_bigsnp <- function(path, gwas_rds, phenotype, top_set, random_sample,
                         markers = NULL, markers2 = NULL,
                         model = c("linear", "logistic"), clump = TRUE){

  gwas_obj <- readRDS(file = file.path(path, gwas_rds))

  if(attr(gwas_obj, "class")[1] == "mhtest"){
    log10p <- predict(gwas_obj)
    # Use bigsnpr predict function to generate p-values
    pre_mash_1 <- as_tibble(cbind(gwas_obj, markers, log10p))

  } else if(attr(gwas_obj, "class")[1] == "tbl_df" &
            "log10p" %in% names(gwas_obj) &
            "CHR" %in% names(gwas_obj) &
            "POS" %in% names(gwas_obj) &
            "bigsnpscore" %in% names(gwas_obj) &
            "estim" %in% names(gwas_obj) &
            "std_err" %in% names(gwas_obj)){
    # CHR POS estim std.err tscore log10p stderr_d effect_d
    pre_mash_1 <- gwas_obj %>% dplyr::select(.data$CHR, .data$POS, .data$estim,
                                             .data$std_err, .data$bigsnpscore,
                                             .data$log10p)

  } else {
    stop(paste0("RDS objects need to be of class 'mhtest' or of class ",
                "'tbl_df' with columns 'CHR', 'POS', 'log10p', 'estim', ",
                "'std_err', and 'bigsnpscore'."))
  }
  standardization <- max(abs(pre_mash_1$estim), na.rm = TRUE)

  pre_mash_strong <- top_set %>%
    dplyr::select(.data$CHR, .data$POS, .data$max_score_log10p) %>%
    dplyr::left_join(pre_mash_1, by = c("CHR", "POS")) %>%
    dplyr::mutate(stderr_d = .data$std_err / standardization,
                  effect_d = .data$estim / standardization,
                  stderr_d = case_when(is.na(.data$stderr_d) ~ 10,
                                       is.nan(.data$stderr_d) ~ 10,
                                       TRUE ~ .data$stderr_d),
                  effect_d = case_when(is.na(.data$effect_d) ~ 0,
                                       is.nan(.data$effect_d) ~ 0,
                                       TRUE ~ .data$effect_d))
  if(clump == FALSE){
    pre_mash_random <- pre_mash_1[random_sample,] %>%
      dplyr::mutate(stderr_d = .data$std_err / standardization,
                    effect_d = .data$estim / standardization,
                    stderr_d = case_when(is.na(.data$stderr_d) ~ 10,
                                         is.nan(.data$stderr_d) ~ 10,
                                         TRUE ~ .data$stderr_d),
                    effect_d = case_when(is.na(.data$effect_d) ~ 0,
                                         is.nan(.data$effect_d) ~ 0,
                                         TRUE ~ .data$effect_d))
  } else {
    markers2 <- markers2 %>%
      dplyr::select(.data$CHR, .data$POS)
    pre_mash_random <- markers2 %>%
      left_join(pre_mash_1, by = c("CHR", "POS"))
    pre_mash_random <- pre_mash_random[random_sample,] %>%
      dplyr::mutate(stderr_d = .data$std_err / standardization,
                    effect_d = .data$estim / standardization,
                    stderr_d = case_when(is.na(.data$stderr_d) ~ 10,
                                         is.nan(.data$stderr_d) ~ 10,
                                         TRUE ~ .data$stderr_d),
                    effect_d = case_when(is.na(.data$effect_d) ~ 0,
                                         is.nan(.data$effect_d) ~ 0,
                                         TRUE ~ .data$effect_d))
  }

  if((model == "linear" & attr(gwas_obj, "class")[1] == "mhtest") |
     attr(gwas_obj, "class")[1] == "tbl_df"){
    # CHR POS max_score_log10p CHRN estim std.err score log10p stderr_d effect_d
    # CHR POS estim std.err tscore log10p stderr_d effect_d

    #names(pre_mash_strong)[3] <- paste0(phenotype, "_maxscore")
    names(pre_mash_strong)[4] <- paste0("Bhat_", phenotype)
    names(pre_mash_strong)[5] <- paste0("Shat_",phenotype)
    names(pre_mash_strong)[6] <- paste0(phenotype, "_tscore")
    names(pre_mash_strong)[7] <- paste0(phenotype, "_log10p")
    names(pre_mash_strong)[8] <- paste0("Stand_Shat", phenotype)
    names(pre_mash_strong)[9] <- paste0("Stand_Bhat", phenotype)

    names(pre_mash_random)[3] <- paste0("Bhat_", phenotype)
    names(pre_mash_random)[4] <- paste0("Shat_",phenotype)
    names(pre_mash_random)[5] <- paste0(phenotype, "_tscore")
    names(pre_mash_random)[6] <- paste0(phenotype, "_log10p")
    names(pre_mash_random)[7] <- paste0("Stand_Shat", phenotype)
    names(pre_mash_random)[8] <- paste0("Stand_Bhat", phenotype)
  }
  if(model == "logistic" & attr(gwas_obj, "class")[1] == "mhtest"){
    # CHR POS max_score_log10p CHRN estim std.err score log10p stderr_d effect_d
    # CHR POS estim std.err niter zscore log10p stderr_d effect_d
    #names(pre_mash_strong)[3] <- paste0(phenotype, "_maxscore")
    names(pre_mash_strong)[4] <- paste0("Bhat_", phenotype)
    names(pre_mash_strong)[5] <- paste0("Shat_",phenotype)
    names(pre_mash_strong)[6] <- paste0(phenotype, "_niter")
    names(pre_mash_strong)[7] <- paste0(phenotype, "_zscore")
    names(pre_mash_strong)[8] <- paste0(phenotype, "_log10p")
    names(pre_mash_strong)[9] <- paste0("Stand_Shat_",phenotype)
    names(pre_mash_strong)[10] <- paste0("Stand_Bhat_",phenotype)

    names(pre_mash_random)[3] <- paste0("Bhat_", phenotype)
    names(pre_mash_random)[4] <- paste0("Shat_",phenotype)
    names(pre_mash_random)[5] <- paste0(phenotype, "_niter")
    names(pre_mash_random)[6] <- paste0(phenotype, "_zscore")
    names(pre_mash_random)[7] <- paste0(phenotype, "_log10p")
    names(pre_mash_random)[8] <- paste0("Stand_Shat_",phenotype)
    names(pre_mash_random)[9] <- paste0("Stand_Bhat_",phenotype)
  }
  # name top_set columns appropriately so that the joined file will have
  # informative column names.

  return(list(strong_df_1 = pre_mash_strong, random_df_1 = pre_mash_random))
}


#' Convert bigsnpr output to mashr input dataframes.
#'
#' @description This function converts bigsnpr output, saved as rds files to
#'    the specified path, to four dataframes used in the R package mashr. It
#'    can clump SNPs based on LD and the maximum -log10(p-value) across all
#'    included GWAS. It can also set the random effect data frames to come from
#'    a subsample of SNPs clumped by MAF and LD.
#'
#' @param path File path to the rds files saved from bigsnpr, a character
#'     string. Defaults to the working directory.
#' @param gwas_rds A character vector of saved gwas rds objects from bigsnpr. If NA, all *.rds files in the path will be used.
#' @param phenotypes A character vector of phenotype names for the GWAS RDS
#'    objects. Must be the same length as gwas_rds, or NA. If NA, these will be
#'    the rds file names.
#' @param clump Logical. Should SNPs be clumped by LD & p-value to
#'    standardize signal strength across different LD blocks? Default is TRUE.
#' @param scaled Logical. Should marker effects in each condition be
#'    scaled to fall between -1 and 1? Default is TRUE.
#' @param snp The "bigSNP" object used to run the gwas; needed if clump is
#'    TRUE. Load with bigsnpr::snp_attach().
#' @param numSNPs The number of most significant SNPs selected from each GWAS.
#'     Ideally this will give 1 million or fewer total cells in the resultant
#'     mash dataframes. Defaults to 1000.
#' @param model Regression used in bigstatsr. One of "logistic" or "linear".
#'     Default is "linear".
#' @param saveoutput Logical. Should the function's output also be saved to RDS
#' files? Default is FALSE.
#'
#' @return A list containing five data frames: the SNPs selected, the B_hat
#'    and S_hat matrices for the strong SNP set and for a random SNP set that
#'    is twice the size.
#'
#' @note To create a vector of phenotype names, use the
#'     \code{\link{get_results_in_folder}} function.
#'
#' @examples
#' \dontrun{cdbn_bigsnp2mashr(path = system.file("inst/extdata"), numSNPs = 20,
#'     model = "linear")}
#' \dontrun{cdbn_bigsnp2mashr(numSNPs = 10000, model = "logistic")}
#' \dontrun{cdbn_bisgnp2mashr(numSNPs = 20000, model = "linear", saveoutput = TRUE)}
#' \dontrun{phenotype_vector <- get_results_in_folder(path = system.file(
#'     "inst/extdata"))
#'     numSNPs <- 1000000 / length(phenotype_vector)^2
#'     cdbn_bigsnp2mashr(phenotypes = phenotype_vector, numSNPs = numSNPs,
#' model = "linear", saveoutput = TRUE)}
#'
#' @import bigsnpr
#' @importFrom dplyr full_join left_join case_when arrange mutate select slice filter group_by n add_tally
#' @importFrom magrittr %>%
#' @importFrom stringr str_sub
#' @importFrom tidyr gather unite
#' @importFrom tidyselect starts_with
#'
#' @export
cdbn_bigsnp2mashr <- function(path = ".", snp = NULL, gwas_rds = NA,
                              phenotypes = NA, clump = TRUE, scaled = TRUE,
                              numSNPs = 1000, model = c("linear", "logistic"),
                              saveoutput = FALSE){
  match.arg(model, c("linear", "logistic"))
  if(is.null(snp)){
    stop(paste0("The bigSNP object that was used to run the GWAS should be ",
                "provided - load this in to R using 'bigsnpr::snp_attach()', ",
                "then specify this value as 'snp'."))
  }
  if(attr(snp, "class") != "bigSNP"){
    stop("snp needs to be a bigSNP object, produced by the bigsnpr package.")
  }
  if(is.na(gwas_rds)[1]){
    phe_col <- get_results_in_folder(path = path, pattern = "GWAS_datatable")
  } else {
    phe_col <- gwas_rds
  }
  if(is.null(phe_col) | is.na(phe_col[1])){
    stop("Can't find any rds files that match 'GWAS_datatable' in this path.")
  }
  if(is.na(phenotypes)[1]){
    message("Phenotypes weren't provided, so using rds file name.")
    phenotypes <- str_sub(phe_col, start = 16, end = -5)
  }
  numSNPs <- as.integer(numSNPs)
  G <- snp$genotypes
  markers <- tibble(CHR = snp$map$chromosome, POS = snp$map$physical.pos) %>%
    mutate(marker = paste(.data$CHR, .data$POS, sep = "_"))

  message(paste0("Starting part one: Making a data frame of all SNPs ",
                 "in the top ", numSNPs, " SNPs
                 by maximum -log10(p-values) for at least one phenotype."))

  top_set <- get_top_effects_log10p(path = path, gwas_rds = phe_col[1],
                                      phenotype = phenotypes[1],
                                      numSNPs = numSNPs, markers = markers)

  for(i in seq_along(phe_col)[-1]){
    message(paste0("Picking top SNPs for GWAS for ", phenotypes[i]," (", i,
                   " of ", length(phe_col), ")."))
    top_set_new <- get_top_effects_log10p(path = path, gwas_rds = phe_col[i],
                                            phenotype = phenotypes[i],
                                            numSNPs = numSNPs,
                                            markers = markers)
    top_set <- top_set %>%
      full_join(top_set_new, by = c("CHR", "POS"))
  }
  top_set <- top_set %>%
    dplyr::arrange(.data$CHR, .data$POS) %>%
    gather(key = "Condition", value = "score_log10p", -(1:2)) %>%
    group_by(.data$CHR, .data$POS) %>%
    filter(!is.na(.data$score_log10p)) %>%
    mutate(max_score_log10p = max(.data$score_log10p)) %>%
    add_tally(name = "n_cond_in_top_set") %>%
    slice(which.max(.data$max_score_log10p)) %>%
    dplyr::select(-.data$score_log10p) %>%
    arrange(.data$CHR, .data$POS) %>%
    mutate(marker = paste(.data$CHR, .data$POS, sep = "_"))
  if(clump == TRUE){
    all_clumps <- snp_clumping(G, infos.chr = markers$CHR, ncores = nb_cores(),
                               infos.pos = markers$POS)
    markers2 <- markers[all_clumps,] # random marker data frame
    top_subset <- which(markers$marker %in% top_set$marker) # subset for ind.col
    rdsfile <- snp_subset(snp, ind.col = top_subset) # make bigSNP subset
    topsnp <- snp_attach(rdsfile) # load bigSNP for top set
    top_clumps <- snp_clumping(topsnp$genotypes, infos.chr = top_set$CHR,
                               S = top_set$max_score_log10p,
                               ncores = nb_cores(), infos.pos = top_set$POS)
    # Outputs a numeric vector corresponding to clumped rows w/ the best p-value
    top_set <- top_set[top_clumps,]
  } else {
    top_set <- top_set# %>%
      #dplyr::select(.data$CHR, .data$POS, .data$max_score_log10p)
  }
  if(saveoutput == TRUE){
    saveRDS(top_set, file = file.path(path, paste0("Part-One-Output_",
                                                   "Top-Effects-", numSNPs,
                                                   "-SNPs.rds")))
  }

  message(paste0("Part One: data frame of SNPs to keep complete.
                 ----------------"))
  message(paste0("Starting Part Two: Creating strong and random dataframes of ",
                 "B_hat and S_hat
                 values for use in mashr."))

  set.seed(1234) # Makes the random data frames reproducible.
  if(clump == TRUE){
    random_sample <- sample(1:nrow(markers2), length(top_clumps)*2) %>%
      sort()
  } else {
    random_sample <- sample(1:nrow(markers), nrow(top_set)*2) %>%
      sort()
  }
  mash_list_1 <- s_hat_bigsnp(path = path, gwas_rds = phe_col[1],
                              phenotype = phenotypes[1],
                              top_set = top_set,
                              random_sample = random_sample, model = model,
                              markers = markers, markers2 = markers2,
                              clump = clump)

  mash_df_strong <- mash_list_1$strong_df_1
  mash_df_random <- mash_list_1$random_df_1

  for(i in seq_along(phe_col)[-1]){
    if(scaled == TRUE){
      message(paste0("Determining standardized B_hat and S_hat values for ",
                     nrow(top_set), " SNPs from ", phenotypes[i],
                     " GWAS (phenotype ", i, " of ", length(phe_col), ")."))
    } else {
      message(paste0("Determining unscaled B_hat and S_hat values for ",
                     nrow(top_set), " SNPs from ", phenotypes[i],
                     " GWAS
                     (phenotype ", i, " of ", length(phe_col), ")."))
    }
    mash_list_1 <- s_hat_bigsnp(path = path, gwas_rds = phe_col[i],
                                phenotype = phenotypes[i],
                                top_set = top_set,
                                random_sample = random_sample, model = model,
                                markers = markers, markers2 = markers2,
                                clump = clump)

    mash_df_strong <- mash_df_strong %>%
      dplyr::left_join(mash_list_1$strong_df_1, by = c("CHR", "POS",
                                                       "max_score_log10p"))
    mash_df_random <- mash_df_random %>%
      dplyr::left_join(mash_list_1$random_df_1, by = c("CHR", "POS"))
  }
  if(scaled == TRUE){
    bhat_strong <- mash_df_strong %>%
      unite(col = "SNP", .data$CHR, .data$POS) %>%
      dplyr::select(.data$SNP, tidyselect::starts_with("Stand_Bhat"))
    shat_strong <- mash_df_strong %>%
      unite(col = "SNP", .data$CHR, .data$POS) %>%
      dplyr::select(.data$SNP, tidyselect::starts_with("Stand_Shat"))
    bhat_random <- mash_df_random %>%
      unite(col = "SNP", .data$CHR, .data$POS) %>%
      dplyr::select(.data$SNP, tidyselect::starts_with("Stand_Bhat"))
    shat_random <- mash_df_random %>%
      unite(col = "SNP", .data$CHR, .data$POS) %>%
      dplyr::select(.data$SNP, tidyselect::starts_with("Stand_Shat"))
  } else {
    bhat_strong <- mash_df_strong %>%
      unite(col = "SNP", .data$CHR, .data$POS) %>%
      dplyr::select(.data$SNP, tidyselect::starts_with("Bhat"))
    shat_strong <- mash_df_strong %>%
      unite(col = "SNP", .data$CHR, .data$POS) %>%
      dplyr::select(.data$SNP, tidyselect::starts_with("Shat"))
    bhat_random <- mash_df_random %>%
      unite(col = "SNP", .data$CHR, .data$POS) %>%
      dplyr::select(.data$SNP, tidyselect::starts_with("Bhat"))
    shat_random <- mash_df_random %>%
      unite(col = "SNP", .data$CHR, .data$POS) %>%
      dplyr::select(.data$SNP, tidyselect::starts_with("Shat"))
    # make a SNP column by uniting CHR and POS
    # make Bhat only and Shat only tables for both strong and random sets.
  }
  log10p_strong <- mash_df_strong %>%
    unite(col = "SNP", .data$CHR, .data$POS, remove = FALSE) %>%
    dplyr::select(.data$SNP, .data$CHR, .data$POS, tidyselect::ends_with("_log10p"))
  log10p_random <- mash_df_random %>%
    unite(col = "SNP", .data$CHR, .data$POS, remove = FALSE) %>%
    dplyr::select(.data$SNP, .data$CHR, .data$POS, tidyselect::ends_with("_log10p"))
  B_hat_random <- data.frame(bhat_random, row.names = "SNP")
  S_hat_random <- data.frame(shat_random, row.names = "SNP")
  B_hat_strong <- data.frame(bhat_strong, row.names = "SNP")
  S_hat_strong <- data.frame(shat_strong, row.names = "SNP")

  if(saveoutput == TRUE){
    saveRDS(B_hat_strong, file = file.path(path, paste0("B_hat_strong_df_",
                                                        numSNPs, "topSNPs.rds")))
    saveRDS(S_hat_strong, file = file.path(path, paste0("S_hat_strong_df_",
                                                        numSNPs, "topSNPs.rds")))
    saveRDS(B_hat_random, file = file.path(path, paste0("B_hat_random_df_",
                                                        numSNPs, "topSNPs.rds")))
    saveRDS(S_hat_random, file = file.path(path, paste0("S_hat_random_df_",
                                                        numSNPs, "topSNPs.rds")))
    saveRDS(log10p_strong, file = file.path(path, paste0("log10p_strong_df_",
                                                         numSNPs, "topSNPs.rds")))
    saveRDS(log10p_random, file = file.path(path, paste0("log10p_random_df_",
                                                         numSNPs, "topSNPs.rds")))
  }
  return(list(top_set = top_set,
              random_sample = random_sample,
              B_hat_strong = B_hat_strong,
              S_hat_strong = S_hat_strong,
              B_hat_random = B_hat_random,
              S_hat_random = S_hat_random))
}
