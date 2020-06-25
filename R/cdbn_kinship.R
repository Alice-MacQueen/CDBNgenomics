#' @title Get current date-time in a filename-appropriate format.
#'
#' @description Converts the current \code{Sys.time()} system time to a format
#'     that is acceptable to include in a filename. Changes punctuation that
#'     won't work in a filename.
#'
#' @return A string containing the current date-time with spaces and colons
#'     replaced with underscores and periods, respectively.
#'
#' @importFrom stringr str_replace_all
get_date_filename <- function(){
  str_replace_all(str_replace_all(Sys.time(), ":", "."), " ", "_")
}


#' @title Find a kinship matrix using the van Raden method.
#'
#' @description Calculate the kinship matrix using code from GAPIT and bigsnpr
#'     and the methods of VanRaden (2009, J. Dairy Sci. 91:4414???C4423). This
#'     code is a modified version of the same code in GAPIT. Note
#'     that this matrix cannot currently be used with the GWAS methods in
#'     bigsnpr; however, this matrix could be used for other analyses, such as
#'     variance components analyses.
#'
#' @param snp A "bigSNP" object; load with bigsnpr::snp_attach().
#' @param ind.row An integer vector of the rows (individuals) to find a
#'     kinship matrix for. Defaults to all rows.
#' @param hasInbred Logical. Does the SNP file contain inbred individuals or
#'     closely related individuals, like siblings? Default is TRUE.
#' @param saveoutput Logical. Should the output be saved to the working
#'     directory?
#'
#' @return A kinship matrix with labeled rows and columns.
#'
#' @import bigsnpr
#' @import bigstatsr
#'
#' @examples
#' \dontrun{
#' K <- get_kinship(snp = snp, saveoutput = TRUE)
#' }
#'
#' @export
get_kinship <- function(snp, ind.row = NA, hasInbred = TRUE,
                          saveoutput = FALSE){
  if(attr(snp, "class") != "bigSNP"){
    stop("snp needs to be a bigSNP object, produced by the bigsnpr package.")
  }
  if(saveoutput == FALSE){
    message(paste0("'saveoutput' is FALSE, so the kinship matrix will not ",
                   "be saved to the working directory."))
  }

  G <- snp$genotypes
  if(!is.na(ind.row[1])){
    nInd <- length(ind.row)
  } else {
    nInd <- snp$genotypes$nrow
    ind.row <- 1:nInd
  }
  # Centered (mean of rows subtracted) transposed crossproduct of snp file.
  K <- big_tcrossprodSelf(G, ind.row = ind.row,
                          fun.scaling = big_scale(center = TRUE,
                                                  scale = FALSE))

  #Extract diagonals
  i = 1:nInd
  j = (i - 1)*nInd
  index = i + j
  d = K[index]
  DL = min(d)
  DU = max(d)
  floor = min(K[i, i])

  K = (K[i, i] - floor)/(DL - floor)
  MD = (DU - floor)/(DL - floor)

  if(!is.na(ind.row[1])){
    rownames(K) <- snp$fam$sample.ID[ind.row]
    colnames(K) <- snp$fam$sample.ID[ind.row]
  } else {
    rownames(K) <- snp$fam$sample.ID
    colnames(K) <- snp$fam$sample.ID
  }

  if(MD > 2){
    K[index] <- K[index]/(MD-1)+1
  }
  #Handler of inbred
  if(MD < 2 & hasInbred){
    K = 2*K / ((DU - floor) / (DL - floor))
  }

  if(saveoutput){
    saveRDS(K, paste0("Kinship_van_Raden_method_", nInd, "_individuals_",
                      ".rds"))
  }
  return(K)
}

#' @title Wrapper for the snp_autoSVD function for the CDBN.
#'
#' @description This is a wrapper to determine population structure for GWAS
#'     for a bigSNP object. Arguments that are recognized by
#'     bigsnpr::snp_autoSVD can also be specified in this function.
#'
#' @param snp A "bigSNP" object; load with bigsnpr::snp_attach().
#' @param k Integer. The number of principal components to find. Default is 10.
#' @param ncores Integer. Number of cores to use. Default is one.
#' @param saveoutput Logical. Should the output be saved to the working
#'     directory?
#' @param ... Other arguments to \code{\link{snp_autoSVD}}.
#'
#' @return A big_SVD object.
#'
#' @import bigsnpr
#' @import bigstatsr
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate case_when
#' @importFrom tibble enframe
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' svd20 <- get_SVD(snp = snp, k = 20, saveoutput = TRUE)
#' }
#'
#' @export
get_SVD <- function(snp, k = 10, ncores = 1, saveoutput = FALSE, ...){
  requireNamespace("dots") # devtools::install_github("lcolladotor/dots")
  fun.scaling <- dots::dots(name = 'fun.scaling', value = snp_scaleBinom(),
                            ...)
  thr.r2 <- dots::dots(name = 'thr.r2', value = 0.2, ...)
  size <- dots::dots(name = 'size', value = 100/thr.r2, ...)
  roll.size <- dots::dots(name = 'roll.size', value = 50, ...)
  int.min.size <- dots::dots(name = 'int.min.size', value = 20, ...)
  #alpha.tukey <- dots::dots(name = 'alpha.tukey', value = 0.05, ...)
  #min.mac <- dots::dots(name = 'min.mac', value = 10, ...)
  #max.iter <- dots::dots(name = 'max.iter', value = 5, ...)
  verbose <- dots::dots(name = 'verbose', value = FALSE, ...)

  # Argument error checks; Set up numeric chromosome data frame.
  if(attr(snp, "class") != "bigSNP"){
    stop("snp needs to be a bigSNP object, produced by the bigsnpr package.")
  }
  G <- snp$genotypes
  CHR <- snp$map$chromosome
  POS <- snp$map$physical.pos
  plants <- snp$fam$sample.ID

  if(saveoutput == FALSE){
    message(paste0("'saveoutput' is FALSE, so the svd will not be saved to ",
                   "the working directory."))
  }

  # Determine population structure
  svd <- snp_autoSVD(G = G, infos.chr = CHR, infos.pos = POS,
                     ncores = ncores, k = k, fun.scaling = fun.scaling,
                     thr.r2 = thr.r2, size = size, roll.size = roll.size,
                     int.min.size = int.min.size, #alpha.tukey = alpha.tukey,
                     #min.mac = min.mac, max.iter = max.iter,
                     verbose = verbose)
  if(saveoutput){
    saveRDS(svd, file = paste0("SVD_", length(plants), "g_", k, "PCs.rds"))
  }
  return(svd)
}


#' @title Wrapper for my standard GWAS functions for the CDBN.
#'
#' @description This is a wrapper to make my standard GWAS simpler to execute.
#'
#' @param snp A "bigSNP" object; load with bigsnpr::snp_attach(). Here, genomic
#'     information for Phaseolus vulgaris from the Common Dry Bean Nursery.
#' @param df Dataframe of phenotypes where the first column is Taxa.
#' @param type Character string. Type of univarate regression to run for GWAS.
#'     Options are "linear" or "logistic".
#' @param covar Optional covariance matrix to include in the regression. You
#'     can generate these using \code{get_SVD()}.
#' @param ncores Number of cores to use. Default detects this with `nb_cores()`.
#' @param lambdagc Default is TRUE - should lambda_GC be used to find the best
#'     population structure correction? Alternatively, you can provide a data
#'     frame containing "NumPCs" and the phenotype names containing lambda_GC
#'     values. This is saved to the output directory by cdbn_standard_gwas and
#'     otherwise found from a GWAS result  using `bigsnpr:::getLambdaGC()`.
#' @param outputdir String or file.path() to the output directory. Default is
#'     the working directory.
#' @param savegwas Logical. Should the gwas output be saved as a rds to the
#'     working directory? These files are typically quite large. Default is
#'     FALSE.
#' @param saveplots Logical. Should Manhattan and QQ-plots be generated and
#'     saved to the working directory? Default is TRUE.
#' @param saveannos Logical. Should annotation tables for top SNPs be generated
#'     and saved to the working directory? Default is FALSE. Can take
#'     additional arguments; requires a txdb.sqlite object used in
#'     AnnotationDbi.
#' @param txdb A txdb object such as 'Pvulgaris_442_v2.1.gene.sqlite'.
#'     Load this into your environment with AnnotationDbi::loadDb.
#' @param minphe Integer. What's the minimum number of phenotyped individuals
#'     to conduct a GWAS on? Default is 200. Use lower values with caution.
#' @param ... Other arguments to \code{\link{get_lambda_GC}} or
#'     \code{\link{get_annotations}}.
#'
#' @return A big_SVD object.
#'
#' @import bigsnpr
#' @import bigstatsr
#' @import ggplot2
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate case_when between summarise_all mutate_at
#' @importFrom tibble enframe as_tibble tibble
#' @importFrom rlang .data
#' @importFrom readr write_rds
#' @importFrom stats p.adjust
#' @importFrom tidyselect all_of
#' @importFrom cowplot save_plot
#'
#' @examples
#' \dontrun{
#' cdbn_standard_gwas(snp, df = phenotypes, type = "linear", covar = svd,
#'     ncores = nb_cores(), lambdagc = TRUE, savegwas = TRUE, saveplots = TRUE,
#'     saveannos = TRUE, txdb = txdb)
#' }
#'
#' @export
cdbn_standard_gwas <- function(snp, df = CDBNgenomics::BLUPs[,c(1,10,15,28)],
                                type = c("linear", "logistic"),
                                ncores = nb_cores(),
                                outputdir = ".", covar = NULL, lambdagc = TRUE,
                                savegwas = FALSE, saveplots = TRUE,
                                saveannos = FALSE, txdb = NULL, minphe = 200,
                                ...){

  if(attr(snp, "class") != "bigSNP"){
    stop("snp needs to be a bigSNP object, produced by the bigsnpr package.")
  }
  if(colnames(df)[1] != "Taxa"){
    stop("First column of phenotype dataframe (df) must be 'Taxa'.")
  }
  stopifnot(type %in% c("linear", "logistic"))
  if(saveannos == TRUE & is.null(txdb)){
    stop(paste0("Need to specify a txdb object created using AnnotationDbi ",
                "in order to generate data frames containing annotated top ",
                "SNPs. If you don't have this, set saveannos = FALSE."))
  }
  if(is.null(colnames(lambdagc))){
    message(paste0("'lambdagc' is TRUE, so lambda_GC will be used to find ",
                   "the best population structure correction using the ",
                   "covariance matrix."))
  } else if(!(colnames(lambdagc)[1] == "NumPCs" &
              colnames(lambdagc)[2] %in% colnames(df))){
    stop(paste0("If lambdagc is a dataframe, the column names must include",
                " NumPCs and the names of the phenotypes to run GWAS on ",
                "(the names of 'df')."))
  }
  if(savegwas == FALSE){
    message(paste0("'savegwas' is FALSE, so the gwas results will not be ",
                   "saved to disk."))
  }

  nSNP_M <- round(snp$genotypes$ncol/1000000, digits = 3)
  nInd <- snp$genotypes$nrow
  if(is.null(covar)){
    message(paste0("Covariance matrix (covar) was not supplied - this will be",
                   " generated using get_SVD()."))
    requireNamespace("dots")
    k <- dots::dots(name = 'k', value = 15, ...)
    covar <- get_SVD(snp, k = k, ncores = ncores, saveoutput = FALSE)
    if(savegwas == TRUE){
      saveRDS(covar, file = file.path(outputdir, paste0("SVD_", nInd, "g_",
                                                        nSNP_M, "M_SNPs_",
                                                        k, "PCs.rds")))
    }
  } else {
    stopifnot(attr(covar, "class") == "big_SVD")
  }

  plants <- snp$fam$sample.ID
  bonferroni <- -log10(0.05/length(snp$map$physical.pos))
  markers <- tibble(CHR = snp$map$chromosome, POS = snp$map$physical.pos)
  df <- plants %>%
    enframe(name = NULL, value = "Taxa") %>%
    left_join(df, by = "Taxa") %>%
    group_by(.data$Taxa) %>%
    summarise_all( ~ mean(., na.rm = TRUE)) %>%
    mutate_at(.vars = vars(-.data$Taxa), .funs = ~replace(., is.nan(.), NA))

  for(i in 2:ncol(df)){

    df1 <- df %>%
      dplyr::select(.data$Taxa, all_of(i))
    phename <- names(df1)[2]
    nPhe <- length(which(!is.na(df1[,2])))
    nLev <- nrow(unique(df1[which(!is.na(df1[,2])),2]))
    # Checks for correct combinations of phenotypes and GWAS types.
    if(nPhe < minphe){
      message(paste0("The phenotype ", phename, " does not have the minimum ",
                     "number of phenotyped Taxa's, (", minphe, ") and so ",
                     "will not be used for GWAS."))
      next
    } else if(nLev < 2){
      message(paste0("The phenotype ", phename, " does not have two or more ",
                     "distinct non-NA values and will not be used for GWAS."))
      next
    } else if(nLev > 2 & type == "logistic"){
      message(paste0("The phenotype ", phename, " has more than two distinct ",
                     "non-NA values and will not be used for GWAS with 'type=",
                     "logistic'."))
      next
    } else if(!(unique(df1[which(!is.na(df1[,2])),2])[1,1] %in% c(0,1)) &
              !(unique(df1[which(!is.na(df1[,2])),2])[2,1] %in% c(0,1)) & type == "logistic"){
      message(paste0("The phenotype ", phename, " has non-NA values that are ",
                     "not 0 or 1 and will not be used for GWAS with 'type=",
                     "logistic'."))
      next
    } else {
      message(paste0("Now starting GWAS pipeline for ", phename, "."))

      if(nLev == 2 & type == "linear"){
        message(paste0("The phenotype ", phename, " has only two distinct non-",
                       "NA values; consider using a logistic model instead.",
                       "(Set type = 'logistic')."))
      }

      if(is.null(colnames(lambdagc))){
        pc_max = ncol(covar$u)
        message(paste0("Now determining lambda_GC for GWAS models with ",
                       pc_max+1, " sets of PCs. This will take some time."))
        lambdagc_df <- get_lambda_GC(df = df1, type = type, snp = snp,
                                       covar = covar, ncores = ncores,
                                       npcs = c(0:pc_max), saveoutput = FALSE)
        if(saveplots == TRUE){
          write_csv(lambdagc_df, path = file.path(outputdir,
                                                  paste0("Lambda_GCs_by_PC_Num_",
                                                         phename, "_", type,
                                                         "_model_", nPhe, "g_",
                                                         nSNP_M, "M_SNPs",
                                                         ".csv")))
        }
        PCdf <- get_best_PC_df(lambdagc_df) # asv_best_PC_df(lambdagc_df)
      } else {
        PC1 <- lambdagc %>%
          dplyr::select(.data$NumPCs, phename)
        PCdf <- get_best_PC_df(PC1) # asv_best_PC_df(PC1)
      }
      PCdf1 <- PCdf[1,]

      # ------ Run GWAS with best pop structure correction -----
      message("Now running GWAS with the best population structure correction.")
      if(PCdf1$NumPCs == 0){
        gwas <- cdbn_gwas(df = df1, type = type, snp = snp, ncores = ncores)
      } else {
        gwas <- cdbn_gwas(df = df1, type = type, snp = snp, covar = covar,
                           ncores = ncores, npcs = PCdf1$NumPCs)
      }

      gwas_data <- tibble(CHR = markers$CHR, POS = markers$POS,
                          estim = gwas$estim, std_err = gwas$std.err,
                          bigsnpscore = gwas$score,
                          pvalue = predict(gwas, log10 = FALSE))
      gwas_data <- gwas_data %>%
        mutate(log10p = -log10(.data$pvalue))
      gwas_data$FDR_adj <- p.adjust(gwas_data$pvalue, method = "BH")
      if(savegwas == TRUE){
        # Save a data.table object with the GWAS results
        write_rds(gwas_data, path = file.path(outputdir,
                                              paste0("GWAS_datatable_", phename,
                                                     "_", type, "_model_", nPhe,
                                                     "g_", nSNP_M, "M_SNPs_",
                                                     PCdf1$NumPCs, "_PCs_",
                                                     "_.rds")), compress = "gz")
      }
      if(saveplots == TRUE){
        message("Now generating and saving Manhattan and QQ plots.")
        # ggplot settings for Manhattans and QQ plots.
        theme_oeco <- theme_classic() +
          theme(axis.title = element_text(size = 10),
                axis.text = element_text(size = 10),
                axis.line.x = element_line(size = 0.35, colour = 'grey50'),
                axis.line.y = element_line(size = 0.35, colour = 'grey50'),
                axis.ticks = element_line(size = 0.25, colour = 'grey50'),
                legend.justification = c(1, 0.75), legend.position = c(1, 0.9),
                legend.key.size = unit(0.35, 'cm'),
                legend.title = element_blank(),
                legend.text = element_text(size = 9),
                legend.text.align = 0, legend.background = element_blank(),
                plot.subtitle = element_text(size = 10, vjust = 0),
                strip.background = element_blank(),
                strip.text = element_text(hjust = 0.5, size = 10 ,vjust = 0),
                strip.placement = 'outside', panel.spacing.x = unit(-0.5, 'cm'))
        # Find 5% FDR threshold
        FDRthreshhi <- gwas_data %>%
          as_tibble() %>%
          filter(between(.data$FDR_adj, 0.10001, 1)) %>%  # 10% FDR threshold
          summarise(thresh = max(.data$log10p))
        FDRthreshlo <- gwas_data %>%
          as_tibble() %>%
          filter(between(.data$FDR_adj, 0, 0.09999)) %>%  # 10% FDR threshold
          summarise(thresh = min(.data$log10p))
        if(FDRthreshhi$thresh[1] > 0 & FDRthreshlo$thresh[1] > 0){
          FDRthreshold = (FDRthreshhi$thresh[1] + FDRthreshlo$thresh[1])/2
        } else if(FDRthreshhi$thresh[1] > 0){
          FDRthreshold = FDRthreshhi$thresh[1]
        } else if(FDRthreshlo$thresh[1] > 0){
          FDRthreshold = FDRthreshlo$thresh[1]
        } else {
          FDRthreshold = NA
        }
        # Save a Manhattan plot with 5% FDR
        ggmanobject1 <- gwas_data %>%
          filter(.data$log10p > 1) %>%
          ggplot(aes(x = .data$POS, y = .data$log10p)) +
          theme_oeco +
          geom_hline(yintercept = c(5, 10), color = "lightgrey") +
          geom_point(aes(color = as.factor(.data$CHR), fill = as.factor(.data$CHR))) +
          geom_hline(yintercept = FDRthreshold, color = "black", linetype = 2,
                     size = 1) +
          facet_wrap(~ .data$CHR, nrow = 1, scales = "free_x",
                     strip.position = "bottom") +
          scale_color_manual(values = rep(c("#1B0C42FF", "grey"), 9),
                             guide = FALSE) +
          theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                panel.background = element_rect(fill=NA),
                legend.position = "none") +
          labs(x = "Chromosome", y = "-log10(p value)") +
          scale_x_continuous(expand = c(0.18, 0.18))

        save_plot(filename = file.path(outputdir,
                                       paste0("Manhattan_", phename, "_", type,
                                              "_model_", nPhe, "g_", nSNP_M,
                                              "M_SNPs_", PCdf1$NumPCs,
                                              "_PCs_5percent_FDR_",
                                              get_date_filename(), ".png")),
                  plot = ggmanobject1, base_asp = 2.5, base_height = 3.5)

        # Save a QQplot
        ggqqplot <- get_qqplot(ps = gwas_data$pvalue, lambdaGC = TRUE)
        save_plot(filename = file.path(outputdir,
                                       paste0("QQplot_", phename, "_", type,
                                              "_model_", nPhe, "g_", nSNP_M,
                                              "M_SNPs_", PCdf1$NumPCs,
                                              "_PCs_FDR_",
                                              get_date_filename(), ".png")),
                  plot = ggqqplot + theme_oeco, base_asp = 1, base_height = 3.5)

        # Save a Manhattan plot with Bonferroni
        ggmanobject2 <- gwas_data %>%
          filter(.data$log10p > 1) %>%
          ggplot(aes(x = .data$POS, y = .data$log10p)) +
          theme_oeco +
          geom_hline(yintercept = c(5, 10), color = "lightgrey") +
          geom_point(aes(color = as.factor(.data$CHR), fill = as.factor(.data$CHR))) +
          geom_hline(yintercept = bonferroni, color = "black", linetype = 2,
                     size = 1) +
          facet_wrap(~ .data$CHR, nrow = 1, scales = "free_x",
                     strip.position = "bottom") +
          scale_color_manual(values = rep(c("#1B0C42FF", "grey"), 9),
                             guide = FALSE) +
          theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                panel.background = element_rect(fill=NA),
                legend.position = "none") +
          labs(x = "Chromosome", y = "-log10(p value)") +
          scale_x_continuous(expand = c(0.18, 0.18))

        save_plot(filename = file.path(outputdir,
                                       paste0("Manhattan_", phename, "_", type,
                                              "_model_", nPhe, "g_", nSNP_M,
                                              "M_SNPs_", PCdf1$NumPCs,
                                              "_PCs_Bonferroni_",
                                              get_date_filename(), ".png")),
                  plot = ggmanobject2, base_asp = 2.5, base_height = 3.5)
      }

      if(saveannos == TRUE){
        message(paste0("Now creating annotation data frames for the top 10 & ",
                       "top 500 SNPs by p-value, and for SNPs above a 10% FDR."))
        ## Save annotation tables for the top associations
        requireNamespace("dots")
        n <- dots::dots(name = 'n', value = c(10, 500), ...)
        FDRalpha <- dots::dots(name = 'FDRalpha', value = NA, ...)
        rangevector <- dots::dots(name = 'rangevector', value = c(0, 50000), ...)
        anno_tables <- get_annotations(df = gwas, type = "bigsnp",
                                       n = n, FDRalpha = FDRalpha,
                                       rangevector = rangevector,
                                       markers = markers,
                                       anno_info = CDBNgenomics::anno_info,
                                       txdb = txdb)
        saveRDS(anno_tables, file.path(outputdir,
                                       paste0("Annotation_tables_", phename,
                                              "_", type, "_model_", nPhe, "g_",
                                              nSNP_M, "M_SNPs_",PCdf1$NumPCs,
                                              "_PCs", ".rds")))
      }
    }
  }
}

#' Return lambda_GC for different numbers of PCs for GWAS on any bigsnp file.
#'
#' @description Given a dataframe of phenotypes associated with Taxa and
#'     output from a PCA to control for population structure, this function will
#'     return a .csv file of the lambda_GC values for the GWAS upon inclusion
#'     of different numbers of PCs. This allows the user to choose a number of
#'     PCs that returns a lambda_GC close to 1, and thus ensure that they have
#'     done adequate correction for population structure.
#'
#' @param df Dataframe of phenotypes where the first column is Taxa and each
#'     Taxa occurs only once in the dataframe.
#' @param type Character string. Type of univarate regression to run for GWAS.
#'     Options are "linear" or "logistic".
#' @param snp Genomic information to include for Phaseolus vulgaris.
#' @param covar Covariance matrix to include in the regression. You
#'     can generate these using \code{bigsnpr::snp_autoSVD()}.
#' @param ncores Number of cores to use. Default is one.
#' @param npcs Integer vector of principle components to use.
#'     Defaults to c(0:10).
#' @param saveoutput Logical. Should output be saved as a csv to the
#'     working directory?
#'
#' @import bigsnpr
#' @import bigstatsr
#' @importFrom dplyr mutate rename case_when mutate_if
#' @importFrom purrr as_vector
#' @importFrom tibble as_tibble enframe
#' @importFrom rlang .data
#' @importFrom readr write_csv
#' @importFrom utils tail
#'
#' @return A dataframe containing the lambda_GC values for each number of PCs
#'     specified. This is also saved as a .csv file in the working directory.
get_lambda_GC <- function(df, type = c("linear", "logistic"), snp,
                            covar = NA, ncores = 1, npcs = c(0:10),
                            saveoutput = FALSE){
  if(attr(snp, "class") != "bigSNP"){
    stop("snp needs to be a bigSNP object, produced by the bigsnpr package.")
  }
  if(colnames(df)[1] != "Taxa"){
    stop("First column of phenotype dataframe (df) must be 'Taxa'.")
  }
  if(length(covar) == 1){
    stop(paste0("Need to specify covariance matrix (covar) and a vector of",
                " PC #'s to test (npcs)."))
  }
  if(saveoutput == FALSE){
    message("saveoutput is FALSE, so lambda_GC values won't be saved to a csv.")
  }

  G <- snp$genotypes

  LambdaGC <- as_tibble(matrix(data = c(npcs,
                                        rep(NA, (ncol(df) - 1)*length(npcs))),
                               nrow = length(npcs), ncol = ncol(df),
                               dimnames = list(npcs, colnames(df))))
  LambdaGC <- LambdaGC %>%
    dplyr::rename("NumPCs" = .data$Taxa) %>%
    mutate_if(is.integer, as.numeric)

  for(i in seq_along(names(df))[-1]){

    for(k in c(1:length(npcs))){

      if(type == "linear"){

        y1 <- as_vector(df[which(!is.na(df[,i])), i])
        ind_y <- which(!is.na(df[,i]))

        if(npcs[k] == 0){

          gwaspc <- big_univLinReg(G, y.train = y1, ind.train = ind_y,
                                   ncores = ncores)
        } else {

          ind_u <- matrix(covar$u[which(!is.na(df[,i])),1:npcs[k]],
                          ncol = npcs[k])
          gwaspc <- big_univLinReg(G, y.train = y1, covar.train = ind_u,
                                   ind.train = ind_y, ncores = ncores)
        }
      } else if(type == "logistic"){
        message(paste0("For logistic models, if convergence is not reached by ",
                       "the main algorithm for some SNPs, the corresponding `niter` element ",
                       "is set to NA, and glm is used instead. If glm can't ",
                       "converge either, those SNP estimations are set to NA."))
        y1 <- as_vector(df[which(!is.na(df[,i])), i])
        ind_y <- which(!is.na(df[,i]))
        if(npcs[k] == 0){
          gwaspc <- suppressMessages(big_univLogReg(G, y01.train = y1,
                                                    ind.train = ind_y,
                                                    ncores = ncores))
        } else {
          ind_u <- matrix(covar$u[which(!is.na(df[,i])),1:npcs[k]],
                          ncol = npcs[k])
          gwaspc <- suppressMessages(big_univLogReg(G, y01.train = y1,
                                                    covar.train = ind_u,
                                                    ind.train = ind_y,
                                                    ncores = ncores))
        }
      }
      ps <- predict(gwaspc, log10 = FALSE)
      LambdaGC[k,i] <- get_lambdagc(ps = ps)
      message(paste0("Finished Lambda_GC calculation for ", names(df)[i],
                     " using ", npcs[k], " PCs."))
    }

    if(saveoutput == TRUE){
      write_csv(LambdaGC, path = paste0("Lambda_GC_", names(df)[i], ".csv"))
    }
    message(paste0("Finished phenotype ", i-1, ": ", names(df)[i]))
  }
  if(saveoutput == TRUE){
    write_csv(LambdaGC, path = paste0("Lambda_GC_", names(df)[2], "_to_",
                                      tail(names(df), n = 1), "_Phenotypes_",
                                      npcs[1], "_to_", tail(npcs, n = 1),
                                      "_PCs.csv"))
    best_LambdaGC <- get_best_PC_df(df = LambdaGC)
    write_csv(best_LambdaGC, path = paste0("Best_Lambda_GC_", names(df)[2],
                                           "_to_", tail(names(df), n = 1),
                                           "_Phenotypes_", npcs[1], "_to_",
                                           tail(npcs, n = 1), "_PCs.csv"))
  }
  return(LambdaGC)
}

#' Return best number of PCs in terms of lambda_GC for the CDBN.
#'
#' @description Given a dataframe created using get_lambda_GC, this function
#'     returns the first lambda_GC less than 1.05, or the smallest lambda_GC,
#'     for each column in the dataframe.
#'
#' @param df Dataframe of phenotypes where the first column is NumPCs and
#'     subsequent column contains lambda_GC values for some phenotype.
#'
#' @importFrom dplyr filter top_n select full_join arrange
#' @importFrom tidyr gather
#' @importFrom rlang .data sym !!
#' @importFrom tidyselect all_of
#'
#' @return A dataframe containing the best lambda_GC value and number of PCs
#'     for each phenotype in the data frame.
get_best_PC_df <- function(df){
  column <- names(df)[ncol(df)]
  bestPCs <- df %>%
    filter(!! sym(column) < 1.05| !! sym(column) == min(!! sym(column))) %>%
    top_n(n = -1, wt = .data$NumPCs) %>%
    select(.data$NumPCs, all_of(column))

  if(ncol(df) > 2){
    for(i in c((ncol(df)-2):1)){
      column <- names(df)[i+1]

      bestPCs <- df %>%
        filter(!! sym(column) < 1.05 | !! sym(column) == min(!! sym(column))) %>%
        top_n(n = -1, wt = .data$NumPCs) %>%
        select(.data$NumPCs, all_of(column)) %>%
        full_join(bestPCs, by = c("NumPCs", (column)))
    }
  }

  bestPCdf <- bestPCs %>%
    arrange(.data$NumPCs) %>%
    gather(key = "trait", value = "lambda_GC", 2:ncol(bestPCs)) %>%
    filter(!is.na(.data$lambda_GC))

  return(bestPCdf)
}


#' Wrapper for bigsnpr for GWAS on Phaseolus vulgaris (the CDBN).
#'
#' @description Given a dataframe of phenotypes associated with Taxa, this
#'     function is a wrapper around bigsnpr functions to conduct linear or
#'     logistic regression on Phaseolus vulgaris. The main advantages of this
#'     function over just using the bigsnpr functions is that it automatically
#'     removes individual genotypes with missing phenotypic data
#'     and that it can run GWAS on multiple phenotypes sequentially.
#'
#' @param df Dataframe of phenotypes where the first column is Taxa.
#' @param type Character string. Type of univarate regression to run for GWAS.
#'     Options are "linear" or "logistic".
#' @param snp Genomic information to include for Phaseolus vulgaris.
#' @param covar Optional covariance matrix to include in the regression. You
#'     can generate these using \code{bigsnpr::snp_autoSVD()}.
#' @param ncores Number of cores to use. Default is one.
#' @param npcs Number of principle components to use. Default is 10.
#' @param saveoutput Logical. Should output be saved as a rds to the
#'     working directory?
#'
#' @import bigsnpr
#' @import bigstatsr
#' @importFrom dplyr mutate rename case_when
#' @importFrom purrr as_vector
#' @importFrom tibble as_tibble enframe
#' @importFrom rlang .data
#'
#' @return The gwas results for the last phenotype in the dataframe. That
#'     phenotype, as well as the remaining phenotypes, are saved as RDS objects
#'     in the working directory.
cdbn_gwas <- function(df, type = c("linear", "logistic"), snp,
                       covar = NA, ncores = 1, npcs = 10, saveoutput = FALSE){
  stopifnot(type %in% c("linear", "logistic"))
  if(attr(snp, "class") != "bigSNP"){
    stop("snp needs to be a bigSNP object, produced by the bigsnpr package.")
  }
  if(colnames(df)[1] != "Taxa"){
    stop("First column of phenotype dataframe (df) must be 'Taxa'.")
  }
  G <- snp$genotypes

  for(i in seq_along(names(df))[-1]){
    y1 <- as_vector(df[which(!is.na(df[,i])), i])
    ind_y <- which(!is.na(df[,i]))

    if(type == "linear"){
      if(!is.na(covar[1])){
        ind_u <- matrix(covar$u[which(!is.na(df[,i])),1:npcs], ncol = npcs)
        gwaspc <- big_univLinReg(G, y.train = y1, covar.train = ind_u,
                                 ind.train = ind_y, ncores = ncores)
      } else {
        gwaspc <- big_univLinReg(G, y.train = y1, ind.train = ind_y,
                                 ncores = ncores)
      }
    } else if(type == "logistic"){
      message(paste0("For logistic models, if convergence is not reached by ",
                     "the main algorithm for any SNP, the corresponding `niter` element ",
                     "is set to NA, and glm is used instead. If glm can't ",
                     "converge either, those SNP estimations are set to NA."))
      if(!is.na(covar[1])){
        ind_u <- matrix(covar$u[which(!is.na(df[,i])),1:npcs], ncol = npcs)
        gwaspc <- suppressMessages(big_univLogReg(G, y01.train = y1,
                                                  covar.train = ind_u,
                                                  ind.train = ind_y,
                                                  ncores = ncores))
      } else {
        gwaspc <- suppressMessages(big_univLogReg(G, y01.train = y1,
                                                  ind.train = ind_y,
                                                  ncores = ncores))
      }
    } else {
      stop(paste0("Type of GWAS not recognized: please choose one of 'linear'",
                  " or 'logistic'"))
    }

    if(saveoutput){
      saveRDS(gwaspc, file = paste0("GWAS_object_", names(df)[i], ".rds"))
    } else {
      print("saveoutput is FALSE so GWAS object will not be saved to disk.")
    }

  }
  return(gwaspc)
}



#' Return a number rounded to some number of digits
#'
#' @description Given some x, return the number rounded to some number of
#'     digits.
#'
#' @param x A number or vector of numbers
#' @param at Numeric. Rounding factor or size of the bin to round to.
#'
#' @return A number or vector of numbers
round2 <- function(x, at) ceiling(x / at) * at

#' Return a dataframe binned into 2-d bins by some x and y.
#'
#' @description Given a dataframe of x and y values (with some optional
#'     confidence intervals surrounding the y values), return only the unique
#'     values of x and y in some set of 2-d bins.
#'
#' @param x Numeric vector. The first vector for binning.
#' @param y Numeric vector. the second vector for binning
#' @param cl Numeric vector. Optional confidence interval for the y vector,
#'     lower bound.
#' @param cu Numeric vector. Optional confidence interval for the y vector,
#'     upper bound.
#' @param roundby Numeric. The amount to round the x and y vectors by for 2d
#'     binning.
#'
#' @return A dataframe containing the 2-d binned values for x and y, and their
#'     confidence intervals.
round_xy <- function(x, y, cl = NA, cu = NA, roundby = 0.001){
  expected <- round2(x, at = roundby)
  observed <- round2(y, at = roundby)
  if(!is.na(cl[1]) & !is.na(cu[1])){
    clower <- round2(cl, at = roundby)
    cupper <- round2(cu, at = roundby)
    tp <- cbind(expected, observed, clower, cupper)
    return(tp[!duplicated(tp),])
  } else {
    tp <- cbind(expected, observed)
    return(tp[!duplicated(tp),])
  }
}

#' Create a quantile-quantile plot with ggplot2.
#'
#' @description Assumptions for this quantile quantile plot:
#'     Expected P values are uniformly distributed.
#'     Confidence intervals assume independence between tests.
#'     We expect deviations past the confidence intervals if the tests are
#'     not independent.
#'     For example, in a genome-wide association study, the genotype at any
#'     position is correlated to nearby positions. Tests of nearby genotypes
#'     will result in similar test statistics.
#'
#' @param ps Numeric vector of p-values.
#' @param ci Numeric. Size of the confidence interval, 0.95 by default.
#' @param lambdaGC Logical. Add the Genomic Control coefficient as subtitle to
#'     the plot?
#'
#' @import ggplot2
#' @importFrom tibble as_tibble
#' @importFrom rlang .data
#' @importFrom stats qbeta ppoints
#' @param tol Numeric. Tolerance for optional Genomic Control coefficient.
#'
#' @return A ggplot2 plot.
get_qqplot <- function(ps, ci = 0.95, lambdaGC = FALSE, tol = 1e-8) {
  ps <- ps[which(!is.na(ps))]
  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  df_round <- round_xy(df$expected, df$observed, cl = df$clower, cu = df$cupper)
  log10Pe <- expression(paste("Expected -log"[10], plain("("), italic(p-value),
                              plain(")")))
  log10Po <- expression(paste("Observed -log"[10], plain("("), italic(p-value),
                              plain(")")))
  p1 <- ggplot(as_tibble(df_round)) +
    geom_point(aes(.data$expected, .data$observed), shape = 1, size = 1) +
    geom_abline(intercept = 0, slope = 1, size = 1.5, color = "red") +
    geom_line(aes(.data$expected, .data$cupper), linetype = 2) +
    geom_line(aes(.data$expected, .data$clower), linetype = 2) +
    xlab(log10Pe) +
    ylab(log10Po)

  if (lambdaGC) {
    lamGC <- get_lambdagc(ps = ps, tol = tol)
    expr <- substitute(expression(lambda[GC] == l), list(l = lamGC))
    p1 + labs(subtitle = eval(expr))
  } else {
    p1
  }
}

#' Find lambda_GC value for non-NA p-values
#'
#' @description Finds the lambda GC value for some vector of p-values.
#'
#' @param ps Numeric vector of p-values. Can have NA's.
#' @param tol Numeric. Tolerance for optional Genomic Control coefficient.
#'
#' @importFrom stats median uniroot
#'
#' @return A lambda GC value (some positive number, ideally ~1)
#'
#' @export
get_lambdagc <- function(ps, tol = 1e-8){
  ps <- ps[which(!is.na(ps))]
  xtr <- log10(ps)
  MEDIAN <- log10(0.5)
  f.opt <- function(x) (x - MEDIAN)
  xtr_p <- median(xtr) / uniroot(f.opt, interval = range(xtr),
                                 check.conv = TRUE,
                                 tol = tol)$root
  lamGC <- signif(xtr_p)
  return(lamGC)
}
