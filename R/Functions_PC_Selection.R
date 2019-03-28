#!/usr/bin/env Rscript

#' Identify trait names for GAPIT results.
#'
#' Creates a vector of phenotype names from GAPIT results in some directory.
#'
#' @param path File path to the csv files that GAPIT has created, a character
#' string or \code{file.path()}.
#' @param model Model type used in GAPIT runs, a character string.
#'
#' @return A vector of phenotype names.
#'
#' @examples
#' gapit_results_in_filepath(path = file.path("inst", "extdata"))
#' gapit_results_in_filepath(path = "./inst/extdata")
#'
#' @export
gapit_results_in_filepath <- function(path = "."){
  result_files <- list.files(path = path, pattern = "*GWAS.Results*")
  phemiddle <- purrr::partial(stringr::str_sub, start = 7, end = -18)
  gapit_phenotypes <- purrr::map(result_files, phemiddle) %>%
    unlist()
  return(gapit_phenotypes)
}

#' Select the best number of PCs using BIC.
#'
#' If you choose \code{Model.Selection = TRUE} in GAPIT, you will get
#' BIC.Model.Selection output files. If you have many traits, it can be tedious
#' to manually go through each output file and select the number of principle
#' components that maximizes BIC. This function automates that process for many
#' traits in some filepath. It returns a data frame containing the trait names
#' in a filepath and the best number of PC's, a.k.a. the number that maximizes
#' the BIC for that trait.
#'
#' @param path File path to the csv Result files that GAPIT has created, a
#' character string or \code{file.path()}.
#' @param saveoutput Should the output of this function be saved to disc?
#' Recommended to be set to TRUE, but the default is FALSE.
#'
#' @return A tbl_df table of traits with the best number of PC's as determined
#' by the maximum BIC.
#'
#' @examples
#' \dontrun{gapit_model_selection(path = "./inst/extdata", saveoutput = TRUE)}
#' gapit_model_selection(path = "./inst/extdata")
#'
#' @export
gapit_model_selection <- function(path, saveoutput = FALSE){
  bic_files <- list.files(path = path, pattern = "*BIC.Model.Selection*")
  phemiddle <- purrr::partial(stringr::str_sub, start = 7, end = -29)
  traits <- purrr::map(bic_files, phemiddle) %>% unlist()
  PCnum <- rep(NA, length(traits))
  for(i in seq_along(traits)){
    filename <- eval(bic_files[i])
    bic_df <- readr::read_csv(file.path(path, filename),
                              col_names = TRUE, col_types = "nnn")
    PCnum[i] <- as_vector(bic_df[which(bic_df[, 2] == max(bic_df[, 2])), 1])
  }
  PCselect <- tibble::as_tibble(list(traits = traits, bestPCnum = PCnum))
  if(saveoutput == TRUE){
    readr::write_csv(PCselect, paste0("PC_selection_via_BIC_", Sys.Date(),
                                      ".csv"))
  }
  return(PCselect)
}

