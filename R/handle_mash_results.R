#' Get number of conditions
#'
#' @param m The mash result
#'
#' @importFrom ashr get_pm
#'
get_ncond = function(m){
  return(ncol(get_pm(m)))
}

#' @title Get column names from a mash object
#'
#' @description This function extracts the column names from the local false
#' sign rate table of a mash object's results. This can tell you the condition
#' names or phenotype names used in the mash object. That can be useful for
#' looking at a subset of these columns, say.
#'
#' @param m An object of type mash
#'
#' @return A vector of phenotype names
#'
#' @examples
#'     \dontrun{get_colnames(m = mash_obj)}
#'
get_colnames <- function(m){
  column_names <- colnames(m$result$lfsr)
  return(column_names)
}


#' From a mash result, get effects that are significant in at least one condition
#'
#' @param m the mash result (from joint or 1by1 analysis)
#' @param thresh indicates the threshold below which to call signals significant
#' @param conditions which conditions to include in check (default to all)
#' @param sig_fn the significance function used to extract significance from mash object; eg could be ashr::get_lfsr or ashr::get_lfdr. (Small values must indicate significant.)
#'
#' @return a vector containing the indices of the significant effects, by order of most significant to least
#'
#' @importFrom ashr get_lfsr
#'
#' @export
get_significant_results = function(m, thresh = 0.05, conditions = NULL,
                                   sig_fn = ashr::get_lfsr) {
  if (is.null(conditions)) {
    conditions = 1:get_ncond(m)
  }
  top = apply(sig_fn(m)[,conditions,drop=FALSE],1,min) # find top effect in each condition
  sig = which(top < thresh)
  ord = order(top[sig],decreasing=FALSE)
  sig[ord]
}


#' @title Get mash marker_df
#'
#' @description Pulls the names of the SNP markers from the mash object.
#'
#' @param m An object of type mash
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom tibble enframe
#' @importFrom dplyr arrange
#'
get_marker_df <- function(m){
  marker_df <- get_significant_results(m, thresh = 1) %>%
    enframe(name = "Marker") %>%
    arrange(.data$value)

  return(marker_df)
}

