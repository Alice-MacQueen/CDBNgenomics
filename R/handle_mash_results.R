
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
#' @export
get_colnames <- function(m){
  column_names <- colnames(m$result$lfsr)
  return(column_names)
}

#' @title Get mash marker_df
#'
#' @description Pulls the names of the SNP markers from the mash object.
#'
#' @param m An object of type mash
#'
#' @export
get_marker_df <- function(m){
  marker_df <- get_significant_results(m, thresh = 1) %>%
    enframe(name = "Marker") %>%
    arrange(value)

  return(marker_df)
}

