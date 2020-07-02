# --- Get Results from Mash Object ---------

#' @title Return the estimated mixture proportions. Use get_estimated_pi to
#'     extract the estimates of the mixture proportions for different types of
#'     covariance matrix. This tells you which covariance matrices have most of
#'     the mass.
#'
#' @param m the mash result
#' @param dimension indicates whether you want the mixture proportions for the
#'     covariances, grid, or all
#'
#' @return a named vector containing the estimated mixture proportions.
#'
#' @details If the fit was done with `usepointmass=TRUE` then the first
#'     element of the returned vector will correspond to the null, and the
#'     remaining elements to the non-null covariance matrices. Suppose the fit
#'     was done with $K$ covariances and a grid of length $L$. If
#'     `dimension=cov` then the returned vector will be of length $K$
#'     (or $K+1$ if `usepointmass=TRUE`).  If `dimension=grid` then
#'     the returned vector will be of length $L$ (or $L+1$).  If
#'     `dimension=all` then the returned vector will be of length $LK$ (or
#'     $LK+1$). The names of the vector will be informative for which
#'     combination each element corresponds to.
#'
#' @importFrom ashr get_fitted_g
#'
get_estimated_pi = function(m, dimension = c("cov","grid","all")){
  dimension = match.arg(dimension)
  if(dimension=="all"){
    get_estimated_pi_no_collapse(m)
  } else {
    g = get_fitted_g(m)
    pihat = g$pi
    pihat_names = NULL
    pi_null = NULL

    if(g$usepointmass){
      pihat_names=c("null",pihat_names)
      pi_null = pihat[1]
      pihat = pihat[-1]
    }

    pihat = matrix(pihat,nrow=length(g$Ulist))
    if(dimension=="cov"){
      pihat = rowSums(pihat)
      pihat_names = c(pihat_names,names(g$Ulist))
    } else if(dimension=="grid"){
      pihat = colSums(pihat)
      pihat_names = c(pihat_names,1:length(g$grid))
    }

    pihat = c(pi_null,pihat)
    names(pihat) = pihat_names
    return(pihat)
  }
}

get_estimated_pi_no_collapse = function(m){
  g = get_fitted_g(m)
  pihat = g$pi
  names(pihat) = names(expand_cov(g$Ulist, g$grid, g$usepointmass))
  pihat
}

#' @title Create expanded list of covariance matrices expanded by
#'   grid, Sigma_{lk} = omega_l U_k
#'
#' @description This is an internal (non-exported) function. This help
#'   page provides additional documentation mainly intended for
#'   developers and expert users.
#'
#' @param Ulist a list of covarance matrices
#'
#' @param grid a grid of scalar values by which the covariance
#'   matrices are to be sc
#'
#' @param usepointmass if TRUE adds a point mass at 0 (null component)
#'   to the list
#'
#' @return This takes the covariance matrices in Ulist and multiplies
#' them by the grid values If usepointmass is TRUE then it adds a null
#' component.
#'
#' @keywords internal
#'
expand_cov = function(Ulist,grid,usepointmass=TRUE){
  scaled_Ulist = scale_cov(Ulist, grid)
  R = nrow(Ulist[[1]])
  if(usepointmass){
    scaled_Ulist = c(list(null=matrix(0,nrow=R,ncol=R)),scaled_Ulist)
  }
  return(scaled_Ulist)
}

#' @title Scale each covariance matrix in list Ulist by a scalar in
#' vector grid
#'
#' @description This is an internal (non-exported) function. This help
#'   page provides additional documentation mainly intended for
#'   developers and expert users.
#'
#' @param Ulist a list of matrices
#'
#' @param grid a vector of scaling factors (standard deviaions)
#'
#' @return a list with length length(Ulist)*length(grid)
#'
#' @keywords internal
#'
scale_cov = function(Ulist, grid){
  orig_names = names(Ulist)
  Ulist = unlist( lapply(grid^2, function(x){multiply_list(Ulist,x)}), recursive=FALSE)
  names(Ulist) = unlist( lapply(1:length(grid), function(x){paste0(orig_names,".",x)}), recursive=FALSE)
  return(Ulist)
}

# Multiply each element of a list by scalar. (In our application each
# element of the list is a matrix.)
multiply_list = function(Ulist, x){lapply(Ulist, function(U){x*U})}


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

#' Return the Bayes Factor for each effect
#'
#' @param m the mash result (from joint or 1by1 analysis); must have been computed using usepointmass=TRUE
#'
#' @return if m was fitted using usepointmass=TRUE then returns a vector of
#'     the log10(bf) values for each effect. That is, the jth element
#'     lbf_j is log10(Pr(Bj | g=ghat-nonnull)/Pr(Bj | g = 0)) where gha
#'     t-nonnull is the non-null part of ghat.  Otherwise returns NULL.
#'
get_log10bf = function(m) {
  if(is.null(m$null_loglik)){
    return(NULL)
  } else {
    return((m$alt_loglik - m$null_loglik)/log(10))
  }
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

#' Get number of conditions
#'
#' @param m The mash result
#'
#' @importFrom ashr get_pm
#'
get_ncond = function(m){
  return(ncol(get_pm(m)))
}

#' Count number of conditions each effect is significant in
#'
#' @param m the mash result (from joint or 1by1 analysis)
#' @param thresh indicates the threshold below which to call signals significant
#' @param conditions which conditions to include in check (default to all)
#' @param sig_fn the significance function used to extract significance from mash object; eg could be ashr::get_lfsr or ashr::get_lfdr
#'
#' @return a vector containing the number of significant conditions
#'
get_n_significant_conditions = function(m, thresh = 0.05, conditions = NULL,
                                        sig_fn = get_lfsr){
  if (is.null(conditions)) {
    conditions = 1:get_ncond(m)
  }
  return(apply(sig_fn(m)[,conditions,drop=FALSE] < thresh, 1, sum))
}

#' Compute the proportion of (significant) signals shared by magnitude in each pair of conditions, based on the poterior mean
#'
#' @param m the mash fit
#' @param factor a number between 0 and 1 - the factor within which effects are
#'     considered to be shared.
#' @param lfsr_thresh the lfsr threshold for including an effect in the
#'     assessment
#' @param FUN a function to be applied to the estimated effect sizes before
#'     assessing sharing. The most obvious choice beside the default
#'     'FUN=identity' would be 'FUN=abs' if you want to ignore the sign of the
#'     effects when assesing sharing.
#' @details For each pair of tissues, first identify the effects that are
#'     significant (by lfsr<lfsr_thresh) in at least one of the two tissues.
#'     Then compute what fraction of these have an estimated (posterior mean)
#'     effect size within a factor `factor` of one another. The results are
#'     returned as an R by R matrix.
#'
#' @examples
#' \dontrun{
#' get_pairwise_sharing(m) # sharing by magnitude (same sign)
#' get_pairwise_sharing(m, factor=0) # sharing by sign
#' get_pairwise_sharing(m, FUN=abs) # sharing by magnitude when sign is ignored
#' }
#'
#' @export
get_pairwise_sharing = function(m, factor=0.5, lfsr_thresh=0.05, FUN= identity){
  R = get_ncond(m)
  lfsr = get_lfsr(m)
  S=matrix(NA,nrow = R, ncol=R)
  for(i in 1:R){
    for(j in i:R){
      sig_i=get_significant_results(m,thresh=lfsr_thresh,conditions = i)
      sig_j=get_significant_results(m,thresh=lfsr_thresh,conditions = j)
      a=union(sig_i,sig_j)
      ratio=FUN(get_pm(m)[a,i])/FUN(get_pm(m)[a,j])##divide effect sizes
      S[i,j]=mean(ratio>factor & ratio<(1/factor))
    }
  }
  S[lower.tri(S, diag = FALSE)] = t(S)[lower.tri(S, diag = FALSE)]
  colnames(S) = row.names(S) = colnames(m$result$PosteriorMean)

  return(S)
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

#' @title Get data frames of types of GxE from a mash object
#'
#' @description Performs set operations to determine pairwise GxE for effects
#'     from a mash object.
#'
#' @param m An object of type mash
#' @param thresh Numeric. The threshold for including an effect in the assessment
#' @param factor a number between 0 and 1. The factor within which effects are
#'     considered to be shared.
#'
#' @return A list containing eight data frames. Those with names that start
#'     "S_" contain significant effects of different types between pairs of
#'     named rows and columns. S_all_pairwise contains all significant effects;
#'     NS_pairwise contains all non-significant effects. S_CN contains effects
#'     significant in only one condition, and effects with a significantly
#'     different magnitude (differential sensitivity). This dataframe is not
#'     conservative using the local false sign rate test - we can't determine
#'     the sign of one of the effects for effects significant in only one
#'     condition - so it's not recommended to use this, but included. S_2_no
#'     contains effects significant in both conditions that do not differ
#'     significantly in magnitude. These effects do not have GxE. S_AP contains
#'     effects significant in both conditions that differ in their sign - and
#'     have antagonistic pleiotropy. S_DS contains effects significant in both
#'     conditions that differ in the magnitude of their effect, but not their
#'     sign - differentially sensitive alleles. S_1_row and S_1_col contain
#'     effects that are significant in just one of the two conditions - the row
#'     or the column, respectively.
#'
#' @importFrom dplyr between mutate filter
#' @importFrom tibble enframe
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @export
get_GxE = function(m, factor = 0.4, thresh = 0.05){
  R = get_ncond(m)                          # Effects to consider

  S_all = matrix(NA, nrow = R, ncol = R, dimnames = list(get_colnames(m),
                                                         get_colnames(m)))
  S_2_no = matrix(NA, nrow = R, ncol = R, dimnames = list(get_colnames(m),
                                                          get_colnames(m)))
  S_CN = matrix(NA, nrow = R, ncol = R, dimnames = list(get_colnames(m),
                                                        get_colnames(m)))
  S_AP = matrix(NA, nrow = R, ncol = R, dimnames = list(get_colnames(m),
                                                        get_colnames(m)))
  S_DS = matrix(NA, nrow = R, ncol = R, dimnames = list(get_colnames(m),
                                                        get_colnames(m)))
  NS_pair = matrix(NA, nrow = R, ncol = R, dimnames = list(get_colnames(m),
                                                          get_colnames(m)))
  S_i = matrix(NA, nrow = R, ncol = R, dimnames = list(get_colnames(m),
                                                           get_colnames(m)))
  S_j = matrix(NA, nrow = R, ncol = R, dimnames = list(get_colnames(m),
                                                           get_colnames(m)))

  for(i in 1:R){
    for(j in 1:R){

      if(i == j){
        S_all[i, j] = length(get_significant_results(m, thresh = thresh,
                                                     conditions = i))
        S_CN[i, j] = 0
        # Not conservative!!
        S_2_no[i, j] = length(get_significant_results(m, thresh = thresh,
                                                      conditions = i))
        S_AP[i, j] = 0
        S_DS[i, j] = 0

        sig_i = get_significant_results(m, thresh = thresh, conditions = i)
        all_i = get_significant_results(m, thresh = 1, conditions = i)
        ns_i = dplyr::setdiff(all_i, sig_i)   # effects that aren't sig in i
        NS_pair[i, j] = length(ns_i)
        S_i[i, j] = 0
        S_j[i, j] = 0

      } else {

        sig_i = get_significant_results(m, thresh = thresh, conditions = i)
        sig_j = get_significant_results(m, thresh = thresh, conditions = j)

        all_i = get_significant_results(m, thresh = 1, conditions = i)
        all_j = get_significant_results(m, thresh = 1, conditions = j)

        ns_i = setdiff(all_i, sig_i)   # effects that aren't sig in i
        ns_j = setdiff(all_j, sig_j)   # effects that aren't sig in j

        # Markers where we aren't sure of the sign in either condition
        # aka most of the effects
        NS_pair[i,j] <- length(intersect(ns_i, ns_j))

        # Markers where we are sure of the sign in just one condition

        # Markers significant in i but not in j
        ms_isigi = intersect(sig_i, ns_j)

        # Markers significant in j but not in i
        ms_jsigj = intersect(sig_j, ns_i)

        # Markers where we are sure of the sign in both conditions
        effects_df <- get_pm(m)[union(sig_i, sig_j), i] %>%
          enframe(name = "Marker", value = "Effect_i") %>%
          mutate(Effect_j = get_pm(m)[union(sig_i, sig_j), j],
                 ratio = .data$Effect_i/.data$Effect_j,  ##divide effect sizes, if this ratio is positive there is not AP
                 APratio = .data$Effect_i/-.data$Effect_j)  ##divide effect sizes, if this ratio is positive there is AP

        ## GxE: we are sure of the sign for two effects, and they are the same sign
        # No GxE in this pair - effects are same sign and same mag
        ms_sig2_noGxE <- effects_df %>%
          filter(between(.data$ratio, factor, 1/factor))

        # DS: we are sure of the sign for two effects, and they are the same sign
        ms_sig2_DS <- effects_df %>%
          filter(.data$ratio > 0 & !between(.data$ratio, factor, 1/factor))

        ## GxE we are sure of the sign for two effects, and they are opposite
        # AP: we are sure of the sign for two effects, and they are opposite
        ms_sig2_AP <- effects_df %>%
          filter(between(.data$APratio, 0, 1E10))

        S_all[i, j] = sum(length(ms_isigi), length(ms_jsigj), nrow(ms_sig2_noGxE),
                          nrow(ms_sig2_DS), nrow(ms_sig2_AP))
        S_CN[i, j] = sum(length(ms_isigi), length(ms_jsigj), nrow(ms_sig2_DS))
        # Not conservative!!
        S_2_no[i, j] = sum(nrow(ms_sig2_noGxE))
        S_AP[i, j] = sum(nrow(ms_sig2_AP))
        S_DS[i, j] = sum(nrow(ms_sig2_DS))
        S_i[i, j] = length(ms_isigi)
        S_j[i, j] = length(ms_jsigj)
      }
    }
  }
  return(list(S_all_pairwise = S_all, S_CN = S_CN, S_2_no = S_2_no, S_AP = S_AP,
              S_DS = S_DS, NS_pairwise = NS_pair, S_1_row = S_i,
              S_1_col = S_j))
}


# --- Plot & Save Plots ---------

#' ggplot of single mash effect
#'
#' @description Creates a plot with point estimates and standard errors for
#'     effects of a single SNP in multiple conditions.
#'
#' @param m An object of type mash
#' @param n Optional. Integer or integer vector. The result number to plot, in
#'     order of significance. 1 would be the top result, for example. Find
#'     these with \code{\link{get_significant_results}}.
#' @param i Optional. Integer or integer vector. The result number to plot, in
#'     the order of the mash object. 1 would be the first marker in the mash
#'     object, for example. Find these with \code{\link{get_marker_df}}.
#' @param marker Optional. Print the marker name on the plot?
#' @param saveoutput Logical. Should the output be saved to the path?
#'
#' @note Specify only one of n or i.
#'
#' @importFrom ashr get_psd
#' @importFrom cowplot save_plot
#' @importFrom tibble enframe
#' @importFrom dplyr mutate case_when
#' @importFrom tidyr separate
#' @import ggplot2
#' @importFrom purrr as_vector
#' @importFrom stringr str_replace str_replace_all
#'
#' @export
mash_plot_effects <- function(m, n = NA, i = NA, marker = TRUE,
                              saveoutput = FALSE){
  stopifnot((typeof(n) %in% c("double", "integer") | typeof(i) %in% c("double", "integer")))
  if(typeof(n) != "logical"){
    i <- get_significant_results(m)[n]
    marker_df <- names(i) %>%
      enframe(name = NULL, value = "Marker") %>%
      separate(.data$Marker, into = c("Chr", "Pos"), sep = "_", convert = TRUE) %>%
      mutate(Mb = round(.data$Pos / 1000000, digits = 1)) %>%
      mutate(Marker = case_when(typeof(.data$Chr) == "integer" &
                                  .data$Chr < 10 ~
                                  paste0("Chr0", .data$Chr, " ",
                                         .data$Mb, " Mb"),
                                typeof(.data$Chr) == "integer" &
                                  .data$Chr >= 10 ~
                                  paste0("Chr", .data$Chr, " ",
                                         .data$Mb, " Mb"),
                                TRUE ~ paste0(Chr, " ", .data$Mb, " Mb")
                                ))
      marker_name <- marker_df$Marker
  } else {
    marker_df <- get_marker_df(m)[i,] %>%
      separate(.data$Marker, into = c("Chr", "Pos"), sep = "_",
               convert = TRUE) %>%
      mutate(Mb = round(.data$Pos / 1000000, digits = 1)) %>%
      mutate(Marker = case_when(typeof(.data$Chr) == "integer"
                                & .data$Chr < 10 ~
                                  paste0("Chr0", .data$Chr, " ", Mb, " Mb"),
                                typeof(.data$Chr) == "integer" &
                                  .data$Chr >= 10 ~
                                  paste0("Chr", .data$Chr, " ", .data$Mb,
                                         " Mb"),
                                TRUE ~ paste0(.data$Chr, " ", .data$Mb, " Mb")
                                ))
    marker_name <- marker_df$Marker
  }
  effectplot <- get_colnames(m) %>%
    enframe(name = NULL, value = "Conditions") %>%
    mutate(mn = get_pm(m)[i,],
           se = get_psd(m)[i,]) %>%
    mutate(Conditions = str_replace(.data$Conditions,
                               "^Stand_Bhat_?",
                               ""),
           Conditions = str_replace(.data$Conditions,
                               "^Bhat_?",
                               ""),
           )

  ggobject <- ggplot(data = effectplot) +
    geom_point(mapping = aes(x = as.factor(.data$Conditions), y = .data$mn)) +
    CDBNgenomics::theme_oeco +
    geom_errorbar(mapping = aes(ymin = .data$mn - .data$se,
                                ymax = .data$mn + .data$se,
                                x = .data$Conditions), width = 0.3) +
    geom_hline(yintercept = 0, lty = 2) +
    labs(x = "Conditions", y = "Effect Size") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  if(marker == TRUE){
    ggobject <- ggobject + ggtitle(label = marker_name)
  }
  if(saveoutput == TRUE){
      save_plot(filename = paste0("Effect_plot_", str_replace_all(marker_name,
                                                              " ", "_"), "_",
                                  get_date_filename(), ".png"),
                plot = ggobject, base_aspect_ratio = 0.9, base_height = 3.5)
  }
  return(list(marker = marker_name, effect_df = effectplot,
              ggobject = ggobject))
}

#' ggplot of covariance matrix masses
#'
#' @description Creates a bar plot using ggplot of the masses that are on each
#'     covariance matrix specified in the mash model.
#'
#' @param m An object of type mash
#' @param saveoutput Logical. Should the output be saved to the path?
#'
#' @importFrom cowplot save_plot
#' @importFrom tibble enframe
#' @importFrom dplyr mutate arrange desc
#' @importFrom stringr str_replace
#' @import ggplot2
#'
#' @note This plot can be useful for seeing the overall patterns of effects in
#'     the data used in mash. Non-significant effects will add mass to the
#'     "no_effects" covariance matrix, while significant effects will add mass
#'     to one of the other covariance matrices. You can use GGally::ggcorr() to
#'     plot the covariance matrix patterns themselves.
#'
#' @export
mash_plot_covar <- function(m, saveoutput = FALSE){
  df <- get_estimated_pi(m)
  df <- enframe(df, name = "Covariance Matrix", value = "Mass") %>%
    mutate(`Covariance Matrix` = str_replace(.data$`Covariance Matrix`,
                                             "^Bhat_?", "single_effect_"),
           `Covariance Matrix` = str_replace(.data$`Covariance Matrix`,
                                             "^null$", "no_effects"),
           `Covariance Matrix` = str_replace(.data$`Covariance Matrix`,
                                             "^Stand_Bhat_?",
                                             "single_effect_")) %>%
    arrange(desc(.data$`Covariance Matrix`))
  df$`Covariance Matrix` <- factor(df$`Covariance Matrix`,
                                   levels = (df$`Covariance Matrix`))
  ggobject <- ggplot(df) +
    geom_bar(aes(x = .data$Mass, y = .data$`Covariance Matrix`),
             stat = "identity") +
    CDBNgenomics::theme_oeco

  if(saveoutput == TRUE){
    save_plot(paste0("Covariance_matrix_mass_plot_", get_date_filename(),
                     ".png"), plot = ggobject, base_aspect_ratio = 1,
              base_height = 4)
  }
  return(list(covar_df = df, ggobject = ggobject))
}


#' @title Manhattan plot in ggplot colored by significant conditions
#'
#' @description Takes a mash object and, for some vector of phenotypes, returns
#'     a Manhattan plot ggplot object (and its dataframe). Each SNP in the plot
#'     is colored by the number of phenotypes it is significant for. Even and
#'     odd chromosomes have different shapes for their SNPs, so that
#'     chromosome identity can be determined.
#'
#' @param m A mash object (outputted by mash).
#' @param cond A vector of phenotypes. Defaults to the names of each
#'     column in the mash object.
#' @param saveoutput Logical. Should the output be saved to the path?
#' @param thresh Numeric. The threshold used for the local false sign rate to
#'     call significance in a condition.
#'
#' @return A \code{tbl_df()} of the data used to make the Manhattan plot, and a
#'     ggplot object containing the Manhattan.
#'
#' @importFrom cowplot save_plot
#' @importFrom dplyr rename select arrange mutate left_join
#' @import ggplot2
#' @importFrom tibble as_tibble rownames_to_column enframe
#' @importFrom tidyr separate
#' @import viridis
#'
#' @examples
#' \dontrun{manhattan_out <- mash_ggman_by_condition(m = m, saveoutput = TRUE)}
#'
#' @export
mash_plot_manhattan_by_condition <- function(m, cond = NA, saveoutput = FALSE,
                                             thresh = 0.05){
  num_sig_in_cond <- c()

  if(is.na(cond)[1]){
    cond <- get_colnames(m = m)
  }

  log10bf_df <- get_log10bf(m = m) %>%
    as.data.frame() %>%
    rownames_to_column(var = "value") %>%
    mutate(value = as.integer(.data$value)) %>%
    as_tibble() %>%
    left_join(get_marker_df(m = m)) %>%
    dplyr::rename(log10BayesFactor = .data$V1) %>%
    dplyr::select(-.data$value)

  ggman_df <- get_n_significant_conditions(m = m, thresh = thresh,
                                           conditions = cond) %>%
    enframe(name = "Marker", value = "Num_Sig_Conditions") %>%
    separate(.data$Marker, into = c("Chr", "Pos"), remove = FALSE, sep = "_",
             extra = "merge", convert = TRUE) %>%
    left_join(log10bf_df, by = "Marker") %>%
    arrange(.data$Chr, .data$Pos)

  log10BF <- expression(paste("log"[10], plain("(Bayes Factor)")))

  ggmanobject <- ggplot(data = ggman_df, aes(x = .data$Pos, y = .data$log10BayesFactor)) +
    CDBNgenomics::theme_oeco +
    geom_point(aes(color = .data$Num_Sig_Conditions, fill = .data$Num_Sig_Conditions,
                   shape = as.factor(.data$Chr))) +
    facet_wrap(~ .data$Chr, nrow = 1, scales = "free_x", strip.position = "bottom") +
    scale_color_viridis(option = "B") + scale_fill_viridis(option = "B") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.background = element_rect(fill=NA)) +
    labs(x = "Chromosome", y = log10BF) +
    scale_x_continuous(expand = c(0.25, 0.25)) +
    scale_shape_manual(values = rep(c(21,22),9), guide = FALSE)

  if(saveoutput == TRUE){
    save_plot(paste0("Manhattan_mash_", get_date_filename(),
                     ".png"), plot = ggmanobject, base_aspect_ratio = 2.4,
              base_height = 3)
  }

  return(list(ggman_df = ggman_df, ggmanobject = ggmanobject))
}


#' @title Create a ggplot of pairwise sharing of mash effects
#'
#' @description Given a correlation matrix, an RDS with a correlation matrix, or
#'     a mash object, create a ggplot of pairwise sharing of mash effects using
#'     \code{\link{get_pairwise_sharing}} and \code{\link{ggcorr}}.
#'
#' @param m An object of type mash
#' @param effectRDS An RDS containing a correlation matrix.
#' @param corrmatrix A correlation matrix
#' @param reorder Logical. Should the columns be reordered by similarity?
#' @param saveoutput Logical. Should the output be saved to the path?
#' @param filename Character string with an output filename. Optional.
#' @param ... Other arguments to \code{\link{get_pairwise_sharing}} or
#'      \code{\link{ggcorr}}.
#'
#' @importFrom GGally ggcorr
#' @import viridis
#' @importFrom stringr str_replace_all
#'
#' @return A list containing a dataframe containing the correlations and a
#'     ggplot2 object containing the correlation plot.
#'
#' @export
mash_plot_pairwise_sharing <- function(m = NULL, effectRDS = NULL,
                                       corrmatrix = NULL, reorder = TRUE,
                                       saveoutput = FALSE, filename = NA, ...){
  # Additional arguments for get_pairwise_sharing, ggcorr, and save_plot
  requireNamespace("dots")
  factor <- dots::dots(name = 'factor', value = 0.5, ...)
  lfsr_thresh <- dots::dots(name = 'lfsr_thresh', value = 0.05, ...)
  FUN <- dots::dots(name = 'FUN', value = identity, ...)
  geom <- dots::dots(name = 'geom', value = 'circle', ...)
  label <- dots::dots(name = 'label', value = FALSE, ...)
  label_alpha <- dots::dots(name = 'label_alpha', value = TRUE, ...)
  label_size <- dots::dots(name = 'label_size', value = 3, ...)
  hjust <- dots::dots(name = 'hjust', value = 0.95, ...)
  vjust <- dots::dots(name = 'vjust', value = 0.3, ...)
  layout.exp <- dots::dots(name = 'layout.exp', value = 9, ...)
  min_size <- dots::dots(name = 'min_size', value = 0, ...)
  max_size <- dots::dots(name = 'max_size', value = 3.5, ...)
  option <- dots::dots(name = 'option', value = 'B', ...)
  dpi <- dots::dots(name = 'dpi', value = 500, ...)

  base_aspect_ratio <- dots::dots(name = 'base_aspect_ratio', value = 1.1, ...)

  if(is.na(filename)[1]){
    filename <- paste0("Mash_pairwise_shared_effects_",
                       get_date_filename(), ".png")
  }

  # look for a shared effects matrix in the path, and if not, generate one
  if(!is.null(effectRDS) && is.null(m) && is.null(corrmatrix)){
    shared_effects <- readRDS(effectRDS)
  } else if(!is.null(corrmatrix) && is.null(effectRDS) && is.null(m)){
    shared_effects <- corrmatrix
  } else if(!is.null(m)){
    shared_effects <- get_pairwise_sharing(m = m, factor = factor,
                                           lfsr_thresh = lfsr_thresh, FUN = FUN)
    rownames(shared_effects) <- str_replace_all(rownames(shared_effects),
                                                "^Stand_Bhat_?", "")
    rownames(shared_effects) <- str_replace_all(rownames(shared_effects),
                                                "^Bhat_?", "")
    colnames(shared_effects) <- str_replace_all(colnames(shared_effects),
                                               "^Stand_Bhat_?", "")
    colnames(shared_effects) <- str_replace_all(colnames(shared_effects),
                    "^Bhat_?", "")
  } else {
    stop(paste0("Please specify one of these: ",
                "1. a mash output object (m), ",
                "2. the path to a effect rds file (mashRDS), ",
                "3.  a correlation matrix (corrmatrix)."))
  }

  base_height <- dots::dots(name = 'base_height',
                            value = nrow(shared_effects)*0.33+1, ...)

  if(reorder == TRUE){
    corrdf <- reorder_cormat(cormat = shared_effects)
    corrplot <- ggcorr(data = NULL, cor_matrix = corrdf, geom = geom,
                       label = label, label_alpha = label_alpha,
                       label_size = label_size, hjust = hjust, vjust = vjust,
                       layout.exp = layout.exp, min_size = min_size,
                       max_size = max_size) +
      scale_color_viridis(option = option)
  } else {
    corrplot <- ggcorr(data = NULL, cor_matrix = shared_effects, geom = geom,
                       label = label, label_alpha = label_alpha,
                       label_size = label_size, hjust = hjust, vjust = vjust,
                       layout.exp = layout.exp, min_size = min_size,
                       max_size = max_size) +
      scale_color_viridis(option = option)
  }

  if(saveoutput == TRUE){
    save_plot(filename = filename, corrplot,
              base_aspect_ratio = base_aspect_ratio, base_height = base_height,
              dpi = dpi)
  }
  return(list(corr_matrix = shared_effects, gg_corr = corrplot))
}

#' @title Significant SNPs per number of conditions
#'
#' @description For some number of columns in a mash object that correspond to
#'     conditions, find the number of SNPs that are significant for that number
#'     of conditions.
#'
#' @param m An object of type mash
#' @param conditions A vector of conditions. Get these with get_colnames(m).
#' @param saveoutput Logical. Save plot output to a file? Default is FALSE.
#' @param thresh What is the threshold to call an effect significant? Default is
#'     0.05.
#'
#' @return A list containing a dataframe of the number of SNPs significant per
#'     number of conditions, and a ggplot object using that dataframe.
#'
#' @import ggplot2
#' @importFrom tibble enframe
#' @importFrom dplyr rename summarise filter group_by n
#'
#' @examples
#'   \dontrun{mash_plot_sig_by_condition(m = mash_obj, saveoutput = TRUE)}
#'
#' @export
mash_plot_sig_by_condition <- function(m, conditions = NA, saveoutput = FALSE,
                                       thresh = 0.05){

  thresh <- as.numeric(thresh)
  num_sig_in_cond <- c()

  if(typeof(conditions) == "logical"){
    cond <- get_colnames(m)
  } else {
    cond <- conditions
  }

  SigHist <- get_n_significant_conditions(m = m, thresh = thresh,
                                          conditions = cond) %>%
    enframe(name = "Marker") %>%
    rename(Number_of_Conditions = .data$value) %>%
    group_by(.data$Number_of_Conditions) %>%
    summarise(Significant_SNPs = n()) %>%
    filter(.data$Number_of_Conditions != 0)

  vis <- ggplot(SigHist, aes(x = .data$Number_of_Conditions, y = .data$Significant_SNPs)) +
    CDBNgenomics::theme_oeco +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = 0, lty = 2) +
    xlab(label = "Number of Conditions") +
    ylab(label = "Number of Significant SNPs")

  if(saveoutput == TRUE){
    ggsave(paste0("SNPs_significant_by_number_of_conditions_",
                  get_date_filename(),
                  ".bmp"), width = 5, height = 3, units = "in", dpi = 400)
  }

  return(list(sighist = SigHist, ggobject = vis))
}


#' @title Reorder correlation matrix
#'
#' @description  Reorder correlation coefficients from a matrix of things
#'     (including NA's) and hierarchically cluster them
#'
#' @param cormat A correlation matrix
#'
#' @importFrom cluster daisy
#' @importFrom stats hclust
#'
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- daisy(cormat, metric = "gower")
  hc <- hclust(dd)
  cormat <- cormat[hc$order, hc$order]
}
