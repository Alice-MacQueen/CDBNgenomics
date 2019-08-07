# --- Get Results from Mash Object ---------

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
#' @param saveoutput Logical. Should the output be saved to the path?
#'
#' @note Specify only one of n or i.
#'
#' @importFrom ashr get_psd
#' @importFrom cowplot save_plot
#' @importFrom tibble enframe
#' @importFrom dplyr mutate
#' @import ggplot2
#' @importFrom purrr as_vector
#'
#' @export
mash_plot_effects <- function(m, n = NA, i = NA, saveoutput = FALSE){
  stopifnot((!is.na(n[1]) | !is.na(i[1])))
  if(is.na(i[1])){
    i <- get_significant_results(m)[n]
  }

  effectplot <- get_colnames(m) %>%
    enframe(name = "Conditions") %>%
    mutate(mn = get_pm(m)[i,],
           se = get_psd(m)[i,])

  ggobject <- ggplot(data = effectplot) +
    geom_point(mapping = aes(x = as.factor(.data$value), y = .data$mn)) +
    geom_errorbar(mapping = aes(ymin = .data$mn - .data$se,
                                ymax = .data$mn + .data$se,
                                x = .data$Conditions), width = 0.3) +
    geom_hline(yintercept = 0, lty = 2) +
    labs(x = "Conditions", y = "Effect Size") +
    scale_x_discrete(labels = as_vector(.data$value)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  if(saveoutput == TRUE){
    if(is.na(n[1])){
      save_plot(filename = paste0("Effect_plot_",
                                  names(get_significant_results(m))[n], ".png"),
                plot = ggobject, base_aspect_ratio = 0.8, base_height = 4.5)
    } else {
      plotname <- get_marker_df(m)[i]
      save_plot(filename = paste0("Effect_plot_", plotname$Marker, ".png"),
                plot = ggobject, base_aspect_ratio = 0.8, base_height = 4.5)

    }
  }
  return(list(marker = i, effect_df = effectplot, ggobject = ggobject))
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
#' @importFrom stringr str_replace_all
#'
#' @examples
#' \dontrun{manhattan_out <- mash_ggman_by_condition(m = m, saveoutput = TRUE)}
#'
#' @export
mash_plot_manhattan_by_condition <- function(m, cond = NA,
                                             saveoutput = FALSE, thresh = 0.05){
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
    enframe(name = "Marker") %>%
    rename(Num_Sig_Conditions = .data$value) %>%
    separate(.data$Marker, into = c("Chr", "Pos"), remove = FALSE, sep = "_",
             extra = "merge") %>%
    mutate(Pos = as.numeric(.data$Pos)) %>%
    left_join(log10bf_df, by = "Marker") %>%
    arrange(.data$Chr, .data$Pos)

  log10BF <- expression(paste("log"[10], plain("(Bayes Factor)")))

  ggmanobject <- ggplot(data = ggman_df, aes(x = .data$Pos, y = .data$log10BayesFactor)) +
    geom_point(aes(color = .data$Num_Sig_Conditions, fill = .data$Num_Sig_Conditions,
                   shape = as.factor(.data$Chr))) +
    facet_wrap(~ .data$Chr, nrow = 1, scales = "free_x", strip.position = "bottom") +
    scale_color_viridis(option = "B") + scale_fill_viridis(option = "B") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.background = element_rect(fill=NA)) +
    labs(x = "Chromosome", y = log10BF) +
    scale_x_continuous(expand = c(0.3, 0.3)) +
    scale_shape_manual(values = rep(c(21,22),9), guide = FALSE)

  if(saveoutput == TRUE){
    save_plot(paste0("Manhattan_mash_", str_replace_all(Sys.time(), ":", "."),
                     ".png"), plot = ggmanobject, base_aspect_ratio = 2.5,
              base_height = 4)
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
                       str_replace_all(Sys.time(), ":", "."), ".png")
  }

  # look for a shared effects matrix in the path, and if not, generate one
  if(!is.null(effectRDS) && is.null(m) && is.null(corrmatrix)){
    shared_effects <- readRDS(effectRDS)
  } else if(!is.null(corrmatrix) && is.null(effectRDS) && is.null(m)){
    shared_effects <- corrmatrix
  } else if(!is.null(m)){
    shared_effects <- get_pairwise_sharing(m = m, factor = factor,
                                           lfsr_thresh = lfsr_thresh, FUN = FUN)
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
#' @param conditions A vector of conditions
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

  if(is.na(conditions)[1]){
    cond <- get_colnames(m = m)
  }

  SigHist <- get_n_significant_conditions(m = m, thresh = thresh,
                                          conditions = cond) %>%
    enframe(name = "Marker") %>%
    rename(Number_of_Conditions = .data$value) %>%
    group_by(.data$Number_of_Conditions) %>%
    summarise(Significant_SNPs = n()) %>%
    filter(.data$Number_of_Conditions != 0)

  vis <- ggplot(SigHist, aes(x = .data$Number_of_Conditions, y = .data$Significant_SNPs)) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = 0, lty = 2) +
    xlab(label = "Number of Conditions") +
    ylab(label = "Number of Significant SNPs")

  if(saveoutput == TRUE){
    ggsave(paste0("SNPs with significant effects in n conditions ",
                  str_replace_all(Sys.time(), ":", "."),
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
