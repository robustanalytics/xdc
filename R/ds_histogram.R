#-----------------------------------------------
#' Getting the Histogram from Observed Effect
#' Sizes and Standard Errors
#'
#' The returned probability vector may not sum
#' to 1. The remaining value represents the
#' probability falling outside the range
#' specified in the boundaries.
#' @param es, a vector of effect sizes
#' @param se, a vector of standard errors (es
#' and se must be equi-length)
#' @return a data frame with two columns:
#' - x: same as input boundaries
#' - p: probability for each bin
#-----------------------------------------------
#' @export
xdc.hist.data = function(es, se) {
  return(get_histogram_from_data(es, se));
}

#-----------------------------------------------
#' Getting the Histogram from a Mixture
#'
#' This function supports an arbitrary number of
#' components. The returned probability vector
#' may not sum to 1. The remaining value
#' represents the probability falling outside
#' the range specified in the boundaries
#' @param boundaries, an array of (upper) bin
#' boundaries. This is usually obtained from
#' the histogram of observed data (xdc.hist.data)
#' @param mxs, a data frame with three columns:
#' (mean, sd, weight), representing the mixture
#' composition
#' @return a data frame with two columns:
#' - x: same as input boundaries
#' - p: probability for each bin
#-----------------------------------------------
#' @export
xdc.hist.mx = function(boundaries, mxs){
  return(get_histogram_from_mixture(boundaries, mxs));
}

#-----------------------------------------------
# get the histogram from a mixture
# input parameters:
#   - boundaries, an array of (upper) bin
#     boundaries
#   - mxs, a data frame with three columns:
#     (mean, sd, weight), hereinafter 3-col
# output:
#   - a data frame with two variables:
#       - x: same as input boundaries
#       - p: probability for each bin
# note:
#   - sum of p should be smaller than or equal
#     to 1. (1-p) represents the probability
#     falling outside the range specified in
#     boundaries (see below for details)
#   - supports an arbitrary number of components
#-----------------------------------------------

get_histogram_from_mixture = function(boundaries, mxs){
  colnames(mxs) = c("mean", "sd", "weight");
  pdf_est = future.apply::future_sapply(seq(1,length(boundaries)), function(x){
    bin_ub = boundaries[x];

    # for the first bin, we set its lower bound following the same width as the second bin
    # theoretically, we could have set bin_lb to -\infty, but doing so would cause a spike
    # on histogram because of the computational optimization we took on shrinking the spread
    # of all bins to as small a range as possible

    bin_lb = ifelse(x==1, boundaries[1]*2-boundaries[2], boundaries[x-1]);
    return(sum((pnorm(bin_ub,mxs$mean,mxs$sd)-pnorm(bin_lb,mxs$mean,mxs$sd))*mxs$weight));
  });
  output = data.frame(boundaries, pdf_est);
  colnames(output) = c("x","p");
  return(output);
}

#-----------------------------------------------
# get the input data structure G_in
# assumed global parameters:
#   - min_w, the minimum
#     weight for a component
#   - eps, the maximum tolerable error
# input parameters:
#   - es, an array of effect sizes
#   - si, an array of standard errors
# output:
#   - the histogram data structure for the
#     observed data (two vars: x, p)
# requirements:
#   - es and si must be equi-length
#-----------------------------------------------

get_histogram_from_data = function(es, si){
  min_w = 0.1; eps = 0.1; input_precision = 0.01;
  # N is the input length;
  N = length(es);
  window_size = round(min_w*N);

  # computing bin_width with sanity check, to avoid unnecessarily large number of bins
  sortes = sort(es);
  bin_width = 3.49*sqrt(pi)*N^(-1/3)*min(future.apply::future_sapply(window_size:N,
                                                       function(x){as.numeric(Lmoments::Lmoments(sortes[(x-window_size+1):x],2)[1,2])}
  ));
  bin_width = max(bin_width, input_precision);

  # limit range of histogram based on min_w, to avoid unnecessarily large number of bins
  h_ub = as.numeric(quantile(es, probs = 1-min_w/4));
  h_lb = as.numeric(quantile(es, probs = min_w/4));
  boundaries = seq(h_lb, h_ub, bin_width);
  mxs = data.frame(es,pmax(eps/3,si),rep(1/N, N));
  return(get_histogram_from_mixture(boundaries, mxs));
}

#-----------------------------------------------
# get the distance between a given data
#   histogram and a given mixture composition
# input parameters:
#   - data_histogram
#   - mxs (3-col format)
#   - style, default to \ell_1-norm, supporting
#       - "L1": \ell_1 norm
#       - "LI": \ell_\infty norm
#       - "SLI": single sided \ell_\infty norm
#         designed to handle cases where there
#         are other unknown mixture components,
#         so only one-side difference matters
# note:
#   - this function is here to be extended
#     later to accomodate more vector norm
#     designs (e.g., L2)
#   - supports an arbitrary number of components
#-----------------------------------------------

get_mixture_distance = function(data_histogram, mxs, style="L1"){
  if (style == "LI") {
    return(max(abs(data_histogram$p - (get_histogram_from_mixture(data_histogram$x, mxs))$p)));
  } else if (style == "SLI") {
    return(max(data_histogram$p - (get_histogram_from_mixture(data_histogram$x, mxs))$p));
  } else if (style == "L1") {
    return(sum(abs(data_histogram$p - get_histogram_from_mixture(data_histogram$x, as.data.frame(mxs))$p)));
  }
}
