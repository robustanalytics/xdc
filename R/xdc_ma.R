#-----------------------------------------------
#' Mixture-based Moderator Estimation
#'
#' Using the latent mixture based method to
#' estimate a moderating effect.
#'
#' @param es, the reported effect sizes
#' @param se, the reported standard errors
#' @param ml, the independent variables
#' @param k, the number of mixture components,
#' default to 2
#' @param method, the method for moderator test,
#'     default to meta-regression, possible
#'     values include:
#'  - RMA: meta-regression,
#'  - GRP: subgroup comparison, using
#'       abs(mean1-mean2)/sqrt(se1^2+se2^2)
#' @return output of moderator analysis
#-----------------------------------------------
# internal parameter:
#   - speed, options:1-5, default to 5
#   - resolution, whether additional resolution
#     is needed (e.g., when most reported es are
#     within a tight interval, e.g., 0-0.5)
#' @export
xdc.ma = function(es, se, ml, k = 2, method="RMA"){

  speed = 5; eps = 0.1; resolution=TRUE;
  # equi_sd is used by mixture decomposition to optimize performance when mdc is provably reduces to regression
  resolution = (resolution & k > 2 & sd(es)/length(es) < eps/10);
  if (resolution | k > 2) {
    equi_sd = FALSE;
  } else {
    vt_model = var.test(es[which(ml==0)], es[which(ml==1)], alternative = "two.sided");
    lm_model = lm(es ~ ml);
    equi_sd = (vt_model$p.value > 0.05 & (speed >= 3 | summary(lm_model)$coefficients[2,4] < 0.05));
  }

  # save the observed histogram as test_histogram
  test_histogram = get_histogram_from_data(es, se);
  output_mx = get_combis(es, test_histogram,
                         k = k,                                               # pass through number of components
                         speed=speed,                                         # pass through speed
                         topx=ifelse(resolution, 10,                          # if resolution is set, leave fine tuning to later
                                     ifelse(speed >= 3, 100,                  # otherwise, test top 100
                                            1000)),                           # test top 1000 (very slow)
                         equisd = equi_sd,                                    # speed optimization when mdc is provably reduces to regression
                         resolution = resolution);  # need fine tuning when values are dense

  # if visualization option is on, display two curves:
  #   - the observed PDF (test_histogram)
  #   - the mixture composition (output_mx)
  # the dotted line represents the observed PDF

  # return output of get_lm
  return(list(
    mixture = output_mx,
    model = get_lm(es, se, ml, output_mx, method = method)
  ));
}
