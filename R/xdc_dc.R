#-----------------------------------------------
#' Mixture Decomposition
#'
#' Using the histogram structure to approximately
#' identify the mixture composition that
#' minimizes the distance with the input
#' distribution.
#'
#' @param es, the reported effect sizes
#' @param se, the reported standard errors
#' @param k, the number of mixture components,
#' default to 2
#' @return output of mixture decomposition, a
#' data frame with three columns: (mean, sd,
#' weight), and k rows representing the k
#' components
#-----------------------------------------------
#' @import dplyr
#' @export
xdc.dc = function(es, se, k = 2){

  speed = 5; eps = 0.1; resolution=TRUE; equi_sd = FALSE;
  resolution = (resolution & k > 2 & sd(es)/length(es) < eps/10);

  return(get_combis(es, get_histogram_from_data(es, se),
                         k = k,                                               # pass through number of components
                         speed=speed,                                         # pass through speed
                         topx=ifelse(resolution, 10,                          # if resolution is set, leave fine tuning to later
                                     ifelse(speed >= 3, 100,                  # otherwise, test top 100
                                            1000)),                           # test top 1000 (very slow)
                         equisd = equi_sd,                                    # speed optimization when mdc is provably reduces to regression
                         resolution = resolution));  # need fine tuning when values are dense
}
