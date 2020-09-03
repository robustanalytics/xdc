#-----------------------------------------------
#' Plot a Mixture
#'
#' Plotting a Gaussian mixture based on its
#' parameter specifications (supporting an
#' arbitrary number of mixture components)
#' @param mx, the mixture specification in 3-col
#' format, i.e., a data frame with three columns
#' (mean, sd, weight)
#' @param hist, optional, may be
#' specified as the output of the xdc.hist.data
#' function. If specified, the histogram of the
#' data will be plotted as a dotted line. The x
#' range of the histogram will also be used for
#' the entire plot
#' @param df, optional, will only be used if
#' mod is also specified, so as to plot multiple
#' lines for different subsets of data. df
#' should be a data frame with at least two
#' columns: es and si. It could also contain a
#' column marking the moderator.
#' @param mod, optional, could be a vector
#' containing the moderator variable.
#' Alternatively, it could specify a column or
#' any expression over the columns in df
#' @param colors, optional, specifying the
#' colors used to plot the lines corresponding
#' to the moderator levels. default to
#' a few high-contrast colors followed by
#' rainbow()
#' @return displays a line chart with the
#' PDF of the mixture distribution and,
#' optionally, the histogram of the data
#-----------------------------------------------
#' @export
xdc.plot.mx = function(mx, hist = c(), df = c(), mod = "", colors=c()) {

  # if the data histogram is specified, then we derive the x range from the histogram
  # and plot the histogram. Otherwise, we compute the range of x from the mixture in mx

  if (length(hist) > 0) {
    plot_x = hist$x;
  } else {
    lb = min(mx[,1]-3*mx[,2]);
    ub = max(mx[,1]+3*mx[,2]);
    plot_x = seq(lb, ub, (ub-lb)/100);
  }

  df_plot = data.frame(es = plot_x);

  mixture_histogram = get_histogram_from_mixture(plot_x,mx);

  df_plot$mixture = mixture_histogram$p;
  line_marker = c("dashed");

  if (length(hist) > 0) {
    df_plot$data = hist$p;
    line_marker = c(line_marker, "dotted");

    gmean = sum(hist$x*hist$p);
    gsd = sqrt(sum(hist$p*(hist$x - gmean)^2));
    gmx = data.frame(t(c(gmean,gsd,1)));
    df_plot$Gaussian = get_histogram_from_mixture(plot_x,gmx)$p;
  }

  # If moderator is specified, we plot the mixture for each moderator level
  if (length(mod) > 1 | mod[1] != "") {
    # get mod if it is specified in terms of an expression of the columns in df
    if (is.character(mod) && length(mod) == 1) mod = eval(substitute(mod),df,parent.frame());

    mod_levels = sort(unique(mod));

    # get the default color palette
    if (length(colors) == 0) colors = c("#4682B4","#AF46B4","#B47846","#4BB446",rainbow(length(mod_levels)));

    # first, compute an n \times k matrix jp_matrix for the n records and k mixture components, marking the joint probabilities
    joint_prob = function(x)(mx[x,3]*dnorm(es,mx[x,1],mx[x,2]));
    jp_matrix = sapply(seq(1,nrow(mx)), function(x){ joint_prob(x)/rowSums(sapply(1:nrow(mx), joint_prob)) });

    for (modv in mod_levels) {
      viz_mx = mx;
      # adjusting the component weights according to the specific moderator level
      viz_mx[,3] = sapply(seq(1,nrow(mx)), function(x) mean(jp_matrix[which(mod == modv),x]));
      score_mixture = get_histogram_from_mixture(plot_x, viz_mx);

      df_plot[[paste0("mod_", modv)]] = score_mixture$p;
      line_marker = c(line_marker, "solid");
    }
  }

  df_plot_long = reshape2::melt(df_plot, id="es");  # convert to long format
  df_plot_long = df_plot_long %>% rename(lines = variable);

  ggplot2::ggplot(data=df_plot_long,
                  ggplot2::aes(x=es, y=value)) +
    ggplot2::geom_line(ggplot2::aes(color=lines)) +
    ggplot2::ylab("Probability") +
    ggplot2::labs(title = "Mixture Histogram");
}
