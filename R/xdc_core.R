#-----------------------------------------------
# get the first cut on the component candidates
# input parameters:
#   - data_histogram, the target histogram to
#     be approximated by the mixture candidates
#   - mean_cand, candidates for mean
#   - sd_cand, candidates for sd
#   - weight_cand, candidates for weight
#   - equisd, whether mdc is reduced to
#     regression, so candidates can be pruned
#     aggressively
#   - speed, the speed preference, default to 5
# output:
#   - a data frame with five variables: mean,
#     sd, weight, mean2, sd2
#-----------------------------------------------

get_first_combis = function(data_histogram, mean_cand, sd_cand, weight_cand, equisd = FALSE, speed = 5, sdes = 100) {

  # this function implements two methods:
  #   - method 1 derives the second component by deducting the first component from the data_histogram
  #     it uses robust statistics to directly infer the parameters for the second component
  #   - method 2 is the brute-force method: simply have the Cartesian product of candidates for each component
  #     as the candidates for the mixture composition
  # In general, Method 1 is more efficient when the number of candidate values in mean_cand and/or sd_cand are
  # very large, because the number of candidates with Method 2 grows quadratically with the size of mean_cand
  # and sd_cand. However, when mean_cand and sd_cand are relatively small, Method 2 could be more efficient,
  # because it saves the computational time of having to deduct each component candidate from the data histogram.
  # Given that R often runs out of memory when we set eps (the tolerable error) below 0.1, it appears that Method
  # 2 is usually the more efficient method. This is particularly pronounced when equisd is TRUE, as the candidates
  # can be aggressively pruned in this degenerated scenario.

  method_id = ifelse(speed >= 3 | equisd, ifelse(length(mean_cand) * length(sd_cand) > 10000 & !equisd, 1, 2), 1);

  if (method_id == 1) {
    rcombis = as.data.frame(expand.grid(mean_cand,sd_cand, weight_cand));
    rcombis2 = as.data.frame(t(apply(rcombis, 1, function(comp){
      comp_hist = get_histogram_from_mixture(data_histogram$x, as.data.frame(t(comp)));
      comp_est = (data_histogram$p - comp_hist$p)/(1-comp[3]);
      return(c(
        data_histogram$x[min(which(cumsum(comp_est)>=min(0.5,max(cumsum(comp_est)))))],
        (data_histogram$x[min(which(cumsum(comp_est)>=min(0.75,max(cumsum(comp_est)))))]-data_histogram$x[min(which(
          cumsum(comp_est)>=min(0.25,max(cumsum(comp_est)))
        ))])/(2*sqrt(2)*erfinv(1/2))
      ));
    })));

    rcombis = bind_cols(rcombis, rcombis2);
    colnames(rcombis) = c("mean","sd","weight","mean2","sd2");
  } else if (method_id == 2) {
    rcombis = as.data.frame(expand.grid(mean_cand,sd_cand,weight_cand,mean_cand,sd_cand));
    colnames(rcombis) = c("mean","sd","weight","mean2","sd2");
    if (equisd) {
      rcombis = rcombis %>%
        filter(sd == sd2) %>%
        filter((mean-mean2)^2+sd^2<sdes^2);
    }
  }

  return(rcombis);
}

#-----------------------------------------------
# select the optimal combination
# assumed global parameters:
#   - eps, the maximum tolerable error
#   - verbose, whether to display debug messages
# input parameters:
#   - cands, the input candidates, a data frame
#     with five variables (mean, sd, weight,
#     mean2, sd2)
#   - es, the observed effect sizes, used here
#     for fast pruning (see below)
#   - data_histogram, the target histogram
#     (for the output components only, excl.
#     the component emean, esd, ew if k = 3)
#   - topx, the target size of pruned candidate
#     set, default to 100
#   - emean, esd, ew: parameters for the
#     external component (used when k = 3)
#-----------------------------------------------

get_final_combis = function(cands, es, data_histogram, topx = 100, emean = 0, esd = 0, ew = 0, k = 2, input_precision = 0.01, eps = 0.1, verbose = verbose) {

  # Additional steps need to be taken if there is an external component emean, esd, ew
  # specifically, we can no longer use the raw effect sizes es for pruning, because they
  # represent both the output composition and the external component. Instead, we need to
  # first regenerate es that represent the output composition only. We do so by sampling from
  # the target histogram data_histogram. Note that the sole purpose of these sampled es is to
  # speed up the pruning process.

  if (ew > 0) {
    # pruning: no need to consider components too close to the external component
    cands = cands %>%
      filter(abs(mean - emean) > 2 * eps & abs(mean2 - emean) > 2 * eps);
    # sample es from data_histogram
    set.seed(1); raw_es = runif(length(es)*20);
    sim_cdf = c(cumsum(data_histogram$p),1);
    es_indices = sapply(raw_es, function(x){min(which(x < sim_cdf))});
    es_indices = es_indices[which(es_indices <= length(data_histogram$x))];
    es = round(data_histogram$x[es_indices]/input_precision)*input_precision;

    # In the extremely unlikely scenario where all remaining density concentrates on the max
    # value, directly return an approximate fit
    if (length(es) == 0) {
      return (c(emean, esd, 0.1, max(data_histogram$x), esd, 0.9));
    }
  }

  # What follows is a two-step process for finding the optimal combination. The reason why a
  # two-step process is needed is because of the significant computational overhead associated
  # with computing get_mixture_distance for each and every candidate mixture composition. To
  # reduce the computational overhead, we first prune from the candidates those mixture
  # compositions that are very unlikely to be the optimal choice. We do so by using a loglihood-
  # like measure lll, and only retain the topx best candidates in the second step of computing
  # get_mixture_distance and finding the best solution.

  cands = cands %>%
    filter((mean != mean2) | (sd != sd2)) %>%
    mutate(weight2 = 1 - weight);

  # In the extremely unlikely scenario where no viable candidates are left, we simply return two
  # equi-weighted components, both representing the input es.
  if (nrow(cands) == 0) {
    return(c(mean(es), sd(es), 0.5, mean(es), sd(es), 0.5));
  }

  cands$lll = rowSums(sapply(es, function(x) log(cands$weight*dnorm(x, cands$mean, cands$sd) + (1-cands$weight)*dnorm(x, cands$mean2, cands$sd2))));

  cands = cands %>%
    top_n(topx, lll) %>%
    rowwise() %>%
    mutate(score = -get_mixture_distance(
      data_histogram,
      as.data.frame((rbind(c(mean, sd, weight),c(mean2, sd2, weight2))))
    )) %>%
    mutate(fscore = ifelse(k == 2, score, floor(ifelse(weight2==1-eps,1/0.9,1)*score/eps*10))) %>%
    ungroup() %>%
    top_n(1, fscore);

  if (nrow(cands) > 3) {
    cands = cands %>% top_n(1, score);
    cands = cands[1,];
  }
  if (nrow(cands) > 1) {
    if (verbose) print(cands);
    cands = cands[order(cands[,1]),];
    if (length(unique(pmin(cands$sd, cands$sd2))) > 1) cands = cands[order(pmin(cands$sd, cands$sd2)),];
    cands = cands[1,];
  }

  return(as.numeric(cands));
}

#-----------------------------------------------
# solve the mixture decomposition problem
# assumed global parameters:
#   - min_w, the minimum weight for a component
#   - max_s, the maximum standard deviation for
#     a component
#   - eps, the maximum tolerable error
#   - beta_t, beta testing indicator
# input parameters:
#   - es. Note that es is here not for the
#     decomposition per se, but to speed up the
#     mixture decomposition process. We do not
#     use es to evaluate the final choice, but
#     use it to prune out unlikely candidates
#   - data_histogram, the observed histogram,
#     which is computed based on both es and si
#   - k, the number of mixture components
#     default to 2
#   - topx, the target size of pruned candidate
#     set, default to 100
#   - speed, default to 5
#   - resolution, default to FALSE (see function
#     mdcma below for an explanation)
#-----------------------------------------------

get_combis = function(es, data_histogram, k = 2, topx = 100, equisd = FALSE, speed = 5, resolution = FALSE){

  #-----------------------------------------------
  # parameter setup:
  # - min_w, the minimum weight for a component
  # - max_s, the maximum possible standard
  #   deivation for a mixture component
  # - input_precision, the precision of input
  #   effect sizes (e.g., 0.01)
  # - eps, the maximum tolerable error for the
  #   parameter estimation of mixture components
  # - verbose, whether to display debug messages
  # - beta_t, the beta testing level. Default to
  #   zero, i.e., no beta feature
  #-----------------------------------------------

  min_w = 0.1;
  max_s = 1;
  input_precision = 0.01;
  eps = 0.1;
  verbose = FALSE;
  beta_t = 0;

  # Constructing the candidate values for mean, sd, and weight of the first mixture component
  # We can further reduce the precision of mean_cand based on the maximum tolerable error
  # Also, we only need to consider weight <= 0.5, given that the order of two components does
  # not matter.

  mean_cand = unique(round(es/eps)*eps);
  sd_cand = seq(eps,max_s,eps);
  weight_cand = seq(min_w, 0.5, eps);

  if (k > 2 & length(mean_cand) > 20) {
    mean_cand = as.numeric(as.character((data.frame(table(round(es*10)/10)) %>% filter(Freq >= length(es)/25))$Var1));
  }

  first_cut = get_first_combis(data_histogram, mean_cand, sd_cand, weight_cand, equisd = equisd, speed = speed, sdes = sd(es));
  final_selection = get_final_combis(first_cut, es, data_histogram, topx = topx, k = k, input_precision = input_precision, eps = eps, verbose = verbose);

  # we are done if k = 2
  if (k == 2) return(as.data.frame(rbind(as.numeric(final_selection[1:3]), as.numeric(final_selection[4:6]))));

  # All subsequent steps are only needed when k = 3 -----------------------------------------------------------
  # A straightforward idea to solve this case is to enumerate all pairs of candidate components before deriving
  # the third for each pair. The problem, of course, is the computational overhead associated with this idea.
  # To speed up the process, we want to derive as much value from the output of the two-component decomposition
  # as possible. Specifically, we note that, when one applies a two-component decomposition algorithm over a
  # mixture of three components, the likely output composites two of the three components into one of the two
  # output components for two-component decomposition. Realizing this, a natural idea to leverage the output of
  # two-component decomposition is to retain the component with the smaller sd, as this component is less
  # likely to be a composite of two separate components. Then, we deduct this component from data_histogram and
  # call upon another round of two-component decomposition over the remaining histogram. The results of the two
  # rounds are then merged to form our final output.

  if (verbose) {
    print("first component:");
    print(final_selection);
  }

  # we want to select the component with the smaller standard deviation
  selected_component = ifelse(final_selection[5] >= final_selection[2], 1, 4);
  final_selection = final_selection[selected_component:(selected_component+2)];

  # the accuracy of weight is crucial, as we need to deduct this component from data_histogram. Thus, we attempt
  # to more accurately pinpoint the weight of the component here by rescaling it
  temp_w_cand = seq(min_w, 1-min_w, eps);
  wscores = sapply(temp_w_cand, function(x){
    final_selection[3] = x;
    return(get_mixture_distance(data_histogram, as.data.frame(t(final_selection))));
  });
  final_selection[3] = temp_w_cand[which(wscores == min(wscores))[1]];
  if (verbose) { print("rescaled weight:"); print(final_selection); }

  # BETA feature level 3: similar to the weight rescaling above, rescale the sd of the first component
  if (beta_t > 3) {
    sd_scores = sapply(sd_cand, function(x){
      final_selection[2] = x;

      index_range_min = min(which(data_histogram$x > (final_selection[2] - 5*min(sd_cand))));
      index_range_max = min(which(data_histogram$x >= min(max(data_histogram$x), final_selection[2] + 5*min(sd_cand))));
      return(get_mixture_distance(data_histogram[index_range_min:index_range_max,], as.data.frame(t(final_selection)), style="SLI"));
    });
    final_selection[2] = sd_cand[which(sd_scores == min(sd_scores))[1]];
    if (verbose) { print("rescaled sd:"); print(final_selection); }
  }

  # Now deduct the first component from data_histogram. A number of stablization operations are in order,
  # as some of the remaining histogram values tend to be negative (due to the sampling error and any
  # estimation error for the first mixture component. Fixing these extreme values is important to make
  # sure the next call of two-component decomposition does not select extreme components to compensate
  # for the extreme values in the input histogram

  new_histogram = data_histogram;
  deduct_selection = final_selection;
  deduct_selection[2] = deduct_selection[2] + eps/2; deduct_selection[3] = deduct_selection[3] * deduct_selection[2]/(deduct_selection[2] - eps/2);
  new_histogram$p = (data_histogram$p - get_histogram_from_mixture(data_histogram$x, as.data.frame(t(deduct_selection)))$p)/(1-as.numeric(deduct_selection[3]));
  new_histogram$p = pmax(0, new_histogram$p);
  if (resolution) new_histogram = new_histogram[round(nrow(new_histogram)/5):round(4*nrow(new_histogram)/5),];

  # run the two-component decomposition again
  next_first_cut = get_first_combis(new_histogram, mean_cand, sd_cand, weight_cand, speed = speed);
  next_final_selection = get_final_combis(next_first_cut, es, new_histogram,
                                          emean = as.numeric(final_selection[1]),
                                          esd = as.numeric(final_selection[2]),
                                          ew = as.numeric(final_selection[3]),
                                          k = 2,
                                          input_precision = input_precision,
                                          eps = eps,
                                          verbose = verbose);
  if (verbose) { print("second and third components:"); print(next_final_selection); }

  # merge all three components by adjusting the weights of the latter two components according to that of the first
  next_final_selection[3] = next_final_selection[3] * (1 - final_selection[3]);
  next_final_selection[6] = next_final_selection[6] * (1 - final_selection[3]);
  cur_outcome = rbind(as.numeric(final_selection[1:3]),
                      as.numeric(next_final_selection[1:3]),
                      as.numeric(next_final_selection[4:6]));
  cur_outcome = cur_outcome[order(cur_outcome[,1]),];

  # BETA feature level 5: rescaling the sds of all components
  if (beta_t > 5) {
    for (comp_id in seq(1,3)) {
      sd_scores = sapply(sd_cand, function(x){
        cur_outcome[comp_id,2] = x;
        index_range_min = min(which(data_histogram$x > (cur_outcome[comp_id,1] - 5*min(sd_cand))));
        index_range_max = min(which(data_histogram$x >= min(max(data_histogram$x), cur_outcome[comp_id,1] + 5*min(sd_cand))));
        return(get_mixture_distance(data_histogram[index_range_min:index_range_max,], as.data.frame(cur_outcome), style="LI"));
      });
      cur_outcome[comp_id,2] = sd_cand[which(sd_scores == min(sd_scores))[1]];
    }
  }

  # BETA feature level 3: rescaling the weights of all components
  if (beta_t > 3) {
    for (comp_id in seq(comp_id,3)) {
      wscores = sapply(temp_w_cand, function(x){
        cur_outcome[comp_id,3] = x;
        return(get_mixture_distance(data_histogram, as.data.frame(cur_outcome)));
      });
      cur_outcome[comp_id,3] = temp_w_cand[which(wscores == min(wscores))[1]];
    }
    cur_outcome[,3] = cur_outcome[,3]/sum(cur_outcome[,3]);
  }

  # Fine tuning for additional resolution, needed when most reported es are within a tight interval, e.g., 0-0.5
  # The basic idea of fine tuning is to test whether adjusting the critical parameters by half of eps would boost
  # the accuracy of the mixture composition. Given that we have eight free parameters, testing 3^8 possible
  # combinations incurs high computational overhead. Thus, we simplify the process into two steps. First, we move
  # simultanously the means of all three components and the sds of the first and third components by half of eps
  # on either direction. Then, for each case, we test all possible combinations of the standard deviation and the
  # weight for the second component. The reason is that the sd and weight of the second component is 1) critical
  # for our subsequent moderator analysis (which only considers the affiliation for the second component), and 2)
  # difficult to estimate accurately by the algorithm, given that it is usually corresponding to the es with the
  # highest probability density in data_histogram.

  if (resolution) {
    mx = cur_outcome;
    mx = mx[order(mx[,1]),];

    row_ids = c(1,1,2,3,3); col_ids = c(1,2,1,1,2);
    backup_mx = mx; temp_mx = mx; opt_mx = mx;
    temp_mx_score = get_mixture_distance(data_histogram, mx);

    for (rev_indx in seq(1, length(row_ids))) {
      temp_mx[row_ids[rev_indx], col_ids[rev_indx]] = temp_mx[row_ids[rev_indx], col_ids[rev_indx]] - eps/2;
    }
    for (sdv in seq(max(mx[2,2]/3,eps/2),mx[2,2]*3,eps/2)) {
      for (wtv in seq(eps/2,1-eps/2,eps/2)) {
        rnd_mx = temp_mx;
        rnd_mx[2,2] = sdv; rnd_mx[2,3] = wtv;
        rnd_mx[1,3] = (1-wtv)*temp_mx[1,3]/(temp_mx[1,3]+temp_mx[3,3]);
        rnd_mx[3,3] = (1-wtv)*temp_mx[1,3]/(temp_mx[3,3]+temp_mx[3,3]);
        rnd_score = get_mixture_distance(data_histogram, rnd_mx);
        if (rnd_score < temp_mx_score) {
          temp_mx_score = rnd_score;
          opt_mx = rnd_mx;
        }
      }
    }

    temp_mx = backup_mx;
    for (rev_indx in seq(1, length(row_ids))) {
      temp_mx[row_ids[rev_indx], col_ids[rev_indx]] = temp_mx[row_ids[rev_indx], col_ids[rev_indx]] + eps/2;
    }
    for (sdv in seq(max(mx[2,2]/3,eps/2),mx[2,2]*3,eps/2)) {
      for (wtv in seq(eps/2,1-eps/2,eps/2)) {
        rnd_mx = temp_mx;
        rnd_mx[2,2] = sdv; rnd_mx[2,3] = wtv;
        rnd_mx[1,3] = (1-wtv)*temp_mx[1,3]/(temp_mx[1,3]+temp_mx[3,3]);
        rnd_mx[3,3] = (1-wtv)*temp_mx[1,3]/(temp_mx[3,3]+temp_mx[3,3]);
        rnd_score = get_mixture_distance(data_histogram, rnd_mx);
        if (rnd_score < temp_mx_score) {
          temp_mx_score = rnd_score;
          opt_mx = rnd_mx;
        }
      }
    }

    cur_outcome = opt_mx;
  }

  return(as.data.frame(cur_outcome));
}
