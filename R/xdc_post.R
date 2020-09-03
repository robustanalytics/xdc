#-----------------------------------------------
# final step: call meta-regression
# supports an arbitrary number of components &
#   an arbitrary set of independent variables
# input parameters:
#   - es, the reported effect sizes
#   - ml, the independent variables
#   - mx, the mixture composition (3-col)
#-----------------------------------------------

get_lm = function(es, si, ml, mx, method = "RMA"){
  # sort mixture components by mean, as we are interested in the second component
  mx = mx[order(mx[,1]),];
  # joint_prob returns a vector containing the joint probability of observing es and component id = x
  joint_prob = function(x)(mx[x,3]*dnorm(es,mx[x,1],mx[x,2]));

  if (method == "GRP") {
    dv = joint_prob(2)/rowSums(future.apply::future_sapply(1:nrow(mx), joint_prob));

    # GRP only functions when ml has two possible values
    tml = rep(NA, length(ml)); tml[which(ml == min(ml))] = 0; tml[which(ml == max(ml))] = 1;
    sd_pair = c(sd(dv[which(tml==0)]), sd(dv[which(tml==1)]));
    se_pair = c(sd(dv[which(tml==0)])/sqrt(sum(1-tml)-1), sd(dv[which(tml==1)])/sqrt(sum(tml)-1));
    mn_pair = c(mean(dv[which(tml==0)]), mean(dv[which(tml==1)]));
    pval = 1-pnorm(abs(mn_pair[1]-mn_pair[2])/sqrt(sum(se_pair^2)));
    diff_lb = abs(mn_pair[1]-mn_pair[2])+qnorm(0.025)*sqrt(se_pair[1]^2+se_pair[2]^2);
    diff_ub = abs(mn_pair[1]-mn_pair[2])+qnorm(0.975)*sqrt(se_pair[1]^2+se_pair[2]^2);

    output = as.data.frame(t(c(mn_pair, sd_pair, se_pair, pval, diff_lb, diff_ub)));
    colnames(output) = c("m1","m2","sd1","sd2","se1","se2","p","diff_lb","diff_ub");
    return(output);
  } else if (method == "RMA") {
    return(lm(joint_prob(2)/rowSums(future.apply::future_sapply(1:nrow(mx), joint_prob)) ~ ml));
  }
}
