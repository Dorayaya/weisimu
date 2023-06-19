# allow mean function to trim above
mean = function(x, trim=0, na.rm=FALSE, trim.upper=FALSE, ...) {

  if (!trim.upper) {
    base::mean(x,trim=trim,na.rm=na.rm, ...)
  } else {
    base::mean(x[x<=quantile(x,probs=1-trim,na.rm=na.rm)], na.rm=na.rm, ...)
  }
}


# simulation
# n = sample size
# shape = shape parameter of weibull
# scale = scale parameter of weibull
# p = trim top 100p% of distirbution for trimmed mean
# S = number of simulations to conduct
simtrim = function(n, shape, scale, p, S) {

  if (length(S) > 1) {
    warning("S has length>1, using only first element")
    S = S[1]
  }

  mu_p_sims = rep(NA,S)
  for (s in 1:S) {
    Y_s = rweibull(n=n, shape=shape, scale=scale)
    mu_p_sims[s] = mean(Y_s, trim=p, trim.upper=TRUE)
  }

  true_mean = scale*gamma(1 + 1/shape)
  true_variance = scale*scale*(gamma(1 + 2/shape) - gamma(1 + 1/shape)^2)
  retlist =
    list(mu_0      = true_mean,
         bias_mu_0 = 0,
         var_mu_0  = true_variance/n,
         MSE_mu_0  = 0^2 + true_variance/n,
         mu_p      = mean(mu_p_sims),
         bias_mu_p = mean(mu_p_sims) - true_mean,
         var_mu_p  = var(mu_p_sims),
         MSE_mu_p  = (mean(mu_p_sims) - true_mean)^2 + var(mu_p_sims)
    )

  return(retlist)
}


# this function plots simulations and also outputs results so you can store them
simtrim_by = function(n, shape, scale, p, S,
                      plot=TRUE, lty=c(2,1), lwd=c(2,2), col=c("black","black"),
                      leg.pos="topright", ...) {
  # check how many params are >1 length
  if ( (length(n)>1) + (length(shape)>1) + (length(scale)>1) + (length(p)>1) > 1 ) {
    stop("only one of (n, shape, scale, p) may be vector of length>1")
  }
  if ( (length(n)==1) & (length(shape)==1) & (length(scale)==1) & (length(p)==1) ) {
    stop("one of (n, shape, scale, p) must be vector of length>1")
  }
  if (length(S) > 1) {
    warning("S has length>1, using only first element")
    S = S[1]
  }

  if (length(n) > 1) {
    allsims = sapply(X=n, FUN=simtrim, shape=shape, scale=scale, p=p, S=S)
    xlab = "n"
  } else if (length(shape) > 1) {
    allsims = sapply(X=shape, FUN=simtrim, n=n, scale=scale, p=p, S=S)
    xlab = "shape"
  } else if (length(scale) > 1) {
    allsims = sapply(X=scale, FUN=simtrim, n=n, shape=shape, p=p, S=S)
    xlab = "scale"
  } else if (length(p) > 1) {
    allsims = sapply(X=p, FUN=simtrim, n=n, shape=shape, scale=scale, S=S)
    xlab = "p"
  }

  if (plot) {
    plot(x=get(xlab),
         y=as.numeric(allsims["MSE_mu_0",]),
         type="l",lty=lty[1],lwd=lwd[1],col=col[1],xlab=xlab,ylab="MSE",
         ylim=c(0,max(as.numeric(allsims["MSE_mu_0",]),as.numeric(allsims["MSE_mu_p",]))),
         ...)
    lines(x=get(xlab),
          y=as.numeric(allsims["MSE_mu_p",]),
          lty=lty[2],lwd=lwd[2],col=col[2])
    legend(leg.pos, lty=lty, lwd=lwd, col=col, legend=c("untrimmed","trimmed"))
  }

  colnames(allsims) = get(xlab)
  return(allsims)
}



