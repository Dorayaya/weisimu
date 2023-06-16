# allow mean function to trim above
mean = function(x, trim=0, na.rm=FALSE, trim.upper=FALSE, ...) {

  if (!trim.upper) {
    base::mean(x,trim=trim,na.rm=na.rm, ...)
  } else {
    base::mean(x[x<=quantile(x,probs=1-trim,na.rm=na.rm)], na.rm=na.rm, ...)
  }
}

# log likelihood function for weibull (stolen from EnvStats package)
neg_ll = function(theta, y) {
  shape = theta[1]
  scale = theta[2]
  sum(-log(shape) - (shape - 1) * log(y) + shape * log(scale) + ((y/scale)^shape))
}


# mle for mean of weibull  (stolen from EnvStats package)
wei_mu_mle = function(y) {
  theta.hat =
    nlminb(start = c(1,1),
           objective = neg_ll,
           lower = c(0,0),
           y = y)$par
  mu_mle = theta.hat[2] * gamma(1 + 1/theta.hat[1])
  return(mu_mle)
}


# simulation
# n = sample size
# shape = shape parameter of weibull
# scale = scale parameter of weibull
# p = trim top 100p% of distirbution for trimmed mean
# S = number of simulations to conduct
simtrim = function(n, shape, scale, p, S) {

  mu_0_sims = rep(NA,S)
  mu_p_sims = rep(NA,S)
  mu_mle_sims = rep(NA,S)
  for (s in 1:S) {
    Y_s = rweibull(n=n, shape=shape, scale=scale)
    mu_0_sims[s] = mean(Y_s)
    mu_p_sims[s] = mean(Y_s, trim=p, trim.upper=TRUE)
    mu_mle_sims[s] = wei_mu_mle(Y_s)
  }

  true_mean = scale*gamma(1 + 1/shape)
  retlist =
    list(true_mean,
         bias=list(mu_0=mean(mu_0_sims) - true_mean,
                   mu_p=mean(mu_p_sims) - true_mean,
                   mu_mle=mean(mu_mle_sims) - true_mean),
         vars=list(mu_0=var(mu_0_sims),
                   mu_p=var(mu_p_sims),
                   mu_mle=var(mu_mle_sims)),
         MSE=list(mu_0=(mean(mu_0_sims) - true_mean)^2 + var(mu_0_sims),
                  mu_p=(mean(mu_p_sims) - true_mean)^2 + var(mu_p_sims),
                  mu_mle=(mean(mu_mle_sims) - true_mean)^2 + var(mu_mle_sims))
    )

  return(retlist)
}

