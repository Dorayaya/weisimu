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

# this function plots the MSE for both mu_0_hat and mu_0_hat given sample size range (type =1), shape range (type =2), and trimming proportion p (type =3)
plotmse <- function(x, y1, y2, type){
  xlabel = c("Sample Size", "Shape", "p")
  xlab = xlabel[type]
  plot (x, y1, ylim = c(min(c(y1, y2)), max(c(y1, y2))), type="l",xlab= xlab, ylab="MSE", main=expression(MSE ~simulation ~plots ~"for" ~"estimators:"), col="red")
  points (x, y2, type="l", col="blue")

  legend("bottomright", c(expression(hat(mu)[0]), expression(hat(mu)[p])), lty=c(7,7), col = c("red", "blue"),
         inset=c(0,1), xpd=TRUE, horiz=TRUE, bty="n")
}

# this function returns 2 MSE vectors along with a pre-set sample size range. The MSEs are respectively for mu_0_hat and mu_0_hat. User can choose shape, scale, p, S (# of simulations), but not sample size range (set to 20:2000).
simtrim_ss = function(n=NA, shape, scale, p, S) {
  min_ss = 20
  max_ss = 2000
  ss = min_ss:max_ss
  true_mean_Vect = c()
  mse_0_vect = c()
  mse_p_vect = c()

  for (i in ss) {
    true_mean = NA
    mu_0_sims = c()
    mu_p_sims = c()

    for (s in 1:S) {
      Y_s = rweibull(n=i, shape=shape, scale=scale)
      mu_0_sims[s] = mean(Y_s)
      mu_p_sims[s] = mean(Y_s, trim=p, trim.upper=TRUE)
    }
    true_mean=scale*gamma(1 + 1/i)
    true_mean_Vect[which(ss==i)] = true_mean
    mse_0_vect [which(ss==i)] = (mean(mu_0_sims) - true_mean)^2 + var(mu_0_sims)
    mse_p_vect[which(ss==i)] = (mean(mu_p_sims) - true_mean)^2 + var(mu_p_sims)
  }
  return(list(ss, mse_0_vect, mse_p_vect))
}

# this function returns 2 MSE vectors along with a pre-set shape range. The MSEs are respectively for mu_0_hat and mu_0_hat. User can choose sample size, scale, p, S (# of simulations), but not shape range.

simtrim_shp = function(n, shape=NA, scale, p, S) {
  min_shape = 0.05
  max_shape = 1
  shape = seq(min_shape,max_shape, 0.01)
  true_mean_Vect = c()
  mse_0_vect = c()
  mse_p_vect = c()

  for (i in shape) {
    true_mean = NA
    mu_0_sims = c()
    mu_p_sims = c()

    for (s in 1:S) {
      Y_s = rweibull(n=n, shape=i, scale=scale)
      mu_0_sims[s] = mean(Y_s)
      mu_p_sims[s] = mean(Y_s, trim=p, trim.upper=TRUE)
    }
    true_mean=scale*gamma(1 + 1/i)
    true_mean_Vect[which(shape==i)] = true_mean
    mse_0_vect [which(shape==i)] = (mean(mu_0_sims) - true_mean)^2 + var(mu_0_sims)
    mse_p_vect[which(shape==i)] = (mean(mu_p_sims) - true_mean)^2 + var(mu_p_sims)
  }
  return(list(shape, mse_0_vect, mse_p_vect))
}

# this function returns 2 MSE vectors along with a pre-set p range. The MSEs are respectively for mu_0_hat and mu_0_hat. User can choose sample size, scale, S (# of simulations), but not p range.

simtrim_p = function(n, shape, scale, p=NA, S) {
  min_p = 0.01
  max_p = 0.20
  p = seq(min_p,max_p, 0.01)
  true_mean_Vect = c()
  mse_0_vect = c()
  mse_p_vect = c()

  for (i in p) {
    true_mean = NA
    mu_0_sims = c()
    mu_p_sims = c()

    for (s in 1:S) {
      Y_s = rweibull(n=n, shape=shape, scale=scale)
      mu_0_sims[s] = mean(Y_s)
      mu_p_sims[s] = mean(Y_s, trim=p, trim.upper=TRUE)
    }
    true_mean=scale*gamma(1 + 1/shape)
    true_mean_Vect[which(p==i)] = true_mean
    mse_0_vect [which(p==i)] = (mean(mu_0_sims) - true_mean)^2 + var(mu_0_sims)
    mse_p_vect[which(p==i)] = (mean(mu_p_sims) - true_mean)^2 + var(mu_p_sims)
  }
  return(list(p, mse_0_vect, mse_p_vect))
}
