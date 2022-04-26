library(mvtnorm)
library(tmvtnorm)
library(pracma)
require(reshape2)
library(knitr)
library(forecast)
library(patchwork)
library(ggfortify)


mlr_gibbs <- function(X_raw, y, mu, tau_2, a, b, iterations=10000) {
  
  X = cbind(1, X_raw)
  
  # calculate needed quantities for sampling
  tau_2_inv = 1 / tau_2
  n = length(y)
  p = ncol(X)
  
  # create matrix to store posterior samples
  posterior_dist = matrix(NA, nrow=iterations, ncol=p+1)
  
  # set starting values for Gibbs sampler to be samples from prior
  beta = t(rmvnorm(1, mean=mu, sigma=tau_2 * eye(p)))
  sig_2_inv = rgamma(1, a, b)
  
  posterior_dist[1, ] = c(beta, sig_2_inv)
  
  # run Gibbs sampling MCMC algo
  for (i in seq(2, iterations)) {
    v = solve(sig_2_inv * t(X) %*% X + tau_2_inv * eye(p))
    m = v %*% (sig_2_inv * t(X) %*% y + tau_2_inv * mu)
    beta = t(rmvnorm(1, m, v))
    
    sig_2_inv = rgamma(1, a + n/2, b + 0.5 * t(y-X %*% beta) %*% (y-X %*% beta))
    
    posterior_dist[i, ] = c(beta, sig_2_inv)
  }
  
  # provide sigma as a parameter in the posterior, remove sig_2_inv
  sigma_column = ncol(posterior_dist)
  posterior_dist[, sigma_column] = 1 / sqrt(posterior_dist[,sigma_column])
  
  colnames(posterior_dist) = c('Intercept', colnames(X_raw), 'sigma')
  
  return(posterior_dist)
}

truncated_gibbs <- function(X_raw, y, mu, tau_2, a, b, iterations=10000, lb=0, ub=Inf){
  X <- cbind(1,X_raw)

  tau_2_inv <- 1 / tau_2
  n <- length(y)
  p <- ncol(X)

  posterior_dist <- matrix(NA, nrow=iterations, ncol=p+1)

  beta <- drop(rtmvnorm(n=1, mean=mu, sigma=tau_2 * eye(p), lower=lb, upper=ub, algorithm="gibbs"))
  sig_2_inv <- rgamma(1, a, b)

  posterior_dist[1, ] <- c(beta, sig_2_inv)

  XtX <- t(X) %*% X
  Xty <- t(X) %*% y
  I <- diag(p)

  for ( i in seq(2, iterations)){
    v <- solve(sig_2_inv * XtX + tau_2_inv * I)
    m <- drop(v %*% (sig_2_inv * Xty + tau_2_inv * mu))
    beta <- drop(rtmvnorm(n=1, mean=m, sigma=v, lower=lb, upper=ub, algorithm="gibbs"))

    sig_2_inv <- rgamma(1, a+n/2, b+ 0.5 * t(y-X %*% beta) %*% (y- X %*% beta))

    posterior_dist[i, ] <- c(beta, sig_2_inv)
  }

  sigma_column <- ncol(posterior_dist)
  posterior_dist[, sigma_column] <- 1/sqrt(posterior_dist[,sigma_column])

  colnames(posterior_dist) <- c('Intercept', colnames(X_raw), 'sigma')

  return(posterior_dist)
}

plot_traces <- function(gibbs_distribution, title='', facets=TRUE) {
  
  plot_df = melt(gibbs_distribution)
  colnames(plot_df) = c('Iteration', 'Parameter', 'Estimate')
  
  if(facets){
    # or plot on different plots
    ggplot(plot_df, aes(Iteration, Estimate)) + geom_line() + 
      facet_grid(Parameter ~ ., scales='free_y') + labs(title=title)
  } else {
    # plot on same grid, each series colored differently -- 
    # good if the series have same scale
    ggplot(plot_df, aes(Iteration, Estimate)) + 
      geom_line(aes(colour = Parameter)) + labs(title=title)
  }
}

acf_plots <- function(dist){
  plot_df <- melt(dist)
  colnames(plot_df) <- c('Iteration','Parameter','Estimate')
  acf.key<-list()
  for(i in 1:length(unique(plot_df$Parameter))){
    acf.key[[i]]<-acf(plot_df$Estimate[plot_df$Parameter==unique(plot_df$Parameter)[[i]]])
  }
  autoplot(acf.key, main="ACF Plots", nbin=5)
}

summarize_dist <- function(distribution, param_names, title='') {
  table = data.frame()
  
  for (i in seq(ncol(distribution))) {
    param_sample = distribution[, i]
    param_summary = c(mean(param_sample), sd(param_sample), 
                      quantile(param_sample, 0.025),
                      quantile(param_sample, 0.975))
    
    table = rbind(table, param_summary)
  }
  
  table = cbind(param_names, table)
  
  kable(table, 
        col.names=c('Parameter', 'Post. Mean', 'Post. Sd', 
                    '95% CI Low', '95% CI High'),
        caption=title
  )
}
