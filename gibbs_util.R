library(mvtnorm)
library(tmvtnorm)
library(pracma)
require(reshape2)
library(knitr)
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
  unique_params = unique(plot_df$Parameter)
  
  acf.key<-list()
  for(i in 1:length(unique_params)){
    acf.key[[i]]<-acf(plot_df$Estimate[plot_df$Parameter==unique_params[[i]]], 
                      plot=FALSE)
  }
  
  names(acf.key) = unique_params
  autoplot(acf.key, nbin=5)
}

summarize_dist <- function(distribution, param_names, title='', round_places=4) {
  table = data.frame()
  
  for (i in seq(ncol(distribution))) {
    param_sample = distribution[, i]
    param_summary = c(mean(param_sample), sd(param_sample), 
                      quantile(param_sample, 0.025),
                      quantile(param_sample, 0.975))
    
    table = rbind(table, param_summary)
  }
  
  table = cbind(param_names, table)
  table = round_df(table, round_places)
  
  kable(table, 
        col.names=c('Parameter', 'Post. Mean', 'Post. Sd', 
                    '95% CI Low', '95% CI High'),
        caption=title
  )
}

round_df <- function(x, digits) {
  # round all numeric variables
  # x: data frame 
  # digits: number of digits to round
  numeric_columns <- sapply(x, mode) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}

niter <- 10000
y <- cost_y
group <- as.factor(burrito$Location)
  
z <- model.matrix(~group-1)

w <- cbind(1, X_proteins, z)

q <- ncol(z)
p <- ncol(w)-q
n <- nrow(w)
  
beta_keep <- matrix(NA, niter, p)
gamma_keep <- matrix(NA, niter, q)
sig2inv_keep <- rep(NA, niter)
kappa2inv_keep <- rep(NA, niter)
  
a1 <- 0.1975
a2 <- 0.44
b1 <- 0.5
b2 <- 0.5
 
lb <- 0
ub <- Inf
#starting values
beta <- rnorm(p, 0, 2)
gamma <- rnorm(q, 0, 1)
theta <- c(beta, gamma)
sig2inv <- rgamma(1, a1, a2)
kappa2inv <- rgamma(1, b1, b2)
    
tau2 <- 4
sigdiag <- c(0, rep(1/tau2, p-1), rep(kappa2inv,1))
  
for(i in 1:niter){
  #update beta
  v <- t(w)%*%w * sig2inv
  sigdiag[(p+1):(p+q)] <- kappa2inv
  diag(v) <- diag(v) + sigdiag
  v <- chol2inv(chol(v))
    
  m <- v %*% (sig2inv*t(w)%*%y)
    
  theta <- drop(m + t(chol(v)) %*% rtmvnorm(n=p+q, mean=0, sigma=1, lower = lb, upper=ub, algorithm = "gibbs"))
    
  #update sig2inv
  sig2inv <- rgamma(1, a1 + n/2, a2 + 0.5*sum((y-w%*%theta)^2))
    
  #update kappa2inv
  kappa2inv <- rgamma(1, b1 + q/2, b2 + sum(theta[(p+1):(p+q)]^2))
    
  #store output
  sig2inv_keep[i] <- sig2inv
  kappa2inv_keep[i] <- kappa2inv
  beta_keep[i,] <- theta[1:p]
  gamma_keep[i,] <- theta[(p+1):(p+q)]
    
}

kap2 <- 1/kappa2inv_keep
sigma1 <- 1/sqrt(sig2inv_keep)

res <- cbind(beta_keep, sigma1, kap2)
colnames(res) <- c("Intercept", "Chicken", "Beef", "Pork", "Shrimp", "Other", "Breakfast", "Sigma", "Kappa2")
dim(res)


dic<-function(x,y,beta,sig2){
  x<-cbind(1,x)
  niter<-ncol(x)
  dbtheta<-0
  xbeta<-x%*%t(beta)
  for(i in 1:niter){
    dbtheta<--2*sum(dnorm(y,xbeta[,i],sig2),log=TRUE)
  }
  dbtheta<-dbtheta/niter
  dtheta_b<--2*sum(dnorm(y,rowMeans(xbeta),mean(sig2)),log=TRUE)
  return(2*dtheta_b-dbtheta)
}
