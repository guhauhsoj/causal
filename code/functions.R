

# Regression outcome

Regression_estimate <- function(z, y, x) {
  x <- scale(x)
  interact <- x * z
  n <- ncol(interact)
  
  new_covars <- as.matrix(cbind(z, interact, x))
  x <- as.matrix(x)
  
  regression_ate <- lm(y ~ new_covars)
  
  coef_treatment <- regression_ate$coefficients[2]
  coef_interact <- regression_ate$coefficients[3:(n + 2)]
  
  
  ate.reg <- coef_treatment + coef_interact %*% colMeans(x)
  
  return(ate.reg)
  
}


Regression_bootstrap <- function(z, y, x, n.boot = 500) {
  point.est  = Regression_estimate(z, y, x)
  
  ## nonparametric bootstrap
  n.sample   = length(y)
  x          = as.matrix(x)
  boot.est   = replicate(n.boot,
                         {
                           id.boot = sample(1:n.sample, n.sample, replace = TRUE)
                           Regression_estimate(z[id.boot], y[id.boot], x[id.boot,])
                         })
  
  
  boot.se    = sd(boot.est)
  
  res = rbind(point.est, boot.se)
  rownames(res) = c("est", "se")
  colnames(res) = c("Lins.est")
  
  return(res)
}


# IPW estimators 


ipw.est = function(z, y, x, truncpscore = c(0, 1))
{
  ## fitted propensity score
  pscore   = glm(z ~ x, family = binomial)$fitted.values
  pscore   = pmax(truncpscore[1], pmin(truncpscore[2], pscore))
  
  ace.ipw0 = mean(z*y/pscore - (1 - z)*y/(1 - pscore))
  ace.ipw  = mean(z*y/pscore)/mean(z/pscore) - 
    mean((1 - z)*y/(1 - pscore))/mean((1 - z)/(1 - pscore))
  
  return(c(ace.ipw0, ace.ipw))     
}


ipw.boot = function(z, y, x, n.boot = 500, truncpscore = c(0, 1))
{
  point.est  = ipw.est(z, y, x, truncpscore)
  
  ## nonparametric bootstrap
  n.sample   = length(z)
  x          = as.matrix(x)
  boot.est   = replicate(n.boot, 
                         {id.boot = sample(1:n.sample, n.sample, replace = TRUE)
                         ipw.est(z[id.boot], y[id.boot], x[id.boot, ], truncpscore)})
  boot.se    = apply(boot.est, 1, sd)
  
  res = rbind(point.est, boot.se)
  rownames(res) = c("est", "se")
  colnames(res) = c("HT", "Hajek")
  
  return(res)
}

# Propensity score stratification


Neyman_SRE = function(z, y, x)
{
  xlevels = unique(x)
  K       = length(xlevels)
  PiK     = rep(0, K)
  TauK    = rep(0, K)
  varK    = rep(0, K)
  for(k in 1:K)
  {
    xk         = xlevels[k]
    zk         = z[x == xk]
    yk         = y[x == xk]
    PiK[k]     = length(zk)/length(z)
    TauK[k]    = mean(yk[zk==1]) - mean(yk[zk==0])
    varK[k]    = var(yk[zk==1])/sum(zk) + 
      var(yk[zk==0])/sum(1 - zk)
  }
  
  return(c(sum(PiK*TauK), sqrt(sum(PiK^2*varK))))
}

# Doubly Robust


OS_est = function(z, y, x, out.family = gaussian, 
                  truncpscore = c(0, 1))
{
  ## fitted propensity score
  pscore   = glm(z ~ x, family = binomial)$fitted.values
  pscore   = pmax(truncpscore[1], pmin(truncpscore[2], pscore))
  
  ## fitted potential outcomes
  outcome1 = glm(y ~ x, weights = z, 
                 family = out.family)$fitted.values
  outcome0 = glm(y ~ x, weights = (1 - z), 
                 family = out.family)$fitted.values
  
  ## regression imputation estimator
  ace.reg  = mean(outcome1 - outcome0) 
  ## IPW estimators
  ace.ipw0 = mean(z*y/pscore - (1 - z)*y/(1 - pscore))
  ace.ipw  = mean(z*y/pscore)/mean(z/pscore) - 
    mean((1 - z)*y/(1 - pscore))/mean((1 - z)/(1 - pscore))
  ## doubly robust estimator
  res1     = y - outcome1
  res0     = y - outcome0
  ace.dr   = ace.reg + mean(z*res1/pscore - (1 - z)*res0/(1 - pscore))
  
  return(c(ace.reg, ace.ipw0, ace.ipw, ace.dr))     
}


OS_ATE = function(z, y, x, n.boot = 2*10^2,
                  out.family = gaussian, truncpscore = c(0, 1))
{
  point.est  = OS_est(z, y, x, out.family, truncpscore)
  
  ## nonparametric bootstrap
  n.sample   = length(z)
  x          = as.matrix(x)
  boot.est   = replicate(n.boot, 
                         {id.boot = sample(1:n.sample, n.sample, replace = TRUE)
                         OS_est(z[id.boot], y[id.boot], x[id.boot, ], 
                                out.family, truncpscore)})
  
  boot.se    = apply(boot.est, 1, sd)
  
  res        = rbind(point.est, boot.se)
  rownames(res) = c("est", "se")
  colnames(res) = c("reg", "HT", "Hajek", "DR")
  
  return(res)
}

# Covariate Balance



cov_balance_plot <- function(title, Bcheck) {
  dat_balance = data.frame(
    est = Bcheck[1,],
    upper = Bcheck[1,] + 1.96 * Bcheck[2,],
    lower = Bcheck[1,] - 1.96 * Bcheck[2,],
    cov = factor(1:ncol(x))
  )
  ggplot(dat_balance) +
    geom_errorbar(aes(x = cov,
                      ymin = lower,
                      ymax = upper),
                  alpha = 0.6) +
    geom_point(aes(x = cov,
                   y = est),
               alpha = 0.6) +
    geom_hline(aes(yintercept = 0),
               alpha = 0.3) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.y = element_blank()
    ) +
    xlab("balance check based on weighting") +
    ggtitle(paste("truncated covariate balance for", title))
}

