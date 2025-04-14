
#--------------------------------------------------------------------------------
# logit and expit functions
#--------------------------------------------------------------------------------

logit <- function(p) {
  out <-  qlogis(p)
  return(out)
}


expit <- function(x) {
  out <-   plogis(x)
  return(out)
}




#--------------------------------------------------------------------------------
# tplogit and tpexpit functions
#--------------------------------------------------------------------------------

tplogit <- function(p, par) {
  
  out <- qtp3(p, mu = 0, par1 = 1, par2 = par, FUN = qlogis, param = "eps")
  
  # if(p < 0.5*(1+par)) out <- (1+par)*log(p/(1-p+par))
  #  else out <-  (1-par)*log((p-par)/(1-p))
  return(out)
}

tpexpit <- function(x, par) {
  
  out <- ptp3(x, mu = 0, par1 = 1, par2 = par, FUN = plogis, param = "eps")
  
  # if(p < 0.5*(1+par)) out <- (1+par)*log(p/(1-p+par))
  #  else out <-  (1-par)*log((p-par)/(1-p))
  return(out)
}

#--------------------------------------------------------------------------------
# Reparameterisations
#--------------------------------------------------------------------------------

# Original parameterisation (alpha,beta) to (mu,phi)
ab2mp <- function(alpha,beta){
  
  mu <- alpha / (alpha + beta)
  phi <- alpha + beta
  
  out <- list(mean = mu, phi = phi)
  return(out)
}


#  (mu,phi) parameterisation to original (alpha,beta) 
mp2ab <- function(mu, phi) {
  alpha <- mu * phi
  beta <- (1 - mu) * phi
  list(a = a, b = b)
  list(alpha = alpha, beta = beta)
}

#--------------------------------------------------------------------------------
# probability density function of the beta distribution 
#--------------------------------------------------------------------------------
# parameterised in terms of the mean (mu) and phi
# mu: mean (0,1)
# phi: positive parameter
# log: logical; if TRUE, probabilities p are given as log(p).

dbetamv <- function(x, mu, phi, log = FALSE) {
  
  logpdf <- lgamma(phi) - lgamma(mu*phi) - lgamma((1-mu)*phi) +
    (mu*phi-1)*log(x) + ((1-mu)*phi-1)*log(1-x)
  
  out <- ifelse(log, logpdf, exp(logpdf))
  return(out)
}


#--------------------------------------------------------------------------------
# MLE for beta regression (logit link currently)
#--------------------------------------------------------------------------------
# y: response variable
# des: design matrix
# init: initial point for maximisation
# method: method for maximisation (nlminb or a method from optim)
# maxit: maximum number of iterations for optimisation

MLEBR <- function(y, des, init, method, control){
  
  y <- as.vector(y)
  des <- as.matrix(des)
  p <- ncol(des)
  init <- as.vector(init)
  
  # negative log-likelihood function
  mloglik <- function(par){
    mu <- plogis(des%*%par[1:p])
    phi <- exp(par[p+1])
    
    mll <- - sum( as.vector(mapply(dbetamv, y, mu, MoreArgs = list(phi = phi, log = TRUE))) )
    
    return(mll)
  }
  
  # Optimisation step
  if (method != "nlminb") {
    OPT <- optim(init, mloglik, control = control)
    DEV <- OPT$value
  }
  if (method == "nlminb") {
    OPT <- nlminb(init, mloglik, control = control)
    DEV <- OPT$objective 
  }
  
  
  
  # output
  loglik <- function(par) -mloglik(par)
  MLE <- c(OPT$par[1:p],exp(OPT$par[p+1]))
  AIC <- 2*DEV + 2*(p+1)
  BIC <- 2*DEV + log(n)*(p+1)
  OUT <- list(loglik = loglik, OPT = OPT, MLE = MLE,
              AIC = AIC, BIC = BIC)
  return(OUT)
}




#--------------------------------------------------------------------------------
# MLE for twopiece beta regression (tplogit link currently)
#--------------------------------------------------------------------------------
# y: response variable
# des: design matrix
# init: initial point for maximisation
# method: method for maximisation (nlminb or a method from optim)
# maxit: maximum number of iterations for optimisation

MLETPBR <- function(y, des, init, method, control){
  
  y <- as.vector(y)
  des <- as.matrix(des)
  p <- ncol(des)
  init <- as.vector(init)
  
  # negative log-likelihood function
  mloglik <- function(par){
    gamma0 <- 2*expit(par[p+2]) - 1
    mu <- tpexpit(des%*%par[1:p], gamma0)
    phi <- exp(par[p+1])
    
    
    mll <- - sum( as.vector(mapply(dbetamv, y, mu, MoreArgs = list(phi = phi, log = TRUE))) )
    
    return(mll)
  }
  
  # Optimisation step
  if (method != "nlminb") {
    OPT <- optim(init, mloglik, control = control)
    DEV <- OPT$value
  }
  if (method == "nlminb") {
    OPT <- nlminb(init, mloglik, control = control)
    DEV <- OPT$objective 
  }
  
  
  
  # output
  loglik <- function(par) -mloglik(par)
  MLE <- c(OPT$par[1:p],exp(OPT$par[p+1]), 2*expit(OPT$par[p+2]) - 1)
  AIC <- 2*DEV + 2*(p+2)
  BIC <- 2*DEV + log(n)*(p+2)
  OUT <- list(loglik = loglik, OPT = OPT, MLE = MLE,
              AIC = AIC, BIC = BIC)
  return(OUT)
}
