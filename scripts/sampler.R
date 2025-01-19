library(tidyverse)
library(Hotelling)

# Define (log) conditional posteriors 
## b0
lncondpost_b0 <- function(y_ps, x1_ps, x2_ps, y_nps, x1_nps, x2_nps, b0, b1, b2, m0, s0, v0, a){
  sum(y_ps*(b0 + b1*x1_ps + b2*x2_ps) - log(1+exp(b0 + b1*x1_ps + b2*x2_ps))) + a*(sum(y_nps*(b0 + b1*x1_nps + b2*x2_nps) - log(1+exp(b0 + b1*x1_nps + b2*x2_nps)))) + (-(v0+1)/2)*log(1+((b0-m0)^2)/v0*s0^2) # baseline prior
}

## b1
lncondpost_b1 <- function(y_ps, x1_ps, x2_ps, y_nps, x1_nps, x2_nps, b0, b1, b2, m0, s0, v0, a){
  sum(y_ps*(b0 + b1*x1_ps + b2*x2_ps) - log(1+exp(b0 + b1*x1_ps + b2*x2_ps))) + a*(sum(y_nps*(b0 + b1*x1_nps + b2*x2_nps) - log(1+exp(b0 + b1*x1_nps + b2*x2_nps)))) + (-(v0+1)/2)*log(1+((b1-m0)^2)/v0*s0^2) # baseline prior
}

## b2
lncondpost_b2 <- function(y_ps, x1_ps, x2_ps, y_nps, x1_nps, x2_nps, b0, b1, b2, m0, s0, v0, a){
  sum(y_ps*(b0 + b1*x1_ps + b2*x2_ps) - log(1+exp(b0 + b1*x1_ps + b2*x2_ps))) + a*(sum(y_nps*(b0 + b1*x1_nps + b2*x2_nps) - log(1+exp(b0 + b1*x1_nps + b2*x2_nps)))) + (-(v0+1)/2)*log(1+((b2-m0)^2)/v0*s0^2) # baseline prior
}

# Normal Approximation
## negatives of conditional posteriors for mode estimation
### -ln(conditional posterior) for mode estimation
neg_lncondpost_b0 <- function(y_ps, x1_ps, x2_ps, y_nps, x1_nps, x2_nps, b0, b1, b2, m0, s0, v0, a){
  -lncondpost_b0(y_ps, x1_ps, x2_ps, y_nps, x1_nps, x2_nps, b0, b1, b2, m0, s0, v0, a)
}

neg_lncondpost_b1 <- function(y_ps, x1_ps, x2_ps, y_nps, x1_nps, x2_nps, b0, b1, b2, m0, s0, v0, a){
  -lncondpost_b1(y_ps, x1_ps, x2_ps, y_nps, x1_nps, x2_nps, b0, b1, b2, m0, s0, v0, a)
}

neg_lncondpost_b2 <- function(y_ps, x1_ps, x2_ps, y_nps, x1_nps, x2_nps, b0, b1, b2, m0, s0, v0, a){
  -lncondpost_b2(y_ps, x1_ps, x2_ps, y_nps, x1_nps, x2_nps, b0, b1, b2, m0, s0, v0, a)
}

# 2. derivatives of conditional posterior (to get fisher information)
## b0
derivcondpost_b0 <- function(y_ps, x1_ps, x2_ps, y_nps, x1_nps, x2_nps, b0, b1, b2, m0, s0, v0, a){
  sum(-(exp(b0 + b1*x1_ps + b2*x2_ps))/(1+exp(b0 + b1*x1_ps + b2*x2_ps))^2) + a*sum(-(exp(b0 + b1*x1_nps + b2*x2_nps))/(1+exp(b0 + b1*x1_nps + b2*x2_nps))^2) + (((v0+1)*(b0^2-2*b0*m0+m0^2-v0*s0^2))/(((b0^2)-2*b0*m0+m0^2+v0*s0^2)))
}

## b1
derivcondpost_b1 <- function(y_ps, x1_ps, x2_ps, y_nps, x1_nps, x2_nps, b0, b1, b2, m0, s0, v0, a){
  sum(-(x1_ps*exp(b0 + b1*x1_ps + b2*x2_ps))/(1+exp(b0 + b1*x1_ps + b2*x2_ps))^2) + a*sum(-(x1_nps*exp(b0 + b1*x1_nps + b2*x2_nps))/(1+exp(b0 + b1*x1_nps + b2*x2_nps))^2)+ (((v0+1)*(b1^2-2*b1*m0+m0^2-v0*s0^2))/(((b1^2)-2*b1*m0+m0^2+v0*s0^2)))
}

## b2
derivcondpost_b2 <- function(y_ps, x1_ps, x2_ps, y_nps, x1_nps, x2_nps, b0, b1, b2, m0, s0, v0, a){
  sum(-(x2_ps*exp(b0 + b1*x1_ps + b2*x2_ps))/(1+exp(b0 + b1*x1_ps + b2*x2_ps))^2) + a*sum(-(x2_nps*exp(b0 + b1*x1_nps + b2*x2_nps))/(1+exp(b0 + b1*x1_nps + b2*x2_nps))^2) + (((v0+1)*(b2^2-2*b2*m0+m0^2-v0*s0^2))/(((b2^2)-2*b2*m0+m0^2+v0*s0^2)))
}


# sampler
MHS <- function(PS, NPS, n.iter, b0, b1, b2, v0, m0, s0, n.chains, seed){
  set.seed(seed)
  #initialise outputs
  b0_est <- numeric(n.iter)
  b1_est <- numeric(n.iter)
  b2_est <- numeric(n.iter)
  b0_acc <- numeric(n.iter)
  b1_acc <- numeric(n.iter)
  b2_acc <- numeric(n.iter)
  
  # save initial starting values to start chains with the same values
  b0.init <- b0
  b1.init <- b1
  b2.init <- b2
  
  y_ps <- PS$Y
  y_nps <- NPS$Y
  x1_ps <- PS$X1
  x1_nps <- NPS$X1
  x2_ps <- PS$X2
  x2_nps <- NPS$X2
  
  # get alpha
  ## get MLE for b_hat
  m_ps <- glm(Y~X1+X2, data = PS, family = "binomial")
  m_nps <- glm(Y~X1+X2, data = NPS, family = "binomial")
  
  htest <- hotelling.test(as.matrix(coef(m_ps)), as.matrix(coef(m_nps)), var.equal = FALSE) # perform hotelling t^2 test
  
  a <- htest$pval
  
  output <- list() # initialise output list
  
  for (c in 1:n.chains) {
    # reset starting values
    b0 <- b0.init
    b1 <- b1.init
    b2 <- b2.init
    
    # start sampling
    for (i in 1:n.iter) {
      
      # b0
      ## use optimise to get the mode
      b0_mode <- optimise(f = neg_lncondpost_b0, interval = c(-100, 100), y_ps = y_ps, y_nps = y_nps, b1 = b1, b2 = b2, x1_ps = x1_ps, x1_nps = x1_nps, x2_ps = x2_ps, x2_nps = x2_nps, m0 = m0, s0 = s0, v0 = v0, a = a)$minimum
      
      b0_I <- derivcondpost_b0(y_ps, x1_ps, x2_ps, y_nps, x1_nps, x2_nps, b0, b1, b2, m0, s0, v0, a)
      
      ## sample proposal value
      b0_prop <- rnorm(1, b0_mode, sqrt(-(1/b0_I)))
      
      ## calculate acceptance ratio
      r <- lncondpost_b0(y_ps, x1_ps, x2_ps, y_nps, x1_nps, x2_nps, b0_prop, b1, b2, m0, s0, v0, a)/lncondpost_b0(y_ps, x1_ps, x2_ps, y_nps, x1_nps, x2_nps, b0, b1, b2, m0, s0, v0, a)
      u <- runif(1, 0, 1)
      
      ## decision rule and save if accepted/rejected
      b0_est[i] <- b0 <- ifelse(r >= u, b0_prop, b0)
      b0_acc[i] <- ifelse(r >= u, 1, 0)
      
      # b1
      b1_mode <- optimise(f = neg_lncondpost_b1, interval = c(-100, 100), y_ps = y_ps, y_nps = y_nps, b0 = b0, b2 = b2, x1_ps = x1_ps, x1_nps = x1_nps, x2_ps = x2_ps, x2_nps = x2_nps, m0 = m0, s0 = s0, v0 = v0, a = a)$minimum
      
      b1_I <- derivcondpost_b1(y_ps, x1_ps, x2_ps, y_nps, x1_nps, x2_nps, b0, b1, b2, m0, s0, v0, a)
      
      b1_prop <- rnorm(1, b1_mode, sqrt(-(1/b1_I)))
      
      r <- lncondpost_b1(y_ps, x1_ps, x2_ps, y_nps, x1_nps, x2_nps, b0, b1_prop, b2, m0, s0, v0, a)/lncondpost_b1(y_ps, x1_ps, x2_ps, y_nps, x1_nps, x2_nps, b0, b1, b2, m0, s0, v0, a)
      u <- runif(1, 0, 1)
      
      b1_est[i] <- b1 <- ifelse(r >= u, b1_prop, b1)
      b1_acc[i] <- ifelse(r >= u, 1, 0)
      
      # b2
      b2_mode <- optimise(f = neg_lncondpost_b2, interval = c(-100, 100), y_ps = y_ps, y_nps = y_nps, b0 = b0, b1 = b1, x1_ps = x1_ps, x1_nps = x1_nps, x2_ps = x2_ps, x2_nps = x2_nps, m0 = m0, s0 = s0, v0 = v0, a = a)$minimum
      
      b2_I <- derivcondpost_b2(y_ps, x1_ps, x2_ps, y_nps, x1_nps, x2_nps, b0, b1, b2, m0, s0, v0, a)
      
      b2_prop <- rnorm(1, b2_mode, sqrt(-(1/b2_I)))
      
      r <- lncondpost_b2(y_ps, x1_ps, x2_ps, y_nps, x1_nps, x2_nps, b0, b1, b2_prop, m0, s0, v0, a)/lncondpost_b2(y_ps, x1_ps, x2_ps, y_nps, x1_nps, x2_nps, b0, b1, b2, m0, s0, v0, a)
      u <- runif(1, 0, 1)
      
      b2_est[i] <- b2 <- ifelse(r >= u, b2_prop, b2)
      b2_acc[i] <- ifelse(r >= u, 1, 0)
    }
    # combine and save in DF for output
    df <- data.frame(b0 = b0_est,
                     b1 = b1_est,
                     b2 = b2_est,
                     b0_acc = b0_acc,
                     b1_acc = b1_acc,
                     b2_acc = b2_acc)
    output[[c]] <- df
  }
  
  
  return(output)
}