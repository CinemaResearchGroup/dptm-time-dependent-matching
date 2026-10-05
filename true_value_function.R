rm(list=ls())

true_value <- function(lambda_z,sample_size) {
  library(survival)
  library(optmatch)
  library(MASS)
  library(coxme) #used for fitting the frailty model
  require(foreach)
  require(doMC)
  library(dplyr)
  library(tableone)
  library(dplyr)
  
  set.seed(123)
  n <- sample_size
  #Each person has his/her X3 records on two time points 
  n_record <- rep(2, n)
  t_f <- rep(0.5*(1:2 -1), n) #X3 checking/changing points: 0, 0.5
  
  #create fixed covariates X1 and X2  
  dt <- data.frame(id=rep(1:n, times=n_record), t_f=t_f)
  X1 <- rnorm(n, mean=1, sd=1) #continous covariate X1
  X2 <- rbinom(n, size=1, prob=0.6) #binary covariate X2
  
  # Add fixed covariate to the time intervals for each patient
  dt$X1 <- rep(X1, times=n_record)  
  dt$X2 <- rep(X2, times=n_record)
  
  # Generate time-varying (piecewise) covariate X3
  alpha0 <- 2
  alpha1 <- 2
  
  #Generate b1, b2 for each individual
  b_sigma <- matrix(c(0.1, 0.05, 0.05, 0.025), byrow=T, nrow=2) 
  b <- mvrnorm(n, mu=c(0, 0), Sigma=b_sigma)
  #Simulate random error for each observation
  sigma2_epsilon <- 0.05
  dt$epsilon <- rnorm(2*n, mean=0, sd=sqrt(sigma2_epsilon)) ####
  
  dt$b1 <- rep(b[,1], times=n_record)
  dt$b2 <- rep(b[,2], times=n_record)
  
  #Create X3 based on linear mixed effect model
  dt$X3 <- with(dt, alpha0 + alpha1*t_f + b1  + b2*t_f + epsilon) 
  
  # stay with a clean data
  dt <- dt[,-c(5:7)]
  
  ##################################
  ## time to treatment initiation ##
  ##################################
  ########################## compute time-to-treatment based on survival probability
  beta_z1 <-  0.3 # decrease this, will help increase the number of controls for time-to-treatment. will try 0.5, 1. 
  beta_z2 <- 0.3 # will try 0.5, 1.
  beta_z3 <- 0.2 #decrease this, or increase beta will help increase T_z
  #lambda_z <- 0.15 #decrease this will help increase T_z
  
  #Simulate "survival" probabilities of being treated
  dt$S_z <- rep(runif(n, 0, 1), times=n_record) #from an uniform distribution
  #Transform data from long to wide format
  dat <- reshape(dt, idvar = "id", timevar = "t_f", v.names = "X3", direction = "wide")
  names(dat)[5:6] <- c("X3_0", "X3_1")
  
  #Based on the formulas of T_Z in pdf, calculate time to treatment for each person
  t1 <- 0.5
  
  # used to identify the time intervals
  # -log(S) in (-infinity,c1)
  cond1 <- -log(dat$S_z)<=(lambda_z*exp(beta_z1*dat$X1 + beta_z2*dat$X2 + beta_z3*dat$X3_0)*t1)
  
  # -log(S) in [c1,infinity)
  cond2 <- -log(dat$S_z)>(lambda_z*exp(beta_z1*dat$X1 + beta_z2*dat$X2 + beta_z3*dat$X3_0)*t1)
  
  # identify the time to treatment initiation
  dat$T_z <- -1
  dat$T_z[cond1] <- (-log(dat$S_z)/(lambda_z*exp(beta_z1*dat$X1 + beta_z2*dat$X2 + beta_z3*dat$X3_0)))[cond1]
  dat$T_z[cond2] <- (( -log(dat$S_z)-lambda_z*exp(beta_z1*dat$X1 + beta_z2*dat$X2 + beta_z3*dat$X3_0)*t1+lambda_z*exp(beta_z1*dat$X1 + beta_z2*dat$X2 + beta_z3*dat$X3_1)*t1)/(lambda_z*exp(beta_z1*dat$X1 + beta_z2*dat$X2 + beta_z3*dat$X3_1)))[cond2]
  
  
  ###################
  ## time to event ##
  ###################
  # dependent on T_z
  dat$group <- ifelse(dat$T_z<t1,"group1","group2")
  
  beta_e1 <- 0.1
  beta_e2 <- 0.1
  beta_e3 <- 0.1
  lambda_e <- 1.5
  phi <- -0.5
  
  #Simulate "survival" probabilities of being treated
  dat$S_e <- runif(nrow(dat), 0, 1)
  
  # case 1: if T_z < t1
  # T_z is actually the tau_4 
  # x4(t)=0, if t<=T_z; x4(t)=1, if t>T_z
  group1 <- dat[which(dat$group=="group1"),]
  cond1e.group1 <- (-log(group1$S_e)<lambda_e*exp(beta_e1*group1$X1 + beta_e2*group1$X2 + beta_e3*group1$X3_0)*group1$T_z)
  cond2e.group1 <- ((-log(group1$S_e)>=lambda_e*exp(beta_e1*group1$X1 + beta_e2*group1$X2 + beta_e3*group1$X3_0)*group1$T_z) & 
                      (-log(group1$S_e)< lambda_e*exp(beta_e1*group1$X1 + beta_e2*group1$X2 + beta_e3*group1$X3_0)*group1$T_z + 
                         lambda_e*exp(beta_e1*group1$X1 + beta_e2*group1$X2 + beta_e3*group1$X3_0 + phi)*(t1-group1$T_z)))
  cond3e.group1 <- (-log(group1$S_e)>=lambda_e*exp(beta_e1*group1$X1 + beta_e2*group1$X2 + beta_e3*group1$X3_0)*group1$T_z + 
                      lambda_e*exp(beta_e1*group1$X1 + beta_e2*group1$X2 + beta_e3*group1$X3_0 + phi)*(t1-group1$T_z))
  
  group1$T_e <- 1
  group1$T_e[cond1e.group1] <- ((-log(group1$S_e))/(lambda_e*exp(beta_e1*group1$X1+beta_e2*group1$X2+beta_e3*group1$X3_0)))[cond1e.group1]
  group1$T_e[cond2e.group1] <- ((-log(group1$S_e) - lambda_e*exp(beta_e1*group1$X1+beta_e2*group1$X2+beta_e3*group1$X3_0)*group1$T_z +
                                   lambda_e*exp(beta_e1*group1$X1+beta_e2*group1$X2+beta_e3*group1$X3_0+phi)*group1$T_z)/(lambda_e*exp(beta_e1*group1$X1+beta_e2*group1$X2+beta_e3*group1$X3_0+phi)))[cond2e.group1]
  group1$T_e[cond3e.group1] <- ((-log(group1$S_e) - lambda_e*exp(beta_e1*group1$X1+beta_e2*group1$X2+beta_e3*group1$X3_0)*group1$T_z -
                                   lambda_e*exp(beta_e1*group1$X1+beta_e2*group1$X2+beta_e3*group1$X3_0+phi)*(t1-group1$T_z) +
                                   lambda_e*exp(beta_e1*group1$X1+beta_e2*group1$X2+beta_e3*group1$X3_1 + phi)*t1)/(lambda_e*exp(beta_e1*group1$X1+beta_e2*group1$X2+beta_e3*group1$X3_1+phi)))[cond3e.group1]
  
  # case 2: if t1 <= T_z 
  group2 <- dat[which(dat$group=="group2"),]
  cond1e.group2 <- (-log(group2$S_e)<lambda_e*exp(beta_e1*group2$X1 + beta_e2*group2$X2 + beta_e3*group2$X3_0)*t1)
  cond2e.group2 <- ((-log(group2$S_e)>=lambda_e*exp(beta_e1*group2$X1 + beta_e2*group2$X2 + beta_e3*group2$X3_0)*t1) & 
                      (-log(group2$S_e)< lambda_e*exp(beta_e1*group2$X1 + beta_e2*group2$X2 + beta_e3*group2$X3_0)*t1 + 
                         lambda_e*exp(beta_e1*group2$X1 + beta_e2*group2$X2 + beta_e3*group2$X3_1)*(group2$T_z-t1)))
  cond3e.group2 <- (-log(group2$S_e)>=lambda_e*exp(beta_e1*group2$X1 + beta_e2*group2$X2 + beta_e3*group2$X3_0)*t1 + 
                      lambda_e*exp(beta_e1*group2$X1 + beta_e2*group2$X2 + beta_e3*group2$X3_1)*(group2$T_z-t1)) 
  
  
  group2$T_e <- 1
  group2$T_e[cond1e.group2] <- ((-log(group2$S_e))/(lambda_e*exp(beta_e1*group2$X1+beta_e2*group2$X2+beta_e3*group2$X3_0)))[cond1e.group2]
  group2$T_e[cond2e.group2] <- ((-log(group2$S_e) - lambda_e*exp(beta_e1*group2$X1+beta_e2*group2$X2+beta_e3*group2$X3_0)*t1 +
                                   lambda_e*exp(beta_e1*group2$X1+beta_e2*group2$X2+beta_e3*group2$X3_1)*t1)/(lambda_e*exp(beta_e1*group2$X1+beta_e2*group2$X2+beta_e3*group2$X3_1)))[cond2e.group2]
  group2$T_e[cond3e.group2] <- ((-log(group2$S_e) - lambda_e*exp(beta_e1*group2$X1+beta_e2*group2$X2+beta_e3*group2$X3_0)*t1 -
                                   lambda_e*exp(beta_e1*group2$X1+beta_e2*group2$X2+beta_e3*group2$X3_1)*(group2$T_z-t1) +
                                   lambda_e*exp(beta_e1*group2$X1+beta_e2*group2$X2+beta_e3*group2$X3_1 + phi)*group2$T_z)/(lambda_e*exp(beta_e1*group2$X1+beta_e2*group2$X2+beta_e3*group2$X3_1+phi)))[cond3e.group2]
  
  
  # combine the groups into a new dataset
  new.dat <- rbind(group1,group2)
  new.dat <- new.dat[order(new.dat$id),]
  
  # generate censoring time
  # non-informative censoring time
  l<- 0.3
  u<- 1
  new.dat$CensorT <- runif(n, l, u)
  
  #####################
  ## identify groups ##
  #####################
  ##### here first define the event indicator based on T_e and censorT; second, define the treatment indicator;
  new.dat$Y_e <- pmin(new.dat$T_e,new.dat$CensorT)
  new.dat$event <- ifelse(new.dat$T_e<new.dat$CensorT,1,0)
  
  new.dat$Y_z <- pmin(new.dat$T_z,pmin(new.dat$T_e,new.dat$CensorT))
  new.dat$treat <- ifelse(new.dat$T_z<pmin(new.dat$T_e,new.dat$CensorT),1,0)
  
  # keep the treatment groups
  newnew.dat <- new.dat[which(new.dat$treat==1),]
  
  #######################################
  ## Counterfactual outcome for treated ##
  #######################################
  # no treated up to T and afterwards
  
  # subset treated individuals
  treated <- newnew.dat 
  
  #Simulate "survival" probabilities of being treated for counterfactual
  treated$S_e_cf <- runif(nrow(treated), 0, 1)
  
  # t < t1 or t >= t1
  cond1_cf <- (-log(treated$S_e_cf) < lambda_e * exp(beta_e1*treated$X1 + beta_e2*treated$X2 + beta_e3*treated$X3_0) * t1)
  cond2_cf <- (-log(treated$S_e_cf) >= lambda_e * exp(beta_e1*treated$X1 + beta_e2*treated$X2 + beta_e3*treated$X3_0) * t1)
  
  # initialize
  treated$T_e_cf <- 1
  treated$T_e_cf[cond1_cf] <- (-log(treated$S_e_cf) / (lambda_e*exp(beta_e1*treated$X1 + beta_e2*treated$X2 + beta_e3*treated$X3_0)))[cond1_cf]
  treated$T_e_cf[cond2_cf] <- ((-log(treated$S_e_cf) - lambda_e*exp(beta_e1*treated$X1 + beta_e2*treated$X2 + beta_e3*treated$X3_0)*t1 +
                                  lambda_e*exp(beta_e1*treated$X1 + beta_e2*treated$X2 + beta_e3*treated$X3_1)*t1) /
                                 (lambda_e*exp(beta_e1*treated$X1 + beta_e2*treated$X2 + beta_e3*treated$X3_1)))[cond2_cf]
  
  # observed time to event for treated group and counterfactual event time
  treated$Y_e_cf <- ifelse(treated$T_e_cf>treated$CensorT, treated$CensorT, treated$T_e_cf)
  
  # event indicator for counterfactual outcomes
  treated$event_cf <- ifelse(treated$T_e_cf<treated$CensorT,1,0) #42795
  
  # keep the subjects when both Y_e_cf and Y_e >= Tz
  new.treated <- treated[which(treated$Y_e>=treated$T_z & treated$Y_e_cf>=treated$T_z),] #42795
  
  # obtain post treatment time to event
  new.treated$Y_pe <- new.treated$Y_e-new.treated$T_z
  new.treated$Y_pe_cf <- new.treated$Y_e_cf-new.treated$T_z
  
  
  # observed post-treatment outcome (Y_PE^1)
  obs <- new.treated %>%
    mutate(arm = "observed",
           treated = 1,
           Y_PE=Y_pe,
           event_E=event)
  
  # counterfactual post-treatment outcome (Y_PE^0)
  cf <- new.treated %>%
    mutate(arm = "counterfactual",
           treated = 0,
           Y_PE=Y_pe_cf,
           event_E=event_cf)
  
  # stack them together
  stacked <- bind_rows(obs, cf)
  
  stacked <- stacked[order(stacked$id),]
  stacked$pair.index<- rep((1:(dim(stacked)[1]/2)), each=2) # assign index for each pair
  
  
  #################### summary of stacking data
  # each treated subject contributes two rows: one observed, one counterfactual
  # then fit a cox frailty model with frailty by id
  
  cox_true <- coxme(Surv(Y_PE, event_E) ~ treated + (1 | pair.index), data = stacked)
  
  re <- summary(cox_true)
  
  treated.ratio <- nrow(newnew.dat)/nrow(new.dat)
  
  dat <- data.frame(true=re$coefficients[1],ratio=treated.ratio)
  print(dat)
}

## example
true_value(1.1,100000) # 1.1 means treatment ratio50%, 100000 is the sample size

## treatment ratio
# TR: 20%; baseline hazard: 0.3; true value: -0.5163079
# TR: 30%; baseline hazard: 0.45; true value: -0.4967315
# TR: 40%; baseline hazard: 0.7; true value: -0.4900718
# TR: 50%; baseline hazard: 1.1; true value: -0.4937465
# TR: 60%; baseline hazard: 1.55; true value: -0.5002375





