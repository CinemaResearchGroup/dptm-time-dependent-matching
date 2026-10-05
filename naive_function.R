rm(list=ls())

#### true_phi is computed from coxme(Surv(Y_PE, event_E) ~ treated + (1 | pair.index), data = stacked) in true_value_function

## treatment ratio
# TR: 20%; baseline hazard: 0.3; true value: -0.5163079
# TR: 30%; baseline hazard: 0.45; true value: -0.4967315
# TR: 40%; baseline hazard: 0.7; true value: -0.4900718

estimated_naive_value <- function(true_phi,sample_size,simulation,lambda_z) {
  library(survival)
  library(optmatch)
  library(MASS)
  library(coxme) #used for fitting the frailty model
  require(foreach)
  require(doMC)
  library(dplyr)
  library(tableone)
  library(ipcwswitch)
  library(MatchIt)
  library(dplyr)
  
  n <- sample_size
  nsim <- simulation
  
  # estimated results
  est1 <- rep(NA, 1000)
  
  # standard error
  se1 <- rep(NA, 1000)
  
  # if in CI, coverage
  ci1 <- rep(NA, 1000)
  
  # mse
  mse_ss1 <- rep(NA, 1000)
  
  # bias
  bias_ss1 <- rep(NA, 1000)
  
  sim_id=c()
  
  sample_size <- list()
  
  
  # set.seed setting
  w=1
  a=5
  
  
  for (sim in 1:nsim){
    set.seed(sim*2021+1207*a+w*100)
    
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
    dt$epsilon <- rnorm(2*n, mean=0, sd=sqrt(sigma2_epsilon)) 
    
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
    
    #table(dat$group)                   
    
    beta_e1 <- 0.1
    beta_e2 <- 0.1
    beta_e3 <- 0.1
    lambda_e <- 1.5
    phi <- -0.5
    
    #Simulate "survival" probabilities of being treated
    dat$S_e <- runif(n, 0, 1)
    
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
    
    # define observed time to event outcome and its indicator
    # Y_e=min(T_e,C)
    new.dat$Y_e <- ifelse(new.dat$T_e>new.dat$CensorT, new.dat$CensorT, new.dat$T_e)
    new.dat$observed.event <- ifelse(new.dat$T_e>new.dat$CensorT, 0, 1) # this is same as "event"
    
    # observed time
    new.dat$Y_z <- ifelse(new.dat$T_z>new.dat$Y_e, new.dat$Y_e, new.dat$T_z)
    new.dat$observed.treat <- ifelse(new.dat$T_z>new.dat$Y_e, 0, 1)  # this is same as "treat"
    
    
    # cut into intervals based on t1
    # tstart: time to treatment
    # tstop: observed time to event (accounting for censoring)
    dat1 <- new.dat[,c("id", "X1", "X2", "X3_0", "X3_1", "T_z", "Y_z", "Y_e", "observed.treat", "observed.event")]
    
    # initialize
    newcgd <- tmerge(dat1, dat1, id = id, tstart = 0, tstop = Y_z)
    
    # Add time interval split at 0.5
    newcgd$cutpt <- t1
    newcgd <- tmerge(newcgd, newcgd, id = id, splitpt = tdc(cutpt)) 
    
    # Add observed time indicator
    newcgd <- tmerge(newcgd, dat1, id = id, treat = event(Y_z, observed.treat))
    
    # arrange X3 to each time interval
    newcgd$X3 <- ifelse(newcgd$tstop<=t1,newcgd$X3_0,newcgd$X3_1)
    
    
    # keep cleaned data
    dat_cut <- newcgd[,c("id", "X1", "X2", "X3", "tstart", "tstop", "treat")]
    
    # build cox model
    cox1 <- coxph(Surv(tstart, tstop, treat) ~ X1 + X2 + X3 + cluster(id), data = dat_cut, control = coxph.control(timefix = FALSE))
    summary(cox1)
    
    beta_hat<-matrix(cox1$coefficients, ncol=1)
    
    # --- NAIVE COX MODEL BEFORE MATCHING ---
    
    # Fit logistic regression for propensity score
    ps_model <- glm(observed.treat ~ X1 + X2 + X3_0, family = binomial(), data = dat1)
    
    # Extract predicted propensity scores
    dat1$ps <- predict(ps_model, type = "response")
    
    # Flip the treatment indicator:
    # original A: treated=1, control=0
    dat1$inv.observed.treat <- 1 - dat1$observed.treat   # now controls = 1, treated = 0
    
    # Nearest neighbor matching with replacement
    m.out <- matchit(inv.observed.treat ~ ps,
                     data = dat1,
                     method = "full")   # allows multiple controls to match same treated
    
    matched_dat1 <- match.data(m.out)
    
    # create group variable
    grp <- cumsum(matched_dat1$observed.treat == 1)
    grp[grp == 0] <- NA
    
    # add Group to Dataset
    new.matched_dat1 <- cbind(matched_dat1,group=grp)
    
    # calculate Y_pe
    new.matched_dat1<- new.matched_dat1 %>%
      group_by(group) %>%
      mutate(new.T_z = first(T_z))
    
    # post time to treatment
    new.matched_dat1$Y_pe <- new.matched_dat1$Y_e-new.matched_dat1$new.T_z
    
    new.matched_dat2 <- new.matched_dat1[-which(new.matched_dat1$Y_pe<0),]
    
    # not all treated subjects have a matched control subject which has event time after time to treatment initiation
    matched_filtered <- new.matched_dat2 %>%
      group_by(group) %>%       # or subclass if using matchit()
      filter(n() > 1) %>%       # keep only groups with 2 or more rows
      ungroup()
    
    
    set.seed(123)# for reproducibility
    
    # from each group: keeping all treated subjects and randomly selecting one control
    matched_random <- matched_filtered %>%
      group_by(group) %>%
      # keep all treated
      filter(inv.observed.treat == 0) %>%
      # add one random control per group
      bind_rows(
        matched_filtered %>%
          group_by(group) %>%
          filter(inv.observed.treat == 1) %>%
          slice_sample(n = 1)
      ) %>%
      arrange(group) %>%
      ungroup()
    
    # build cluster cox model
    naive_cox <- coxph(Surv(Y_pe, observed.event) ~ observed.treat + cluster(group), data=matched_random)
    
    est1[sim] <- naive_cox$coefficients
    
    se1[sim] <- sqrt(diag(vcov(naive_cox)))
    
    # Compute z critical value (for 95% CI)
    z_crit <- qnorm(0.975)  # = 1.96
    
    # Manual CI on log(HR) scale
    lower_log1 <- est1[sim] - z_crit * se1[sim]
    upper_log1 <- est1[sim] + z_crit * se1[sim]
    
    ci1[sim]<- ifelse(lower_log1<true_phi&upper_log1>true_phi, 1, 0)
    
    # bias
    bias_ss1[sim] <- est1[sim]-true_phi
    
    # mse
    mse_ss1[sim] <- bias_ss1[sim]^2
    
    
    
  }
  return(list(naive_est=est1,naive_se=se1,naive_ci=ci1,naive_mse=mse_ss1))
}


######################## example
estimated <- estimated_naive_value(true_phi,1000,1000,0.7)
# here, true_phi is the true value, 1000 is the sample size, 1000 is the simulation times, 0.7 is value of lambda_z which means treatment ratio is 40%

###############################################
## ASE
naive_ase <- mean(na.omit(estimated$naive_se))

## bias
naive_bias <- mean(na.omit(estimated$naive_est))-true_phi

## MAE
naive_absolutebias2 <- mean(abs(na.omit(estimated$naive_est)-na.omit(true_phi)))

## MCSD
naive_mcsd <- sd(na.omit(estimated$naive_est))

## MSE
naive_mse <- mean((na.omit(estimated$naive_est)-true_phi)^2)

