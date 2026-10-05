rm(list=ls())

ipcw.new <- function(data, id, tstart, tstop, cens, arm, bas.cov, conf, trunc=NULL, type = "kaplan-meier"){
  
  tempcall <- match.call()
  
  data$weights <- vector("numeric", length = nrow(data))
  
  levArms <- levels(data[, arm])
  
  # ~ baseline covariates in numerator
  num <- paste(c(bas.cov), collapse = "+")
  # ~ baseline covariates and time-dependent confoundersin the denominator
  denom <- paste(c(bas.cov, conf), collapse = "+")
  
  if (length(levArms) > 2) {
    stop("Method not implemented yet for more than 2 arms!")
  } else if (length(levArms) < 2) {
    # if numerator
    if (length(bas.cov) != 0) {
      # Cox model for stabilized weights
      #    fit.cox.num <- coxph(formula = eval(parse(text = paste("Surv(",
      #                                                          deparse(tempcall$tstart), ", ",
      #                                                         deparse(tempcall$tstop), ", ",
      #                                                        deparse(tempcall$cens), ") ~  ",
      #                                                       num, sep = "" ))),
      #                    data = data)
      fit.cox.denom <- coxph(eval(parse(text = paste("Surv(",
                                                     deparse(tempcall$tstart), ", ",
                                                     deparse(tempcall$tstop), ", ",
                                                     deparse(tempcall$cens), ") ~  ",
                                                     denom, sep = ""))),
                             data = data)
      #  fit.cox.denom<- coxph(Surv(tstart,tstop,cens)~X1+X2+X3, data=data)
      idpat <- unique(data[, id])
      for (i in idpat) {
        datai <- data[data[, id] == i, ]
        km.num<- survfit(Surv(data$tstart, data$tstop, data$cens)~1, data=data, type=type)
        #    Z_cox <- coxph(Surv(tstart, tstop_Z, Z) ~ X1 + X2 + X3, data = dat_cut, cluster=id, control = coxph.control(timefix = FALSE)) 
        
        #   km.num <- survfit(fit.cox.num, newdata = datai, individual = TRUE, type = type)
        km.denom <- survfit(fit.cox.denom, newdata = datai, individual = TRUE, type = type)
        # data$weights[data[, id] == i] <-
        #  1 / summary(km.denom,times = datai[, as.character(tempcall$tstart)])$surv
        data$weights[data[, id] == i] <-
          summary(km.num,times = datai[, 'tstart'])$surv / summary(km.denom,times = datai[, 'tstart'])$surv
        
      }
    } else{
      # if denominator only
      fit.cox.denom <- coxph(eval(parse(text = paste("Surv(", deparse(tempcall$tstart), ", ",
                                                     deparse(tempcall$tstop), ", ",
                                                     deparse(tempcall$cens), ") ~  ",
                                                     denom, sep = ""))), data = data)
      idpat <- unique(data[, id])
      for (i in idpat) {
        datai <- data[data[, id] == i, ]
        km.denom <- survfit(fit.cox.denom, newdata = datai, individual = TRUE, type = type)
        data$weights[data[, id] == i] <-
          1 / summary(km.denom, times = datai[, as.character(tempcall$tstart)])$surv
      }
    }
  } else{
    # if numerator
    if (length(bas.cov) != 0) {
      # For the first arm
      data1 <- data[data[, arm] == levArms[1],]
      #   fit.cox.num.1 <- coxph(formula = eval(parse(text = paste("Surv(",
      #                                                             deparse(tempcall$tstart), ", ",
      #                                                             deparse(tempcall$tstop), ", ",
      #                                                             deparse(tempcall$cens), ") ~  ",
      #                                                             num, sep = ""))),
      #                           data = data1)
      fit.cox.denom.1 <- coxph(eval(parse(text = paste("Surv(",
                                                       deparse(tempcall$tstart), ", ",
                                                       deparse(tempcall$tstop), ", ",
                                                       deparse(tempcall$cens), ") ~  ",
                                                       denom, sep = ""))),
                               data = data1)
      
      idpat1 <- unique(data1[, id])
      for (i in idpat1) {
        datai <- data1[data1[, id] == i, ]
        km.num<- survfit(Surv(data1$tstart, data1$tstop, data1$cens)~1, data=datai, type=type)
        #    Z_cox <- coxph(Surv(tstart, tstop_Z, Z) ~ X1 + X2 + X3, data = dat_cut, cluster=id, control = coxph.control(timefix = FALSE)) 
        
        #   km.num <- survfit(fit.cox.num, newdata = datai, individual = TRUE, type = type)
        # data$weights[data[, id] == i] <-
        #  1 / summary(km.denom,times = datai[, as.character(tempcall$tstart)])$surv
        # data$weights[data[, id] == i] <-
        #  summary(km.num,times = datai[, 'tstart'])$surv / summary(km.denom,times = datai[, 'tstart'])$surv
        # km.num <- survfit(fit.cox.num.1, newdata = datai, individual = TRUE, type = type)
        km.denom <- survfit(fit.cox.denom.1, newdata = datai,type = type, id=id)
        #km.denom<- survfit(fit.cox.denom.1, newdata = datai,type = type)
        denom_surv<- summary(km.denom,
                             times = datai[, as.character(tempcall$tstart)],extend=TRUE)$surv
        num_surv<-  summary(km.num, times = datai[, as.character(tempcall$tstart)], extend=TRUE)$surv
        data1$weights[data1[, id] == i] <- 1 / denom_surv 
        
        
      }
      
      # For the second arm
      data2 <- data[data[, arm] == levArms[2],]
      #    fit.cox.num.2 <- coxph(eval(parse(text = paste("Surv(",
      #                                                  deparse(tempcall$tstart), ", ",
      #                                                 deparse(tempcall$tstop), ", ",
      #                                                deparse(tempcall$cens), ") ~  ",
      #                                               num, sep = ""))), 
      #                      data = data2)
      fit.cox.denom.2<- coxph(Surv(tstart,tstop,cens)~X1+X2+X3, data=data2)
      
      fit.cox.denom.2 <- coxph(eval(parse(text = paste("Surv(",
                                                       deparse(tempcall$tstart), ", ",
                                                       deparse(tempcall$tstop), ", ",
                                                       deparse(tempcall$cens), ") ~  ",
                                                       denom, sep = ""))),
                               data = data2)
      
      idpat2 <- unique(data2[, id])
      for (i in idpat2) {
        datai <- data2[data2[, id] == i, ]
        km.num<- survfit(Surv(data2$tstart, data2$tstop, data2$cens)~1, data=data2, type=type)
        
        # km.num <- survfit(fit.cox.num.2, newdata = datai, individual = TRUE, type = type)
        km.denom <- survfit(fit.cox.denom.2, newdata = datai, individual = TRUE,id=id, type = type)
        denom_surv<- summary(km.denom,
                             times = datai[, as.character(tempcall$tstart)], extend=TRUE)$surv
        num_denom<- summary(km.num, times = datai[, as.character(tempcall$tstart)], extend=TRUE)$surv
        data2$weights[data2[, id] == i] <- 1/denom_surv
        
        
      }
      rows.1 <- which(data[, arm] == levArms[1])
      data$weights[rows.1] <- data1$weights
      rows.2 <- which(data[, arm] == levArms[2])
      data$weights[rows.2] <- data2$weights
    } else{
      # if denominator only
      # arm 1
      data1 <- data[data[, arm] == levArms[1],]
      fit.cox.denom.1 <- coxph(eval(parse(text = paste("Surv(",
                                                       deparse(tempcall$tstart), ", ",
                                                       deparse(tempcall$tstop), ", ",
                                                       deparse(tempcall$cens), ") ~  ",
                                                       denom, sep = ""))),
                               data = data1)
      
      idpat1 <- unique(data1[, id])
      for (i in idpat1) {
        datai <- data1[data1[, id] == i, ]
        km.num<- survfit(Surv(data1$tstart, data1$tstop, data1$cens)~1, data=data1, type=type)
        
        km.denom <- survfit(fit.cox.denom.1, newdata = datai, individual = TRUE, type = type)
        denom_surv<- summary(km.denom,times = datai[, as.character(tempcall$tstart)],extend=TRUE)$surv
        data1$weights[data1[, id] == i] <-  1 / denom_surv
        #  data1$weights[data1[, id] == i] <-
        #   1 / summary(km.denom, times = datai[, as.character(tempcall$tstart)])$surv
      }
      
      # arm 2
      data2 <- data[data[, arm] == levArms[2],]
      fit.cox.denom.2 <- coxph(eval(parse(text = paste("Surv(",
                                                       deparse(tempcall$tstart), ", ",
                                                       deparse(tempcall$tstop), ", ",
                                                       deparse(tempcall$cens), ") ~  ",
                                                       denom, sep = ""))),
                               data = data2)
      
      idpat2 <- unique(data2[, id])
      for (i in idpat2) {
        datai <- data2[data2[, id] == i, ]
        km.num<- survfit(Surv(data2$tstart, data2$tstop, data2$cens)~1, data=data2, type=type)
        km.denom <- survfit(fit.cox.denom.2, newdata = datai, individual = TRUE, type = type)
        denom_surv<- summary(km.denom, times = datai[, as.character(tempcall$tstart)],extend=TRUE)$surv
        data2$weights[data2[, id] == i] <- 1 / denom_surv
        #     data2$weights[data2[, id] == i] <-
        #      1 / summary(km.denom, times = datai[, as.character(tempcall$tstart)])$surv
      }
      # rows with 1st level of arm
      rows.1 <- which(data[, arm] == levArms[1])
      data$weights[rows.1] <- data1$weights
      # rows with 2nd level of arm
      rows.2 <- which(data[, arm] == levArms[2])
      data$weights[rows.2] <- data2$weights
    }
    
  }
  
  # Truncated weights (optional)
  if (!(is.null(tempcall$trunc))) {
    data$weights.trunc <- data$weights
    data$weights.trunc[data$weights <= quantile(data$weights,
                                                0 + trunc)] <-
      quantile(data$weights, 0 + trunc)
    data$weights.trunc[data$weights > quantile(data$weights,
                                               1 - trunc)] <-
      quantile(data$weights, 1 - trunc)
  }
  return(data)
}

#### true_phi is computed from coxme(Surv(Y_PE, event_E) ~ treated + (1 | pair.index), data = stacked) in true_value_function

## treatment ratio
# TR: 20%; baseline hazard: 0.3; true value: -0.5163079
# TR: 30%; baseline hazard: 0.45; true value: -0.4967315
# TR: 40%; baseline hazard: 0.7; true value: -0.4900718
# TR: 50%; baseline hazard: 1.1; true value: -0.4937465
# TR: 60%; baseline hazard: 1.55; true value: -0.5002375

estimated_purecontrol_value <- function(true_phi, sample_size, simulation, lambda_z) {
  library(survival)
  library(optmatch)
  library(MASS)
  library(coxme) #used for fitting the frailty model
  require(foreach)
  require(doMC)
  library(dplyr)
  library(tableone)
  library(ipcwswitch)
  
  # sample size
  n <- sample_size
  nsim <- simulation
  
  #treated size after matching
  true.treat.after_matched <- rep(NA, 1000)
  
  #treated ratio before matching
  treat.ratio_before_matched <- rep(NA, 1000)
  
  # treated moved to control
  treated.in.control <- rep(NA, 1000)
  
  # estimated results
  est1 <- rep(NA, 1000)
  est2 <- rep(NA, 1000)
  
  # standard error
  se1 <- rep(NA, 1000)
  se2 <- rep(NA, 1000)
  
  # if in 95% CI, coverage
  ci1 <- rep(NA, 1000)
  ci2 <- rep(NA, 1000)
  
  # mse
  mse_ss1 <- rep(NA, 1000)
  mse_ss2 <- rep(NA, 1000)
  
  # bias
  bias_ss1 <- rep(NA, 1000)
  bias_ss2 <- rep(NA, 1000)
  
  sim_id=c()
  
  # set.seed setting
  w=1
  a=5
  
  
  # start with an empty matched data
  #DtMatched <- data.frame()
  
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
    #lambda_z <- 0.7 #decrease this will help increase T_z
    
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
    dat1 <- new.dat[,c("id", "X1", "X2", "X3_0", "X3_1", "T_e", "T_z", "CensorT", "Y_z", "Y_e", "observed.treat", "observed.event")]
    
    # define a new competing variable
    # censor: 0; T_z:1; T_e:2
    dat1$cr <- ifelse(dat1$CensorT<pmin(dat1$T_z,dat1$T_e),0,ifelse(dat1$T_z<pmin(dat1$T_e,dat1$CensorT),1,2)) 
    
    # initialize
    newcgd <- tmerge(dat1, dat1, id = id, tstart = 0, tstop = Y_z)
    
    # Add time interval split at 0.5
    newcgd$cutpt <- t1
    newcgd <- tmerge(newcgd, newcgd, id = id, splitpt = tdc(cutpt)) 
    
    # Add observed time indicator
    newcgd <- tmerge(newcgd, dat1, id = id, competing = event(Y_z, cr))
    
    # arrange X3 to each time interval
    newcgd$X3 <- ifelse(newcgd$tstop<=t1,newcgd$X3_0,newcgd$X3_1)
    
    
    # keep cleaned data
    dat_cut <- newcgd[,c("id", "X1", "X2", "X3", "tstart", "tstop", "competing")]
    
    # build cox model
    # 0.3,0.3,0.2
    cox1 <- coxph(Surv(tstart, tstop, competing==1) ~ X1 + X2 + X3 + cluster(id), data = dat_cut, control = coxph.control(timefix = FALSE))
    summary(cox1)
    # get the estimated coefficients
    beta_hat<-matrix(cox1$coefficients, ncol=1)
    
    
    ##############
    ## matching ##
    ##############
    
    # take the time to treatment out from the data
    times <- dat_cut$tstop[which(dat_cut$competing==1)]
    times <- sort(unique(times))
    
    DtMatched<- NULL
    j<- 0 
    res1<- NULL
    dat_res<- new.dat
    dat_long<- dat_cut
    
    # calculate propensity score
    dat_long$ps<- exp((as.matrix(dat_long[,c('X1', 'X2', 'X3')]))%*%beta_hat) # decided by time-varying X3
    
    for (i in 1:length(times)) {
      time<- times[i]
      time.chara<- as.character(time)
      time_cum<- c(0, times[1:i])
      ############## find the risk set and event set
      risk.set<- dat_res[dat_res$observed.treat==0&dat_res$Y_z>time,] # who doesn't get treatment at the matching time and afterward
      event.set<- dat_res[as.character(dat_res$Y_z)==time.chara,] # who has the treatment at the matching time "times"
      # there will be only one event at each matching point
      
      if (nrow(risk.set) > 0 & nrow(event.set) > 0) { 
        
        # get the data with long format
        risk.set_long <- dat_long[dat_long$id%in%risk.set$id,]
        event.set_long<- dat_long[dat_long$id%in%event.set$id,]
        dist_sq<- 0
        
        if(time<=0.5){ # it means that event set only has one row (>0 and <t1)
          tmp_dat.pt<- data.frame(id=c(event.set_long$id, risk.set_long$id), time=rep(time,nrow(event.set_long)+nrow(risk.set_long)),
                                  zStatus=c(1, rep(0, nrow(risk.set_long))), hz=c(event.set_long$ps, risk.set_long$ps))
          dist.pt<- match_on(zStatus~hz, data=tmp_dat.pt, method='euclidean') # Computes pairwise Euclidean distances between treated (eStatus == 1) and control (eStatus == 0) units based on the hz variable
          dist_sq.pt<- dist.pt^2
          # Accumulates weighted squared distance, here we set weight to be 1
          dist_sq<- dist_sq+(time_cum[2]-time_cum[1])*dist_sq.pt
          #dist_sq<- dist_sq+dist_sq.pt
        } 
        
        if(time>0.5){# do a for loop for the event set have two time intervals
          
          for(k in 1:nrow(event.set_long)){
            
            tmp_dat.pt<- data.frame(id=c(event.set$id, risk.set$id), 
                                    time=rep(time,length(event.set$id)+length(risk.set$id)),
                                    zStatus=c(1, rep(0,length(risk.set$id))), 
                                    hz=c(event.set_long$ps[k], 
                                         risk.set_long$ps[seq(k, nrow(risk.set_long), 
                                                              by=nrow(event.set_long))]))
            
            rownames(tmp_dat.pt)<- tmp_dat.pt$id
            # calculate the distance metric between control and treated subjects
            dist.pt<- match_on(zStatus~hz, data=tmp_dat.pt, method='euclidean')
            dist_sq.pt<- dist.pt^2
            
            dist_sq<- dist_sq + (event.set_long$tstop[k] - event.set_long$tstart[k])*dist_sq.pt # calculate summation of distance based on nrow(event.set)
          }
        }
        
        # choose the pair with closest distance
        mahal.match <- pairmatch(dist_sq, data = tmp_dat.pt, controls = 1)
        DTwithGrp <- cbind(tmp_dat.pt, matches = mahal.match)
        dtMatched <- DTwithGrp[!is.na(DTwithGrp$matches),] # find the matched pairs
        dat_res <- dat_res[!(dat_res$id %in% dtMatched$id), ] # remove the id in the matched pairs for the original data, and then keep finding new matching sets.
        DtMatched <- rbind(DtMatched, dtMatched)
      } else {
        j <- j + 1
        res1 <- c(res1, time)
        next
      }
      print(i)
    }
    
    DtMatched$id<- as.factor(DtMatched$id)
    new.dat$id<- as.factor(new.dat$id)
    
    # final data after matching
    DtMatched<- plyr::join(DtMatched, new.dat, by='id', match='all')
    
    # assign index for each pair
    DtMatched$pair.index<- rep((1:(dim(DtMatched)[1]/2)), each=2)
    
    ######### note: Y_e=min(T_e,C); Y_z=min(T_z,Y_e)
    # identify post-treatment time
    
    DtMatched$Y_pe <- -1
    
    ### when zstatus=1, the subject is treated 
    # T_z is the matching time and also treated time (time=T_z)
    # the post-treatment event time is L=Y_e-time
    DtMatched$Y_pe[DtMatched$zStatus==1] <- DtMatched$Y_e[DtMatched$zStatus==1]-DtMatched$time[DtMatched$zStatus==1]
    
    ### when zstatus=0, the subject is control
    # T_z is the matching time, but not the real treated time
    # there are three scenarios
    
    ### scenario 1:
    # observed.treat=1, and the subject is treated after the event time
    # the post-treatment event time is L=Y_e-time
    ### scenario 2:
    # observed.treat=1, and the subject is treated before the event time
    # the post-treatment events are censored by the treatment initiation
    # the censored event time is L=T_z-time
    DtMatched$Y_pe[DtMatched$zStatus==0&DtMatched$observed.treat==1] <- ifelse(DtMatched$T_z[DtMatched$zStatus==0&DtMatched$observed.treat==1]>DtMatched$Y_e[DtMatched$zStatus==0&DtMatched$observed.treat==1],
                                                                               DtMatched$Y_e[DtMatched$zStatus==0&DtMatched$observed.treat==1]-DtMatched$time[DtMatched$zStatus==0&DtMatched$observed.treat==1],
                                                                               DtMatched$T_z[DtMatched$zStatus==0&DtMatched$observed.treat==1]-DtMatched$time[DtMatched$zStatus==0&DtMatched$observed.treat==1])
    
    ### scenario 3:
    # observed.treat=0, the subject is never treated
    # the post-treatment event time is L=Y_e-time
    DtMatched$Y_pe[DtMatched$zStatus==0&DtMatched$observed.treat==0] <-DtMatched$Y_e[DtMatched$zStatus==0&DtMatched$observed.treat==0]-DtMatched$time[DtMatched$zStatus==0&DtMatched$observed.treat==0]
    
    
    # identify the event indicator
    DtMatched$delta<- -1
    
    ### when zstatus=1, we just need to identify censoring or event occurs first
    DtMatched$delta[DtMatched$zStatus==1] <- ifelse(DtMatched$T_e[DtMatched$zStatus==1]<DtMatched$CensorT[DtMatched$zStatus==1],1, 0) #treatment time before censored time
    
    ### when zstatus=0 and observed.treat=1, the later treatment can happen after or before the event
    # event<censor<T_z
    # censor<event<T_z
    # censor<T_z<event
    # event<T_z<censor
    # T_z<censor/event
    DtMatched$delta[DtMatched$zStatus==0&DtMatched$observed.treat==1]<- ifelse(((DtMatched$Y_e[DtMatched$zStatus==0&DtMatched$observed.treat==1]<
                                                                                   DtMatched$T_z[DtMatched$zStatus==0&DtMatched$observed.treat==1]) & (DtMatched$T_e[DtMatched$zStatus==0&DtMatched$observed.treat==1]<
                                                                                                                                                         DtMatched$CensorT[DtMatched$zStatus==0&DtMatched$observed.treat==1])),1, 
                                                                               ifelse(((DtMatched$Y_e[DtMatched$zStatus==0&DtMatched$observed.treat==1]<
                                                                                          DtMatched$T_z[DtMatched$zStatus==0&DtMatched$observed.treat==1]) & (DtMatched$T_e[DtMatched$zStatus==0&DtMatched$observed.treat==1]>
                                                                                                                                                                DtMatched$CensorT[DtMatched$zStatus==0&DtMatched$observed.treat==1])),0,
                                                                                      ifelse(((DtMatched$CensorT[DtMatched$zStatus==0&DtMatched$observed.treat==1]<DtMatched$T_z[DtMatched$zStatus==0&DtMatched$observed.treat==1]) & (DtMatched$T_e[DtMatched$zStatus==0&DtMatched$observed.treat==1] > DtMatched$T_z[DtMatched$zStatus==0&DtMatched$observed.treat==1])),0,
                                                                                             ifelse(((DtMatched$T_e[DtMatched$zStatus==0&DtMatched$observed.treat==1] < DtMatched$T_z[DtMatched$zStatus==0&DtMatched$observed.treat==1]) & (DtMatched$CensorT[DtMatched$zStatus==0&DtMatched$observed.treat==1] > DtMatched$T_z[DtMatched$zStatus==0&DtMatched$observed.treat==1])),1,0))))
    
    ### when zstatus=0 and observed.treat=0, there is no treatment observed, so we just need to identify censoring or event occurs first  
    DtMatched$delta[DtMatched$zStatus==0&DtMatched$observed.treat==0]<- ifelse(DtMatched$T_e[DtMatched$zStatus==0&DtMatched$observed.treat==0]<
                                                                                 DtMatched$CensorT[DtMatched$zStatus==0&DtMatched$observed.treat==0],1, 0)
    
    ## redefine tstop 
    DtMatched$new.Y_e <- NA
    DtMatched$new.Y_e[DtMatched$zStatus==1] <- DtMatched$Y_e[DtMatched$zStatus==1]
    DtMatched$new.Y_e[DtMatched$zStatus==0&DtMatched$observed.treat==1] <- ifelse(DtMatched$T_z[DtMatched$zStatus==0&DtMatched$observed.treat==1]>DtMatched$Y_e[DtMatched$zStatus==0&DtMatched$observed.treat==1],
                                                                                  DtMatched$Y_e[DtMatched$zStatus==0&DtMatched$observed.treat==1],
                                                                                  DtMatched$T_z[DtMatched$zStatus==0&DtMatched$observed.treat==1])
    DtMatched$new.Y_e[DtMatched$zStatus==0&DtMatched$observed.treat==0] <-DtMatched$Y_e[DtMatched$zStatus==0&DtMatched$observed.treat==0]
    
    
    
    
    
    # --- FRAILTY COX MODEL AFTER MATCHING ---
    test_frailty <- coxme(Surv(Y_pe, delta) ~ zStatus+(1|pair.index), data=DtMatched)
    
    # estimated log hazard ratio for treatment vs control post-initiation.
    est1[sim] <- test_frailty$coefficients
    
    se1[sim] <- sqrt(diag(vcov(test_frailty)))
    
    z_crit <- qnorm(0.975)  # = 1.96
    
    # Manual CI on log(HR) scale
    lower_log1 <- est1[sim] - z_crit * se1[sim]
    upper_log1 <- est1[sim] + z_crit * se1[sim]
    
    ci1[sim]<- ifelse(lower_log1<true_phi&upper_log1>true_phi, 1, 0)
    
    bias_ss1[sim] <- est1[sim]-true_phi
    mse_ss1[sim] <- bias_ss1[sim]^2
    
    # --- COX MODEL AFTER MATCHING ---
    test_normal <- coxph(Surv(Y_pe, delta) ~ zStatus, data=DtMatched)
    
    # estimated log hazard ratio for treatment vs control post-initiation.
    est2[sim] <- test_normal$coefficients
    
    se2[sim] <- sqrt(diag(vcov(test_normal)))
    
    z_crit <- qnorm(0.975)  # = 1.96
    
    # Manual CI on log(HR) scale
    lower_log2 <- est2[sim] - z_crit * se2[sim]
    upper_log2 <- est2[sim] + z_crit * se2[sim]
    
    ci2[sim]<- ifelse(lower_log2<true_phi&upper_log4>true_phi, 1, 0)
    
    bias_ss2[sim] <- est2[sim]-true_phi
    mse_ss2[sim] <- bias_ss2[sim]^2
    
    
    sim_id<- c(sim_id, sim) #collect simulation times
    
    # check treated and control groups after matching
    true.treat.after_matched[sim] <- length(which(DtMatched$zStatus==1))
    
    # check treated and control groups before matching
    treat.ratio_before_matched[sim] <- length(which(new.dat$observed.treat==1))
    
    # check who are treated before matching are in control after matching
    treated.control.after <- DtMatched[which(DtMatched$zStatus==0),] # control subjects after matching
    treated.treated.after <- DtMatched[which(DtMatched$zStatus==1),] # treated subjects after matching
    treated.treated.before <- new.dat[which(new.dat$observed.treat==1),] # treated subjects before matching
    treated.in.control[sim] <- length(which(treated.treated.before$id %in% treated.control.after$id))
    
    sim <- sim+1
    
  }
  
  return(list(frailty_est=est1,frailty_se=se1,frailty_ci=ci1,frailty_mse=mse_ss1,
              cox_est=est2,cox_se=se2,cox_ci=ci2,cox_mse=mse_ss2,
              true.treated_after=true.treat.after_matched,
              treated_before=treat.ratio_before_matched, treated.to.control=treated.in.control))
}


######################## example
true_phi <- -0.4900718

estimated <- estimated_purecontrol_value(true_phi3,true_phi4,500,1000,0.7)
# here, true_phi3 is the true value, 500 is the sample size, 1000 is the simulation times, 0.7 is value of lambda_z which means treatment ratio is 40%

## ASE
frailty_ase <- mean(na.omit(estimated$frailty_se))
ipcw_ase <- mean(na.omit(estimated$ipcw_se))

## bias
frailty_bias <- mean(na.omit(estimated$frailty_est))-true_phi
ipcw_bias <- mean(na.omit(estimated$ipcw_est))-true_phi

## MAB (mean absolute bias)
frailty_absolutebias <- mean(abs(na.omit(estimated$frailty_est)-na.omit(true_phi)))
ipcw_absolutebias <- mean(abs(na.omit(estimated$ipcw_est)-na.omit(true_phi)))

## coverage rate
frailty_cover <- length(which(estimated$frailty_ci==1))/length(na.omit(estimated$frailty_ci))
ipcw_cover <- length(which(estimated$ipcw_ci==1))/length(na.omit(estimated$ipcw_ci))

## MCSD
frailty_mcsd <- sd(na.omit(estimated$frailty_est))
ipcw_mcsd <- sd(na.omit(estimated$ipcw_est))

## MSE
frailty_mse <- mean((na.omit(estimated$frailty_est)-true_phi)^2)
ipcw_mse <- mean((na.omit(estimated$ipcw_est)-true_phi)^2)




