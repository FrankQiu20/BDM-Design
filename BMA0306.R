get.oc.bma.oc<-function(pT.true, pE.true, rho, targetT=0.35,targetE=0.25,ncohort=20,cohortsize=3,
                        startdose=1,cutoff.eli.T=0.9, alpha, deltaE = 0.9,Cp=0.43,
                        cutoff.eli.E=0.90,ntrial=1000){
  library(rjags)
  library(coda)
  library(runjags)
  
  library(rjags)
  
  BMA<-function(N,YE){
    H<-length(N[N!=0])
    
    marginal_likelihood_estimates<-rep(0,H)
    p_samples_list <- list()
    for (h in 1:H) {
      if(h==1){
        n<-sum(N[h:H])
        y<-sum(YE[h:H])
        J=length(n)
        # Data list for JAGS
        dataList <- list(
          y = y,
          n = n,
          J = J
        )
        # JAGS model string
        modelString <- "
    model {
    for (r in 1:J) {
        # Priors for beta
        beta[r] ~ dbeta(1, 1)
    }
    for (j in 1:J) {
        # likelihood
        y[j] ~ dbin(p[j], n[j])
        # link function for the probability using a product of beta terms
        p[j] <- 1 - prod(1 - beta[1:j])
    }
    }
    "
        initValues <- function() {
          list(beta = runif(J))
        }
        
        # Parameters to monitor
        params <- c("p")
        # Setting up the model
        jagsModel <- jags.model(textConnection(modelString), data = dataList, inits = initValues, n.chains = 1,quiet=T)
        
        # Burn-in
        update(jagsModel, n.iter = 1000,progress.bar="none")  # adjust the burn-in period as needed
        
        # MCMC
        samples <- coda.samples(jagsModel, variable.names = params, n.iter = 10000,progress.bar="none")
        p_samples <- as.matrix(samples)
        
        ####
        likelihood <- function(params) {
          # This function needs to calculate the log-likelihood for a given set of beta parameters
          # params is a vector of the beta parameters
          # You will need to adapt this function based on your specific model
          ll <- 1
          for (i in 1:length(y)) {
            #p <- 1 - prod(1 - params[1:i])
            ll <- ll * dbinom(y[i], size = n[i], prob = params[i])/(choose(n[i],y[i]))
          }
          return(ll)
        }
        ####
        
        likelihoods <- apply(p_samples, 1, function(params) likelihood(params))
        
        # Compute the harmonic mean of the likelihoods
        harmonic_mean_likelihood <- (1/length(likelihoods)) * sum(1 / likelihoods)
        
        # Convert to the estimate of the marginal likelihood
        marginal_likelihood_estimate <- 1 / harmonic_mean_likelihood
        marginal_likelihood_estimates[h]<-marginal_likelihood_estimate
        replicate_indicator <- matrix(1, ncol = H)
        p_samples_list[[h]] <- p_samples %*% replicate_indicator
      }else{
        n<-c(N[1:(h-1)],sum(N[h:H]))
        y<-c(YE[1:(h-1)],sum(YE[h:H]))
        J=length(n)
        # Data list for JAGS
        dataList <- list(
          y = y,
          n = n,
          J = J
        )
        modelString <- "
      model {
    for (r in 1:J) {
        # Priors for beta
        beta[r] ~ dbeta(1, 1)
    }
    for (j in 1:J) {
        # likelihood
        y[j] ~ dbin(p[j], n[j])
        # link function for the probability using a product of beta terms
        p[j] <- 1 - prod(1 - beta[1:j])
    }
    }
      "
        initValues <- function() {
          list(beta = runif(J))
        }
        
        # Parameters to monitor
        params <- c("p")
        # Setting up the model
        jagsModel <- jags.model(textConnection(modelString), data = dataList, inits = initValues, n.chains = 1,quiet=T)
        
        # Burn-in
        update(jagsModel, n.iter = 1000,progress.bar="none")  # adjust the burn-in period as needed
        
        # MCMC
        samples <- coda.samples(jagsModel, variable.names = params, n.iter = 10000,progress.bar="none")
        p_samples <- as.matrix(samples[,c(1:h)])
        
        replicate_indicator <- matrix(1, ncol = length(h:H))
        p_samples_list[[h]] <- cbind(p_samples[,1:h-1],p_samples[,h:ncol(p_samples)] %*% replicate_indicator)
        ####
        likelihood <- function(params) {
          # This function needs to calculate the log-likelihood for a given set of beta parameters
          # params is a vector of the beta parameters
          # You will need to adapt this function based on your specific model
          ll <- 1
          for (i in 1:length(y)) {
            #p <- 1 - prod(1 - params[1:i])
            ll <- ll * dbinom(y[i], size = n[i], prob = params[i])/(choose(n[i],y[i]))
          }
          return(ll)
        }
        ####
        
        likelihoods <- apply(p_samples, 1, function(params) likelihood(params))
        harmonic_mean_likelihood <- (1/length(likelihoods)) * sum(1 / likelihoods)
        marginal_likelihood_estimate <- 1 / harmonic_mean_likelihood
        marginal_likelihood_estimates[h]<-marginal_likelihood_estimate
        
      }
      
    }
    return(list(marginal_likelihood_estimates = marginal_likelihood_estimates, 
                p_samples = p_samples_list))
  }
  
  
  Tox_bayes<-function(N,YT){
    H<-length(N[N!=0])
    n<-N[1:H]
    y<-YT[1:H]
    J=length(n)
    # Data list for JAGS
    dataList <- list(
      y = y,
      n = n,
      J = J
    )
    # JAGS model string
    modelString <- "
    model {
    for (r in 1:J) {
        # Priors for beta
        beta[r] ~ dbeta(1, 1)
    }
    for (j in 1:J) {
        # likelihood
        y[j] ~ dbin(p[j], n[j])
        # link function for the probability using a product of beta terms
        p[j] <- 1 - prod(1 - beta[1:j])
    }
    }
    "
    initValues <- function() {
      list(beta = runif(J))
    }
    # Parameters to monitor
    params <- c("p")
    # Setting up the model
    jagsModel <- jags.model(textConnection(modelString), data = dataList, inits = initValues, n.chains = 1,quiet=T)
    
    # Burn-in
    update(jagsModel, n.iter = 1000,progress.bar="none")  # adjust the burn-in period as needed
    
    # MCMC
    samples <- coda.samples(jagsModel, variable.names = params, n.iter = 10000,progress.bar="none")
    p_samples <- as.matrix(samples)
    return(tox_post = p_samples)
  }
  
  prob_T<-function(x){
    mean(x>targetT)
  }
  
  joint.p=function(ttox,teff,c){
    f1<-function(x,bn.m1,bn.m2,rho){
      ff<-dnorm(x,bn.m1,1)*(1-pnorm(0,bn.m2+rho*(x-bn.m1),sqrt(1-rho^2)))
      return(ff)
    }
    ## ttox: marginal toxicity probability
    ## teff: marginal efficacy probability
    ## c: association parameter between tox and eff, c>0 indicates a positive correlation
    ndose=length(ttox) ## dose level
    out=matrix(rep(0,4*ndose),nrow=4) ## joint tox-eff probability matrix; row is dose level, column=(tox=0,eff=0; tox=0,eff=1; tox=1,eff=0; tox=1,eff=1)
    for (i in 1: ndose){
      out[4,i]=integrate(f1,bn.m1=qnorm(ttox[i]),bn.m2=qnorm(teff[i]),rho=c,lower=0,upper=Inf)$value
      out[2,i]=teff[i]-out[4,i]
      out[3,i]=ttox[i]-out[4,i]
      out[1,i]=1-out[4,i]-out[2,i]-out[3,i]
    }
    return(out)
  }
  
  ### Get outcome for each cohort
  pts.outcome <- function(d,cohortsize,jointp){
    n = rep(1,cohortsize)
    yT = rep(0,cohortsize)
    yE = rep(0,cohortsize)
    for (pts in 1:cohortsize) {
      res=rmultinom(1,1,jointp[,d])
      yT[pts] = res[3]+res[4]
      yE[pts] = res[2]+res[4]
    }
    return(list(
      yT=sum(yT),
      yE=sum(yE)))
  }
  model_post_E<-function(alpha,high,result){
    prior<-c(1:high)^alpha/sum(c(1:high)^alpha)
    model_post<-result$marginal_likelihood_estimates*prior/sum(result$marginal_likelihood_estimates*prior)
    return(model_post)
  }
  
  prob_E<-function(T_H,list_p,phiE){
    matrix_prob<-matrix(0,nrow = T_H,ncol = T_H)
    for (i in 1:T_H) {
      for (j in 1:T_H) {
        matrix_prob[i,j]<-mean(list_p[[i]][,j]<phiE)
      }
    }
    return(matrix_prob)
  }
  
  get_plateau_prob<-function(T_H,list_p,delta){
    matrix_prob<-matrix(0,nrow = T_H,ncol = T_H)
    for (i in 1:T_H) {
      for (j in 1:T_H) {
        matrix_prob[i,j]<-mean(list_p[[i]][,j]>delta*list_p[[i]][,T_H])
      }
    }
    return(matrix_prob)
  }
  
  get_eff_est<-function(T_H,list_p){
    matrix_eff<-matrix(0,nrow = T_H,ncol = T_H)
    for (i in 1:T_H) {
      for (j in 1:T_H) {
        matrix_eff[i,j]<-mean(list_p[[i]][,j])
      }
    }
    return(matrix_eff)
  }
  
  joinp<-joint.p(pT.true,pE.true,c=rho)
  dselect = rep(0, ntrial)
  ndose = length(pE.true)
  N = matrix(rep(0, ndose * ntrial), ncol = ndose)
  YTOX = matrix(rep(0, ndose * ntrial), ncol = ndose)
  YEFF = matrix(rep(0, ndose * ntrial), ncol = ndose)
  
  for (trial in 1:ntrial) {
    d=startdose
    cohortsize = cohortsize
    ytox = rep(0,ndose)
    yeff = rep(0,ndose)
    n = rep(0,ndose)
    for (i in 1:ncohort) {
      if(d==0){
        break
      }
      outcome = pts.outcome(d=d,cohortsize = cohortsize,joinp)
      n[d] = n[d] + cohortsize
      ytox[d] = ytox[d] + outcome$yT
      yeff[d] = yeff[d] + outcome$yE
      h_tried = length(n[n!=0])
      
      bma_E_sample<-BMA(N=n,YE=yeff)
      T_sample<-Tox_bayes(N=n,YT=ytox)
      Adm_tox<-which(apply(T_sample, 2, prob_T)<cutoff.eli.T)
      
      
      model_post_Eff<-model_post_E(alpha = alpha,high = h_tried,result = bma_E_sample)
      Adm_eff<-which(model_post_Eff %*% prob_E(h_tried,bma_E_sample$p_samples,targetE)<cutoff.eli.E)
      
      #plateau_prob = model_post_Eff %*% get_plateau_prob(T_H = h_tried,list_p = bma_E_sample$p_samples,delta = deltaE)
      eff_est = model_post_Eff %*% get_eff_est(h_tried,bma_E_sample$p_samples)
      
      A = intersect(Adm_eff,Adm_tox)
      if (length(A)==0){
        d = 0;
        break
      } else {
        if (h_tried %in% Adm_tox & h_tried < ndose){
          d = d + 1
        }else{
          B = c(1:max(min(d+1,ndose),h_tried))
          adm_next = intersect(A,B)
          platuau = adm_next[which(eff_est[,adm_next]>deltaE*eff_est[,h_tried])]
          if(length(platuau)==0){
            #d = max(adm_next)
            lobd = max(adm_next)
            if(lobd>d){
              d = d + 1
            }else if(lobd<d){
              d = d - 1
            }else{
              d = d
            }
          }else{
            lobd = min(platuau)
            if(lobd>d){
              d = d + 1
            }else if(lobd<d){
              d = d - 1
            }else{
              d = d
            }
          }
        }
      }
    }
    
    YTOX[trial, ] = ytox
    YEFF[trial, ] = yeff
    N[trial,]=n;
    
    if(d == 0){
      dselect[trial] = 0 
    }else{
      h_tried = length(n[n!=0])
      bma_E_sample<-BMA(N=n,YE=yeff)
      T_sample<-Tox_bayes(N=n,YT=ytox)
      Adm_tox<-which(apply(T_sample, 2, prob_T)<cutoff.eli.T)
      
      model_post_Eff<-model_post_E(alpha = alpha,high = h_tried,result = bma_E_sample)
      Adm_eff<-which(model_post_Eff %*% prob_E(h_tried,bma_E_sample$p_samples,targetE)<cutoff.eli.E)
      
      plateau_prob = model_post_Eff %*% get_plateau_prob(T_H = h_tried,list_p = bma_E_sample$p_samples,delta = deltaE)
      
      A = intersect(Adm_eff,Adm_tox)
      
      if(length(A)==0){
        obd = 0
      } else{
        MTD = which.min(abs(apply(T_sample, 2, mean)-targetT))
        A = A[A<=MTD]
       # eff_est_final = model_post_Eff %*% get_eff_est(h_tried,bma_E_sample$p_samples)
        
      #  platuau = A[which(eff_est_final[A]>deltaE*eff_est_final[h_tried])]
        platuau = A[which(plateau_prob[A]>Cp)]
        
        if(length(platuau)==0){
          obd = max(A)
        }else{
          obd = min(platuau)
        }
      }
      
      
      dselect[trial] =obd
    }
  }
  
  selpercent = rep(0, ndose + 1)
  patpercent = matrix(rep(0, ntrial * ndose), ncol = ntrial, nrow = ndose)
  efficacy = rep(0, ntrial)
  toxicity = rep(0, ntrial)
  f <- function(x) {
    x[i] / sum(x)
  }
  ## Summarize results
  for (i in 0:ndose) {
    selpercent[(i + 1)] = sum(dselect == i) / ntrial * 100
  }
  print("selection probablity")
  cat(formatC(selpercent, digits = 1, format = "f"), sep = " ", "\n")
  for (i in 1:ndose) {
    patpercent[i, ] = apply(N, 1, f)
  }
  print("average percent of patients")
  cat(formatC(
    apply(patpercent, 1, mean) * 100,
    digits = 1,
    format = "f"
  ),
  sep = " ", "\n")
  print("average number of patients")
  cat(formatC(c(apply(N, 2, mean), sum(apply(
    N, 2, mean
  ))), digits = 1, format = "f"),
  sep = " ", "\n")
  print("average number of patients response to efficacy")
  cat(formatC(c(apply(YEFF, 2, mean), sum(
    apply(YEFF, 2, mean)
  )), digits = 1, format = "f"),
  sep = " ", "\n")
  for (i in 1:ntrial) {
    efficacy[i] = sum(YEFF[i, ]) / sum(N[i, ])
  }
  print("average percent of efficacy")
  cat(formatC(mean(efficacy) * 100, digits = 1, format = "f"),
      sep = " ", "\n")
  print("average number of patients response to toxicity")
  cat(formatC(c(apply(YTOX, 2, mean), sum(
    apply(YTOX, 2, mean)
  )), digits = 1, format = "f"),
  sep = " ", "\n")
  for (i in 1:ntrial) {
    toxicity[i] = sum(YTOX[i, ]) / sum(N[i, ])
  }
  print("average percent of toxicity")
  cat(formatC(mean(toxicity) * 100, digits = 1, format = "f"),
      sep = " ", "\n")
}


get.oc.bma.oc(c(0.01,0.02,0.05,0.07,0.1,0.14),c(0.1,0.2,0.4,0.55,0.78,0.8),0,0.3,0.2,alpha = 1.2,ntrial = 100,Cp=0.5)
