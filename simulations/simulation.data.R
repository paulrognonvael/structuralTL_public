##############################################################################
##############################################################################
############################# Simulate data ##################################
##############################################################################
##############################################################################
path = "~/structuralTL/simulations/sim.data"

data.sim = function(n.t, p.t, s.t, sd.t, n.s, p.s, sd.s,trans.inf){
  #Input:
  # - n.t: target sample size
  # - p.t: target dimension
  # - s.t: target number of nonzero coefficients
  # - sd.t: error standard deviation in target linear model
  # - n.s: a vector of source datasets sample size
  # - p.s: a vector of source datasets dimension
  # - sd.s: a vector of source datasets error standard deviation (linear model)
  # - trans.inf: a vector of source datasets informativeness (each one of 
  #   'v.informative', 'm.informative', 'adversarial')
  #
  # Output: a list where
  # - target.data.diffbetmin: target dataset where the small signal is different across blocks
  # - target.data.samebetmin: target dataset where the small signal is the same across blocks
  # - source.data.covshift: a list with source datasets
  # - beta.t: target data generating coefficients
  # - beta.s: source data generating coefficients
  
  
  # Set up objects to store output
  k = length(n.s)
  source.data.covshift = list()
  source.data.nocovshift = list()
  store.beta.s = data.frame(matrix(NA, ncol = p.s, nrow = k))
  
  for(i in 1:k){
    n.s.i <- n.s[i]
    p.s.i <- p.s
    trans.inf.i <- trans.inf[i]
    
    # Define source support
    s.s.i <- s.t
    if(trans.inf.i=='v.informative'){
      S.sintS.t = 10 #round(0.9*s.t)
      S.s.i <- c(rep(1,S.sintS.t), rep(0,s.t-S.sintS.t), 
                 rep(1,s.s.i-S.sintS.t), rep(0,p.s.i-s.s.i-(s.t-S.sintS.t)))
    }
    
    if(trans.inf.i=='m.informative'){
      S.sintS.t = 8 #round(2*s.t/3)
      S.s.i <- c(rep(0,s.t-S.sintS.t), rep(1,S.sintS.t), 
               rep(1,s.s.i-S.sintS.t), rep(0,p.s.i-s.s.i-(s.t-S.sintS.t)))
    }
    
    if(trans.inf.i=='adversarial'){
      S.s.i <- c(rep(0,s.t), rep(1,s.s.i), rep(0,p.s.i-(s.s.i+s.t)))
    }
    
    ## Simulate source coefficients
    beta.s.i <- rep(0,p.s)
    beta.s.i[as.logical(S.s.i)] <- sample(c(-1,1),size=s.s.i,replace = TRUE)*runif(s.s.i, 1/3, 1.5)
    
    ## save source parameters
    store.beta.s[i,] = beta.s.i
    
    ## Generate source predictors and response
    ### No covariate shift
    # cov.s <- cov.t
    # X.nocov.s <- mvtnorm::rmvnorm(n.s.i, sigma = cov.s)
    # X.nocov.s <- X.nocov.s %*% diag(1/sqrt(diag(t(X.nocov.s) %*% X.nocov.s)/n.s.i))
    
    ### With covariate shift
    cov.s <- outer(1:p.t, 1:p.t, function(x,y){0.1^(abs(x-y))})
    X.cov.s <- mvtnorm::rmvnorm(n.s.i, sigma = cov.s)
    X.cov.s <- X.cov.s %*% diag(1/sqrt(diag(t(X.cov.s) %*% X.cov.s)/n.s.i))
    
    error = rnorm(n.s.i, sd = sd.s[[i]])
    y.cov.s <- X.cov.s %*% beta.s.i + error
    colnames(X.cov.s) = sprintf('var%s',1:ncol(X.cov.s))
    source.data.covshift[[i]] = list(y=y.cov.s, x=X.cov.s)
  }
  
  ## Simulate target coefficients
  beta.t.diffbetmin <- rep(0, p.t)
  beta.t.samebetmin <- rep(0, p.t)
  
  beta.t.diffbetmin[1:s.t] <- c(1/3,sample(c(-1,1),size=s.t-2,replace = TRUE)*runif(s.t-2, 1/2, 1.5),1/2)
  beta.t.samebetmin[1:s.t] <- c(1/3,sample(c(-1,1),size=s.t-2,replace = TRUE)*runif(s.t-2, 1/2, 1.5),1/3)
  
  ## Generate target predictors
  cov.t <- matrix(rep(0.5,p.t*p.t), ncol=p.t)
  diag(cov.t) <- rep(1,p.t)
  X.t <- mvtnorm::rmvnorm(n.t, sigma = cov.t)
  X.t <- X.t %*% diag(1/sqrt(diag(t(X.t) %*% X.t)/n.t))
  colnames(X.t) = sprintf('var%s',1:ncol(X.t))
  
  ## Simulate target response
  y.t.diffbetmin <- X.t %*% beta.t.diffbetmin + rnorm(n.t, sd = sd.t)
  target.data.diffbetmin = list(y=y.t.diffbetmin, x=X.t)
  
  y.t.samebetmin <- X.t %*% beta.t.samebetmin + rnorm(n.t, sd = sd.t)
  target.data.samebetmin = list(y=y.t.samebetmin, x=X.t)
  
  return(list(target.data.diffbetmin=target.data.diffbetmin,
              target.data.samebetmin=target.data.samebetmin,
              source.data.covshift = source.data.covshift,
              beta.t = data.frame(samebetmin = beta.t.samebetmin, diffbetmin = beta.t.diffbetmin),
              beta.s = store.beta.s))
}

set.seed(55)

for(n.t in c(50,seq(100,400,100),700)){
  cat('   simulating data for n_t=',n.t,'\n')
  dir.create(file.path(path, paste0('nt', n.t)), showWarnings = FALSE)
  path.nt = file.path(path, paste0('nt', n.t))
  
  s.t = 12 #2*round(log(n.t))
  p.t = 2*n.t
  args.sim = data.frame(n.t=c(n.t,NA,NA), p.t=c(p.t,NA,NA), s.t=c(s.t,NA,NA),
                        sd.t=c(1,NA,NA), n.s=c(1000,1000,1000), p.s = c(p.t,NA,NA),
                        sd.s=c(1,1,1), trans.inf = c('v.informative','m.informative', 'adversarial'))

  write.csv(args.sim,paste0(path.nt,'/args.sim.csv'),row.names = FALSE)
  n.sim = 100
  for(i in 1:n.sim){
    sim.data.res = data.sim(n.t=args.sim$n.t[1], p.t=args.sim$p.t[1], s.t=args.sim$s.t[1], sd.t=args.sim$sd.t[1],
                            n.s=args.sim$n.s, p.s = args.sim$p.s[1], sd.s=args.sim$sd.s,
                            trans.inf = args.sim$trans.inf)
    saveRDS(sim.data.res, paste0(path.nt,"/sim_",i,".rds"))
  }
}
