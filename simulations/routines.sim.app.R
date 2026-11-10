
##########################################################################
# . MISCELLANEOUS HELPER FUNCTIONS
##########################################################################

## Converts a string of variable indices separated by commas into a vector
model.char2vec <- function(model.char){
  return(as.numeric(strsplit(model.char,",")[[1]]))
}

## Gets the indices corresponding to nonzero coeff lasso estimates 
get.model <- function(glmnet.beta){
  return(paste(unname(which(glmnet.beta!=0)),collapse=","))
}

##########################################################################
# . FUNCTIONS TO COMPUTE BLOCK PARATION AND l0 SCORES
##########################################################################

## Outputs a block partition from a list of subsets selected on source 
fromSelToBlocks = function(source_sel,cols_target_X){
  #Input:
  # source_sel: a list of vectors of names of target covariates found relevant in source datasets
  # col_target_X: names of target covariates
  #Output:
  # blocks: a list where each element is a vector of covariates corresponding to a block
  
  
  # converting source selection subsets into indicator value for each target covariate
  for(k in 1:length(source_sel)){
    source_sel[[k]] <- cols_target_X %in% source_sel[[k]]
  }
  # Binding indicator vectors in a single matrix
  source_sel.concat = data.frame(do.call('cbind',source_sel)) 
  colnames(source_sel.concat) = sprintf('source%s',1:length(source_sel))
  # Get a source selection value for each target covariate e.g 1,0,1
  comb = apply(source_sel.concat, MARGIN = 1, FUN=paste, collapse=',')
  # Partition target covariates by source selection value (e.g 1,0,1)
  source_sel.concat$block = factor(comb)
  source_sel.concat$var = 1:nrow(source_sel.concat)
  blocks = lapply(levels(source_sel.concat$block),FUN =function (l) which(source_sel.concat$block==l))
  return(blocks)
}

## Returns the l0 scores for a given model and block partition
l0score <- function(model, y, X.design, blocks, s.j.hat){
  # Input:
  # - model: vector of indices of target covariates
  # - y: target outcome
  # - X.design: target covariates
  # - block: list where each element is a vector of indices of target covariates forming a block
  # - s.j.hat: block estimates of number of nonzero coefficients
  # Output:
  # a vector with model name and l0 scores: S^ATL (kappa_j^A), S^ATL (kappa_j^AEB), S^asd (kappa^A), S^asd (kappa^AEB), and selected by EBIC
  
  # initialize quantities
  beta.hat <- rep(0,ncol(X.design))
  k <- 0
  k.j <-rep(0,length(blocks))
  n <- length(y)
  
  # compute number of variables in each block in the model
  if (sum(is.na(model))!=length(model)){
    model <- model[!is.na(model)]
    X.model <- X.design[,model]
    beta.hat[model] <- solve(t(X.model)%*%X.model)%*%t(X.model) %*% y
    k.j <- sapply(blocks, function(b){sum(beta.hat[b] !=0)})
    k <- sum(beta.hat !=0)
  }
  
  # compute estimate of total number of nonzero coeff and prediction by model
  s.hat <- sum(s.j.hat)
  pred <- X.design%*%beta.hat
  
  # score EBIC
  sc.ebic <- 0.5*n*log(sum((y-pred)^2)/n) + 0.5*k*log(n) + log(choose(ncol(X.design),k))
  
  # score S^asd (kappa^AEB)
  sc.A.sd <- 0.5*n*log(sum((y-pred)^2)/n) + 0.5*k*log(n) + k*log(ncol(X.design)/s.hat-1)

  # score S^asd (kappa^A)
  sc.Abis.sd <- 0.5*n*log(sum((y-pred)^2)/n) + 0.5*k*log(n) + k*log(ncol(X.design)-s.hat)
  
  # score S^ATL (kappa_j^AEB)
  penAEB.j <- rep(NA,length(blocks))
  for(j in 1:length(blocks)){
    penAEB.j[j] = k.j[j]*log(ncol(X.design[,blocks[[j]]])/s.j.hat[j]-1)
  }
  sc.A.TL <- 0.5*n*log(sum((y-pred)^2)/n) + 0.5*k*log(n) + sum(penAEB.j)
  
  # score S^ATL (kappa_j^A)
  penA.j <- rep(NA,length(blocks))
  for(j in 1:length(blocks)){
    penA.j[j] = k.j[j]*log(ncol(X.design[,blocks[[j]]])-s.j.hat[j])
  }
  sc.Abis.TL <- 0.5*n*log(sum((y-pred)^2)/n) + 0.5*k*log(n) + sum(penA.j)
  
  # model short name
  models.val <- paste(model, collapse = "," )
  
  return(c(models.val, sc.ebic,
           sc.A.sd, sc.Abis.sd,
           sc.Abis.TL, sc.A.TL))
}

## Returns the l0 scores for all models in a given list
compute.l0scores<-function(models, y, X.design, blocks, s.j.hat){
  l0score.df<- apply(models,1,l0score,y=y,X.design=X.design,
                     blocks=blocks, s.j.hat=s.j.hat)
  l0score.df <-data.frame(t(l0score.df))
  l0score.df[,-1] <- sapply(l0score.df[, -1], as.numeric)
  colnames(l0score.df)<-c('model','EBIC','S.A.sd','S.Abis.sd',
                          'S.A.TL', 'S.Abis.TL')
  return(l0score.df)
}

## For a given block partition
## Uses MCMC algorithm to score models by initial l0 penalty kappa^0 = 1/2 log(n) + log(p)
## estimates block number of nonzero coefficients by sum of pseudo posterior probabilities
## scores models visited by MCMC by the penalties S^ATL (kappa_j^A), S^ATL (kappa_j^AEB), S^asd (kappa^A), S^asd (kappa^AEB), and EBIC 
## Returns S^ATL (kappa_j^A), S^ATL (kappa_j^AEB), S^asd (kappa^A), S^asd (kappa^AEB), and selected by EBIC for a given block partition 
selectionl0 <-function(y, X.design, blocks){
  # Input:
  # y: target outcome
  # X.design: target covariates
  # blocks: block partition
  #
  # Output:
  # S^ATL (kappa_j^A), S^ATL (kappa_j^AEB), S^asd (kappa^A), S^asd (kappa^AEB), and selected by EBIC
  # l0 scores of models visited by MCMC for all the penalties
  require(modelSelection)
  t0 <- Sys.time()
  
  # Get top models with kappa.o
  best.kappa.o.res <- modelSelection::bestIC(y, X.design, penalty= log(nrow(X.design))+2*log(ncol(X.design)), 
                              verbose=FALSE, maxvars= nrow(X.design)
                              #vmaxvars=nrow(X.design)*0.9,
                              # deltaini=ini,
  )
  sel.kappa.o <- best.kappa.o.res$models$modelid[1]
  top.kappa.o.list <- distinct(best.kappa.o.res$models[,1])
  
  # Add models visited by lasso
  glmnet.res<- cv.glmnet(X.design,y,intercept = FALSE)
  models.glmnet <- data.frame(modelid=apply(glmnet.res$glmnet.fit$beta,2,get.model))
  top.kappa.o.list<- distinct(rbind(top.kappa.o.list,models.glmnet))
  
  # Format list of models
  max.model.size <- max(sapply(str_split(top.kappa.o.list$modelid, ","),length))
  top.kappa.o.fmted <- str_split_fixed(top.kappa.o.list$modelid, ",",max.model.size)
  top.kappa.o.fmted[top.kappa.o.fmted == ""] <- NA
  class(top.kappa.o.fmted) <- "numeric"
  
  # sanity check can't have models with more than n variables
  top.kappa.o.fmted <- top.kappa.o.fmted[,1:min(nrow(X.design)-1,ncol(top.kappa.o.fmted))]
  
  # Estimate s, s.hat, s.j.hat
  margprob <- best.kappa.o.res$msfit$margpp
  s.hat <- sum(margprob,na.rm=TRUE)
  s.j.hat <- sapply(blocks, function(b){sum(margprob[b],na.rm=TRUE)})
  
  # Cleaning
  rm(best.kappa.o.res)
  
  # Compute l0 scores for 
  l0score.df <- compute.l0scores(top.kappa.o.fmted, y, X.design, blocks, s.j.hat)
  
  t1 <- Sys.time()
  
  time.secs <- round(difftime(t1, t0, units = "secs"),3)
  
  sel.ebic <- l0score.df[which.min(l0score.df$EBIC),"model"]
  sel.S.A.sd <- l0score.df[which.min(l0score.df$S.A.sd),"model"]
  sel.S.Abis.sd <- l0score.df[which.min(l0score.df$S.Abis.sd),"model"]
  sel.S.A.TL <- l0score.df[which.min(l0score.df$S.A.TL),"model"]
  sel.S.Abis.TL   <- l0score.df[which.min(l0score.df$S.Abis.TL),"model"]
  
  
  ## formating results 
  method.vec <- c(
    'EBIC','S.A.sd','S.Abis.sd',
    'S.A.TL', 'S.Abis.TL'
  )
  
  sel.model.res <-c(
    sel.ebic,
    sel.S.A.sd,
    sel.S.Abis.sd,
    sel.S.A.TL,
    sel.S.Abis.TL
  )
  
  res.all <- c(sel.model.res,
               s.j.hat,
               s.hat,
               time.secs)
  
  df <- data.frame(t(res.all))
  colnames(df) <- c(method.vec,
                    sprintf('s%s.hat',1:length(blocks)),
                    's.hat',
                    'time.l0')
  
  gc()
  
  return(list(sel=df,scores= l0score.df, blocks=blocks))
}

##########################################################################
# . WRAPPER FUNCTIONS RETURNING MODEL SELECTED by ASTLO, standard l0 and EBIC
##########################################################################

# Gets the block partition and returns S^ATL (kappa_j^A), S^ATL (kappa_j^AEB), S^asd (kappa^A), S^asd (kappa^AEB), and selected by EBIC
ATL0andOtherl0 <- function(target_y, target_X, source_sel=NULL, source_data=NULL){
  require(modelSelection)
  
  # Construction of the block partition 
  source_sel_all = source_sel
  if(!is.null(source_data)){
    cat('   Selecting variables in the source datasets\n')
    source.data_sel <- vector(mode = "list", length = length(source_data))
    
    ## Select subset in each source dataset provided
    for(k in 1:length(source_data)){
      source.X = source_data[[k]]$x
      source.y = source_data[[k]]$y[,1]
      bestIC.k <- modelSelection::bestIC(source.y, source.X, 
                                penalty= log(nrow(source.X))+2*log(ncol(source.X)), 
                                verbose=FALSE, maxvars= nrow(source.X)
                                #vmaxvars=nrow(X.design)*0.9,# deltaini=ini,
      )
      ### retrieving variable names of selected subset
      source.data_sel[[k]] <- colnames(source.X)[model.char2vec(bestIC.k$models$modelid[1])]
    }
    ## Adding subsets selected on source data
    source_sel_all = append(source_sel_all,source.data_sel)
  }  

  ## Get block partition
  blocks = fromSelToBlocks(source_sel_all,colnames(target_X))
  
  # Run selectionl0 for the block partition
  cat('   Selecting variables in the target dataset\n')
  res <- selectionl0(y = target_y, X.design = target_X, blocks= blocks)
  names(source_sel_all) = 'source_sel_all'
  res <- append(res, source_sel_all)
  return(res)
}


##########################################################################
# . FUNCTIONS TO PROCESS SIMULATION OUTPUTS
##########################################################################

## Returns MSE in coef estimates for a selected model
mse.est <- function(model,y, X.design, beta_star){
  
  beta.est <- rep(0,length(beta_star))
  X.design.df <- data.frame(X.design)
  model.vec <- model.char2vec(model)
  if(sum(is.na(model.vec))>0){
    mse <- NA
  }else{
    data<-cbind(data.frame(y=y),X.design.df[,model.vec,drop=FALSE])
    lin.reg <- lm(y~.+0,data=data)
    beta.est[model.vec] <- coef(lin.reg)
    mse <- sum((beta.est-beta_star)^2)/length(beta_star) 
  }
  return(mse)
}

## Returns CV prediction MSE with l0 penalties.
cvmse.l0.comp <- function(y, X.design, blocks, K, mc.cores=1){
  require(parallel)
  if (K > nrow(X.design)) stop("The number of folds cannot be larger than nrow(x)")
  subset <- rep(1:K,ceiling(nrow(X.design)/K))[1:nrow(X.design)]
  subset <- sample(subset,size=nrow(X.design),replace=FALSE)
  f <- function(k,...) {
    sel <- subset==k
    
    sel.l0 <- selectionl0(y=y[!sel],X.design[!sel,,drop=FALSE], 
                               blocks)$sel[,c('EBIC','S.A.sd','S.Abis.sd',
                                              'S.A.TL', 'S.Abis.TL')]
    
    sel.list <- lapply(sel.l0, model.char2vec)
    
    pred.f <- function(model,y, X.design, sel){
      X.design.df <- data.frame(X.design)
      if(sum(is.na(model))==0){
        data<-cbind(data.frame(y=y[!sel]),X.design.df[!sel,model,drop=FALSE])
        lin.reg <- lm(y~.+0,data=data)
        pred <- predict(lin.reg, newdata= X.design.df[sel,model,drop=FALSE])
      }else{
        pred <- rep(mean(y[!sel]),nrow(X.design.df[sel,]))
      }
      return(pred)
    }
    
    pred.subset <- sapply(sel.list, pred.f, y, X.design, sel)
    
   }
  
  if (mc.cores > 1) {
    if ("parallel" %in% loadedNamespaces())  {
      allpred <- parallel::mclapply(1:K, f, mc.preschedule=FALSE)
    } else {
      stop("Did not find mclapply. Please load parallel package")
    }
  } else {
    allpred <- lapply(1:K, f)
  }
    
  pred= data.frame(matrix(NA,nrow=nrow(X.design), ncol= length(c('EBIC','S.A.sd','S.Abis.sd',
                                                                 'S.A.TL', 'S.Abis.TL'))))
  for (k in 1:K) {
      pred[subset==k,] <- allpred[[k]]
  }
    
  res.cvmse <- colMeans((pred-y)^2)
  df.cvmse <- data.frame(t(res.cvmse))

  res.rsquared <- (apply(pred,2,cor,y))^2
  df.rsquared <- data.frame(t(res.rsquared))
  
  colnames(df.cvmse) <- c('EBIC','S.A.sd','S.Abis.sd',
                          'S.A.TL', 'S.Abis.TL')
  colnames(df.rsquared) <- c('EBIC','S.A.sd','S.Abis.sd',
                             'S.A.TL', 'S.Abis.TL')
  return(list(rsquared=df.rsquared,cvmse=df.cvmse))
  }


## Returns model selected, estimation MSE and CV prediction MSE with CV LASSO, SCAD and Trans-GLM
lasso.scad.transglm.cv <- function(y, X.design, trans.glm = TRUE, source.data=NULL, beta_star=NULL){
  require(glmnet)
  require(ncvreg)
  cat('   Fitting lasso','\n')
  ## CV-LASSO 
  t0lasso <- Sys.time()
  glmnet.res<- cv.glmnet(X.design,y,intercept = FALSE)
  t1lasso <- Sys.time()
  sel.lasso.cv <- paste(which(coef(glmnet.res, s=glmnet.res$lambda.min)[-1]!=0),collapse=",")
  cv.mse.lasso <- glmnet.res$cvm[which(glmnet.res$lambda==glmnet.res$lambda.min)]
  time.secs.lasso <- round(difftime(t1lasso, t0lasso, units = "secs"),3)
  
  cat('   Fitting SCAD','\n')
  ## SCAD
  t0scad <- Sys.time()
  cvscad <- cv.ncvreg(X=X.design,y=y,family="gaussian",penalty="SCAD",nfolds=10,dfmax=1000,max.iter=10^4)
  t1scad <- Sys.time()
  sel.scad.cv <- paste(which(cvscad$fit$beta[-1,cvscad$min]!=0),collapse=",")
  cv.mse.scad <- cvscad$cve[which(cvscad$lambda==cvscad$lambda.min)]
  time.secs.scad <- round(difftime(t1scad, t0scad, units = "secs"),3)
  
  ## Trans-glm
  if(trans.glm){
    cat('   Fitting Trans-glm','\n')
    require(glmtrans)
    target.data = list(y=y, x=X.design)
    t0transglm <- Sys.time()
    res.glmtrans = glmtrans(target = target.data, source = source.data, 
                             family = 'gaussian', intercept = FALSE, alpha=1)
    t1transglm <- Sys.time()
    sel.glmtrans <- paste(which(res.glmtrans$beta[-1]!=0),collapse=",")
    soursel.glmtrans = ifelse(length(res.glmtrans$transfer.source.id)==0,'none',paste0(res.glmtrans$transfer.source.id,collapse=','))
    cv.mse.glmtrans = res.glmtrans$fitting.list$target.valid.loss
    time.secs.glmtrans <- round(difftime(t1transglm , t0transglm, units = "secs"),3)
    
    ## Gathering selection and cv prediction mse results
    df.sel <- data.frame(t(c(sel.lasso.cv, sel.scad.cv, sel.glmtrans, time.secs.lasso,
                             time.secs.scad, time.secs.glmtrans))) 
    df.cvmse <- data.frame(t(c(cv.mse.lasso, cv.mse.scad, cv.mse.glmtrans)))
    colnames(df.sel) <- c("lasso.cv","scad.cv","trans.glm","time.lasso.cv","time.scad.cv","time.trans.glm")
    colnames(df.cvmse) <- c("lasso.cv","scad.cv","trans.glm")
    
    ## Estimation mse results
    if(!is.null(beta_star)){
      est.mse.cvlasso <- sum((coef(glmnet.res, s=glmnet.res$lambda.min)[-1]-beta_star)^2)/length(beta_star)
      est.mse.cvscad <- sum((cvscad$fit$beta[-1,cvscad$min]-beta_star)^2)/length(beta_star)
      est.mse.glmtrans = sum((res.glmtrans$beta[-1]-beta_star)^2)/length(beta_star)
      df.est.mse <- data.frame(t(c(est.mse.cvlasso, est.mse.cvscad,est.mse.glmtrans)))
      colnames(df.est.mse) <- c("lasso.cv","scad.cv","trans.glm")
      return(list(sel=df.sel,cvmse=df.cvmse,est.mse = df.est.mse, soursel.glmtrans = soursel.glmtrans))
    }
    return(list(sel=df.sel,cvmse=df.cvmse, soursel.glmtrans = soursel.glmtrans))
    
  } else {
    ## Gathering selection and cv prediction mse results
    df.sel <- data.frame(t(c(sel.lasso.cv, sel.scad.cv, time.secs.lasso,
                             time.secs.scad))) 
    df.cvmse <- data.frame(t(c(cv.mse.lasso, cv.mse.scad)))
    colnames(df.sel) <- c("lasso.cv","scad.cv","time.lasso.cv","time.scad.cv")
    colnames(df.cvmse) <- c("lasso.cv","scad.cv")
    
    ## Estimation mse results
    if(!is.null(beta_star)){
      est.mse.cvlasso <- sum((coef(glmnet.res, s=glmnet.res$lambda.min)[-1]-beta_star)^2)/length(beta_star)
      est.mse.cvscad <- sum((cvscad$fit$beta[-1,cvscad$min]-beta_star)^2)/length(beta_star)
      df.est.mse <- data.frame(t(c(est.mse.cvlasso, est.mse.cvscad)))
      colnames(df.est.mse) <- c("lasso.cv","scad.cv")
      return(list(sel=df.sel,cvmse=df.cvmse,est.mse = df.est.mse))
    }
    return(list(sel=df.sel,cvmse=df.cvmse))
  }
  

}

## Returns the number of false positive and false negatives
nb.error<-function(sel.vect, true.beta, type, block=NULL){
  truth <- true.beta != 0
  sel.model <- rep(FALSE, length(true.beta))
  sel.model[sel.vect] <- TRUE
  
  if(is.null(block)){
    block=1:length(true.beta)
  }
  
  if(type=='I'){
    tI <- truth==FALSE & sel.model==TRUE
    return(sum(tI[block]))
  } else if (type=='II'){
    tII <- truth==TRUE & sel.model==FALSE
    return(sum(tII[block]))
  } else{
    stop("Invalid type of error")
  }
}

## Returns false positives, false negatives, FDR, power from a simulation output
postprocess.sel.sim <- function(sim.result.df,beta_star,blocks=NULL){
  sim.result.df['recovery']<-  sim.result.df$sel.model==paste(as.character(which(beta_star!=0)),collapse=",")
  sel.list <- lapply(sim.result.df$sel.model, model.char2vec)
  if(!is.null(blocks)){
    for(i in 1:length(blocks)){
      t1.blocki = paste0('nb.tI.b',j)
      t2.blocki = paste0('nb.tII.b',j)
      sim.result.df[t1.blocki] <- sapply(sel.list, nb.error, true.beta=beta_star, type='I', block=blocks[[j]])
      sim.result.df[t2.blocki] <- sapply(sel.list, nb.error, true.beta=beta_star, type='II', block=blocks[[j]])
    }
  }
  sim.result.df['nb.tI'] <- sapply(sel.list, nb.error, true.beta=beta_star, type='I') 
  sim.result.df['nb.tII'] <- sapply(sel.list, nb.error, true.beta=beta_star, type='II') 
  sim.result.df['est.size'] <- sapply(sel.list,length)
  sim.result.df <- sim.result.df %>% mutate(FP=nb.tI, TP = sum(beta_star!=0)-nb.tII, 
  FDR = FP/(FP+TP), power = TP/sum(beta_star!=0))
  return(sim.result.df)
}
