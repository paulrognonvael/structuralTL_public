################################################################################
################################### ASTL0 ######################################
################################################################################

ASTL0<- function(target_y, target_X, source_sel=NULL, source_data=NULL){
  # Input
  # - target_y: a vector with target outcome 
  # - target_y: a matrix with target covariates 
  # - source_sel: a list of vectors of indices of target covariates relevant in the source
  # - source_data: a list of source datasets, each source dataset is itself a list
  #             where the outcome is a vector named 'y' and the covariates are in a
  #             matrix named 'x'
  # One should provide either source_sel or source_data or both 
  #
  # Output
  # - sel: a data.frame with indices of selected target covariates, estimates block
  #   block sparsities, and runnning time
  # - scores: the scores of models last visited by MCMC
  # - blocks: the block partition fitted from source_sel and source_data
  # - source_sel_all: a list of indices of target covariates found relevant in the source data 
  #   (source_sel and those selected in source_data if non empty)
  require(modelSelection)
  require(dplyr)
  require(glmnet)
  source_sel_all = source_sel
  if(!is.null(source_data)){
    cat('   Selecting variables in the source datasets\n')
    source.data_sel <- vector(mode = "list", length = length(source_data))
    
    for(k in 1:length(source_data)){
      # Selecting subset in each source 
      source.X = source_data[[k]]$x
      source.y = source_data[[k]]$y[,1]
      bestIC.k <- modelSelection::bestIC(source.y, source.X, 
                                penalty= log(nrow(source.X))+2*log(ncol(source.X)), 
                                verbose=FALSE, maxvars= nrow(source.X)
                                #vmaxvars=nrow(X.design)*0.9,# deltaini=ini,
      )
      # retrieving variable names of selected subset
      source.data_sel[[k]] <- colnames(source.X)[model.char2vec(bestIC.k$models$modelid[1])]
    }
    ## Adding subsets selected on source data
    source_sel_all = append(source_sel_all,source.data_sel)
  }  

  ## Getting block partition
  blocks = fromSelToBlocks(source_sel_all,colnames(target_X))
  
  cat('   Selecting variables in the target dataset\n')
  res <- selectionl0(y = target_y, X.design = target_X, blocks= blocks)
  names(source_sel_all) = 'source_sel_all'
  res <- append(res, source_sel_all)
  return(res)
}

################################################################################
############################ Helper functions ##################################
################################################################################

# Converts a string of variable indices separated by commas into a vector
model.char2vec <- function(model.char){
  return(as.numeric(strsplit(model.char,",")[[1]]))
}

# Gets the indices corresponding to nonzero coeff estimates 
get.model <- function(glmnet.beta){
  return(paste(unname(which(glmnet.beta!=0)),collapse=","))
}


# Outputs a block partition from a list of subsets selected on source 
fromSelToBlocks = function(source_sel,cols_target_X){
  # converting source selection subsets into 0,1 vector 
  for(k in 1:length(source_sel)){
    source_sel[[k]] <- cols_target_X %in% source_sel[[k]]
  }
  source_sel.concat = data.frame(do.call('cbind',source_sel)) 
  colnames(source_sel.concat) = sprintf('source%s',1:length(source_sel))
  comb = apply(source_sel.concat, MARGIN = 1, FUN=paste, collapse=',')
  source_sel.concat$block = factor(comb)
  source_sel.concat$var = 1:nrow(source_sel.concat)
  blocks = lapply(levels(source_sel.concat$block),FUN =function (l) which(source_sel.concat$block==l))
  return(blocks)
}



## Returns the l0 score for a given model
l0score <- function(model, y, X.design, blocks,s.j.hat){
  beta.hat <- rep(0,ncol(X.design))
  k <- 0
  k.j <-rep(0,length(blocks))
  n <- length(y)
  
  if (sum(is.na(model))!=length(model)){
    model <- model[!is.na(model)]
    X.model <- X.design[,model]
    beta.hat[model] <- solve(t(X.model)%*%X.model)%*%t(X.model) %*% y
    k.j <- sapply(blocks, function(b){sum(beta.hat[b] !=0)})
    k <- sum(beta.hat !=0)
  }
  
  s.hat <- sum(s.j.hat)
  
  pred <- X.design%*%beta.hat
  
  # S^ATL0
  penA.j <- rep(NA,length(blocks))
  for(j in 1:length(blocks)){
    penA.j[j] = k.j[j]*log(ncol(X.design[,blocks[[j]]])-s.j.hat[j])
  }
  sc.A.TL <- 0.5*n*log(sum((y-pred)^2)/n) + 0.5*k*log(n) + sum(penA.j)
  
  models.val <- paste(model, collapse = "," )
  
  return(c(models.val, sc.A.TL))
}

## Returns the l0 score for all models in a given list
compute.l0scores<-function(models, y, X.design, blocks, s.j.hat){
  l0score.df<- apply(models,1,l0score,y=y,X.design=X.design,
                     blocks=blocks, s.j.hat=s.j.hat)
  l0score.df <-data.frame(t(l0score.df))
  l0score.df[,-1] <- sapply(l0score.df[, -1], as.numeric)
  colnames(l0score.df)<-c('model', 'sc.A.TL')
  return(l0score.df)
}


## Returns S^ATL for a given block partition 
selectionl0 <-function(y, X.design, blocks){
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
  
  # Estimate s, s.hat, sj.hat
  margprob <- best.kappa.o.res$msfit$margpp
  s.hat <- sum(margprob,na.rm=TRUE)
  s.j.hat <- sapply(blocks, function(b){sum(margprob[b],na.rm=TRUE)})
  
  # Cleaning
  rm(best.kappa.o.res)
  
  ## selected models
  
  l0score.df <- compute.l0scores(top.kappa.o.fmted, y, X.design, blocks, s.j.hat)
  
  t1 <- Sys.time()
  
  time.secs <- round(difftime(t1, t0, units = "secs"),3)
  
  sel.ASTL0 <- l0score.df[which.min(l0score.df$sc.A.TL),"model"]
  
  
  ## formating results 
  method.vec <- c(
    'sel.ASTL0'
  )
  
  sel.model.res <-c(
    sel.ASTL0
  )
  
  res.all <- c(sel.model.res,
               s.j.hat,
               s.hat,
               time.secs)
  print(res.all)
  df <- data.frame(t(res.all))
  colnames(df) <- c(method.vec,
                    sprintf('s%s.hat',1:length(blocks)),
                    's.hat',
                    'time.l0')
  
  gc()
  
  return(list(sel=df,scores= l0score.df, blocks=blocks))
}
