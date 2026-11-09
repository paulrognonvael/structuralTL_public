### Run the simulations ###
path = "~/structuralTL/simulations/"

TL.sim <- function(scenario, n.val=NULL, path, nb.sim=100){
  require(mvtnorm)
  require(tidyverse)
  require(ggplot2)
  require(modelSelection)
  require(stringr)
  require(glmnet)
  require(ncvreg)
  require(parallel)
  require(glmtrans)
  
  # Working directory
  dir.create(paste0(path,'scenario',scenario,'/'), showWarnings = FALSE)
  path.scen = paste0(path,'scenario',scenario,'/')
  
  # Helper functions
  source(paste0(path,"routines.sim.app.R"))
  
  # Scenario

  # Lists of methods
  l0method.vec <- c('EBIC','S.A.sd','S.Abis.sd',
                    'S.A.TL', 'S.Abis.TL')
  
  method.vec <- c(l0method.vec,"lasso.cv", "scad.cv",'trans.glm')


  ###################################################
  
  if(is.null(n.val)){
    values.n <- c(50,seq(100,700,100))
  }else{
    values.n <- n.val
  }
  
  for (n  in values.n){
    cat('n=',n,'\n')
    t0 <- Sys.time()
    
    sim.result.sel <- data.frame()
    sim.result.est.mse <- data.frame()
    sim.result.cvmse <- data.frame()
    
    for(i in 1:nb.sim){
      cat('\n','   iteration ',i,'\n')
      sim.data = readRDS(paste0(path,'sim.data/nt',n,'/sim_',i,'.rds'))
      
      
      if(scenario %in% c(1,2,3,4)){
        target.data = sim.data$target.data.diffbetmin
        beta_star = sim.data$beta.t$diffbetmin
      }
      
      if(scenario %in% c(5,6)){
        target.data = sim.data$target.data.samebetmin
        beta_star = sim.data$beta.t$samebetmin
      }
      
      if(scenario %in% c(1,5)){
        source.data = list(sim.data$source.data.covshift[[1]])
      }
      
      if(scenario %in% c(2,6)){
        source.data = list(sim.data$source.data.covshift[[3]])
      }
      
      if(scenario %in% c(3)){
        source.data = list(sim.data$source.data.covshift[[1]],
                           sim.data$source.data.covshift[[2]])
      }
      
      if(scenario %in% c(4)){
        source.data = list(sim.data$source.data.covshift[[1]],
                           sim.data$source.data.covshift[[3]])
      }
      
      res.sel.l0 <- ATL0andOtherl0(target_y = target.data$y[,1], target_X = target.data$x, 
                                   source_data = source.data, source_sel=NULL)
      res.est.mse.l0 <- sapply(res.sel.l0$sel[,l0method.vec], mse.est, y=target.data$y[,1], 
                               X.design=target.data$x, beta_star= beta_star)
      res.est.mse.l0 <- data.frame(t(res.est.mse.l0))
      colnames(res.est.mse.l0) <- l0method.vec
      res.lasso.scad.transglm <- lasso.scad.transglm.cv(y=target.data$y[,1], X.design=target.data$x, 
                                                        beta_star= beta_star, source.data = source.data)
      res.sel <- data.frame(res.sel.l0$sel,res.lasso.scad.transglm$sel,sim=i,n=n, soursel.glmtrans=res.lasso.scad.transglm$soursel.glmtrans)
      sim.result.sel <- bind_rows(sim.result.sel,res.sel)
      
      res.est.mse <- data.frame(res.est.mse.l0,res.lasso.scad.transglm$est.mse,sim=i,n=n)
      sim.result.est.mse <- bind_rows(sim.result.est.mse,res.est.mse)
      
      
      res.mse.l0 <- cvmse.l0.comp(target.data$y[,1], target.data$x, res.sel.l0$blocks, K=5, mc.cores=4)
      res.mse <- data.frame(res.mse.l0$cvmse,res.lasso.scad.transglm$cvmse,sim=i,n=n)
      sim.result.cvmse <- bind_rows(sim.result.cvmse,res.mse)
    }
    
    
    sim.result.df <- pivot_longer(sim.result.sel,cols=method.vec,names_to = 'method',values_to = 'sel.model')
  
    sim.result.est.mse.df <- pivot_longer(sim.result.est.mse,cols=method.vec,names_to = 'method',values_to = 'est.mse')
    sim.result.df <- merge(sim.result.df, sim.result.est.mse.df, by = c('method','sim','n'))
      
    sim.result.cvmse.df <- pivot_longer(sim.result.cvmse,cols=method.vec,names_to = 'method',values_to = 'cv.mse')
    sim.result.df <- merge(sim.result.df, sim.result.cvmse.df, by = c('method','sim','n'))
    
    sim.result.df <- postprocess.sel.sim(sim.result.df=sim.result.df, beta_star= beta_star)
    
    write.csv(sim.result.df,paste0(path.scen,"sim.result.scenario",scenario,".n",n,".csv"),row.names = FALSE)
    
    t1 <- Sys.time() 
    cat('Running time n=',n,':', round(difftime(t1, t0, units = "mins"),3),'minutes,')
    cat('Time:', format(Sys.time(), "%X")); cat('\n');
    
  }
}

set.seed(1278)
for(scenario in 5){
  cat('--------------------------------------------------------------','\n')
  cat('SCENARIO',scenario,'\n')
  cat('--------------------------------------------------------------','\n')
  TL.sim(scenario,n.val=#c(50,seq(100,400,100),700)
           400,path=path, nb.sim=50)
}
set.seed(1278)
for(scenario in 5:6){
  cat('--------------------------------------------------------------','\n')
  cat('SCENARIO',scenario,'\n')
  cat('--------------------------------------------------------------','\n')
  TL.sim(scenario,n.val=#c(50,seq(100,400,100),700)
           700,path=path, nb.sim=50)
}
