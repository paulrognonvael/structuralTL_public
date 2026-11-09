library(tidyverse)
library(modelSelection)
library(stringr)
library(glmnet)
library(ncvreg)
library(parallel)
#BiocManager::install("hgu133plus2.db")
library(hgu133plus2.db)


##### Paths #####
path = "~/structuralTL/"
#path = "C:/Users/Usuario/Downloads/"
source(paste0(path,"simulations/routines.sim.app.R"))
#source(paste0(path,"colon cancer/routines.sim.app.R"))

##### Load data #####
x = read.table(paste0(path,"colon cancer/tgfb.txt"), header=TRUE)
shortlist= as.character(read.table(paste0(path,"colon cancer/mouse_shortlist.txt"), header=TRUE)[,1])
y = as.matrix(x[,1])
X.design = as.matrix(x[,-1])
mouselist = which(colnames(X.design) %in% shortlist)
S1 = colnames(X.design)[colnames(X.design) %in% shortlist]

# Lists of methods
l0method.vec <- c('EBIC','S.A.sd','S.Abis.sd',
                  'S.A.TL', 'S.Abis.TL')

method.vec <- c(l0method.vec,"lasso.cv", "scad.cv")

set.seed(954)

##### Gene selection #####
res.sel.l0 <- ATL0andOtherl0(target_y=y, target_X=X.design, source_sel = list(S1))
res.lasso.scad <- lasso.scad.transglm.cv(y,X.design=X.design,trans.glm = FALSE)
res.sel <- data.frame(res.sel.l0$sel,res.lasso.scad$sel)

res.cvmse.l0 <- cvmse.l0.comp(y, X.design, blocks=res.sel.l0$blocks, K=10, mc.cores=1)
res.cvmse <- data.frame(res.cvmse.l0$cvmse,res.lasso.scad$cvmse)

result.cvmse.df <- pivot_longer(res.cvmse,cols=method.vec,names_to = 'method',values_to = 'cv.mse')
result.sel.df <- pivot_longer(res.sel,cols=method.vec,names_to = 'method',values_to = 'sel.model')
result.df <- merge(result.sel.df, result.cvmse.df, by = c('method'))
write.csv(result.df, paste0(path,"colon cancer/selected.cvmse.csv"), row.names = FALSE)

##### Retrieving gene names #####
get_annotation_by_block <- function(index.sel.ids,all.ids,mouselist){
  # In mouse list
  index <- index.sel.ids[index.sel.ids%in%mouselist]
  b <- all.ids[index]
  probe_ids= sub("X", "", b)
  annotations= AnnotationDbi::select(hgu133plus2.db, keys = probe_ids, columns = c("SYMBOL", "GENENAME"), keytype = "PROBEID",multiVals='first')
  tab1= cbind(index,annotations)
  rownames(tab1)= NULL
  cat('In mouse list\n')
  print(tab1)
  
  # Not in mouse list 
  index <- index.sel.ids[!index.sel.ids%in%mouselist]
  b <- all.ids[index]
  probe_ids= sub("X", "", b)
  annotations= AnnotationDbi::select(hgu133plus2.db, keys = probe_ids, columns = c("SYMBOL", "GENENAME"), keytype = "PROBEID",multiVals='first')
  annotations <- annotations %>% distinct(PROBEID,.keep_all = TRUE)
  
  tab2= cbind(index,annotations)
  rownames(tab2)= NULL
  cat('Not in mouse list\n') 
  print(tab2)
  res <- rbind(c('Mouse list',NA,NA,NA),tab1,c('Not Mouse list',NA,NA,NA),tab2)
  return(res)
} 

cat('Selected by EBIC\n')
write.csv(get_annotation_by_block(model.char2vec(res.sel$EBIC),colnames(X.design),mouselist), 
          paste0(path,"colon cancer/selected.genes.EBIC.annot.csv"), row.names = FALSE)

cat('Selected by S.A.sd\n')
write.csv(get_annotation_by_block(model.char2vec(res.sel$S.A.sd),colnames(X.design),mouselist), 
          paste0(path,"colon cancer/selected.genes.S.A.sd.annot.csv"), row.names = FALSE)

cat('Selected by S.A.TL\n')
write.csv(get_annotation_by_block(model.char2vec(res.sel$S.A.TL),colnames(X.design),mouselist), 
          paste0(path,"colon cancer/selected.genes.S.A.TL.annot.csv"), row.names = FALSE)

cat('Selected by S.Abis.sd\n')
write.csv(get_annotation_by_block(model.char2vec(res.sel$S.Abis.sd),colnames(X.design),mouselist), 
          paste0(path,"colon cancer/selected.genes.S.Abis.sd.annot.csv"), row.names = FALSE)

cat('Selected by S.Abis.TL\n')
write.csv(get_annotation_by_block(model.char2vec(res.sel$S.Abis.TL),colnames(X.design),mouselist), 
          paste0(path,"colon cancer/selected.genes.S.Abis.TL.annot.csv"), row.names = FALSE)

cat('Selected by LASSO')
write.csv(get_annotation_by_block(model.char2vec(res.sel$lasso.cv),colnames(X.design),mouselist), 
          paste0(path,"colon cancer/selected.genes.lasso.annot.csv"), row.names = FALSE)

cat('Selected by SCAD')
write.csv(get_annotation_by_block(model.char2vec(res.sel$scad.cv),colnames(X.design),mouselist), 
          paste0(path,"colon cancer/selected.genes.SCAD.annot.csv"), row.names = FALSE)


##### Pseudo posterior probabilities #####

l0.norm.scores = data.frame(apply(res.sel.l0$scores[,-1], MARGIN = 2 , FUN =function(x) exp(-x)/sum(exp(-x))))
l0.norm.scores$model = res.sel.l0$scores[,1]

max.model.size <- max(sapply(str_split(l0.norm.scores$model, ","),length))
models.cols.fmted <- str_split_fixed(l0.norm.scores$model, ",",max.model.size)
models.cols.fmted[models.cols.fmted == ""] <- NA
l0.norm.scores$is.ESM1.in = unlist(apply(models.cols.fmted, MARGIN =1, 
                                         FUN = function(x) '11' %in% x))
l0.norm.scores$is.GAS1.in = unlist(apply(models.cols.fmted, MARGIN =1, 
                                         FUN = function(x) '55' %in% x))
l0.norm.scores$is.HIC1.in = unlist(apply(models.cols.fmted, MARGIN =1, 
                                         FUN = function(x) '68' %in% x))
l0.norm.scores$is.CILP.in = unlist(apply(models.cols.fmted, MARGIN =1, 
                                         FUN = function(x) '155' %in% x))

incl.prob = data.frame(method = c('EBIC', 'S.Abis.TL'), ESM1 = c(NA, NA), HIC1 = c(NA, NA), GAS1 = c(NA, NA), 
                       CILP = c(NA, NA))
incl.prob$ESM1 = l0.norm.scores %>% filter(is.ESM1.in) %>% dplyr::select('EBIC', 'S.Abis.TL') %>% colSums() %>% round(digits = 3)
incl.prob$GAS1 = l0.norm.scores %>% filter(is.GAS1.in) %>% dplyr::select('EBIC', 'S.Abis.TL') %>% colSums() %>% round(digits = 3)
incl.prob$HIC1 = l0.norm.scores %>% filter(is.HIC1.in) %>% dplyr::select('EBIC', 'S.Abis.TL') %>% colSums() %>% round(digits = 3)
incl.prob$CILP = l0.norm.scores %>% filter(is.CILP.in) %>% dplyr::select('EBIC', 'S.Abis.TL') %>% colSums() %>% round(digits = 3)
write.csv(incl.prob, paste0(path,'colon cancer/inclusion.prob.EBIC.S.Abis.TL.csv'), row.names = FALSE)
