library(glmnet)
library(ggplot2)
library(cowplot)
library(foreach)
library(doMC)

getPopulationOffset = function(y){
  if(!is.factor(y))
    y=factor(y)
  if(length(levels(y))!=2)
    stop("y must be a two-level factor")
  off = sum(y==levels(y)[2])/length(y)
  off = log(off/(1-off))
  return(rep(off,length(y)))
}


#' Do the OvR fit for every variable.  This just does a simple CV selection of regularisation amount.  Far from ideal, but should be good enough for the main conclusions.
multinomialFitCV = function(x,y,nParallel=1,...){
  fits = list()
  if(nParallel>1)
    registerDoMC(cores=nParallel)
  #Do them in order of size
  marks = names(sort(table(as.character(y))))
  for(mark in marks){
    message(sprintf("Fitting model for variable %s",mark))
    fac = factor(y==mark)
    #The two main modes of failure are too few positives and errors constructing lambda.  These should be handled semi-gracefully
    fits[[mark]] = tryCatch(
      cv.glmnet(x,fac,offset=getPopulationOffset(fac),family='binomial',intercept=FALSE,alpha=0.99,nfolds=10,type.measure='class',parallel=nParallel>1,...),
      error = function(e) {
        tryCatch(
          cv.glmnet(x,fac,offset=getPopulationOffset(fac),family='binomial',intercept=FALSE,alpha=0.99,nfolds=10,type.measure='class',parallel=nParallel>1,lambda=exp(seq(-10,-3,length.out=100)),...),
          error = function(e) {
            warning(sprintf("Could not fit model for variable %s",mark))
            return(NULL)
          })
      })
  }
  return(fits)
}

#' Load training data
loadTrainingData = function(dirs){
  mDat = NULL
  toc = NULL
  for(dir in dirs){
    metadata = readRDS(file.path(dir,'metadata.RDS'))
    ttoc = readRDS(file.path(dir,'tableOfCounts.RDS'))
    rownames(metadata$meta.data) = paste0(metadata$meta.data$Sanger_study_ID,'___',gsub('^[0-9]+_','',rownames(metadata$meta.data)))
    colnames(ttoc) = rownames(metadata$meta.data)
    if(is.null(mDat)){
      mDat = metadata$meta.data
    }else{
      mDat = rbind(mDat,metadata$meta.data)
    }
    if(is.null(toc)){
      toc = ttoc
    }else{
      toc = cbind(toc,ttoc)
    }
  }
  return(list(toc=toc,mDat=mDat))
}
