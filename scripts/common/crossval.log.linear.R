## crossval.log.linear.R
## David M. Budden
## 07/04/2014
##
## from 'Exploring redundancy in mammalian transcriptional regulation by predictive modeling of gene expression'
## perform k-fold cross validation of log-linear regression model
## 
## www.davidbudden.com/contact

## PACKAGES AND PARAMETERS =========================================

require(parallel)
require(foreach)
require(doParallel)

## seed random number generators
seed <- 1337
## number of cores for parallel processing
num.cores <- detectCores()

## CROSS VALIDATE MODEL ============================================

crossval.log.linear <- function(data, data.range, sigma, num.folds) {  
  ## create folds
  set.seed(seed)
  flds <- createFolds(data$Expression, k=num.folds)  
  ## set up parallel processing
  cl <- makeCluster(num.cores)
  registerDoParallel(cl)
  ## calculate R^2 for each fold
  rsq.per.fld <- foreach(i=1:num.folds, .combine=rbind, .export=c('log.linear.model', 'find.best.sigma'), .packages=c('preprocessCore')) %dopar% {
    ## train on k-1 folds
    training.flds <- setdiff(1:nrow(data), flds[[i]])
    validation <- log.linear.model(data[training.flds,], data.range)
    validation.model <- validation$model
    validation.sigma <- validation$sigma      
    ## log transform and quantile normalise data for validation
    prepared.data <- data[flds[[i]], data.range]
    prepared.data <- normalize.quantiles(as.matrix(prepared.data))
    prepared.data <- data.frame(log(prepared.data + sigma))
    prepared.data <- cbind(data$Expression[flds[[i]]], prepared.data)
    colnames(prepared.data) <- c('Expression', colnames(data)[data.range])    
    ## evaluate on remaining fold  
    validation.y <- log(data$Expression[flds[[i]]] + validation.sigma)
    validation.y.hat <- as.numeric(predict(validation.model, prepared.data))    
    ## calculate adjusted R^2 from PCC
    1-(1-cor(validation.y, validation.y.hat)^2)*((length(validation.y)-1)/(length(validation.y)-ncol(prepared.data)-1))
  }  
  ## stop parallel processing
  stopCluster(cl)
  ## return R^2 values for each fold
  return(rsq.per.fld)
}
