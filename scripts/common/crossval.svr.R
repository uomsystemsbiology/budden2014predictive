## crossval.svr.R
## David M. Budden
## 07/04/2014
##
## from 'Exploring redundancy in mammalian transcriptional regulation by predictive modeling of gene expression'
## perform k-fold cross validation of support vector regression model
## 
## www.davidbudden.com/contact

## TO-DO: Track down bug with parallelising SVR model (doesn't affect results)

## PACKAGES AND PARAMETERS =========================================

#require(parallel)
#require(foreach)
#require(doParallel)

## seed random number generators
seed <- 1337

## CROSS VALIDATE MODEL ============================================

crossval.svr <- function(prepared.data, y, num.folds) {
  ## create folds
  set.seed(seed)
  flds <- createFolds(data$Expression, k=num.folds)
  
  ## set up parallel processing
  #cl <- makeCluster(detectCores())
  #registerDoParallel(cl)
  
  rsq.per.fld <- numeric(num.folds)  
  #rsq.per.fld <- foreach(i=1:num.folds, .combine=rbind, .packages=c('preprocessCore','e1071','caret')) %dopar% {  
  for (i in 1:num.folds) {
    ## training data
    training.flds <- setdiff(1:nrow(data), flds[[i]])
    training.data <- prepared.data[training.flds,]
    training.y <- y[training.flds]
    ## validation data
    validation.flds <- flds[[i]]
    validation.data <- prepared.data[validation.flds,]
    validation.y <- y[validation.flds]
    ## train on k-1 folds
    svr.model <- svm(training.data, training.y)
    ## evaluate on remaining fold 
    y.hat <- predict(svr.model, validation.data)
    ## calculate adjusted R^2 from PCC
    rsq.per.fld[i] <- 1-(1-cor(validation.y, y.hat)^2)*((length(validation.y)-1)/(length(validation.y) - ncol(prepared.data) - 1))
  }
  
  ## stop parallel processing
  #stopCluster(cl)
  ## return R^2 values for each fold
  return(rsq.per.fld)
}