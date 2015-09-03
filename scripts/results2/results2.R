## results2.R
## David M. Budden
## 12/05/2014
##
## from 'Predictive modeling of gene expression from transcriptional regulatory elements'
## generate results for 'Prediction accuracy'
## 
## www.davidbudden.com/contact

## PACKAGES AND PARAMETERS =========================================

require(preprocessCore)
require(e1071)
require(caret)

source('../common/load.all.data.R')
source('../common/load.all.data.human.R')
source('../common/log.linear.model.R')
source('../common/find.best.sigma.R')
source('../common/crossval.log.linear.R')
source('../common/crossval.svr.R')

## folds for cross validation
num.folds <- 10

# introduction
print('executing results2.R from "Predictive modeling of gene expression from transcriptional regulatory elements"')
print('WARNING: this script may take several hours to complete')

## REPEAT FOR MOUSE AND HUMAN ======================================

for (mouse in c(TRUE,FALSE)) {  
  ## path to data files
  if (mouse)  {data.path <- '../../data/mouse esc eb/'
  } else      {data.path <- '../../data/gm12878/'}
  
  ## path for output files
  if (mouse)  {output.path <- '../../output/results2/mouse/'
  } else      {output.path <- '../../output/results2/gm12878/'}  
  
  ## range of columns in processed data file
  if (mouse) {
    tfas.range <- 7:18
    histone.range <- 19:26
    all.range <- 7:26
  } else {
    tfas.range <- 4:14
    histone.range <- 15:21
    all.range <- 4:21
  }
  
  ## LOAD DATA =======================================================
  
  ## load all data files
  print('Loading data... ')
  
  if (mouse)  {data <- load.all.data(data.path)
  } else      {data <- load.all.data.human(data.path)}
  
  ## LOG-LINEAR REGRESSION ===========================================
  
  ## TFAS -----------------------------------------
  
  print('Evaluating TFAS log-linear regression...')
  ## train model 
  tfas.regression <- log.linear.model(data, tfas.range)
  tfas.regression.model <- tfas.regression$model
  tfas.sigma <- tfas.regression$sigma
  tfas.data <- tfas.regression$subset
  ## prepare data (for plots)
  tfas.y <- log(tfas.data$Expression + tfas.sigma)
  ## predict on training data (for plots)
  tfas.y.hat <- as.numeric(fitted(tfas.regression.model))
  ## 10-fold cross validate for R^2 performance
  tfas.rsq.per.fld <- as.numeric(crossval.log.linear(tfas.data, tfas.range, tfas.sigma, num.folds))
  ## write measured v.s. predicted to file
  write.csv(data.frame(measured=tfas.y, predicted=tfas.y.hat), file=sprintf('%stfas.loglinear.results.csv', output.path), row.names=FALSE)
  
  ## Histone+DNase --------------------------------
  
  print('Evaluating Histone+DNase log-linear regression...')
  ## train model
  histone.regression <- log.linear.model(data, histone.range)
  histone.regression.model <- histone.regression$model
  histone.sigma <- histone.regression$sigma
  histone.data <- histone.regression$subset
  ## prepare data (for plots)
  histone.y <- log(histone.data$Expression + histone.sigma)
  ## predict on training data (for plots)
  histone.y.hat <- as.numeric(fitted(histone.regression.model))
  ## 10-fold cross validate for R^2 performance
  histone.rsq.per.fld <- as.numeric(crossval.log.linear(histone.data, histone.range, histone.sigma, num.folds))
  ## write measured v.s. predicted to file
  write.csv(data.frame(measured=histone.y, predicted=histone.y.hat), file=sprintf('%shistone.loglinear.results.csv', output.path), row.names=FALSE)
  
  ## TFAS+Histone+DNase ---------------------------
  
  print('Evaluating TFAS+Histone+DNase log-linear regression...')
  ## train model
  all.regression <- log.linear.model(data, all.range)
  all.regression.model <- all.regression$model
  all.sigma <- all.regression$sigma
  all.data <- all.regression$subset
  ## prepare data (for plots)
  all.y <- log(all.data$Expression + all.sigma)
  ## predict on training data (for plots)
  all.y.hat <- as.numeric(fitted(all.regression.model))
  ## 10-fold cross validate for R^2 performance
  all.rsq.per.fld <- as.numeric(crossval.log.linear(all.data, all.range, all.sigma, num.folds))
  ## write measured v.s. predicted to file
  write.csv(data.frame(measured=all.y, predicted=all.y.hat), file=sprintf('%sall.loglinear.results.csv', output.path), row.names=FALSE)
  
  ## SUPPORT VECTOR REGRESSION =======================================
  
  ## TFAS -----------------------------------------
  
  print('Evaluating TFAS SV regression...')
  ## prepare data
  tfas.prepared.data <- log(normalize.quantiles(as.matrix(data[,tfas.range])) + tfas.sigma)
  ## train model
  tfas.svr.model <- svm(tfas.prepared.data, tfas.y)
  ## predict on training data (for plots)
  tfas.svr.y.hat <- predict(tfas.svr.model, tfas.prepared.data)
  ## 10-fold cross validate for R^2 performance
  tfas.svr.rsq.per.fold <- crossval.svr(tfas.prepared.data, tfas.y, num.folds)
  ## write measured v.s. predicted to file
  write.csv(data.frame(measured=tfas.y, predicted=tfas.svr.y.hat), file=sprintf('%stfas.svr.results.csv', output.path), row.names=FALSE)
  
  ## Histone+DNase --------------------------------
  
  print('Evaluating Histone+DNase SV regression...')
  ## prepare data
  histone.prepared.data <- log(normalize.quantiles(as.matrix(data[,histone.range])) + histone.sigma)
  ## train model
  histone.svr.model <- svm(histone.prepared.data, histone.y)
  ## predict on training data (for plots)
  histone.svr.y.hat <- predict(histone.svr.model, histone.prepared.data)
  ## 10-fold cross validate for R^2 performance
  histone.svr.rsq.per.fold <- crossval.svr(histone.prepared.data, histone.y, num.folds)
  ## write measured v.s. predicted to file
  write.csv(data.frame(measured=histone.y, predicted=histone.svr.y.hat), file=sprintf('%shistone.svr.results.csv', output.path), row.names=FALSE)
  
  ## TFAS+Histone+DNase ---------------------------
  
  print('Evaluating TFAS+Histone+DNase SV regression...')
  ## prepare data
  all.prepared.data <- log(normalize.quantiles(as.matrix(data[,all.range])) + all.sigma)
  ## train model
  all.svr.model <- svm(all.prepared.data, all.y)
  ## predict on training data (for plots)
  all.svr.y.hat <- predict(all.svr.model, all.prepared.data)
  ## 10-fold cross validate for R^2 performance
  all.svr.rsq.per.fold <- crossval.svr(all.prepared.data, all.y, num.folds)
  ## write measured v.s. predicted to file
  write.csv(data.frame(measured=all.y, predicted=all.svr.y.hat), file=sprintf('%sall.svr.results.csv', output.path), row.names=FALSE)
  
  ## COMPLETE ========================================================
  save(tfas.rsq.per.fld, histone.rsq.per.fld, all.rsq.per.fld, file=sprintf('%slog.linear.crossval.RData', output.path))
  save(tfas.svr.rsq.per.fold, histone.svr.rsq.per.fold, all.svr.rsq.per.fold, file=sprintf('%ssvr.crossval.RData', output.path))
}
print(sprintf('Complete'))
