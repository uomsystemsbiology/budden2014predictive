## results1.R
## David M. Budden
## 23/05/2014
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
print('executing results1.R from "Predictive modeling of gene expression from transcriptional regulatory elements"')
print('WARNING: this script may take several hours to complete')

## REPEAT FOR MOUSE AND HUMAN ======================================

for (tfas.source in c('chipseq','pwm','pwm+histone')) {
  for (mouse in c(TRUE,FALSE)) {  
    ## path to data files
    if (mouse)  {data.path <- '../../data/mouse esc eb/'
    } else      {data.path <- '../../data/gm12878/'}
    
    ## path for output files
    if (mouse)  {output.path <- '../../output/results1/mouse/'
    } else      {output.path <- '../../output/results1/gm12878/'}  
    
    ## range of columns in processed data file
    if (mouse) {tfas.range <- 7:18
    } else {tfas.range <- 4:14}
    
    ## LOAD DATA =======================================================
    
    ## load all data files
    print('Loading data... ')
    
    if (mouse)  {data <- load.all.data(data.path, sprintf('tfas mouse %s.csv', tfas.source))
    } else      {data <- load.all.data.human(data.path, sprintf('tfas human %s.csv', tfas.source))}
    
    ## LOG-LINEAR REGRESSION ===========================================
    
    print(sprintf('Evaluating TFAS log-linear regression (%s)...', tfas.source))
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
    write.csv(data.frame(measured=tfas.y, predicted=tfas.y.hat), file=sprintf('%s%s.loglinear.results.csv', output.path, tfas.source), row.names=FALSE)
    
    
    ## SUPPORT VECTOR REGRESSION =======================================
    
    print(sprintf('Evaluating TFAS support vector regression (%s)...', tfas.source))
    ## prepare data
    tfas.prepared.data <- log(normalize.quantiles(as.matrix(data[,tfas.range])) + tfas.sigma)
    ## train model
    tfas.svr.model <- svm(tfas.prepared.data, tfas.y)
    ## predict on training data (for plots)
    tfas.svr.y.hat <- predict(tfas.svr.model, tfas.prepared.data)
    ## 10-fold cross validate for R^2 performance
    tfas.svr.rsq.per.fold <- crossval.svr(tfas.prepared.data, tfas.y, num.folds)
    ## write measured v.s. predicted to file
    write.csv(data.frame(measured=tfas.y, predicted=tfas.svr.y.hat), file=sprintf('%s%s.svr.results.csv', output.path, tfas.source), row.names=FALSE)
    
   
    ## COMPLETE ========================================================
    save(tfas.rsq.per.fld, file=sprintf('%s%s.loglinear.crossval.RData', output.path, tfas.source))
    save(tfas.svr.rsq.per.fold, file=sprintf('%s%s.svr.crossval.RData', output.path, tfas.source))
  }
}
print(sprintf('Complete'))
