## results3.R
## David M. Budden
## 6/06/2014
##
## from 'Predictive modeling of gene expression from transcriptional regulatory elements'
## generate results for 'Experimental verification'
## 
## www.davidbudden.com/contact

## PACKAGES AND PARAMETERS =========================================

require(preprocessCore)
require(e1071)
require(caret)
library(boot);
library(mixtools);

source('../common/load.all.data.R')
source('../common/load.all.data.human.R')
source('../common/log.linear.model.R')
source('../common/find.best.sigma.R')
source('../common/crossval.log.linear.R')
source('../common/crossval.svr.R')

## folds for cross validation
num.folds <- 10

# introduction
print('executing results3.R from "Predictive modeling of gene expression from transcriptional regulatory elements"')
print('WARNING: this script may take several hours to complete')

## REPEAT FOR MOUSE AND HUMAN ======================================

for (mouse in c(TRUE,FALSE)) {  
  ## path to data files
  if (mouse)  {data.path <- '../../data/mouse esc eb/'
  } else      {data.path <- '../../data/gm12878/'}
  
  ## path for output files
  if (mouse)  {output.path <- '../../output/results3/mouse/'
  } else      {output.path <- '../../output/results3/gm12878/'}  
  
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
  
  print('Evaluating TFAS+Histone+DNase log-linear regression...')
  ## train model
  all.regression <- log.linear.model(data, all.range)
  ##all.regression.model <- all.regression$model
  all.sigma <- all.regression$sigma
  all.data <- all.regression$subset
  ## prepare data (for plots)
  all.y <- log(all.data$Expression + all.sigma)
  ## predict on training data (for plots)
  ##all.y.hat <- as.numeric(fitted(all.regression.model))
  
  ## PRINCIPAL COMPONENT ANALYSIS ==================================

  ## massage TF and HM data into consistent form
  tmp <- all.data
  if (mouse) {
    tmp[,histone.range] <- log(normalize.quantiles(as.matrix(tmp[,histone.range]))+all.sigma)
  } else {
    tmp[,all.range] <- log(normalize.quantiles(as.matrix(tmp[,all.range]))+all.sigma)
  }
  non.pca.data <- cbind(all.y, tmp[,all.range])
  colnames(non.pca.data)[1] <- 'Expression'
  
  ## singular value decomposition
  svd.data <- svd(non.pca.data[,2:length(non.pca.data)])
  u.matrix <- svd.data$u
  sigma.matrix <- diag(svd.data$d)
  v.matrix <- svd.data$v

  ## construct eigengene data for each PC
  pca.data <- non.pca.data
  pca.data[,2:length(pca.data)] <- u.matrix %*% sigma.matrix
  colnames(pca.data)[2:length(pca.data)] <- paste('PC', c(1:(length(pca.data)-1)), sep='')
  
  ## train regression model for each PC and calculate adj. R^2
  pca.rsq <- vector()
  pca.tval <- vector()
  for (i in 2:length(pca.data)) {
    pca.regression.model <- lm(pca.data$Expression ~ pca.data[,i])
    pca.summary <- summary(pca.regression.model)
    pca.rsq[i-1] <- pca.summary$adj.r.squared
    pca.tval[i-1] <- pca.summary$coefficients[2,3]
    
    ## manual extraction of adj.r.squared (required for SVR but can validate here)    
    #pca.y <- pca.data$Expression
    #pca.y.hat <- as.numeric(fitted(pca.regression.model))
    #pca.rsq[i-1] <- print(1-(1-cor(pca.y, pca.y.hat)^2)*((length(pca.y)-1)/(length(pca.y) - 2)))
  }
  names(pca.rsq) <- paste('PC', c(1:(length(pca.data)-1)), sep='')
  names(pca.tval) <- paste('PC', c(1:(length(pca.data)-1)), sep='')

  ## weight loadings by adj. R^2 for each PC
  weighted.loading <- data.frame(tmp=numeric(length(pca.rsq)))
  for (i in 1:length(pca.rsq)) {
    weighted.loading[,i] <- v.matrix[,i]*pca.rsq[i]*sign(pca.tval[i])
  }
  colnames(weighted.loading) <- names(pca.rsq)
  rownames(weighted.loading) <- colnames(non.pca.data[2:length(non.pca.data)])

  ## SUPPORT VECTOR REGRESSION =====================================  

  print('Evaluating TFAS+Histone+DNase SV regression...')

  ## train SVR model for each PC and calculate adj. R^2
  svr.pca.rsq <- vector()
  for (i in 2:length(pca.data)) {
    #svr.pca.regression.model <- lm(pca.data$Expression ~ pca.data[,i])
    pca.y <- pca.data$Expression
    svr.pca.model <- svm(pca.data[,i], pca.y)
    pca.y.hat <- predict(svr.pca.model, pca.data[,i])
    
    svr.pca.rsq[i-1] <- (1-(1-cor(pca.y, pca.y.hat)^2)*((length(pca.y)-1)/(length(pca.y) - 2)))
    ##print(svr.pca.rsq[i-1]) ## DEBUG
  }
  names(svr.pca.rsq) <- paste('PC', c(1:(length(pca.data)-1)), sep='')

  ## weight loadings by adj. R^2 for each PC
  svr.weighted.loading <- data.frame(tmp=numeric(length(svr.pca.rsq)))
  for (i in 1:length(svr.pca.rsq)) {
    svr.weighted.loading[,i] <- v.matrix[,i]*svr.pca.rsq[i]*sign(pca.tval[i]) ## CHECK
  }
  colnames(svr.weighted.loading) <- names(svr.pca.rsq)
  rownames(svr.weighted.loading) <- colnames(non.pca.data[2:length(non.pca.data)])

  ## COMPLETE ========================================================
  
  write.csv(weighted.loading, file=sprintf('%sweighted.loading.loglinear.csv', output.path), row.names=TRUE)
  save(pca.rsq, file=sprintf('%spca.rsq.loglinear.RData', output.path))  
  write.csv(svr.weighted.loading, file=sprintf('%ssvr.weighted.loading.loglinear.csv', output.path), row.names=TRUE)
  save(svr.pca.rsq, file=sprintf('%ssvr.pca.rsq.loglinear.RData', output.path))
  
}
print(sprintf('Complete'))
