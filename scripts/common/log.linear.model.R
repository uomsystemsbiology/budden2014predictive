## log.linear.model.R
## David M. Budden 
## 07/04/2014
##
## from 'Exploring redundancy in mammalian transcriptional regulation by predictive modeling of gene expression'
## construct log-linear regression model
## 
## www.davidbudden.com/contact

log.linear.model <- function(all.data, data.range) {   
  ## extract 20% held-out data set (for calculating sigma)
  subset.ind <- sample(nrow(all.data), floor(0.2*nrow(all.data)))
  ## NOTE: HELD-ASIDE SET NOT USED TO ALLOW GENE-LEVEL ENRICHMENT ANALYSIS
  #data.subset <- all.data[subset.ind,]  
  #data.remainder <- all.data[-subset.ind,]
  data.subset <- all.data
  data.remainder <- all.data  
  ## define model equation from data frame column names (Y ~ A + B + ...)
  model.features <- colnames(all.data)[data.range]
  equation <- as.formula(paste('Expression ~ ', paste(model.features, collapse='+')))  
  ## fit best sigma value (log(Y + sigma) = ...) to prevent log(0)  
  sigma <- find.best.sigma(data.range, data.subset, equation)    
  ## prepare data for regression
  y <- log(data.remainder$Expression + sigma)
  tmp.remainder <- data.remainder
  tmp.remainder[,data.range] <- normalize.quantiles(as.matrix(tmp.remainder[,data.range]))
  tmp.remainder[,data.range] <- log(tmp.remainder[,c(data.range)] + sigma)  
  prepared.data <- cbind(y, tmp.remainder)
  colnames(prepared.data)[1] <- 'Expression'  
  ## fit and return linear model  
  model <- lm(equation, data=prepared.data)  
  return(list(model=model, sigma=sigma, subset=data.remainder))
}
