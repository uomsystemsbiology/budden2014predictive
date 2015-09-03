## find.best.sigma.R
## David M. Budden 
## 07/04/2014
##
## from ' '
## calculate sigma value for log-linear regression from held-aside dataset
## 
## www.davidbudden.com/contact

find.best.sigma <- function(data.range, data.subset, formula) {
  rsq.values <- c(0,0)
  for (exponent in -10:10) {
    tmp.subset <- data.subset
    ## perform quantile normalisation and log transformation    
    tmp.subset[,data.range] <- normalize.quantiles(as.matrix(tmp.subset[,data.range]))        
    tmp.subset[,data.range] <- log(tmp.subset[,data.range] + 10^exponent)    
    tmp.subset$Expression <- log(data.subset$Expression+10^exponent)    
    ## perform linear regression to determine R^2 value    
    model <- lm(formula, data=tmp.subset)
    rsq.values <- rbind(rsq.values, c(10^exponent, summary(model)$r.squared))
  }
  ## return sigma value that yields highest R^2 performance
  as.numeric(rsq.values[which(max(rsq.values[,2]) == rsq.values[,2]),1]);
}
