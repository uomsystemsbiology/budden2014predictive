## figure5.R
## David M. Budden 
## 08/06/2014
##
## from 'Predictive modeling of gene expression from transcriptional regulatory elements'
## generate figure 5
## 
## www.davidbudden.com/contact

## PACKAGES AND PARAMETERS =========================================

require(ggplot2)
require(gridExtra)

## path to data files
data.path <- '../../output/results3/mouse/'

## path for output files
output.path <- '../../output/figures/'

## colour definitions
black <- rgb(0,0,0)
red <- rgb(1,0,0)

## LOAD DATA =======================================================

load(sprintf('%spca.rsq.loglinear.RData', data.path))
load(sprintf('%ssvr.pca.rsq.loglinear.RData', data.path))

weighted.loading <- read.csv(sprintf('%sweighted.loading.loglinear.csv', data.path))
svr.weighted.loading <- read.csv(sprintf('%ssvr.weighted.loading.loglinear.csv', data.path))

## GENERATE PLOTS ==================================================

weighted.loading$X <- factor(weighted.loading$X, levels=weighted.loading$X)
## plot weighted loadings for most informative PC (highest adj. R^2)
fig5 <- ggplot(weighted.loading, aes(weighted.loading$X, weighted.loading[,which.max(pca.rsq)+1])) +
  geom_bar(alpha=0.2, color=black, stat='identity') +
  ylab('weighted loading') +
  theme_bw(base_size=14) +
  theme(axis.text.x=element_text(angle=90, size=16, vjust=0.5, hjust=1), axis.title.x=element_blank()) +
  theme(axis.text.y=element_text(size=14), axis.title.y=element_text(size=20))
plot(fig5)

## NOTE: export at 6.6 by 10" (convert to 84mm)
ggsave(fig5, file=sprintf('%sfigure5.pdf', output.path), height=6.6, width=14)

## COMPLETE ========================================================

print(sprintf('Complete'))
