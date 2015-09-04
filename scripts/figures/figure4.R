## figure4.R
## David M. Budden 
## 23/05/2014
##
## from 'Predictive modeling of gene expression from transcriptional regulatory elements'
## generate figure 4
## 
## www.davidbudden.com/contact

## PACKAGES AND PARAMETERS =========================================

require(ggplot2)
require(gridExtra)
require(grid)

## path to data files
data.path <- '../../output/results2/gm12878/'

## path for output files
output.path <- '../../output/figures/'

## colour definitions
black <- rgb(0,0,0)
red <- rgb(1,0,0)

## LOAD DATA =======================================================

tfas.data <- read.csv(sprintf('%stfas.loglinear.results.csv', data.path), row.names=NULL)
tfas.svr.data <- read.csv(sprintf('%stfas.svr.results.csv', data.path), row.names=NULL)
histone.data <- read.csv(sprintf('%shistone.loglinear.results.csv', data.path), row.names=NULL)
histone.svr.data <- read.csv(sprintf('%shistone.svr.results.csv', data.path), row.names=NULL)
all.data <- read.csv(sprintf('%sall.loglinear.results.csv', data.path), row.names=NULL)
all.svr.data <- read.csv(sprintf('%sall.svr.results.csv', data.path), row.names=NULL)

## read adj R^2 performance ---------------------
load(sprintf('%slog.linear.crossval.RData', data.path))
load(sprintf('%ssvr.crossval.RData', data.path))

tfas.log.linear.rsq <- mean(tfas.rsq.per.fld)
tfas.svr.rsq <- mean(tfas.svr.rsq.per.fold)
histone.log.linear.rsq <- mean(histone.rsq.per.fld)
histone.svr.rsq <- mean(histone.svr.rsq.per.fold)
all.log.linear.rsq <- mean(all.rsq.per.fld)
all.svr.rsq <- mean(all.svr.rsq.per.fold)

## plot labels ----------------------------------
tfas.log.linear.label <- paste(sprintf('adj.~R^2 == %.2f', tfas.log.linear.rsq))
tfas.svr.label <- paste(sprintf('adj.~R^2 == %.2f', tfas.svr.rsq))
histone.log.linear.label <- paste(sprintf('adj.~R^2 == %.2f', histone.log.linear.rsq))
histone.svr.label <- paste(sprintf('adj.~R^2 == %.2f', histone.svr.rsq))
all.log.linear.label <- paste(sprintf('adj.~R^2 == %.2f', all.log.linear.rsq))
all.svr.label <- paste(sprintf('adj.~R^2 == %.2f', all.svr.rsq))

## GENERATE FIGURE =================================================

## create individual scatter plots --------------
scatter.plot.tfas.linear <- ggplot(tfas.data, aes(tfas.data$measured, tfas.data$predicted)) +
  geom_point(alpha=0.2, color=black) + 
  geom_smooth(method=lm, color=red, size=1) +
  scale_x_continuous(limits=c(-5,10), breaks=seq(from=-6, to=10, by=2)) +
  scale_y_continuous(limits=c(-8,8), breaks=seq(from=-8, to=8, by=2)) +
  ylab('log(predicted expression)') +
  xlab('log(RNA-seq expression)') +
  theme_bw(base_size=14) +
  annotate('text', x=-4.5, y=7.5, label='(a)', size=6) +
  annotate('text', x=7, y=-7, label=tfas.log.linear.label, parse=TRUE, size=6)

scatter.plot.tfas.svr <- ggplot(tfas.svr.data, aes(tfas.svr.data$measured, tfas.svr.data$predicted)) +
  geom_point(alpha=0.2, color=black) + 
  geom_smooth(method=lm, color=red, size=1) +
  scale_x_continuous(limits=c(-5,10), breaks=seq(from=-6, to=10, by=2)) +
  scale_y_continuous(limits=c(-8,8), breaks=seq(from=-8, to=8, by=2)) +
  ylab('log(predicted expression)') +
  xlab('log(RNA-seq expression)') +
  theme_bw(base_size=14) +
  annotate('text', x=-4.5, y=7.5, label='(d)', size=6) +
  annotate('text', x=7, y=-7, label=tfas.svr.label, parse=TRUE, size=6)

scatter.plot.histone.linear <- ggplot(histone.data, aes(histone.data$measured, histone.data$predicted)) +
  geom_point(alpha=0.2, color=black) + 
  geom_smooth(method=lm, color=red, size=1) +
  scale_x_continuous(limits=c(-5,10), breaks=seq(from=-6, to=10, by=2)) +
  scale_y_continuous(limits=c(-8,8), breaks=seq(from=-8, to=8, by=2)) +
  ylab('log(predicted expression)') +
  xlab('log(RNA-seq expression)') +
  theme_bw(base_size=14) +
  annotate('text', x=-4.5, y=7.5, label='(b)', size=6) +
  annotate('text', x=7, y=-7, label=histone.log.linear.label, parse=TRUE, size=6)

scatter.plot.histone.svr <- ggplot(histone.svr.data, aes(histone.svr.data$measured, histone.svr.data$predicted)) +
  geom_point(alpha=0.2, color=black) + 
  geom_smooth(method=lm, color=red, size=1) +
  scale_x_continuous(limits=c(-5,10), breaks=seq(from=-6, to=10, by=2)) +
  scale_y_continuous(limits=c(-8,8), breaks=seq(from=-8, to=8, by=2)) +
  ylab('log(predicted expression)') +
  xlab('log(RNA-seq expression)') +
  theme_bw(base_size=14) +
  annotate('text', x=-4.5, y=7.5, label='(e)', size=6) +
  annotate('text', x=7, y=-7, label=histone.svr.label, parse=TRUE, size=6)

scatter.plot.all.linear <- ggplot(all.data, aes(all.data$measured, all.data$predicted)) +
  geom_point(alpha=0.2, color=black) + 
  geom_smooth(method=lm, color=red, size=1) +
  scale_x_continuous(limits=c(-5,10), breaks=seq(from=-6, to=10, by=2)) +
  scale_y_continuous(limits=c(-8,8), breaks=seq(from=-8, to=8, by=2)) +
  ylab('log(predicted expression)') +
  xlab('log(RNA-seq expression)') +
  theme_bw(base_size=14) +
  annotate('text', x=-4.5, y=7.5, label='(c)', size=6) +
  annotate('text', x=7, y=-7, label=all.log.linear.label, parse=TRUE, size=6)

scatter.plot.all.svr <- ggplot(all.svr.data, aes(all.svr.data$measured, all.svr.data$predicted)) +
  geom_point(alpha=0.2, color=black) + 
  geom_smooth(method=lm, color=red, size=1) +
  scale_x_continuous(limits=c(-5,10), breaks=seq(from=-6, to=10, by=2)) +
  scale_y_continuous(limits=c(-8,8), breaks=seq(from=-8, to=8, by=2)) +
  ylab('log(predicted expression)') +
  xlab('log(RNA-seq expression)') +
  theme_bw(base_size=14) +
  annotate('text', x=-4.5, y=7.5, label='(f)', size=6) +
  annotate('text', x=7, y=-7, label=all.svr.label, parse=TRUE, size=6)

## arrange plots --------------------------------

## portrait orientation
##grid.arrange(scatter.plot.linear, scatter.plot.svr, ncol=1, nrow=2)
## landscape orientation
#grid.arrange(scatter.plot.linear, scatter.plot.svr, ncol=2, nrow=1)
fig4 <- arrangeGrob(scatter.plot.tfas.linear, scatter.plot.histone.linear, scatter.plot.all.linear, scatter.plot.tfas.svr, scatter.plot.histone.svr, scatter.plot.all.svr, ncol=3, nrow=2)
grid.draw(fig4)

## removing class check from ggsave as a workaround.  This will break at some point in the future
ggsave <- ggplot2::ggsave; body(ggsave) <- body(ggplot2::ggsave)[-2]

## NOTE: export at 6.6 by 10" (convert to 84mm)
ggsave(fig4, file=sprintf('%sfigure4.pdf', output.path), height=6.6, width=14)

## COMPLETE ========================================================

print(sprintf('Complete'))
