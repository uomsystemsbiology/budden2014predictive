## figure2.R
## David M. Budden 
## 24/05/2014
##
## from 'Predictive modeling of gene expression from transcriptional regulatory elements'
## generate figure 2
## 
## www.davidbudden.com/contact

## PACKAGES AND PARAMETERS =========================================

require(ggplot2)
require(gridExtra)
require(grid)

## path to data files
data.path <- '../../output/results1/gm12878/'

## path for output files
output.path <- '../../output/figures/'

## colour definitions
black <- rgb(0,0,0)
red <- rgb(1,0,0)

## LOAD DATA =======================================================

#tfas.data <- read.csv(sprintf('%stfas.loglinear.results.csv', data.path), row.names=NULL)
#tfas.svr.data <- read.csv(sprintf('%stfas.svr.results.csv', data.path), row.names=NULL)
#histone.data <- read.csv(sprintf('%shistone.loglinear.results.csv', data.path), row.names=NULL)
#histone.svr.data <- read.csv(sprintf('%shistone.svr.results.csv', data.path), row.names=NULL)
#all.data <- read.csv(sprintf('%sall.loglinear.results.csv', data.path), row.names=NULL)
#all.svr.data <- read.csv(sprintf('%sall.svr.results.csv', data.path), row.names=NULL)

pwm.loglinear.data <- read.csv(sprintf('%spwm.loglinear.results.csv', data.path), row.names=NULL)
pwm.histone.loglinear.data <- read.csv(sprintf('%spwm+histone.loglinear.results.csv', data.path), row.names=NULL)
chipseq.loglinear.data <- read.csv(sprintf('%schipseq.loglinear.results.csv', data.path), row.names=NULL)

pwm.svr.data <- read.csv(sprintf('%spwm.svr.results.csv', data.path), row.names=NULL)
pwm.histone.svr.data <- read.csv(sprintf('%spwm+histone.svr.results.csv', data.path), row.names=NULL)
chipseq.svr.data <- read.csv(sprintf('%schipseq.svr.results.csv', data.path), row.names=NULL)

## read adj R^2 performance ---------------------
#load(sprintf('%slog.linear.crossval.RData', data.path))
#load(sprintf('%ssvr.crossval.RData', data.path))

#tfas.log.linear.rsq <- mean(tfas.rsq.per.fld)
#tfas.svr.rsq <- mean(tfas.svr.rsq.per.fold)
#histone.log.linear.rsq <- mean(histone.rsq.per.fld)
#histone.svr.rsq <- mean(histone.svr.rsq.per.fold)
#all.log.linear.rsq <- mean(all.rsq.per.fld)
#all.svr.rsq <- mean(all.svr.rsq.per.fold)

load(sprintf('%spwm.loglinear.crossval.RData', data.path))
pwm.loglinear.rsq <- mean(tfas.rsq.per.fld)
load(sprintf('%spwm+histone.loglinear.crossval.RData', data.path))
pwm.histone.loglinear.rsq <- mean(tfas.rsq.per.fld)
load(sprintf('%schipseq.loglinear.crossval.RData', data.path))
chipseq.loglinear.rsq <- mean(tfas.rsq.per.fld)

load(sprintf('%spwm.svr.crossval.RData', data.path))
pwm.svr.rsq <- mean(tfas.svr.rsq.per.fold)
load(sprintf('%spwm+histone.svr.crossval.RData', data.path))
pwm.histone.svr.rsq <- mean(tfas.svr.rsq.per.fold)
load(sprintf('%schipseq.svr.crossval.RData', data.path))
chipseq.svr.rsq <- mean(tfas.svr.rsq.per.fold)

## plot labels ----------------------------------
#tfas.log.linear.label <- paste(sprintf('adj.~R^2 == %.2f', tfas.log.linear.rsq))
#tfas.svr.label <- paste(sprintf('adj.~R^2 == %.2f', tfas.svr.rsq))
#histone.log.linear.label <- paste(sprintf('adj.~R^2 == %.2f', histone.log.linear.rsq))
#histone.svr.label <- paste(sprintf('adj.~R^2 == %.2f', histone.svr.rsq))
#all.log.linear.label <- paste(sprintf('adj.~R^2 == %.2f', all.log.linear.rsq))
#all.svr.label <- paste(sprintf('adj.~R^2 == %.2f', all.svr.rsq))

pwm.loglinear.label <- paste(sprintf('adj.~R^2 == %.2f', pwm.loglinear.rsq))
pwm.histone.loglinear.label <- paste(sprintf('adj.~R^2 == %.2f', pwm.histone.loglinear.rsq))
chipseq.loglinear.label <- paste(sprintf('adj.~R^2 == %.2f', chipseq.loglinear.rsq))

pwm.svr.label <- paste(sprintf('adj.~R^2 == %.2f', pwm.svr.rsq))
pwm.histone.svr.label <- paste(sprintf('adj.~R^2 == %.2f', pwm.histone.svr.rsq))
chipseq.svr.label <- paste(sprintf('adj.~R^2 == %.2f', chipseq.svr.rsq))

## GENERATE FIGURE =================================================

## create individual scatter plots --------------
scatter.plot.pwm.loglinear <- ggplot(pwm.loglinear.data, aes(pwm.loglinear.data$measured, pwm.loglinear.data$predicted)) +
  geom_point(alpha=0.2, color=black) + 
  geom_smooth(method=lm, color=red, size=1) +
  scale_x_continuous(limits=c(-8,10), breaks=seq(from=-8, to=10, by=2)) +
  scale_y_continuous(limits=c(-20,8), breaks=seq(from=-20, to=8, by=2)) +
  ylab('log(predicted expression)') +
  xlab('log(RNA-seq expression)') +
  theme_bw(base_size=14) +
  annotate('text', x=-7.5, y=7.5, label='(a)', size=6) +
  annotate('text', x=7, y=-7, label=pwm.loglinear.label, parse=TRUE, size=6)

scatter.plot.pwm.histone.loglinear <- ggplot(pwm.histone.loglinear.data, aes(pwm.histone.loglinear.data$measured, pwm.histone.loglinear.data$predicted)) +
  geom_point(alpha=0.2, color=black) + 
  geom_smooth(method=lm, color=red, size=1) +
  scale_x_continuous(limits=c(-3,10), breaks=seq(from=-2, to=10, by=2)) +
  scale_y_continuous(limits=c(-4,4), breaks=seq(from=-4, to=4, by=2)) +
  ylab('log(predicted expression)') +
  xlab('log(RNA-seq expression)') +
  theme_bw(base_size=14) +
  annotate('text', x=-2.5, y=3.5, label='(b)', size=6) +
  annotate('text', x=7, y=-0, label=pwm.histone.loglinear.label, parse=TRUE, size=6)

scatter.plot.chipseq.loglinear <- ggplot(chipseq.loglinear.data, aes(chipseq.loglinear.data$measured, chipseq.loglinear.data$predicted)) +
  geom_point(alpha=0.2, color=black) + 
  geom_smooth(method=lm, color=red, size=1) +
  scale_x_continuous(limits=c(-5,10), breaks=seq(from=-6, to=10, by=2)) +
  scale_y_continuous(limits=c(-6,4), breaks=seq(from=-6, to=4, by=2)) +
  ylab('log(predicted expression)') +
  xlab('log(RNA-seq expression)') +
  theme_bw(base_size=14) +
  annotate('text', x=-4.5, y=3.5, label='(c)', size=6) +
  annotate('text', x=7, y=-3, label=chipseq.loglinear.label, parse=TRUE, size=6)

scatter.plot.pwm.svr <- ggplot(pwm.svr.data, aes(pwm.svr.data$measured, pwm.svr.data$predicted)) +
  geom_point(alpha=0.2, color=black) + 
  geom_smooth(method=lm, color=red, size=1) +
  scale_x_continuous(limits=c(-8,10), breaks=seq(from=-8, to=10, by=2)) +
  scale_y_continuous(limits=c(-22,8), breaks=seq(from=-22, to=8, by=2)) +
  ylab('log(predicted expression)') +
  xlab('log(RNA-seq expression)') +
  theme_bw(base_size=14) +
  annotate('text', x=-7.5, y=7.5, label='(d)', size=6) +
  annotate('text', x=7, y=-7, label=pwm.svr.label, parse=TRUE, size=6)

scatter.plot.pwm.histone.svr <- ggplot(pwm.histone.svr.data, aes(pwm.histone.svr.data$measured, pwm.histone.svr.data$predicted)) +
  geom_point(alpha=0.2, color=black) + 
  geom_smooth(method=lm, color=red, size=1) +
  scale_x_continuous(limits=c(-3,10), breaks=seq(from=-2, to=10, by=2)) +
  scale_y_continuous(limits=c(-4,8), breaks=seq(from=-4, to=8, by=2)) +
  ylab('log(predicted expression)') +
  xlab('log(RNA-seq expression)') +
  theme_bw(base_size=14) +
  annotate('text', x=-2.5, y=7.5, label='(e)', size=6) +
  annotate('text', x=7, y=-0, label=pwm.histone.svr.label, parse=TRUE, size=6)

scatter.plot.chipseq.svr <- ggplot(chipseq.svr.data, aes(chipseq.svr.data$measured, chipseq.svr.data$predicted)) +
  geom_point(alpha=0.2, color=black) + 
  geom_smooth(method=lm, color=red, size=1) +
  scale_x_continuous(limits=c(-5,10), breaks=seq(from=-6, to=10, by=2)) +
  scale_y_continuous(limits=c(-6,8), breaks=seq(from=-6, to=8, by=2)) +
  ylab('log(predicted expression)') +
  xlab('log(RNA-seq expression)') +
  theme_bw(base_size=14) +
  annotate('text', x=-4.5, y=7.5, label='(f)', size=6) +
  annotate('text', x=7, y=-3, label=chipseq.svr.label, parse=TRUE, size=6)

## arrange plots --------------------------------

fig2 <- arrangeGrob(scatter.plot.pwm.loglinear, scatter.plot.pwm.histone.loglinear, scatter.plot.chipseq.loglinear, scatter.plot.pwm.svr, scatter.plot.pwm.histone.svr, scatter.plot.chipseq.svr, ncol=3, nrow=2)
grid.draw(fig2)

## removing class check from ggsave as a workaround.  This will break at some point in the future
ggsave <- ggplot2::ggsave; body(ggsave) <- body(ggplot2::ggsave)[-2]

## NOTE: export at 6.6 by 10" (convert to 84mm)
ggsave(fig2, file=sprintf('%sfigure2.pdf', output.path), height=6.6, width=14)

## COMPLETE ========================================================

print(sprintf('Complete'))
