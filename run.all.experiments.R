## run.all.experiments.R
## David M. Budden
## 24/06/2014
##
## from 'Predictive modeling of gene expression from transcriptional regulatory elements'
## run all experiments
## 
## www.davidbudden.com/contact

print('Running all experiments from "Predictive modeling of gene expression from transcriptional regulatory elements"')
setwd('./scripts/results1/')
source('results1.R')
setwd('../results2/')
source('results2.R')
setwd('../results3/')
source('results3.R')
print('All experiments complete!')
setwd('../../')
