## generate.all.figures.R
## David M. Budden
## 24/06/2014
##
## from 'Predictive modeling of gene expression from transcriptional regulatory elements'
## run all experiments
## 
## www.davidbudden.com/contact

print('Generating all figures from "Predictive modeling of gene expression from transcriptional regulatory elements"')
setwd('./scripts/figures/')
source('figure1.R')
source('figure2.R')
source('figure3.R')
source('figure4.R')
source('figure5.R')
source('figure6.R')
print('All figures generated!')
setwd('../../')
