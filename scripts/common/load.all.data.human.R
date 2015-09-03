## load.all.data.R
## David M. Budden 
## 07/04/2014
##
## from 'Exploring redundancy in mammalian transcriptional regulation by predictive modeling of gene expression'
## load expression, TF-binding, HM and DNase-I hypersensitivity data from human GM12878 data
## 
## www.davidbudden.com/contact

## NOTE: Correct col ranges for output features
## tfas.range = 3:13
## histone.range = 14:20
## all.range = 3:20

load.all.data.human <- function(data.path = '', tfas.source = 'tfas human chipseq.csv') {
  
  ## read in gene mapping
  gene.mapping <- read.csv('../../data/gm12878/ensembl gene mapping.csv')
  colnames(gene.mapping) <- c('Gene','DB_Object_Symbol')
    
  ## read in gene expression data 
  gene.expression <- read.csv(sprintf('%sgene expression human.csv', data.path))
  colnames(gene.expression) <- c('Gene', 'Expression')
   
  ## read in transcription factor association strengths
  tfas.scores <- read.csv(sprintf('%s%s', data.path, tfas.source))
  #colnames(tfas.scores)[1] <- 'entrezgene'  
  
  ## read in histone and DNase scores
  histone.scores <- read.csv(sprintf('%shistone score human.csv', data.path))
  
  ## merge data
  merged.data <- merge(gene.expression, gene.mapping, by='Gene')
  merged.data <- merge(merged.data, tfas.scores, by='Gene')
  merged.data <- merge(merged.data, histone.scores, by='Gene')
}
