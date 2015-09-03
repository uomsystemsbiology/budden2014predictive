## load.all.data.R
## David M. Budden 
## 07/04/2014
##
## from 'Exploring redundancy in mammalian transcriptional regulation by predictive modeling of gene expression'
## load expression, TF-binding, HM and DNase-I hypersensitivity data from Ouyang09 and McLeay12
## 
## www.davidbudden.com/contact

## NOTE: Correct col ranges for output features
## tfas.range = 7:18
## histone.range = 19:26
## all.range = 7:26

load.all.data <- function(data.path = '', tfas.source = 'tfas mouse chipseq.csv') {
  ## read in entrez -> ensemble -> mgi gene ID mapping
  gene.mapping <- read.csv(sprintf('%sgenerated gene mapping.csv', data.path))
  gene.mapping[,1] <- NULL  
  ## read in gene expression data (Ouyang09, McLeay12)  
  gene.expression <- read.csv(sprintf('%sgene expression mouse.csv', data.path))
  colnames(gene.expression)[1] <- 'entrezgene'  
  ## read in transcription factor association strengths
  tfas.scores <- read.csv(sprintf('%s%s', data.path, tfas.source))
  colnames(tfas.scores)[1] <- 'entrezgene'  
  ## read in histone and DNase scores
  histone.scores <- read.csv(sprintf('%shistone score mouse.csv', data.path))
  colnames(histone.scores)[1] <- 'ensembl_gene_id'
  ## merge and sort
  merged.data <- merge(gene.mapping, gene.expression, by='entrezgene')
  merged.data <- merge(merged.data, tfas.scores, by=c('entrezgene','Gene'))
  merged.data <- merge(merged.data, histone.scores, by='ensembl_gene_id')
  merged.data <- merged.data[with(merged.data, order(Gene, entrezgene, ensembl_gene_id, mgi_id)),]
  colnames(merged.data)[5] <- 'Expression'
  ## return merged data
  return(merged.data)
}
