#Import Libraries and Sources
library(SummarizedExperiment)
library(dplyr)
source("~/unmc/func4.0.R")

#Import Gene expression Data and Meta Data

#Gene expression data
gene_exp <- input_gene_exp("GSE235874_QC_BQC_Q3Norm_matrixonly.csv")
class(gene_exp)
dim(gene_exp) 

#Meta data
meta_data <- input_meta_data("GSE235874_segmentSummary_readTypes.csv")
class(meta_data)
dim(meta_data) 

#Reorder meta data to the gene exp data order
meta_data <-reorder_metadata(gene_exp,meta_data)

# Create SummarizedExperiment object (logcounts - Normalised Data)
se_object <- SummarizedExperiment(
  assays = list(Q3norm  = gene_exp),  
  colData = meta_data
)


p<-colData(se_object)
p<- assay(se_object, "Q3norm")  
head(assays(se_object)$Q3norm[, 1:5])
se_object[, se_object$Sample == "trt"]


