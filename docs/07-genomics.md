
# Genomics


## Single Cell Analysis

### Minimising the size of a Seurat Object

Single cell analyses are ofetn collated in Seurat objects which can be huge. We reduced the size of our Seurat object by pulling the relevant sections using the code below reducing the object to <20% of its original size. This works by only including the sparse matrices. This is not a workable example, but the PBMC dataset could be used if necessary. 



```r
#library(Seurat)
#library(dplyr)


#Load in the data

pathtoglobaldata <-"mac_shiny/data/Global_sc_object"
object_global<-readRDS(pathtoglobaldata)


# Get sparse assay data
sizetest_global<-GetAssayData(object = object_global, slot = "data") # counts don’t work, scaled breaks my pc with memory errors, not sure it matters that much
object_global_small<- CreateSeuratObject(sizetest_global)

# Add in classifications. 

object_global_small@reductions$tsne<-object_global@reductions$tsne # reduction embeddings for tsnegraph
object_global_small$Phenotype<- object_global$Phenotype

object_global_small@meta.data$final_classification<- object_global@meta.data$final_classification
object_global_small$final_classificaiton<- object_global$final_classification #cell classifications
Idents(object = object_global_small) <- "final_classificaiton" #set forever


#same
saveRDS(object_global_small,"mac_shiny/data/global_object_sml")
```
