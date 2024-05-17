###################################################################
## 1. Prepare the count input file
###################################################################
DefaultAssay(data_annotated) <- "RNA"
exp_rawdata <- as.matrix(GetAssayData(data_annotated, slot="counts")) 

Non_epithelial_clusters <- c(0, 1, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 16, 17, 18, 19, 20)
genome <- "hg20"

###################################################################
## 2. Running copykat
###################################################################
samples <- unique(data_annotated@meta.data$orig.ident)

# Run one sample at a time.
for(s in 1:length(samples)) {
  print(samples[s])
  cur_ind <- which(data_annotated@meta.data$orig.ident == samples[s])
  
  ref_cells <-  which(data_annotated@meta.data$seurat_clusters %in% Non_epithelial_clusters)
  
  cur_ind_ref <- cur_ind[cur_ind %in% ref_cells]
  copykat.test <- copykat(rawmat=exp_rawdata[,cur_ind], id.type="S", cell.line="no", ngene.chr=5, win.size=25, norm.cell.names=colnames(exp_rawdata)[cur_ind_ref], KS.cut=0.1, sam.name=samples[s], distance="euclidean", output.seg="FALSE", plot.genes="TRUE", genome=genome, n.cores=1)
}
