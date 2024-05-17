# Data was dowloaded from GSE141445, from Supplementary file: "GSM4203181_data.raw.matrix.txt.gz".

# Phenotypic data was stored in Supplementary Table 1 (Chen, S., Zhu, G., Yang, Y. et al. Single-cell analysis reveals transcriptomic remodellings in distinct cell types that contribute to human prostate cancer progression. Nat Cell Biol 23, 87–98 (2021). https://doi.org/10.1038/s41556-020-00613-6). 

#Load in packages
library(Seurat) 
library(DoubletFinder) 
library(ggplot2) 

####################################################################
## Set-up the Seurat Object
####################################################################
seurat_obj <- CreateSeuratObject(counts = raw_data_matrix, names.delim = ".", project = "Chen_all")
pheno_data <- read.delim(paste("RAW data/pheno_data.txt", sep ="")) # phenotypic data was stored in the meta.data

# Add a pheno column to metadata with the info of the contrasts for DEGs. 
seurat_obj[["pheno"]] <- NA
seurat_obj$pheno[grep("SC153", seurat_obj$GEO_sample)] <- "P_N"
seurat_obj$pheno[grep("SC154", seurat_obj$GEO_sample)] <- "P"
seurat_obj$pheno[grep("SC155", seurat_obj$GEO_sample)] <- "P"
seurat_obj$pheno[grep("SC156", seurat_obj$GEO_sample)] <- "P"
seurat_obj$pheno[grep("SC159", seurat_obj$GEO_sample)] <- "P"
seurat_obj$pheno[grep("SC162", seurat_obj$GEO_sample)] <- "P_N"
seurat_obj$pheno[grep("SC171", seurat_obj$GEO_sample)] <- "P_B"
seurat_obj$pheno[grep("SC174", seurat_obj$GEO_sample)] <- "P"
seurat_obj$pheno[grep("SC175", seurat_obj$GEO_sample)] <- "P_B"
seurat_obj$pheno[grep("SC176", seurat_obj$GEO_sample)] <- "P"
seurat_obj$pheno[grep("SC177", seurat_obj$GEO_sample)] <- "P_M"
seurat_obj$pheno[grep("SC172", seurat_obj$GEO_sample)] <- "LN_M"
seurat_obj$pheno[grep("SC173", seurat_obj$GEO_sample)] <- "P_T"

####################################################################
## Standard pre-processing workflow - QC and filtering
####################################################################
Idents(seurat_obj) <- seurat_obj$orig.ident

# Mitochondrial Ratio
seurat_obj[["percent_mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
# Ribosomal Ratio
seurat_obj[["percent_rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")

# Calculate the number of genes per UMI for each cell
seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)

# no need of filtering (alredy pre-processed)

##################################################
## DoubletFinder
##################################################
split_data <- SplitObject(seurat_obj, split.by = "orig.ident")
dim_v <- c(15, 13, 18, 12, 17, 13, 14, 15, 16, 12, 17, 15, 14)

for (i in 1:length(split_data)) {
  
  x <- SCTransform(split_data[[i]], vars.to.regress = c("percent_mt"))
  sct_param <- TRUE
  
  x <- RunPCA(x)
  ElbowPlot(x, ndims =  60)
  
  x <- FindNeighbors(object = x, dims = 1:dim_v[i])
  x <- FindClusters(object = x)
  
  x <- RunUMAP(x, dims = 1:dim_v[i], reduction = "pca", min.dist = 0.6)
  sweep_res_x <- paramSweep_v3(x, PCs = 1:dim_v[i], sct = sct_param)
  
  sweep_x <- summarizeSweep(sweep_res_x, GT = FALSE)
  bcmvn_x <- find.pK(sweep_x)
  ggplot(bcmvn_x, aes(pK, BCmetric, group = 1)) +
    geom_point() +
    geom_line()
  
  pK <- bcmvn_x %>%
    filter(BCmetric == max(BCmetric)) %>%
    dplyr::select(pK) 
  pK <- as.numeric(as.character(pK[[1]]))
  
  annotations <- x@meta.data$seurat_clusters
  homotypic_prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.075*nrow(x@meta.data))
  nExp_poi_adj <- round(nExp_poi*(1-homotypic_prop))
  
  x <- doubletFinder_v3(x, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi_adj, reuse.pANN = FALSE, sct = sct_param)
  
  # visualize doublets
  colnames(x@meta.data)
  colname <- colnames(x@meta.data[grep("DF.classifications", colnames(x@meta.data))])
  DimPlot(x, reduction = 'umap', group.by = colname)
}
