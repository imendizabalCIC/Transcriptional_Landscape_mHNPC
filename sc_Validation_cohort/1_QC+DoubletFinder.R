samples_v <- c("S_03", "S_04", "S_05", "S_06")

####################################################################
## Set-up the Seurat Object
####################################################################
for (sample in samples_v){
  path_to_sample <- paste("CellRangerCount", sample, "outs/filtered_feature_bc_matrix/", sep = "/")
  
  seurat_data <- Read10X(data.dir = path_to_sample) 
  
  seurat_obj <- CreateSeuratObject(counts = seurat_data, project = sample, min.cells = 0, min.features = 0)
  assign(sample, seurat_obj)
}

if (length(samples_v) > 1){
  seurat_obj <- merge(x = get(samples_v[1]), y = sapply(samples_v[-1], get), add.cell.ids = samples_v)
}

seurat_obj[["pheno"]] <- NA
seurat_obj$pheno[which(seurat_obj$orig.ident == "S_03")] <- "LPCa"
seurat_obj$pheno[which(seurat_obj$orig.ident == "S_04")] <- "LPCa"
seurat_obj$pheno[which(seurat_obj$orig.ident == "S_05")] <- "mHNPC"
seurat_obj$pheno[which(seurat_obj$orig.ident == "S_06")] <- "mHNPC"

####################################################################
## Standard pre-processing workflow - QC and filtering
####################################################################
# Mitochondrial Ratio
seurat_obj[["percent_mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
# Ribosomal Ratio
seurat_obj[["percent_rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")

# Calculate the number of genes per UMI for each cell
seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)

# Filtering
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & log10GenesPerUMI > 0.8 & percent_mt < 20)

counts <- GetAssayData(object = seurat_obj, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
seurat_obj <- seurat_obj[keep_genes, ]

##################################################
## DoubletFinder
##################################################
split_data <- SplitObject(seurat_obj, split.by = "orig.ident")
dim_v <- c(14, 17, 16, 13)

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
