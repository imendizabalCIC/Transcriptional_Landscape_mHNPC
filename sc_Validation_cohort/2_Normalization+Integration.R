###################################################################
## Normalization
###################################################################
split_data <- SplitObject(data_combined, split.by = "orig.ident")

##################################################
## Norm_Feature_Scale
##################################################
# normalize and identify variable features for each dataset independently
split_data <- lapply(X = split_data, FUN = function(x) {

  # Normalizing the data
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # Identification of highly variable features (feature selection)
  gene_info_distribution <- summary(Matrix::colSums(x@assays$RNA@counts[,]>0))
  hvg_number <- round(gene_info_distribution[4]+100)

  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = hvg_number)
})

##################################################
## SCTransform
##################################################
for (l in 1:length(split_data)){
  DefaultAssay(split_data[[l]]) <- "SCT"
}

split_data <- lapply(X = split_data, FUN = function(x) {
  x@meta.data$orig.ident[1]
  x <- SCTransform(x, vars.to.regress = c("percent_mt"), vst.flavor = "v2")
  return(x)
})

###################################################################
## Integration
###################################################################
integ_features <- SelectIntegrationFeatures(object.list = split_data, nfeatures = 3000) 

for (l in 1:length(split_data)){
  DefaultAssay(split_data[[l]]) <- "SCT"
}

split_data <- PrepSCTIntegration(object.list = split_data, anchor.features = integ_features)
split_data <- lapply(X = split_data, FUN = RunPCA, features = integ_features)
integ_anchors <- FindIntegrationAnchors(object.list = split_data, normalization.method = "SCT", anchor.features = integ_features, reduction = "CCA")
data_integrated <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT")

####################################################################
## Linear dimensional reduction
####################################################################
DefaultAssay(data_integrated) <- "integrated"

# Run PCA
data_integrated <- RunPCA(object = data_integrated)

# Plot PCA
PCAPlot(data_integrated, group.by = "orig.ident", split.by = "orig.ident")  

# Run UMAP
data_integrated <- RunUMAP(data_integrated, dims = 1:19, reduction = "pca", min.dist = 0.6)

# VizDimLoadings
VizDimLoadings(data_integrated, dims = 1:2, reduction = "pca")

# Plot UMAP
DimPlot(data_integrated, group.by = "orig.ident")

####################################################################
## Clustering
####################################################################
# Plot the elbow plot
ElbowPlot(object = data_integrated, ndims = 60)

##################################################
## Cluster the cells 
##################################################
data_integrated <- FindNeighbors(object = data_integrated, dims = 1:19)

# Determine the clusters                               
data_integrated <- FindClusters(object = data_integrated, resolution = 0.6)
resolution_find_clusters <- "integrated_snn_res.0.6"
Idents(object = data_integrated) <- resolution_find_clusters
data_integrated$seurat_clusters <- data_integrated$integrated_snn_res.0.6

# Plot the UMAP
DimPlot(data_integrated, reduction = "umap", label = TRUE, label.size = 6)