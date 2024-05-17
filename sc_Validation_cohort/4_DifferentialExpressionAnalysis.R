###################################################################
## Differential expression testing - pseudobulk DESeq2
###################################################################
data_annotated$pheno <- as.factor(data_annotated$pheno)
data_annotated$orig.ident <- as.factor(data_annotated$orig.ident)

data_singlet <- data_annotated[, data_annotated@meta.data[ ,"DF"] == "Singlet"]

# Cluster-sample cell-counts and identify any possible cluster_individual combination with 0 cells
# While DE analysis is typically used for comparison between cell-types, and may struggle with rare subpopulations, DS analysis compares cluster-sample instances that are likely to be much smaller. Thus, DS analysis may only be applicable to more prominent populations. It is therefore recommended to check cluster-sample cell-counts, and to possibly exclude small instances from downstream analyses. In our example, we might consider, for instance, omitting DS analysis on the "Megakaryocytes" and "Dendritic cells" clusters, as these contain less than 30 cells across almost all samples.
table_singlet <- as.data.frame(table(data_singlet$seurat_clusters, data_singlet$orig.ident))
ind_missing <- which(table_singlet$Freq == 0)
cluster_missing <- as.character(table_singlet$Var1[ind_missing])
pheno_missing <- as.character(table_singlet$Var2[ind_missing])

groups <- data_singlet@meta.data[, c("seurat_clusters", "pheno", "orig.ident")]

#try adding a new fake cell to the count matrix
if (length(ind_missing) >= 1) {
  # grep PTEN exactly 
    pten_number <- grep("\\PTEN\\b", rownames(data_singlet))

  fake_cell <- rep(0,dim(data_singlet@assays$SCT@counts)[1])
  for (i in 1:length(pheno_missing)){
    fake_cell[pten_number] <- 1
    if (i == 1){
      matrix_counts <- cbind.data.frame(data_singlet@assays$SCT@counts, fake_cell)
    } else {
      matrix_counts <- cbind.data.frame(matrix_counts, fake_cell)
    }
    colnames(matrix_counts)[dim(matrix_counts)[2]] <-  paste(pheno_missing[i], "_", pheno_missing[i], "_FAKE-", cluster_missing[i],sep='')
    # For aggregation, we use Matrix.utils's aggregate. Matrix function with the colData columns "cluster_id" and "pheno" as groupings. Note that aggregate.Matrix initially yields a matrix of dimensions #(cluster-sample-instances) ? #(genes), which we re-split and transform to obtain, for each cluster, and matrix of dimensions #(genes) ? #(samples):
    # aggregate by cluster-sample
    group_missing <- as.character(groups$pheno[grep(pheno_missing[i],groups$orig.ident)[1]])
    new_row <- c(cluster_missing[i],group_missing,pheno_missing[i])
    groups <- rbind.data.frame(groups,new_row)
    row.names(groups)[dim(groups)[1]] <- paste(colnames(matrix_counts)[dim(matrix_counts)[2]], "_", i, sep ="")
  }
  pb <- aggregate.Matrix(t(matrix_counts),  groupings = groups, fun = "sum")
} else {
  matrix_counts <- data_singlet@assays$SCT@counts
  pb <- aggregate.Matrix(t(data_singlet@assays$SCT@counts), groupings = groups, fun = "sum")
}

table(groups$seurat_clusters, groups$orig.ident) # check if the 0 have been converted to 1

cluster_names <- length(clust <- purrr::set_names(levels(groups$seurat_clusters)))
pheno_names <- length(phen <- purrr::set_names(levels(groups$pheno)))
orig_names <- length(orig <- purrr::set_names(levels(groups$orig.ident)))

# Finally, we compile a table that summarizes the experimental design:
n_cells <- as.numeric(table(groups$orig.ident))
m <- match(orig, groups$orig.ident)
ei <- data.frame(groups[m, ], n_cells, row.names = NULL) %>% select(-"seurat_clusters")

# The approach used here for aggregation presents but one of many ways, and was selected due to its simplicity and efficiency. Alternatively, one could, for example, split the character vector of cells (colnames(sce)) by cluster-pheno, and use rowSums(counts(sce[, ...])) for aggregation. However, this approach is comparatively more verbose and less efficient:
# split cell by cluster-pheno
cs_by_ks <- as.data.frame(groups) %>% 
  rownames_to_column("cells") %>% setDT %>% 
  split(flatten = FALSE, sorted = TRUE,
        by = c("seurat_clusters", "orig.ident")) %>% 
  map_depth(2, "cells")

#lets make a list of clusters, where each element is a matrix of genesx individuals (after all cells are sum into a single value per individual)
cum_list <- vector(mode='list', length=length(clust))
for(clu in 1:length(cs_by_ks)) {
  sub_list <-  cs_by_ks[[clu]]
  clust_list <- vector(mode='list', length=length(orig))
  mat_clust <- as.data.frame(matrix(nrow=dim(matrix_counts)[1],ncol=length(orig)))
  for(ind in 1:length(sub_list)) {
    cluind_cells <- sub_list[[ind]]
    fake_counts <- grep("FAKE",cluind_cells)
    if(length(fake_counts) >0) {     cluind_cells <- gsub("_\\d*$","",cluind_cells[grep("FAKE",cluind_cells)]) }
    mat_clust[,ind] <- apply(as.data.frame(matrix_counts[,match(cluind_cells,colnames(matrix_counts))]),1,sum)
  } 
  colnames(mat_clust) <- orig
  row.names(mat_clust) <- row.names(matrix_counts)
  cum_list[[clu]] <- mat_clust
}
pb2 <- cum_list
names(pb2) <- clust

# Get sample names for each of the cell type clusters
#First, we will create a vector of sample names combined for each of the cell type clusters.
# prep. data.frame for plotting
get_sample_ids <- function(x){  pb2[[x]] %>%    colnames() }
de_samples <- map(1:length(cluster_names), get_sample_ids) %>%  unlist()

#Then we can get the cluster IDs corresponding to each of the samples in the vector.
# Get cluster IDs for each of the samples
samples_list <- map(1:length(clust), get_sample_ids)
get_cluster_ids <- function(x){   rep(names(pb2)[x], each = length(samples_list[[x]])) }
de_cluster_ids <- map(1:length(clust), get_cluster_ids) %>%   unlist()

#Finally, let's create a data frame with the cluster IDs and the corresponding sample IDs. We will merge together the condition information.
# Create a data frame with the sample IDs, cluster IDs and condition
gg_df <- data.frame(cluster_id = de_cluster_ids, sample_id = de_samples)
colnames(gg_df) <- c("seurat_clusters","orig.ident")
gg_df <- left_join(gg_df, ei[, c("orig.ident", "pheno")]) 
metadata <- gg_df %>%   dplyr::select(seurat_clusters, orig.ident, pheno) 

metadata

#########################################################################
#CORRECT the fake count in PTEN
if (length(ind_missing) >= 1) {
    for(m in 1:length(cluster_missing)) {
      pb2[[cluster_missing[m]]][grep("\\PTEN\\b",row.names(pb2[[cluster_missing[m]]])), which(colnames(pb2[[cluster_missing[m]]]) == pheno_missing[m])] <- 0
    }
} 

#########################################################################
#Subsetting dataset to cluster(s) of interest and contrasts of interest
table(metadata$pheno)
de_metadata <- metadata[which(metadata$pheno %in% c("Adeno","NDM")),]
clusters <- levels(as.factor(de_metadata$seurat_clusters))

for(i in 1:length(clusters)) {
  print(paste("Running current cluster:  ", clusters[i], sep=''))
  
  ################################
  #A- subset the cluster of interest
  ################################
  cluster_metadata <- de_metadata[which(de_metadata$seurat_clusters == clusters[i]), ]
  head(cluster_metadata)
  # Assign the rownames of the metadata to be the sample IDs
  rownames(cluster_metadata) <- cluster_metadata$orig.ident
  head(cluster_metadata)
  counts <- pb2[[clusters[i]]]
  cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])
  # Check that all of the row names of the metadata are the same and in the same order as the column names of the counts in order to use as input to DESeq2
  all(rownames(cluster_metadata) == colnames(cluster_counts))  
  
  ################################
  #B-Create DESeq2 object
  ################################
  dim(cluster_counts)
  function_zero <-  function(x){ return(sum(x==0))}
  zeros <- apply(cluster_counts,1,function_zero)
  if(length(zeros)>0) {   cluster_counts <- cluster_counts[-which(zeros == dim(cluster_counts)[2]),] }
  dim(cluster_counts)
  
  cluster_counts <- cluster_counts + 1
  dds <- DESeqDataSetFromMatrix(cluster_counts, colData = cluster_metadata, design = ~ pheno)
  
  # Transform counts for data visualization
  rld <- rlog(dds, blind=TRUE)
  # Plot PCA
  DESeq2::plotPCA(rld, intgroup = "pheno")
  # Hierarchical clustering: Extract the rlog matrix from the object and compute pairwise correlation values
  rld_mat <- assay(rld)
  rld_cor <- cor(rld_mat)
  # Plot heatmap
  pheatmap(rld_cor, annotation = cluster_metadata[, c("pheno"), drop=F])
  
  ################################
  #C: Running DESeq2
  ################################
  dds <- DESeq(dds)
  #We can check the fit of the model to our data by looking at the plot of dispersion estimates.
  plotDispEsts(dds)
  
  #Results
  # /////////////////////////////////////////////////////////////////
  # Output results of Wald test for contrast
  levels(cluster_metadata$pheno)[2]
  levels(cluster_metadata$pheno)[1]
  
  contrast <- c("pheno", levels(cluster_metadata$pheno)[1], levels(cluster_metadata$pheno)[2])
  # /////////////////////////////////////////////////////////////////
  
  resultsNames(dds)
  res <- results(dds,  contrast = contrast, alpha = p_thres)
  
  ################################
  #D- Generate tables
  ################################
  #Table of results for significant genes
  # Turn the results object into a tibble for use with tidyverse functions
  res_tbl <- res %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>%
    as_tibble()
  
  # Check results output
  res_tbl
}
