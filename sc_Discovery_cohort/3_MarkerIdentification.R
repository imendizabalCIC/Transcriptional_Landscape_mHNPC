DefaultAssay(data_integrated) <- "SCT"
data_integrated <- PrepSCTFindMarkers(data_integrated, assay = "SCT")

###################################################################
## Mannual annotation
###################################################################
data_markers <- FindAllMarkers(data_integrated, slot = "data", assay = "SCT", min.pct = 0.1, logfc.threshold = 0.25)

top15_markers <- as.data.frame(data_markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC))
top100_markers <- as.data.frame(data_markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC))

features_to_plot <- c("CD3D", "CD3E", "CD3G", "PTPRC", 
                      "CD79A", "CD79B", "IGKC", "MS4A1",
                      "CD14", "CD68", "CSF1R", "FCGR3A", "LYZ",
                      "KIT", "MS4A2", "TPSAB1", "TPSB2",
                      "COL1A1", "COL1A2", "COL3A1", "DCN",
                      "RGS5", "ACTA2",
                      "CDH5", "ENG", "PECAM1", "VWF",
                      "AR", "EPCAM", "KRT5", "KRT8", "KRT14",
                      "BIRC5", "CENPF")

DefaultAssay(data_integrated) <- "RNA"

for (feature in features_to_plot){
  FeaturePlot(data_integrated,
              reduction = "umap",
              features = feature,
              order = TRUE,
              min.cutoff = 'q10',
              label = TRUE,
              keep.scale = "all")
  
  VlnPlot(data_integrated, features = feature, pt.size = 0)
}

new.cluster.ids <- c("0" = "0-Luminal",
                     "1" = "1-Luminal", 
                     "2" = "2-Luminal",
                     "3" = "3-Luminal",
                     "4" = "4-Luminal",
                     "5" = "5-Luminal",
                     "6" = "6-Luminal",
                     "7" = "7-Luminal",
                     "8" = "8-Intermediate",
                     "9" = "9-Mast cell",
                     "10" = "10-Macrophage",
                     "11" = "11-T cell",
                     "12" = "12-Luminal",
                     "13" = "13-Endothelial",
                     "14" = "14-Fibroblast",
                     "15" = "15-Luminal",
                     "16" = "16-Endothelial",
                     "17" = "17-T cell", 
                     "18" = "18-Fibroblast", 
                     "19" = "19-B Cell", 
                     "20" = "20-Luminal cycling")

data_annotated <- data_integrated
names(new.cluster.ids) <- levels(data_annotated)
data_annotated <- RenameIdents(data_annotated, new.cluster.ids)
data_annotated$annotation <- data_annotated@active.ident

DotPlot(data_annotated, features = features_to_plot) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))

DimPlot(data_annotated, 
        reduction = "umap", 
        label = TRUE, 
        pt.size = 0.8) + 
  guides(color=guide_legend(ncol =1))