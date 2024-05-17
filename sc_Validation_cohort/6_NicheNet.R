
#version of nichenetr package is : 1.1.1

#Differential NicheNet  uses the differential expression between the conditions/niches of the ligand-receptor pairs for prioritization in addition to the ligand activities. The classic NicheNet pipeline on the contrary uses only ligand activity for prioritization (and shows differential expression only in visualizations).

#################################################################################################
#CODE extracted from:
#https://github.com/saeyslab/nichenetr/blob/master/vignettes/differential_nichenet_pEMT.md
#################################################################################################
#The goal of Differential NicheNet is to predict ligand-receptors pairs that are both differentially expressed and active between different conditions  of interest.

# "wilcox", test used for differential expression of L-R pairs between pairwise sender cell-types (or receiver) across niches
# "SCT, for Wilcox DEG only: need to use SCT
# 0.25, recommended for 10x as min_lfc cutoff. We recommend using a cutoff of 0.15 if you have > 2 receiver cells/niches to #compare and use the min_lfc as specificity score. If you have only 2 receivers/niche, we recommend using a higher threshold (such as using 0.25). If you have single-cell data like Smart-seq2 with high sequencing depth, we recommend to also use higher threshold.
# 0.1, typically with Wilcox will be 0.1

####################################################################
## Read in the expression data of interest, and the NicheNet ligand-receptor network and ligand-target matrix
####################################################################
#Load in packages
library(nichenetr)
library(RColorBrewer)
library(tidyverse)
library(Seurat) 

####################################################################
## Data preparation
####################################################################
#choose the receiver cell-types one or several

receiver_v <- c("0-Endothelial", "1-Fibroblast", "4-T cell", "5-Pericyte", "6-T cell", "7-Macrophage", "8-Fibroblast", "10-Endothelial", "11-B Cell", "13-Mast cell", "14-Fibroblast", "16-Endothelial", "17-Fibroblast",  "18-T cell",  "19-Endothelial",  "20-Schwann cell")

for (r in 1:length(receiver_v)){
  
  tag_niche <- paste("12_Luminal_to_",receiver_v[r],sep='')
  niches = list( # user adaptation required on own dataset
  "mHNPC_niche" = list(
   "sender" = c("12_Luminal_mHNPC"),
   "receiver" = paste(receiver_v[r],"_mHNPC",sep="")),
  "LPCa_niche" = list(
   "sender" = c("12_Luminal_LPCa"),
   "receiver" = paste(receiver_v[r],"_LPCa",sep=""))) 
     
  ####################################################################
  ## Read in the expression data
  ####################################################################
  annotated <- readRDS("Annotation_Manual_Validation.rds")

  vector <- as.character(annotated$annotation)
  vector[which(vector =="0-Endothelial")] <- "0_Endothelial"
  vector[which(vector =="1-Fibroblast")] <- "1_Fibroblast"
  vector[which(vector =="4-T cell")] <- "4_T_cell"
  vector[which(vector =="5-Pericyte")] <- "5_Pericyte"
  vector[which(vector =="6-T cell")] <- "6_T_cell"
  vector[which(vector =="7-Macrophage")] <- "7_Macrophage"
  vector[which(vector =="8-Fibroblast")] <- "8_Fibroblast"
  vector[which(vector =="10-Endothelial")] <- "10_Endothelial"
  vector[which(vector =="11-B Cell")] <- "11_B_Cell"
  vector[which(vector =="13-Mast cell")] <- "13_Mast_cell"
  vector[which(vector =="14-Fibroblast")] <- "14_Fibroblast"
  vector[which(vector =="16-Endothelial")] <- "16_Endothelial"
  vector[which(vector =="17-Fibroblast")] <- "17_Fibroblast"
  vector[which(vector =="18-T cell")] <- "18_T_cell"
  vector[which(vector =="19-Endothelial")] <- "19_Endothelial"
  vector[which(vector =="20-Schwann cell")] <- "20_Schwann_cell"
  annotated$annotation <- vector
  
  # I will remove the cells belonging to other categories other than localized and mHNPC
  seurat_obj <- subset(x = annotated, subset = pheno %in% c("LPCa","mHNPC"))
  
  table(seurat_obj$pheno)
  
  DimPlot(seurat_obj, group.by = "annotation") # "celltype" in the tutorial
  DimPlot(seurat_obj, group.by = "pheno") #"pEMT" in the dataset

  # We will now also check the number of cells per cell type condition combination
  table(seurat_obj@meta.data$annotation, seurat_obj@meta.data$pheno) # cell types vs conditions 
  
  # For the Differential NicheNet, we need to compare at least 2 niches or conditions to each other. In this case, the 2 niches are the pEMT-high-niche and the pEMT-low-niche. We will adapt the names of the cell types based on their niche of origin.
  seurat_obj@meta.data$celltype_aggregate = paste(seurat_obj@meta.data$annotation, seurat_obj@meta.data$pheno, sep = "_") # user adaptation required on own dataset
  DimPlot(seurat_obj, group.by = "celltype_aggregate")
  
  seurat_obj@meta.data$celltype_aggregate %>% 
    table() %>% 
    sort(decreasing = TRUE)
  
  celltype_id = "celltype_aggregate" # metadata column name of the cell type of interest
  seurat_obj = SetIdent(seurat_obj, value = seurat_obj[[celltype_id]])
  
  if(sum(unlist(niches) %in% seurat_obj@meta.data$celltype_aggregate)/length(unlist(niches)) != 1) {
    print("check niche names!")
    unlist(niches)[-which(unlist(niches) %in% seurat_obj@meta.data$celltype_aggregate)]
  }
  
  #Read in the NicheNet ligand-receptor network and ligand-target matrix
  ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
  ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns
  
  lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
  
  lr_network = lr_network %>% 
    mutate(bonafide = ! database %in% c("ppi_prediction","ppi_prediction_go"))
  
  lr_network = lr_network %>% 
    dplyr::rename(ligand = from, receptor = to) %>% 
    distinct(ligand, receptor, bonafide)
  head(lr_network)
  
  ####################################################################
  ## Define the niches/microenvironments of interest
  ####################################################################
  #Each niche should have at least one “sender/niche” cell population and one “receiver/target” cell population (present in your expression data)
  #! Important: your receiver cell type should consist of 1 cluster!
  table(seurat_obj@meta.data$celltype_aggregate)
  
  ####################################################################
  ## Calculate differential expression between the niches
  ####################################################################
  # In this step, we will determine DE between the different niches for both senders and receivers to define the DE of L-R pairs.
  
  # Calculate DE; The method to calculate the differential expression is here the standard Seurat Wilcoxon test, but this can be replaced if wanted by the user (only requirement: output tables DE_sender_processed and DE_receiver_processed should be in the same format as shown here).
  #DE will be calculated for each pairwise sender (or receiver) cell type comparison between the niches (so across niches, not within niche). In our case study, this means that DE of myofibroblast_High ligands will be calculated by DE analysis of myofibroblast_High vs myofibroblast_Low; myofibroblast_High vs Endothelial_Low; and myofibroblast_High vs CAF_Low. We split the cells per cell type instead of merging all cells from the other niche to avoid that the DE analysis will be driven by the most abundant cell types.
  
  DE_sender = calculate_niche_de(seurat_obj = seurat_obj %>% 
                                   subset(features = lr_network$ligand %>% 
                                                                      unique()), niches = niches, type = "sender", assay_oi = "SCT") # only ligands important for sender cell types
  DE_receiver = calculate_niche_de(seurat_obj = seurat_obj %>% 
                                     subset(features = lr_network$receptor %>% 
                                                                        unique()), niches = niches, type = "receiver", assay_oi = "SCT") # only receptors now, later on: DE analysis to find targets
  
  DE_sender = DE_sender %>% 
    mutate(avg_log2FC = ifelse(avg_log2FC == Inf, max(avg_log2FC[is.finite(avg_log2FC)]), ifelse(avg_log2FC == -Inf, min(avg_log2FC[is.finite(avg_log2FC)]), avg_log2FC)))
  DE_receiver = DE_receiver %>% 
    mutate(avg_log2FC = ifelse(avg_log2FC == Inf, max(avg_log2FC[is.finite(avg_log2FC)]), ifelse(avg_log2FC == -Inf, min(avg_log2FC[is.finite(avg_log2FC)]), avg_log2FC)))
  
  # process DE results
  DE_sender_processed = process_niche_de(DE_table = DE_sender, niches = niches, expression_pct = 0.1, type = "sender")
  DE_receiver_processed = process_niche_de(DE_table = DE_receiver, niches = niches, expression_pct = 0.1, type = "receiver")
  
  # Combine sender-receiver DE based on L-R pairs:
  #  As mentioned above, DE of ligands from one sender cell type is determined be calculating DE between that cell type, and all the sender cell types of the other niche. To summarize the DE of ligands of that cell type we have several options: we could take the average LFC, but also the minimum LFC compared to the other niche. We recommend using the minimum LFC, because this is the strongest specificity measure of ligand expression, because a high min LFC means that a ligand is more strongly expressed in the cell type of niche 1 compared to all cell types of niche 2 (in contrast to a high average LFC, which does not exclude that one or more cell types in niche 2 also strongly express that ligand).
  specificity_score_LR_pairs = "min_lfc"
  DE_sender_receiver = combine_sender_receiver_de(DE_sender_processed, DE_receiver_processed, lr_network, specificity_score = specificity_score_LR_pairs)
  
  ####################################################################
  ## Calculate ligand activities and infer active ligand-target links
  ####################################################################
  #In this step, we will predict ligand activities of each ligand for each of the receiver cell types across the different niches. This is similar to the ligand activity analysis done in the normal NicheNet pipeline.
  #To calculate ligand activities, we first need to define a geneset of interest for each niche. In this case study, the geneset of interest for the pEMT-high niche are the genes upregulated in pEMT-high tumors compared to pEMT-low tumors, and vice versa.Note that you can also define these geneset of interest in a different way! (eg pathway-based geneset etc). Ligand-target links are inferred in the same way as described in the basic NicheNet vignettes.
  
  specificity_score_targets = "min_lfc"
  
  DE_receiver_targets = calculate_niche_de_targets(seurat_obj = seurat_obj, niches = niches, lfc_cutoff = 0.25, expression_pct = 0.1, assay_oi = "SCT") 
  DE_receiver_processed_targets = process_receiver_target_de(DE_receiver_targets = DE_receiver_targets, niches = niches, expression_pct = 0.1, specificity_score = specificity_score_targets)
  
  background = DE_receiver_processed_targets  %>% 
    pull(target) %>% 
    unique()
  
  geneset_niche1 = DE_receiver_processed_targets %>% 
    filter(receiver == niches[[1]]$receiver & target_score >= 0.25 & target_significant == 1 & target_present == 1) %>% 
    pull(target) %>% 
    unique()
  
  geneset_niche2 = DE_receiver_processed_targets %>% 
    filter(receiver == niches[[2]]$receiver & target_score >= 0.25 & target_significant == 1 & target_present == 1) %>% 
    pull(target) %>% 
    unique()
  
  # Good idea to check which genes will be left out of the ligand activity analysis (=when not present in the rownames of the ligand-target matrix).
  # If many genes are left out, this might point to some issue in the gene naming (eg gene aliases and old gene symbols, bad human-mouse mapping)
  geneset_niche1 %>% setdiff(rownames(ligand_target_matrix))
  geneset_niche2 %>% setdiff(rownames(ligand_target_matrix))
  
  print(c("!!!!!!!!!!!!!!!!!!!!!!!!!!"))
  print("PAY ATTENTION TO THE NUMBER OF GENES INCLUDED IN THE ANALYSIS")
  length(geneset_niche1)
  print(c("We recommend having between 20 and 1000 genes in the geneset of interest"))
  length(geneset_niche2)
  print(c("We recommend having between 20 and 1000 genes in the geneset of interest"))
  length(background)
  print(c("We recommend having background of at least 5000 genes"))
  print(c("!!!!!!!!!!!!!!!!!!!!!!!!!!"))
  
  # We recommend having between 20 and 1000 genes in the geneset of interest, and a background of at least 5000 genes for a proper ligand activity analysis. If you retrieve too many DE genes, it is recommended to use a higher lfc_cutoff threshold. We recommend using a cutoff of 0.15 if you have > 2 receiver cells/niches to compare and use the min_lfc as specificity score. If you have only 2 receivers/niche, we recommend using a higher threshold (such as using 0.25). If you have single-cell data like Smart-seq2 with high sequencing depth, we recommend to also use higher threshold.
  
  top_n_target = 250
  niche_geneset_list = list(
  "mHNPC_niche" = list(
  "receiver" = niches[[1]]$receiver,
  "geneset" = geneset_niche1,
  "background" = background),
  "LPCa_niche" = list(
  "receiver" = niches[[2]]$receiver,
  "geneset" = geneset_niche2 ,
  "background" = background))
  ligand_activities_targets = get_ligand_activities_targets(niche_geneset_list = niche_geneset_list, ligand_target_matrix = ligand_target_matrix, top_n_target = top_n_target)
  
  ####################################################################
  ## Calculate (scaled) expression of ligands, receptors and targets across cell types of interest (log expression values and expression fractions).
  ####################################################################
  # In this step, we will calculate average (scaled) expression, and fraction of expression, of ligands, receptors, and target genes across all cell types of interest. Now this is here demonstrated via the DotPlot function of Seurat, but this can also be done via other ways of course.
  features_oi = union(lr_network$ligand, lr_network$receptor) %>% 
    union(ligand_activities_targets$target) %>% 
    setdiff(NA)
  
  dotplot = suppressWarnings(Seurat::DotPlot(seurat_obj %>% 
                                               subset(idents = niches %>% 
                                                        unlist() %>% 
                                                        unique()), features = features_oi, assay = "SCT"))
  
  exprs_tbl = dotplot$data %>% 
    as_tibble()
  
  exprs_tbl= exprs_tbl %>% 
    dplyr::rename(celltype = id, gene = features.plot, expression = avg.exp, expression_scaled = avg.exp.scaled, fraction = pct.exp)%>% 
    mutate(fraction = fraction/100) %>% 
    as_tibble() %>% 
    select(celltype, gene, expression, expression_scaled, fraction) %>% 
    distinct() %>% 
    arrange(gene) %>% 
    mutate(gene = as.character(gene))
  
  exprs_tbl_ligand = exprs_tbl %>% 
    filter(gene %in% lr_network$ligand) %>% 
    dplyr::rename(sender = celltype, ligand = gene, ligand_expression = expression, ligand_expression_scaled = expression_scaled, ligand_fraction = fraction) 
  
  exprs_tbl_receptor = exprs_tbl %>% 
    filter(gene %in% lr_network$receptor) %>% 
    dplyr::rename(receiver = celltype, receptor = gene, receptor_expression = expression, receptor_expression_scaled = expression_scaled, receptor_fraction = fraction)
  
  exprs_tbl_target = exprs_tbl %>% 
    filter(gene %in% ligand_activities_targets$target) %>% 
    dplyr::rename(receiver = celltype, target = gene, target_expression = expression, target_expression_scaled = expression_scaled, target_fraction = fraction)
  
  exprs_tbl_ligand = exprs_tbl_ligand %>%  
    mutate(scaled_ligand_expression_scaled = scale_quantile_adapted(ligand_expression_scaled)) %>% 
    mutate(ligand_fraction_adapted = ligand_fraction) %>% 
    mutate_cond(ligand_fraction >= 0.1, ligand_fraction_adapted = 0.1)  %>% 
    mutate(scaled_ligand_fraction_adapted = scale_quantile_adapted(ligand_fraction_adapted))
  
  exprs_tbl_receptor = exprs_tbl_receptor %>% 
    mutate(scaled_receptor_expression_scaled = scale_quantile_adapted(receptor_expression_scaled))  %>% 
    mutate(receptor_fraction_adapted = receptor_fraction) %>% 
    mutate_cond(receptor_fraction >= 0.1, receptor_fraction_adapted = 0.1)  %>% 
    mutate(scaled_receptor_fraction_adapted = scale_quantile_adapted(receptor_fraction_adapted))
  
  ####################################################################
  ## Expression fraction and receptor
  ####################################################################
  #In this step, we will score ligand-receptor interactions based on expression strength of the receptor, in such a way that we give higher scores to the most strongly expressed receptor of a certain ligand, in a certain celltype. This will not effect the rank of individual ligands later on, but will help in prioritizing the most important receptors per ligand (next to other factors regarding the receptor - see later).
  exprs_sender_receiver = lr_network %>% 
  inner_join(exprs_tbl_ligand, by = c("ligand")) %>% 
  inner_join(exprs_tbl_receptor, by = c("receptor")) %>% 
    inner_join(DE_sender_receiver %>% 
                 distinct(niche, sender, receiver))
  
  ligand_scaled_receptor_expression_fraction_df = exprs_sender_receiver %>% 
    group_by(ligand, receiver) %>% 
    mutate(rank_receptor_expression = dense_rank(receptor_expression), rank_receptor_fraction = dense_rank(receptor_fraction)) %>% 
    mutate(ligand_scaled_receptor_expression_fraction = 0.5*((rank_receptor_fraction/max(rank_receptor_fraction)) + ((rank_receptor_expression/max(rank_receptor_expression))))) %>% 
    distinct(ligand, receptor, receiver, ligand_scaled_receptor_expression_fraction, bonafide) %>% 
    distinct() %>% 
    ungroup() 
  
  ####################################################################
  ## Prioritization of ligand-receptor and ligand-target links
  ####################################################################
  #In this step, we will combine all the above calculated information to prioritize ligand-receptor-target links. We scale every property of interest between 0 and 1, and the final prioritization score is a weighted sum of the scaled scores of all the properties of interest.
  
  #Note: these settings will give substantially more weight to DE ligand-receptor pairs compared to activity. Users can change this if wanted, just like other settings can be changed if that would be better to tackle the specific biological question you want to address.
  
  prioritizing_weights = c("scaled_ligand_score" = 5,
                       "scaled_ligand_expression_scaled" = 1,
                       "ligand_fraction" = 1,
                       "scaled_ligand_score_spatial" = 0, 
                       "scaled_receptor_score" = 0.5,
                       "scaled_receptor_expression_scaled" = 0.5,
                       "receptor_fraction" = 1, 
                       "ligand_scaled_receptor_expression_fraction" = 1,
                       "scaled_receptor_score_spatial" = 0,
                       "scaled_activity" = 0,
                       "scaled_activity_normalized" = 1,
                       "bona_fide" = 1)
  
  output = list(DE_sender_receiver = DE_sender_receiver, ligand_scaled_receptor_expression_fraction_df = ligand_scaled_receptor_expression_fraction_df, sender_spatial_DE_processed = sender_spatial_DE_processed, receiver_spatial_DE_processed = receiver_spatial_DE_processed, ligand_activities_targets = ligand_activities_targets, DE_receiver_processed_targets = DE_receiver_processed_targets, exprs_tbl_ligand = exprs_tbl_ligand, exprs_tbl_receptor = exprs_tbl_receptor, exprs_tbl_target = exprs_tbl_target)
  prioritization_tables = get_prioritization_tables(output, prioritizing_weights)
  
  prioritization_tables$prioritization_tbl_ligand_receptor %>% 
    filter(receiver == niches[[1]]$receiver) %>% 
    head(10)
  prioritization_tables$prioritization_tbl_ligand_target %>% 
    filter(receiver == niches[[1]]$receiver) %>% 
    head(10)
  prioritization_tables$prioritization_tbl_ligand_receptor %>% 
    filter(receiver == niches[[2]]$receiver) %>% 
    head(10)
  prioritization_tables$prioritization_tbl_ligand_target %>% 
    filter(receiver == niches[[2]]$receiver) %>% 
    head(10)
}

# Differential expression of ligand and expression. 
# Before visualization, we need to define the most important ligand-receptor pairs per niche. We will do this by first determining for which niche the highest score is found for each ligand/ligand-receptor pair. And then getting the top x ligands per niche.
top_ligand_niche_df = prioritization_tables$prioritization_tbl_ligand_receptor %>% 
  select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% 
  group_by(ligand) %>% 
  top_n(1, prioritization_score) %>% 
  ungroup() %>% 
  select(ligand, receptor, niche) %>% 
  dplyr::rename(top_niche = niche)

top_ligand_receptor_niche_df = prioritization_tables$prioritization_tbl_ligand_receptor %>% 
  select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% 
  group_by(ligand, receptor) %>% 
  top_n(1, prioritization_score) %>% 
  ungroup() %>% 
  select(ligand, receptor, niche) %>% 
  dplyr::rename(top_niche = niche)

ligand_prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% 
  select(niche, sender, receiver, ligand, prioritization_score) %>% 
  group_by(ligand, niche) %>% 
  top_n(1, prioritization_score) %>% 
  ungroup() %>% 
  distinct() %>% 
  inner_join(top_ligand_niche_df) %>% 
  filter(niche == top_niche) %>% 
  group_by(niche) %>% 
  top_n(100, prioritization_score) %>% 
  ungroup() # get the top50 ligands per niche



