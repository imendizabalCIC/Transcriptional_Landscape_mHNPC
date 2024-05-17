###################################################################
# Load RAW count file
###################################################################
raw.df <- read.delim(counts.file)
rawCounts <- as.matrix(raw.df)

# Ordenamos coldata igual que la matriz de datos
rawCounts <- rawCounts[,which(colnames(rawCounts) %in% row.names(coldata))]
coldata = coldata[(match(colnames(rawCounts), coldata$Name)),]
coldata$Group <- factor(coldata$Group)
if(sum(row.names(coldata) == colnames(rawCounts)) != dim(rawCounts)[2]) { print("ERROR:REORDER COUNT MATRIX!!")}

###################################################################
# Prefiltering for low count genes
###################################################################
rawCounts <- rawCounts[rowSums(rawCounts >= 10) >= min(table(coldata$Group)),]
print(paste("Total individuals:",dim(coldata)[1]))
print(paste("Total Genes:",nrow(rawCounts)))

###################################################################
# Analysis.
###################################################################
coldata$Age <- as.numeric((as.numeric(coldata$Age)))
coldata$DV200 <-as.numeric((as.numeric(coldata$DV200)))
coldata$Group <- as.factor(coldata$Group)

dds <- DESeqDataSetFromMatrix(countData = rawCounts, colData = coldata, design = ~Age+DV200+Group) 
resultsNames(dds)

anotacion <- read.table(paste(indir,length.file,sep=''), header=T, sep="\t", stringsAsFactors=FALSE)
anotacionFiltrada <- anotacion[anotacion$GeneID %in% rownames(counts(dds)), ]
geneNames <- data.frame(GeneID=anotacionFiltrada$GeneID, Symbol=anotacionFiltrada$gene_name)

###################################################################
# STABILIZING COUNT VARIANCE
# For large datasets, apply VST transformation. It does not use the design to remove variation in the data.
dds_transformed <- vst(dds, blind=var.blind) 

###################################################################
# DEGs
###################################################################
dds <- DESeq(dds)
resultsNames(dds) 

#DEGs: EXTRACT info per contrast
mod_mat <- model.matrix(design(dds), colData(dds))
contrast.ndm <- colMeans(mod_mat[dds$Group == "mHNPC", ])
contrast.loc <- colMeans(mod_mat[dds$Group %in% c("LPCa"),])

resDESeq2 <- results(dds, contrast = contrast.ndm - contrast.loc,alpha=0.05)
resDESeq2.shrunken <- lfcShrink(dds, res=resDESeq2,contrast= contrast.ndm - contrast.loc,type="ashr")

#visualize the DEG results
  mcols(resDESeq2, use.names = TRUE)
  summary(resDESeq2)
  table(resDESeq2$padj < 0.05)
