
library(WGCNA)

# CHOOSE THE INPUT FILE OF GENE EXPRESSION. Use vst counts, as recommended
mhnpc.lpca <- read.table(paste("vst_matrix.txt",sep=''),header=T)

exp.mat <- (mhnpc.lpca[,grep("vst_",colnames(mhnpc.lpca))])
# make it logarithmic
# exp.mat <- log2(exp.mat+1)
print(dim(mhnpc.lpca))

###################################################################
## from WGCNA https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/
###################################################################
###################################################################
## Data input and cleaning
###################################################################
###################################################################
#  Check individuals with high missingness data
###################################################################
datExpr0 <- t(exp.mat)
rownames(datExpr0) <- colnames(exp.mat)
colnames(datExpr0) <- mhnpc.lpca[,2]

gsg <-  goodSamplesGenes(datExpr0, verbose = 3);
if(gsg$allOK ==TRUE) { print("DATA OK, no excessive missing values, nothing to filter out")}
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes] }

###################################################################
# Next we cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious outliers
###################################################################
sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
pdf(file = paste(tag,"_SampleClustering.pdf",sep=''), width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
dev.off()

###################################################################
#Remove outlier.One can remove it by hand, or use an automatic approach.
#Choose a height cut that will remove the oending sample, say 15 (the red line in the plot), and use a branch cut at
#that height.
###################################################################
# Plot a line to show the cut
#abline(h = 15, col = "red");
# Determine cluster under the line
#clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
#table(clust)
# clust 1 contains the samples we want to keep.
#keepSamples = (clust==1)
#datExpr = datExpr0[keepSamples, ]
#nGenes = ncol(datExpr)
#nSamples = nrow(datExpr)

datExpr <- datExpr0

###################################################################
#Load clinical/covariate trait data
###################################################################
# remove columns that hold information we do not need.
allTraits <- cbind(covariates[,which(colnames(covariates) %in% c("FASTQ.Sample.Name","Groups2","age","RIN","DV200","TumorPurity","PSA","Gleason"))])
colnames(allTraits)[grep("FASTQ.Sample.Name",colnames(allTraits))] <- "sample"
colnames(allTraits)[grep("Groups2",colnames(allTraits))] <- "group"
dim(allTraits)
names(allTraits)

trait.df <- allTraits
if(analysis =="SingleNetwork") {
  
#Transform data. Put group ID as numeric (mHNPC vs LPCa) 
new.v <- rep(1,dim(trait.df)[1])
new.v[which(allTraits$group =="LPCa")] <- 2
}
if(analysis !="SingleNetwork") {  new.v <- rep(1,dim(trait.df)[1]) }
trait.df$group <- new.v
Gleason <- trait.df$Gleason
Gleason <- gsub("3",3,Gleason,fixed=T)
Gleason <- gsub("3+4",7,Gleason,fixed=T)
Gleason <- gsub("4+3",7,Gleason,fixed=T)
Gleason <- gsub("3+3",6,Gleason,fixed=T)
Gleason <- gsub("4+4",8,Gleason,fixed=T)
Gleason <- gsub("4+5",9,Gleason,fixed=T)
Gleason <- gsub("5+4",9,Gleason,fixed=T)
Gleason <- gsub("5+5",10,Gleason,fixed=T)
Gleason[which(Gleason =="7 bilateral")] <- 7
Gleason <- as.numeric(Gleason)
trait.df$Gleason <- Gleason
trait.df$RIN[which(is.na(trait.df$RIN))] <- mean(na.omit(trait.df$RIN))

# Form a data frame analogous to expression data that will hold the clinical traits.
exp.sample.v <- gsub("vst_","",rownames(datExpr));
rownames(datExpr) <- exp.sample.v
traitRows <- match(exp.sample.v, trait.df$sample);
datTraits = trait.df[traitRows,-c(grep("sample",colnames(trait.df))) ]; #definitive reordered trait dataframe
rownames(datTraits) <- trait.df$sample[traitRows];
collectGarbage()

###################################################################
# Plot cluster with covariate info
###################################################################
# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
quant <- datTraits
traitColors = numbers2colors(quant, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
pdf(file = paste(tag,"_SampleClustering_Covariate.pdf",sep=''), width = 12, height = 9);
plotDendroAndColors(sampleTree2, traitColors, groupLabels = names(quant),  main = "Sample dendrogram and trait heatmap")
dev.off()

save(datExpr, datTraits, exp.mat,mhnpc.lpca,file = paste(tag,"_input_WGCNA.RData",sep=''))

###################################################################
#from: https://github.com/mochar/wgcna
#docker pull mochar/wgcna
#docker images
library(WGCNA)
allowWGCNAThreads()


###################################################################
# Network construction and module detection
###################################################################
# my settings
comp <- c("1_mHNPCvsLPCa")  #mHNPC is 1 and Localized is 2
analysis <- "SingleNetwork"
method <- c("joint") #or splitted
filtered.genes <- c("DESeq2_adjusted_filtered")
tag <- paste(analysis,comp,method,filtered.genes,sep="_")

###################################################################
# Load the data saved in the first part
###################################################################
options(stringsAsFactors = FALSE);
enableWGCNAThreads()
#lnames <- load(file = paste(indir, comp,"_input_WGCNA.RData",sep=''));
#lnames
#make sure data looks good
is.numeric(datTraits$group)		
#is.numeric(datTraits$RIN)		
is.numeric(datTraits$age)		
is.numeric(datTraits$DV200)
dim(datExpr) 
#complains about the data not having the correct format:
#checkSets(datExpr, checkStructure = FALSE)
# Get the number of sets in the multiExpr structure.
#nSets <- checkSets(datExpr, checkStructure = TRUE)$nSets
tag <- paste(comp,"_",filtered.genes,"_signed_deepSplit2",sep='')

###################################################################
##Powers analysis
#typical values for signed 11-16 and for unsigned is 6-8
###################################################################
powers = c(seq(2,30,2))
sft=pickSoftThreshold(datExpr,powerVector=powers,verbose = 5, blockSize=30000, networkType = "signed",RsquaredCut = 0.85) 

# now plot:
pdf(paste(tag,"_SoftThresholdingPower.pdf",sep=''),width=6,height=4)
par(mfrow = c(1,2), mar=c(5.1,5.1,4.1,2.1));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n",main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red");
abline(h=0.5,col="red"); abline(h=0.8,col="blue");abline(h=0.9,col="black")
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

###################################################################
#Automatic, one-step network construction and module detection
###################################################################
PWR=sft$powerEstimate
net = blockwiseModules(datExpr,
corType="bicor",
maxBlockSize = 30000,
networkType="signed",
#minCoreKME = 0.4, 
#minKMEtoStay = 0.5,
power=PWR, 
#checkMissingData = TRUE,
#minModuleSize=25,
nThreads=28,
saveTOMs=TRUE,
saveTOMFileBase=tag,
#TOMType = "signed",
#TOMDenom = "mean",
deepSplit=2,
verbose=1,
#mergeCutHeight=0.10,
#reassignThreshold = 1e-10,
numericLabels=TRUE)

str(net)
table(net$color)

###################################################################
#The dendrogram can be displayed together with the color assignment using the following code
###################################################################
# open a graphics window
pdf(paste(tag,"_DendroClusterColors.pdf",sep=''))
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], 
                    "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

###################################################################
#Save the module information necessary for subsequent analysis.
###################################################################
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
pdf(paste(tag,"_barplot_ClusterColors.pdf"))
par(mar=c(10,5,2,2))
barplot(table(moduleColors),las=2,col=names(table(moduleColors)),ylab="Frequency",cex.names=0.5,main=tag)
dev.off()
save.image(file = paste(tag,"_02_networkConstruction-auto_singleblock.RData",sep=''))
