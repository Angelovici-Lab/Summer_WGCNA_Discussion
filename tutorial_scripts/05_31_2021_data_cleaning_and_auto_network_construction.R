rm(list=ls())

library(jpeg)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)

library(foreach)
library(iterators)
library(parallel)
library(doParallel)

library(WGCNA)
# library(KEGGREST)
# library(biomaRt)

set.seed(1)

# Enable WGCNA threads to speed up calculations
enableWGCNAThreads()


##################################################
# Constants/Variables
##################################################
softPower <- 6
minModuleSize <- 30

# Eigengenes clustering tree cutting threshold
MEDissThres <- 0.25


##################################################
# Output folder
##################################################
output_path <- file.path("/home/ycth8/data/projects/05_30_2021_summer_WGCNA/tutorial_output/05_31_2021_data_cleaning_and_auto_network_construction")

if(!dir.exists(output_path)){
  dir.create(output_path, showWarnings=FALSE, recursive=TRUE)
  if(!dir.exists(output_path)){
    quit(status=1)
  }
}


##################################################
# Read in input file
##################################################

folder_path = file.path("/home/ycth8/data/projects/05_30_2021_summer_WGCNA/tutorial_data")

femData = read.csv(
  file = file.path(folder_path, "LiverFemale3600.csv"),
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if(ncol(femData) > 6) { print(head(femData[, 1:6])) } else { print(head(femData)) }
if(ncol(femData) > 6) { print(tail(femData[, 1:6])) } else { print(tail(femData)) }
print(dim(femData))


##################################################
# Pre-process input file
##################################################

datExpr0 = as.data.frame(t(femData[, -c(2:8)]), stringsAsFactors = FALSE)
colnames(datExpr0) = datExpr0[1,]
datExpr0 <- datExpr0[-1,]
for (i in 1:ncol(datExpr0)) {
  datExpr0[,i] = as.numeric(as.character(datExpr0[,i]))
}


cat(rep("\n", 2))
if(ncol(datExpr0) > 6) { print(head(datExpr0[, 1:6])) } else { print(head(datExpr0)) }
if(ncol(datExpr0) > 6) { print(tail(datExpr0[, 1:6])) } else { print(tail(datExpr0)) }
print(dim(datExpr0))


##################################################
# Quality check
##################################################

gsg = goodSamplesGenes(datExpr0, verbose = 3)

cat(rep("\n", 2))
print(gsg$allOK)

# If gsg$allOK is FALSE, remove genes and samples
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}


##################################################
# Cluster the data
##################################################

sampleTree = hclust(dist(datExpr0), method = "average")

# Plot tree
cat(rep("\n", 2))
jpeg(filename = file.path(output_path, "sampleTree.jpeg"),  width = 1920, height = 720)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
abline(h = 15.2, col = "red")
dev.off()

# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 15.2, minSize = 10)

print(table(clust))

# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# Collect garbage
collectGarbage()


##################################################
# Save processed female expression data
##################################################
write.csv(
  x = datExpr,
  file = file.path(output_path, "..", "fem_datExpr.csv"),
  na = "",
  quote = FALSE
)


#######################################################################
## Performs principal components analysis (PCA)
#######################################################################
pca_df <- datExpr
for (i in 1:ncol(pca_df)) {
  pca_df[is.na(pca_df[,i]),i] = mean(pca_df[,i], na.rm = TRUE)
}
pca <- prcomp(pca_df, center = TRUE, scale. = TRUE)

pca_variance <- pca$sdev^2
pca_variance_percentage <- round(pca_variance/sum(pca_variance, na.rm = TRUE)*100, digits=2)

x <- pca$x %>% as.data.frame(stringsAsFactors=TRUE)

## Plot PC1 vs PC2 with GGplot2
p <- ggplot(x, aes(x=PC1, y=PC2)) +
  geom_text(label=rownames(x)) +
  labs(
    x=paste0("PC1 - ", pca_variance_percentage[1], "%"),
    y=paste0("PC2 - ", pca_variance_percentage[2], "%")
  )

ggsave(
  filename="PCA_Plot.jpg",
  plot=p,
  path=output_path
)


##################################################
# Choose a set of soft-thresholding powers
##################################################

powers = 1:60
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot tree
cat(rep("\n", 2))
jpeg(filename = file.path(output_path, "softThreshold.jpeg"), width = 1920, height = 480)
par(mfrow = c(1,2))
cex1 = 0.9

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=sft$fitIndices[sft$fitIndices$Power == softPower, "SFT.R.sq"], col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
abline(h=sft$fitIndices[sft$fitIndices$Power == softPower, "mean.k."], col="red")
dev.off()


# Collect garbage
collectGarbage()


##################################################
# Create blockwise modules
##################################################

net = blockwiseModules(
  datExpr,
  power = softPower,
  TOMType = "unsigned",
  minModuleSize = minModuleSize,
  reassignThreshold = 0,
  mergeCutHeight = MEDissThres,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = file.path(output_path, "femaleMouseTOM"),
  verbose = 3
)

# The numbers on top are the color labels
print(table(net$colors))


##################################################
# Plot the dendrogram and colors underneath
##################################################

# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)

cat(rep("\n", 2))
jpeg(filename = file.path(output_path, "networkConstruction_merged_colors.jpeg"), width = 1920, height = 720)
plotDendroAndColors(
  net$dendrograms[[1]],
  mergedColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)
dev.off()


##################################################
# Save important variables as RData
##################################################

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(
  MEs, moduleLabels, moduleColors, geneTree,
  file = file.path(output_path, "FemaleLiver-02-networkConstruction-auto.RData")
)
