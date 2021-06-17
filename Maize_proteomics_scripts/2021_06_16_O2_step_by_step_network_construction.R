rm(list=ls())

library(jpeg)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)

library(argparse)

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


#######################################################################
# Constants/Variables
#######################################################################
selected_genotype <- "O2"
softPower <- 22
minModuleSize <- 5
mergeCutHeight <- 0.0000001


##################################################
# Output folder
##################################################
output_path <- file.path(
  paste0(
    "/home/ycth8/data/projects/2021_05_30_summer_WGCNA/Maize_proteomics_output/",
    paste0("2021_06_10_", selected_genotype, "_step_by_step_network_construction")
  )
)

if(!dir.exists(output_path)){
  dir.create(output_path, showWarnings=FALSE, recursive=TRUE)
  if(!dir.exists(output_path)){
    quit(status=1)
  }
}


##################################################
# Read in input file
##################################################

folder_path = file.path("/home/ycth8/data/projects/2021_05_30_summer_WGCNA/Maize_proteomics_output/")

datExpr = read.csv(
  file = file.path(folder_path, "datExpr.csv"),
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

datExpr = datExpr[startsWith(rownames(datExpr), selected_genotype),]


##################################################
# Quality check
##################################################

gsg = goodSamplesGenes(datExpr, verbose = 3)

cat(rep("\n", 2))
print(gsg$allOK)

# If gsg$allOK is FALSE, remove genes and samples
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")))

  # Remove the offending genes and samples from the data:
  datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
}


##################################################
# Choose a set of soft-thresholding powers
##################################################

powers = 1:40
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
# Make a gene tree and identify modules
##################################################
adjacency = adjacency(datExpr, power = softPower)

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average")

# Plot the resulting clustering tree (dendrogram)
cat(rep("\n", 2))
jpeg(filename = file.path(output_path, "geneClustering.jpeg"), width = 1920, height = 720)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()

jpeg(filename = file.path(output_path, "geneClustering_with_gene_id.jpeg"), width = 1920, height = 720)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = colnames(datExpr), hang = 0.04)
dev.off()


# We like large modules, so we set the minimum module size relatively high
# Check minModuleSize in the constants/variables section

# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(
  dendro = geneTree,
  distM = dissTOM,
  deepSplit = 2,
  pamRespectsDendro = FALSE,
  minClusterSize = minModuleSize
)
print(table(dynamicMods))


# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)

print(table(dynamicColors))

# Plot the dendrogram and colors underneath
cat(rep("\n", 2))
jpeg(filename = file.path(output_path, "networkConstruction_dynamic_colors.jpeg"), width = 1920, height = 720)
plotDendroAndColors(
  geneTree,
  dynamicColors,
  "Dynamic Tree Cut",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Gene dendrogram and module colors"
)
dev.off()


##################################################
# Modules merging by clustering module eigengenes
##################################################
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes

# Create expression table
expr <- MEs %>%
  rownames_to_column(var = "Sample") %>%
  pivot_longer(!Sample, names_to = "Module", values_to = "Measurement") %>%
  separate(Sample, c("Sample", "Time")) %>%
  as.data.frame(stringsAsFactors = FALSE)

expr$Sample <- factor(expr$Sample, levels = unique(expr$Sample))

expr$Module <- sub("^ME", "", expr$Module)
expr$Module <- factor(expr$Module, levels = unique(expr$Module))

write.csv(
  x = expr,
  file = file.path(output_path, "expr.csv"),
  na = "",
  quote = FALSE
)

# Plot expression plot
p <- ggplot(data = expr, mapping = aes(x = Time, y = Measurement, group=1)) +
  geom_line() +
  facet_grid(Module ~ Sample, scales = "free") +
  labs(x = "Time Point", y = "Expression Value")


ggsave(
  filename = "expression.png",
  plot = p,
  path = output_path,
  width = ifelse(length(unique(expr$Sample))==3, 21, 14),
  height = 35
)

# Create average expression table
averageExpr <- MEList$averageExpr %>%
  rownames_to_column(var = "Sample") %>%
  pivot_longer(!Sample, names_to = "Module", values_to = "Measurement") %>%
  separate(Sample, c("Sample", "Time")) %>%
  as.data.frame(stringsAsFactors = FALSE)

averageExpr$Sample <- factor(averageExpr$Sample, levels = unique(averageExpr$Sample))

averageExpr$Module <- sub("^AE", "", averageExpr$Module)
averageExpr$Module <- factor(averageExpr$Module, levels = unique(averageExpr$Module))

write.csv(
  x = averageExpr,
  file = file.path(output_path, "averageExpr.csv"),
  na = "",
  quote = FALSE
)

# Plot average expression plot
p <- ggplot(data = averageExpr, mapping = aes(x = Time, y = Measurement, group=1)) +
  geom_line() +
  facet_grid(Module ~ Sample, scales = "free") +
  labs(x = "Time Point", y = "Average Expression Value")


ggsave(
  filename = "average_expression.png",
  plot = p,
  path = output_path,
  width = ifelse(length(unique(averageExpr$Sample))==3, 21, 14),
  height = 35
)


# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)

# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")

# Plot the result
cat(rep("\n", 2))
jpeg(filename = file.path(output_path, "moduleEigengenesClustering.jpeg"), width = 1440, height = 720)
plot(
  METree,
  main = "Clustering of module eigengenes",
  xlab = "",
  sub = ""
)

# Plot the cut line into the dendrogram
abline(h=mergeCutHeight, col = "red")
dev.off()


# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = mergeCutHeight, verbose = 3)

# The merged module colors
mergedColors = merge$colors

print(table(mergedColors))

# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

cat(rep("\n", 2))
jpeg(filename = file.path(output_path, "networkConstruction_dynamic_and_merged_colors.jpeg"), width = 1920, height = 720)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

cat(rep("\n", 2))
jpeg(filename = file.path(output_path, "networkConstruction_merged_colors.jpeg"), width = 1920, height = 720)
plotDendroAndColors(geneTree, mergedColors,
                    "Merged dynamic",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()


##################################################
# Plot merged eigengenes
##################################################
expr <- mergedMEs %>%
  rownames_to_column(var = "Sample") %>%
  pivot_longer(!Sample, names_to = "Module", values_to = "Measurement") %>%
  separate(Sample, c("Sample", "Time")) %>%
  arrange(Module) %>%
  as.data.frame(stringsAsFactors = FALSE)

expr$Sample <- factor(expr$Sample, levels = unique(expr$Sample))

expr$Module <- sub("^ME", "", expr$Module)
expr$Module <- factor(expr$Module, levels = unique(expr$Module))

write.csv(
  x = expr,
  file = file.path(output_path, "mergedExpr.csv"),
  na = "",
  quote = FALSE
)

# Plot average expression plot
p <- ggplot(data = expr, mapping = aes(x = Time, y = Measurement, group=1)) +
  geom_line() +
  facet_grid(Module ~ Sample, scales = "free") +
  labs(x = "Time Point", y = "Expression Value")


ggsave(
  filename = "merged_expression.png",
  plot = p,
  path = output_path,
  width = ifelse(length(unique(expr$Sample))==3, 21, 14),
  height = 35
)


##################################################
# Save genes and colors
##################################################
genes_colors_df <- data.frame(
  "Gene" = colnames(datExpr),
  "Color" = mergedColors,
  stringsAsFactors = FALSE
)

write.csv(
  x = genes_colors_df,
  file = file.path(output_path, "genes_colors_df.csv"),
  na = "",
  quote = FALSE,
  row.names = FALSE
)


genes_colors_summary_df <- genes_colors_df %>%
  group_by(Color) %>%
  summarize(Count = n()) %>%
  arrange(Count) %>%
  as.data.frame(stringsAsFactors = FALSE)

genes_colors_summary_df$Color <- factor(genes_colors_summary_df$Color, levels = unique(genes_colors_summary_df$Color))

p <- ggplot(data=genes_colors_summary_df, aes(x = Color, y=Count)) +
  geom_bar(mapping = aes(fill = Color), stat="identity") +
  geom_text(aes(label=Count), hjust=1, color="white", size=3.5) +
  coord_flip() +
  scale_fill_manual(values = levels(genes_colors_summary_df$Color))

ggsave(
  filename = "genes_colors_summary.png",
  plot = p,
  path = output_path
)


##################################################
# Save as RData
##################################################

# Rename to moduleColors
moduleColors = mergedColors

# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1

MEs = mergedMEs

save(
  datExpr, MEs, moduleLabels, moduleColors, geneTree,
  file = file.path(output_path, paste0(selected_genotype, "-networkConstruction-stepByStep.RData"))
)
