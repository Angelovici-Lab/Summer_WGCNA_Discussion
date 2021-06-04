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
library(KEGGREST)
# library(biomaRt)

set.seed(1)

# Enable WGCNA threads to speed up calculations
enableWGCNAThreads()


##################################################
# Constants/Variables
##################################################
softPower <- 30
minModuleSize <- 30

# Eigengenes clustering tree cutting threshold
MEDissThres <- 0.1


##################################################
# Output folder
##################################################
output_path <- file.path("/storage/htc/joshilab/yenc/projects/03_27_2021_aWGCNA/WGCNA_Col_output")

if(!dir.exists(output_path)){
  dir.create(output_path, showWarnings=FALSE, recursive=TRUE)
  if(!dir.exists(output_path)){
    quit(status=1)
  }
}


##################################################
# Read in input files
##################################################
folder_path <- file.path("/storage/htc/joshilab/yenc/projects/03_27_2021_aWGCNA")

dat <- read.csv(
  file=file.path(folder_path, "processed_inputs", "STAR_out_counts_cpm_log_edgeR_filtered_genes.csv"),
  header=TRUE,
  row.names = 1,
  comment.char="#",
  check.names=FALSE,
  stringsAsFactors = FALSE
)

dat <- dat[(startsWith(rownames(dat), "Col")),]

cat(rep("\n", 2))
print((dat[1:6, 1:6]))
print(dim(dat))

cat(rep("\n", 2))
print(rownames(dat))
print(length(rownames(dat)))


#######################################################################
## Get data from KEGG using KEGGREST
#######################################################################
ath_gene <- keggList("ath")
ath_gene_df <- as.data.frame(ath_gene, row.names = names(ath_gene), stringsAsFactors = FALSE) %>%
  rownames_to_column(var = "ath_id") %>%
  separate(ath_gene, c("ath_gene", "ath_gene_description"), sep = "; ", extra = "drop", fill = "left") %>%
  separate_rows(ath_gene, sep = ", ", convert = TRUE) %>%
  as.data.frame(stringsAsFactors = FALSE)
ath_gene_df$ath_id <- gsub("ath:", "", ath_gene_df$ath_id)
print(head(ath_gene_df))

ath_pathway_id <- keggLink("pathway", "ath")
ath_pathway_id_df <- as.data.frame(ath_pathway_id, row.names = names(ath_pathway_id), stringsAsFactors = FALSE) %>%
  rownames_to_column(var = "ath_id")
ath_pathway_id_df$ath_id <- gsub("ath:", "", ath_pathway_id_df$ath_id)
print(head(ath_pathway_id_df))

ath_pathway <- keggList("pathway", "ath")
ath_pathway_df <- as.data.frame(ath_pathway, row.names = names(ath_pathway), stringsAsFactors = FALSE) %>%
  rownames_to_column(var = "ath_pathway_id")
ath_pathway_df$ath_pathway <- gsub(" \\- Arabidopsis thaliana \\(thale cress)", "", ath_pathway_df$ath_pathway)
print(head(ath_pathway_df))


# ##################################################
# # Get ensembl info
# ##################################################
# # ensembl <- useMart("ensembl", dataset="athaliana_eg_gene")
# ensembl <- useMart(biomart="plants_mart", host="plants.ensembl.org", dataset="athaliana_eg_gene")
# ath_gene_ids  <- unique(colnames(dat))
#
# ensembl_dat <- getBM(attributes=c('external_gene_name','ensembl_gene_id'),
#              filters = "ensembl_gene_id",
#              values = ath_gene_ids,
#              mart = ensembl)
#
# print(head(ensembl_dat))
# print(tail(ensembl_dat))
# print(dim(ensembl_dat))


##################################################
# Check if the genes are good using WGCNA::goodSamplesGenes
##################################################
gsg = goodSamplesGenes(dat, verbose = 3)

print(gsg$allOK)
# print(gsg$goodGenes)
# print(gsg$goodSamples)

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(dat)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(dat)[!gsg$goodSamples], collapse = ", ")))
  # Remove the offending genes and samples from the data:
  dat = dat[gsg$goodSamples, gsg$goodGenes]
}


#######################################################################
## Performs a principal components analysis
#######################################################################
pca_df <- dat
pca <- prcomp(pca_df, center = TRUE, scale. = TRUE)

pca_variance <- pca$sdev^2
pca_variance_percentage <- round(pca_variance/sum(pca_variance)*100, digits=2)

x <- pca$x %>% as.data.frame(stringsAsFactors=TRUE)


#######################################################################
## Plot PC1 vs PC2 with GGplot2
#######################################################################
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
# Cluster the data
##################################################
sampleTree = hclust(dist(dat), method = "average")

# Plot tree
cat(rep("\n", 2))
jpeg(filename = file.path(output_path, "sampleTree.jpeg"))
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
# abline(h = 300, col = "red")
dev.off()


# Collect garbage
collectGarbage()


##################################################
# Choose a set of soft-thresholding powers
##################################################
powers = 1:60
# Call the network topology analysis function
sft = pickSoftThreshold(dat, powerVector = powers, verbose = 5)

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
adjacency = adjacency(dat, power = softPower)

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
     labels = colnames(dat), hang = 0.04)
dev.off()


# We like large modules, so we set the minimum module size relatively high
# Check minModuleSize in the constants/variables section

# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
print(table(dynamicMods))


# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)

print(table(dynamicColors))

# Plot the dendrogram and colors underneath
cat(rep("\n", 2))
jpeg(filename = file.path(output_path, "networkConstruction_dynamic_colors.jpeg"), width = 1920, height = 720)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()


##################################################
# Modules merging by clustering module eigengenes
##################################################
# Calculate eigengenes
MEList = moduleEigengenes(dat, colors = dynamicColors)
MEs = MEList$eigengenes

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
abline(h=MEDissThres, col = "red")
dev.off()




# Call an automatic merging function
merge = mergeCloseModules(dat, dynamicColors, cutHeight = MEDissThres, verbose = 3)

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
  "Gene" = colnames(dat),
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


genes_colors_df <- genes_colors_df %>%
  left_join(ath_gene_df, by = c("Gene" = "ath_id")) %>%
  distinct() %>%
  left_join(ath_pathway_id_df, by = c("Gene" = "ath_id")) %>%
  distinct() %>%
  left_join(ath_pathway_df, by = "ath_pathway_id") %>%
  distinct() %>%
  group_by(Gene, Color, ath_gene, ath_gene_description) %>%
  mutate(
    ath_pathway_id = paste(ath_pathway_id, collapse = "; "),
    ath_pathway = paste(ath_pathway, collapse = "; ")
  ) %>%
  distinct(Gene, .keep_all = TRUE) %>%
  as.data.frame(stringsAsFactors = FALSE)

genes_colors_df$ath_pathway_id[is.na(genes_colors_df$ath_pathway_id) | genes_colors_df$ath_pathway_id == "NA"] = NA
genes_colors_df$ath_pathway[is.na(genes_colors_df$ath_pathway) | genes_colors_df$ath_pathway == "NA"] = NA

write.table(
  x = genes_colors_df,
  file = file.path(output_path, "genes_colors_df.txt"),
  sep = "\t",
  na = "",
  quote = FALSE,
  row.names = FALSE
)


##################################################
# Save as RData
##################################################
datExpr = dat

# Rename to moduleColors
moduleColors = mergedColors

# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1

MEs = mergedMEs

save(datExpr, MEs, moduleLabels, moduleColors, geneTree,
     file = file.path(output_path, "Col-networkConstruction-stepByStep.RData"))
