#!/usr/bin/Rscript --vanilla
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
output_path <- file.path("/home/ycth8/data/projects/05_30_2021_summer_WGCNA/Maize_proteomics_output/2021_06_10_data_cleaning")

if(!dir.exists(output_path)){
  dir.create(output_path, showWarnings=FALSE, recursive=TRUE)
  if(!dir.exists(output_path)){
    quit(status=1)
  }
}


##################################################
# Read in input file
##################################################

folder_path = file.path("/home/ycth8/data/projects/05_30_2021_summer_WGCNA/Maize_proteomics_data")

dat = read.csv(
  file = file.path(folder_path, "3plus_2fold_genelist.csv"),
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

datTraits = read.csv(
  file = file.path(folder_path, "PBAA_B73o2_Abs_Raw.csv"),
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


##################################################
# Process the input file
##################################################

datExpr0 = dat %>%
  pivot_longer(!c(Accession, Description), names_to = "Genotype", values_to = "Measurement") %>%
  separate(Genotype, c("Genotype", "Repetition"), sep = "\\_(?=[^\\_]+$)", extra = "drop", fill = "right") %>%
  group_by(Accession, Description, Genotype) %>%
  summarise_at(vars(Measurement), list(Measurement = mean)) %>%
  ungroup() %>%
  select(-Description) %>%
  pivot_wider(names_from = Accession, values_from = Measurement, values_fill = NA) %>%
  as.data.frame(stringsAsFactors = FALSE)

# Copy values in column 1 into row names
rownames(datExpr0) <- datExpr0[,1]

# Remove column 1
datExpr0 <- datExpr0[,-1]

print(head(datExpr0))
print(tail(datExpr0))
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
# abline(h = 15, col = "red")
dev.off()

# # Determine cluster under the line
# clust = cutreeStatic(sampleTree, cutHeight = 15.2, minSize = 10)
#
# print(table(clust))
#
# # Remove outlier sample
# keepSamples = (clust==1)
# datExpr = datExpr0[keepSamples, ]

datExpr = datExpr0
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# Collect garbage
collectGarbage()


#######################################################################
## Performs principal components analysis (PCA)
#######################################################################
pca_df <- datExpr
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
# Save processed data
##################################################
write.csv(
  x = datExpr,
  file = file.path(output_path, "..", "datExpr.csv"),
  na = "",
  quote = FALSE
)


##################################################
# Process the trait data
##################################################

datTraits[datTraits[,1] == "o2", 1] <- "O2"

datTraits = datTraits %>%
  pivot_longer(!c(Genotype, DAP, Replicaiton), names_to = "Trait", values_to = "Measurement") %>%
  group_by(Genotype, DAP, Trait) %>%
  summarise_at(vars(Measurement), list(Measurement = mean)) %>%
  ungroup() %>%
  unite("Genotype", Genotype:DAP, sep = "_") %>%
  pivot_wider(names_from = Trait, values_from = Measurement, values_fill = NA) %>%
  as.data.frame(stringsAsFactors = FALSE)

# Copy values in column 1 into row names
rownames(datTraits) <- datTraits[,1]

# Remove column 1
datTraits <- datTraits[,-1]

print(head(datTraits))
print(tail(datTraits))
print(dim(datTraits))


##################################################
# Save processed trait data
##################################################
write.csv(
  x = datTraits,
  file = file.path(output_path, "..", "datTraits.csv"),
  na = "",
  quote = FALSE
)


##################################################
# Construct heatmap to visualize trait data
##################################################
# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")

# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE)

cat(rep("\n", 2))
jpeg(filename = file.path(output_path, "tree_trait_heatmap.jpeg"), width = 1920, height = 720)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors, groupLabels = names(datTraits), main = "Sample dendrogram and trait heatmap")
dev.off()
