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


#######################################################################
# Constants/Variables
#######################################################################
selected_genotype <- "O2"
softPower <- 22

# Number to top hub genes
nTop = 30


##################################################
# Output folder
##################################################
output_path <- file.path(
  paste0(
    "/home/ycth8/data/projects/2021_05_30_summer_WGCNA/Maize_proteomics_output/",
    paste0("2021_06_16_", selected_genotype, "_traits_heatmap")
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

folder_path = file.path("/home/ycth8/data/projects/2021_05_30_summer_WGCNA/Maize_proteomics_output")

datTraits = read.csv(
  file = file.path(folder_path, "datTraits.csv"),
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

datTraits = datTraits[startsWith(rownames(datTraits), selected_genotype),]


lnames = load(
  file = file.path(
    folder_path,
    paste0("2021_06_10_", selected_genotype, "_step_by_step_network_construction"),
    paste0(selected_genotype, "-networkConstruction-stepByStep.RData")
  )
)
print(lnames)


# Define numbers of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


##################################################
# Module-trait correlation matrix
##################################################
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)


textMatrix = paste0(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")");
dim(textMatrix) = dim(moduleTraitCor)

modTotals = c()
for(i in 1:length(names(MEs))){
  modTotals = c(modTotals, length(moduleColors[moduleColors == substring(names(MEs)[i], 3)]))
}

cat(rep("\n", 2))
jpeg(filename = file.path(output_path, paste0(selected_genotype, "_PBAA_modules_absolute_traits_correlation.jpeg")),
     pointsize = 15, quality =95, height = 480*1.5, width = 480*1.7)
par(mar = c(5, 10, 3, 3))
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = names(datTraits),
  yLabels = names(MEs),
  ySymbols = paste0(names(MEs), ": ", modTotals),
  colorLabels = FALSE,
  colors = greenWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.5,
  zlim = c(-1,1),
  main = paste(selected_genotype, " module-absolute_trait relationships")
)
dev.off()


##################################################
# Get top hub gene using WGCNA chooseTopHubInEachModule
##################################################
top_hub <- chooseTopHubInEachModule(datExpr, moduleColors, omitColors = "")

top_hub_df <- data.frame(
  "Module" = names(top_hub),
  "Hub_Gene" = as.vector(top_hub),
  stringsAsFactors = FALSE
)

write.csv(
  x = top_hub_df,
  file = file.path(output_path, paste0(selected_genotype, "_top_hub_df.csv")),
  na = "",
  quote = FALSE,
  row.names = FALSE
)


##################################################
# Calculate kME
##################################################
signed_KME_df <- as.data.frame(signedKME(datExpr, MEs), stringsAsFactors = FALSE) %>%
  rownames_to_column(var = "Gene") %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  pivot_longer(!Gene, names_to = "Color", values_to = "kME") %>%
  mutate(Color = gsub("^kME", "", Color)) %>%
  arrange(desc(kME)) %>%
  arrange(desc(abs(kME))) %>%
  as.data.frame(stringsAsFactors = FALSE)

write.csv(
  x = signed_KME_df,
  file = file.path(output_path, paste0(selected_genotype, "_kME_df.csv")),
  na = "",
  quote = FALSE,
  row.names = FALSE
)


##################################################
# Calculate soft connectivity
##################################################
IMConn = softConnectivity(datExpr)

IMConn_df <- data.frame(
  "Gene" = names(datExpr),
  "SoftConnectivity" = IMConn,
  stringsAsFactors = FALSE
)

IMConn_df = IMConn_df %>%
  arrange(desc(SoftConnectivity)) %>%
  as.data.frame(stringsAsFactors = FALSE)

write.csv(
  x = IMConn_df,
  file = file.path(output_path, paste0(selected_genotype, "_softConnectivity.csv")),
  na = "",
  quote = FALSE,
  row.names = FALSE
)

top = (rank(-IMConn) <= nTop)


##################################################
# Export nodes and edges
##################################################
# Select modules
# modules = c("brown", "red")
modules = unique(moduleColors)

# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(datExpr, power = softPower)

# Select module probes
probes = names(datExpr)
inModule = is.finite(match(moduleColors, modules))
modProbes = probes[inModule]
modGenes = modProbes

# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]

dimnames(modTOM) = list(modProbes, modProbes)

# print(head(modTOM))

# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(
  modTOM,
  edgeFile = file.path(output_path, paste0(selected_genotype, "_CytoscapeInput-edges-", paste(modules, collapse = "-"), ".txt")),
  nodeFile = file.path(output_path, paste0(selected_genotype, "_CytoscapeInput-nodes-", paste(modules, collapse = "-"), ".txt")),
  weighted = TRUE,
  threshold = 0.02,
  nodeNames = modProbes,
  altNodeNames = modGenes,
  nodeAttr = moduleColors[inModule]
)
