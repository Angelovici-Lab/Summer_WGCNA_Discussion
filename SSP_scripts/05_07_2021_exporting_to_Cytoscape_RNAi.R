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
softPower <- 30

# Number to top hub genes
nTop = 30


##################################################
# Output folder
##################################################
output_path <- file.path("/storage/htc/joshilab/yenc/projects/03_27_2021_aWGCNA/WGCNA_RNAi_output/Networks")

if(!dir.exists(output_path)){
  dir.create(output_path, showWarnings=FALSE, recursive=TRUE)
  if(!dir.exists(output_path)){
    quit(status=1)
  }
}


##################################################
# Read in input files
##################################################
lnames = load(file = file.path("/storage/htc/joshilab/yenc/projects/03_27_2021_aWGCNA/WGCNA_RNAi_output/RNAi-networkConstruction-stepByStep.RData"))
print(lnames)


##################################################
# Calculate soft connectivity
##################################################
print(moduleColors)
print(length(moduleColors))

IMConn = softConnectivity(datExpr)

IMConn_df <- data.frame(
  "Gene" = names(datExpr),
  "SoftConnectivity" = IMConn,
  stringsAsFactors = FALSE
)
IMConn_df = IMConn_df %>% arrange(desc(SoftConnectivity)) %>% as.data.frame(stringsAsFactors = FALSE)
write.csv(
  x = IMConn_df,
  file = file.path(output_path, "softConnectivity.csv"),
  na = "",
  quote = FALSE,
  row.names = FALSE
)

top = (rank(-IMConn) <= nTop)

# print(IMConn[top])
# print(names(datExpr)[top])

print(length(IMConn[top]))
print(length(names(datExpr)[top]))


##################################################
# Export nodes and edges
##################################################
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

print(head(modTOM))

# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(
  modTOM,
  edgeFile = file.path(output_path, paste0("CytoscapeInput-edges-", paste(modules, collapse = "-"), ".txt")),
  nodeFile = file.path(output_path, paste0("CytoscapeInput-nodes-", paste(modules, collapse = "-"), ".txt")),
  weighted = TRUE,
  threshold = 0.02,
  nodeNames = modProbes,
  altNodeNames = modGenes,
  nodeAttr = moduleColors[inModule]
)
