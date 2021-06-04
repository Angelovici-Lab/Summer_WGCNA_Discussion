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


##################################################
# Output folder
##################################################
cru_output_path <- file.path("/storage/htc/joshilab/yenc/projects/03_27_2021_aWGCNA/WGCNA_CRU_output")
col_output_path <- file.path("/storage/htc/joshilab/yenc/projects/03_27_2021_aWGCNA/WGCNA_Col_output")
rnai_output_path <- file.path("/storage/htc/joshilab/yenc/projects/03_27_2021_aWGCNA/WGCNA_RNAi_output")


##################################################
# Read in input files
##################################################
lnames = load(file = file.path("/storage/htc/joshilab/yenc/projects/03_27_2021_aWGCNA/WGCNA_Col_output/Col-networkConstruction-stepByStep.RData"))
print(lnames)

col_datExpr <- datExpr
col_MEs <- MEs
col_moduleLabels <- moduleLabels
col_moduleColors <- moduleColors
col_geneTree <- geneTree


lnames = load(file = file.path("/storage/htc/joshilab/yenc/projects/03_27_2021_aWGCNA/WGCNA_CRU_output/CRU-networkConstruction-stepByStep.RData"))
print(lnames)

cru_datExpr <- datExpr
cru_MEs <- MEs
cru_moduleLabels <- moduleLabels
cru_moduleColors <- moduleColors
cru_geneTree <- geneTree


lnames = load(file = file.path("/storage/htc/joshilab/yenc/projects/03_27_2021_aWGCNA/WGCNA_RNAi_output/RNAi-networkConstruction-stepByStep.RData"))
print(lnames)

rnai_datExpr <- datExpr
rnai_MEs <- MEs
rnai_moduleLabels <- moduleLabels
rnai_moduleColors <- moduleColors
rnai_geneTree <- geneTree


col_top_hub <- chooseTopHubInEachModule(col_datExpr, col_moduleColors, omitColors = "")
cru_top_hub <- chooseTopHubInEachModule(cru_datExpr, cru_moduleColors, omitColors = "")
rnai_top_hub <- chooseTopHubInEachModule(rnai_datExpr, rnai_moduleColors, omitColors = "")


col_top_hub_df <- data.frame(
  "Module" = names(col_top_hub),
  "Hub_Gene" = as.vector(col_top_hub),
  stringsAsFactors = FALSE
)

cru_top_hub_df <- data.frame(
  "Module" = names(cru_top_hub),
  "Hub_Gene" = as.vector(cru_top_hub),
  stringsAsFactors = FALSE
)

rnai_top_hub_df <- data.frame(
  "Module" = names(rnai_top_hub),
  "Hub_Gene" = as.vector(rnai_top_hub),
  stringsAsFactors = FALSE
)

print(col_top_hub_df)
print(cru_top_hub_df)
print(rnai_top_hub_df)


write.csv(
  x = col_top_hub_df,
  file = file.path(col_output_path, "col_top_hub_df.csv"),
  na = "",
  quote = FALSE,
  row.names = FALSE
)

write.csv(
  x = cru_top_hub_df,
  file = file.path(cru_output_path, "cru_top_hub_df.csv"),
  na = "",
  quote = FALSE,
  row.names = FALSE
)

write.csv(
  x = rnai_top_hub_df,
  file = file.path(rnai_output_path, "rnai_top_hub_df.csv"),
  na = "",
  quote = FALSE,
  row.names = FALSE
)
