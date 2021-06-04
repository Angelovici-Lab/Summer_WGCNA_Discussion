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
output_path <- file.path("/storage/htc/joshilab/yenc/projects/03_27_2021_aWGCNA/WGCNA_Col_CRU_RNAi_output")

if(!dir.exists(output_path)){
  dir.create(output_path, showWarnings=FALSE, recursive=TRUE)
  if(!dir.exists(output_path)){
    quit(status=1)
  }
}


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

col_index <- c("MEsienna3", "MEsteelblue", "MEwhite", "MEdarkorange", "MElightsteelblue1", "MEbrown4", "MEdarkturquoise", "MEcyan", "MEdarkred",
               "MEpurple", "MEdarkslateblue", "MEsalmon4", "MEivory", "MEskyblue")
col_index <- substring(col_index, 3)

col_genes_colors <- data.frame(
  "Gene" = colnames(col_datExpr),
  "Color" = col_moduleColors,
  stringsAsFactors = FALSE
)
col_genes_colors <- col_genes_colors[col_genes_colors$Color %in% col_index,]

print(head(col_genes_colors))
print(dim(col_genes_colors))


lnames = load(file = file.path("/storage/htc/joshilab/yenc/projects/03_27_2021_aWGCNA/WGCNA_CRU_output/CRU-networkConstruction-stepByStep.RData"))
print(lnames)

cru_datExpr <- datExpr
cru_MEs <- MEs
cru_moduleLabels <- moduleLabels
cru_moduleColors <- moduleColors
cru_geneTree <- geneTree

cru_index <- c("MEcoral1", "MEantiquewhite4", "MEdarkseagreen4", "MElightcyan1", "MEcoral2", "MEfloralwhite", "MEplum2", "MEskyblue1",
               "MEbisque4", "MEdarkolivegreen", "MElightpink4")
cru_index <- substring(cru_index, 3)

cru_genes_colors <- data.frame(
  "Gene" = colnames(cru_datExpr),
  "Color" = cru_moduleColors,
  stringsAsFactors = FALSE
)
cru_genes_colors <- cru_genes_colors[cru_genes_colors$Color %in% cru_index,]

print(head(cru_genes_colors))
print(dim(cru_genes_colors))


lnames = load(file = file.path("/storage/htc/joshilab/yenc/projects/03_27_2021_aWGCNA/WGCNA_RNAi_output/RNAi-networkConstruction-stepByStep.RData"))
print(lnames)

rnai_datExpr <- datExpr
rnai_MEs <- MEs
rnai_moduleLabels <- moduleLabels
rnai_moduleColors <- moduleColors
rnai_geneTree <- geneTree

rnai_index <- c(
  "MEdarkturquoise", "MEpaleturquoise", "MElightyellow", "MEmediumpurple3", "MEdarkred", "MEbrown", "MEyellowgreen", "MEbrown4",
  "MEdarkorange2", "MEdarkgrey", "MEorangered4", "MEcyan", "MEdarkolivegreen", "MEplum1"
)
rnai_index <- substring(rnai_index, 3)

rnai_genes_colors <- data.frame(
  "Gene" = colnames(rnai_datExpr),
  "Color" = rnai_moduleColors,
  stringsAsFactors = FALSE
)
rnai_genes_colors <- rnai_genes_colors[rnai_genes_colors$Color %in% rnai_index,]

print(head(rnai_genes_colors))
print(dim(rnai_genes_colors))


##################################################
# Intersect genes
##################################################

col_genes <- col_genes_colors$Gene
cru_genes <- cru_genes_colors$Gene
rnai_genes <- rnai_genes_colors$Gene


all_intersect_genes <- col_genes[(col_genes %in% cru_genes) & (col_genes %in% rnai_genes)]

print(head(all_intersect_genes))
print(length(all_intersect_genes))

col_cru_genes <- col_genes[(col_genes %in% cru_genes) & !(col_genes %in% all_intersect_genes)]
cru_rnai_genes <- cru_genes[(cru_genes %in% rnai_genes) & !(cru_genes %in% all_intersect_genes)]
rnai_col_genes <- rnai_genes[(rnai_genes %in% col_genes) & !(rnai_genes %in% all_intersect_genes)]

print(length(col_cru_genes))
print(length(cru_rnai_genes))
print(length(rnai_col_genes))

col_unique_genes <- col_genes[!(col_genes %in% col_cru_genes) & !(col_genes %in% rnai_col_genes) & !(col_genes %in% all_intersect_genes)]
cru_unique_genes <- cru_genes[!(cru_genes %in% col_cru_genes) & !(cru_genes %in% cru_rnai_genes) & !(cru_genes %in% all_intersect_genes)]
rnai_unique_genes <- rnai_genes[!(rnai_genes %in% cru_rnai_genes) & !(rnai_genes %in% rnai_col_genes) & !(rnai_genes %in% all_intersect_genes)]

print(length(col_unique_genes))
print(length(cru_unique_genes))
print(length(rnai_unique_genes))
