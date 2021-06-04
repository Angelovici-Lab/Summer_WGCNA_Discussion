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
output_path <- file.path("/storage/htc/joshilab/yenc/projects/03_27_2021_aWGCNA/WGCNA_CRU_output")

if(!dir.exists(output_path)){
  dir.create(output_path, showWarnings=FALSE, recursive=TRUE)
  if(!dir.exists(output_path)){
    quit(status=1)
  }
}


##################################################
# Read in input files
##################################################
lnames = load(file = file.path("/storage/htc/joshilab/yenc/projects/03_27_2021_aWGCNA/WGCNA_CRU_output/CRU-networkConstruction-stepByStep.RData"))
print(lnames)

datTraits <- read.csv(
  file = file.path("/storage/htc/joshilab/yenc/projects/03_27_2021_aWGCNA/processed_inputs/TAA_SSP_Developmental_R1R2R3_filtered_time_points.csv"),
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

datTraits2 <- read.csv(
  file = file.path("/storage/htc/joshilab/yenc/projects/03_27_2021_aWGCNA/processed_inputs/TAA_SSP_Developmental_R1R2R3_filtered_time_points_absolute_vs_total_dat.csv"),
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

datTraits3 <- read.csv(
  file = file.path("/storage/htc/joshilab/yenc/projects/03_27_2021_aWGCNA/processed_inputs/FAA_SSP_Developmental_R1R2R3_filtered_time_points.csv"),
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

datTraits4 <- read.csv(
  file = file.path("/storage/htc/joshilab/yenc/projects/03_27_2021_aWGCNA/processed_inputs/FAA_SSP_Developmental_R1R2R3_filtered_time_points_absolute_vs_total_dat.csv"),
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

datTraits <- datTraits[startsWith(rownames(datTraits), "CRU"),]
datTraits2 <- datTraits2[startsWith(rownames(datTraits2), "CRU"),]
datTraits3 <- datTraits3[startsWith(rownames(datTraits3), "CRU"),]
datTraits4 <- datTraits4[startsWith(rownames(datTraits4), "CRU"),]

print(head(datTraits))
print(tail(datTraits))
print(dim(datTraits))

print(head(datTraits2))
print(tail(datTraits2))
print(dim(datTraits2))

print(head(datTraits3))
print(tail(datTraits3))
print(dim(datTraits3))

print(head(datTraits4))
print(tail(datTraits4))
print(dim(datTraits4))

# Define numbers of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

##################################################
# Quantifying module–trait associations (absolute)
##################################################
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

textMatrix = paste0(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")");
dim(textMatrix) = dim(moduleTraitCor)

modTotals = c()
for(i in 1:length(names(MEs))){
  modTotals = c(modTotals, length(moduleColors[moduleColors == substring(names(MEs)[i], 3)]))
}

cat(rep("\n", 2))
jpeg(filename = file.path(output_path, "TAA_modules_absolute_traits_correlation.jpeg"),
     pointsize = 15, quality =95, height = 480*1.5, width = 480*1.7)
par(mar = c(5, 10, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = paste0(names(MEs), ": ", modTotals),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-absolute_trait relationships"))
dev.off()


##################################################
# Quantifying module–trait associations (compositional)
##################################################
moduleTraitCor = cor(MEs, datTraits2, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

textMatrix = paste0(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")");
dim(textMatrix) = dim(moduleTraitCor)

modTotals = c()
for(i in 1:length(names(MEs))){
  modTotals = c(modTotals, length(moduleColors[moduleColors == substring(names(MEs)[i], 3)]))
}

cat(rep("\n", 2))
jpeg(filename = file.path(output_path, "TAA_modules_compositional_traits_correlation.jpeg"),
     pointsize = 15, quality =95, height = 480*1.5, width = 480*1.7)
par(mar = c(5, 10, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits2),
               yLabels = names(MEs),
               ySymbols = paste0(names(MEs), ": ", modTotals),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-compositional_trait relationships"))
dev.off()


##################################################
# Quantifying module–trait associations (absolute)
##################################################
moduleTraitCor = cor(MEs, datTraits3, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

textMatrix = paste0(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")");
dim(textMatrix) = dim(moduleTraitCor)

modTotals = c()
for(i in 1:length(names(MEs))){
  modTotals = c(modTotals, length(moduleColors[moduleColors == substring(names(MEs)[i], 3)]))
}

cat(rep("\n", 2))
jpeg(filename = file.path(output_path, "FAA_modules_absolute_traits_correlation.jpeg"),
     pointsize = 15, quality =95, height = 480*1.5, width = 480*1.7)
par(mar = c(5, 10, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits3),
               yLabels = names(MEs),
               ySymbols = paste0(names(MEs), ": ", modTotals),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-absolute_trait relationships"))
dev.off()


##################################################
# Quantifying module–trait associations (compositional)
##################################################
moduleTraitCor = cor(MEs, datTraits4, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

textMatrix = paste0(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")");
dim(textMatrix) = dim(moduleTraitCor)

modTotals = c()
for(i in 1:length(names(MEs))){
  modTotals = c(modTotals, length(moduleColors[moduleColors == substring(names(MEs)[i], 3)]))
}

cat(rep("\n", 2))
jpeg(filename = file.path(output_path, "FAA_modules_compositional_traits_correlation.jpeg"),
     pointsize = 15, quality =95, height = 480*1.5, width = 480*1.7)
par(mar = c(5, 10, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits4),
               yLabels = names(MEs),
               ySymbols = paste0(names(MEs), ": ", modTotals),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-compositional_trait relationships"))
dev.off()
