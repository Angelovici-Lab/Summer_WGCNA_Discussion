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


#######################################################################
# Constants/Variables
#######################################################################


##################################################
# Output folder
##################################################
output_path <- file.path(
  "/home/ycth8/data/projects/2021_05_30_summer_WGCNA/Maize_proteomics_output/2021_06_16_overlap_B73_O2_modules"
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

selected_genotype <- "B73"
lnames = load(
  file = file.path(
    folder_path,
    paste0("2021_06_10_", selected_genotype, "_step_by_step_network_construction"),
    paste0(selected_genotype, "-networkConstruction-stepByStep.RData")
  )
)
print(lnames)

# Replace B73 realted variable names
b73Expr = datExpr

index <- match(sub("^ME", "", colnames(MEs)), moduleColors)

colnames(MEs) <- paste0("ME", moduleLabels[index])

b73MEs = orderMEs(MEs, greyName = "ME0")
b73Labels = moduleLabels
b73Colors = moduleColors
b73Tree = geneTree

cat(rep("\n", 2))
print(head(b73Labels))
print(head(b73Colors))
print(head(b73MEs))


selected_genotype <- "O2"
lnames = load(
  file = file.path(
    folder_path,
    paste0("2021_06_10_", selected_genotype, "_step_by_step_network_construction"),
    paste0(selected_genotype, "-networkConstruction-stepByStep.RData")
  )
)
print(lnames)


# Replace O2 realted variable names
o2Expr = datExpr

index <- match(sub("^ME", "", colnames(MEs)), moduleColors)

colnames(MEs) <- paste0("ME", moduleLabels[index])

o2MEs = orderMEs(MEs, greyName = "ME0")
o2Labels = moduleLabels
o2Colors = moduleColors
o2Tree = geneTree

cat(rep("\n", 2))
print(head(o2Labels))
print(head(o2Colors))
print(head(o2MEs))


# Isolate the module labels in the order they appear in ordered module eigengenes
b73ModuleLabels = substring(names(b73MEs), 3)
o2ModuleLabels = substring(names(o2MEs), 3)


# Convert the numeric module labels to color labels
b73Modules = labels2colors(as.numeric(b73ModuleLabels))
o2Modules = labels2colors(as.numeric(o2ModuleLabels))

# Numbers of B73 and O2 modules
nB73Mods = length(b73Modules)
nO2Mods = length(o2Modules)

# Initialize tables of p-values and of the corresponding counts
pTable = matrix(0, nrow = nB73Mods, ncol = nO2Mods)
CountTbl = matrix(0, nrow = nB73Mods, ncol = nO2Mods)


for (i in 1:nB73Mods) {
  for (cmod in 1:nO2Mods) {
    b73Members = (b73Colors == b73Modules[i])
    o2Members = (o2Colors == o2Modules[cmod])
    pTable[i, cmod] = -log10(fisher.test(b73Members, o2Members, alternative = "greater")$p.value)
    CountTbl[i, cmod] = sum(b73Colors == b73Modules[i] & o2Colors == o2Modules[cmod])
  }
}


# Truncate p values smaller than 10^{-50} to 10^{-50}
pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)])
pTable[pTable>50 ] = 50

# Marginal counts (really module sizes)
b73ModTotals = apply(CountTbl, 1, sum)
o2ModTotals = apply(CountTbl, 2, sum)


cat(rep("\n", 2))
jpeg(file = file.path(output_path, "B73_vs_O2.jpg"), height = 480*2, width = 480*2)
par(mfrow=c(1,1))
par(cex = 1.0)
par(mar=c(10, 12, 3, 1))
labeledHeatmap(Matrix = pTable,
               xLabels = paste(" ", o2Modules),
               yLabels = paste(" ", b73Modules),
               colorLabels = TRUE,
               xSymbols = paste0("O2 ", o2Modules, ": ", o2ModTotals),
               ySymbols = paste0("B73 ", b73Modules, ": ", b73ModTotals),
               textMatrix = CountTbl,
               colors = greenWhiteRed(100)[50:100],
               main = "Correspondence of B73 set-specific and O2 set-specific modules",
               cex.text = 1.0, cex.lab = 1.0, setStdMargins = FALSE)

dev.off()
