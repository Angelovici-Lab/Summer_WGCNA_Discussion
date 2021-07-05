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

b73_selected_module <- "yellow"
o2_selected_module <- "turquoise"


##################################################
# Output folder
##################################################
output_path <- file.path(
  "/home/ycth8/data/projects/2021_05_30_summer_WGCNA/Maize_proteomics_output/2021_07_05_plot_intersected_genes_line_plot"
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


##################################################
# Overlap modules
##################################################

b73_genes_colors_df <- data.frame(
  "Gene" = colnames(b73Expr),
  "Color" = b73Colors,
  stringsAsFactors = FALSE
)

o2_genes_colors_df <- data.frame(
  "Gene" = colnames(o2Expr),
  "Color" = o2Colors,
  stringsAsFactors = FALSE
)


b73_genes_colors_df <- b73_genes_colors_df %>%
  filter(Color == b73_selected_module) %>%
  as.data.frame(stringsAsFactors = FALSE)

o2_genes_colors_df <- o2_genes_colors_df %>%
  filter(Color == o2_selected_module) %>%
  as.data.frame(stringsAsFactors = FALSE)

print(head(b73_genes_colors_df))
print(head(o2_genes_colors_df))


genes_colors_df <- b73_genes_colors_df %>%
  inner_join(o2_genes_colors_df, by = "Gene") %>%
  as.data.frame(stringsAsFactors = FALSE)

print(head(genes_colors_df))


b73Expr <- b73Expr[, genes_colors_df$Gene] %>%
  rownames_to_column(var = "Sample") %>%
  separate(Sample, c("Sample", "Timepoint"), sep = "\\_(?=[^\\_]+$)", extra = "drop", fill = "right") %>%
  pivot_longer(!c("Sample", "Timepoint"), names_to = "Genotype", values_to = "Measurement") %>%
  as.data.frame(stringsAsFactors = FALSE)
o2Expr <- o2Expr[, genes_colors_df$Gene] %>%
  rownames_to_column(var = "Sample") %>%
  separate(Sample, c("Sample", "Timepoint"), sep = "\\_(?=[^\\_]+$)", extra = "drop", fill = "right") %>%
  pivot_longer(!c("Sample", "Timepoint"), names_to = "Genotype", values_to = "Measurement") %>%
  as.data.frame(stringsAsFactors = FALSE)


expr <- rbind(b73Expr, o2Expr)

expr$Sample_Genotype <- paste(expr$Sample, expr$Genotype, sep = "_")

print(head(expr))
print(tail(expr))


##################################################
# Plot expression of genes
##################################################

expr$Timepoint <- factor(expr$Timepoint, levels = unique(expr$Timepoint))
expr$Genotype <- factor(expr$Genotype, levels = unique(expr$Genotype))
expr$Sample <- factor(expr$Sample, levels = unique(expr$Sample))
expr$Sample_Genotype <- factor(expr$Sample_Genotype, levels = unique(expr$Sample_Genotype))

p <- ggplot(data=expr, mapping=aes(x=Timepoint, y=Measurement, group=Sample_Genotype, color=Genotype, linetype=Sample, shape=Sample)) +
  geom_line()+
  geom_point()

ggsave(
  filename=paste0(b73_selected_module, "_", o2_selected_module, "_line_plot.jpg"),
  plot=p,
  path=output_path
)
