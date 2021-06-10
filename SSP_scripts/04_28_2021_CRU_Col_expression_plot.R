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
output_path <- file.path("/storage/htc/joshilab/yenc/projects/03_27_2021_aWGCNA/WGCNA_CRU_Col_output/CRU_Col_expression_plots")

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
col_dat = col_MEs[, col_index] %>%
  rownames_to_column(var = "Time_Point") %>%
  pivot_longer(!Time_Point, names_to = "Module", values_to = "Measurement") %>%
  separate(Time_Point, c("Genotype", "Time_Point")) %>%
  as.data.frame(stringsAsFactors = FALSE)
col_dat$Module <- substring(col_dat$Module, 3)
col_index <- substring(col_index, 3)


lnames = load(file = file.path("/storage/htc/joshilab/yenc/projects/03_27_2021_aWGCNA/WGCNA_CRU_output/CRU-networkConstruction-stepByStep.RData"))
print(lnames)

cru_datExpr <- datExpr
cru_MEs <- MEs
cru_moduleLabels <- moduleLabels
cru_moduleColors <- moduleColors
cru_geneTree <- geneTree

cru_index <- c("MEcoral1", "MEantiquewhite4", "MEdarkseagreen4", "MElightcyan1", "MEcoral2", "MEfloralwhite", "MEplum2", "MEskyblue1",
               "MEbisque4", "MEdarkolivegreen", "MElightpink4")
cru_dat = cru_MEs[, cru_index] %>%
  rownames_to_column(var = "Time_Point") %>%
  pivot_longer(!Time_Point, names_to = "Module", values_to = "Measurement") %>%
  separate(Time_Point, c("Genotype", "Time_Point")) %>%
  as.data.frame(stringsAsFactors = FALSE)
cru_dat$Module <- substring(cru_dat$Module, 3)
cru_index <- substring(cru_index, 3)


##################################################
# Plot expression plots
##################################################

print(head(col_dat))
print(head(cru_dat))

for(i in cru_index){
  for(j in col_index){

    temp_cru_dat = cru_dat[cru_dat$Module == i, ]
    temp_col_dat = col_dat[col_dat$Module == j, ]

    temp_dat = rbind(temp_cru_dat, temp_col_dat)

    # temp_cru_dat$Module <- factor(temp_cru_dat$Module, levels = unique(temp_cru_dat$Module))
    # temp_col_dat$Module <- factor(temp_col_dat$Module, levels = unique(temp_col_dat$Module))
    # temp_dat$Module <- factor(temp_dat$Module, levels = unique(temp_dat$Module))

    print(head(temp_dat))
    print(tail(temp_dat))
    print(dim(temp_dat))

    p <- ggplot(data = temp_dat, mapping = aes(x = Time_Point, y = Measurement, group = Module)) +
      geom_line(mapping = aes(color = Module)) +
      scale_color_manual(values=unique(temp_dat$Module))
      # geom_line(data = temp_cru_dat, mapping = aes(x = Time_Point, y = Measurement, group = 1, color = Module)) +
      # geom_line(data = temp_col_dat, mapping = aes(x = Time_Point, y = Measurement, group = 1, color = Module))

    ggsave(
      filename = paste0("CRU_Col_", i, "_", j, ".png"),
      plot = p,
      path = output_path
    )

  }
}


print(warnings())
