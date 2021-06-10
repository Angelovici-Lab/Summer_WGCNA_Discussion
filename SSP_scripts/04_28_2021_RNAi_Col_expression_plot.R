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
output_path <- file.path("/storage/htc/joshilab/yenc/projects/03_27_2021_aWGCNA/WGCNA_RNAi_Col_output/RNAi_Col_expression_plots")

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
rnai_dat = rnai_MEs[, rnai_index] %>%
  rownames_to_column(var = "Time_Point") %>%
  pivot_longer(!Time_Point, names_to = "Module", values_to = "Measurement") %>%
  separate(Time_Point, c("Genotype", "Time_Point")) %>%
  as.data.frame(stringsAsFactors = FALSE)
rnai_dat$Module <- substring(rnai_dat$Module, 3)
rnai_index <- substring(rnai_index, 3)


##################################################
# Plot expression plots
##################################################

print(head(col_dat))
print(head(rnai_dat))

for(i in rnai_index){
  for(j in col_index){

    temp_rnai_dat = rnai_dat[rnai_dat$Module == i, ]
    temp_col_dat = col_dat[col_dat$Module == j, ]

    temp_dat = rbind(temp_rnai_dat, temp_col_dat)

    print(head(temp_dat))
    print(tail(temp_dat))
    print(dim(temp_dat))

    p <- ggplot(data = temp_dat, mapping = aes(x = Time_Point, y = Measurement, group = Module)) +
      geom_line(mapping = aes(color = Module)) +
      scale_color_manual(values=unique(temp_dat$Module))

    ggsave(
      filename = paste0("RNAi_Col_", i, "_", j, ".png"),
      plot = p,
      path = output_path
    )

  }
}


print(warnings())
