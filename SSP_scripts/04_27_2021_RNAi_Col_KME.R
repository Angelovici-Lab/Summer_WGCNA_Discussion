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
output_path <- file.path("/storage/htc/joshilab/yenc/projects/03_27_2021_aWGCNA/WGCNA_RNAi_Col_output")

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


##################################################
# Determines significant overlap between modules
# in two networks based on kME tables
##################################################

full_overlap_genes_df <- data.frame(
  "Gene" = NA,
  "RNAi_Color" = NA,
  "Col_Color" = NA,
  stringsAsFactors = FALSE
)
full_overlap_genes_df <- full_overlap_genes_df[!is.na(full_overlap_genes_df$Gene),]

overlop_KME_results <- overlapTableUsingKME(
  rnai_datExpr, col_datExpr,
  rnai_moduleColors, col_moduleColors,
  MEs1 = rnai_MEs, MEs2 = col_MEs,
  name1 = "MM1", name2 = "MM2",
  cutoffMethod = "assigned", cutoff = 0.5,
  omitGrey = TRUE, datIsExpression = TRUE
)

write.csv(
  x = overlop_KME_results$PvaluesHypergeo,
  file = file.path(output_path, "RNAi_Col_KME_PvaluesHypergeo.csv"),
  na = "",
  quote = FALSE
)

for(i in rnai_index){
  for(j in col_index){

    element <- paste0("MM1_", i, "_MM2_", j)

    overlap_genes <- overlop_KME_results$OverlappingGenes[[element]]

    if(length(overlap_genes) > 0){
      overlap_genes_df <- data.frame(
        "Gene" = overlap_genes,
        "RNAi_Color" = i,
        "Col_Color" = j,
        stringsAsFactors = FALSE
      )

      full_overlap_genes_df <- rbind(full_overlap_genes_df, overlap_genes_df)

    }

  }
}
print(head(full_overlap_genes_df))
print(tail(full_overlap_genes_df))


col_signed_KME <- as.data.frame(signedKME(col_datExpr, col_MEs), stringsAsFactors = FALSE) %>%
  rownames_to_column(var = "Gene") %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  pivot_longer(!Gene, names_to = "Col_Color", values_to = "Col_kME") %>%
  arrange(desc(Col_kME)) %>%
  arrange(desc(abs(Col_kME))) %>%
  as.data.frame(stringsAsFactors = FALSE)

rnai_signed_KME <- as.data.frame(signedKME(rnai_datExpr, rnai_MEs), stringsAsFactors = FALSE) %>%
  rownames_to_column(var = "Gene") %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  pivot_longer(!Gene, names_to = "RNAi_Color", values_to = "RNAi_kME") %>%
  arrange(desc(RNAi_kME)) %>%
  arrange(desc(abs(RNAi_kME))) %>%
  as.data.frame(stringsAsFactors = FALSE)

col_signed_KME$Col_Color <- substring(col_signed_KME$Col_Color, 4)
rnai_signed_KME$RNAi_Color <- substring(rnai_signed_KME$RNAi_Color, 4)

print(head(rnai_signed_KME))
print(tail(rnai_signed_KME))
print(head(col_signed_KME))
print(tail(col_signed_KME))


full_overlap_genes_df <- full_overlap_genes_df %>%
  left_join(rnai_signed_KME, by = c("Gene", "RNAi_Color")) %>%
  left_join(col_signed_KME, by = c("Gene", "Col_Color")) %>%
  group_by(RNAi_Color, Col_Color) %>%
  arrange(desc(RNAi_kME), desc(Col_kME), .by_group = TRUE) %>%
  arrange(desc(abs(RNAi_kME)), desc(abs(Col_kME)), .by_group = TRUE) %>%
  ungroup() %>%
  as.data.frame(stringsAsFactors = FALSE)

print(head(full_overlap_genes_df))
print(tail(full_overlap_genes_df))


write.csv(
  x = full_overlap_genes_df,
  file = file.path(output_path, "RNAi_Col_KME.csv"),
  na = "",
  quote = FALSE,
  row.names = FALSE
)
