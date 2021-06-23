#!/usr/bin/Rscript --vanilla
rm(list=ls())

library(jpeg)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)

set.seed(1)


##################################################
# Constants/Variables
##################################################


##################################################
# Output folder
##################################################
output_path <- file.path("/home/ycth8/data/projects/2021_05_30_summer_WGCNA/Maize_proteomics_output/2021_06_23_B73_O2_average_expression_plot")

if(!dir.exists(output_path)){
  dir.create(output_path, showWarnings=FALSE, recursive=TRUE)
  if(!dir.exists(output_path)){
    quit(status=1)
  }
}


##################################################
# Read in input files
##################################################
# b73_dat = read.csv(
#   file = file.path("/home/ycth8/data/projects/2021_05_30_summer_WGCNA/Maize_proteomics_output/2021_06_22_B73_average_expression_plot/averageExpr.csv"),
#   header = TRUE,
#   row.names = 1,
#   check.names = FALSE,
#   stringsAsFactors = FALSE
# )

b73_dat = read.csv(
  file = file.path("/home/ycth8/data/projects/2021_05_30_summer_WGCNA/Maize_proteomics_output/2021_06_10_B73_step_by_step_network_construction/mergedAverageExpr.csv"),
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# o2_dat = read.csv(
#   file = file.path("/home/ycth8/data/projects/2021_05_30_summer_WGCNA/Maize_proteomics_output/2021_06_22_O2_average_expression_plot/averageExpr.csv"),
#   header = TRUE,
#   row.names = 1,
#   check.names = FALSE,
#   stringsAsFactors = FALSE
# )

o2_dat = read.csv(
  file = file.path("/home/ycth8/data/projects/2021_05_30_summer_WGCNA/Maize_proteomics_output/2021_06_10_O2_step_by_step_network_construction/mergedAverageExpr.csv"),
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


print(head(b73_dat))
print(tail(b73_dat))
print(dim(b73_dat))

print(head(o2_dat))
print(tail(o2_dat))
print(dim(o2_dat))


dat <- b73_dat %>%
  dplyr::left_join(o2_dat, by = "Time", suffix = c(".B73", ".O2")) %>%
  as.data.frame(stringsAsFactors = FALSE)

print(head(dat))
print(tail(dat))
print(dim(dat))


##################################################
# Plot average expression plot
##################################################
p <- ggplot(data = dat, mapping = aes(x = Time, y = Measurement.B73)) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = Time, y = Measurement.B73, group = Module.B73, color = Module.B73)) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = Time, y = Measurement.O2, group = Module.O2, color = Module.O2)) +
  facet_grid(Sample.B73 + Module.B73 ~ Sample.O2 + Module.O2, scales = "free") +
  labs(x = "Time Point", y = "Average Expression Value", color="Color") +
  scale_color_manual(values=unique(c(dat$Module.B73, dat$Module.O2)))


ggsave(
  filename = "average_expression.png",
  plot = p,
  path = output_path,
  width = 20,
  height = 20
)
