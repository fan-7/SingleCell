Sys.setenv(PATH = paste("~/miniconda3/envs/astro/bin/", Sys.getenv("PATH"), sep = ":"))
core_packages <- c("Seurat", "dplyr", "ggplot2", "Signac", "RColorBrewer", "scSHC", 
                   "tidyr", "data.table", "stringr", "tidyverse", "CytoTRACE2", 
                   "viridis", "org.Mm.eg.db", "pheatmap", "ggridges", "gridExtra", 
                   "harmony", "slingshot", "monocle")
suppressPackageStartupMessages({lapply(core_packages, library, character.only = TRUE)})
