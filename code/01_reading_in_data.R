# Set up packages for scRNA analysis using Seurat ----

## Install the required packages ----

install.packages("BiocManager")
install.packages('Seurat')
install.packages("devtools")
remotes::install_github("mojaveazure/seurat-disk")
BiocManager::install("SingleR")
BiocManager::install("celldex")
BiocManager::install("scran")
BiocManager::install("limma")
BiocManager::install("org.Mm.eg.db")
install.packages("vegan")
install.packages("pheatmap")
install.packages("tidyverse")
install.packages("enrichR")

sessionInfo()
# these libraries were picked from the BIOS 289 Spring class 
# and from the Seurat Parse tutorial 
# https://support.parsebiosciences.com/hc/en-us/articles/360053078092-Seurat-Tutorial-65k-PBMCs

## Load installed packages ----

library(Seurat)
library(SeuratDisk)
library(SingleR)
library(celldex)
library(dplyr)
library(Matrix)
library(scran)
library(vegan)
library(pheatmap)
library(ggplot2)
library(tidyverse)
library(enrichR)
library(ggpubr)


# Location paths and convenience functions --------------------------------

## Setting Location paths ----

data_path <- "~/Documents/bhutani_lab/PGDH_scRNA_Mamta/analysis/PGDH_chondro_project/data/"
fig_path <- "~/Documents/bhutani_lab/PGDH_scRNA_Mamta/analysis/PGDH_chondro_project/figures/"

## Convenience functions ----

SaveFigure <- function(plots, name, type = "png", width, height, res){
  if(type == "png") {
    png(paste0(fig_path, name, ".", type),
        width = width, height = height, units = "in", res = 200)
  } else {
    pdf(paste0(fig_path, name, ".", type),
        width = width, height = height)
  }
  print(plots)
  dev.off()
}

SaveObject <- function(object, name){
  saveRDS(object, paste0(data_path, name, ".RDS"))
}

ReadObject <- function(name){
  readRDS(paste0(data_path, name, ".RDS"))
}


# Reading in data ---------------------------------------------------------

mat_path <- "~/Documents/bhutani_lab/PGDH_scRNA_Mamta/analysis/PGDH_chondro_project/data/DGE_filtered"
mat <- ReadParseBio(mat_path)

## Check to see if empty gene names are present, add name if so. 

table(rownames(mat) == "")

# No need to run if FALSE
rownames(mat)[rownames(mat) == ""] <- "unknown"

## Read in cell meta data ----

cell_meta <- read.csv(paste0(mat_path, "/cell_metadata.csv"), row.names = 1)

## Create Seurat object ----

pgdh.combined <-  CreateSeuratObject(mat, 
                                     min.features = 100, 
                                     min.cells = 100,
                                     names.field = 0, 
                                     meta.data = cell_meta)

# Setting our initial cell class to a single type, this will change
# after clustering. 

pgdh.combined@meta.data$orig.ident <- factor(rep("knee", 
                                                 nrow(pgdh.combined@meta.data)))

# Changing the active Identity to this new class name 

Idents(pgdh.combined) <- pgdh.combined@meta.data$orig.ident

## Save the object up to this point so it can be loaded quickly ----

SaveObject(pgdh.combined, "01_seurat_obj_before_QC")

## remove large files no longer needed

rm(mat)
rm(cell_meta)
