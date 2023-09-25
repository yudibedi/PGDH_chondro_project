# Read in data ----
pgdh.filtered <- ReadObject("04_seurat_obj_clustered")
# Determine cluster markers ----
# check currently loaded packages
sessionInfo()

# load the missing packages
library(Seurat)
library(tidyverse)

## Finding cluster markers based on diff gene expression ----
# NOTE grouping of clusters here is by top 2 genes only whereas below it is using 5
all_markers <- FindAllMarkers(pgdh.filtered, 
                              min.pct = 0.25, 
                              logfc.threshold = 0.25)

all_markers

## top 2 for each cluster
all_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)


## Visualizing top n genes per cluster in a single figure ----
# DoHeatmap is a good tool for few thousand cell but loses resolution as number of cells increase
# DotPlot addresses this by displaying average expression for all cells in a cluster

top5 <- all_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
to_plot <- unique(top5$gene)

plot <- DotPlot(pgdh.filtered, 
                features = to_plot, 
                group.by = "tree.ident") + 
  coord_flip()
SaveFigure(plot, "dplot_top5", width = 9, height = 20)

# same as above but for top 2 genes
top2_genes <- all_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
to_plot <- unique(top2_genes$gene)

plot <- DotPlot(pgdh.filtered, 
                features = to_plot, 
                group.by = "tree.ident")
coord_flip()
SaveFigure(plot, "dplot_top2", width = 9, height = 10)

# Heatmap plot of same

plot <- DoHeatmap(pgdh.filtered, features = top5$gene, ) + 
  FontSize(y.text = 8) +
  NoLegend()
SaveFigure(plot, "heatmap_top5", width = 25, height = 20)

# Violin plot of top 2 cluster markers

plot <- VlnPlot(pgdh.filtered, 
                features = to_plot, 
                group.by = "tree.ident", 
                stack = TRUE) + NoLegend()
SaveFigure(plot, "vln_exp4", width = 16, height = 8)

# Looking at some known chondrocyte and immune cell markers
plot<- VlnPlot(pgdh.filtered, 
        features = c("Acan","Col2a1","Mmp13","Sox9","Nt5e","Cd200","Ptprc","Pax5"), 
        stack = TRUE, 
        flip = TRUE) +
  NoLegend()
SaveFigure(plot, 'random_markers', width = 5, height =5)

FeaturePlot(pgdh.filtered, features = c("Dcn","Prg4","Cox5a",))
FeaturePlot(pgdh.filtered, features = "Cd200", split.by = 'age_inj', ncol = 2)


## exporting top 5 markers per cluster as a csv file ----
write.csv(top5,'data/top5_markers.csv')


# Adding in metadata for age and SW treatment ----

Idents(pgdh.filtered)

factor(pgdh.filtered$sample)

## set the identitites to the sample names

Idents(pgdh.filtered) <- pgdh.filtered$sample

## check identities

Idents(pgdh.filtered)

## create new ids vectors for age and inj status ----

age_ids <- c(
  "C10M1" = "Young",
  "C10M2" = "Young",
  "C10M3" = "Young",
  "C10M4" = "Young",
  "C1M1" = "Aged",
  "C1M2" = "Aged",
  "C2M1" = "Aged",
  "C3M1" = "Aged",
  "C4M2" = "Aged",
  "C4M5" = "Aged",
  "C5M2" = "Aged",
  "C5M5" = "Aged",
  "C9M1" = "Young",
  "C9M2" = "Young",
  "C9M3" = "Young",
  "C9M4" = "Young"
)

inj_ids <- c(
  "C10M1" = "Veh",
  "C10M2" = "Veh",
  "C10M3" = "SW",
  "C10M4" = "SW",
  "C1M1" = "SW",
  "C1M2" = "Veh",
  "C2M1" = "Veh",
  "C3M1" = "SW",
  "C4M2" = "SW",
  "C4M5" = "SW",
  "C5M2" = "Veh",
  "C5M5" = "Veh",
  "C9M1" = "Veh",
  "C9M2" = "Veh",
  "C9M3" = "SW",
  "C9M4" = "SW"
)

## Now Rename identities to the age vector 
## and make this into a new metadata column called 'age'

pgdh.filtered <- RenameIdents(pgdh.filtered, age_ids)
pgdh.filtered$age <- Idents(pgdh.filtered) 

plot <- DimPlot(pgdh.filtered)
SaveFigure(plot, 'umap_by_age', width = 7, height = 6)

## Re-run the code to set identities as sample names

## Now Rename identities to the inj vector 
## and make this into a new metadata column called 'inj'

pgdh.filtered <- RenameIdents(pgdh.filtered, inj_ids)
pgdh.filtered$inj <- Idents(pgdh.filtered) 

plot <- DimPlot(pgdh.filtered, cols = c("lightgreen","darkorchid2"))
SaveFigure(plot, 'umap_by_inj', width = 7, height = 6)

## creating metadata column from combined info ----
## Based off of these two metadata, we can create a 3rd column combining both info
pgdh.filtered$age_inj <- paste(pgdh.filtered$age, 
                               pgdh.filtered$inj, 
                               sep = "_")
## Re-set the identities to tree.ident
Idents(pgdh.filtered) <- pgdh.filtered$tree.ident
DimPlot(pgdh.filtered,
        label = TRUE)
DimPlot(pgdh.filtered,
        group.by = 'age',
        split.by = 'inj')

## check the levels of the factors of each metadata
factor(pgdh.filtered$age)
factor(pgdh.filtered$inj)
factor(pgdh.filtered$age_inj)

## looks like the levels for age_inj need re-rodering
pgdh.filtered$age_inj <- factor(pgdh.filtered$age_inj, 
                                levels = c("Young_Veh", "Young_SW", "Aged_Veh", "Aged_SW"))

## plotting of all clusters split by groups ----
Idents(pgdh.filtered) <- pgdh.filtered$tree.ident
plot1 <- DimPlot(pgdh.filtered,
                 group.by = 'tree.ident',
                 split.by = 'age_inj',
                 ncol = 2) 
plot2 <- DimPlot(pgdh.filtered,
                 group.by = 'tree.ident',
                 label = TRUE) + 
  NoLegend()
SaveFigure((plot1+plot2),
           'split_by_age_inj_with_labels', 
           width = 8, 
           height = 15)
table(pgdh.filtered$age, pgdh.filtered$inj)

DimPlot(pgdh.filtered,
        group.by = "tree.ident")

# Identifying cluster cell types manually ----
library(scCustomize)
library(Nebulosa)
library(viridis)
library(qs)
library(patchwork)

Idents(pgdh.filtered) <- pgdh.filtered$tree.ident
plot <- Clustered_DotPlot(pgdh.filtered,
                          features = top5$gene,
                          k=13,
                          flip = T,
                          column_label_size = 10,
                          row_label_size =10)
SaveFigure (plot[[2]], 
            'clustered_top5_all_cells',
            width = 20,
            height = 5)
FeaturePlot(pgdh.filtered, features = "Lepr")

## Meniscus markers from CellMarker2.0 [Cd105 = Eng, Cd34, Cd73 = Nt5e, Cd90 = Thy1, Sca-1= Ly6a]
meniscus_markers <- c("Eng", "Cd34", "Cd44", "Nt5e", "Thy1", "Ly6a")
FeaturePlot(pgdh.filtered,
            features = meniscus_markers,
            ncol = 3)
VlnPlot(pgdh.filtered,
        group.by = 'tree.ident',
        features = meniscus_markers,
        stack = TRUE,
        flip = TRUE)

VlnPlot(pgdh.filtered,
        group.by = 'tree.ident',
        features = c("Tie1","Apod","Tagln","Myl9"),
        stack = TRUE,
        flip = TRUE)

Plot_Density_Joint_Only(pgdh.filtered,
                        features = chondrocyte_markers)

## Lets make marker gene lists based on the Loots paper

immune_markers <- c("Ptprc")
blood_markers <- c("Hemgn")
chondrocyte_markers <- c("Sox9","Col2a1","Acan")
cellcycle_markers <- c("Mki67","Cdk1","Stmn1","Top2a","Cenpa")
syno_subintfibro_markers <- c("Cxcl12","Col3a1","Col14a1")
syno_intfibro_markers <- c("Prg4","Has1","Htra1")
mesen_fibro_markers <-c("Pdgfra","Pdpn","Clec3b","Abi3bp","S100a4","Thy1")
cytokines <- c("Ccl2","Ccl7","C3","C4b")
osteoblast_markers <- c("Col1a1","Col1a2","Bglap","Alpl")
endo_markers <- c("Pecam1","Ptprb","Cdh5")
pericyte_markers <- c("Rgs5","Myh11","Mcam")

## Combine marker genes from all lists
loots_genes <- c(immune_markers, 
                 blood_markers,
                 chondrocyte_markers,
                 cellcycle_markers,
                 syno_intfibro_markers,
                 syno_subintfibro_markers,
                 mesen_fibro_markers,
                 cytokines,
                 osteoblast_markers,
                 endo_markers,
                 pericyte_markers)

## Create a metadata column indicating the list name for each gene
## - this DID NOT WORK
#list_assignment <- character(nrow(pgdh.filtered))
#names(list_assignment) <- rownames(pgdh.filtered)

#list_assignment[rownames(pgdh.filtered) %in% immune_markers] <- "Immune"
#list_assignment[rownames(pgdh.filtered) %in% blood_markers] <- "Blood"
#list_assignment[rownames(pgdh.filtered) %in% chondrocyte_markers] <- "Chondrocyte"
#list_assignment[rownames(pgdh.filtered) %in% cellcycle_markers] <- "Cell Cycle"
#list_assignment[rownames(pgdh.filtered) %in% syno_intfibro_markers] <- "Synovial Intimal Fibroblast"
#list_assignment[rownames(pgdh.filtered) %in% syno_subintfibro_markers] <- "Synovial Subintimal Fibroblast"
#list_assignment[rownames(pgdh.filtered) %in% mesen_fibro_markers] <- "Mesenchymal or Fibroblast"
#list_assignment[rownames(pgdh.filtered) %in% osteoblast_markers] <- "Osteoblast"
#list_assignment[rownames(pgdh.filtered) %in% endo_markers] <- "Endothelial"
#list_assignment[rownames(pgdh.filtered) %in% pericyte_markers] <- "Pericyte"
#list_assignment[!rownames(pgdh.filtered) %in% loots_genes] <- "Other"

## Assign metadata to object
#pgdh.filtered$marker <- list_assignment


## plotting the genes listed above

plot <- VlnPlot(pgdh.filtered,
                features = loots_genes,
                stack = TRUE,
                flip = TRUE,
                group.by = 'tree.ident',
                ) + NoLegend()
SaveFigure(plot, 'Loots_markers', width = 6, height = 10)

# other misc plots below

VlnPlot(pgdh.filtered, features = "Acta2")

VlnPlot(pgdh.filtered,
        idents = c('9','10','19','20','21','22'),
        group.by = 'tree.ident',
        split.by = 'age_inj',
        features = c("Cd200", "Hpgd", "Sox9"),
        stack = TRUE,
        flip = TRUE)

VlnPlot(pgdh.filtered,
        group.by = 'tree.ident',
        features = c(chondrocyte_markers,osteoblast_markers),
        stack = TRUE,
        flip = TRUE) + scale_fill_manual(values = c("red","blue"))

write.csv(table(pgdh.filtered$age_inj,pgdh.filtered$tree.ident),'data/group_by_cluster.csv')

FeaturePlot(pgdh.filtered, features = "Hpgd")

VlnPlot(pgdh.filtered,
        idents = c('9','10','19','20','21','22'),
        group.by = 'tree.ident',
        split.by = 'age_inj',
        features = c("Ndufa4","Atp5a1","Ndufs7","Cox5a","mt-Nd5", "Htra1"),
        stack = TRUE,
        flip = TRUE)

FeaturePlot(pgdh.filtered,
            features = c("Dcn", "Creb5", "Htra1","Prg4"))

# Hematopoeitic SC markers in mouse
FeaturePlot(pgdh.filtered, features = c("Mmrn1","Tie1","Trap1","Ctsk","Mmp9"))
VlnPlot(pgdh.filtered,
        group.by = 'tree.ident',
        features = c("Slamf1","Mmrn1","Tie1","Ly6a","Cd27","Meis1",
                     "Mpl","Vwf","Hlf","Gfi1","Hoxa9"),
        stack = TRUE,
        flip = TRUE)
# Neutrophil markers
neutrophil_markers <- c("S100a8","S100a9","Ltf","Itgam")
plot <- FeaturePlot(pgdh.filtered, 
            features = neutrophil_markers) 
  
VlnPlot(pgdh.filtered,
        group.by = 'tree.ident',
        features = neutrophil_markers,
        stack = TRUE,
        flip = TRUE) + 
  NoLegend() +
  labs(title = "Neutrophil Markers")

RidgePlot(pgdh.filtered,
          features = neutrophil_markers,
          ncol = 2,
          group.by = 'tree.ident')

RidgePlot(pgdh.filtered,
          features = neutrophil_markers,
          group.by = 'tree.ident',
          same.y.lims = T,
          combine = T,
          log = T, 
          ncol = 4,
          sort = 'decreasing')

FeaturePlot(chondo,
            features = "Ptges")

## Density plot for Hpgd in individual age_inj groups
Plot_Density_Custom(subset(x = chondo, age_inj == 'Young_Veh'),
            features = c("Hpgd"))
Plot_Density_Custom(subset(x = chondo, age_inj == 'Young_SW'),
                      features = c("Hpgd"))
Plot_Density_Custom(subset(x = chondo, age_inj == 'Aged_Veh'),
                      features = c("Hpgd"))
Plot_Density_Custom(subset(x = chondo, age_inj == 'Aged_SW'),
                      features = c("Hpgd"))

### Lets loop this density plot function
density_plot_list <- list()
age_inj_levels <- levels(chondo$age_inj)
goi <- "Thbs1"
for (i in age_inj_levels) {  i
  density_plot_list[[i]] <- Plot_Density_Custom(subset(x = chondo, age_inj == i),
                                                features = goi) + 
    ggtitle(label = goi, subtitle = i)
}

cowplot::plot_grid(plotlist = density_plot_list, labels = "AUTO",
                     label_size = 18, align = "vh") 
ggsave(filename = paste0("figures/density_plots/",goi,"_densityplot_by_group.png"), 
         width = 10, height = 10, limitsize = FALSE)
  
## an alternative is to run facet_grid
Plot_Density_Custom(chondo, goi) + 
  facet_grid(.~chondo$age_inj)

# Skeletal stem and progenitor cell markers
SSPC_markers <- c("Pdgfra","Pdgfrb","Cd200",
                  "Ly6a", "Lepr","Itgav",
                  "Grem1","Cxcl12","Eng")
plot <- FeaturePlot(pgdh.filtered, 
                    features = SSPC_markers,
                    ncol = 3) 
SaveFigure(plot, 'SSPC_markers_featureplot', width = 8, height = 8)

plot <- VlnPlot(pgdh.filtered,
                group.by = 'tree.ident',
                features = SSPC_markers,
                stack = TRUE,
                flip = TRUE) + 
          NoLegend() +
          labs(title = "SSPC Markers")
SaveFigure(plot, 'SSPC_markers_violinplot', width = 5, height = 5)

sublining_fibroblast_markers <- c("Cadm1","Mfap2","Dkk3","Ogn")
stromal_cell_markers <- c("Lepr","Kitl","Ebf3")
osteoclast_markers <- c("Acp5","Atp6v0d2","Ctsk","Mmp9","Ppargc1b")

FeaturePlot(pgdh.filtered,
            features = stromal_cell_markers,
            ncol = 3)

VlnPlot(pgdh.filtered, 
        group.by = 'tree.ident',
        stack = TRUE,
        flip = TRUE,
        features = osteoclast_markers)

DimPlot(pgdh.filtered, split.by = 'tree.ident')
view(as.data.frame(sigmarkers10))
VlnPlot(pgdh.filtered, 
        group.by = 'tree.ident',
        stack = TRUE,
        flip = TRUE,
        features = c("Kcna1","Mcam","Lyz2","Prg4","Pdgfra","Rspo2")) +
  NoLegend()

# Examining mitochondrial DNA genes

FeaturePlot(pgdh.filtered, 
            features = c("mt-Nd5"),
            split.by = 'age_inj')

DotPlot(pgdh.filtered,
        features = c("mt-Nd1","Il31ra","mt-Nd2","mt-Nd4","mt-Nd5", "Ndufa4"),
        group.by = 'tree.ident') + coord_flip()

table(pgdh.filtered$tree.ident == '16', pgdh.filtered$age_inj )
FeatureScatter(pgdh.filtered, 
               group.by = 'tree.ident',jitter = TRUE,
               feature1 = "nFeature_RNA", 
               feature2 = "percent.mt")
