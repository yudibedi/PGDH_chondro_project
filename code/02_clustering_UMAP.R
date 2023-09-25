
# Read data in from code 01 ----------------------------------------------------

pgdh.combined <- ReadObject("01_seurat_obj_before_QC")

set.seed(42)

# Cell quality control ----------------------------------------------------

## create a new data column for percent of features from 
## mitochondrial DNA ----

pgdh.combined[["percent.mt"]] <- PercentageFeatureSet(pgdh.combined, 
                                                      pattern = "^mt-")
plot <- VlnPlot(pgdh.combined, 
                pt.size = 0.10, 
                features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                ncol = 3)

SaveFigure(plot, "vln_pre_QC", width = 12, height = 6)

plot1 <- FeatureScatter(pgdh.combined, 
                        feature1 = "nCount_RNA", 
                        feature2 = "percent.mt") + 
  NoLegend()
plot2 <- FeatureScatter(pgdh.combined, 
                        feature1 = "nCount_RNA", 
                        feature2 = "nFeature_RNA") + 
  NoLegend()
plot3 <- FeatureScatter(pgdh.combined, 
                        feature1 = "nFeature_RNA", 
                        feature2 = "percent.mt") + 
  NoLegend()
SaveFigure((plot1 + plot2 + plot3),
           "scatter_QC",
           width = 18,
           height = 6,
           res = 200)

rm(list = c("plot","plot1","plot2","plot3"))

## Summary stats for nCount, nFeature and percent.mt
summary(pgdh.combined$nCount_RNA)
summary(pgdh.combined$nFeature_RNA)
summary(pgdh.combined$percent.mt)

## Perform the filtering for outliers ----

pgdh.filtered <- subset(pgdh.combined, 
                        subset = nFeature_RNA < 5000 & 
                          nCount_RNA < 20000 & 
                          percent.mt < 10)

## Save the object up to this point so it can be loaded quickly

SaveObject(pgdh.filtered, "02_seurat_obj_post_QC")
pgdh.filtered <- ReadObject("02_seurat_obj_post_QC")

rm(pgdh.combined)
sessionInfo()

# Clustering the data using Seurat -------------------------------

## Normalizing the data ----

?NormalizeData
pgdh.filtered <- NormalizeData(pgdh.filtered, 
                               normalization.method = "LogNormalize", 
                               scale.factor = 10000)

## Identification of highly variable features ----

pgdh.filtered <- FindVariableFeatures(pgdh.filtered,
                                      selection.method = "vst",
                                      nfeatures = 2000)
# encountered a strange error here - uninstalled all seurat packages 
# in the packages side bar, reinstalled, then loaded and problem fixed itself

## Identify the 20 most highly variable genes

top20 <- head(VariableFeatures(pgdh.filtered), 20)

## Plot variable features with labels

plot1 <- VariableFeaturePlot(pgdh.filtered)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
SaveFigure(plot2, "top20_var_features", width = 7, height = 10)

## Scaling the data to apply linear transformation ----

pgdh.filtered <- ScaleData(pgdh.filtered)

## perform linear dimension reduction

pgdh.filtered <- RunPCA(pgdh.filtered)

## save the object up to this point so it can be loaded quickly

SaveObject(pgdh.filtered, "03_seurat_obj_after_PCA")

pgdh.filtered <- ReadObject("03_seurat_obj_after_PCA")

## Examine PCA results ----
## Examine and visualize PCA results a few different ways - 

#1 VizDimloading

print(pgdh.filtered[["pca"]], dims = 1:5, nfeatures = 5)

plot <- VizDimLoadings(pgdh.filtered, dims = 1:2, reduction = "pca")
SaveFigure(plot, "viz_PCA_loadings", width = 10, height = 8)

#2 DimPlot

plot <- DimPlot(pgdh.filtered, 
                reduction = "pca", 
                group.by = "orig.ident")
SaveFigure(plot, "pc1_2_scatter", width = 8, height = 6)

#3 DimHeatmap - Image doesn't save as png unless fast = FALSE

plot <- DimHeatmap(pgdh.filtered, 
                   dims = 1, 
                   cells = 500, 
                   balanced = TRUE, 
                   fast = FALSE)
SaveFigure(plot, "dim_heatmap1", width = 8, height = 6)

plot <- DimHeatmap(pgdh.filtered, 
                   dims = 1:15, 
                   cells = 500, 
                   balanced = TRUE, 
                   fast = FALSE)
SaveFigure(plot, "dim_heatmap1_15", width = 12, height = 18)

## Determine the dimensionality of the dataset ----
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time - I skipped Jackstraw
pgdh.filtered <- JackStraw(pgdh.filtered, num.replicate = 100)
pgdh.filtered <- ScoreJackStraw(pgdh.filtered, dims = 1:20)
plot <- JackStrawPlot(pgdh.filtered, dims = 1:20)
SaveFigure(plot, "PC_jackstraw_plot", width = 10, height = 8)

## alternate method is to generate elbow plot
?ElbowPlot
plot <- ElbowPlot(pgdh.filtered,ndims = 50)
SaveFigure(plot, "PC_elbow_plot", width = 8, height = 10)

## Cluster the cells ----
 
## these values were chosen to make 25 PCs and 
## tried 0.6 resolution copied from Stanford tutorial
pgdh.filtered <- FindNeighbors(pgdh.filtered, dims = 1:25)
pgdh.filtered <- FindClusters(pgdh.filtered, resolution = 0.60)

plot <- DimPlot(pgdh.filtered)
SaveFigure(plot,'dimplot_pre_cluster_tree', width = 6, height = 5)

# reordering clusters according to their similarity - 
# done for ease in merging cell clusters by reassigning each cluster position in phylogenic tree
install.packages("ape")
library(ape)

pgdh.filtered <- BuildClusterTree(pgdh.filtered, 
                                  reorder = TRUE, 
                                  reorder.numeric = TRUE)
DimPlot(pgdh.filtered)

## Non-linear dimensional reduction UMAP/tSNE ----
pgdh.filtered <- RunUMAP(pgdh.filtered, dims = 1:25)
?DimPlot
plot <- DimPlot(pgdh.filtered, 
                reduction = "umap", 
                label = TRUE) + NoLegend()
SaveFigure(plot, "dimplot_umap", width = 8, height = 8)

## Sanity check
pgdh.filtered

# number of cells per samples
per_sample_counts <- table(pgdh.filtered$sample)
write.csv(per_sample_counts,'data/sample_counts.csv')
# visualizing clustered expression of some genes
plot <- FeaturePlot(pgdh.filtered, 
                    features = c('Acan','Col2a1','Sox9','Nt5e'))
SaveFigure(plot,'chondrocyte_markers',width = 8, height = 7)
# Save output object ----
## Saving the object at this point so it can be easily loaded back in without 
## computationally intensive steps or to share with collaborators
SaveObject(pgdh.filtered, "04_seurat_obj_clustered")
