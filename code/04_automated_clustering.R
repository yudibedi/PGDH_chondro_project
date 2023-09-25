# Automated classification using  singleR ## ----
## Predict using the MouseRNAseq data from celldex ----

BiocManager::install("celldex")
browseVignettes("celldex")
library(celldex)
library(SingleR)

mouseref <- celldex::MouseRNAseqData()
mouseref$label.main

# What is the structure of mouseref?
mouseref
View(mouseref)
View(as.data.frame(mouseref@colData))

### Predict using the Broad labels ----
sceP <- as.SingleCellExperiment(pgdh.filtered)
pred.mref <- SingleR(test = sceP, 
                     ref = mouseref, 
                     labels = mouseref$label.main)
table(pgdh.filtered$tree.ident)
?plotScoreHeatmap

plot <- plotScoreHeatmap(pred.mref,
                         max.labels = 25,
                         clusters = pgdh.filtered$seurat_clusters,
                         order.by = "clusters",
                         show_colnames = FALSE)
SaveFigure(plot, "broad_label_markers", width = 10, height = 10)

### Predict using the Fine-grained labels ----
pred.mref.fine <- SingleR(test = sceP, 
                          ref = mouseref, 
                          labels = mouseref$label.fine)
plot <- plotScoreHeatmap(pred.mref.fine,
                         max.labels = 30,
                         clusters = pgdh.filtered$seurat_clusters,
                         order.by = "clusters",
                         show_colnames = FALSE)
SaveFigure(plot, "fine_label_markers_expanded", width = 10, height = 10)

### Combine the labels for improved accuracy ----
pred.comb <- combineCommonResults(
  list(
    "Broad" = pred.mref,
    "Fine" = pred.mref.fine
  )
)

plotScoreHeatmap(pred.comb,
                 max.labels = 30,
                 clusters = pgdh.filtered$seurat_clusters,
                 order.by = "clusters",
                 show_colnames = FALSE)
?plotScoreHeatmap
table(pgdh.filtered$seurat_clusters,pgdh.filtered$tree.ident)
## umap plot of data depicting seurat_clusters as labels were generated with
## reference to these clusters
plot <- DimPlot(pgdh.filtered, 
                group.by = 'seurat_clusters', 
                label = TRUE) + 
  NoLegend()
SaveFigure(plot,'umap_by_seurat', width = 6, height = 6)

plot2 <- DimPlot(pgdh.filtered,
                 label = TRUE) +
  NoLegend()
SaveFigure((plot + plot2),'umap_seurat_tree_side', width = 12, height = 6)

## Add the predicted labels and compare to clusters ----
# You are looking for clusters which have unambiguous labels
pgdh.filtered$predicted_id <- pred.comb$pruned.labels
table(pgdh.filtered$predicted_id, pgdh.filtered$seurat_clusters)

pgdh.filtered$predicted_broad_id <- pred.mref$pruned.labels
pgdh.filtered$predicted_fine_id <- pred.mref.fine$pruned.labels

table(pgdh.filtered$predicted_fine_id, pgdh.filtered@active.ident)
table(pgdh.filtered$predicted_broad_id, pgdh.filtered@active.ident)



