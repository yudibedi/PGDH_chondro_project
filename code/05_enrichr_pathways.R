# Cluster Marker Analysis  ----

library(ggpubr)
library(tidyverse)

## First, re-assign idents as the clusters
Idents(pgdh.filtered) <- pgdh.filtered$tree.ident

## Part 1: Pathway enrichment using enrichr ----
## can be done online in enrichr: https://maayanlab.cloud/Enrichr/ 

### Find all positive markers only for all clusters ----
all_pos_markers <- FindAllMarkers(pgdh.filtered, 
                                  only.pos = TRUE, 
                                  min.pct = 0.25, 
                                  logfc.threshold = 0.5)

sig_markers <- all_pos_markers %>% 
  filter(p_val_adj < .05)

write_csv(sig_markers, file = "data/sig_markers.csv")

library(enrichR)
dbs <- listEnrichrDbs()
to_check_cell_type <- c("Mouse_Gene_Atlas",
                        "PanglaoDB_Augmented_2021",
                        "CellMarker_Augmented_2021",
                        "Tabula_Muris")
to_check_pathways <- c("GO_Biological_Process_2023",
                       "KEGG_2019_Mouse",
                       "MSigDB_Hallmark_2020",
                       "Reactome_2022")
to_check_aging <- c("GTEx_Aging_Signatures_2021")
to_check_chip <- c("ChEA_2022",
                   "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")

## How many sig genes in each cluster
table(sig_markers$cluster)

### Plotting function for cluster enrichment results ----

plot_eres <- function(eres_name, eres_list, n = 10) {
  eres <- eres_list[[eres_name]]
  eres %>%
    top_n(n = n, wt = -log10(Adjusted.P.value)) %>%
    arrange(-log10(Adjusted.P.value)) %>%
    mutate(Term = factor(Term, levels = Term)) %>%
    ggplot(mapping = aes(x = Term, y = -log10(Adjusted.P.value), fill = Combined.Score)) +
    geom_bar(stat = "identity") +
    ggpubr::rotate() +
    theme_bw(base_size = 16) +
    rremove("ylab") +
    labs(title = eres_name)
}
### Using FOR loops to run analyses on the significant markers for each cluster ----
levels(factor(pgdh.filtered$tree.ident))
### pull out markers for each cluster ----
for(i in levels(factor(pgdh.filtered$tree.ident))) {
  assign(paste0("sigmarkers_",i), (sig_markers %>%
    filter(cluster == i) %>%
    pull(gene)))
}
## put them in a list 
sigmarkers <- list()
for(i in levels(factor(pgdh.filtered$tree.ident))) {
  sigmarkers[[i]] <- (sig_markers %>%
                        filter(cluster == i) %>%
                        pull(gene))
}

### Cell type analysis ----
celltype = list()
for(i in levels(factor(pgdh.filtered$tree.ident))) {
  celltype[[i]] <- enrichr(sigmarkers[[i]], databases = to_check_cell_type)
}
## Get the top hits from each and plot them
plotlist_cell_type <- list()
for (i in levels(factor(pgdh.filtered$tree.ident))) {
  plotlist_cell_type[[i]] <- lapply(names(celltype[[i]]), 
                                    plot_eres, 
                                    eres_list = celltype[[i]])
}
for (i in levels(factor(pgdh.filtered$tree.ident))) {
cowplot::plot_grid(plotlist = plotlist_cell_type[[i]], labels = "AUTO",
                   label_size = 18, align = "vh") +
  patchwork::plot_annotation(title = paste0("Cluster ",i," cell type"), 
                             theme = theme(title = element_text(size = 20)))
  ggsave(filename = paste0("figures/enrichR/cell_type/cluster",i,"_cell_type.png"), 
         width = 25, height = 10, limitsize = FALSE)
}
### Pathway analysis on chondrocyte clusters ----
#### Pull the Cluster 19 genes as a vector - should have 213 genes ----
cluster19_genes <- sig_markers %>%
  filter(cluster == "19") %>%
  pull(gene)
## Run through enrichr
cluster19_eresList <- enrichr(cluster19_genes, databases = to_check_pathways)
## Get the top hits from each and plot them
plotList <- lapply(names(cluster19_eresList), plot_eres, eres_list = cluster19_eresList)
cowplot::plot_grid(plotlist = plotList, labels = "AUTO",
                   label_size = 18, align = "vh") +
  patchwork::plot_annotation(title = "Cluster 19 Pathway Analysis", 
                             theme = theme(title = element_text(size = 20))) +
  ggsave(filename = "cluster19_enrichr_pathways.png", width = 25, height = 10, limitsize = FALSE)

#### Pull the Cluster 20 genes as a vector - should have 487 genes ----
cluster20_genes <- sig_markers %>%
  filter(cluster == "20") %>%
  pull(gene)
## Run through enrichr
cluster20_eresList <- enrichr(cluster20_genes, databases = to_check_pathways)
## Get the top hits from each and plot them
plotList <- lapply(names(cluster20_eresList), 
                   plot_eres, 
                   eres_list = cluster20_eresList)
cowplot::plot_grid(plotlist = plotList, labels = "AUTO",
                   label_size = 18, align = "vh") +
  patchwork::plot_annotation(title = "Cluster 20 Pathway Analysis", 
                             theme = theme(title = element_text(size = 20))) +
  ggsave(filename = "cluster20_enrichr_pathways.png", width = 25, height = 10, limitsize = FALSE)

#### Pull the Cluster 21 genes as a vector - should have 432 genes ----
cluster21_genes <- sig_markers %>%
  filter(cluster == "21") %>%
  pull(gene)
## Run through enrichr
cluster21_eresList <- enrichr(cluster21_genes, databases = to_check_pathways)
## Get the top hits from each and plot them
plotList <- lapply(names(cluster21_eresList), 
                   plot_eres, 
                   eres_list = cluster21_eresList)
cowplot::plot_grid(plotlist = plotList, labels = "AUTO",
                   label_size = 18, align = "vh") +
  patchwork::plot_annotation(title = "Cluster 21 Pathway Analysis", 
                             theme = theme(title = element_text(size = 20))) +
  ggsave(filename = "cluster21_enrichr_pathways.png", width = 25, height = 10, limitsize = FALSE)

#### Pull the Cluster 22 genes as a vector - should have 257 genes ----
cluster22_genes <- sig_markers %>%
  filter(cluster == "22") %>%
  pull(gene)
## Run through enrichr
cluster22_eresList <- enrichr(cluster22_genes, databases = to_check_pathways)
## Get the top hits from each and plot them
plotList <- lapply(names(cluster22_eresList), 
                   plot_eres, 
                   eres_list = cluster22_eresList)
cowplot::plot_grid(plotlist = plotList, labels = "AUTO",
                   label_size = 18, align = "vh") +
  patchwork::plot_annotation(title = "Cluster 22 Pathway Analysis", 
                             theme = theme(title = element_text(size = 20))) +
  ggsave(filename = "cluster22_enrichr_pathways.png", width = 25, height = 10, limitsize = FALSE)

## Part 2: Finding differences using genes ----
# not ideal method

# Find differences in marker genes
cluster_gene_list <- lapply(c("19", "20"), function(cluster_now) {
  sig_markers %>%
    filter(cluster == cluster_now) %>%
    pull(gene)
})
names(cluster_gene_list) <- c("Cluster: 19", "Cluster: 20")

# Use Venn Diagram to Compare
library(VennDiagram)
venn.diagram(cluster_gene_list)
venn.diagram(cluster_gene_list, filename = "figures/venn_clusters_19_20_compared.png", 
             fill = c("firebrick", "skyblue"), margin = .05)
overlap_list <- calculate.overlap(cluster_gene_list)

## Part 3: Use a differential model to find differences in biology ----
# Find markers in Cluster x VS Cluster y
cluster7vs14.markers <- FindMarkers(pgdh.filtered,
                                     ident.1 = 7, 
                                     ident.2 = 14,
                                     logfc.threshold = .25,
                                     min.pct = 0.25)

# Pull the Cluster x vs y genes as a vector
upreg_7v14 <- cluster7vs14.markers %>%
  filter(avg_log2FC > 0) %>%
  rownames_to_column("gene") %>%
  pull(gene)
down_7v14 <- cluster7vs14.markers %>%
  filter(avg_log2FC < 0) %>%
  rownames_to_column("gene") %>%
  pull(gene)
# Run through enrichr
upreg_7v14_list <- enrichr(upreg_7v14, databases = to_check_cell_type)
down_7v14_list <- enrichr(down_7v14, databases = to_check_cell_type)

# Get the top hits from each and plot them
plotList <- lapply(names(upreg_7v14_list), plot_eres, eres_list = upreg_7v14_list)
cowplot::plot_grid(plotlist = plotList, labels = "AUTO",
                   label_size = 18, align = "vh") +
  patchwork::plot_annotation(title = "Upregulated 7 vs 14", 
                             theme = theme(title = element_text(size = 20)))
  ggsave(filename = "figures/enrichR/cell_type/up_cluster7v14_enrichr_pathways.png", 
         width = 25, 
         height = 30)

plotList <- lapply(names(down_7v14_list), plot_eres, eres_list = down_7v14_list)
cowplot::plot_grid(plotlist = plotList, labels = "AUTO",
                   label_size = 18, align = "vh") +
  patchwork::plot_annotation(title = "Downregulated 7 vs 14", 
                             theme = theme(title = element_text(size = 20)))
  ggsave(filename = "figures/enrichR/cell_type/down_cluster7v14_enrichr_pathways.png", 
         width = 25, 
         height = 30)

