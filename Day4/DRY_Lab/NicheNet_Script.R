    ##### CELL-CELL INTERACTION: NicheNet ##### 
### Universit� degli Studi di Milano | Single cell analysis boot camp
### 14/09/2024

#### Prepare NicheNet analysis ####
### Load Packages (tidyverse allow you to write code in a smart way)
library(Seurat)
library(nichenetr)
library(tidyverse)
library(openxlsx)
library(circlize)

### Set the working directory
setwd("S:/NicheNet/Data")

### Read and visualize the Seurat objects
## MAIT cells
MAIT = readRDS("Mait.RDS")
head(MAIT@meta.data)

Idents(MAIT) = 'seurat_cluster'
DimPlot(MAIT, reduction = "umap", pt.size = 1.5, label = T, label.size = 5, label.box = T) + NoLegend()

## B cells
B_P3 = readRDS("B_P3.RDS")
head(B_P3@meta.data)

Idents(B_P3) = 'seurat_cluster'
DimPlot(B_P3, reduction = "umap", pt.size = 1.5, label = T, label.size = 5, label.box = T) + NoLegend()

### Read in NicheNet's ligand-target prior model and ligand-receptor network
## Ligand-receptor network: a set of potential ligands and receptors
lr_network =  readRDS("lr_network.rds")    
head(lr_network)[1:2]

## Ligand-target matrix: target genes in rows, ligands in columns
ligand_target_matrix =  readRDS("ligand_target_matrix.rds")   
head(ligand_target_matrix)[1:5, 1:5]

#### Perform NicheNet analysis ####
### Define the sender/receiver cell population and determine which genes are expressed in both populations
## Sender: MAIT cells
Idents(MAIT) = 'time'
expressed_genes_sender = get_expressed_genes("P3", MAIT, pct = 0.10)
expressed_genes_sender =  expressed_genes_sender %>% .[. %in% rownames(ligand_target_matrix)]
head(expressed_genes_sender)

## Receiver: B cells
  # Since there are several B cell clusters, to study the ligand-target interactions, each cluster has to be considered singularly.
  #   In this tutorial, we will focus on cluster 2:
B_P3_c2 = subset(B_P3, idents = '2')

Idents(B_P3_c2) = 'time'
expressed_genes_receiver =  get_expressed_genes("P3", B_P3_c2, pct = 0.10)
background_expressed_genes =  expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
head(background_expressed_genes)

### Define a gene set of interest
  # These are the genes in the receiver/target� cell population that are differently expressed between the baseline (P0/pre-vaccine)
  #   and the condition of interest (P3/three days after the vaccine boost).
  # We want to know whether these genes are potentially affected by the ligands expressed by the interacting cells (MAIT cells).
  # Differentially expressed genes (DEGs) were computed by applying the Seurat Wilcoxon statistic test (FindMarkers function)

DEG = read.xlsx("DEGs Bc2_P3vsP0.xlsx", rowNames = TRUE, colNames = TRUE)
head(DEG)
length(DEG$names)

## Only significant DEGs have to be considered in the analysis:
DEG_filter = filter(DEG, pvals_adj <= 0.05 & abs(logfoldchanges) >= 0.25  & pct_nz_group > 0.1)
head(DEG_filter)
length(DEG_filter$names)

rownames(DEG_filter) = DEG_filter$names
geneset_oi = rownames(DEG_filter)
head(geneset_oi)
length(geneset_oi)

## Among the DEGs, we have to select only the once in the ligand_target_matrix
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

### Define a gene set of potential ligands
  # These are ligands that are expressed by the sender cell population (MAIT cells) and bind a putative receptor expressed
  #   by the receiver population (B cells)
ligands =  lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

expressed_ligands = intersect(ligands, expressed_genes_sender)
expressed_receptors =  intersect(receptors, expressed_genes_receiver)

potential_ligands =  lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
head(potential_ligands)

#### Perform NicheNet ligand activity analysis ####
### Rank the potential ligands based on the presence of their target genes in the gene set of interest
  # This Pearson correlation indicates the ability of each ligand to predict differential expressed genes P3 vs P0,
  #   and better predictive ligands are thus ranked higher.

ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
ligand_activities[1:20, ]

  # The ligand activity measures how well a ligand can predict the observed differentially expressed genes compared to the background of expressed genes.
  #   Since it is a predictive analysis, we will focus on the top 20 pearson-ranked ligands.

best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
list(best_upstream_ligands)

ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)
vis_ligand_pearson = ligand_pearson_matrix[best_upstream_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("pearson")

vis_ligand_pearson %>%
  make_heatmap_ggplot("Prioritized ligands", "Ligand activity", color = "darkorange", 
                      legend_position = "top", x_axis_position = "bottom", legend_title = "Pearson coefficient")  +
  theme(axis.text.y = element_text(face = "italic"), axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0),
        legend.position = "bottom", legend.box = "horizontal") + 
  theme(text = element_text(size = 14)) +
  scale_fill_gradient2(low = "whitesmoke",  high = "darkorange")


### Infer the top-predicted target genes of ligands that are top-ranked in the ligand activity analysis
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links, geneset = geneset_oi,
                                                                 ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()
head(active_ligand_target_links_df)

## To visualize as a heatmap matrix:
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

vis_ligand_target %>%
  make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "red",legend_position = "top",
                      x_axis_position = "top", legend_title = "Regulatory potential") +
  theme(axis.text.x = element_text(face = "italic")) + theme(text = element_text(size=14)) + scale_fill_gradient2(low = "whitesmoke",  high = "purple")

## To visualize the  data as a chord diagram/circos plot:
  # To avoid making circos plots with too many ligand-target links, it's possible to show only links with a weight higher than a
  #   predefined cutoff (e.g., links belonging to the 40% (quantile(0.4)) of lowest scores were removed)

cutoff_include_all_ligands = active_ligand_target_links_df$weight %>% quantile(0.4)
active_ligand_target_links_df_circos = active_ligand_target_links_df %>% filter(weight > cutoff_include_all_ligands)

ligands_to_remove = setdiff(active_ligand_target_links_df$ligand %>% unique(), active_ligand_target_links_df_circos$ligand %>% unique())
targets_to_remove = setdiff(active_ligand_target_links_df$target %>% unique(), active_ligand_target_links_df_circos$target %>% unique())

circos_links = active_ligand_target_links_df %>% filter(!target %in% targets_to_remove &!ligand %in% ligands_to_remove)

links_circle = circos_links %>% select(ligand, target, weight)
target_order = circos_links$target %>% unique()
ligand_order = circos_links$ligand %>% unique()
order = c(ligand_order, target_order)

circos.par(gap.after = c(rep(2, length(unique(links_circle[[1]]))-1), 10, 
                         rep(2, length(unique(links_circle[[2]]))-1), 10))
chordDiagram(links_circle, directional = 1, order = order, link.sort = TRUE, link.decreasing = FALSE,
             transparency = 0, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) 








