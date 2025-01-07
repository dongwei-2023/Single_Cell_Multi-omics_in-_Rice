library(Seurat)
library(monocle3)
library(dplyr)
library(ggplot2)

tissue_order <- c("Root", "ST", "leaf", "Flag", "SAM", "Bud", "SP", "Seed")
tissue.colors <- scales::hue_pal()(8)
names(tissue.colors) <- c("Bud", "Flag", "leaf", "Root", "SAM", "Seed", "SP", "ST")

fd <-
    data.frame(
        gene_short_name = rownames(counts),
        row.names = rownames(counts)
    )

cds <-
    new_cell_data_set(counts, cell_metadata = pd, gene_metadata = fd)
cds <- preprocess_cds(cds, num_dim = 100)
cds <- align_cds(cds, alignment_group = "tissue")
cds <-
    reduce_dimension(cds, umap.min_dist = 0.5, reduction_method = "UMAP", umap.n_neighbors = 100, spread = 1.5)
cds <- learn_graph(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds, close_loop = F, learn_graph_control = list(nn.k = 50, ncenter = 600))
cds <- order_cells(cds)
plot_cells(cds,
    color_cells_by = "pseudotime",
    cell_size = 0.5, trajectory_graph_segment_size = 1,
    label_branch_points = F, label_leaves = F, label_roots = T
)
plot_cells(cds,
    color_cells_by = "tissue",
    group_label_size = 6,
    label_branch_points = F, label_leaves = F, label_roots = F
)

pdf("Vascular.monocle3.pdf", 6.3, 5.5)
plot_cells(cds,
    color_cells_by = "tissue",
    label_branch_points = F,
    label_principal_points = F, rasterize = T,
    label_leaves = F, label_roots = F
) +
    scale_color_manual(values = tissue.colors) + theme(legend.position = "right")
plot_cells(cds,
    color_cells_by = "pseudotime",
    cell_size = 0.5, trajectory_graph_segment_size = 1, rasterize = T,
    label_branch_points = F, label_leaves = F, label_roots = F
) +
    scale_color_gradient(high = "red", low = "yellow")
dev.off()

save(list = ls(), file = "monocle3.Rdata")
savehistory(file = "monocle3.R")
