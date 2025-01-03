library(Seurat)
library(ggplot2)
library(dplyr)
library(harmony)
library(ggplot2)
library(ggrastr)
library(ggrepel)
library(ggpubr)

load("rice.sub.Rdata")
obj <- rice.sub

tissue_order <- c("Root", "ST", "leaf", "Flag", "SAM", "Bud", "SP", "Seed")
tissue.colors <- scales::hue_pal()(8)
names(tissue.colors) <- c("Bud", "Flag", "leaf", "Root", "SAM", "Seed", "SP", "ST")


DefaultAssay(obj) <- "ATAC"
obj <-
    RunHarmony(
        object = obj,
        group.by.vars = "tissue",
        reduction.use = "lsi",
        assay.use = "ATAC",
        project.dim = FALSE
    )

obj[["integrated_atac"]] <-
    CreateDimReducObject(
        embeddings = obj[["harmony"]]@cell.embeddings,
        key = "integratedATAC_",
        assay = "ATAC"
    )
obj[["harmony"]] <- NULL

DefaultAssay(obj) <- "RNA"
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj, features = features)
obj <- RunPCA(obj, features = features, npcs = 50)
obj <- RunHarmony(
    object = obj, group.by.vars = "tissue",
    dims.use = 1:50,
    sigma = 10, nclust = 16, max_iter = 30,
    theta = 2, lambda = NULL,
    project.dim = FALSE
)
obj[["integrated_rna"]] <- CreateDimReducObject(
    embeddings = obj[["harmony"]]@cell.embeddings,
    key = "integratedRNA_",
    assay = "RNA"
)

obj <- FindMultiModalNeighbors(
    object = obj,
    reduction.list = list("integrated_rna", "integrated_atac"),
    dims.list = list(1:30, 2:30),
    modality.weight.name = "RNA.weight",
    verbose = TRUE
)
### Build a joint UMAP visualization

obj <- RunUMAP(
    object = obj,
    nn.name = "weighted.nn",
    assay = "RNA",
    verbose = TRUE,
    reduction.key = "wsnnUMAP_",
    reduction.name = "wsnnumap",
    n.neighbors = 100,
    n.epochs = 200,
    min.dist = 0.1,
    spread = 10
)
pdf("harmony.pdf", 12, 9)
DimPlot(obj,
    group.by = "cluster_names",
    label = T,
    repel = T
) + scale_color_manual(values = clusters_colors)
DimPlot(obj, group.by = "cluster_names") + scale_color_manual(values = clusters_colors)
dev.off()

obj$final_id2 <-
    data.frame(cluster = obj$cluster_names, tissue = obj$tissue) %>%
    mutate(scluster = paste(tissue, cluster, sep = ":")) %>%
    mutate(cluster_id = dense_rank(scluster)) %>%
    pull(cluster_id)
dd <- data.frame(
    obj@reductions$wsnnumap@cell.embeddings,
    tissue = obj$tissue,
    final_id = obj$final_id2,
    name_id = obj$cluster_names
)
sid.cor <- dd %>%
    group_by(final_id) %>%
    summarise(
        wsnnUMAP_1 = median(wsnnUMAP_1),
        wsnnUMAP_2 = median(wsnnUMAP_2)
    )
clusters_colors[grepl("unknown", names(clusters_colors))] <- "grey"

pdf("harmony.tissue.pdf", 6.5, 5.5)
DimPlot(obj, reduction = "wsnnumap", group.by = "tissue") + scale_color_manual(values = tissue.colors)
dev.off()
pdf("harmony.number.pdf", 12, 9)
ggplot(dd, aes(wsnnUMAP_1, wsnnUMAP_2)) +
    geom_point(aes(color = name_id), size = 0.1) +
    geom_text_repel(data = sid.cor, aes(label = final_id)) +
    theme_pubr(legend = "right") +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    scale_color_manual(values = clusters_colors)
dev.off()



pdf("new.harmony4.number.pdf", 12, 9)
ggplot(dd, aes(wsnnUMAP_1, wsnnUMAP_2)) +
    rasterise(geom_point(aes(color = name_id), size = 0.1, alpha = 0.4)) +
    geom_text_repel(data = sid.cor, aes(label = final_id)) +
    theme_pubr(legend = "right") +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    scale_color_manual(values = clusters_colors) +
    theme(aspect.ratio = 1)
ggplot(dd, aes(wsnnUMAP_1, wsnnUMAP_2)) +
    rasterise(geom_point(aes(color = name_id), size = 0.1, alpha = 0.2)) +
    geom_text_repel(data = sid.cor, aes(label = final_id)) +
    theme_pubr(legend = "right") +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    scale_color_manual(values = clusters_colors) +
    theme(aspect.ratio = 1)
dev.off()

data.frame(
    cell = colnames(obj), obj@reductions$wsnnumap@cell.embeddings[, 1:2],
    sample = obj$sample, tissue = obj$tissue, cluster_name = obj$cluster_names
) %>%
    left_join(data.frame(clusters_colors, cluster_name = names(clusters_colors))) %>%
    write.csv("cell.embedding.infos.csv", row.names = F)

dd %>%
    dplyr::select(tissue, final_id, name_id) %>%
    distinct() %>%
    arrange(final_id) %>%
    write.csv("new.harmony.number.csv", row.names = F)
