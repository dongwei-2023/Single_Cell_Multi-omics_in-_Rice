library(Seurat)
library(monocle)
library(dplyr)

load("leaf_Flag.counts.Rdata")

combined <- CreateSeuratObject(counts, meta.data = meta.data)
combined$group <- gsub("_.+", "", combined$orig.ident)
combined <- NormalizeData(combined)
table(combined$cluster)
Idents(combined) <- combined$cluster
sb <- subset(combined, idents = "Epidermis")
sb <- ScaleData(sb)
sb <- FindVariableFeatures(sb)
sb <- RunPCA(sb, npcs = 50)
sb <- RunUMAP(sb, dims = 1:50)

sb <- FindNeighbors(sb, dims = 1:50)
sb <- FindClusters(sb, resolution = 0.5)
DimPlot(sb, group.by = "orig.ident")

ifnb.list <- SplitObject(sb, split.by = "orig.ident")
ifnb.list <- lapply(X = ifnb.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
leaf_Flag.anchors <- FindIntegrationAnchors(
    object.list = ifnb.list, normalization.method = "SCT",
    anchor.features = features
)
leaf_Flag.combined.sct <- IntegrateData(anchorset = leaf_Flag.anchors, normalization.method = "SCT")
leaf_Flag.combined.sct <- RunPCA(leaf_Flag.combined.sct, verbose = FALSE)
leaf_Flag.combined.sct <- RunUMAP(leaf_Flag.combined.sct, reduction = "pca", dims = 1:30)
leaf_Flag.combined.sct <- FindNeighbors(leaf_Flag.combined.sct, dims = 1:30)
leaf_Flag.combined.sct <- FindClusters(leaf_Flag.combined.sct, resolution = 0.3)

DimPlot(leaf_Flag.combined.sct)
DimPlot(leaf_Flag.combined.sct, group.by = "orig.ident")

DefaultAssay(leaf_Flag.combined.sct) <- "RNA"
leaf_Flag.combined.sct <- NormalizeData(leaf_Flag.combined.sct)
leaf_Flag.combined.sct <- ScaleData(leaf_Flag.combined.sct)
allmarkers <- FindAllMarkers(leaf_Flag.combined.sct, only.pos = T, min.pct = 0.5)
pdf("epidermis.combined.UMAP.pdf", 7, 6.5)
DimPlot(leaf_Flag.combined.sct, label = T, label.size = 8)
DimPlot(leaf_Flag.combined.sct, group.by = "orig.ident")
dev.off()

source("/enrichment.R")

for (i in unique(allmarkers$cluster)) {
    anno <- get_annotation(allmarkers %>% filter(cluster %in% i) %>% pull(gene))
    pdf(paste0("c", i, ".annotation.pdf"))
    for (xx in names(anno)) {
        pdata <- anno[[xx]] %>%
            mutate(Description = gsub("(.{1,90})(\\s|$)", "\\1\n", Description)) %>%
            mutate(ratio = Count / as.numeric(gsub(".+/", "", GeneRatio))) %>%
            mutate(Description = forcats::fct_reorder(Description, ratio)) %>%
            top_n(10, ratio)
        p <- ggplot(pdata, aes(Description, ratio, color = -log10(pvalue), size = Count)) +
            geom_point() +
            coord_flip() +
            theme_bw() +
            ggtitle(xx) +
            scale_color_gradient(low = "grey", high = "red")
        print(p)
    }
    dev.off()
}

sel.gene <- unique(allmarkers$gene)

cds <- as.CellDataSet(sb)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- monocle::setOrderingFilter(cds, sel.gene)
cds <- monocle::reduceDimension(cds, method = "DDRTree", )
cds <- monocle::orderCells(cds)
pData(cds)$cluster <- Idents(leaf_Flag.combined.sct)
pData(cds)$group <- gsub("_.+", "", pData(cds)$orig.ident)

pdf("pseudotime.pdf")
monocle::plot_cell_trajectory(cds, color_by = "cluster")
monocle::plot_cell_trajectory(cds, color_by = "State")
monocle::plot_cell_trajectory(cds, color_by = "group")
dev.off()


leaf_Flag.combined.sct$State <- pData(cds)$State
DimPlot(leaf_Flag.combined.sct, group.by = "State")
Idents(leaf_Flag.combined.sct) <- leaf_Flag.combined.sct$State
State.markers <- FindAllMarkers(leaf_Flag.combined.sct, only.pos = T, min.pct = 0.5)

for (i in unique(State.markers$cluster)) {
    anno <- get_annotation(State.markers %>% filter(cluster %in% i) %>% pull(gene))
    pdf(paste0("State.", i, ".annotation.pdf"))
    for (xx in names(anno)) {
        pdata <- anno[[xx]] %>%
            mutate(Description = gsub("(.{1,90})(\\s|$)", "\\1\n", Description)) %>%
            mutate(ratio = Count / as.numeric(gsub(".+/", "", GeneRatio))) %>%
            mutate(Description = forcats::fct_reorder(Description, ratio)) %>%
            top_n(10, ratio)
        p <- ggplot(pdata, aes(Description, ratio, color = -log10(pvalue), size = Count)) +
            geom_point() +
            coord_flip() +
            theme_bw() +
            ggtitle(xx) +
            scale_color_gradient(low = "grey", high = "red")
        print(p)
    }
    dev.off()
}
save(cds, leaf_Flag.combined.sct, sb, file = "epidermis.RData")



load("/alldata.Rdata")
Flag <- alldata$Flag
leaf <- alldata$leaf

ifnb.list <- list(Flag, leaf)
ifnb.list <- lapply(X = ifnb.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
leaf_Flag.anchors <- FindIntegrationAnchors(
    object.list = ifnb.list, normalization.method = "SCT",
    anchor.features = features
)
leaf_Flag.combined.sct <- IntegrateData(anchorset = leaf_Flag.anchors, normalization.method = "SCT")
leaf_Flag.combined.sct <- RunPCA(leaf_Flag.combined.sct, verbose = FALSE)
leaf_Flag.combined.sct <- RunUMAP(leaf_Flag.combined.sct, reduction = "pca", dims = 1:30)
DimPlot(leaf_Flag.combined.sct)
leaf_Flag.combined.sct <- FindNeighbors(leaf_Flag.combined.sct, dims = 1:30)
leaf_Flag.combined.sct <- FindClusters(leaf_Flag.combined.sct, resolution = 0.3)
DimPlot(leaf_Flag.combined.sct)
DimPlot(leaf_Flag.combined.sct, group.by = "orig.ident")

DefaultAssay(leaf_Flag.combined.sct) <- "RNA"
leaf_Flag.combined.sct <- NormalizeData(leaf_Flag.combined.sct)
leaf_Flag.combined.sct <- ScaleData(leaf_Flag.combined.sct, features = rownames(leaf_Flag.combined.sct))
Idents(leaf_Flag.combined.sct) <- leaf_Flag.combined.sct$cluster
allmarkers <- parallel::mclapply(sort(unique(leaf_Flag.combined.sct$cluster)), function(x) {
    diff <- FindMarkers(leaf_Flag.combined.sct, ident.1 = x, min.pct = 0.3, only.pos = T)
    return(data.frame(diff, cluster = x, gene = rownames(diff)))
}, mc.cores = 10) %>% do.call(rbind, .)


for (i in unique(allmarkers$cluster)) {
    anno <- get_annotation(allmarkers %>% filter(cluster %in% i) %>% pull(gene))
    do.call(rbind, anno) %>% write.csv(paste0("leaf_Flag.", i, ".annotation.csv"), row.names = F)
    pdf(paste0("leaf_Flag.", i, ".annotation.pdf"), 8, 6)
    print(plot_anno(anno, type = "BP", style = "bar"))
    print(plot_anno(anno, type = "MF", style = "bar"))
    print(plot_anno(anno, type = "kegg", style = "bar"))
    dev.off()
}

table(leaf_Flag.combined.sct$cluster, leaf_Flag.combined.sct$tissue) %>%
    data.frame() %>%
    group_by(Var2) %>%
    mutate(percent = 100 * Freq / sum(Freq)) %>%
    write.csv("celltype.percent.csv", row.names = F)
p <- table(leaf_Flag.combined.sct$cluster, leaf_Flag.combined.sct$tissue) %>%
    data.frame() %>%
    group_by(Var2) %>%
    mutate(percent = 100 * Freq / sum(Freq)) %>%
    ggplot(aes(Var2, percent, fill = Var1)) +
    geom_histogram(stat = "identity") +
    theme_classic() +
    scale_fill_npg() +
    BoldTheme()
ggsave(p, filename = "celltype.percent.pdf", width = 8, height = 7)

pdf("leaf_Flag.combined.UMAP.pdf", 7, 6.5)
DimPlot(leaf_Flag.combined.sct, label = T, label.size = 6, repel = T) + NoLegend()
DimPlot(leaf_Flag.combined.sct, group.by = "tissue")
DimPlot(leaf_Flag.combined.sct, group.by = "orig.ident")
dev.off()

MS.diff <- FindMarkers(leaf_Flag.combined.sct, ident.1 = "Flag", ident.2 = "leaf", subset.ident = "Mesophyll", group.by = "tissue")
VS.diff <- FindMarkers(leaf_Flag.combined.sct, ident.1 = "Flag", ident.2 = "leaf", subset.ident = "Vascular_cylinder", group.by = "tissue")
Epi.diff <- FindMarkers(leaf_Flag.combined.sct, ident.1 = "Flag", ident.2 = "leaf", subset.ident = "Epidermis", group.by = "tissue")
MS.diff <- MS.diff %>% data.frame(., gene = rownames(.))
VS.diff <- VS.diff %>% data.frame(., gene = rownames(.))
Epi.diff <- Epi.diff %>% data.frame(., gene = rownames(.))

MS.up.ann <- get_annotation(MS.diff %>% filter(avg_log2FC > 0) %>% pull(gene))
MS.down.ann <- get_annotation(MS.diff %>% filter(avg_log2FC < 0) %>% pull(gene))
VS.up.ann <- get_annotation(VS.diff %>% filter(avg_log2FC > 0) %>% pull(gene))
VS.down.ann <- get_annotation(VS.diff %>% filter(avg_log2FC < 0) %>% pull(gene))
Epi.up.ann <- get_annotation(Epi.diff %>% filter(avg_log2FC > 0) %>% pull(gene))
Epi.down.ann <- get_annotation(Epi.diff %>% filter(avg_log2FC < 0) %>% pull(gene))

pdf("Epi.diff.ann.pdf", 9, 10)
plot_compare_anno(
    Epi.up.ann, Epi.down.ann,
    lab.1 = "Flag.up", lab.2 = "Flag.down",
    type = "BP", stype = "bar"
)
plot_compare_anno(
    Epi.up.ann, Epi.down.ann,
    lab.1 = "Flag.up", lab.2 = "Flag.down",
    type = "MF", stype = "bar"
)
plot_compare_anno(
    Epi.up.ann, Epi.down.ann,
    lab.1 = "Flag.up", lab.2 = "Flag.down",
    type = "kegg", stype = "bar"
)
dev.off()

pdf("MS.diff.ann.pdf", 9, 10)
plot_compare_anno(
    MS.up.ann, MS.down.ann,
    lab.1 = "Flag.up", lab.2 = "Flag.down",
    type = "BP", stype = "bar"
)
plot_compare_anno(
    MS.up.ann, MS.down.ann,
    lab.1 = "Flag.up", lab.2 = "Flag.down",
    type = "MF", stype = "bar"
)
plot_compare_anno(
    MS.up.ann, MS.down.ann,
    lab.1 = "Flag.up", lab.2 = "Flag.down",
    type = "kegg", stype = "bar"
)
dev.off()

pdf("VS.diff.ann.pdf", 9, 10)
plot_compare_anno(
    VS.up.ann, VS.down.ann,
    lab.1 = "Flag.up", lab.2 = "Flag.down",
    type = "BP", stype = "bar"
)
plot_compare_anno(
    VS.up.ann, VS.down.ann,
    lab.1 = "Flag.up", lab.2 = "Flag.down",
    type = "MF", stype = "bar"
)
plot_compare_anno(
    VS.up.ann, VS.down.ann,
    lab.1 = "Flag.up", lab.2 = "Flag.down",
    type = "kegg", stype = "bar"
)
dev.off()


################### Mesophyll

ms.data <- subset(leaf_Flag.combined.sct, idents = "Mesophyll")

ifnb.list <- SplitObject(ms.data, split.by = "orig.ident")
ifnb.list <- lapply(X = ifnb.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
leaf_Flag.anchors <- FindIntegrationAnchors(
    object.list = ifnb.list, normalization.method = "SCT",
    anchor.features = features
)
ms.sct <- IntegrateData(anchorset = leaf_Flag.anchors, normalization.method = "SCT")
ms.sct <- RunPCA(ms.sct, verbose = FALSE, npcs = 50)
ms.sct <- RunUMAP(ms.sct, reduction = "pca", dims = 1:50)
DimPlot(ms.sct)
ms.sct <- FindNeighbors(ms.sct, dims = 1:50)
ms.sct <- FindClusters(ms.sct, resolution = 0.3)

DefaultAssay(ms.sct) <- "RNA"
ms.sct <- ms.sct %>%
    NormalizeData() %>%
    ScaleData(., features = rownames(.))
ms.sct.markers <- FindAllMarkers(ms.sct, only.pos = T)

ms.sct.markers %>% write.csv("Mesophyll.markers.csv", row.names = F)

pdf("Mesophyll.umap.pdf")
DimPlot(ms.sct, label = T)
DimPlot(ms.sct, label = F)
dev.off()
ms.sct.markers.ann <- lapply(unique(ms.sct.markers$cluster), function(x) {
    return(get_annotation(ms.sct.markers %>% filter(cluster %in% x) %>% pull(gene)))
})
names(ms.sct.markers.ann) <- unique(ms.sct.markers$cluster)
pdf("Mesophyll.cluster.BP.annotation.pdf", 9, 6)
for (i in names(ms.sct.markers.ann)) {
    p <- plot_anno(ms.sct.markers.ann[[i]], type = "BP", style = "bar", title = i)
    print(p)
}
dev.off()

pdf("Mesophyll.cluster.MF.annotation.pdf", 9, 6)
for (i in names(ms.sct.markers.ann)) {
    p <- plot_anno(ms.sct.markers.ann[[i]], type = "MF", style = "bar", title = i)
    print(p)
}
dev.off()

pdf("Mesophyll.cluster.KEGG.annotation.pdf", 9, 6)
for (i in names(ms.sct.markers.ann)) {
    p <- plot_anno(ms.sct.markers.ann[[i]], type = "kegg", style = "bar", title = i)
    print(p)
}
dev.off()
# allbp <- lapply(ms.sct.markers.ann,plot_anno,type = "BP",style = "bar")

###########################################

setwd("Mesophyll/")
counts <- ms.sct@assays$RNA@counts
### run monocle.R
load("Mesophyll.monocle.pd.Rdata")
ms.sct$State <- pd$State
DimPlot(ms.sct, group.by = "State")
Idents(ms.sct) <- ms.sct$State
State.markers <- FindAllMarkers(ms.sct, only.pos = T, min.pct = 0.2)
table(State.markers$cluster)
for (i in unique(State.markers$cluster)) {
    anno <- get_annotation(State.markers %>% filter(cluster %in% i) %>% pull(gene))
    pdf(paste0("Mesophyll.", "State.", i, ".annotation.pdf"))
    for (xx in names(anno)) {
        p <- plot_anno(anno, type = xx, style = "bar", title = i)
        print(p)
    }
    dev.off()
}
save(State.markers, file = "State.markers.Rdata")

DefaultAssay(leaf_Flag.combined.sct) <- "integrated"
leaf_Flag.combined.sct$tissue_cluster <- paste(leaf_Flag.combined.sct$tissue, leaf_Flag.combined.sct$cluster, sep = "-")
texp <- AverageExpression(leaf_Flag.combined.sct, assays = "RNA", slot = "data", group.by = "tissue_cluster")

p <- plot_triangle(
    cor(texp$RNA[VariableFeatures(leaf_Flag.combined.sct), -9], method = "spearman"),
    limits = c(0, 0.9)
) +
    scale_fill_gradient(low = "yellow", high = "red")
ggsave(p, filename = "leaf_flag.tissue.correlation.spearman.pdf", width = 12, height = 7)



clu <- unique(leaf_Flag.combined.sct$cluster)[1:5]

allconserved.markers <- parallel::mclapply(clu, function(x) {
    marker1 <- FindConservedMarkers(leaf_Flag.combined.sct, ident.1 = x, grouping.var = "tissue", verbose = FALSE)
    return(data.frame(marker1, cluster = x, gene = rownames(markers1)))
}, mc.cores = 10) %>% do.call(rbind, .)

tissue_order <- c("Root", "ST", "leaf", "Flag", "SAM", "Bud", "SP", "Seed")
tissue.colors <- scales::hue_pal()(8)
names(tissue.colors) <- c("Bud", "Flag", "leaf", "Root", "SAM", "Seed", "SP", "ST")
DimPlot(leaf_Flag.combined.sct, group.by = "tissue") + scale_color_manual(values = tissue.colors[unique(leaf_Flag.combined.sct$tissue)])
