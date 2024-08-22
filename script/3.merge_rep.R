########## set up environment
library(Seurat)
library(dplyr)
library(Seurat)
library(Signac)
library(ggplot2)
library(cowplot)


UMI <- 1000
samples <- c("leaf_1", "leaf_2", "Root_1", "Root_2", "Bud_2", "Flag_1", "Flag_2", "Seed_31", "Seed_32", "SAM_1", "SAM_2", "Bud_21", "Bud_22")

all.list <- lapply(samples, function(sampleName) {
    load(paste0(sampleName, "/data.", "1.9", ".Rdata"))
    p <- DimPlot(pbmc.rna, label = T)
    ggsave(p, filename = paste0(sampleName, ".UMAP.pdf"), width = 8, height = 7)
    pbmc.markers <- FindAllMarkers(pbmc.rna, min.pct = 0.5, min.diff.pct = 0.2, only.pos = T)
    pbmc.markers <- data.frame(pbmc.markers, g1 = rownames(pbmc.markers)) %>% left_join(genelist, by = "g1")
    write.csv(pbmc.markers, file = "markers.csv")
    pbmc.markers %>%
        group_by(cluster) %>%
        top_n(n = 3, wt = avg_log2FC) -> top5
    p <- DotPlot(pbmc.rna, features = unique(top5$gene), cols = c("grey", "red")) + NoLegend() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(p, filename = "top3.DotPlot.pdf", width = 8, height = 5)
    pdf("Feature.UMAP.pdf")
    for (gene in top5$gene) {
        p <- FeaturePlot(pbmc.rna, features = gene, cols = c("grey", "red"), max.cutoff = "q90") + ggtitle(paste(gene, gName(gene, invert = T), sep = "/"))
        print(p)
    }
    return(pbmc.rna)
})
names(all.list) <- samples
save(all.list, file = "all.list.Rdata")

Flag <- merge(all.list$Flag_1, all.list$Flag_2)
leaf <- merge(all.list$leaf_1, all.list$leaf_2)
Root <- merge(all.list$Root_1, all.list$Root_2)
SAM <- merge(all.list$SAM_1, all.list$SAM_2)
SP <- merge(all.list$SP_1, all.list$SP_2)
ST <- merge(all.list$ST_1, all.list$ST_2)
Seed <- merge(all.list$Seed_31, all.list$Seed_32)
Bud <- merge(all.list$Bud_21, all.list$Bud_22)
tissues <- c("Bud", "Flag", "leaf", "Root", "SAM", "Seed", "SP", "ST")
alldata <- lapply(tissues, function(tis) {
    obj <- get(tis)
    obj <- ScaleData(obj, features = rownames(obj))
    obj <- FindVariableFeatures(obj)
    obj <- RunPCA(obj, npcs = 50)
    obj <- RunUMAP(obj, dims = 1:30)
    return(obj)
})
names(alldata) <- tissues
save(alldata, file = "alldata.Rdata")
