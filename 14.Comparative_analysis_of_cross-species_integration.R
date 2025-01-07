library(dplyr)
library(Seurat)
load("/GSE225118_Maize_Sorghum_Setaria_Cells_Nuclei_expression_matrix_raw_counts.RData")
Orthogroups <- read.delim("/Orthogroups.tsv")
colnames(Orthogroups) <- gsub("\\..+", "", colnames(Orthogroups))
Orthogroups <- Orthogroups %>% filter(!(Oryza_sativa %in% ""), !(Setaria_viridis %in% ""), !(Tair10 %in% ""), !(Zea_mays %in% ""))

allorthologs <- lapply(1:nrow(Orthogroups), function(x) {
    a1 <- unique(gsub("\\..+", "", unlist(strsplit(Orthogroups$Tair10[x], ", "))))
    a2 <- unique(gsub("_.+", "", unlist(strsplit(Orthogroups$Zea_mays[x], ", "))))
    a3 <- unique(gsub("t", "g", gsub("-.+", "", unlist(strsplit(Orthogroups$Oryza_sativa[x], ", ")))))
    a4 <- unique(gsub("-transcript.+", "", unlist(strsplit(Orthogroups$Sorghum_bicolor[x], ", "))))
    a5 <- unique(gsub("-transcript.+", "", unlist(strsplit(Orthogroups$Setaria_viridis[x], ", "))))
    return(list(tair = a1, zea_mays = a2, ory_sat = a3, sor_bic = a4, set_vir = a5))
})
names(allorthologs) <- Orthogroups$Orthogroup


maize.data <- Maize_Sorghum_Setaria_Nuclei_Cells_SCT_OrthologousgenesGoodname_Sept2022@assays$RNA@counts[, Maize_Sorghum_Setaria_Nuclei_Cells_SCT_OrthologousgenesGoodname_Sept2022$Species %in% "Maize"] %>% as.matrix()
Setaria.data <- Maize_Sorghum_Setaria_Nuclei_Cells_SCT_OrthologousgenesGoodname_Sept2022@assays$RNA@counts[, Maize_Sorghum_Setaria_Nuclei_Cells_SCT_OrthologousgenesGoodname_Sept2022$Species %in% "Setaria"] %>% as.matrix()
Sorghum.data <- Maize_Sorghum_Setaria_Nuclei_Cells_SCT_OrthologousgenesGoodname_Sept2022@assays$RNA@counts[, Maize_Sorghum_Setaria_Nuclei_Cells_SCT_OrthologousgenesGoodname_Sept2022$Species %in% "Sorghum"] %>% as.matrix()


maize.counts <- lapply(names(allorthologs), function(x) {
    ids <- intersect(rownames(maize.data), allorthologs[[x]]$zea_mays)
    if (length(ids) > 1) {
        return(apply(maize.data[ids, ], 2, sum))
    }
    if (length(ids) == 1) {
        return(maize.data[ids, ])
    }
    if (length(ids) == 0) {
        return(rep(0, ncol(maize.data)))
    }
})
maize.data <- do.call(rbind, maize.counts)
rownames(maize.data) <- names(allorthologs)

Setaria.counts <- lapply(names(allorthologs), function(x) {
    ids <- intersect(rownames(Setaria.data), allorthologs[[x]]$set_vir)
    if (length(ids) > 1) {
        return(apply(Setaria.data[ids, ], 2, sum))
    }
    if (length(ids) == 1) {
        return(Setaria.data[ids, ])
    }
    if (length(ids) == 0) {
        return(rep(0, ncol(Setaria.data)))
    }
})
Setaria.data <- do.call(rbind, Setaria.counts)
rownames(Setaria.data) <- names(allorthologs)

Sorghum.counts <- lapply(names(allorthologs), function(x) {
    ids <- intersect(rownames(Sorghum.data), allorthologs[[x]]$sor_bic)
    if (length(ids) > 1) {
        return(apply(Sorghum.data[ids, ], 2, sum))
    }
    if (length(ids) == 1) {
        return(Sorghum.data[ids, ])
    }
    if (length(ids) == 0) {
        return(rep(0, ncol(Sorghum.data)))
    }
})
Sorghum.data <- do.call(rbind, Sorghum.counts)
rownames(Sorghum.data) <- names(allorthologs)

load("/alldata.Rdata")
Rice.data <- as.matrix(alldata$Root@assays$RNA@counts)
Rice.counts <- lapply(names(allorthologs), function(x) {
    ids <- intersect(rownames(Rice.data), allorthologs[[x]]$ory_sat)
    if (length(ids) > 1) {
        return(apply(Rice.data[ids, ], 2, sum))
    }
    if (length(ids) == 1) {
        return(Rice.data[ids, ])
    }
    if (length(ids) == 0) {
        return(rep(0, ncol(Rice.data)))
    }
})
Rice.data <- do.call(rbind, Rice.counts)
rownames(Rice.data) <- names(allorthologs)

maize <- CreateSeuratObject(maize.data, project = "maize")
Rice <- CreateSeuratObject(Rice.data, project = "Rice")
Setaria <- CreateSeuratObject(Setaria.data, project = "Setaria")
Sorghum <- CreateSeuratObject(Sorghum.data, project = "Sorghum")
Tair <- CreateSeuratObject(Tair.data, project = "Tair")

maize.data <- Maize_Sorghum_Setaria_Nuclei_Cells_SCT_OrthologousgenesGoodname_Sept2022@assays$RNA@counts[, Maize_Sorghum_Setaria_Nuclei_Cells_SCT_OrthologousgenesGoodname_Sept2022$Species %in% "Maize"] %>% as.matrix()
maize.counts <- lapply(names(allorthologs), function(x) {
    ids <- intersect(rownames(maize.data), allorthologs[[x]]$zea_mays)
    if (length(ids) > 1) {
        return(apply(maize.data[ids, ], 2, sum))
    }
    if (length(ids) == 1) {
        return(maize.data[ids, ])
    }
    if (length(ids) == 0) {
        return(rep(0, ncol(maize.data)))
    }
})
maize.data <- do.call(rbind, maize.counts)
rownames(maize.data) <- names(allorthologs)
maize <- CreateSeuratObject(maize.data, project = "maize")


adata <- merge(Rice, c(maize, Setaria, Sorghum, Tair))
ifnb.list <- SplitObject(adata, split.by = "orig.ident")
ifnb.list <- lapply(X = ifnb.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 2000)
ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "SCT", anchor.features = features)
immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")
immune.combined.sct <- RunPCA(immune.combined.sct, verbose = FALSE)
immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:30, a = 1.2, b = 0.9, n.neighbors = 50)
DimPlot(immune.combined.sct)

DimPlot(immune.combined.sct, group.by = "cluster", split.by = "orig.ident", label = T, ncol = 3, repel = T) + NoLegend()

save(immune.combined.sct, file = "5species.merge.Rdata")


prop.table(table(subdata$orig.ident, subdata$cluster3), margin = 1) %>%
    data.frame() %>%
    ggplot(aes(Var2, Freq, fill = Var1)) +
    geom_histogram(stat = "identity", position = position_dodge2()) +
    scale_fill_npg() +
    theme_bw() +
    labs(x = "") +
    BoldTheme()


cn.list <- read.delim("/C-N.list")
cn.list$term %>% unique()
# [1] "Nitrogen_absorption"      "Amino_acid_synthesis"     "Photosynthesis"
# [4] "Sugar_metabolism-related" "Transport-related"        "TCA_cycle"
for (i in cn.list$term %>% unique()) {
    ids <- lapply(cn.list %>% filter(term %in% i) %>% pull(gene), function(x) {
        id <- ory.orth %>%
            filter(grepl(x, ory)) %>%
            pull(id)
        if (length(id) > 0) {
            return(data.frame(os = x, id = id))
        }
    })
    ids <- ids %>% do.call(rbind, .)
    tmpdata <- subdata.avg$RNA[ids$id, ]
    rownames(tmpdata) <- ids$os
    pheatmap(tmpdata,
        scale = "row", cluster_cols = F,
        breaks = seq(-2, 2, length.out = 100), filename = paste0(i, ".CN.pdf"), height = 3 + length(ids) / 10
    )
}

Endo.markers <- FindConservedMarkers(subdata, ident.1 = "Endodermis", grouping.var = "orig.ident", verbose = FALSE)
Epi.markers <- FindConservedMarkers(subdata, ident.1 = "Epidermis", grouping.var = "orig.ident", verbose = FALSE)
Cortex.markers <- FindConservedMarkers(subdata, ident.1 = "Cortex", grouping.var = "orig.ident", verbose = FALSE)
Vas.markers <- FindConservedMarkers(subdata, ident.1 = "Vascularcylinder", grouping.var = "orig.ident", verbose = FALSE)

Vas.ids <- Vas.markers %>%
    filter_at(dplyr::vars(ends_with("log2FC")), dplyr::all_vars(. > 0)) %>%
    rownames()
Endo.ids <- Endo.markers %>%
    filter_at(dplyr::vars(ends_with("log2FC")), dplyr::all_vars(. > 0)) %>%
    rownames()
Epi.ids <- Epi.markers %>%
    filter_at(dplyr::vars(ends_with("log2FC")), dplyr::all_vars(. > 0)) %>%
    rownames()
Cortex.ids <- Cortex.markers %>%
    filter_at(dplyr::vars(ends_with("log2FC")), dplyr::all_vars(. > 0)) %>%
    rownames()

rand_cells <- data.frame(cell = colnames(subdata), cluster = subdata$cluster3, sample = subdata$orig.ident) %>%
    group_by(sample, cluster) %>%
    sample_n(300)

pdf("Conserved.heatmap.pdf", 12, 8)
DoHeatmap(subset(subdata, cells = rand_cells %>% filter(sample %in% "Rice") %>% pull(cell)), features = c(Endo.ids, Epi.ids, Vas.ids, Cortex.ids)) + scale_fill_gradient2(low = "blue", mid = "white", high = "#D73027")
DoHeatmap(subset(subdata, cells = rand_cells %>% filter(sample %in% "Maize") %>% pull(cell)), features = c(Endo.ids, Epi.ids, Vas.ids, Cortex.ids)) + scale_fill_gradient2(low = "blue", mid = "white", high = "#D73027")
DoHeatmap(subset(subdata, cells = rand_cells %>% filter(sample %in% "Set") %>% pull(cell)), features = c(Endo.ids, Epi.ids, Vas.ids, Cortex.ids)) + scale_fill_gradient2(low = "blue", mid = "white", high = "#D73027")
DoHeatmap(subset(subdata, cells = rand_cells %>% filter(sample %in% "Sorghum") %>% pull(cell)), features = c(Endo.ids, Epi.ids, Vas.ids, Cortex.ids)) + scale_fill_gradient2(low = "blue", mid = "white", high = "#D73027")
DoHeatmap(subset(subdata, cells = rand_cells %>% filter(sample %in% "Tair") %>% pull(cell)), features = c(Endo.ids, Epi.ids, Vas.ids, Cortex.ids)) + scale_fill_gradient2(low = "blue", mid = "white", high = "#D73027")
dev.off()
for (i in unique(subdata$orig.ident)) {
    hdata <- subset(subdata, cells = rand_cells %>% filter(sample %in% i) %>% pull(cell))
    pheatmap(hdata@assays$RNA@data[c(Cortex.ids, Endo.ids, Epi.ids, Vas.ids), order(hdata$orig.ident)],
        scale = "row", show_colnames = F,
        breaks = seq(-2, 2, length.out = 100), cluster_cols = F, cluster_rows = F,
        annotation_col = data.frame(cluster = hdata$cluster3, row.names = colnames(hdata)),
        # color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="RdBu")))(100),
        color = viridis::viridis_pal()(100),
        gaps_col = cumsum(table(hdata$cluster3)),
        gaps_row = cumsum(sapply(list(Cortex.ids, Endo.ids, Epi.ids, Vas.ids), length)),
        filename = paste0(i, ".conserved.pdf"),
        width = 12, height = 8
    )
}

Endo.xids <- setdiff(Endo.markers %>% filter(Rice_avg_log2FC > 0) %>% rownames(), Endo.ids)
Epi.xids <- setdiff(Epi.markers %>% filter(Rice_avg_log2FC > 0) %>% rownames(), Epi.ids)
Cortex.xids <- setdiff(Cortex.markers %>% filter(Rice_avg_log2FC > 0) %>% rownames(), Cortex.ids)
Vas.xids <- setdiff(Vas.markers %>% filter(Rice_avg_log2FC > 0) %>% rownames(), Vas.ids)

for (i in unique(subdata$orig.ident)) {
    hdata <- subset(subdata, cells = rand_cells %>% filter(sample %in% i) %>% pull(cell))
    pheatmap(hdata@assays$RNA@data[c(Cortex.xids, Endo.xids, Epi.xids, Vas.xids), order(hdata$orig.ident)],
        scale = "row", show_colnames = F,
        breaks = seq(-2, 2, length.out = 100), cluster_cols = F, cluster_rows = F,
        annotation_col = data.frame(cluster = hdata$cluster3, row.names = colnames(hdata)),
        # color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="RdBu")))(100),
        color = viridis::viridis_pal()(100),
        gaps_col = cumsum(table(hdata$cluster3)),
        gaps_row = cumsum(sapply(list(Cortex.xids, Endo.xids, Epi.xids, Vas.xids), length)),
        filename = paste0(i, ".unconserved.pdf"),
        width = 12, height = 8
    )
}


rice_maize.diff <- parallel::mclapply(subdata$cluster3 %>% unique(), function(x) {
    rice.markers <- FindMarkers(subdata, ident.1 = "Rice", ident.2 = "Maize", subset.ident = x, group.by = "orig.ident", logfc.threshold = 0.01, test.use = "negbinom")
    return(data.frame(rice.markers, cluster = x, gene = rownames(rice.markers)))
}, mc.cores = 5)
p <- rice_maize.diff %>%
    do.call(rbind, .) %>%
    mutate(sig = ifelse(p_val > 0.01, "nosig", ifelse(avg_log2FC > 0, "up", "down"))) %>%
    ggplot(aes(avg_log2FC, -log10(p_val), color = cluster)) +
    geom_point(size = 0.3) +
    coord_cartesian(xlim = c(-4, 4)) +
    theme_bw() +
    scale_color_manual(values = clusters_colors[unique(subdata$cluster3)]) +
    facet_wrap(~cluster, ncol = 2) +
    labs(x = "Rice vs Maize log2FC")
ggsave(p, filename = "Rice_vs_maize.log2FC.pdf", width = 7, height = 5.5)

xx <- get_annotation(allorthologs[Vas.ids] %>% unlist() %>% grep("Os", ., value = T) %>% as.character())
save_anno(xx, filename = "Vascular.conserved.GO_KEGG.csv")
xx <- get_annotation(allorthologs[Epi.ids] %>% unlist() %>% grep("Os", ., value = T) %>% as.character())
save_anno(xx, filename = "Epidermis.conserved.GO_KEGG.csv")
xx <- get_annotation(allorthologs[Endo.ids] %>% unlist() %>% grep("Os", ., value = T) %>% as.character())
save_anno(xx, filename = "Endodermis.conserved.GO_KEGG.csv")
xx <- get_annotation(allorthologs[Cortex.ids] %>% unlist() %>% grep("Os", ., value = T) %>% as.character())
save_anno(xx, filename = "Cortex.conserved.GO_KEGG.csv")
for (i in c("Vas", "Epi", "Endo", "Cortex")) {
    aa <- get(paste0(i, ".ids"))
    xx <- get_annotation(allorthologs[aa] %>% unlist() %>% grep("Os", ., value = T) %>% as.character())
    pdf(paste0(i, ".conserved.anno.pdf"))
    print(plot_anno(xx))
    print(plot_anno(xx, type = "MF"))
    print(plot_anno(xx, type = "kegg"))
    dev.off()
}

gg_color_hue <- function(n) {
    hues <- seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}
ccolor <- gg_color_hue(24)
names(ccolor) <- unique(immune.combined.sct$cluster3)
ccolor[intersect(names(ccolor), names(clusters_colors))] <- clusters_colors[intersect(names(ccolor), names(clusters_colors))]
DimPlot(immune.combined.sct, label = T, repel = T) + scale_color_manual(values = ccolor)

rice_maize.vas.upID <- rice_maize.diff[[3]] %>%
    filter(avg_log2FC > 0) %>%
    rownames()
rice_maize.vas.downID <- rice_maize.diff[[3]] %>%
    filter(avg_log2FC < 0) %>%
    rownames()
xx1 <- get_annotation(allorthologs[rice_maize.vas.upID] %>% unlist() %>% grep("Os", ., value = T) %>% as.character())
xx2 <- get_annotation(allorthologs[rice_maize.vas.downID] %>% unlist() %>% grep("Os", ., value = T) %>% as.character())
plot_compare_anno(xx1, xx2, lab.1 = "Rice_up", lab.2 = "Maize_up")
save_anno(xx1, filename = "Rice_Maize.vascular.up.anno.csv")
save_anno(xx2, filename = "Rice_Maize.vascular.down.anno.csv")

rice_tair.diff <- parallel::mclapply(subdata$cluster3 %>% unique(), function(x) {
    rice.markers <- FindMarkers(subdata, ident.1 = "Rice", ident.2 = "Tair", subset.ident = x, group.by = "orig.ident", logfc.threshold = 0.01, test.use = "negbinom")
    return(data.frame(rice.markers, cluster = x, gene = rownames(rice.markers)))
}, mc.cores = 5)
p <- rice_tair.diff %>%
    do.call(rbind, .) %>%
    mutate(sig = ifelse(p_val > 0.01, "nosig", ifelse(avg_log2FC > 0, "up", "down"))) %>%
    ggplot(aes(avg_log2FC, -log10(p_val), color = cluster)) +
    geom_point(size = 0.3) +
    coord_cartesian(xlim = c(-4, 4)) +
    theme_bw() +
    scale_color_manual(values = clusters_colors[unique(subdata$cluster3)]) +
    facet_wrap(~cluster, ncol = 2) +
    labs(x = "Rice vs tair log2FC")
ggsave(p, filename = "Rice_vs_tair.log2FC.pdf", width = 7, height = 5.5)
rice_tair.vas.upID <- rice_tair.diff[[3]] %>%
    filter(avg_log2FC > 0) %>%
    rownames()
rice_tair.vas.downID <- rice_tair.diff[[3]] %>%
    filter(avg_log2FC < 0) %>%
    rownames()
xx1 <- get_annotation(allorthologs[rice_tair.vas.upID] %>% unlist() %>% grep("Os", ., value = T) %>% as.character())
xx2 <- get_annotation(allorthologs[rice_tair.vas.downID] %>% unlist() %>% grep("Os", ., value = T) %>% as.character())
plot_compare_anno(xx1, xx2, lab.1 = "Rice_up", lab.2 = "tair_up")
save_anno(xx1, filename = "Rice_tair.vascular.up.anno.csv")
save_anno(xx2, filename = "Rice_tair.vascular.down.anno.csv")

rice_Set.diff <- parallel::mclapply(subdata$cluster3 %>% unique(), function(x) {
    rice.markers <- FindMarkers(subdata, ident.1 = "Rice", ident.2 = "Set", subset.ident = x, group.by = "orig.ident", logfc.threshold = 0.01, test.use = "negbinom")
    return(data.frame(rice.markers, cluster = x, gene = rownames(rice.markers)))
}, mc.cores = 5)
p <- rice_Set.diff %>%
    do.call(rbind, .) %>%
    mutate(sig = ifelse(p_val > 0.01, "nosig", ifelse(avg_log2FC > 0, "up", "down"))) %>%
    ggplot(aes(avg_log2FC, -log10(p_val), color = cluster)) +
    geom_point(size = 0.3) +
    coord_cartesian(xlim = c(-4, 4)) +
    theme_bw() +
    scale_color_manual(values = clusters_colors[unique(subdata$cluster3)]) +
    facet_wrap(~cluster, ncol = 2) +
    labs(x = "Rice vs Set log2FC")
ggsave(p, filename = "Rice_vs_Set.log2FC.pdf", width = 7, height = 5.5)
rice_Set.vas.upID <- rice_Set.diff[[3]] %>%
    filter(avg_log2FC > 0) %>%
    rownames()
rice_Set.vas.downID <- rice_Set.diff[[3]] %>%
    filter(avg_log2FC < 0) %>%
    rownames()
xx1 <- get_annotation(allorthologs[rice_Set.vas.upID] %>% unlist() %>% grep("Os", ., value = T) %>% as.character())
xx2 <- get_annotation(allorthologs[rice_Set.vas.downID] %>% unlist() %>% grep("Os", ., value = T) %>% as.character())
plot_compare_anno(xx1, xx2, lab.1 = "Rice_up", lab.2 = "Set_up")
save_anno(xx1, filename = "Rice_Set.vascular.up.anno.csv")
save_anno(xx2, filename = "Rice_Set.vascular.down.anno.csv")

rice_Sorghum.diff <- parallel::mclapply(subdata$cluster3 %>% unique(), function(x) {
    rice.markers <- FindMarkers(subdata, ident.1 = "Rice", ident.2 = "Sorghum", subset.ident = x, group.by = "orig.ident", logfc.threshold = 0.01, test.use = "negbinom")
    return(data.frame(rice.markers, cluster = x, gene = rownames(rice.markers)))
}, mc.cores = 5)
p <- rice_Sorghum.diff %>%
    do.call(rbind, .) %>%
    mutate(sig = ifelse(p_val > 0.01, "nosig", ifelse(avg_log2FC > 0, "up", "down"))) %>%
    ggplot(aes(avg_log2FC, -log10(p_val), color = cluster)) +
    geom_point(size = 0.3) +
    coord_cartesian(xlim = c(-4, 4)) +
    theme_bw() +
    scale_color_manual(values = clusters_colors[unique(subdata$cluster3)]) +
    facet_wrap(~cluster, ncol = 2) +
    labs(x = "Rice vs Sorghum log2FC")
ggsave(p, filename = "Rice_vs_Sorghum.log2FC.pdf", width = 7, height = 5.5)
rice_Sorghum.vas.upID <- rice_Sorghum.diff[[3]] %>%
    filter(avg_log2FC > 0) %>%
    rownames()
rice_Sorghum.vas.downID <- rice_Sorghum.diff[[3]] %>%
    filter(avg_log2FC < 0) %>%
    rownames()
xx1 <- get_annotation(allorthologs[rice_Sorghum.vas.upID] %>% unlist() %>% grep("Os", ., value = T) %>% as.character())
xx2 <- get_annotation(allorthologs[rice_Sorghum.vas.downID] %>% unlist() %>% grep("Os", ., value = T) %>% as.character())
plot_compare_anno(xx1, xx2, lab.1 = "Rice_up", lab.2 = "Sorghum_up")
save_anno(xx1, filename = "Rice_Sorghum.vascular.up.anno.csv")
save_anno(xx2, filename = "Rice_Sorghum.vascular.down.anno.csv")




Vas.ids1 <- Vas.markers %>%
    dplyr::select(-Tair_avg_log2FC) %>%
    filter_at(dplyr::vars(ends_with("log2FC")), dplyr::all_vars(. > 0)) %>%
    rownames()
Endo.ids1 <- Endo.markers %>%
    dplyr::select(-Tair_avg_log2FC) %>%
    filter_at(dplyr::vars(ends_with("log2FC")), dplyr::all_vars(. > 0)) %>%
    rownames()
Epi.ids1 <- Epi.markers %>%
    dplyr::select(-Tair_avg_log2FC) %>%
    filter_at(dplyr::vars(ends_with("log2FC")), dplyr::all_vars(. > 0)) %>%
    rownames()
Cortex.ids1 <- Cortex.markers %>%
    dplyr::select(-Tair_avg_log2FC) %>%
    filter_at(dplyr::vars(ends_with("log2FC")), dplyr::all_vars(. > 0)) %>%
    rownames()

rand_cells <- data.frame(cell = colnames(subdata), cluster = subdata$cluster3, sample = subdata$orig.ident) %>%
    group_by(sample, cluster) %>%
    sample_n(300)

pdf("Conserved.heatmap.pdf", 12, 8)
DoHeatmap(subset(subdata, cells = rand_cells %>% filter(sample %in% "Rice") %>% pull(cell)), features = c(Endo.ids, Epi.ids, Vas.ids, Cortex.ids)) + scale_fill_gradient2(low = "blue", mid = "white", high = "#D73027")
DoHeatmap(subset(subdata, cells = rand_cells %>% filter(sample %in% "Maize") %>% pull(cell)), features = c(Endo.ids, Epi.ids, Vas.ids, Cortex.ids)) + scale_fill_gradient2(low = "blue", mid = "white", high = "#D73027")
DoHeatmap(subset(subdata, cells = rand_cells %>% filter(sample %in% "Set") %>% pull(cell)), features = c(Endo.ids, Epi.ids, Vas.ids, Cortex.ids)) + scale_fill_gradient2(low = "blue", mid = "white", high = "#D73027")
DoHeatmap(subset(subdata, cells = rand_cells %>% filter(sample %in% "Sorghum") %>% pull(cell)), features = c(Endo.ids, Epi.ids, Vas.ids, Cortex.ids)) + scale_fill_gradient2(low = "blue", mid = "white", high = "#D73027")
DoHeatmap(subset(subdata, cells = rand_cells %>% filter(sample %in% "Tair") %>% pull(cell)), features = c(Endo.ids, Epi.ids, Vas.ids, Cortex.ids)) + scale_fill_gradient2(low = "blue", mid = "white", high = "#D73027")
dev.off()
for (i in unique(subdata$orig.ident)[1:4]) {
    hdata <- subset(subdata, cells = rand_cells %>% filter(sample %in% i) %>% pull(cell))
    pheatmap(hdata@assays$RNA@data[c(Cortex.ids1, Endo.ids1, Epi.ids1, Vas.ids1), order(hdata$orig.ident)],
        scale = "row", show_colnames = F,
        breaks = seq(-2, 2, length.out = 100), cluster_cols = F, cluster_rows = F,
        annotation_col = data.frame(cluster = hdata$cluster3, row.names = colnames(hdata)),
        # color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="RdBu")))(100),
        color = viridis::viridis_pal()(100),
        gaps_col = cumsum(table(hdata$cluster3)),
        gaps_row = cumsum(sapply(list(Cortex.ids1, Endo.ids1, Epi.ids1, Vas.ids1), length)),
        filename = paste0(i, ".conserved.noTair.pdf"),
        width = 12, height = 8
    )
}

lapply(grep("diff$", ls(), value = T), function(x) {
    get(x) %>%
        do.call(rbind, .) %>%
        write.csv(paste0(x, ".csv"), row.names = F)
})

allorthologs.infos <- lapply(names(allorthologs), function(x) {
    return(data.frame(orthoID = x, lapply(allorthologs[[x]], paste, collapse = ";")))
}) %>% do.call(rbind, .)
allorthologs.infos %>% write.csv("allorthologs.infos.csv", row.names = F)
