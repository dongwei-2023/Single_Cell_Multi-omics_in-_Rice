library(Seurat)
library(WGCNA)
load("/alldata.Rdata")

library(hdWGCNA)
library(igraph)
library(cowplot)
library(patchwork)
library(WGCNA)
enableWGCNAThreads(nThreads = 28)

DefaultAssay(alldata) <- "RNA"
seurat_obj <- SetupForWGCNA(
    alldata,
    gene_select = "fraction", # the gene selection approach
    fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
    wgcna_name = "alldata" # the name of the hdWGCNA experiment
)
seurat_obj$group <- gsub("_.+", "", seurat_obj$orig.ident)
seurat_obj <- MetacellsByGroups(
    seurat_obj = seurat_obj,
    group.by = c("tissue"), # specify the columns in seurat_obj@meta.data to group by
    k = 25, # nearest-neighbors parameter
    max_shared = 100, # maximum number of shared cells between two metacells
    ident.group = "tissue" # set the Idents of the metacell seurat object
)
seurat_obj <- NormalizeMetacells(seurat_obj)

DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- SetDatExpr(
    seurat_obj,
    group_name = c(seurat_obj$tissue %>% unique()), # the name of the group of interest in the group.by column
    group.by = "tissue", # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
    assay = "RNA", # using RNA assay
    slot = "data" # using normalized data
)
seurat_obj <- TestSoftPowers(
    seurat_obj,
    networkType = "signed" # you can also use "unsigned" or "signed hybrid"
)
plot_list <- PlotSoftPowers(seurat_obj)
wrap_plots(plot_list, ncol = 2)

power_table <- GetPowerTable(seurat_obj)
seurat_obj <- ConstructNetwork(
    seurat_obj,
    soft_power = 8,
    setDatExpr = FALSE,
    tom_name = "alldata" # name of the topoligical overlap matrix written to disk
)
PlotDendrogram(seurat_obj, main = "hdWGCNA Dendrogram")
TOM <- GetTOM(seurat_obj)

seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))
seurat_obj <- ModuleEigengenes(seurat_obj, group.by.vars = "group")
# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
    seurat_obj,
    group.by = "tissue", group_name = c(seurat_obj$tissue %>% unique()), sparse = F
)
PlotKMEs(seurat_obj, ncol = 5)
seurat_obj <- ResetModuleNames(
    seurat_obj,
    new_name = "alldata-M"
)
select <- dplyr::select

plot_list <- ModuleFeaturePlot(seurat_obj, features = "hMEs", order = TRUE, reduction = "wsnnumap", raster = T)
pdf("allModules.new.pdf")
lapply(plot_list, print)
dev.off()

ModuleFeaturePlot(
    module_names = "alldata-M2",
    seurat_obj,
    features = "hMEs", # plot the hMEs
    order = TRUE # order so the points with highest hMEs are on top
) + scale_color_gradient2(low = "grey75", mid = "grey95", high = "red") + NoLegend()

modules <- GetModules(seurat_obj)
modules %>%
    filter(module %in% "alldata-M1") %>%
    pull(gene_name) %>%
    head()

source("/enrichment.R")

all.module.anno <- lapply(modules$module %>% unique(), function(x) {
    anno <- get_annotation(modules %>% filter(module %in% x) %>% pull(gene_name))
    return(anno)
})
names(all.module.anno) <- modules$module %>% unique()

pdf("alldata.module.annotation.pdf", 9, 8)
lapply(names(all.module.anno), function(x) {
    # plot_anno(all.module.anno[[x]],type = c("BP","kegg"),style = "bar",title = x)
    do.call(rbind, all.module.anno[[x]]) %>% write.csv(paste0(x, ".anno.csv"), row.names = F)
})
dev.off()

plot_anno(M1.anno, type = c("BP", "kegg"), style = "bar", title = "alldata-M1")
M2.anno <- get_annotation(modules %>% filter(module %in% "alldata-M2") %>% pull(gene_name))
plot_anno(M2.anno, type = "BP", style = "bar")
M3.anno <- get_annotation(modules %>% filter(module %in% "alldata-M3") %>% pull(gene_name))
plot_anno(M3.anno, type = "BP", style = "bar")
M4.anno <- get_annotation(modules %>% filter(module %in% "alldata-M4") %>% pull(gene_name))
plot_anno(M4.anno, type = "BP", style = "bar")
M5.anno <- get_annotation(modules %>% filter(module %in% "alldata-M5") %>% pull(gene_name))
plot_anno(M5.anno, type = "BP", style = "bar")
M6.anno <- get_annotation(modules %>% filter(module %in% "alldata-M6") %>% pull(gene_name))
plot_anno(M6.anno, type = "BP", style = "bar")


M4.genes <- module.genes %>%
    filter(module %in% "alldata-M4") %>%
    top_n(500, `kME_alldata-M4`)
mm.genes <- rbind(M4.genes, module.genes %>% filter(!(module %in% "alldata-M4"))) %>%
    filter(!(color %in% "grey")) %>%
    group_by(color) %>%
    dplyr::sample_frac(0.2)
mm.cor <- cor(t(as.matrix(seurat_obj@assays$RNA@data[mm.genes$gene_name, sample(1:ncol(seurat_obj), 3000)])))
pheatmap::pheatmap(mm.cor, scale = "none", show_rownames = F, show_colnames = F, annotation_col = data.frame(module = mm.genes$module, row.names = mm.genes$gene_name), breaks = seq(-1, 1, length.out = 100), color = colorRampPalette(scales::div_gradient_pal("blue", "grey", "red")(seq(0, 1, length.out = 25)))(100), cluster_cols = F, cluster_rows = F, filename = "module.correlation.heatmap.pdf", width = 9, height = 7)
