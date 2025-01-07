########## set up environment
library(Seurat)
library(dplyr)
library(Seurat)
library(Signac)
library(ggplot2)
library(cowplot)

get_doublet <- function(obj) {
    require(DoubletFinder)
    dcell.rate <- read.delim("double_rate.txt")
    sweep.res.list_kidney <- paramSweep_v3(obj, PCs = 1:30, sct = FALSE)
    sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
    bcmvn_kidney <- find.pK(sweep.stats_kidney)
    pk <- bcmvn_kidney %>%
        filter(BCmetric == max(BCmetric)) %>%
        pull(pK) %>%
        as.character() %>%
        as.numeric()
    homotypic.prop <- modelHomotypic(Idents(obj))
    drate <- dcell.rate %>%
        filter(cells > ncol(obj)) %>%
        dplyr::slice(1) %>%
        pull(rate)
    drate <- ifelse(isEmpty(drate), 7.6, drate)
    nExp_poi <- round(drate * nrow(obj@meta.data) / 100)
    nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop / 100))
    pANN_name <- paste("pANN", "0.25", pk, nExp_poi, sep = "_")
    pANN_name2 <- paste("DF.classifications", "0.25", pk, nExp_poi.adj, sep = "_")
    seu_kidney <- doubletFinder_v3(obj, PCs = 1:30, pN = 0.25, pK = pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
    seu_kidney <- doubletFinder_v3(seu_kidney, PCs = 1:30, pN = 0.25, pK = pk, nExp = nExp_poi.adj, reuse.pANN = pANN_name, sct = FALSE)
    obj$doublet <- seu_kidney[[pANN_name2]]
    return(obj)
}

############### create seurat obj
samples <- c("leaf_1", "leaf_2", "Root_1", "Root_2", "Bud_2", "Flag_1", "Flag_2", "Seed_31", "Seed_32", "SAM_1", "SAM_2", "Bud_21", "Bud_22")
for (sampleName in samples) {
    a <- Read10X(paste0("../", sampleName, "/outs/filtered_feature_bc_matrix/"))
    print(paste0(ana.dir, sampleName))
    dir.create(paste0(ana.dir, sampleName))
    setwd(paste0(ana.dir, sampleName))

    ###################################

    pbmc.rna <- CreateSeuratObject(a$`Gene Expression`, project = sampleName, min.cells = 50, min.features = 200)
    p1 <- VlnPlot(pbmc.rna, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
    ggsave(p1, filename = "QC1_violin.pdf")
    p2 <- FeatureScatter(pbmc.rna, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    ggsave(p2, filename = "QC1_scatter.pdf")

    ########## remove doublet cells
    # for(i in seq(0,5000,1000)){
    pbmc <- pbmc.rna
    lapply(seq(0, 5000, 1000), function(umi) {
        dir.create(paste0(ana.dir, sampleName, "/", umi))
        setwd(paste0(ana.dir, sampleName, "/", umi))

        pbmc.rna <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & nCount_RNA > umi)
        dim(pbmc.rna)
        p1 <- VlnPlot(pbmc.rna, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
        ggsave(p1, filename = "QC_violin.pdf")
        p2 <- FeatureScatter(pbmc.rna, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
        ggsave(p2, filename = "QC_scatter.pdf")

        pbmc.rna <- NormalizeData(pbmc.rna) %>%
            FindVariableFeatures() %>%
            ScaleData() %>%
            RunPCA()
        pbmc.rna <- get_doublet(pbmc.rna)
        p <- DimPlot(pbmc.rna, group.by = "doublet")
        ggsave(p, filename = "QC_doublet.pdf")
        table(pbmc.rna$doublet) %>%
            data.frame() %>%
            write.csv("doublet.csv", row.names = F)
        pbmc.rna <- subset(pbmc.rna, subset = `doublet` == "Singlet")

        ############# filter low UMI cells
        pbmc.rna <- NormalizeData(pbmc.rna)
        pbmc.rna <- FindVariableFeatures(pbmc.rna)
        pbmc.rna <- ScaleData(pbmc.rna)
        pbmc.rna <- RunPCA(pbmc.rna)
        pbmc.rna <- RunUMAP(pbmc.rna, dims = 1:50, n.neighbors = 30, min.dist = 0.001)
        pbmc.rna <- FindNeighbors(pbmc.rna, dims = 1:50)
        for (i in seq(0.1, 2, 0.2)) {
            pbmc.rna <- FindClusters(pbmc.rna, resolution = i)
        }
        library(clustree)
        p <- clustree(pbmc.rna, prefix = "RNA_snn_res.", node_label = "nFeature_RNA", node_label_aggr = "median") +
            guides(edge_colour = "none", edge_alpha = "none") +
            scale_color_manual(values = colorRampPalette(pal_npg()(9))(20)) +
            theme(legend.position = "bottom") + ggtitle(label = paste0("cell number:", ncol(pbmc.rna), "\n", "UMI cutoff:", umi))
        # p
        ggsave(p, filename = "QC_clustree.pdf", height = 9, width = 7)
        save(pbmc.rna, file = paste0("data.", i, ".Rdata"))
    })
}
############## stop to select suitable resolution
