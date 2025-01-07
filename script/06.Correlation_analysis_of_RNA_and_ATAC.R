library(Seurat)
library(Signac)
library(ggplot2)
library(dplyr)
library(harmony)
library(ggplot2)
library(ggrastr)
library(ggrepel)
library(ggpubr)
library(corrplot)


load("rice.sub.Rdata")
DefaultAssay(rice.sub) <- "RNA"

rice.sub.avg <- AverageExpression(
    rice.sub,
    assays = "RNA",
    features = VarableFeatures(rice.sub)
)
############## RNA correlation
celltypes <- unique(colnames(rice.sub.avg$RNA))
corM <- cor(rice.sub.avg$RNA, method = "pearson")
pdf("RNA.correaltion.pdf", width = 9, height = 9)
corrplot(
    corM[celltypes, celltypes],
    method = "square",
    type = "upper",
    tl.col = "black",
    tl.cex = 0.6,
    is.corr = F,
    col = rev(COL2("RdBu", 100)),
    order = "original", col.lim = c(-1, 1)
)
dev.off()
########  ATAC correlation

DefaultAssay(rice.sub) <- "ATAC"
rice.sub.atac.markers <-
    parallel::mclapply(unique(rice.sub$tissue_cluster), function(x) {
        xx <- FindMarkers(
            rice.sub,
            ident.1 = x,
            only.pos = T,
            test.use = "LR",
            max.cells.per.ident = 300,
            latent.vars = "nCount_ATAC"
        )
        return(data.frame(xx, gene = rownames(xx), cluster = x))
    }, mc.cores = 20)
peaks <-
    lapply(rice.sub.atac.markers, function(x) {
        if (!is.null(nrow(x))) {
            return(x)
        }
    })
peaksID <-
    do.call(rbind, peaks) %>%
    group_by(cluster) %>%
    top_n(500, avg_log2FC) %>%
    pull(gene) %>%
    unique()

rice.sub.atac.avg <- AverageExpression(rice.sub,
    assays = "ATAC",
    features = peaksID
)
corATAC <- cor(rice.sub.atac.avg$ATAC, method = "pearson")
pdf("ATAC.correaltion.pdf", width = 9, height = 9)
corrplot(
    corATAC[celltypes, celltypes],
    method = "square",
    type = "lower",
    tl.col = "black",
    tl.cex = 0.6,
    is.corr = F,
    col = rev(COL2("RdBu", 100)),
    order = "original", col.lim = c(-1, 1)
)
dev.off()

############## RNA and ATAC correlation

DefaultAssay(rice.sub) <- "ATAC"
gene.activities <- GeneActivity(rice.sub)
rice.sub[["atacRNA"]] <- CreateAssayObject(counts = gene.activities)
rice.sub <- NormalizeData(
    object = rice.sub,
    assay = "atacRNA",
    normalization.method = "LogNormalize",
    scale.factor = median(rice.sub$nCount_atacRNA)
)

atacRNA.avg <- AverageExpression(
    rice.sub,
    assays = "atacRNA",
    features = VarableFeatures(rice.sub)
)
colnames(atacRNA.avg$atacRNA) <-
    paste(colnames(atacRNA.avg$atacRNA), "ATAC", sep = "-")
atac_RNA.cor <-
    cor(atacRNA.avg$atacRNA, rice.sub.avg$RNA[rownames(atacRNA.avg$atacRNA), ], method = "spearman")

pheatmap::pheatmap(
    atac_RNA.cor,
    cluster_cols = F,
    cluster_rows = F,
    filename = "RNA_ATAC.correlation.pdf",
    height = 9,
    width = 11
)
atac_RNA.cor %>%
    reshape2::melt() %>%
    write.csv("atac_RNA.cor.csv", row.names = F)
