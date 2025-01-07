library(Seurat)
library(Signac)
library(ggplot2)
library(dplyr)
library(motifmatchr)
library(TFBSTools)
library(EnsDb.Osativa.v55)
library(BSgenome.O.sativa.Ensembl.MSU7)


load("rice.sub.Rdata")
pfm <-
    readRDS("rice.motif.pfmList2.Rds")
motif.family <-
    lapply(pfm, function(x) {
        return(x@tags$family)
    }) %>%
    unlist() %>%
    data.frame(family = ., motif = names(.))
DefaultAssay(rice.sub) <- "ATAC"
rice.sub <- AddMotifs(
    object = rice.sub,
    genome = BSgenome.O.sativa.Ensembl.MSU7,
    pfm = pfm
)
rice.sub <- RunChromVAR(
    object = rice.sub,
    genome = BSgenome.O.sativa.Ensembl.MSU7
)

rice.sub.atac.markers <-
    parallel::mclapply(unique(rice.sub$tissue_cluster), function(x) {
        print(x)
        xx <-
            FindMarkers(
                rice.sub,
                ident.1 = x,
                only.pos = T,
                test.use = "LR",
                logfc.threshold = 0.05,
                max.cells.per.ident = 300,
                latent.vars = "nCount_ATAC"
            )

        if (nrow(xx) > 0) {
            return(data.frame(xx, gene = rownames(xx), cluster = x))
        }
    }, mc.cores = 10)

names(rice.sub.atac.markers) <- unique(rice.sub$tissue_cluster)
rice.sub.atac.markers <-
    rice.sub.atac.markers[!sapply(rice.sub.atac.markers, is.null)]

get_motif <- function(rice.sub, cluster = "") {
    top.da.peak <-
        rice.sub.atac.markers[[cluster]] %>%
        dplyr::filter(p_val < 0.05, avg_log2FC >
            0.5) %>%
        pull(gene)
    if (length(top.da.peak) < 10) {
        top.da.peak <-
            rice.sub.atac.markers[[cluster]] %>%
            top_n(300, avg_log2FC) %>%
            pull(gene)
    }
    open.peaks <- AccessiblePeaks(rice.sub, idents = cluster)
    # match the overall GC content in the peak set
    meta.feature <-
        GetAssayData(rice.sub, assay = "ATAC", slot = "meta.features")
    peaks.matched <- MatchRegionStats(
        meta.feature = meta.feature[open.peaks, ],
        query.feature = meta.feature[top.da.peak, ],
        n = 50000
    )
    enriched.motifs <- FindMotifs(
        object = rice.sub,
        features = top.da.peak,
        background = peaks.matched
    ) %>% mutate(family = gsub(":.+", "", motif))
    enriched.motifs$cluster <- cluster
    return(enriched.motifs)
}

allmotifs <-
    lapply(names(rice.sub.atac.markers), function(x) {
        get_motif(rice.sub, cluster = x)
    })
allmotifs <- do.call(rbind, allmotifs)
DefaultAssay(rice.sub) <- "chromvar"

DefaultAssay(rice.sub) <- "ATAC"

pdf("Dof.CoveragePlot.pdf", 8, 6)
for (xx in unique(rice.sub$tissue)) {
    p <- CoveragePlot(
        object = rice.sub,
        region = "Os03g0821200",
        # Dof:LOC-Os03g60630
        features = "Os03g0821200",
        expression.assay = "RNA",
        extend.upstream = 5000,
        extend.downstream = 5000,
        idents = grep(xx, unique(rice.sub$tissue_cluster), value = T),
        annotation = T
    )
    print(cowplot::plot_grid(p, labels = xx))
}
dev.off()
