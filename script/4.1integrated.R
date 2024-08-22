library(Signac)
library(Seurat)
library(dplyr)
library(EnsDb.Osativa.v55)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Osativa.v55)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "ory"
alldata <-
    Read10X(
        "/aggr/all/outs/filtered_feature_bc_matrix/"
    )

rice <- CreateSeuratObject(alldata$`Gene Expression`)

atac_counts <- alldata$Peaks
grange.counts <-
    StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <-
    seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
frag.file <-
    "/aggr/all/outs/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c(":", "-"),
    fragments = frag.file,
    min.cells = 10,
    annotation = annotations
)
rice[["ATAC"]] <- chrom_assay

###### transfer cluster label to integrated data
load("alldata.Rdata")
cluster_names <-
    lapply(names(allmerge.list), function(x) {
        return(
            data.frame(
                cluster = allmerge.list[[x]]$cluster_name,
                cells = paste(gsub("-.+", "", colnames(allmerge.list[[x]])), allmerge.list[[x]]$orig.ident, sep = "-"),
                group = x,
                row.names = NULL
            )
        )
    })
names(cluster_names) <- names(allmerge.list)
cluster_names <- cluster_names %>% do.call(rbind, .)

rice.sub <- subset(rice, subset = cellnames %in% cluster_names$cells)
rice.sub$cluster <-
    cluster_names$cluster[match(colnames(rice.sub), cluster_names$cells)]

run_wnn <- function(obj) {
    DefaultAssay(obj) <- "RNA"
    obj <-
        NormalizeData(obj) %>%
        FindVariableFeatures() %>%
        ScaleData() %>%
        RunPCA(npcs = 50) %>%
        RunUMAP(
            dims = 1:50,
            reduction.name = "umap.rna",
            reduction.key = "rnaUMAP_"
        )
    DefaultAssay(obj) <- "ATAC"
    obj <- RunTFIDF(obj)
    obj <- FindTopFeatures(obj, min.cutoff = "q0")
    obj <- RunSVD(obj)
    obj <- RunUMAP(
        obj,
        reduction = "lsi",
        dims = 2:50,
        reduction.name = "umap.atac",
        reduction.key = "atacUMAP_"
    )
    obj <-
        FindMultiModalNeighbors(obj,
            reduction.list = list("pca", "lsi"),
            dims.list = list(1:50, 2:50)
        )
    obj <-
        RunUMAP(
            obj,
            nn.name = "weighted.nn",
            reduction.name = "wnn.umap",
            reduction.key = "wnnUMAP_"
        )
    obj <-
        FindClusters(
            obj,
            graph.name = "wsnn",
            algorithm = 3,
            resolution = 0.1
        )
    DefaultAssay(obj) <- "ATAC"
    obj <- NucleosomeSignal(obj)
    obj <- TSSEnrichment(obj)
    return(obj)
}
rice.sub <- run_wnn(rice.sub)
source(
    "DensityScatter.R"
)
p <-
    DensityScatter(
        rice.sub,
        x = "nCount_ATAC",
        y = "TSS.enrichment",
        log_x = TRUE,
        quantiles = TRUE
    )
rice.sub <- RunUMAP(
    rice.sub,
    nn.name = "weighted.nn",
    reduction.name = "wnn.umap",
    reduction.key = "wnnUMAP_",
    dims = 1:50,
    n.neighbors = 30,
    n.epochs = 200,
    min.dist = 0.5,
    a = 0.03,
    b = 1.5
)
DimPlot(rice.sub, reduction = "wnn.umap", group.by = c("sample"))
save(rice.sub,file="rice.sub.Rdata")