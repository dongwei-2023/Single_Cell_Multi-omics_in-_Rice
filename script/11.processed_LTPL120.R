rm(list=ls())
gc()
options(stringsAsFactors = F) 
setwd("/data2/lidongwei/work/work_scMultiome/part_mut_LTPL120/results/R/")
library(Seurat)
library(future)
library(harmony)
library(readr)
library(DoubletFinder)
library(ggplot2)
library(dplyr)
library(clustree)
library(ggsci)
options(future.globals.maxSize = 100 * 1024^3) 
plan("multicore", workers=20)

colors <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#1F78B4", "#33A02C")
#############################################

library(rtracklayer)
rice_gtf <- import("/data2/lidongwei/work/work_scMultiome/refGenome/genes.gtf")
mt_gtf <- rice_gtf[which(rice_gtf@seqnames == "chrMt")]
rice_mt_gene <- unique(mt_gtf$gene_id)

paramSweep_v3_Seurat5<-function (seu, PCs = 1:10, sct = FALSE, num.cores = 1) 
{
  require(Seurat)
  require(fields)
  pK <- c(5e-04, 0.001, 0.005, seq(0.01, 0.3, by = 0.01))
  pN <- seq(0.05, 0.3, by = 0.05)
  min.cells <- round(nrow(seu@meta.data)/(1 - 0.05) - nrow(seu@meta.data))
  pK.test <- round(pK * min.cells)
  pK <- pK[which(pK.test >= 1)]
  orig.commands <- seu@commands
  if (nrow(seu@meta.data) > 10000) {
    real.cells <- rownames(seu@meta.data)[sample(1:nrow(seu@meta.data), 
                                                 10000, replace = FALSE)]
    data <- seu@assays$RNA$counts[, real.cells]
    n.real.cells <- ncol(data)
  }
  if (nrow(seu@meta.data) <= 10000) {
    real.cells <- rownames(seu@meta.data)
    data <- seu@assays$RNA$counts
    n.real.cells <- ncol(data)
  }
  if (num.cores > 1) {
    require(parallel)
    cl <- makeCluster(num.cores)
    output2 <- mclapply(as.list(1:length(pN)), FUN = parallel_paramSweep_v3, 
                        n.real.cells, real.cells, pK, pN, data, orig.commands, 
                        PCs, sct, mc.cores = num.cores)
    stopCluster(cl)
  }
  else {
    output2 <- lapply(as.list(1:length(pN)), FUN = parallel_paramSweep_v3, 
                      n.real.cells, real.cells, pK, pN, data, orig.commands, 
                      PCs, sct)
  }
  sweep.res.list <- list()
  list.ind <- 0
  for (i in 1:length(output2)) {
    for (j in 1:length(output2[[i]])) {
      list.ind <- list.ind + 1
      sweep.res.list[[list.ind]] <- output2[[i]][[j]]
    }
  }
  name.vec <- NULL
  for (j in 1:length(pN)) {
    name.vec <- c(name.vec, paste("pN", pN[j], "pK", pK, 
                                  sep = "_"))
  }
  names(sweep.res.list) <- name.vec
  return(sweep.res.list)
}



doubletFinder_v3_v5 <- function (seu, PCs, pN = 0.25, pK, nExp, reuse.pANN = FALSE, 
                                 sct = FALSE, annotations = NULL) 
{
  require(Seurat)
  require(fields)
  require(KernSmooth)
  if (reuse.pANN != FALSE) {
    pANN.old <- seu@meta.data[, reuse.pANN]
    classifications <- rep("Singlet", length(pANN.old))
    classifications[order(pANN.old, decreasing = TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("DF.classifications", pN, pK, nExp, 
                          sep = "_")] <- classifications
    return(seu)
  }
  if (reuse.pANN == FALSE) {
    real.cells <- rownames(seu@meta.data)
    data <- seu[["RNA"]]$counts[, real.cells]
    n_real.cells <- length(real.cells)
    n_doublets <- round(n_real.cells/(1 - pN) - n_real.cells)
    print(paste("Creating", n_doublets, "artificial doublets...", 
                sep = " "))
    real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
    real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
    doublets <- (data[, real.cells1] + data[, real.cells2])/2
    colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
    data_wdoublets <- cbind(data, doublets)
    if (!is.null(annotations)) {
      stopifnot(typeof(annotations) == "character")
      stopifnot(length(annotations) == length(Cells(seu)))
      stopifnot(!any(is.na(annotations)))
      annotations <- factor(annotations)
      names(annotations) <- Cells(seu)
      doublet_types1 <- annotations[real.cells1]
      doublet_types2 <- annotations[real.cells2]
    }
    orig.commands <- seu@commands
    if (sct == FALSE) {
      print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
      print("Normalizing Seurat object...")
      seu_wdoublets <- NormalizeData(seu_wdoublets, normalization.method = orig.commands$NormalizeData.RNA@params$normalization.method, 
                                     scale.factor = orig.commands$NormalizeData.RNA@params$scale.factor, 
                                     margin = orig.commands$NormalizeData.RNA@params$margin)
      print("Finding variable genes...")
      seu_wdoublets <- FindVariableFeatures(seu_wdoublets, 
                                            selection.method = orig.commands$FindVariableFeatures.RNA$selection.method, 
                                            loess.span = orig.commands$FindVariableFeatures.RNA$loess.span, 
                                            clip.max = orig.commands$FindVariableFeatures.RNA$clip.max, 
                                            mean.function = orig.commands$FindVariableFeatures.RNA$mean.function, 
                                            dispersion.function = orig.commands$FindVariableFeatures.RNA$dispersion.function, 
                                            num.bin = orig.commands$FindVariableFeatures.RNA$num.bin, 
                                            binning.method = orig.commands$FindVariableFeatures.RNA$binning.method, 
                                            nfeatures = orig.commands$FindVariableFeatures.RNA$nfeatures, 
                                            mean.cutoff = orig.commands$FindVariableFeatures.RNA$mean.cutoff, 
                                            dispersion.cutoff = orig.commands$FindVariableFeatures.RNA$dispersion.cutoff)
      print("Scaling data...")
      seu_wdoublets <- ScaleData(seu_wdoublets, features = orig.commands$ScaleData.RNA$features, 
                                 model.use = orig.commands$ScaleData.RNA$model.use, 
                                 do.scale = orig.commands$ScaleData.RNA$do.scale, 
                                 do.center = orig.commands$ScaleData.RNA$do.center, 
                                 scale.max = orig.commands$ScaleData.RNA$scale.max, 
                                 block.size = orig.commands$ScaleData.RNA$block.size, 
                                 min.cells.to.block = orig.commands$ScaleData.RNA$min.cells.to.block)
      print("Running PCA...")
      seu_wdoublets <- RunPCA(seu_wdoublets, features = orig.commands$ScaleData.RNA$features, 
                              npcs = length(PCs), rev.pca = orig.commands$RunPCA.RNA$rev.pca, 
                              weight.by.var = orig.commands$RunPCA.RNA$weight.by.var, 
                              verbose = FALSE)
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[, 
                                                                PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets)
      gc()
    }
    if (sct == TRUE) {
      require(sctransform)
      print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
      print("Running SCTransform...")
      seu_wdoublets <- SCTransform(seu_wdoublets)
      print("Running PCA...")
      seu_wdoublets <- RunPCA(seu_wdoublets, npcs = length(PCs))
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[, 
                                                                PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets)
      gc()
    }
    print("Calculating PC distance matrix...")
    dist.mat <- fields::rdist(pca.coord)
    print("Computing pANN...")
    pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, 
                                 ncol = 1))
    if (!is.null(annotations)) {
      neighbor_types <- as.data.frame(matrix(0L, nrow = n_real.cells, 
                                             ncol = length(levels(doublet_types1))))
    }
    rownames(pANN) <- real.cells
    colnames(pANN) <- "pANN"
    k <- round(nCells * pK)
    for (i in 1:n_real.cells) {
      neighbors <- order(dist.mat[, i])
      neighbors <- neighbors[2:(k + 1)]
      pANN$pANN[i] <- length(which(neighbors > n_real.cells))/k
      if (!is.null(annotations)) {
        for (ct in unique(annotations)) {
          neighbors_that_are_doublets = neighbors[neighbors > 
                                                    n_real.cells]
          if (length(neighbors_that_are_doublets) > 0) {
            neighbor_types[i, ] <- table(doublet_types1[neighbors_that_are_doublets - 
                                                          n_real.cells]) + table(doublet_types2[neighbors_that_are_doublets - 
                                                                                                  n_real.cells])
            neighbor_types[i, ] <- neighbor_types[i, 
            ]/sum(neighbor_types[i, ])
          }
          else {
            neighbor_types[i, ] <- NA
          }
        }
      }
    }
    print("Classifying doublets..")
    classifications <- rep("Singlet", n_real.cells)
    classifications[order(pANN$pANN[1:n_real.cells], decreasing = TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("pANN", pN, pK, nExp, sep = "_")] <- pANN[rownames(seu@meta.data), 
                                                                    1]
    seu@meta.data[, paste("DF.classifications", pN, pK, nExp, 
                          sep = "_")] <- classifications
    if (!is.null(annotations)) {
      colnames(neighbor_types) = levels(doublet_types1)
      for (ct in levels(doublet_types1)) {
        seu@meta.data[, paste("DF.doublet.contributors", 
                              pN, pK, nExp, ct, sep = "_")] <- neighbor_types[, 
                                                                              ct]
      }
    }
    return(seu)
  }
}


rm_doublet_v5 <- function(EC = EC,dim.usage=30,project_names=NULL,n_cores = 1) {
  
  EC <- NormalizeData(EC, normalization.method = "LogNormalize", scale.factor = 10000)
  EC <- FindVariableFeatures(EC, selection.method = "vst", nfeatures = 2000)
  EC <- ScaleData(EC,vars.to.regress = c("percent.mt"))
  EC <- RunPCA(EC, features = VariableFeatures(object = EC))
  EC <- FindNeighbors(EC, dims = 1:dim.usage)
  EC <- RunUMAP(EC, dims = 1:dim.usage)
  EC <- FindClusters(EC, resolution = 0.8)
  
 
  Find_doublet <- function(data,dim.usage_1=30){
    sweep.res.list <- paramSweep_v3_Seurat5(data, PCs = 1:dim.usage_1, sct=FALSE,num.cores = n_cores)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    p<-as.numeric(as.vector(bcmvn[bcmvn$MeanBC==max(bcmvn$MeanBC),]$pK))
    
    homotypic.prop <- modelHomotypic(data@meta.data$seurat_clusters) 
    Doubletrate <- 0.05
    pN_value <- 0.25
    nExp_poi <- round(Doubletrate*ncol(data))
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    
    data <- doubletFinder_v3_v5(data, PCs = 1:dim.usage_1, pN = pN_value, pK = p, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = F)  # sct：表示是否使用了SCTransform方法进行标准化
    colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_info"
    data
  }
  EC<-Find_doublet(EC)
  EC<-subset(EC,subset=doublet_info=="Singlet")
  c <- grep("pANN_",colnames(EC@meta.data))
  EC@meta.data <- EC@meta.data[,-c]
  print(paste0(" cells: ", length(EC@meta.data$orig.ident)," is read!"," Time:",format(Sys.time(), "%Y%m%d %X"),sep = " "))
  return(EC)
}


run_QC <- function(matrix_path=NULL,nFeature_RNA_NUM=500,nCount_RNA_num = 1000,mt_genes= rice_mt_gene,percent_mt= 10, out_path="./s1_QC/",sample_name = "seurat"){
  
  seurat_count <- Read10X(data.dir = matrix_path) 
  seurat_obj <- CreateSeuratObject(counts = seurat_count, project = sample_name, min.cells=3, min.features = 200) 
  rm(seurat_count)
  mt_genes_use_seurat_obj <- intersect(rownames(seurat_obj), mt_genes)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, features=mt_genes_use_seurat_obj)
  p1 <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(filename=paste0(out_path,sample_name,"_","Vlnplot_before_QC.pdf"),plot=p1, width = 12, height = 8, units = "cm")
  
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > nFeature_RNA_NUM  & percent.mt < percent_mt & nCount_RNA > nCount_RNA_num)
  p2 <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  ggsave(filename=paste0(out_path,sample_name,"_","Vlnplot_after_QC.pdf"),plot=p2, width = 12, height = 8, units = "cm")
  
  seurat_obj <- rm_doublet_v5(seurat_obj,project_names= sample_name,n_cores=30)
  write_rds(seurat_obj,paste0(out_path,sample_name,"_","after_rmDOublet.rds"))
  return(seurat_obj)
}


###################################################################################################################
##################################  s1 QC #############################################
###################################################################################################################
### Create output directory if it doesn't exist
out_dir <- "./s1_QC"
if(!file.exists(out_dir)){
  dir.create(out_dir, showWarnings = TRUE, recursive=TRUE)
}

# LTP_G_1
LTP_G_1 <- run_QC(matrix_path = "/data2/lidongwei/work/work_scMultiome/part_mut_LTPL120/results/LTP_G_1/outs/filtered_feature_bc_matrix/",
                  nFeature_RNA_NUM=500,
                  nCount_RNA_num = 1000,
                  mt_genes= rice_mt_gene,
                  percent_mt= 10, 
                  out_path="./s1_QC/",
                  sample_name = "LTP_G_1")
# LTP_G_2
LTP_G_2 <- run_QC(matrix_path = "/data2/lidongwei/work/work_scMultiome/part_mut_LTPL120/results/LTP_G_2/outs/filtered_feature_bc_matrix/",
                  nFeature_RNA_NUM=500,
                  nCount_RNA_num = 1000,
                  mt_genes= rice_mt_gene,
                  percent_mt= 10, 
                  out_path="./s1_QC/",
                  sample_name = "LTP_G_2")

# NIP_G_1
NIP_G_1 <- run_QC(matrix_path = "/data2/lidongwei/work/work_scMultiome/part_mut_LTPL120/results/NIP_G_1/outs/filtered_feature_bc_matrix/",
                  nFeature_RNA_NUM=500,
                  nCount_RNA_num = 1000,
                  mt_genes= rice_mt_gene,
                  percent_mt= 10, 
                  out_path="./s1_QC/",
                  sample_name = "NIP_G_1")

# NIP_G_2
NIP_G_2 <- run_QC(matrix_path = "/data2/lidongwei/work/work_scMultiome/part_mut_LTPL120/results/NIP_G_2/outs/filtered_feature_bc_matrix/",
                  nFeature_RNA_NUM=500,
                  nCount_RNA_num = 1000,
                  mt_genes= rice_mt_gene,
                  percent_mt= 10, 
                  out_path="./s1_QC/",
                  sample_name = "NIP_G_2")






LTP_NIP_list <- list("LTP_G_1" = LTP_G_1,
                     "LTP_G_2" = LTP_G_2,
                     "NIP_G_1" = NIP_G_1,
                     "NIP_G_2" = NIP_G_2)
do.call(rbind,lapply(LTP_NIP_list, dim))
LTP_NIP <- merge(LTP_NIP_list[[1]], y = LTP_NIP_list[2:length(LTP_NIP_list)],
             add.cell.ids = names(LTP_NIP_list),project="LTP_NIP")
names(LTP_NIP@assays$RNA@layers)

LTP_NIP$condition <- ifelse(str_detect(LTP_NIP@meta.data$orig.ident, "^NIP"),
                            "WT","LTP")


write_rds(LTP_NIP,"./s1_QC/LTP_NIP_merge_list.rds")

LTP_NIP[["RNA"]]$counts 
# Alternate accessor function with the same result
LayerData(LTP_NIP, assay = "RNA", layer = "counts")
LTP_NIP <- JoinLayers(LTP_NIP)

write_rds(LTP_NIP,"./s1_QC/LTP_NIP_merge_list_JoinLayers.rds")


LTP_NIP@meta.data %>%
  ggplot(aes(x=orig.ident, fill=orig.ident)) +
  geom_bar(color="black") +
  stat_count(geom = "text", colour = "black", size = 3.5,
             aes(label = ..count..),
             position=position_stack(vjust=0.5))+
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Number of Cells per Sample")


###################################################################################################################
##################################  s2 harmony #############################################
###################################################################################################################

out_dir <- "./s2_harmony"
if(!file.exists(out_dir)){
  dir.create(out_dir, showWarnings = TRUE, recursive=TRUE)
}


LTP_NIP <- NormalizeData(LTP_NIP, 
                              normalization.method = "LogNormalize",
                              scale.factor = 1e4) 

LTP_NIP <- FindVariableFeatures(LTP_NIP)
p4 <- VariableFeaturePlot(LTP_NIP) 
p4

LTP_NIP <- ScaleData(LTP_NIP)

LTP_NIP <- RunPCA(LTP_NIP, features = VariableFeatures(object = LTP_NIP))

seuratObj <- RunHarmony(LTP_NIP, "orig.ident")

seuratObj <- RunUMAP(seuratObj,  dims = 1:20, reduction = "harmony")
DimPlot(seuratObj,reduction = "umap",label=F ) 
seuratObj <- RunTSNE(seuratObj, dims = 1:20, reduction = "harmony")
DimPlot(seuratObj,reduction = "tsne",label=F ) 
seuratObj <- FindNeighbors(seuratObj, reduction = "harmony", dims = 1:20) 

seuratObj <- FindClusters(object = seuratObj , resolution = seq(from = 0.1, to = 1.0,  by = 0.1))

write_rds(seuratObj,"./s2_harmony/LTP_NIP_after_harmony_FindClusters.rds")

DimPlot(object = seuratObj, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
DimPlot(object = seuratObj, reduction = "umap", pt.size = .1, group.by = "orig.ident")

p2_tree=clustree(seuratObj@meta.data, prefix = "RNA_snn_res.")
ggsave(plot=p2_tree, filename="./s2_harmony/Tree_diff_resolution.pdf")



###################################################################################################################
##################################  s3 marker gene #############################################
###################################################################################################################

out_dir <- "./s3_marker_gene"
if(!file.exists(out_dir)){
  dir.create(out_dir, showWarnings = TRUE, recursive=TRUE)
}

sel.clust = "RNA_snn_res.0.5"
seuratObj <- SetIdent(seuratObj, value = sel.clust)

# find markers for every cluster compared to all remaining cells, report only the positive ones
seuratObj.markers <- FindAllMarkers(seuratObj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

write.csv(seuratObj.markers, file="./s3_marker_gene/LTP_NIP_all_res_0.5_clusters.csv")
seuratObj.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)

DimPlot(seuratObj, reduction = "umap", 
        group.by = "RNA_snn_res.0.5", pt.size = .3,label = TRUE,label.size = 7) +
  scale_color_npg()+
  NoLegend()


seurat_inte_label <- seuratObj
seurat_inte_label@meta.data$seurat_clusters <- seurat_inte_label@meta.data$RNA_snn_res.0.5

cluster2celltype <- c("0"="cortex",
                      "1"="cortex", 
                      "2"="root cap", 
                      "3"= "vascular", 
                      "4"= "meristem", 
                      "5"= "exodermis",
                      "6"= "epidermis (neat root hair)", 
                      "7"= "sclerenchyma", 
                      "8"= "vascular", 
                      "9"= "epidermis", 
                      "10"= "epidermis", 
                      "11"= "endodermis", 
                      "12"="exodermis", 
                      "13"= "vascular", 
                      "14"= "root cap", 
                      "15"= "vascular",
                      "16"= "exodermis", 
                      "17"= "exodermis", 
                      "18"= "vascular", 
                      "19"= "epidermis (neat root hair)", 
                      "20"= "epidermis (neat root hair)")


seurat_inte_label[['cell_type']] = unname(cluster2celltype[seurat_inte_label@meta.data$seurat_clusters])
Idents(seurat_inte_label) <- "cell_type"

write_rds(seurat_inte_label,"./s3_marker_gene/LTP_NIP_label_240622.rds")

sel.clust = "cell_type"
seurat_inte_label <- SetIdent(seurat_inte_label, value = sel.clust)

# find markers for every cluster compared to all remaining cells, report only the positive ones
seuratObj.markers <- FindAllMarkers(seurat_inte_label, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

write.csv(seuratObj.markers, file="./s3_marker_gene/LTP_NIP_all_res_cell_type.csv")
seuratObj.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)


###################################################################################################################
##################################  s4_compare_LTP_NIP #############################################
###################################################################################################################
col2 = pal_igv("default")(33)
custom_colors <- col2


out_dir <- "./s4_compare_LTP_NIP_cluster"
if(!file.exists(out_dir)){
  dir.create(out_dir, showWarnings = TRUE, recursive=TRUE)
}


p1=DimPlot(seurat_inte_label, split.by = "condition")

ct_stat2 = as.data.frame(table(seurat_inte_label$cell_type, seurat_inte_label$condition))
sums = rep(c(sum(ct_stat2$Freq[1:9]),sum(ct_stat2$Freq[10:18])),each=9)
ct_stat2$percentage = ct_stat2$Freq/sums

ct_stat2$Var2 = factor(ct_stat2$Var2, levels = c("LTP","WT"))
library(ggplot2)
p2=ggplot(ct_stat2, aes(x=Var1, fill=Var2, y = percentage)) +
  geom_bar( position=position_dodge(width=0.8), width=0.6,
            stat="identity",color = "black")  + 
  scale_fill_manual(values= c("#DC143C","#1E90FF")) +
  theme_classic() + ylab("Proportion of total nuclei(%)") + 
  theme(axis.title.x = element_blank(), legend.title = element_blank())

seurat_inte_label$compare = paste(seurat_inte_label$condition, seurat_inte_label$cell_type, sep = "_")
ct = levels(seurat_inte_label@active.ident)
all_markers = lapply(ct, function(x){
  print(x)
  markers <- FindMarkers(seurat_inte_label, group.by = "compare",
                         logfc.threshold = 0.1,
                         ident.1 = paste("LTP",x,sep = "_"),
                         ident.2 = paste("WT",x,sep = "_"))
  return(markers)
})

lapply(all_markers,nrow)

for (i in 1:length(all_markers)){
  all_markers[[i]]$gene <- rownames(all_markers[[i]])
}



saveRDS(all_markers,"./s4_compare_LTP_NIP_cluster/all_markers_compare.rds")

library(xlsx)
write.xlsx(all_markers, "s4_compare_LTP_NIP_cluster/all_markers_compare.xlsx", rowNames = FALSE)



all_markers_sig = lapply(all_markers, function(x){
  markers_sig <- subset(x, p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)
})
sapply(all_markers_sig,nrow)
marker_stat = as.data.frame(lapply(all_markers_sig,function(x){
  Up = sum(x$avg_log2FC>0)
  Down = sum(x$avg_log2FC<0)
  Total = length(x$avg_log2FC)
  return(c(Up, Down, Total))
}))
rownames(marker_stat) = c("Up","Down","Total")

sce <- subset(x = seurat_inte_label, features = all_markers_sig$cortex$gene) 
nCount = colSums(x = sce, layer = "counts")  
nFeature = colSums(x = GetAssayData(object = sce, layer  = "counts") > 0) 
sce$nFeature_RNA <- nFeature
sce$nCount_RNA <- nCount

FeaturePlot(sce,features="nFeature_RNA",reduction = 'umap',cols = c("#EAEFD7" ,"#17645C"),raster=FALSE)
ggsave("./s4_compare_LTP_NIP_cluster/cortex_gene_num_umap.pdf")

