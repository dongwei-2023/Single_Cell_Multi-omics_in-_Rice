# Single Cell Multiomics in Rice analysis code

**Accompanying code repository for the paper entitled: "A single-cell multiomics atlas of rice".**

## **Abstract**
Cell functions across eucaryotes are driven by specific gene expression programs that are dependent on chromatin structure. Here, we report a single-cell multiomics atlas of rice, one of the major crop plants. We simultaneously profiled chromatin accessibility and RNA expression in 116,564 cells from eight organs. We annotated the majority of cell types and identified cell-type expression programs, and found a high correlation between RNA and chromatin accessibility across cells. We constructed single-cell gene regulatory networks and co-expression networks, which allowed us to identify the cell type-specific regulators RSR1, F3H and LTPL120 in rice development. Furthermore, our analysis revealed a correlation between cell type and agronomic traits of rice, and the conserved and divergent functions of cell type during evolution. In summary, this study not only offers a valuable single-cell multiomics resource, but also enriches our understanding of the intricate roles and molecular underpinnings of individual cell types in rice.

## **Schematic overview of the workflow**
![image](https://github.com/dongwei-2023/Single_Cell_Multiomics_in_Rice/blob/master/img/pipeline.png)

## Raw and processed data 
* The raw and processed data of single-cell multiomics, scRNA-seq and RNA-seq have been deposited into NCBI GEO with accession numbers: [GSE232863](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE232863), [GSE273875](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE273875), [GSE245410](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245410). Processed single-cell multiomics seurat objects are available from [GSE232863_rice_scRNA_scATAC_all.Rdata.gz](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE232863&format=file&file=GSE232863%5Frice%5FscRNA%5FscATAC%5Fall%2ERdata%2Egz) and processed LTPL120 seurat objects are available from [GSE273875_LTP_NIP_label_240622.rds.gz](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE273875&format=file&file=GSE273875%5FLTP%5FNIP%5Flabel%5F240622%2Erds%2Egz).The raw and processed data of single-cell multiomics of Panicle have been deposited into NCBI GEO with accession numbers:[GSE285639](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE285639).

* Single Cell Multi-omics in Rice chromatin accessibility and RNA expression zonation profiles can be browsed via our web-app: https://www.elabcaas.cn/scmr/

* Processed zonation tables are provided as supplementary tables in the manuscript. 

# Scripts for single-cell multiomics analyses

This repository contains a series of scripts designed to handle various aspects of single-cell multiomics data analysis for rice. The scripts cover the entire workflow, including quality control, doublet removal, data integration, clustering, cell type annotation, differential expression analysis, motif and regulatory network analysis, and visualization. The analysis leverages tools like Seurat, Signac, DoubletFinder, Harmony, WGCNA, cellOracle and Monocle to process and interpret single-cell RNA-seq and ATAC-seq data. The outputs include UMAP visualizations, marker gene identification, pseudotime trajectory analysis, motif enrichment results, and cross-species comparisons, providing a comprehensive pipeline for single-cell multiomics studies in rice.

## Overview of Scripts

### 1.preprocess.sh

- **Purpose**:This script processes single-cell multiome data (RNA-seq + ATAC-seq) using **Cell Ranger ARC** to generate count matrices for multiple samples.

- **Inputs**:
  - **Cell Ranger ARC executable**
  - **Reference genome**
  - **Sample names**
  - **Raw FASTQ files**

- **Outputs**:
  - **Count matrices** for gene expression and chromatin accessibility
  - **Summary metrics** and analysis reports

### 2.QC.R

- **Purpose**: This script processes single-cell RNA-seq data to perform quality control, doublet detection, and clustering analysis using **Seurat** and **DoubletFinder**.

- **Inputs**:
  - **Raw single-cell RNA-seq data** (filtered feature-barcode matrices) for multiple samples.
  - **Doublet rate information** from a file (`double_rate.txt`).

- **Outputs**:
  - **Quality control plots** (violin plots, scatter plots) for cell filtering.
  - **Doublet detection results** and filtered single-cell data.
  - **Clustering results** and visualizations (UMAP, clustree plots).
  - **Processed Seurat objects** saved as `.Rdata` files for downstream analysis.

### 3.merge_rep.R

- **Purpose**: This script processes and analyzes single-cell RNA-seq data across multiple samples and tissues using **Seurat**. It performs clustering, marker gene identification, and visualization, followed by data integration and dimensionality reduction.

- **Inputs**:
  - **Processed Seurat objects** for each sample (loaded from `.Rdata` files).
  - **Gene list** for marker gene annotation.

- **Outputs**:
  - **UMAP plots** for each sample.
  - **Marker gene lists** saved as CSV files.
  - **Dot plots** and **feature plots** for top marker genes.
  - **Integrated Seurat objects** for each tissue type.
  - **Processed data** saved as `.Rdata` files for downstream analysis.




### 4.1 integrated.R

- **Purpose**: This script processes single-cell multiome data (RNA-seq + ATAC-seq) using **Signac** and **Seurat**. It integrates gene expression and chromatin accessibility data, performs clustering, and generates visualizations.

- **Inputs**:
  - **Filtered feature-barcode matrices** from 10X multiome data.
  - **ATAC fragment file** and genome annotations.
  - **Pre-processed cluster labels** from `alldata.Rdata`.

- **Outputs**:
  - **Integrated Seurat object** with RNA-seq and ATAC-seq data.
  - **UMAP visualizations** for RNA-seq, ATAC-seq, and weighted nearest neighbor (WNN) integration.
  - **Quality control metrics** (nucleosome signal, TSS enrichment).
  - **Processed data** saved as `rice.sub.Rdata` for downstream analysis.


### 4.2 integrated.R

- **Purpose**: This script integrates single-cell RNA-seq and ATAC-seq data using **Harmony** and **Seurat**, performs joint clustering, and generates visualizations for tissue and cluster annotations.

- **Inputs**:
  - **Pre-processed Seurat object** (`rice.sub.Rdata`) containing RNA-seq and ATAC-seq data.
  - **Tissue and cluster annotations**.

- **Outputs**:
  - **Integrated UMAP visualizations** for RNA-seq, ATAC-seq, and joint WNN analysis.
  - **Cluster and tissue-specific plots** (PDF files).
  - **CSV files** containing cell embeddings, cluster IDs, and tissue annotations.


### 5.correlation.R

- **Purpose**: This script performs correlation analysis on single-cell RNA-seq and ATAC-seq data to explore relationships between gene expression and chromatin accessibility.

- **Inputs**:
  - **Pre-processed Seurat object** (`rice.sub.Rdata`) containing RNA-seq and ATAC-seq data.
  - **Variable features** for RNA-seq and ATAC-seq data.

- **Outputs**:
  - **Correlation heatmaps** for RNA-seq, ATAC-seq, and RNA-ATAC integration (PDF files).
  - **CSV files** containing correlation values between RNA-seq and ATAC-seq data.


### 6.1_motif_TF.R

- **Purpose**: This script performs motif analysis and visualization for single-cell ATAC-seq data using **Signac** and **motifmatchr**. It identifies enriched transcription factor motifs, computes chromatin accessibility profiles, and generates coverage plots for specific genes.

- **Inputs**:
  - **Pre-processed Seurat object** (`rice.sub.Rdata`) containing ATAC-seq data.
  - **Position frequency matrices (PFMs)** for motif analysis (`rice.motif.pfmList2.Rds`).

- **Outputs**:
  - **Enriched motifs** for each tissue or cluster.
  - **Coverage plots** for specific genes (e.g., `Os03g0821200` and `Os05g0121600`) across tissues (PDF files).


### 6.2_run_rice_create_cistarget.sh

- **Purpose**: This script creates a **cisTarget motif database** for rice using the `create_cisTarget_databases` pipeline. It processes genomic regions and motif data to generate a database for regulatory network analysis.

- **Inputs**:
  - **FASTA file** containing genomic regions (e.g., gene TSS ± 3 kb).
  - **Motif files** in `cb` format and a list of motifs.

- **Outputs**:
  - **cisTarget motif database** files for use in regulatory network analysis (e.g., SCENIC).


### 6.3_runPySCENIC.sh
- **Purpose**: This script runs the **pySCENIC** pipeline to infer gene regulatory networks (GRNs) from single-cell RNA-seq data. It performs three main steps: GRN inference, motif enrichment analysis, and regulon activity scoring.

- **Inputs**:
  - **Single-cell RNA-seq data** in loom format.
  - **Transcription factor list** and **motif database** for rice.
  - **cisTarget motif rankings** for motif enrichment.

- **Outputs**:
  - **Gene regulatory network** (GRN) adjacency list.
  - **Regulons** with enriched motifs.
  - **Regulon activity scores** in loom format.



### 6.4_RSR1_gene_KO_simulation_with_rice_root.ipynb
- **Purpose**: This script performs **in silico gene perturbation analysis** using **CellOracle** to simulate the effects of gene knockout (e.g., `Os05g0121600`) on cell identity and gene regulatory networks. It visualizes the predicted cell state transitions and vector fields.

- **Inputs**:
  - **Pre-processed CellOracle object** (`rice_root.celloracle.oracle`) containing single-cell RNA-seq data and GRN information.
  - **Gene regulatory network links** (`rice_root_links.celloracle.links`).

- **Outputs**:
  - **Quiver plots** and **vector field graphs** showing predicted cell state transitions after gene perturbation.
  - **Visualizations** of gene expression changes and cell identity shifts. 
 
### 7.leaf_Flag.R
- **Purpose**: This script performs single-cell RNA-seq data analysis using **Seurat** and **Monocle** to study cell type-specific gene expression, pseudotime trajectories, and differential expression between tissues (e.g., leaf and flag). It also integrates data from multiple samples and performs functional enrichment analysis.

- **Inputs**:
  - **Single-cell RNA-seq data** for leaf and flag tissues.
  - **Pre-processed Seurat object** (`alldata.Rdata`) containing integrated data.

- **Outputs**:
  - **UMAP visualizations** for cell clusters and tissue-specific expression.
  - **Differential expression analysis** results and functional enrichment plots.
  - **Pseudotime trajectory analysis** using Monocle.
  - **Integrated Seurat object** and processed data saved for downstream analysis. 
 
### 8.hdWGCNA.R
- **Purpose**: This script performs **Weighted Gene Co-expression Network Analysis (WGCNA)** using the **hdWGCNA** package to identify co-expression modules and analyze their functional relevance in single-cell RNA-seq data.

- **Inputs**:
  - **Pre-processed Seurat object** (`alldata.Rdata`) containing single-cell RNA-seq data.
  - **Gene expression data** for constructing co-expression networks.

- **Outputs**:
  - **Co-expression modules** and their associated genes.
  - **Module eigengene visualizations** and functional enrichment analysis.
  - **Heatmaps** showing module correlations and gene expression patterns. 
 
### 9.Epidermis_Vascular_monocle3.R
- **Purpose**: This script performs **pseudotime trajectory analysis** using **Monocle3** to study cell differentiation and developmental trajectories in single-cell RNA-seq data.

- **Inputs**:
  - **Single-cell RNA-seq count data** and metadata.
  - **Tissue-specific annotations** and color mappings.

- **Outputs**:
  - **UMAP visualizations** of pseudotime trajectories and tissue-specific cell states.
  - **Pseudotime plots** showing differentiation trajectories.
  - **Processed Monocle3 object** saved for further analysis. 
  
### 10.multi-species.R
- **Purpose**: This script integrates and analyzes single-cell RNA-seq data across multiple species (rice, maize, Setaria, Sorghum, and Arabidopsis) using **Seurat**. It identifies conserved and species-specific gene expression patterns, performs differential expression analysis, and visualizes results.

- **Inputs**:
  - **Single-cell RNA-seq data** for multiple species.
  - **Ortholog information** to map genes across species.

- **Outputs**:
  - **Integrated Seurat object** for cross-species analysis.
  - **Differential expression results** and functional enrichment analysis.
  - **Heatmaps** and **UMAP visualizations** of conserved and species-specific gene expression patterns. 
 
### 11.processed_LTPL120.R
- **Purpose**: This script performs single-cell RNA-seq data analysis using **Seurat** and **DoubletFinder** to process and analyze data from multiple samples (LTP and NIP). It includes quality control, doublet removal, data integration, clustering, and differential expression analysis.

- **Inputs**:
  - **Single-cell RNA-seq data** for LTP and NIP samples.
  - **Mitochondrial gene list** for quality control.

- **Outputs**:
  - **Quality control plots** and filtered Seurat objects.
  - **Integrated Seurat object** after harmony batch correction.
  - **Cluster-specific marker genes** and cell type annotations.
  - **Differential expression analysis** between LTP and NIP conditions. 

### DensityScatter.R
- **Purpose**: This script defines a function `DensityScatter` to create a density scatter plot for visualizing the relationship between two variables in a Seurat object. It includes options for log transformation and quantile-based thresholding.

- **Inputs**:
  - **Seurat object** containing metadata with the specified variables.
  - **Optional parameters** for log transformation (`log_x`, `log_y`) and quantile-based thresholding (`quantiles`).

- **Outputs**:
  - **Density scatter plot** with color-coded density values.
  - **Optional quantile lines** for thresholding.



## Citation 
For usage of the code and associated manuscript,  If you use our codes, please cite our paper [A single-cell multiomics atlas of rice](https://github.com/dongwei-2023/Single_Cell_Multiomics_in_Rice).

## Questions and errors
If you have a question, error, bug, or problem, please use the [Github issue page](https://github.com/dongwei-2023/Single_Cell_Multiomics_in_Rice/issues).

## Contact  
  - **Zhe Liang**: [liangzhe@caas.cn](mailto:liangzhe@caas.cn)  
  - **Xiaofeng Gu**: [guxiaofeng@caas.cn](mailto:guxiaofeng@caas.cn)  
  Biotechnology Research Institute, Chinese Academy of Agricultural Sciences
