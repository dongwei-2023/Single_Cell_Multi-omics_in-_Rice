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

### 01.preprocess.sh

- **Purpose**: This script processes single-cell multiomics data (RNA-seq + ATAC-seq) using **Cell Ranger ARC** to generate count matrices for multiple rice samples. It includes data preparation, quality control, and analysis of gene expression and chromatin accessibility.

- **Inputs**:
  - **Reference genome**: Path to the rice reference genome (Oryza sativa MSU7`).
  - **Sample names**: List of sample names (e.g., `leaf_1`,  etc.).
  - **Raw data paths**: Paths to the raw data directories for each sample, containing FASTQ files for RNA-seq and ATAC-seq (e.g., `/rawdata/leaf_1`).

- **Outputs**:
  - **CSV files**: Specifying input data and sample information for each sample (e.g., `leaf_1.csv`).
  - **Cell Ranger ARC output directories**:
    - **Filtered matrices**: Gene expression and chromatin accessibility matrices.
    - **Clustering and visualization results**: UMAP plots and other visualizations.
    - **Analysis summaries**: Web reports (`web_summary.html`) and metrics (`metrics_summary.csv`).

### 02.QC.R

- **Purpose**: This script processes single-cell RNA-seq data using **Seurat** and **DoubletFinder** to perform quality control, doublet removal, clustering, and visualization for multiple rice tissue samples. It includes steps for filtering low-quality cells, identifying doublets, and analyzing gene expression patterns (Figure S1).

- **Inputs**:
  - **Single-cell RNA-seq data**: Filtered feature-barcode matrices for multiple samples (e.g., `leaf_1`, `Root_1`, etc.).
  - **Doublet rate file**: A file (`double_rate.txt`) containing expected doublet rates for each sample.

- **Outputs**:
  - **Quality control plots** (violin plots, scatter plots) for cell filtering.
  - **Doublet detection results** and filtered single-cell data.
  - **Clustering results** and visualizations (UMAP, clustree plots).
  - **Processed Seurat objects** saved as `.Rdata` files for downstream analysis.

### 03.merge_rep.R

- **Purpose**: This script processes single-cell RNA-seq data for multiple rice tissue samples using **Seurat**. It performs clustering, marker gene identification, and visualization, followed by data integration and dimensionality reduction for each tissue type.

- **Inputs**:
  - **Pre-processed Seurat objects**: Loaded from saved RData files for each sample (e.g., `leaf_1/data.Rdata`).

- **Outputs**:
  - **UMAP plots**: Visualizations of cell clusters for each sample.
  - **Marker gene lists**: CSV files containing marker genes for each cluster.
  - **Dot plots**: Visualizations of top marker genes for each cluster.
  - **Feature plots**: UMAP plots showing expression of top marker genes.
  - **Integrated Seurat objects**: Merged and processed data for each tissue type, saved as `alldata.Rdata`.

### 04.wnn_integrated.R

- **Purpose**: This script processes single-cell multiomics data (RNA-seq + ATAC-seq) using **Signac** and **Seurat** to integrate gene expression and chromatin accessibility data. It performs clustering, visualization, and quality control for rice samples (Figure 1).

- **Inputs**:
  - **Filtered feature-barcode matrices**: Combined RNA-seq and ATAC-seq data from `/aggr/all/outs/filtered_feature_bc_matrix/`.
  - **ATAC fragment file**: Path to the ATAC fragment file (`/aggr/all/outs/atac_fragments.tsv.gz`).
  - **Genome annotations**: Rice genome annotations from `EnsDb.Osativa.v55`.
  - **Pre-processed cluster labels**: Loaded from `alldata.Rdata`.

- **Outputs**:
  - **Integrated Seurat object**: Combines RNA-seq and ATAC-seq data into a single Seurat object (`rice.sub`).
  - **UMAP visualizations**: For RNA-seq (`umap.rna`), ATAC-seq (`umap.atac`), and weighted nearest neighbor (WNN) integration (`wnn.umap`).
  - **Quality control metrics**: Nucleosome signal and TSS enrichment scores.
  - **Saved Seurat object**: Processed data saved as `rice.sub.Rdata`.


### 05.Harmony_integrated.R

- **Purpose**: This script integrates single-cell RNA-seq and ATAC-seq data using **Harmony** and **Seurat** to perform batch correction, clustering, and visualization. It focuses on analyzing rice tissue samples and generating UMAP visualizations for integrated data.

- **Inputs**:
  - **Pre-processed Seurat object**: Loaded from `rice.sub.Rdata`, containing RNA-seq and ATAC-seq data.
  - **Tissue and cluster annotations**: Pre-defined tissue and cluster labels.

- **Outputs**:
  - **Integrated UMAP visualizations**: For RNA-seq, ATAC-seq, and weighted nearest neighbor (WNN) integration.
  - **Cluster and tissue-specific plots**: PDF files showing UMAP visualizations for clusters and tissues.
  - **CSV files**: Containing cell embeddings, cluster IDs, and tissue annotations.
  - **Quality control metrics**: Nucleosome signal and TSS enrichment scores.


### 06.correlation.R

- **Purpose**: This script performs correlation analysis on single-cell RNA-seq and ATAC-seq data to explore relationships between gene expression and chromatin accessibility. It generates correlation heatmaps and integrates RNA-seq and ATAC-seq data for rice samples(Figure 2a).

- **Inputs**:
  - **Pre-processed Seurat object**: Loaded from `rice.sub.Rdata`, containing RNA-seq and ATAC-seq data.
  - **Variable features**: Used for calculating average expression and accessibility.

- **Outputs**:
  - **Correlation heatmaps**:
    - RNA-seq, ATAC-seq, RNA-ATAC correlation heatmap.
  - **CSV files**:
    - Correlation values between RNA-seq and ATAC-seq data (`atac_RNA.cor.csv`).


### 07.motif_TF.R

- **Purpose**: This script performs motif analysis and visualization for single-cell ATAC-seq data using **Signac** and **motifmatchr**. It identifies enriched transcription factor motifs, computes chromatin accessibility profiles, and generates coverage plots for specific genes.

- **Inputs**:
  - **Pre-processed Seurat object**: Loaded from `rice.sub.Rdata`, containing ATAC-seq data.
  - **Position frequency matrices (PFMs)**: Loaded from  rice motif pfm list for motif analysis.

- **Outputs**:
  - **Enriched motifs**: For each tissue or cluster, saved as a data frame.
  - **Coverage plots**: For specific genes (e.g., `Os03g0821200` and `Os05g0121600`) across tissues.


### 08.run_rice_create_cistarget.sh

- **Purpose**: This script creates a **cisTarget motif database** for rice using the [create_cisTarget_databases] pipeline(https://github.com/aertslab/create_cisTarget_databases). It processes genomic regions and motif data to generate a database for regulatory network analysis.

- **Inputs**:
  - **FASTA file**: Contains genomic regions (e.g., gene TSS ± 3 kb) for motif enrichment analysis (`Rice_gene_tss_updown3k.fasta`).
  - **Motif files**: Position frequency matrices (PFMs) in `cb` format use [cluster-buster](https://github.com/ghuls/cluster-buster/), located in `motifs_dir`.
  - **Motif list**: A text file listing the motifs to include in the database (`Rice_cisbp_motifs_list.txt`).

- **Outputs**:
  - **cisTarget motif database**: A set of files prefixed with `db_prefix`, containing motif enrichment scores for the provided genomic regions. These files can be used as input for **SCENIC** or other regulatory network analysis tools.


### 09.runPySCENIC.sh
- **Purpose**: This script runs the **pySCENIC** pipeline to infer gene regulatory networks (GRNs) from single-cell RNA-seq data. It performs three main steps: GRN inference, motif enrichment analysis, and regulon activity scoring(Figure S2 h-k).

- **Inputs**:
  - **Single-cell RNA-seq data**: Loom file containing gene expression data (`s1_rice_scRNA_count_cells_genes.loom`).
  - **Transcription factor list**: A list of transcription factors for rice (`Rice_plantTFDB_TF_list.txt`).
  - **Motif database**: A motif database file in TBL format (`rice_plantTFDB_motif.tbl`).
  - **cisTarget motif rankings**: Feather files containing motif rankings for rice (`Rice_gene_tss*.regions_vs_motifs.rankings.feather`).

- **Outputs**:
  - **Gene regulatory network (GRN)**: Adjacency list of regulatory interactions (`s2_rice_scRNA_count_cells_genes.adj.tsv`).
  - **Regulons**: List of enriched motifs and their associated regulons (`s2_rice_scRNA_count_cells_genes.reg.tsv`).
  - **Regulon activity scores**: Loom file containing AUCell scores for regulons (`s2_rice_scRNA_count_cells_genes.pyscenic.loom`).



### 10.RSR1_gene_KO_simulation_with_rice_root.ipynb
- **Purpose**: This script performs **in silico gene perturbation analysis** using [CellOracle](https://github.com/morris-lab/CellOracle) to simulate the effects of gene knockout (e.g., `Os05g0121600`) on cell identity and gene regulatory networks. It visualizes predicted cell state transitions and vector fields(Figure 2f-h).

- **Inputs**:
  - **Pre-processed CellOracle object** (`rice_root.celloracle.oracle`) containing single-cell RNA-seq data and GRN information.
  - **Gene regulatory network links** (`rice_root_links.celloracle.links`).

- **Outputs**:
  - **Quiver plots** and **vector field graphs** showing predicted cell state transitions after gene perturbation.
  - **Visualizations** of gene expression changes and cell identity shifts. 
 
### 11.leaf_Flag.R
- **Purpose**: This script performs single-cell RNA-seq data analysis using **Seurat** and **Monocle** to study cell type-specific gene expression, pseudotime trajectories, and differential expression between tissues (e.g., leaf and flag). It also integrates data from multiple samples and performs functional enrichment analysis.

- **Inputs**:
  - **Single-cell RNA-seq data** for leaf and flag tissues.
  - **Pre-processed Seurat object** (`alldata.Rdata`) containing integrated data.

- **Outputs**:
  - **UMAP visualizations** for cell clusters and tissue-specific expression.
  - **Differential expression analysis** results and functional enrichment plots.
  - **Pseudotime trajectory analysis** using Monocle.
  - **Integrated Seurat object** and processed data saved for downstream analysis. 
 
### 12.hdWGCNA.R
- **Purpose**: This script performs **Weighted Gene Co-expression Network Analysis (WGCNA)** using the **hdWGCNA** package to identify co-expression modules and analyze their functional relevance in single-cell RNA-seq data.

- **Inputs**:
  - **Pre-processed Seurat object** (`alldata.Rdata`) containing single-cell RNA-seq data.
  - **Gene expression data** for constructing co-expression networks.

- **Outputs**:
  - **Co-expression modules** and their associated genes.
  - **Module eigengene visualizations** and functional enrichment analysis.
  - **Heatmaps** showing module correlations and gene expression patterns. 
 
### 13.Epidermis_Vascular_monocle3.R
- **Purpose**: This script performs **pseudotime trajectory analysis** using **Monocle3** to study cell differentiation and developmental trajectories in single-cell RNA-seq data.

- **Inputs**:
  - **Single-cell RNA-seq count data** and metadata.
  - **Tissue-specific annotations** and color mappings.

- **Outputs**:
  - **UMAP visualizations** of pseudotime trajectories and tissue-specific cell states.
  - **Pseudotime plots** showing differentiation trajectories.
  - **Processed Monocle3 object** saved for further analysis. 
  
### 14.multi-species.R
- **Purpose**: This script integrates and analyzes single-cell RNA-seq data across multiple species (rice, maize, Setaria, Sorghum, and Arabidopsis) using **Seurat**. It identifies conserved and species-specific gene expression patterns, performs differential expression analysis, and visualizes results.

- **Inputs**:
  - **Single-cell RNA-seq data** for multiple species.
  - **Ortholog information** to map genes across species.

- **Outputs**:
  - **Integrated Seurat object** for cross-species analysis.
  - **Differential expression results** and functional enrichment analysis.
  - **Heatmaps** and **UMAP visualizations** of conserved and species-specific gene expression patterns. 
 
### 15.processed_LTPL120.R
- **Purpose**: This script performs single-cell RNA-seq data analysis using **Seurat** and **DoubletFinder** to process and analyze data from multiple samples (LTP and NIP). It includes quality control, doublet removal, data integration, clustering, and differential expression analysis.

- **Inputs**:
  - **Single-cell RNA-seq data** for LTP and NIP samples.
  - **Mitochondrial gene list** for quality control.

- **Outputs**:
  - **Quality control plots** and filtered Seurat objects.
  - **Integrated Seurat object** after harmony batch correction.
  - **Cluster-specific marker genes** and cell type annotations.
  - **Differential expression analysis** between LTP and NIP conditions. 

### 16.DensityScatter.R
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
