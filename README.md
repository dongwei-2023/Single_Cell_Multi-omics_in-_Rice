# Single Cell Multiomics in Rice analysis code

**Accompanying code repository for the paper entitled: "A single-cell multiomics atlas of rice".**

## **Abstract**
Cell functions across eucaryotes are driven by specific gene expression programs that are dependent on chromatin structure. Here, we report a single-cell multiomics atlas of rice, one of the major crop plants. We simultaneously profiled chromatin accessibility and RNA expression in 116,564 cells from eight organs. We annotated the majority of cell types and identified cell-type expression programs, and found a high correlation between RNA and chromatin accessibility across cells. We constructed single-cell gene regulatory networks and co-expression networks, which allowed us to identify the cell type-specific regulators RSR1, F3H and LTPL120 in rice development. Furthermore, our analysis revealed a correlation between cell type and agronomic traits of rice, and the conserved and divergent functions of cell type during evolution. In summary, this study not only offers a valuable single-cell multiomics resource, but also enriches our understanding of the intricate roles and molecular underpinnings of individual cell types in rice.

## **Schematic overview of the workflow**
![image](https://github.com/dongwei-2023/Single_Cell_Multiomics_in_Rice/blob/master/img/pipeline.png)

## Raw and processed data 
* The raw and processed data of single-cell multiomics, scRNA-seq and RNA-seq have been deposited into NCBI GEO with accession numbers: [GSE232863](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE232863), [GSE273875](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE273875), [GSE245410](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245410). Processed single-cell multiomics seurat objects are available from [GSE232863_rice_scRNA_scATAC_all.Rdata.gz](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE232863&format=file&file=GSE232863%5Frice%5FscRNA%5FscATAC%5Fall%2ERdata%2Egz) and processed LTPL120 seurat objects are available from [GSE273875_LTP_NIP_label_240622.rds.gz](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE273875&format=file&file=GSE273875%5FLTP%5FNIP%5Flabel%5F240622%2Erds%2Egz).

* Single Cell Multi-omics in Rice chromatin accessibility and RNA expression zonation profiles can be browsed via our web-app: https://www.elabcaas.cn/scmr/

* Processed zonation tables are provided as supplementary tables in the manuscript. 


## Citation 
For usage of the code and associated manuscript,  If you use our codes, please cite our paper [A single-cell multiomics atlas of rice](https://github.com/dongwei-2023/Single_Cell_Multiomics_in_Rice).

## Questions and errors
If you have a question, error, bug, or problem, please use the [Github issue page](https://github.com/dongwei-2023/Single_Cell_Multiomics_in_Rice/issues).

## Notes

Contact: Zhe Liang & Xiaofeng Gu at Biotechnology Research Institute, Chinese Academy of Agricultural Sciences
