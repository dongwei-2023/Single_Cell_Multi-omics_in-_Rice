## inputs
f_loom_path_scenic=/data2/lidongwei/work/work_scMultiome/part_scenic/allTissue/s1_rice_scRNA_count_cells_genes.loom

## outputs
grn_output=/data2/lidongwei/work/work_scMultiome/part_scenic/allTissue/results/s2_rice_scRNA_count_cells_genes.adj.tsv
ctx_output=/data2/lidongwei/work/work_scMultiome/part_scenic/allTissue/results/s2_rice_scRNA_count_cells_genes.reg.tsv
f_pyscenic_output=/data2/lidongwei/work/work_scMultiome/part_scenic/allTissue/results/s2_rice_scRNA_count_cells_genes.pyscenic.loom

## reference
f_tfs=/data2/lidongwei/work/work_scMultiome/part_scenic/allTissue/planTFDB_motif/input/motif_data/Rice_plantTFDB_TF_list.txt
f_motif_path=/data2/lidongwei/work/work_scMultiome/part_scenic/allTissue/planTFDB_motif/input/motif_data/rice_plantTFDB_motif.tbl
f_db_names=`find /data2/lidongwei/work/work_scMultiome/part_scenic/cisTarget_databases/rice/planTFDB_motif/output/ -name "Rice_gene_tss*.regions_vs_motifs.rankings.feather"`


echo "################# runing  grn ###################### "
pyscenic grn \
    --num_workers 30 \
    --output $grn_output \
    --method grnboost2 \
    $f_loom_path_scenic \
    $f_tfs


echo "################# runing  ctx ###################### "
pyscenic ctx \
    $grn_output \
    $f_db_names \
    --annotations_fname $f_motif_path \
    --expression_mtx_fname $f_loom_path_scenic \
    --output $ctx_output \
    --num_workers 30 \
    --no_pruning \
    --mode "dask_multiprocessing"


echo "################# runing  aucell ###################### "
pyscenic aucell \
    $f_loom_path_scenic \
    $ctx_output \
    --output $f_pyscenic_output \
    --num_workers 30 \
    --seed 777
