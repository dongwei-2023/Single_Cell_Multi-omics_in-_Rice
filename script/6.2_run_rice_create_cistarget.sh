conda activate create_cistarget_databases
create_cistarget_databases_dir="/data2/lidongwei/work/work_scMultiome/part_scenic/Github/create_cisTarget_databases-master"
motifs_dir="/data2/lidongwei/work/work_scMultiome/part_scenic/cisTarget_databases/rice/planTFDB_motif/input/motif_data/motifs_cb_format"
motifs_list_filename="/data2/lidongwei/work/work_scMultiome/part_scenic/cisTarget_databases/rice/planTFDB_motif/input/motif_data/Rice_cisbp_motifs_list.txt"
nbr_threads=30

fasta_filename="/data2/lidongwei/work/work_scMultiome/part_scenic/cisTarget_databases/rice/planTFDB_motif/input/fasta/Rice_gene_tss_updown3k.fasta"
db_prefix="/data2/lidongwei/work/work_scMultiome/part_scenic/cisTarget_databases/rice/planTFDB_motif/output/Rice_gene_tss_updown3k"

"${create_cistarget_databases_dir}/create_cistarget_motif_databases.py" \
    -f "${fasta_filename}" \
    -M "${motifs_dir}" \
    -m "${motifs_list_filename}" \
    -o "${db_prefix}" \
    -t "${nbr_threads}"
