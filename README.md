### QC/mapping scripts
generate_index.sh ## is used to generate the salmon index
run_QC.sh ## runs quality control and trimming programs on the fastq files
map_reads.sh ## maps the reads to the mouse transcriptome

### gene expression analysis
all_gene_expression_time_series_BBvsW1.r ## generates the DEG table deg_voomWQW.tsv used by all the fig2 scripts. Also generates figS2b. 


### scripts used to generate supplementary figures
PCA_supfig.r ## figS2a
upset_plots_supfig.r ## figS2d
plot_ABC_gene_expr.r ## was used to answer a comment from a referee 

### scripts used to generate figure 2
gsea_fig2.r ## fig2d
heatmap_fig2.r ## fig2b
upset_plots_fig2.r ## fig2c

##### Informations
The necessary data to run the scripts are in the folder "data". all the results will be output in the directory "ouput"
The scripts from the section "QC/mapping scripts" cannot be run without the data available in the SRA accession. Some intermediary files
to run the other scripts are therefore provided to the form of expression matrices
