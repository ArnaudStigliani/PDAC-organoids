### QC/mapping scripts
*Located in `mapping/`*  
`generate_index.sh` is used to generate the salmon index  
`run_QC.sh`  runs quality control and trimming programs on the fastq files  
`map_reads.sh` maps the reads to the mouse transcriptome  

### Informations
The necessary data to run the scripts are in the folder `data/`. all the results will be output in the directory `fig*/results`  
The scripts from the section "QC/mapping scripts" cannot be run without the data available in the SRA accession. Some intermediary files 
to run the other scripts are therefore provided to the form of expression matrices.  
In addition, the scripts needed to determine the organoid surface from the brightfield pictures have been adeed in the folder 
`orgaquant_organoids_detection/`.  


### gene expression analysis
`figS2/limmaDEG_figS2C.r` generates the DEG table deg_voomWQW.tsv used by all the `fig2/*` scripts.  
Also generates `figS2/figS2C`.  

### scripts used to generate supplementary figures
`figS2/PCA_figS2a.r` figS2a  
`figS2/upset_S2E.r` figS2d  

`figS3/figS3_boxplots_barplots_figS3C.r` boxplots figS3b, barplots figS3f,g, part of dotplots figS3c,d.  

`figS4/heatmap_figS4df.r` heatmap plot figS4 df

`figS5/plot_ABC_gene_expr_figS5e.r` figS5e, was used to answer a comment from a referee  

### scripts used to generate main figures
`fig1/fig1_boxplots_anova.r` fig1e  

`fig2/gsea_dotplot_fig3D.r` fig2d  
`fig2/heatmap_fig2B.r` fig2b  
`fig2/barplots_fig2C.r` fig2c  

`fig3/fig3_boxplots_barplots_figS3C.r` boxplots fig3b, barplots fig3c,d, part of dotplots figS3c,d.  

`fig4/permutation_fig4bc.r` permutation plot fig4 bc  

`fig6/03_subclone_analysis.qmd` contains the code to generate figure 6  
