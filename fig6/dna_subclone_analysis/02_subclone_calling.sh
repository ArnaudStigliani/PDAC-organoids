#!/bin/bash

run_thread=60
export working_dir=WORKING_DIR
export data=$working_dir/data
this_script=$working_dir/script/02_subclone_calling.sh
export results=$working_dir/results
export TMP_DIR=$working_dir/tmp
export phylowgs_folder=PATH_TO/phylowgs/2205be1

cat $this_script

########################################
### 1. PhyloWGS
## https://github.com/morrislab/phylowgs
## https://github.com/morrislab/phylowgs/tree/master/parser
## before everything, run code in 
## parse your cnv data to cnv_data.txt 
########################################

pipeline_phyloWGS () {

    module load phylowgs
    local tumor=$1
    local normal=$2
    local input_file=$3
    local output_dir=$4
    local out_prfix=$output_dir/"$tumor"vs"$normal"

    mkdir -p $output_dir $output_dir/results $output_dir/results/"$tumor"vs"$normal" $output_dir/reports/"$tumor"vs"$normal"

    python $phylowgs_folder/parser/create_phylowgs_inputs.py \
        --vcf-type sample1=mutect_smchet \
        --regions all \
        --tumor-sample $tumor \
        --cnvs sample1=PATH_TO/"$tumor"vs"$normal".cnv.txt \
        --output-cnvs $out_prfix.cnv_data.txt \
        --output-variants $out_prfix.ssm_data.txt \
        --output-params $out_prfix.param.json \
        --nonsubsampled-variants $out_prfix.nonsubsampled.variants \
        --nonsubsampled-variants-cnvs $out_prfix.nonsubsampled.cnvs \
        sample1=$input_file &&
    

    python $phylowgs_folder/multievolve.py \
        --num-chains 4 \
        --ssms $out_prfix.ssm_data.txt \
        --cnvs $out_prfix.cnv_data.txt \
        --output-dir $out_prfix/chains &&

    python $phylowgs_folder/write_results.py \
        --include-ssm-names \
        --include-multiprimary \
        "$tumor"vs"$normal" \
        $out_prfix/chains/trees.zip \
        $output_dir/results/"$tumor"vs"$normal"/"$tumor"vs"$normal".summ.json.gz \
        $output_dir/results/"$tumor"vs"$normal"/"$tumor"vs"$normal".muts.json.gz \
        $output_dir/results/"$tumor"vs"$normal"/"$tumor"vs"$normal".mutass.zip &&
   

    echo "$(date +"%Y_%m_%d %H:%M") ..... finish pipeline_phyloWGS" 
}


pipeline_phyloWGS_blacklist_exclude () {
    local tumor=$1
    local normal=$2
    local input_file=$3

    pipeline_phyloWGS \
        $tumor \
        $normal \
        $results/mutect_pairs_tumor_normal/$tumor.filtered.blacklist_exclude.vcf.PASS.vcf \
        $results/subclone_phylowgs_exclude_blacklist

}

########################################
### To run this script
########################################

export -f pipeline_phyloWGS pipeline_phyloWGS_blacklist_exclude 

(   
    parallel --dry-run -N 2 pipeline_phyloWGS_blacklist_exclude :::: paired_tumor_normal_samples.csv 
) | parallel --jobs 80% --memfree 100g --tag --lb &&

echo "===============ALL FINISHED $(date +"%Y_%m_%d_%H_%M")================" 

exit 0
