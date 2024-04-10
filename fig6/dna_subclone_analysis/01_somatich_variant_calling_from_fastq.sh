#!/bin/bash

run_thread=60
working_dir=WORKING_DIR
export data=$working_dir/data
export raw_fa_dir=$data/raw_files
export bwa_reference=$data/DIR_TO/GRCm38.primary_assembly.genome.fa
this_script=$working_dir/script/somatich_variant_calling_from_fastq.sh
export results=$working_dir/results
export PICARD=$working_dir/softwares/picard.jar
export TMP_DIR=$working_dir/tmp
export recal_out=$results/recalibrated_bam
export mutect_normal_samples_for_PoN=$results/mutect_normal_samples_for_PoN
export pon_db=$results/pon_db
export pon_vcf=$results/pon_vcf
export mutect_germline_db=$results/mutect_germline_db
export mutect_out_tumor_only=$results/mutect_tumor_only_mode
export mutect_out_pairs=$results/mutect_pairs_tumor_normal
export mutect_out_germline=$results/mutect_germline
export cnv_called=$results/cnv_called
export interval=$data/contig_mm10/mm10_contigs.interval.list
export CNV_FACET=$working_dir/softwares/cnv_facets/0.16.0/cnv_facets.R
export contig_regions=$data/contig_mm10/mm10_contigs.interval.list
export blacklist_regions=$data/blacklist_mm10/mm10-blacklist.interval.list

mkdir -p $results $results/bwa_mapped_bam $recal_out $mutect_out_tumor_only $TMP_DIR \
    $mutect_normal_samples_for_PoN $pon_db $pon_vcf $mutect_out_tumor_only \
    $mutect_out_pairs $results/qc/fastqc $results/qc/bwa_qc/ $mutect_out_germline \
    $mutect_germline_db "$mutect_out_germline"/combined $cnv_called $results/cnv_varscan_called \
    $phylowgs_out

cat $this_script

###################################### 
### 2. Load modules 
######################################

#
module load gcc/11.2.0
module load R/4.2.2 
module load R-packages
module load openjdk/20.0.0
module load bwa
module load samtools
module load gatk
module load perl
module load htslib
module load fastqc
module load multiqc
module load cnv_facets
module load bedtools 


###################################### 
### 3. Build bwa reference and mouse variation db
######################################

# {wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz && \
# gzip -dk GRCm38.primary_assembly.genome.fa.gz && \
# samtools faidx GRCm38.primary_assembly.genome.fa
# nohup bwa index GRCm38.primary_assembly.genome.fa } &


## create fasta dict file for GATK
# java -jar /home/gkc251/project_folder/202305_organoids_dna_seq_somatic/softwares/picard.jar CreateSequenceDictionary \
#       -R GRCm38.primary_assembly.genome.fa \
#       -O GRCm38.primary_assembly.genome.dict

## download snp all dbsnp142 from genome browser
## mgpV5MergedSNPsAlldbSNP142.vcf.gz
## get only the first columns
## zcat mgpV5MergedSNPsAlldbSNP142.vcf.gz | awk '{ printf "%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5; } ' > mgpV5MergedSNPsAlldbSNP142.tsv &

## download mm10 contigs and convert to interval list
# awk '{ printf "%s:%s-%s\n",$1,$2,$3; } ' mm10_contigs.bed | awk '$1~/^chr([1-9]|1[0-9]|X|Y):[0-9]*-[0-9]*$/' > mm10_contigs.interval.list

## download ENCODE blacklist regions and convert to interval list
# awk '{ printf "%s:%s-%s\n",$1,$2,$3; } ' mm10-blacklist.v2.0to1.bed | awk '$1~/^chr([1-9]|1[0-9]|X|Y):[0-9]*-[0-9]*$/' > mm10-blacklist.interval.list


###################################### 
### 3.(2) Quality check 
######################################


quality_check_fastq () {
    fastqc $raw_fa_dir/*/*/*fq.gz -t 32 -o $results/qc/fastqc
}

export -f quality_check_fastq

###################################### 
### 4. bwa alignment + sort
### tutorial: https://hbctraining.github.io/In-depth-NGS-Data-Analysis-Course/sessionVI/lessons/01_alignment.html
### https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/
######################################



pipeline_bwa_alignment () {
    local sample=$1
    echo "$(date +"%Y_%m_%d %H:%M") ..... bwa align $sample" 
    r1=$(find $raw_fa_dir/$sample/ -type f -maxdepth 5 -name *_1.fq.gz) &&
    r2=$(find $raw_fa_dir/$sample/ -type f -maxdepth 5 -name *_2.fq.gz) &&

    bwa mem -t 8 -T 0 \
        -R $(echo "@RG\tID:$sample\tSM:$sample\tLB:$sample\tPL:ILLUMINA") \
        $bwa_reference \
        $r1 \
        $r2 |
        samtools sort \
        -@8 \
        -o $results/bwa_mapped_bam/$sample.bam - &&

    echo "$(date +"%Y_%m_%d %H:%M") ..... bwa align finish $sample" &&

    samtools stats -@ 4 $results/bwa_mapped_bam/$sample.bam > $results/qc/bwa_qc/$sample.stats &&

    ## calculate average depth
    samtools depth $results/bwa_mapped_bam/$sample.bam  |  \
        awk '{sum+=$3} END { print "Average = ",sum/NR}' > \
        $results/qc/bwa_qc/$sample.depth 
}

export -f pipeline_bwa_alignment

###################################### 
### 5. preprocess the bam files
### tutorial: https://github.com/broadinstitute/gatk-docs/blob/master/gatk3-methods-and-algorithms/Base_Quality_Score_Recalibration_(BQSR).md
######################################


pipeline_recalibration () {
   local sample=$1
   echo "=========processing $sample ============="

   echo "$(date +"%Y_%m_%d %H:%M") ..... mark duplicates for $sample" &&
   # 4.2. Mark duplicates
   java -Xmx8g -jar $PICARD MarkDuplicates \
        -CREATE_INDEX true \
        -VALIDATION_STRINGENCY SILENT \
        -TMP_DIR $TMP_DIR \
        -I $results/bwa_mapped_bam/$sample.bam \
        -O "$recal_out"/"$sample"_bwa.marked_dup.bam \
        -M "$recal_out"/"$sample"_marked_dup_metrics.txt &&

   ## 4.4. BaseRecalibrator
   gatk BaseRecalibrator \
      --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms4g -Xmx8g" \
      -tmp-dir $TMP_DIR \
      --use-original-qualities \
      -I "$recal_out"/"$sample"_bwa.marked_dup.bam \
      -R $bwa_reference \
      --known-sites "$data/variation_mouse_db/mgpV5MergedSNPsAlldbSNP142.vcf.gz" \
      -O "$recal_out"/"$sample"_recal_data.table1 &&

   ## 4.5. ApplyBQSR
   gatk ApplyBQSR \
      --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms3g -Xmx8g" \
      -tmp-dir $TMP_DIR \
      --use-original-qualities \
      --add-output-sam-program-record \
      -R "$bwa_reference" \
      -I "$recal_out"/"$sample"_bwa.marked_dup.bam \
      --bqsr-recal-file "$recal_out"/"$sample"_recal_data.table1 \
      -O "$recal_out"/"$sample"_bwa.marked_dup.recalibrated.bam &&

   echo "$(date +"%Y_%m_%d %H:%M") ..... baseRecalibrator 2nd round for $sample" &&
   ## 4.6 BaseRecalibrator 2nd
   gatk BaseRecalibrator \
      --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms4g -Xmx8g" \
      -tmp-dir $TMP_DIR \
      --use-original-qualities \
      -I "$recal_out"/"$sample"_bwa.marked_dup.recalibrated.bam \
      -R $bwa_reference \
      --known-sites "$data/variation_mouse_db/mgpV5MergedSNPsAlldbSNP142.vcf.gz" \
      -O "$recal_out"/"$sample"_recal_data.table2 &&

   gatk AnalyzeCovariates \
      -before "$recal_out"/"$sample"_recal_data.table1 \
      -after "$recal_out"/"$sample"_recal_data.table2 \
      -plots "$recal_out"/"$sample"_AnalyzeCovariates.pdf &&

   echo "$(date +"%Y_%m_%d %H:%M") ..... recalibration finish $sample" 

}

export -f pipeline_recalibration

########################################
### 6. Mutect calling 
### (1) Tumor-only mode to create panel of normal
### (2) Create Panel of Normal
### (3) Tumor-Normal pair mode to call somatic mutations
### (4) Germline variants calling 
### https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2
########################################

########################################
### 6.(1) Tumor-only mode to create panel of normal
### need max-mnp-distance 0 : https://gatk.broadinstitute.org/hc/en-us/community/posts/360071895952-Mutect2-s-support-for-MNP-in-GATK-4-1-8-1-
########################################

pipeline_mutect_normal_samples_for_PoN () {
    # Tumor-only mode
   local sample=$1
   echo "$(date +"%Y_%m_%d %H:%M") ..... pipeline_mutect_normal_samples_for_PoN $sample" &&
   gatk --java-options "-Xms6g -Xmx8g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
      Mutect2 \
      -tmp-dir $TMP_DIR \
      -R $bwa_reference \
      -L $interval \
      --max-mnp-distance 0 \
      -I "$recal_out"/"$sample"_bwa.marked_dup.recalibrated.bam \
      -O "$mutect_normal_samples_for_PoN"/"$sample".vcf \
      -dont-use-soft-clipped-bases true &&

   echo "$(date +"%Y_%m_%d %H:%M") ..... pipeline_mutect_normal_samples_for_PoN finish $sample" 
}

########################################
### 6.(2) Create Panel of Normal (to remove technical noise)
### tutorial: https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2
########################################

pipeline_create_panel_of_normal () {

   gatk --java-options "-Xms6g -Xmx8g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
    GenomicsDBImport \
        -tmp-dir $TMP_DIR \
        -R $bwa_reference \
        -L $interval \
        --genomicsdb-workspace-path $pon_db/pon_normal_w1_samples \
        -V "$mutect_normal_samples_for_PoN"/D2301171049.vcf \
        -V "$mutect_normal_samples_for_PoN"/D2301171051.vcf \
        -V "$mutect_normal_samples_for_PoN"/D2301171053.vcf &&

   gatk --java-options "-Xms6g -Xmx8g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
    CreateSomaticPanelOfNormals \
      -tmp-dir $TMP_DIR \
      -R $bwa_reference \
      -V gendb://$pon_db/pon_normal_w1_samples \
      -O $pon_vcf/pon.vcf.gz &&

   echo "$(date +"%Y_%m_%d %H:%M") ..... pipeline_create_panel_of_normal finish" 
}

########################################
### 6.(3) Tumor-normal pair mode to call somatic mutations 
########################################


pipeline_mutect_tumor_normal_pairs () {
   local tumor=$1
   local normal=$2
   echo "=========pipeline_mutect_tumor_normal_pairs $tumor vs $normal ============="
   gatk --java-options "-Xms6g -Xmx8g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
        Mutect2 \
        -tmp-dir $TMP_DIR \
        -R $bwa_reference \
        -L $interval \
        -pon $pon_vcf/pon.vcf.gz \
        -I "$recal_out"/"$tumor"_bwa.marked_dup.recalibrated.bam \
        -I "$recal_out"/"$normal"_bwa.marked_dup.recalibrated.bam \
        -normal $normal \
        -O "$mutect_out_pairs"/"$tumor".vcf \
        -dont-use-soft-clipped-bases true &&
   
   gatk --java-options "-Xms6g -Xmx8g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
   	FilterMutectCalls \
        -tmp-dir $TMP_DIR \
        -R $bwa_reference \
        -L $interval \
        -V "$mutect_out_pairs"/"$tumor".vcf \
        -O "$mutect_out_pairs"/"$tumor".filtered.vcf &&

    gatk --java-options "-Xms6g -Xmx8g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
   	FilterMutectCalls \
        -tmp-dir $TMP_DIR \
        -R $bwa_reference \
        -L $interval \
        --exclude-intervals $blacklist_regions \
        -V "$mutect_out_pairs"/"$tumor".vcf \
        -O "$mutect_out_pairs"/"$tumor".filtered.blacklist_exclude.vcf &&

    helper_exract_PASS_rows_from_vcf \
        "$mutect_out_pairs"/"$tumor".filtered.vcf &&

    helper_exract_PASS_rows_from_vcf \
        "$mutect_out_pairs"/"$tumor".filtered.blacklist_exclude.vcf &&

   echo "$(date +"%Y_%m_%d %H:%M") ..... pipeline_mutect_tumor_normal_pairs finish $tumor vs $normal" 

}


########################################
### 6.(4) Germline calling
### Haplotype caller + GenotypeGVCFs
### https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller
### https://gatk.broadinstitute.org/hc/en-us/articles/360037057852-GenotypeGVCFs
########################################


pipeline_mutect_haplotypecall () {
    local sample=$1
    echo "$(date +"%Y_%m_%d %H:%M") ..... pipeline_haplotypecall $1" 

    gatk --java-options "-Xms6g -Xmx8g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
        HaplotypeCaller  \
        -tmp-dir $TMP_DIR \
        -R $bwa_reference \
        -L $interval \
        -I "$recal_out"/"$sample"_bwa.marked_dup.recalibrated.bam \
        -O "$mutect_out_germline"/"$sample".g.vcf.gz \
        -dont-use-soft-clipped-bases true \
        -ERC GVCF &&

    echo "$(date +"%Y_%m_%d %H:%M") ..... pipeline_haplotypecall finished $1" 
}

pipeline_mutect_germline_call () {

   echo "$(date +"%Y_%m_%d %H:%M") ..... pipeline_mutect_germline_call" 

    gatk --java-options "-Xms6g -Xmx8g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
        CombineGVCFs \
        -tmp-dir $TMP_DIR \
        -R $bwa_reference \
        -L $interval \
        -V "$mutect_out_germline"/D2301171049.g.vcf.gz \
        -V "$mutect_out_germline"/D2301171051.g.vcf.gz \
        -V "$mutect_out_germline"/D2301171053.g.vcf.gz \
        -O "$mutect_out_germline"/combined/germline_combined.g.vcf.gz &&

    gatk --java-options "-Xms6g -Xmx8g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
        GenotypeGVCFs \
        -tmp-dir $TMP_DIR \
        -R $bwa_reference \
        -L $interval \
        -V "$mutect_out_germline"/combined/germline_combined.g.vcf.gz \
        -O "$mutect_out_germline"/combined/germline_w1.vcf.gz &&


    # https://gatk.broadinstitute.org/hc/en-us/articles/360037055952-SelectVariants
    # https://gatk.broadinstitute.org/hc/en-us/articles/360035891011-JEXL-filtering-expressions

    gatk --java-options "-Xms6g -Xmx8g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
        SelectVariants \
        -tmp-dir $TMP_DIR \
        -R $bwa_reference \
        -V "$mutect_out_germline"/combined/germline_w1.vcf.gz \
        --select-type-to-include SNP \
        -L $interval \
        --exclude-intervals $blacklist_regions \
        --intervals $contig_regions \
        -select "MQ>=40.0 && QD>=2.0 && FS<=60.0 && SOR<=3.0 && MQRankSum>=-12.5 && ReadPosRankSum>=-8.0" \
        -O "$mutect_out_germline"/combined/germline_w1.snp.hardfiltered.blacklist_exclude.vcf.gz && 

    
    echo "$(date +"%Y_%m_%d %H:%M") ..... pipeline_mutect_germline_call finished" 

}


##########################################
### 7 Helper functions
##########################################

pipeline_snpeff_annotate () {
   local vcffile=$1
   local samplename=$2
   echo "=========processing $vcffile ============="

    ## 6.1 SnpEff 
    java -Xmx8g -jar $SNPEFF GRCm38.99 \
        $vcffile \
        > "$vcffile".snpeff.vcf &&

    ## 6.2 Convert to maf
    perl $VCF2MAF \
        --input-vcf "$vcffile".snpeff.vcf \
        --output-maf "$vcffile".snpeff.maf \
        --tmp-dir $TMP_DIR \
        --tumor-id $samplename \
        --inhibit-vep \
        --ref-fasta $bwa_reference \
        --species mus_musculus \
        --ncbi-build GRCm38 &&

   echo "$(date +"%Y_%m_%d %H:%M") ..... finish snpeff $sample" 

}

helper_exract_PASS_rows_from_vcf () {
    local vcffile=$1
    awk -F '\t' '{if($0 ~ /\#/) print; else if($7 == "PASS") print}' $vcffile > $vcffile.PASS.vcf
}

export -f pipeline_mutect_normal_samples_for_PoN \
        pipeline_create_panel_of_normal \
        pipeline_mutect_tumor_normal_pairs \
        pipeline_snpeff_annotate \
        helper_exract_PASS_rows_from_vcf \
        pipeline_mutect_haplotypecall \
        pipeline_mutect_germline_call


########################################
### 8. CNA calling 
### https://github.com/dariober/cnv_facets
########################################


pipeline_cna_calling () {
    local tumor=$1
    local normal=$2
    local het_snps=$3
    local out_prefix=$4
    mkdir -p $cnv_called/"$out_prefix"_WGS
    tumor_bam="$recal_out"/"$tumor"_bwa.marked_dup.recalibrated.bam
    normal_bam="$recal_out"/"$normal"_bwa.marked_dup.recalibrated.bam
    pileup_bam=$cnv_called/"$out_prefix"_WGS/"$tumor"vs"$normal"_$out_prefix.csv.gz
    min_depth=10
    max_depth=100
    min_cval=25
    max_cval=400
    nbhd_snp=500
    snp_mapq=40
    param_file=$cnv_called/"$out_prefix"_WGS/"$tumor"vs"$normal"_$out_prefix.param

    echo "Parameter of cnv_facet" > $param_file
    echo "depth $min_depth $max_depth" >> $param_file
    echo "cval $min_cval $max_cval" >> $param_file
    echo "nbhd_snp $nbhd_snp" >> $param_file
    echo "snp_mapq $snp_mapq" >> $param_file
    echo "snp sites $het_snps" >> $param_file

    ## parameter for WGS\
    ## if no pile up file is available 
    $CNV_FACET \
        --snp-nprocs 2 \
        --gbuild mm10 \
        --depth $min_depth $max_depth \
        --cval $min_cval $max_cval \
        --snp-mapq $snp_mapq \
        --nbhd-snp $nbhd_snp \
        -t $tumor_bam \
        -n $normal_bam \
        -vcf $het_snps \
        -o $cnv_called/"$out_prefix"_WGS/"$tumor"vs"$normal"_$out_prefix &&

    # if pile up file is available
    # $CNV_FACET \
    #     --gbuild mm10 \
    #     --depth 10 100 \
    #     --cval 25 400 \
    #     --snp-mapq 40 \
    #     --nbhd-snp 500 \
    #     --rnd-seed 12345 \
    #     -p $pileup_bam \
    #     -o $cnv_called/"$out_prefix"_WGS/"$tumor"vs"$normal"_$out_prefix &&

    echo "$(date +"%Y_%m_%d %H:%M") ..... finish pipeline_cna_calling $tumor vs $normal" 

}

pipeline_cna_calling_snp () {
    local tumor=$1
    local normal=$2

    pipeline_cna_calling \
        $tumor \
        $normal \
        $data/variation_mouse_db/mgpV5MergedSNPsAlldbSNP142.vcf.gz \
        snp_20230608_snp142
}


export -f pipeline_cna_calling pipeline_cna_calling_snp pipeline_cna_calling_pseudo_for_sanity


########################################
### To run this script parallely
########################################

parallel --jobs 80% --memfree 100g -N 1 --lb --tag pipeline_bwa_alignment :::: all_samples.csv  &&

parallel --jobs 80% --memfree 100g -N 1 --lb --tag pipeline_recalibration :::: all_samples.csv &&

parallel --jobs 80% --memfree 100g -N 1 --lb --tag pipeline_mutect_normal_samples_for_PoN :::: normal_w1_samples.csv &&

pipeline_create_panel_of_normal &&

parallel --jobs 80% --memfree 100g -N 2 --lb --tag pipeline_mutect_tumor_normal_pairs :::: paired_tumor_normal_samples.csv 

parallel --jobs 80% --memfree 100g -N 2 --lb --tag pipeline_cna_calling_snp :::: paired_tumor_normal_samples.csv 

wait

echo "===============ALL FINISHED $(date +"%Y_%m_%d_%H_%M")================" 

exit 0

