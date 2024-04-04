#!/bin/bash

data=../data
shared_data1=/isdata/sandelin/projects/pH_Pancreatic_Organoids/pdac_organiods_1_54/F20FTSEUHT0248_MUSasdE/Clean
shared_data2=/isdata/sandelin/projects/pH_Pancreatic_Organoids/pdac_organoids_55_90/F20FTSEUHT0248_MUSyogE/Clean
shared_data3=/isdata/sandelin/projects/pH_Pancreatic_Organoids/pdac_organoids_91_108/F20FTSEUHT0248-02_MUSbmtE/Clean
shared_data4=/isdata/sandelin/projects/pH_Pancreatic_Organoids/pdac_organoids_mN_andp53R273H/raw_data/F21FTSEUHT1004_MOUrelbR/Filter_SOAPnuke/Clean
shared_data5=/isdata/sandelin/projects/pH_Pancreatic_Organoids/pdac_organoids_drug_adaptation/F22FTSEUET0027_MUSydrtR/Clean

results=../results/run_QC

results1=$results/results_batch_1
results2=$results/results_batch_2
results3=$results/results_batch_3
results4=$results/results_batch_4
results5=$results/results_batch_5
multiqc1=$results/multiqc_1
multiqc2=$results/multiqc_2
multiqc3=$results/multiqc_3
multiqc4=$results/multiqc_4
multiqc5=$results/multiqc_5


### run fastqc 
mkdir -p $results1 $results2 $results3 $results4 $multiqc1 $multiqc2 $multiqc $multiqc4

ls  $shared_data1 | xargs -I {} --max-procs=20 fastqc $shared_data1/{}/{}_1.fq.gz  $shared_data1/{}/{}_2.fq.gz -o $results1/
ls  $shared_data2 | xargs -I {} --max-procs=20 fastqc $shared_data2/{}/{}_1.fq.gz  $shared_data2/{}/{}_2.fq.gz -o $results2/ 
ls  $shared_data3 | xargs -I {} --max-procs=40 fastqc $shared_data3/{}/{}_1.fq.gz  $shared_data3/{}/{}_2.fq.gz -o $results3/
ls  $shared_data4 | xargs -I {} --max-procs=40 fastqc $shared_data4/{}/{}_1.fq.gz  $shared_data4/{}/{}_2.fq.gz -o $results4/ 
ls  $shared_data5 | xargs -I {} --max-procs=40 fastqc $shared_data5/{}/{}_1.fq.gz  $shared_data5/{}/{}_2.fq.gz -o $results5/ 

### run multiQC ####

multiqc $results1/*  -o $multiqc1 
multiqc $results2/*  -o $multiqc2
multiqc $results3/*  -o $multiqc3
multiqc $results4/*  -o $multiqc4
multiqc $results5/*  -o $multiqc5


###  triming stuff ###

##### trim_galore ######

trim_out1=$results/trim_out1
trim_out2=$results/trim_out2
trim_out3=$results/trim_out3
trim_out4=$results/trim_out4
trim_out5=$results/trim_out5

multiqc_trim1=$results/multiqc_trim1
multiqc_trim2=$results/multiqc_trim2
multiqc_trim3=$results/multiqc_trim3
multiqc_trim4=$results/multiqc_trim4
multiqc_trim5=$results/multiqc_trim5


mkdir -p $trim_out3  $trim_out1 $trim_out2 $trim_out4 $trim_out5 $multiqc_trim1 $multiqc_trim2 $multiqc_trim3 $multiqc_trim4 $multiqc_trim5

ls  $shared_data1 | xargs -I {} --max-procs=40 trim_galore --fastqc --paired  --clip_R1 13 --clip_R2 13  $shared_data1/{}/{}_1.fq.gz $shared_data1/{}/{}_2.fq.gz -o $trim_out1  
ls  $shared_data2 | xargs -I {} --max-procs=40 trim_galore --fastqc --paired  --clip_R1 13 --clip_R2 13  $shared_data2/{}/{}_1.fq.gz $shared_data2/{}/{}_2.fq.gz -o $trim_out2
ls  $shared_data3 | xargs -I {} --max-procs=40 trim_galore --fastqc --paired  --clip_R1 13 --clip_R2 13  $shared_data3/{}/{}_1.fq.gz $shared_data3/{}/{}_2.fq.gz -o $trim_out3
ls  $shared_data4 | xargs -I {} --max-procs=20 trim_galore --fastqc --paired  --clip_R1 9 --clip_R2 9  $shared_data4/{}/{}_1.fq.gz $shared_data4/{}/{}_2.fq.gz -o $trim_out4
ls  $shared_data5 | xargs -I {} --max-procs=20 trim_galore --fastqc --paired  --clip_R1 9 --clip_R2 9  $shared_data5/{}/{}_1.fq.gz $shared_data5/{}/{}_2.fq.gz -o $trim_out5


multiqc $trim_out1/*  -o $multiqc_trim1
multiqc $trim_out2/*  -o $multiqc_trim2
multiqc $trim_out3/*  -o $multiqc_trim3
multiqc $trim_out4/*  -o $multiqc_trim4
multiqc $trim_out5/*  -o $multiqc_trim5


exit 0
