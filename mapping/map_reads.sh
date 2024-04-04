#!/bin/bash

mouse_index=../../shared_data/mouse_index/mouse_salmon_index
data1=../results/run_QC/trim_out1
data2=../results/run_QC/trim_out2
data3=../results/run_QC/trim_out3
data4=../results/run_QC/trim_out4
data5=../results/run_QC/trim_out5

results1=../results/map1
results2=../results/map2
results3=../results/map3
results4=../results/map4
results5=../results/map5


mkdir -p $results3  $results2 $results1 $results4 $results5
map_recap=../results/map_recap_validate.txt


ls $data1  | grep 1.fq.gz$ | sed 's/_.*//' | xargs -I {} --max-procs=5 bash -c  "salmon quant -i $mouse_index -l IU -1 <(zcat $data1/{}_1_val_1.fq.gz) -2 <(zcat $data1/{}_2_val_2.fq.gz) -o $results1/{} --validateMappings"

ls $data2 | grep 1.fq.gz$  | sed 's/_.*//' | xargs -I {} --max-procs=5 bash -c  "salmon quant -i $mouse_index -l IU -1 <(zcat $data2/{}_1_val_1.fq.gz) -2 <(zcat $data2/{}_2_val_2.fq.gz) -o $results2/{} --validateMappings"

ls $data3  | grep 1.fq.gz$ | sed 's/_.*//' | xargs -I {} --max-procs=5 bash -c  "salmon quant -i $mouse_index -l IU -1 <(zcat $data3/{}_1_val_1.fq.gz) -2 <(zcat $data3/{}_2_val_2.fq.gz)  -o $results3/{} --validateMappings"

ls $data4  | grep 1.fq.gz$ | sed 's/_.*//' | xargs -I {} --max-procs=5 bash -c  "salmon quant -i $mouse_index -l IU -1 <(zcat $data4/{}_1_val_1.fq.gz) -2 <(zcat $data4/{}_2_val_2.fq.gz)  -o $results4/{} --validateMappings"

ls $data5  | grep 1.fq.gz$ | sed 's/_.*//' | xargs -I {} --max-procs=5 bash -c  "salmon quant -i $mouse_index -l IU -1 <(zcat $data5/{}_1_val_1.fq.gz) -2 <(zcat $data5/{}_2_val_2.fq.gz)  -o $results5/{} --validateMappings"


grep -R " rate" "$results1" > $map_recap
grep -R " rate" "$results2" >> $map_recap
grep -R " rate" "$results3" >> $map_recap
grep -R " rate" "$results4" >> $map_recap
grep -R " rate" "$results5" >> $map_recap

exit 0
