#!/bin/bash



mouse_index=../../shared_data/mouse_salmon_index
mkdir -p  $mouse_index

### index and prepare  mouse transcriptome

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.transcripts.fa.gz -P $mouse_index # transcriptome
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz -P $mouse_index # genome

grep "^>" <(gunzip -c $mouse_index/GRCm38.primary_assembly.genome.fa.gz) | cut -d " " -f 1 > $mouse_index/decoys.txt
sed -i.bak -e 's/>//g' $mouse_index/decoys.txt

cat $mouse_index/gencode.vM25.transcripts.fa.gz $mouse_index/GRCm38.primary_assembly.genome.fa.gz > $mouse_index/gentrome.fa.gz
salmon index -t $mouse_index/gentrome.fa.gz -d $mouse_index/decoys.txt -p 12 -i $mouse_index/mouse_salmon_index --gencode

exit 0
