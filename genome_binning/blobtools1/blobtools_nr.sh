#running blastn for scaffolds vs nt database
/hpcshare/appsunit/HusnikU/ncbi-blast-2.11.0+/bin/blastn -query scaffolds.fasta -db /flash/HusnikU/databases/nt/nt -outfmt '6 qseqid staxids bitscore std' -max_target_seqs 1 -max_hsps 1 -evalue 1e-25 -num_threads 64 > scaffolds_vs_blastn.out

#combining the output from nt database to create the json file which will have the data to write plots/tables
/apps/unit/HusnikU/blobtools/blobtools create -i scaffolds.fasta -y spades -t scaffolds_vs_blastn.out --db /home/y/yong-phua/blobtaxid/nodesDB.txt 

#plotting blobplot at different taxonomical rank 
/apps/unit/HusnikU/blobtools/blobtools plot -i blobDB.json
/apps/unit/HusnikU/blobtools/blobtools plot -i blobDB.json -r order

#writing tables with order and family ranks
/apps/unit/HusnikU/blobtools/blobtools view -i blobDB.json -r phylum -r order -r family