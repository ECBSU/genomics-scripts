#!/bin/bash
#SBATCH -p compute
#SBATCH -t 0-3
#SBATCH --mem=15G
#SBATCH -c 32
#SBATCH -n 1
#SBATCH --job-name=AvP

#activating conda environment
source /home/y/yong-phua/conda_ini.sh 
conda activate avp

#Diamond blastp hits file, output format should remain as written here
diamond blastp -q proteins.faa -d /flash/HusnikU/databases/AvP_swissprot/uniprot_sprot.fasta.dmnd --evalue 1e-5 --max-target-seqs 500 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids --out similarity.out

#AvP pipeline change out_dir to path to output folder
/home/y/yong-phua/AvP/aux_scripts/calculate_ai.py -i similarity.out -x groups.yaml
/home/y/yong-phua/AvP/avp prepare -a similarity.out_ai.out -o out_dir -f proteins.faa -b similarity.out -x groups.yaml -c config.yaml
/home/y/yong-phua/AvP/avp detect -i out_dir/mafftgroups/ -o out_dir -g out_dir/groups.tsv -t out_dir/tmp/taxonomy_nexus.txt -c config.yaml
/home/y/yong-phua/AvP/avp evaluate -i out_dir/mafftgroups/ -t out_dir/fasttree_tree_results.txt -o out_dir -c config.yaml