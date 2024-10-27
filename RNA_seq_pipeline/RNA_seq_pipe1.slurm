#!/bin/bash
#SBATCH -p compute
#SBATCH -t 3-0
#SBATCH --mem=300G
#SBATCH -c 96
#SBATCH -n 1
#SBATCH --job-name=RNASPAdes
#SBATCH --output=RNAspades.log

out_dir=/path/SPAdes
in_path=/path/to/transcript

#set database for BUSCO
busco_db=/bucket/.deigo/HusnikU/Unit_Members/Phua_internship_files/Sino_prj/transptome/clean_transcriptome/Trinity/sino_transcriptome_SL1234/busco_downloads/lineages/alveolata_odb10

#set organisms type for targetp
organism_type='pl'

#Run SPAdes
ml python/3.11.4
/apps/unit/HusnikU/SPAdes-4.0.0-Linux/bin/spades.py -o $out_dir --pe1-1 $in_path/LOY1_1_trimed.fastq.gz --pe1-2 $in_path/LOY1_2_trimed.fastq.gz --pe2-1 $in_path/LOY2_1_trimed.fastq.gz --pe2-2 $in_path/LOY2_2_trimed.fastq.gz --pe3-1 $in_path/LOY3_1_trimed.fastq.gz --pe3-2 $in_path/LOY3_2_trimed.fastq.gz --rna -t 96 -m 300

#run busco
module load BUSCO
busco -m transcriptome -i $out_dir/transcripts.fasta -c 2 -o busco -l $busco_db --offline -f
module unload BUSCO

#Run transdecoder
/apps/unit/HusnikU/TransDecoder-TransDecoder-v5.5.0/TransDecoder.LongOrfs -t $out_dir/transcripts.fasta
/apps/unit/HusnikU/TransDecoder-TransDecoder-v5.5.0/TransDecoder.Predict -t $out_dir/transcripts.fasta

#Run eggnog
source /home/a/arno-hagenbeek/miniconda3/etc/profile.d/conda.sh
conda activate eggnog-mapper
emapper.py --data_dir /bucket/HusnikU/Databases/Eggnog_db/Eggnog_Diamond_db/ --cpu 96 -i $out_dir/transcripts.fasta.transdecoder.pep -o transcripts

#Run Targetp
/hpcshare/appsunit/HusnikU/targetp-2.0/bin/targetp -fasta $out_dir/transcripts.fasta.transdecoder.pep -org $organism_type