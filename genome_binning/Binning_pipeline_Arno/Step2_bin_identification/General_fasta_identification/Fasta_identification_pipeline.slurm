#!/bin/bash

#SBATCH -p compute
#SBATCH --output=bin_identification_pipe%a.out
#SBATCH -t 2-0
#SBATCH --mem=48G
#SBATCH -c 32
#SBATCH --array=0-4

#Performs megan analysis, and extracts and identifies rRNA genes from the provided fasta files

ml bioinfo-ugrp-modules
ml ncbi-blast
ml DB/blastDB/ncbi/238
ml DebianMed
ml barrnap

#######################################
#INPUT VARIABLES !!CHANGE AS REQUIRED!!
#######################################
#Directory containing the fasta files
fasta_dir=/flash/HusnikU/Arno/Bin_pipeline

#Output directory
id_out_dir=/flash/HusnikU/Arno/Bin_identification

#####################################
#SSU EXTRACTION AND IDENTIFICATION
####################################
#Run barrnap on all fastas, set to bacteria, archea, and eukaryotes
mkdir ${id_out_dir}/barrnap
cd ${fasta_dir}
bins=*.fa

echo "-----------------------------------------------"
echo "Running pipeline on:"
echo ${bins[${SLURM_ARRAY_TASK_ID}]}
echo "-----------------------------------------------"
echo "-----------------------------------------------"
echo "Running Barrnap"
echo "-----------------------------------------------"
barrnap ${bins[${SLURM_ARRAY_TASK_ID}]} --kingdom bac --outseq ${id_out_dir}/barrnap/${bins[${SLURM_ARRAY_TASK_ID}]}_bac.fasta
barrnap ${bins[${SLURM_ARRAY_TASK_ID}]} --kingdom arc --outseq ${id_out_dir}/barrnap/${bins[${SLURM_ARRAY_TASK_ID}]}_arc.fasta
barrnap ${bins[${SLURM_ARRAY_TASK_ID}]} --kingdom euk --outseq ${id_out_dir}/barrnap/${bins[${SLURM_ARRAY_TASK_ID}]}_euk.fasta
cat ${id_out_dir}/barrnap/${bins[${SLURM_ARRAY_TASK_ID}]}_bac.fasta ${id_out_dir}/barrnap/${bins[${SLURM_ARRAY_TASK_ID}]}_arc.fasta ${id_out_dir}/barrnap/${bins[${SLURM_ARRAY_TASK_ID}]}_euk.fasta > ${id_out_dir}/barrnap/${bins[${SLURM_ARRAY_TASK_ID}]}_total.fasta

#Remove the intermediate barrnap outputs
rm ${id_out_dir}/barrnap/${bins[${SLURM_ARRAY_TASK_ID}]}_bac.fasta
rm ${id_out_dir}/barrnap/${bins[${SLURM_ARRAY_TASK_ID}]}_arc.fasta
rm ${id_out_dir}/barrnap/${bins[${SLURM_ARRAY_TASK_ID}]}_euk.fasta

echo "-----------------------------------------------"
echo "BLASTing SSU sequences"
echo "-----------------------------------------------"

#blast extracted 16S sequences against ncbi NT
blastn -db nt -outfmt 6 -query ${id_out_dir}/barrnap/${bins[${SLURM_ARRAY_TASK_ID}]}_total.fasta -out ${id_out_dir}/barrnap/${bins[${SLURM_ARRAY_TASK_ID}]}_SSU.blast -num_threads 32

##################################
#MEGAN ANALYSES
##################################
echo "-----------------------------------------------"
echo "BLASTing bins against NT for MEGAN"
echo "-----------------------------------------------"

#run BLASTn on the fasta files
mkdir ${id_out_dir}/MEGAN
blastn -db nt -query ${fasta_dir}/${bins[${SLURM_ARRAY_TASK_ID}]} -out ${id_out_dir}/MEGAN/${bins[${SLURM_ARRAY_TASK_ID}]}.blast -num_threads 32 -outfmt 5

echo "-----------------------------------------------"
echo "Performing LCA on the BLAST results, and MEGANizing the output"
echo "-----------------------------------------------"

#convert blast to .rma file for reading in MEGAN
/apps/unit/HusnikU/megan/tools/blast2rma -i ${id_out_dir}/MEGAN/${bins[${SLURM_ARRAY_TASK_ID}]}.blast  -o ${id_out_dir}/MEGAN/${bins[${SLURM_ARRAY_TASK_ID}]}.rma  -f BlastXML -bm BlastN -c

#Run MEGAN BLAST2LCA script, and summarize LCA using python script
/apps/unit/HusnikU/megan/tools/blast2lca -i ${id_out_dir}/MEGAN/${bins[${SLURM_ARRAY_TASK_ID}]}.blast -f BlastXML -m BlastN -o ${id_out_dir}/MEGAN/${bins[${SLURM_ARRAY_TASK_ID}]}.lca
python3 /bucket/HusnikU/Unit_Members/Arno/Scripts/MEGAN_LCA_summary.py ${id_out_dir}/MEGAN/${bins[${SLURM_ARRAY_TASK_ID}]}.lca ${id_out_dir}/MEGAN/${bins[${SLURM_ARRAY_TASK_ID}]}_LCA_summary 5





