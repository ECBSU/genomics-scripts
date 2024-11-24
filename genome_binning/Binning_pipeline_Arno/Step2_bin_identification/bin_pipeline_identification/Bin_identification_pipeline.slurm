#!/bin/bash

#SBATCH -p compute
#SBATCH --output=bin_identification_pipe%a.out
#SBATCH -t 4-0
#SBATCH --mem=48G
#SBATCH -c 32
#SBATCH --array=0-4

#Performs megan analysis, and extracts and identifies rRNA genes from the bins

ml bioinfo-ugrp-modules
ml ncbi-blast
ml DB/blastDB/ncbi/238
ml DebianMed
ml barrnap

#######################################
#INPUT VARIABLES !!CHANGE AS REQUIRED!!
#######################################
#Directory containing the binning pipeline output
bin_pipeline_out_dir=/path/to/binning/pipeline/output/

#Output directory
id_out_dir=/path/to/identification/pipeline/output

######################################
#INPUT PREPARATION
######################################
#Creates a single directory containing all successful bins for the sample
cd ${bin_pipeline_out_dir}
bin_dir=(*_concoct)
mkdir ${id_out_dir}/
mkdir ${id_out_dir}/${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}
cp ${bin_dir[${SLURM_ARRAY_TASK_ID}]}/*.fa ${id_out_dir}/${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}
cp ${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_semibin/*.fa ${id_out_dir}/${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}/
cp ${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_metabat/*.fa ${id_out_dir}/${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}/
cp ${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_rosella/*.fa ${id_out_dir}/${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}/
cp ${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_DASTool/*.fa ${id_out_dir}/${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}/

echo "-----------------------------------------------"
echo "Running pipeline on:"
echo ${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}
echo "-----------------------------------------------"

#####################################
#SSU EXTRACTION AND IDENTIFICATION
####################################

echo "-----------------------------------------------"
echo "Running Barrnap"
echo "-----------------------------------------------"

#Run barrnap on all bins, set to bacteria, archea, and eukaryotes
mkdir ${id_out_dir}/${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_barrnap
cd ${id_out_dir}/${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}/
bins=*.fa
for i in ${bins}; do barrnap ${i} --kingdom bac --outseq ${id_out_dir}/${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_barrnap/${i}_bac.fasta; done
for i in ${bins}; do barrnap ${i} --kingdom arc --outseq ${id_out_dir}/${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_barrnap/${i}_arc.fasta; done
for i in ${bins}; do barrnap ${i} --kingdom euk --outseq ${id_out_dir}/${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_barrnap/${i}_euk.fasta; done
for i in ${bins}; do cat ${id_out_dir}/${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_barrnap/${i}_bac.fasta ${id_out_dir}/${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_barrnap/${i}_arc.fasta ${id_out_dir}/${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_barrnap/${i}_euk.fasta > ${id_out_dir}/${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_barrnap/${i}_total.fasta; done

#Remove the intermediate barrnap outputs
rm ${id_out_dir}/${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_barrnap/*_bac.fasta
rm ${id_out_dir}/${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_barrnap/*_arc.fasta
rm ${id_out_dir}/${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_barrnap/*_euk.fasta

echo "-----------------------------------------------"
echo "BLASTing SSU sequences"
echo "-----------------------------------------------"

#blast extracted 16S sequences against ncbi NT
for i in ${bins}; do blastn -db nt -outfmt "6 qseqid sseqid sscinames pident length mismatch gapopen qstart qend sstart send evalue bitscore" -query ${id_out_dir}/${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_barrnap/${i}_total.fasta -out ${id_out_dir}/${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_barrnap/${i}_SSU_blast -num_threads 32; done

##################################
#MEGAN ANALYSES
##################################

echo "-----------------------------------------------"
echo "BLASTing bins against NT for MEGAN"
echo "-----------------------------------------------"

#run BLASTn on the bin files
mkdir ${id_out_dir}/${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_MEGAN
for i in $bins; do blastn -db nt -query ${id_out_dir}/${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}/${i} -out ${id_out_dir}/${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_MEGAN/${i}.blast -num_threads 32 -outfmt 5; done

echo "-----------------------------------------------"
echo "Performing LCA on the BLAST results, and MEGANizing the output"
echo "-----------------------------------------------"

#convert blast to .rma file for reading in MEGAN
for i in $bins; do /apps/unit/HusnikU/megan/tools/blast2rma -i ${id_out_dir}/${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_MEGAN/${i}.blast  -o ${id_out_dir}/${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_MEGAN/${i}.rma  -f BlastXML -bm BlastN -c; done

#Run MEGAN BLAST2LCA script, and summarize LCA using python script
for i in $bins; do /apps/unit/HusnikU/megan/tools/blast2lca -i ${id_out_dir}/${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_MEGAN/${i}.blast -f BlastXML -m BlastN -o ${id_out_dir}/${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_MEGAN/${i}.lca; done
for i in $bins; do python3 /bucket/HusnikU/Unit_Members/Arno/Scripts/MEGAN_LCA_summary.py ${id_out_dir}/${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_MEGAN/${i}.lca ${id_out_dir}/${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_MEGAN/${i}_lca_summary 5; done

