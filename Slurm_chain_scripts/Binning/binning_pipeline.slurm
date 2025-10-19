#!/bin/bash

#SBATCH -p compute
#SBATCH --output=bin_pipe_%a.out
#SBATCH -t 4-0
#SBATCH --mem=64G
#SBATCH -c 32


#Load modules
ml bioinfo-ugrp-modules
ml DebianMed
ml fastp
ml samtools
ml bowtie2
ml metabat
ml DB/diamondDB/ncbi

#############################################
#Insert directories containing only the assemblies and reads for binning
assembly_dir=${2}/spades
reads_dir=${2}/fastp
#Output directory
out_dir=${2}/binning
#Conda installation
conda_ex=$3

#INPUT READING
#############################################
#Creates array of assemblies to iterate over
cd ${reads_dir}
assemblies=(*R1_001.fastq.gz) #Change this to the extension of your assemblies
#create and move to out_dir
mkdir $out_dir
cd ${out_dir}
mkdir ${assemblies[${SLURM_ARRAY_TASK_ID}]}_intermediate_files
#creates arrays of forward and reverse reads. 
cd $reads_dir
reads1=(*R1_001.fastq.gz) 
reads2=(*R2_001.fastq.gz) 
echo "Finished finding reads"


cd ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_intermediate_files

#READ ALIGNMENT AND BAM FILE GENERATION
############################################
#Run bowtie for alignment
echo "-----------------------------------------------"
echo "Running Bowtie2"
echo "-----------------------------------------------"
bowtie2-build $assembly_dir/${assemblies[${SLURM_ARRAY_TASK_ID}]}/scaffolds.fasta  ./${assemblies[${SLURM_ARRAY_TASK_ID}]}
bowtie2 -x ./${assemblies[${SLURM_ARRAY_TASK_ID}]} -1 ${reads_dir}/${reads1[${SLURM_ARRAY_TASK_ID}]} -2 ${reads_dir}/${reads2[${SLURM_ARRAY_TASK_ID}]} -S ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_intermediate_files/${assemblies[${SLURM_ARRAY_TASK_ID}]}.sam -p 16

#Dont use the trimmed reads anymore, so delete to make space
rm ${reads_dir}/${reads1[${SLURM_ARRAY_TASK_ID}]}
rm ${reads_dir}/${reads2[${SLURM_ARRAY_TASK_ID}]}

#Sort and compress sam file into bam file
echo "-----------------------------------------------"
echo "Running samtools"
echo "-----------------------------------------------"
samtools sort -@16 -O BAM -o ./${assemblies[${SLURM_ARRAY_TASK_ID}]}.bam ./${assemblies[${SLURM_ARRAY_TASK_ID}]}.sam
samtools index ./${assemblies[${SLURM_ARRAY_TASK_ID}]}.bam ./${assemblies[${SLURM_ARRAY_TASK_ID}]}.bai

#remove sam file and bowtie indices to save space
rm ./${assemblies[${SLURM_ARRAY_TASK_ID}]}.sam
rm ./${assemblies[${SLURM_ARRAY_TASK_ID}]}*.bt2


#RUN METABAT
#############################################
#Create folder for metabat output (per sample)
cd $out_dir
mkdir ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_metabat
mkdir ${assemblies[${SLURM_ARRAY_TASK_ID}]}_intermediate_files/Metabat
#Create metabat coverage input, and run metabat
jgi_summarize_bam_contig_depths --outputDepth ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_intermediate_files/Metabat/${assemblies[${SLURM_ARRAY_TASK_ID}]}_depth.txt ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_intermediate_files/${assemblies[${SLURM_ARRAY_TASK_ID}]}.bam
echo "-----------------------------------------------"
echo "Running Metabat2"
echo "-----------------------------------------------"
metabat2 -i ${assembly_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]} -a ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_intermediate_files/Metabat/${assemblies[${SLURM_ARRAY_TASK_ID}]}_depth.txt -o ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_metabat/metabatbin -t 16

#RUN CONCOCT
#############################################
#Activate concoct conda env
source ${conda_ex}
conda activate /apps/unit/HusnikU/Conda-envs/concoct
concoct -v
mkdir ${assemblies[${SLURM_ARRAY_TASK_ID}]}_concoct
mkdir ${assemblies[${SLURM_ARRAY_TASK_ID}]}_intermediate_files/Concoct
#Run concoct and its associated python scripts
echo "-----------------------------------------------"
echo "Running Concoct and its scripts"
echo "-----------------------------------------------"
cut_up_fasta.py $assembly_dir/${assemblies[${SLURM_ARRAY_TASK_ID}]} -c 10000 -o 0 --merge_last -b ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_intermediate_files/Concoct/${assemblies[${SLURM_ARRAY_TASK_ID}]}_10k.bed > ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_intermediate_files/Concoct/${assemblies[${SLURM_ARRAY_TASK_ID}]}_contigs_10k.fa
concoct_coverage_table.py ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_intermediate_files/Concoct/${assemblies[${SLURM_ARRAY_TASK_ID}]}_10k.bed ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_intermediate_files/${assemblies[${SLURM_ARRAY_TASK_ID}]}.bam > ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_intermediate_files/Concoct/${assemblies[${SLURM_ARRAY_TASK_ID}]}_coverage_table.tsv
concoct --composition_file ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_intermediate_files/Concoct/${assemblies[${SLURM_ARRAY_TASK_ID}]}_contigs_10k.fa --coverage_file ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_intermediate_files/Concoct/${assemblies[${SLURM_ARRAY_TASK_ID}]}_coverage_table.tsv -b ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_intermediate_files/Concoct/${assemblies[${SLURM_ARRAY_TASK_ID}]}_concoct -t 16
merge_cutup_clustering.py ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_intermediate_files/Concoct/${assemblies[${SLURM_ARRAY_TASK_ID}]}_concoct_clustering_gt1000.csv > ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_intermediate_files/Concoct/${assemblies[${SLURM_ARRAY_TASK_ID}]}_concoct_clustering_merged.csv
extract_fasta_bins.py $assembly_dir/${assemblies[${SLURM_ARRAY_TASK_ID}]} ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_intermediate_files/Concoct/${assemblies[${SLURM_ARRAY_TASK_ID}]}_concoct_clustering_merged.csv --output_path ${assemblies[${SLURM_ARRAY_TASK_ID}]}_concoct

#RUN SEMIBIN
#############################################
#Activate semibin conda env
conda deactivate
conda activate /apps/unit/HusnikU/Conda-envs/semibin
#Run semibin
echo "-----------------------------------------------"
echo "Running Semibin"
echo "-----------------------------------------------"
SemiBin2 -v
mkdir ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_semibin/
mkdir ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_intermediate_files/Semibin
SemiBin2 single_easy_bin -i $assembly_dir/${assemblies[${SLURM_ARRAY_TASK_ID}]} -b ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_intermediate_files/${assemblies[${SLURM_ARRAY_TASK_ID}]}.bam --environment global -o ${assemblies[${SLURM_ARRAY_TASK_ID}]}_intermediate_files/Semibin 
#Semibin was run in the intermediate file folder since it generates a bunch of files. Move the bins to the appropriate folder 
mv ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_intermediate_files/Semibin/output_bins/*.gz ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_semibin/
#Unzip the bin fasta files
gunzip ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_semibin/*.gz
#Remove the now empty output bin folder
rm -r ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_intermediate_files/Semibin/output_bins/

#RUN ROSELLA
#############################################
#Activate rosella conda env
conda deactivate 
conda activate /bucket/HusnikU/Conda-envs/rosella/
#Run Rosella
echo "-----------------------------------------------"
echo "Running Rosella"
echo "-----------------------------------------------"
mkdir ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_rosella/
mkdir ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_intermediate_files/Rosella
rosella recover -r $assembly_dir/${assemblies[${SLURM_ARRAY_TASK_ID}]} -o ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_intermediate_files/Rosella -b ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_intermediate_files/${assemblies[${SLURM_ARRAY_TASK_ID}]}.bam
#Rosella generates non-bin fasta files, with naming similar to the bins.
#To avoid these fastas from interfering downstream, they'll be moved to a separate folder 
mkdir ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_intermediate_files/Rosella/unbinned_sequences
mv ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_intermediate_files/Rosella/*unbinned.fna ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_intermediate_files/Rosella/unbinned_sequences
mv ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_intermediate_files/Rosella/rosella_*.fna ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_rosella/

#Rosella produces both bins and refined bins. The refined bins appear to be higher quality and both will be included in the analysis

#RUN DASTOOL
#############################################
#Activate Dastool conda environment 
conda deactivate
conda activate /apps/unit/HusnikU/Conda-envs/das-tool
mkdir ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_intermediate_files/DASTool/
#Generate contig2bin inputs .tsv files for dastool for the tools that don't generate their own
for i in ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_metabat/*.fa; do grep "^>" $i | sed 's/>.*/&\t'${i##*/}' /' | sed 's/>//' >> ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_intermediate_files/DASTool/metabat_bins.tsv; done
for i in ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_concoct/*.fa; do grep "^>" $i | sed 's/>.*/&\t'${i##*/}' /' | sed 's/>//' >> ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_intermediate_files/DASTool/concoct_bins.tsv; done
for i in ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_rosella/rosella_*.fna; do grep "^>" $i | sed 's/>.*/&\t'${i##*/}' /' | sed 's/>//' >> ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_intermediate_files/DASTool/rosella_bins.tsv; done

#Semibin generates its own contig2bin .tsv files, but it contains headers, which DASTools doesnt accept.
#Remove the header files 
tail -n +2 ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_intermediate_files/Semibin/contig_bins.tsv > ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_intermediate_files/DASTool/semibin_bins.tsv

cd ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_intermediate_files/DASTool
#Run DASTool
echo "-----------------------------------------------"
echo "Running DASTool"
echo "-----------------------------------------------"
DAS_Tool -i metabat_bins.tsv,concoct_bins.tsv,rosella_bins.tsv,semibin_bins.tsv -c $assembly_dir/${assemblies[${SLURM_ARRAY_TASK_ID}]} -o ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_intermediate_files/DASTool/${assemblies[${SLURM_ARRAY_TASK_ID}]} --write_bin_evals --write_bins -t 16
#Move refined bins to the DASTool bin folder
mv ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_intermediate_files/DASTool/${assemblies[${SLURM_ARRAY_TASK_ID}]}_DASTool_bins ${out_dir}/${assemblies[${SLURM_ARRAY_TASK_ID}]}_DASTool
