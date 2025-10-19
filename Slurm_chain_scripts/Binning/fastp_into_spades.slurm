#!/bin/bash

#SBATCH -p compute
#SBATCH --output=fastp_into_spades%a.out
#SBATCH -t 4-0
#SBATCH --mem=480G
#SBATCH -c 64

#Read in the fastq files
cd $1
R1=(*R1_001.fastq.gz)
R2=(*R2_001.fastq.gz)
cd $2

#Call the fastp module
ml bioinfo-ugrp-modules
ml DebianMed
ml fastp
#Run fastp to trim the fastq files
mkdir ${2}/fastp
echo ${R1[${SLURM_ARRAY_TASK_ID}]}
echo ${R2[${SLURM_ARRAY_TASK_ID}]}
fastp -i ${1}/${R1[${SLURM_ARRAY_TASK_ID}]} -I ${1}/${R2[${SLURM_ARRAY_TASK_ID}]} -o ${2}/fastp/${R1[${SLURM_ARRAY_TASK_ID}]} -O ${2}/fastp/${R2[${SLURM_ARRAY_TASK_ID}]} -w 64 -z 9
#Run spades
ml purge
mkdir ${2}/spades
/apps/unit/HusnikU/SPAdes-4.0.0-Linux/bin/spades.py -o ${2}/spades/${R1[${SLURM_ARRAY_TASK_ID}]} -1 ${2}/fastp/${R1[${SLURM_ARRAY_TASK_ID}]} -2 ${2}/fastp/${R2[${SLURM_ARRAY_TASK_ID}]} -t 64 -m 480


