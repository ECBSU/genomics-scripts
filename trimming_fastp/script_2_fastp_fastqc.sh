#!/bin/bash
#SBATCH -p compute
#SBATCH -t 01:00:00
#SBATCH -c 64
#SBATCH --mem=128G

#metadata_version: 1.0.0
#script_version: 1.0.0
#software: [fastp: 0.23.2, fastqc: 0.12.1]
#environment: hpc
#author: Vasyl Vaskivskyi

#script_info: 
# Performs adapter trimming and fastqc before and after trimming.
# https://github.com/OpenGene/fastp
# https://www.bioinformatics.babraham.ac.uk/projects/fastqc/


#inputs:
# - paths to reads
# - output directory


#outputs:
# - trimmed reads
# - fastp reports
# - fastqc reports
# - fastqc archives

eval "$(conda shell.bash hook)"
conda activate smartseq2 

r1=$1
r2=$2
outdir=$3

out_r1=$outdir/"tr_"$(basename $r1)
out_r2=$outdir/"tr_"$(basename $r2)

mkdir -p $outdir


echo "FASTQC before trimming"
fastqc $r1 $r2 -o $outdir -t 64 --memory 10000

echo "TRIMMING of adapters"
fastp -w 16 -i $r1 -o $out_r1 -I $r2 -O $out_r2 --json $outdir/fastp.json --html $outdir/fastp.html --adapter_fasta /flash/HusnikU/databases/smartseq2_oligos.fasta

echo "FASTQC after trimming"
fastqc $out_r1  $out_r2 -o $outdir -t 64 --memory 10000
