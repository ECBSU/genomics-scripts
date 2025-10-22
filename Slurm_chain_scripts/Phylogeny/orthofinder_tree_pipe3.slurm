#!/bin/bash

#SBATCH -p compute
#SBATCH --output=orthofinder_tree_part2.out
#SBATCH -t 4-0
#SBATCH --mem=256G
#SBATCH -c 128

indir=$1/trimal_alignment
outdir=$1/tree
ml bioinfo-ugrp-modules
ml DebianMed
ml iqtree
mkdir $outdir
cd $outdir
cp ${indir}/concatenated_aligned_trimmed_SGC.fasta . 
iqtree2 -s ./concatenated_aligned_trimmed_SGC.fasta  -bb 1000 -T 128 
