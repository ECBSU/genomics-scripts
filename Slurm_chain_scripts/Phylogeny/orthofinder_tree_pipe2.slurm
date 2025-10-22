#!/bin/bash

#SBATCH -p compute
#SBATCH --output=orthofinder_tree_part2.out
#SBATCH -t 0-5
#SBATCH --mem=32G
#SBATCH -c 32

ml mafft
indir=$1/orthofinder/*/Single_Copy_Orthologue_Sequences 
outdir=$1/trimal_alignment
#Empty the outdir if it exists, to prevent adding genes from previous runs
rm -r $outdir
mkdir $outdir
cd $indir
for i in *.fa; do mafft --localpair --maxiterate 1000 --thread 32 ${i} > ${outdir}/${i}_aligned; /home/a/arno-hagenbeek/trimal/trimal/source/trimal -in ${outdir}/${i}_aligned -out ${outdir}/${i}_trimmed -automated1; done
python orthofinder_tree_pipe_gene_concat.py $outdir ${outdir}/concatenated_aligned_trimmed_SGC.fasta
