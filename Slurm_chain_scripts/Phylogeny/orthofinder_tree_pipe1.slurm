#!/bin/bash

#SBATCH -p compute
#SBATCH --output=orthofinder_tree_part1.out
#SBATCH -t 2-0
#SBATCH --mem=128G
#SBATCH -c 128


indir=$1
ortho_dir=$2/orthofinder
prokka_dir=$2/prokka
extension=$3
source $4
conda activate $5
cd $indir
mkdir $2
mkdir $prokka_dir
#Generate the .faa files using prokka
for i in *.${extension}; do prokka --force --outdir $prokka_dir --prefix $i --cpus 128 --compliant $i ; done
#Rename the protein codes to sample names, so the samples show in the tree downstream
python orthofinder_tree_pipe_gene_rename.py $prokka_dir
#Run orthofinder
ml ncbi-blast
rm -r $ortho_dir
ulimit -n 14500 #Set open file limit higher so larger jobs can still run
/apps/unit/HusnikU/OrthoFinder/orthofinder -f $prokka_dir -o $ortho_dir -S blast -t 128
