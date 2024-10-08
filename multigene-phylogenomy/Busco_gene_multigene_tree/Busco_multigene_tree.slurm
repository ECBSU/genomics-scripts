#!/bin/bash

#SBATCH -p compute
#SBATCH --output=Busco_tree_build.out
#SBATCH -t 0-12 #This number will need to be increased for larger aligments
#SBATCH --mem=40G
#SBATCH -c 64 #!!!If you lower this value, also lower in the tool commands found below

#This script finds shared busco genes in busco output, and aligns, trims, and 
#concatenates them for input into iqtree. Then runs iqtree to generate a phylogenetic tree
##########################################################################################################################################
ml mafft
ml bioinfo-ugrp-modules
ml DebianMed
ml iqtree

#########################
#INPUT VARIABLES, CHANGE AS NEEDED!!!
#########################
input_dir=/bucket/HusnikU/Unit_Members/Pradeep/Multigene_Phylogeny/Busco/Facultative_bacterial_endosymbionts/  #Directory containing your busco outputs
output_dir= #Directory to which to write the output
fraction=0.8 #fraction of organisms a gene needs to be present in for it to be included in the analysis
MAFFT_cmd="mafft --globalpair --maxiterate 1000 --thread 64 {} > {}" #Command used to execute mafft. Arguments can be modified or added freely. The {} denote input and output respectively, please leave unmodified.
trimal_cmd="/home/a/arno-hagenbeek/trimal/trimal/source/trimal -in {} -out {} -automated1" #Command used to execute trimal. Arguments can be modified or added freely. The {} denote input and output respectively, please leave unmodified.
iqtree_cmd="iqtree2 -s {} --prefix {} -bb 1000 -T 64" #Command used to execute iqtree. Arguments can be modified or added freely. The {} denote input and output respectively, please leave unmodified.

#########################
#Run script
#########################
python /bucket/HusnikU/Unit_Members/Arno/Scripts/Busco_multigene_tree.py $input_dir $output_dir $fraction $MAFFT_cmd $trimal_cmd $iqtree_cmd
