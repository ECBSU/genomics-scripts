#!/bin/bash

#SBATCH -p compute
#SBATCH --output=Busco_tree_build.out
#SBATCH -t 0-12 #This number will need to be increased for larger aligments
#SBATCH --mem=40G
#SBATCH -c 64 #!!!Don't change, as currently this number is hardcoded for the incorporated tools. 

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

#########################
#Run script
#########################
python /bucket/HusnikU/Unit_Members/Arno/Scripts/Busco_multigene_tree.py $input_dir $output_dir $fraction
