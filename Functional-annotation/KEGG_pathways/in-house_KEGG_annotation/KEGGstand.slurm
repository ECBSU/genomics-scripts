#!/bin/bash

#SBATCH -p compute
#SBATCH --output=KEGG_stand_%a.out
#SBATCH -t 0-12
#SBATCH --mem=40G
#SBATCH -c 32
#SBATCH --array=0-0   #CHANGE THIS TO MATCH THE NUMBER OF FASTA FILES 0-0=1 0-2=3 0-6=7 ETC.

#This script runs eggnog on the fasta files in the provided directory, and uses the output to calculate the completion of each KEGG module.
##########################################################################################################################################

#########################
#INPUT VARIABLES, CHANGE AS NEEDED!!!
#########################
input_dir=/path/to/dir  #Directory containing your fasta files
output_dir=/path/to/dir #Directory to which to write the output
extension=.fasta #Extension used in your fasta files

#########################
#Read input
#########################
cd $input_dir
sf=(*$extension)
mkdir $output_dir
cd $output_dir

#########################
#Run Eggnog
#########################
source /home/a/arno-hagenbeek/miniconda3/etc/profile.d/conda.sh #This and the following line initiate the conda env
conda activate /home/a/arno-hagenbeek/miniconda3/envs/eggnog-mapper
emapper.py --data_dir /bucket/HusnikU/Databases/Eggnog_db/Eggnog_Diamond_db/ -i ${input_dir}/${sf[${SLURM_ARRAY_TASK_ID}]} --itype genome --translate -o ${output_dir}/${sf[${SLURM_ARRAY_TASK_ID}]} --cpu 32

#########################
#Invoke black coding magic to run Arno's script of horrible spaghetti code
#########################
python /bucket/HusnikU/Unit_Members/Arno/Scripts/Eggnog_KEGG_completion_checker.py ${output_dir}/${sf[${SLURM_ARRAY_TASK_ID}]}.emapper.annotations
