#!/bin/bash
#SBATCH -p compute
#SBATCH -t 0-2
#SBATCH --mem=5G
#SBATCH -c 4
#SBATCH -n 1
#SBATCH --job-name=blob2
#SBATCH --output=blob2_setup.log

work_dir=/path/to/work_dir

#activate conda environment 
source /home/y/yong-phua/conda_ini.sh 
conda activate btk

#create database
blobtools create --fasta $work_dir/scaffolds.fasta $work_dir/blob2/

#add hits, diamond outputs needs to be formatted, (taxyfy breaks format)
blobtools add --hits $work_dir/scaffolds_v_sprot.out --hits $work_dir/scaffolds_v_nt.out --taxrule bestsumorder --taxdump /home/y/yong-phua/blobtaxid/taxdump $work_dir/blob2

#add coverage
blobtools add --cov $work_dir/sorted.mapped.bam --threads 4 $work_dir/blob2
