#!/bin/bash
#SBATCH -p compute
#SBATCH -t 4-0
#SBATCH --mem=64G
#SBATCH -c 64
#SBATCH -n 1
#SBATCH --job-name=deeploc_array
#SBATCH --output=/path/to/logs/%A_%a.out
#SBATCH --array=XX-XX%XX

files=(path/to/folder/of/fasta/*)
echo "slurm: " ${files[${SLURM_ARRAY_TASK_ID}]}

source /home/y/yong-phua/conda_ini.sh 
conda activate deeploc
deeploc2 -f ${files[${SLURM_ARRAY_TASK_ID}]} -o /flash/HusnikU/Phua_2/sino_transcriptome/deeploc_full -m Accurate
conda deactivate