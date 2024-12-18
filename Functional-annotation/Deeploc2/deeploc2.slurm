#!/bin/bash

#SBATCH -p gpu
#SBATCH --output=deeploc2.out
#SBATCH -t 1-0
#SBATCH --mem=46G
#SBATCH -c 4
#SBATCH --gres=gpu:1
#SBATCH --job-name=deeploc

#INPUT VARIABLES !!CHANGE AS NEEDED!!
input_fasta=/path/to/fasta
out_dir=/path/to/dir
conda_ex=/path/to/activate/your/conda #OR use mine /work/HusnikU/phua/conda_ini.sh 

#Activate modules and conda enviroment
source ${conda_ex}
conda activate /home/y/yong-phua/miniforge3/envs/deeploc
ml cuda

#Run Deeploc2.1
deeploc2 -f $input_fasta -o $out_dir -m Accurate -d cuda
