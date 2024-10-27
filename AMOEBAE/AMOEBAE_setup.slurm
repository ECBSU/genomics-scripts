#!/bin/bash
#SBATCH -p compute
#SBATCH -t 2-0
#SBATCH --mem=2G
#SBATCH -c 32
#SBATCH -n 1
#SBATCH --job-name=amoebae
#SBATCH --output=amoebae_setup.log

cd /path/to/your/amoebae_dir

#activate conda
source /home/y/yong-phua/conda_ini.sh
conda activate /hpcshare/appsunit/HusnikU/mamba/snakemake

snakemake get_ref_seqs -j 100 --use-conda
