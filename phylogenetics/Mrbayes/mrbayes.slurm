#!/bin/bash
#SBATCH -p compute
#SBATCH -t 4-0
#SBATCH --mem=64G
#SBATCH -c 16
#SBATCH -n 1
#SBATCH --job-name=bayesproteins
#SBATCH --output=bayesproteins.log

nex=/flash/HusnikU/Petra/iqtree_Dioscoreales/conc_20g_Diosc.nex

ml load MrBayes.mpi
mpirun mb /$nex > bayesproteins.log

