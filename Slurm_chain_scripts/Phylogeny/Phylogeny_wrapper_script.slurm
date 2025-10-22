#!/bin/bash

##########################
#Input your variables here
#           ||
#           ||
#           ||
#           \/
##########################
#Put the path of the directory containing your genomes here after the "=" symbol 
genome_dir=/bucket/HusnikU/Unit_Members/Arno/Gustavo_bobtail_analysis/Central_analysis/All_squid_fastq_temp
#Put the desired path of the output directory here after the "=" symbol . Most likely /flash/HusnikU/your_name/this_job
outdir=/flash/HusnikU/Arno/chain_pipe_test/
#Put your fasta extension here (excluding the dot)
extension=fa

###########################
#Other variables
###########################
#Conda executable. Unless you have specific requirements, changing this shouldnt be necessary
conda_ex=/home/a/arno-hagenbeek/miniforge3/etc/profile.d/conda.sh
prokka_path=/bucket/HusnikU/Conda-envs/prokka/
############################
mkdir $outdir
#Run the script for prokka and orthofinder
#Use the "jobstring" construction to create a chain job, automatically running the next job after the previous finishes
jobstring=$(sbatch orthofinder_tree_pipe1.slurm $genome_dir $outdir $extension $conda_ex $prokka_path)
jobid=${jobstring##* }

#Run the script for mafft and trimal
jobstring=$(sbatch --dependency=afterany:${jobid} orthofinder_tree_pipe2.slurm $outdir)
jobid=${jobstring##* }

#Run the script for iqtree 
sbatch --dependency=afterany:${jobid} orthofinder_tree_pipe2.slurm $outdir
