#!/bin/bash

##########################
#Input your variables here
#           ||
#           ||
#           ||
#           \/
##########################
#Put the path of the directory containing your fastq files here after the "=" symbol 
fastq_dir=/bucket/HusnikU/Unit_Members/Arno/Gustavo_bobtail_analysis/Central_analysis/All_squid_fastq_temp
#Put the desired path of the output directory here after the "=" symbol . Most likely /flash/HusnikU/your_name/this_job
outdir=/flash/HusnikU/Arno/chain_pipe_test/
#Last index of the array
#This equals your number of samples/fastq pairs minus one (one sample = 0, three samples = 2 etc.)
filenumber=0

###########################
#Other variables
###########################
#Conda executable. Unless you have specific requirements, changing this shouldnt be necessary
conda_ex=/home/a/arno-hagenbeek/miniforge3/etc/profile.d/conda.sh
############################
mkdir $outdir
#Run the script running fastp and spades
#Use the "jobstring" construction to create a chain job, automatically running the next job after the previous finishes
jobstring=$(sbatch --array=0-${filenumber} fastp_into_spades.slurm $fastq_dir $outdir)
jobid=${jobstring##* }

#Run the binning pipeline
jobstring=$(sbatch --dependency=afterany:${jobid} --array=0-${filenumber} binning_pipeline.slurm $fastq_dir $outdir $conda_ex)
jobid=${jobstring##* }

#Finally run checkm2 and gtdb-tk to assess bin quality, and find bins of interest.
#This job also includes a python script that summarizes the MAG stats and groups bins by assigned taxonomy 
sbatch --dependency=afterany:${jobid} --array=0-${filenumber} bin_stats.slurm $outdir $conda_ex
