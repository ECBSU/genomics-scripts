# Binning pipeline using Metabat2, CONCOCT, Semibin, Rosella and DASTool
## Info
This pipeline will align provided reads against the provided assemblies, and use the provided assemblies
and the created .bam files to perform binning. Binning will be performed using Metabat2, CONCOCT, Semibin, and Rosella,
and subsequently refined using DASTool.

This pipeline is a .slurm script that (after giving your input) can be run directly on the OIST Deigo server 
using "sbatch Bin_identification_pipeline.slurm"

## Input
The pipeline requires reads as paired .fastq files, as well as assembled fasta files.
Both need to be placed in their own directories, and named similarly (or at least in the same alphabetical order) to ensure they occupy the same index in their respective arrays.


The input is hardcoded, so needs to be adjusted in the script itself.
Please change the following lines to suit your case:

Line 8: #SBATCH --array=x-x <-- Change the array size to fit your number of samples (0-0 = 1 sample, 0-1 = 2 samples etc)

Line 20: assembly_dir=/Path/To/X <-- Change the path to the path of your directory containing the assemblies

Line 21: rawreads_dir=/Path/To/X <-- Change the path to the path of your directory containing the reads

Line 23: out_dir=/Path/To/Output <-- Change the path to your desired output folder (probably somewhere on /flash)

Line 25: conda_ex=/Path/To/Conda <-- Change the path to the path of your conda installation executable (if you have no conda, you can give Yong Heng's conda, 
which he kindly made public: /home/y/yong-phua/conda_ini.sh)

Line 31: assemblies=(\*.extension) <-- Change the extension to whatever the extension of your assemblies files are

Line 36: reads1=(/${rawreads_dir}/\*.extension) <-- Change the extension to whatever the extension of your read files are

Line 36: reads2=(/${rawreads_dir}/\*.extension) <-- Change the extension to whatever the extension of your read files are


If all variables are correct, just run the script using "sbatch Binning_pipeline.slurm"

## Output
The output will all be generated in the specified folder. Due to the use of multiple tools, many files will be generated. However, the most important ones are the bin.fasta files, which are located in the AssemblyName_BinningSoftware directories. 

Note: I am aware that CONCOCT generates a huge amount of loop errors. The tool functions correctly regardless and since I am not sure how to fix it, please just ignore.  
