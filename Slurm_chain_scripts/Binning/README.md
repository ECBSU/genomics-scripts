# Binning chain script
The various scripts forming a pipeline for binning. The pipeline runs various analyses split over 3 (slurm) jobs which are performed successively:

### Job 1
1. Trimming (fastp)
2. Assembly (spades)
### Job 2
3. Binning (Metabat2, CONCOCT, Semibin, Rosella, and refinement by DASTool)
### Job 3
4. Quality check (Checkm2)
5. Taxonomic assignment (GTDB-TK)
6. Results summarizing (python)

The jobs will be executed in series, with the next one running after the previous finishes. Additionally, the script also allows for array jobs, allowing multiple samples to be analyzed in parallel. 

## How to run
First ensure you have appropriate input. The pipeline takes a directory as input. This directory should be filled with paired (short-read) fastq files. The script by default finds fastq files based on the standard extensions received from the SQC at OIST:
"R1_001.fastq.gz" for forward reads and "R2_001.fastq.gz" for reverse reads. 

With the input directory prepared copy the scripts to a directory writeable by deigo compute, either your home folder or one of your directories on flash. Note that the logs of the script will be outputted in the same location as where the scripts are stored and executed.

Next, open the Pipeline_wrapper_script.slurm file in a text editor (in nano for example). In line 11 input the path to the directory containing the fastq files. In line 13 put the path of the directory where you want the binning output. 
Lastly, if you want to run an array job (like if you have multiple fastq pairs), change the filenumber of line 16. IMPORTANT: This number is an index, which usually start counting from 0. In short, put in 0 for one sample, and 99 for a hundred samples etc. 

Finally, move to the directory containing the scripts, and execute the wrapper script using the command: "./Pipeline_wrapper_script.slurm"

