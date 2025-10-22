# Phylogeny chain script
The various scripts forming a pipeline for binning. The pipeline runs various analyses split over 3 (slurm) jobs which are performed successively:

### Job 1
1. Gene calling (prokka)
2. Ortholog finding (orthofinder)
### Job 2
3. Ortholog alignment (MAFFT)
4. Ortholog trimming (trimal)
### Job 3
4. Phylogenetic tree reconstruction (iqtree)

The jobs will be executed in series, with the next one running after the previous finishes.

## How to run
First ensure you have appropriate input. The pipeline takes a directory as input. This directory should be filled with (prokaryote) genome files in fasta format. The extension of the files can be specified in the wrapper script

With the input directory prepared copy the scripts in this repository (everything but the README.md) to a directory writeable by deigo compute (either your home folder or one of your directories on flash). Note that the logs of the script will be written in the same location as where the scripts are stored and executed.

Next, open the Phylogeny_wrapper_script.slurm file in a text editor (in nano for example). In line 11 input the path to the directory containing the genome files. In line 13 put the path of the directory where you want the pipeline output. You can set the appropriate extension of your genome fasta files on line 15.

Finally, move to the directory containing the scripts, and execute the wrapper script using the command: "./Phylogeny_wrapper_script.slurm"

## How to use the output
There are four main outputs for this script. \
\

Firstly is are the gene sequences generate by prokka, found in: \
"/specified_output_folder/prokka/". \
Secondly is the output of an orthofinder run, found in:\
"/specified_output_folder/orthofinder/"\
Thirdly, are the aligned and trimmed sequences of the orthologs found in:\
"/specified_output_folder/trimal_alignment/"\
Lastly, is the phylogenetic tree found in:\
"/specified_output_folder/tree/"\

## Notes
+ This pipeline is not suited for eukaryotes, only including prokaryote gene calling
+ Most scripts have been extensively tested in isolation, but the chain script is mostly untested. Please let me know if you encounter any bugs.
+ For ease of use, most tool settings are hardcoded, but after downloading, these settings can be changed freely in their relevant scripts if you have different use cases. 


