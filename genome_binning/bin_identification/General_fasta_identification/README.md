# Fasta identification pipeline
## Info
This pipeline identifies .fasta files by annotating with prokka, extracting the 16S sequences from the prokka annotation,
and BLASTing the 16S genes against the ncbi nt database.

It also runs BLASTn and MEGAN on each fasta, and creates a summary of the MEGAN output (so the megan application doesnt
have to be opened, and the .rma doesn't have to be imported)

## Input
This pipeline only requires a directory containing .fasta files to identify


The input is hardcoded, so needs to be adjusted in the script itself.
Please change the following lines to suit your case:

Line 8: #SBATCH --array=x-x <-- Change the array size to fit your number of samples (0-0 = 1 sample, 0-1 = 2 samples etc)

Line 22: fasta_dir=/Path/To/X <-- Change the path to the path of your directory containing the pipeline output

Line 25: id_out_dir=/Path/To/X <-- Change the path to your desired output folder (probably somewhere on /flash)

If all variables are correct, just run the script using "sbatch Fasta_identification_pipeline.slurm"

## Output
Per fasta creates a directory containing all generated bins, a directory containing the .ffn prokka outputs for each bin,
a directory containing the nt BLAST and MEGAN results. In addition to a .rma file used in the MEGAN application, a python script is used to generate a MEGAN LCA summary which lists the most abundant hits per taxonomy level. 

Within the output directory creates (again per sample) a .fasta file containing all 16S rRNA genes found among the bins, and a .blast file which contains the BLAST results for all the found 16S rRNA genes. 
