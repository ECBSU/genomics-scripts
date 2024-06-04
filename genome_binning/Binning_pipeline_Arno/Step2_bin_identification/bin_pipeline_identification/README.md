# Bin identification pipeline
## Info
This pipeline identifies the .fasta files created by the binning pipeline by annotating with barrnap, extracting the SSU and 
BLASTing the 16S genes against the ncbi nt database.

It also runs BLASTn and MEGAN on each fasta, and creates a summary of the MEGAN output (so the megan application doesnt
have to be opened, and the .rma doesn't have to be imported)

As with the binning pipeline, this is designed to work directly on the OIST Deigo server.

## Input
This pipeline is designed to work on the output generated by the DASTool-based binning pipeline. For general fasta's check the 
bin_identification_pipeline_general_fasta directory

The input is hardcoded, so needs to be adjusted in the script itself.
Please change the following lines to suit your case:

Line 8: #SBATCH --array=x-x <-- Change the array size to fit your number of samples (0-0 = 1 sample, 0-1 = 2 samples etc)

Line 22: bin_pipeline_out_dir=/Path/To/X <-- Change the path to the path of your directory containing the pipeline output

Line 25: id_out_dir=/Path/To/X <-- Change the path to your desired output folder (probably somewhere on /flash)

If all variables are correct, just run the script using "sbatch Bin_identification_pipeline.slurm"

## Output
Per sample creates a directory containing all generated bins, a directory containing the total SSU annotation of barrnap, as well as the BLAST results. Additionally, creates a directory containing the nt BLAST (of the entire bin) and MEGAN results. In addition to a .rma file used in the MEGAN application, a python script is used to generate a MEGAN LCA summary which lists the most abundant hits per taxonomy level. 