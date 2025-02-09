# Bin identification summary script
## Info
This is a python script that makes a short, (somewhat) human readable summary of the bin identification pipeline. It will read the SSU BLAST results, 
and the MEGAN output per bin, and give this information per bin per sample, in a single text file. 

As this script is a python script, not a slurm script, you will need to run it using python. 
Usage:
Python Pipeline_total_taxonomy_summary.py PIPELINE/OUTPUT/DIRECTORY/PATH OUTPUT_NAME.tsv

## Input
Unlike the previous scripts, this script requires no hardcoded input, and instead takes its input as above:
Python Pipeline_total_taxonomy_summary.py PIPELINE/OUTPUT/DIRECTORY/PATH OUTPUT_NAME.tsv

Where PIPELINE/OUTPUT/DIRECTORY/PATH is the path to the identification pipeline's output directory and 
OUTPUT_NAME.tsv is the name (and path) of the desired output file.

## Output
This script generates an tab-delimited text file.

The file lists each bin per sample, and lists per bin:

The best (based on percentage identity times match length) 5 BLAST results for the barrnap output in the format:
BLAST_match identity    match_length
Per taxonomic level, the number of LCA matches found by MEGAN.

## Summary interpretation
This is how I tend to judge the results of the binning:

A good bin tends to have one or several full length 16S or 18S hit(s). Whether the found 16S or 18S makes sense should be informed by the type of sample. 
Multiple 16S/18S hits are expected and fine if the organisms in the hits are related. If there are multiple unrelated 16S or 18S genes, this is likely an indication
of contamination of the bin.

MEGAN results are usually quite biased, with a tendency to give many matches in group with many reference genomes. If the MEGAN results give many hits for the
same organisms as the found 16S/18S, this is again a good sign. However, given the biases in MEGAN, a lack of relevant matches does not mean the bin is poor, especially 
in organisms that are poorly studied. 

I tend to use the MEGAN results to roughly judge contamination. For example, a large number of human hits in a non-animal bin would be a bad sign. Similarly,
a large number of Propionibacteria or similar skin contamintants in unrelated bacterial bins are typically a sign of poor binning.

For bins that seem to contain organisms of interest, I typically assess its completion with BUSCO/CheckM and compare bin length, GC content and coding density to the those
known from related organisms. If bin length (given the degree of completion), GC content, and coding density are in line with related organisms, I move on to annotation using
for example bakta,eggnog or BLASTKoala. 
