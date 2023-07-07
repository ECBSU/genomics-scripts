# Multigene phylogeny for plastid genome

Code by Dewi, uses Geneious (locally) and splace, mafft and iqtree on Deigo

## Create fasta files with Geneious

1- import all plastid genomes
2- Tools/Extract annotations with the following settings:
    2.1 name does not contain fragment
    2.2 type is gene
    2.3 names does not contain orf
    2.4 name does not contain rpo
    2.5 name does not contain |

3- export all the sequences list to fasta file (Export/Multiple files/FASTA sequences-alignment)

## install Splace
``` bash
git clone https://github.com/reinator/splace.git

```

## rename fasta files headers
for splace to work, fasta headers should be at the >gene_accession  or >gene formats
By default, geneious produces headers with >accession - gene

``` bash
cd /flash/HusnikU/Dewi/Phylogeny/Bacillariales/

for f in *\ *; do mv "$f" "${f//  annotations/}"; done  # removes "annotations" at the end of the file
for f in *\ *; do mv "$f" "${f// annotations/}"; done # removes "annotations" at the end of the file
for f in *\ *; do mv "$f" "${f// /_}"; done # removes "annotations" at the end of the file

for i in *.fasta; do sed 's|\ gene.*||;s/.*- />/' ${i} > rename/${i}; done # reformat the fasta file header

# sed 's|\ gene.*||;s/.*- />/' *.fasta > rename/*.fasta
```
## rename fasta files headers
Now that files are prepared, run splace (requires mafft for alignment)

``` bash
mkdir /flash/HusnikU/Dewi/Phylogeny/Bacillariales/splace

find /flash/HusnikU/Dewi/Phylogeny/Bacillariales/rename/*.fasta -type f > /flash/HusnikU/Dewi/Phylogeny/Bacillariales/splace/Bac_files.txt ### creates a text file with the list of fasta files to align

srun -p compute -t 0-8 --mem=32G -c 64 --pty bash  ###Run interactive code 

ml load mafft
cd /flash/HusnikU/Dewi/Phylogeny/Bacillariales/splace

python /home/d/dewi-langlet/splace/Dockerfile/SPLACE_v2.py Bac_files.txt 64 Bacillariales #Runs Splace to the Bacillariales output using 64 cores ()


sed 's/.*\//>/' Bacillariales_result/final_concat_aligned_Bacillariales.fasta > Bacillariales_result/final_concat_aligned_Bacillariales_simple.fasta ### removes path name from the sequences headers
```

## Run IQ-tree
Run IQ-tree on the sequence alignment produced by splace
``` bash
ml bioinfo-ugrp-modules Other/iqtree

mkdir /flash/HusnikU/Dewi/Phylogeny/Bacillariales/splace/iqtree
cd /flash/HusnikU/Dewi/Phylogeny/Bacillariales/splace/iqtree

# it might be a good idea to request more ressources to run iqtree here
iqtree -s /flash/HusnikU/Dewi/Phylogeny/Bacillariales/splace/Bacillariales_result/final_concat_aligned_Bacillariales_simple.fasta -m TEST -bb 1000 -alrt 1000


```

