# Small python script to extract KEGG k terms from EggNOG output
## Info
Among the EggNOG output are the K terms associated with the genes. This script extracts these K terms and outputs them in a format similar 
to blastKOALA. This way KEGG pathways can be assessed using KEGG reconstruct (https://www.genome.jp/kegg/mapper/reconstruct.html) but without
individually processing fasta files on the blastKOALA webserver.

Note: EggNOG and blastKOALA results often differ. Expect EggNOG to yield slighlty more pathways. On average EggNOG finds 3% more k terms, which often results in at least several more completed modules. 
In my experience, these additional k terms often complete almost complete pathways, suggesting they might be true positives. However, there is no way to confirm this. 

## Input
The input of the script is the prefix.emapper.annotations file generated by EggNOG provided in the following format:

"python Eggnog_KEGG_KO_extracter.py input.emapper.annotations"

## Output
By default the output is a text file named based on input (prefix.emapper.annotations_only_KEGG), which contains the KEGG K terms in a format
similar to what is received when downloading the blastKOALA results (a tab-delimited file with column 1 being the gene, and column 2 being the K term).

Note: EggNOG sometimes gives multiple K terms per gene. In these cases BlastKOALA only gives 1 K term. To prevent changes in format clashing with downstream
applications, this script only outputs the first K term given by EggNOG per gene. 
