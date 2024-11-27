# Python scripts to predict KEGG module completion and gene presence/absence
## Info 
These scripts parse KEGG k-terms from EggNOG output. They compare these k-terms against databases obtained from KEGG API to estimate KEGG module completion, and gene counts per BRITE category or pathway.

## Module checker script
### Input
The script takes 3 inputs: 
- The .emapper.annotations output generated by EggNOG
- The desired name/path of the script's prefix for generation output files
- The KEGG module database, which can be generated with the corresponding script from:  
https://github.com/ECBSU/genomics-scripts/tree/main/Functional-annotation/KEGG_pathways/in-house_KEGG_annotation/Database_generation

The script is run with the following command:
```
(python) KEGGstand_module_checker.py input.emapper.annotations output_prefix KEGG_module_database
```
### Output
Writes 2 tab-delimited text files: 
- *prefix*_KEGG_completion.tsv: lists the found completion for all KEGG modules. 
- *prefix*_KEGG_complete_modules.tsv: reports only the modules found to be 100% complete. 

The format is the same for the two files, and looks as follows:
```
#Entry	Name	Highest_completion	Most_complete_pathway	Non-essential_genes_found
M00001	Glycolysis (Embden-Meyerhof pathway), glucose => pyruvate	0.9	K01810,K00850,K01624,K01803,K00134,K00927,K15633,K01689,K00873	None
M00002	Glycolysis, core module involving three-carbon compounds	1.0	K01803,K00134,K00927,K15633,K01689,K00873	None
M00003	Gluconeogenesis, oxaloacetate => fructose-6P	0.875	K01610,K01689,K15633,K00927,K00134,K01803,K01624	None
M00004	Pentose phosphate pathway (Pentose phosphate cycle)	0.7142857142857143	K01783,K01807,K00615,K00616,K01810	None
M00005	PRPP biosynthesis, ribose 5P => PRPP	1.0	K00948	None
```
Where the columns are:
- First column: The module entry code
- Second column: Full name of the module
- Third column: highest completion found for this module as a fraction. 
- Fourth column: A comma-delimited list of the genes that comprised the most complete pathway found. 
- Fifth column: A list of module genes found to be present, but not essential to the functioning of the module.

## BRITE checker script
### Input
The script takes 3 inputs: 
- The .emapper.annotations output generated by EggNOG
- The desired name/path of the script output file
- The KEGG gene database, which can be generated with the corresponding script from:  
https://github.com/ECBSU/genomics-scripts/tree/main/Functional-annotation/KEGG_pathways/in-house_KEGG_annotation/Database_generation

The script is run with the following command:
```
(python) KEGGstand_BRITE_checker.py input.emapper.annotations output_file KEGG_k_term_database
```
### Output
Writes a single tab-delimited text file with a specified name in the following format:
```
Ribosome biogenesis [BR:ko03009]	33/544
-Prokaryotic type	32/95
--rRNA modification factors	12/42
---23S rRNA modification factors	6/29
---16S rRNA modification factors	6/14
--Other maturation factors	6/17
--rRNA transcription factors	4/4
--ribosomal protein modification factors	2/10
--Sno RNPs	0/6
---Box C/D snoRNPs	0/3
---Box H/ACA snoRNPs	0/3
--GTPases	4/7
--rRNA processing factors	3/5
--rRNA helicases	1/4
```

The first column of the file denotes the pathway or BRITE category. The second column shows the number of genes found as a fraction of the total number of genes that exist in the category: found/total. 
The categories are listed in hierarchical fashion, with the dashes preceding the entry name denoting subcategories. 