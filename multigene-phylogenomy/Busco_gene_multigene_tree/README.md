# Multigene tree reconstruction based on shared BUSCO genes
## Info
This pipeline will build a multigene tree based on BUSCO genes shared by provided organisms. The script will read BUSCO output of multiple organisms and make 
multifasta files for shared genes. These multifasta's are then aligned using MAFFT, and trimmed using trimal, and concatenated. Concatenated output is then 
reconstructed into a phylogenetic tree with iqtree (using modelfinder). 

## Input
The pipeline requires BUSCO output directories (with unmodified structure) as input. All directories to be analyzed need to be placed into a single directory. 
(e.g. input_dir/organism1_busco, input_dir/organism2_busco), which is provided as input to the script. Note that all BUSCO output should be generated with the same BUSCO database to prevent weird results. 

The script also requires one or multiple "fractions" as input. BUSCOs are not necessarily shared by all organisms, particularly with high numbers of input organisms. 
Therefore, only considering genes shared by all provided organisms often results in trees built based on very few genes, or outright failure. To prevent this, 
a for which the gene needs to be present can be provided. 
For example, 0.8 analyzes genes shared by at least 80% of provided organisms. For the organisms missing this gene, it will be a gap in the final output. 
A fraction of 0 will result in analyzing every gene (at the cost of generating many gaps), and a fraction of 1 will result in only analyzing genes shared by all organisms (often resulting in few genes). To more easily find the optimal fraction, multiple can be provided to the script in a comma delimited format (e.g. 0.77,0.8,0.9,0.999999)

Any input is hardcoded, so needs to be adjusted in the script itself.
Please change the following lines to suit your case:

Line 5: #SBATCH -t 0-12 <-- Change to the time allotment of the script (the left number denotes days, the right denotes hours). Times should be longer, especially for multiple fractions or large inputs. 12h tends to be sufficient for a single run for 50 organisms. This should be multiplied per fraction. 

Line 20: input_dir=/Path/To/X <-- Change the path to the path of your directory containing the BUSCO output directories

Line 21: output_dir=/Path/To/X <-- Change the path to your desired output folder (probably somewhere on /flash)

Line 22: fraction=number <-- Change number to desired fraction (see above for explanation. 

If all variables are correct, just run the script using "sbatch Busco_multigene_tree.slurm"

## Output
The output will all be generated in the specified folder. 

Multifasta output:
Output includes protein sequence multifasta files of all BUSCO genes found in the outputs, where each line in the fasta represents an input organism (if the gene was found for it). If a gene was included in of any of the provided fractions, MAFFT alignments and trimal trimmed alignments of the corresponding multifasta's will also be generated. 

Concatenated output:
For each provided fraction, a multifasta file of the concatenated alignment of the matching fraction is generated. These are the input of iqtree for the tree reconstruction

Iqtree output:
For each fraction, iqtree output is generated. The .treefile can be visualized into a phylogenetic tree using figtree or other software. Since the pipeline doesnt specify a model, it makes iqtree run
its modelfinder function. To find the optimal model iqtree found, check the .log file. 

## Notes on tool use
For now, the tool commands are hardcoded in the python script. Please let me know if you would like a version with altered commands.

MAFFT is run using the following command:
mafft --globalpair --maxiterate 1000 --thread 64 input > output

Trimal is run using the command:
/home/a/arno-hagenbeek/trimal/trimal/source/trimal -in input -out output -automated1

Iqtree is run using the command:
iqtree2 -s input --prefix output -bb 1000 -T 64


