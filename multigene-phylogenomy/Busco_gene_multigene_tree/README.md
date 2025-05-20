# Multigene tree reconstruction based on shared BUSCO genes
## Note
Akito rewrote this script to be neater and more modular. I recommend using that instead: 
https://github.com/ASUQ/OIST_Scripts/tree/main/Python_codes/Busco_multigene_tree

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
A fraction of 0 will result in analyzing every gene (at the cost of generating many gaps), and a fraction of 1 will result in only analyzing genes shared by all organisms (often resulting in few genes). To more easily find the optimal fraction, multiple can be provided to the script in a comma delimited format (e.g. 0.77,0.8,0.9,0.999999). 

If you wish, the script also allows you to modify the MAFFT, trimal, and iqtree commands to match your needs. These commands contain two {} values. Python substitutes these with the appropriate input and output values. Please leave these parts unchanged.  

Note on time allotment: Deigo slurm scripts require a hardcoded time allotment. Longer is better, but will also make it harder for the server to schedule your job (resulting in longer queue times). The default of the script is 12h, which is sufficient to run a single fraction for 50 organisms. Adding fractions will aproximately multiply the required time.

Any input is hardcoded, so needs to be adjusted in the script itself.
Please change the following lines to suit your case:

Line 5: #SBATCH -t 0-12 <-- Change to the time allotment of the script (the left number denotes days, the right denotes hours). 12h tends to be sufficient for a single run/fraction for 50 organisms.

Line 20: input_dir=/Path/To/X <-- Change the path to the path of your directory containing the BUSCO output directories

Line 21: output_dir=/Path/To/X <-- Change the path to your desired output folder (probably somewhere on /flash)

Line 22: fraction=number <-- Change number to desired fraction (see above for explanation). 

Line 23-25 <- Tool commands that can be changed to match your desired commands. 

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


