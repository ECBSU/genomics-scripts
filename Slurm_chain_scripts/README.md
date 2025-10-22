# Chain scripts for slurm
These are a collection of bash scripts for submitting chain jobs to Deigo. Their main advantage is consistency and ease of use, with a single script submitting multiple jobs to deigo, which are then performed in series. 

To run any of the scripts, copy them over to your Deigo home folder (or any other folder writeable from the compute partition). Then write your input into the appropriate lines of the wrapper script, and then execute the wrapper script. 

## Contained scripts
### Binning
fastp -> spades -> metabat2 -> concoct -> semibin2 -> rosella -> Dastool -> checkm2 -> gtdb-tk
### Phylogeny
prokka -> orthofinder -> MAFFT -> trimal -> iqtree2
