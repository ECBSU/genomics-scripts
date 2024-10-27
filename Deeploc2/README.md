For DeepLoc2, If you large fasta file, it will take an extremely long time to run.
I will first split the fasta file (split_fasta.py) then do an array slurm script

split_fasta.py -i <input fasta file> -o <output dir> -n <entries per split fasta file>

change #SBATCH --array= to fit the number of entries from split_fasta.py
For more details on array slurm files, please check deigo documentation. 
