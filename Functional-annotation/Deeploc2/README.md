On Saion, 
Change the input variables in the slurm script
Then run deeploc2.slurm 

---

On Deigo if you large fasta file, it will take an extremely long time to run.
You will also need to change line 22 `cuda` -> `cpu` 

I recommend to first split the fasta file (split_fasta.py) then do an array slurm script

```
source /home/y/yong-phua/conda_ini.sh 
conda activate deeploc
split_fasta.py -i <input fasta file> -o <output dir> -n <entries per split fasta file>

```

change #SBATCH --array= to fit the number of entries from split_fasta.py
For more details on array slurm files, please check deigo documentation. 
