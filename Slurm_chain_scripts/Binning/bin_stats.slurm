#!/bin/bash

#SBATCH -p compute
#SBATCH --output=bin_stats%a.out
#SBATCH -t 1-0
#SBATCH --mem=128G
#SBATCH -c 32

#move all bins of the same sample to the same folder
bin_dir=(${1}/binning/*_concoct)
mkdir ${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_combined_bins
mv ${bin_dir[${SLURM_ARRAY_TASK_ID}]}/*.fa ${1}/${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_combined_bins/
mv ${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_semibin/*.fa ${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_combined_bins/
mv ${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_metabat/*.fa ${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_combined_bins/
mv ${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_rosella/*.fna ${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_combined_bins/
mv ${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_DASTool/*.fa ${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_combined_bins/
rm -r ${bin_dir[${SLURM_ARRAY_TASK_ID}]}/
rm -r ${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_semibin/
rm -r ${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_metabat/
rm -r ${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_rosella/
rm -r ${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_DASTool/

#Change the Rosella bin extensions to match
for i in ${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_combined_bins/*.fna; do mv ${i} ${i%.*}.fa; done

#Run checkm2 on the bins 
source $2
conda activate '/bucket/HusnikU/Conda-envs/checkm2'
/apps/unit/HusnikU/checkm2/bin/checkm2 predict --threads 32 --input ${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_combined_bins/  --output-directory ${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_bin_stats/checkm2 --database_path '/bucket/HusnikU/Databases/checkm2/CheckM2_database/uniref100.KO.1.dmnd' -x fa --force

#Run gtdb-tk on the bins
conda activate /bucket/HusnikU/Conda-envs/gtdb-tk/
gtdbtk classify_wf --genome_dir ${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_combined_bins/ --out_dir ${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_bin_stats/gtdbtk --skip_ani_screen --extension fa --cpus 32
python ./Pipeline_bin_stats.py ${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_bin_stats/checkm2/quality_report.tsv ${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_bin_stats/gtdbtk/gtdbtk.bac120.summary.tsv ${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_bin_stats/gtdbtk/gtdbtk.ar53.summary.tsv ${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_bin_stats/bin_stats_summary.txt

#Copy the best bin over for each taxonomy
python ./Pipeline_bin_copier ${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_bin_stats/bin_stats_summary.txt ${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_combined_bins/ ${bin_dir[${SLURM_ARRAY_TASK_ID}]%%_concoct}_best_bins/