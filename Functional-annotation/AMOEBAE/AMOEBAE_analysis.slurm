#!/bin/bash
#SBATCH -p compute
#SBATCH -t 2-0
#SBATCH --mem=2G
#SBATCH -c 32
#SBATCH -n 1
#SBATCH --job-name=amoebae
#SBATCH --output=amoebae_analysis.log

cd /path/to/your/amoebae

#activate conda
source /home/y/yong-phua/conda_ini.sh
conda activate /hpcshare/appsunit/HusnikU/mamba/snakemake

#copy results files to config files
cp results/Ref_seqs_1_manual_selections.csv config/Ref_seqs_1_manual_selections.csv 
cp results/coulson_plot_organization.csv config/coulson_plot_organization.csv
cp results/db_list.txt config/output_plot_row_order.txt

snakemake -j 100 --use-conda

#results may fail to plot
echo 'retrying plot results'

python3 amoebae plot results/fwd_srchs_1_rev_srch_1_interp_with_ali_col_nonredun.csv --complex_info config/coulson_plot_organization.csv --row_order config/output_plot_row_order.txt --out_pdf results/plot.pdf