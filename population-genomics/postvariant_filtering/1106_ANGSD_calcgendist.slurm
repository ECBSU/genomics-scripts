#!/bin/bash
#SBATCH --partition compute
#SBATCH --time 4-0
#SBATCH --mem=200G
#SBATCH --cpus-per-task=128
#SBATCH --job-name=ANGSD_calcgendist
#SBATCH --output=/flash/HusnikU/vivienne/Palythoa/Jobs/202411/1106_ANGSD_calcgendist.log
#SBATCH --error=/flash/HusnikU/vivienne/Palythoa/Jobs/202411/1106_ANGSD_calcgendist.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=vivienne.gao@oist.jp

set -euo pipefail
set -x # enable debugging

#===============================================================================================#
#-------------Softwares--------------#
ANGSD_DIR='/home/v/vivienne-gao/angsd'
#------------------------------------#

#---------------Inputs---------------#
RUN_DIR='/flash/HusnikU/vivienne/Palythoa/Jobs/202411/ANGSD'
INPUT_DIR='/flash/HusnikU/vivienne/Palythoa/Jobs/202410/bowtie2'
reference='/flash/HusnikU/vivienne/Palythoa/Palythoa_Final_assembly.fasta'
#------------------------------------#
#===============================================================================================#

# set software paths!
set_software_paths() {

        ANGSD='/home/v/vivienne-gao/angsd/angsd'

        mkdir -p ${RUN_DIR}
        cd ${RUN_DIR} || exit

        # create a temporary directory for the list of paths to all BAM files
        all_bam_list=$(mktemp)

        for bam_file in "${INPUT_DIR}"/*.bam; do
                echo "${bam_file}" >> "${all_bam_list}"
        done


	return 0
}


#### 02.call.genotypes.GL1.slurm
# call genotypes using the samtools model (GL 1) - genotype likelihoods are the likelihood of the data given a genotype

call_genotypes() {

        # set filters (quality same as above), discard sites where there is data from >10 ind
        FILTERS="-minMapQ 30 -minQ 20 -uniqueOnly 1 -remove_bads 1 -baq 1 -only_proper_pairs 1 -trim 0 -C 50 -maxDepth 500 -minInd 10 -setMinDepth -setMaxDepth"

        # find major/minor alleles, calculate per site freq, only biallelic sites, remove sites with pvalue>0.01, call genotypes, calculate posterior prob, 
        ${ANGSD} -b "${all_bam_list}" -ref ${reference} ${FILTERS} -GL 1 -doCounts 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-3 -doGeno 32 -doPost 1 -out "./all_genotypes"

        return 0
}

#### 02.genetic.dist.GL1.slurm
#calculate genetic distances using the samtools model (GL 1)

calculate_gen_dist() {

        FILTERS_3="-minMapQ 30 -minQ 20 -uniqueOnly 1 -remove_bads 1 -baq 1 -only_proper_pairs 1 -trim 0 -C 50 -maxDepth 500"
        ${ANGSD} -b "${all_bam_list}" -ref ${reference} ${FILTERS_3} -GL 1 -doCounts 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-3 -doGeno 8 -doPost 1 -out ./all_gene_dist

        return 0
}


set_software_paths
calculate_gen_dist
