#!/bin/bash
#SBATCH --partition compute
#SBATCH --time 4-0
#SBATCH --mem=200G
#SBATCH --cpus-per-task=32
#SBATCH --job-name=ANGSD_admixture
#SBATCH --output=/flash/HusnikU/vivienne/Palythoa/Jobs/202411/1106_ANGSD_admixture.log
#SBATCH --error=/flash/HusnikU/vivienne/Palythoa/Jobs/202411/1106_ANGSD_admixture.err
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


#### 03.admixture.GL1.slurm
#generate output file in BEAGLE format to estimate admixture proportions from genotype likelihoods using NGSadmix (later)
#GL 1 (samtools model)

admixture_GL1() {

        FILTERS="-minMapQ 30 -minQ 20 -uniqueOnly 1 -remove_bads 1 -baq 1 -only_proper_pairs 1 -trim 0 -C 50 -maxDepth 500"
        ${ANGSD} -b "${all_bam_list}" -ref ${reference} ${FILTERS} -GL 1 -doCounts 1 -doMajorMinor 4 -doMaf 1 -skipTriallelic 1 -SNP_pval 1e-6 -doGlf 2 -out ./all_admixture

        return 0
}


# run functions
set_software_paths
admixture_GL1
