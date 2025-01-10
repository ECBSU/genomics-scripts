#!/bin/bash
#SBATCH --partition compute
#SBATCH --time 4-0
#SBATCH --mem=200G
#SBATCH --cpus-per-task=128
#SBATCH --job-name=ANGSD_qc
#SBATCH --output=/flash/HusnikU/vivienne/Palythoa/Jobs/202411/1106_ANGSD_qc.log
#SBATCH --error=/flash/HusnikU/vivienne/Palythoa/Jobs/202411/1106_ANGSD_qc.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=vivienne.gao@oist.jp

set -euo pipefail
set -x # enable debugging

#===============================================================================================#

#-------------Softwares--------------#
ANGSD_DIR='/home/v/vivienne-gao/angsd'
SAMTOOLS_DIR='/apps/unit/HusnikU/samtools-1.15.1' 
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
        SAMTOOLS_INDEX='/apps/unit/HusnikU/samtools-1.15.1/samtools index'

        mkdir -p ${RUN_DIR}
        cd ${RUN_DIR} || exit

        all_bam_list=$(mktemp)

        for bam_file in "${INPUT_DIR}"/*.sorted.bam; do
                echo "${bam_file}" >> ${all_bam_list}
        done

	return 0
}


#### 01.qc.filtering.slurm
#assessing base qualities and coverage depth

# remove reads with low mapping quality and sites with low quality or covered by few reads (low depth)
qc_filtering() {

        # set filters min mapping quality, min base quality score. removes reads with multiple best hits, remove sup. align/duplicates/fails, perform BAQ, only include read pairs that both mapped correctly. no trimming, reduce effect of reads with excessive mismatches, s$
        FILTERS="-minMapQ 30 -minQ 20 -uniqueOnly 1 -remove_bads 1 -baq 1 -only_proper_pairs 1 -trim 0 -C 50 -maxDepth 500"

        # outputs qscores, sequencing depth, count of ACGT, seq depth per sample
        ${ANGSD} -b "${all_bam_list}" -ref ${reference} ${FILTERS} -doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2 -out "${RUN_DIR}/ANGSD_qc"

        return 0
}


set_software_paths

qc_filtering
