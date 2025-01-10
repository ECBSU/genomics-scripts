#!/bin/bash
#SBATCH --partition compute
#SBATCH --time 4-0
#SBATCH --mem=200G
#SBATCH --cpus-per-task 128
#SBATCH --nodes 1
#SBATCH --job-name=1031_bcftools_palythoa
#SBATCH --output=/flash/HusnikU/vivienne/Palythoa/Jobs/202410/1031_bcftools_palythoa.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=vivienne.gao@oist.jp

set -euo pipefail


#===============================================================================================#

#-------------Softwares--------------#
BCFTOOLS='/apps/unit/HusnikU/bcftools'
#------------------------------------#

#---------------Inputs---------------#
RUN_DIR='/flash/HusnikU/vivienne/Palythoa/Jobs/202410'
INPUT_DIR='/flash/HusnikU/vivienne/Palythoa/Jobs/202410/bowtie2'
OUTPUT_DIR='/flash/HusnikU/vivienne/Palythoa/Jobs/202410/bcftools'
reference='/flash/HusnikU/vivienne/Palythoa/Palythoa_Final_assembly.fasta'
#------------------------------------#

#===============================================================================================#

# set software paths
set_software_paths() {

        MPILEUP="${BCFTOOLS}/bcftools mpileup"
        CALL="${BCFTOOLS}/bcftools call"
        return 0
}

#main function

main() {

        # make directory if it does not exist
        mkdir -p "${OUTPUT_DIR}"
        cd "${RUN_DIR}" || exit

        # set software paths
        set_software_paths

        echo "----------------------------------------"
        echo "Running bcftools..."
        echo "----------------------------------------"

        # Loop over each .sorted.bam file in the input directory
        for bam_file in "${INPUT_DIR}"/*.sorted.bam; do # starts a for loop to analyze all bam files in the directory

                # Get the base filename without path and extension
                base_name=$(basename "${bam_file}" .sorted.bam)

                # Step 1: Convert BAM to VCF
                echo "Processing ${base_name}"

                if [ ! -f "$OUTPUT_DIR/${base_name}_variants.vcf" ]; then # if ! runs command if files do not exist, [] is for evaluating test expressions (one file)
                        echo 'Running bcftools mpileup and bcftools call  for' ${base_name}
                        ${MPILEUP} -Ou -f "${reference}" "${bam_file}" | ${CALL} -mv -Ov -o "$OUTPUT_DIR/${base_name}_variants.vcf" # bcftools mpileup output piped to bcftools call
                        # -Ou specifies uncompressed BCF output (ideal for piping to other commands), -mv detect snps and indels (-m) and output only variant sites (-v)

                else
                    	echo "${base_name} vcf file already exists. Skipping vcf file generation."
                fi

                echo "Finished processing ${base_name}."
        done

	echo "----------------------------"
        echo "ALL SAMPLES COMPLETED"
        echo "----------------------------"

        return 0

}

# run main function
main
