## Various commands used to call SNPs and generate genotype frequency files

#### 1031_bcftools_palythoa.slurm
- this script generates individual vcf files for each bam file
- sample names in VCF file are included as entire path to file so it looks kind of ugly - can change this with the bcftools reheader command
- no filtering at all - output all multiallelic variant sites


Tools used:
- BCFtools
- VCFtools
- ANGSD
- bgzip
- tabix




ANGSD pipeline

Link to ANGSD github repository: https://github.com/ANGSD/angsd

ANGSD website: https://www.popgen.dk/angsd/index.php/Main_Page#Overview
