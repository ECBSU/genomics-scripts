## Pre-variant filtering tools and alignment quality checks

#### 1101_samtools_coverage_palythoa.slurm
- calculates average coverage for each sample

#### 1106_ANGSD_qc.slurm
- generates 6 files with various quality metrics
- counts.gz: sequencing depth at every position in the genome for each individual
- depthGlobal, depthSample: Column1 in the .depthSample,.depthGlobal contains the number of sites with sequencing depth of 0. Column2 is the number of sites with a sequencing depth of 1, etc. The .depthSample contains depth per sample. Line one corresponds to individual 1. Column2 corresponds to individual 2 etc.
- pos.gz: sum of reads covering a sites for all individuals
- .qs: count of q-score values within minQ maxQ thresholds


Tools used
- SAMtools
- VCFtools
- ANGSD
