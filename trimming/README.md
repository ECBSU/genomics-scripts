# Adapter trimming with `fastp`

### Info

[Fastp](https://github.com/OpenGene/fastp) - tool for preprocessing for FastQ files allows to do adapter trimming, filtering reads by quality and many more.


Most important parameters are:

`-i`                 input read 1

`-o`                 output read 1

`-I`                 input read 2

`-O`                 output read 2

`--adapter_fasta`    provide a `fasta` file with sequences of Illumina adapters and oligo-dT primers
#### Example usage
`fastp -w 16 -i r1 -o out_r1 -I r2 -O out_r2 --json fastp.json --html fastp.html --adapter_fasta smartseq2_oligos.fasta`

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) - quality control of reads. 
FastQC Output from multiple datasets can also be visualised in [MultiQC](https://multiqc.info/).
#### Example usage
`fastqc read1.fastq.gz read2.fastq.gz read3.fastq.gz readN.fastq.gz -o /output/directory/path/`

### Requirements and installation
`fastp fastqc`
You can create a conda environment with these packages
`conda create -n trim "fastp==0.23.2" "fastqc==0.12.1"`
