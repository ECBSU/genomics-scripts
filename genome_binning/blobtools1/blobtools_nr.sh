# Load necessary modules
ml bioinfo-ugrp-modules
ml DB/blastDB/ncbi/2022-07-nt

# Activate conda
source /home/y/yong-phua/conda_ini.sh
conda activate blobtools

# Run BLASTn
/hpcshare/appsunit/HusnikU/ncbi-blast-2.11.0+/bin/blastn -query scaffolds.fasta -db /bucket/HusnikU/Databases/nt/nt -outfmt '6 qseqid staxids bitscore std' \
    -max_target_seqs 1 -max_hsps 1 -evalue 1e-25 -num_threads 64 > scaffolds_vs_blastn.out

# Create BlobTools database
/apps/unit/HusnikU/blobtools/blobtools create -i scaffolds.fasta -y spades -t scaffolds_vs_blastn.out --db /bucket/.deigo/HusnikU/Databases/blobtaxid/nodesDB.txt 

# Generate BlobPlots
/apps/unit/HusnikU/blobtools/blobtools plot -i blobDB.json

/apps/unit/HusnikU/blobtools/blobtools plot -i blobDB.json -r order

# Generate taxonomic tables
/apps/unit/HusnikU/blobtools/blobtools view -i blobDB.json -r phylum -r order -r family
