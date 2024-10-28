
#diamond blastx using database from EukprotV3 querying scaffolde from spades.  
/hpcshare/appsunit/HusnikU/diamond blastx --db "/flash/HusnikU/databases/EukProt_v3/proteins/databse.dmnd" -q "/flash/HusnikU/Phua_2/Sino_genome/SPAdes_all_seq/scaffolds.fasta" --outfmt 6 --sensitive --max-target-seqs 1 --evalue 1e-25 -p 64 > scaffolds_v_eukprot.out