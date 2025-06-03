
/hpcshare/appsunit/HusnikU/diamond blastx --db "/flash/HusnikU/databases/swissprot/uniprot_sprot.fasta" -q "/flash/HusnikU/Phua_2/Sino_genome/SPAdes_all_seq/scaffolds.fasta" --outfmt 6 --sensitive --max-target-seqs 1 --evalue 1e-25 -p 64 > scaffolds_v_sprot.diablastx
sed 's/|/\t/g' "/flash/HusnikU/Phua_2/Sino_genome/SPAdes_all_seq/scaffolds_v_sprot.diablastx" | awk '{$2=$4=""; print $0}' OFS="\t" > scaffolds_v_sprot.diablastx
export PYTHONPATH="/apps/unit/HusnikU/pysam/lib64/python3.6/site-packages/:$PYTHONPATH"
/apps/unit/HusnikU/blobtools/blobtools taxify -f scaffolds_v_sprot.diablastx -m "/flash/HusnikU/databases/reference_proteomes.taxid_map" -s 0 -t 2
/apps/unit/HusnikU/blobtools/blobtools create -i scaffolds.fasta -y spades -t scaffolds_v_sprot.diablastx2.taxified.out --db /bucket/.deigo/HusnikU/Databases/blobtaxid/nodesDB.txt 
/apps/unit/HusnikU/blobtools/blobtools plot -i blobDB.json
/apps/unit/HusnikU/blobtools/blobtools plot -i blobDB.json -r order
/apps/unit/HusnikU/blobtools/blobtools view -i blobDB.json -r phylum -r order -r family 
