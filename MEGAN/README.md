MEGAN - analysis for metagenomes or in my case for identifying HGT candidates in a transcriptome.
#diamond blastp or blastx can be used like this for megan anaysis 

```
input=/flash/HusnikU/Phua_2/sino_transcriptome/SPAdes_cln1234/transcript_process/transcripts.fasta.transdecoder.pep

/apps/unit/HusnikU/diamond blastp -p 64 -f 100 --db /flash/HusnikU/databases/NR_data_diamond/data.dmnd --sensitive -q $input -o diamond_blast_out.daa

/apps/unit/HusnikU/megan/tools/daa-meganizer -i diamond_blast_out.daa -mdb /bucket/.deigo/HusnikU/Databases/megandb/megan-map-Feb2022.db -t 64

/apps/unit/HusnikU/megan/tools/blast2lca -i 'diamond_blast_out.daa' -o '/flash/HusnikU/Phua_2/legend/transcripts_taxa.txt' -f DAA -mdb /bucket/.deigo/HusnikU/Databases/megandb/megan-map-Feb2022.db
```

For visualisation, I use the MEGAN programme installed on my laptop. when you open the DAA file, you should get the reads linked to taxon visualised in a phylogenetic tree. 
