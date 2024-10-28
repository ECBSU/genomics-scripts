MEGAN - analysis for metagenomes or in my case for identifying HGT candidates in a transcriptome.
#diamond blastp or blastx can be used like this for megan anaysis 

```
input=/flash/HusnikU/Phua_2/sino_transcriptome/SPAdes_cln1234/transcript_process/transcripts.fasta.transdecoder.pep

/apps/unit/HusnikU/diamond blastp -p 64 -f 100 --db /flash/HusnikU/databases/NR_data_diamond/data.dmnd --sensitive -q $input -o transcripts_pep.daa
```

I run the remaining steps on my own laptop, you will need the MEGAN programme,link can be found here: https://software-ab.cs.uni-tuebingen.de/download/megan6/welcome.html

You will also need this database megan-map-Feb2022.db.zip, to link the .daa output into MEGAN visualisation format.

#In the megan GUI, there is an option to meganise a file, you will need the fasta file and database file 
