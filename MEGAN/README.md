MEGAN - analysis for metagenomes or in my case for identifying HGT candidates in a transcriptome.
#diamond blastp for megan anaysis 
#the database is --db /flash/HusnikU/databases/NR_data_diamond/data.dmnd
```
infile=/flash/HusnikU/Phua_2/sino_transcriptome/SPAdes_cln1234/transcript_process/transcripts.fasta.transdecoder.pep
/apps/unit/HusnikU/diamond blastp -p 64 -f 100 --db /flash/HusnikU/databases/NR_data_diamond/data.dmnd --sensitive -q $infile -o transcripts_pep.daa
```

#After running diamond, you will need to 'meganise' the file for Megan to give you the detailed taxonomical info of each read
#that uses this database megan-map-Feb2022.db.zip , link can be found here: https://software-ab.cs.uni-tuebingen.de/download/megan6/welcome.html
#In the megan GUI, there is an option to meganise a file, you will need the fasta file and database file 
