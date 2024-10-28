You first have to download the AMOEBAE code.
I placed this in my flash folder. You can put it in yours.

```
git clone https://github.com/laelbarlow/amoebae.git
```

After which you should do the edits to the config files. refer to https://github.com/laelbarlow/amoebae/blob/master/documentation/workflow_protocol.md
- genomes.csv
- queries.csv
- reference_db_list.txt

*Direct download from NCBI via AMOEBAE tends to fail. It is alot smoother to get the files and place them in resources/local_db_files and resources/local_query_files

Once you have that set up, run the setup script (AMOEBAE_setup.slurm). #make sure to change line 10 to your AMOEBAE directory.

Next is the manual curation of the hits against your database. 
When you are done with the manual curation, save the file in the results folder as "Ref_seqs_1_manual_selections.csv"

Once you have that set up, run the analysis script (AMOEBAE_analysis.slurm). #make sure to change line 10 to your AMOEBAE directory.
