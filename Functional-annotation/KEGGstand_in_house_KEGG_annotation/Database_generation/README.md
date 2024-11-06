# Python scripts, used to generate the KEGG databases used by KEGGstand
## Info
The KEGGstand scripts compare the EggNOG-generated K terms against KEGG resources. To obtain these, KEGG API downloads are performed using
the biopython KEGG REST package. The two scripts here generate the KEGG module definition and the KEGG gene entry databases used by the 
KEGGstand scripts. Running these scripts allows for updating existing databases to the current version of KEGG, or to generate new databases. 

Running the module database generation will take several minutes, generating a database of approx. 6kb. The KEGG gene entry database
is significantly larger, and running the corresponding script will take several hours. The generated database will be approx. 27Mb.

## Important notes
Both scripts query the KEGG server. KEGG monitors their traffic and will break connection if excessive traffic is generated. I thus recommend
to not combine running these scripts with other KEGG-related activities. 

The scripts include several measures to minimize and handle disconnects. However, if KEGG disallows connection for a prolonged time (10 consecutive minutes), the script will exit. 
The script will store its current progress, and if called again using the same output file, will check it's progress and continue where it left off. This functionality also allows
the scripts to update existing databases to the current version of KEGG.

## Input
The scripts both take only a single argument, which is the name/path of the database that it will generate.

Example command:
```
python KEGG_kterm_db_generate.py /db_dir/KEGGstand_kterm_db
```

# Current issues and workaround
The module completion calculation script was written under the assumption that commas could only exist within brackets in a KEGG definition. This turned out to not be the case, and module M00134 and M00745
do have commas outside of brackets. This causes the script to incorrectly read these definitions. 

**To get these modules to work correctly, please change them manually in the database:**

M00134: K01476,K27545 K01581 -> (K01476,K27545) K01581

M00745: (K18072 K18073),(K07644 K07665),K18297 K18093 -> ((K18072 K18073),(K07644 K07665),K18297) K18093
