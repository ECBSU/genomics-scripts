#!/usr/bin/env python3
"""
Reads all the EggNOG .annotations file and extracts the KEGG ko number. Outputs a file similar to the downloaded
outputs from BLASTkoala. 

Usage: (python) Eggnog_KEGG_KO_extracter.py input.emapper.annotations output_name (optional)
"""
import sys
import os

#Generator function that reads file line by line, without reading the full file into memory
def gen_line_reader(file_path):
    for line in open(file_path, "r"):
        yield line

#Obtain inputs
eggnog_file = sys.argv[1]
if len(sys.argv) > 2:
    outname = sys.argv[2]
else:
    outname = eggnog_file + "_only_KEGG"
out = open(outname, "w")

#Read over EggNOG file, and write the KEGG KO terms per gene to the output
for line in gen_line_reader(eggnog_file):
    if line.startswith("#"):
        continue
    gene = line.split("\t")[0]
    ko = line.split("\t")[11]
    #!!!!A comma means there is multiple ko terms. BLASTkoala appears to only save the first one, so will do the same here
    if ko == "-":
        ko = ""
    elif "," in ko: 
        ko = ko.partition(",")[0]
    ko = ko.replace("ko:", "")
    out.write("{}\t{}\n".format(gene, ko))
    
                    
