#!/usr/bin/env python3
"""
Parses checkm2 and gtdb-tk output and summarizes them together, listing the bins found to match each taxonomy,
sorted by their checkm2 quality. 

EG
Gtdbtk_bin_extracter.py checkm2_quality_report gtdb_tk_bac_summary gtdb_tk_arc_summary outfile
"""
###################
#Import statements
###################  

import sys
import os
import shutil

###################
#File handling functions
###################  
def gen_line_reader(file_path):
    for line in open(file_path, "r"):
        yield line

checkm2_file = sys.argv[1]
gtdbtk_bac_file = sys.argv[2]
gtdbtk_arc_file = sys.argv[3]
outfile = sys.argv[4]

#find the MAG stats and make a tuple containing the MAG name, completion, contamination and size. 
checkm_list = []
for line in gen_line_reader(checkm2_file):        
    if line.startswith("Name\tCompleteness"):
        continue
    if line.strip():
        s = line.split("\t")
        checkm_list.append((s[0], s[1], s[2], s[8]))

#Read all gtdb-tk entries (bacteria and archea). Make a dictionary to group classified MAGs per taxonomy. 
#The assigned taxonomy is the key, and a list of tuples (each containing stats of a MAG) as the value
tax_dict = {}
for line in gen_line_reader(gtdbtk_bac_file):
    if line.startswith("user_genome"):
        continue
    if line.strip():
        s = line.split("\t")
        if s[1].startswith("Unclassified"):
            continue
        if s[1].strip() in tax_dict:
            for i in checkm_list:
                if i[0] == s[0].strip():
                    tax_dict[s[1]].append(i)
        else:
            tax_dict[s[1]] = []
            for i in checkm_list:
                if i[0] == s[0].strip():
                    tax_dict[s[1]].append(i)
for line in gen_line_reader(gtdbtk_arc_file):
    if line.startswith("user_genome"):
        continue
    if line.strip():
        s = line.split("\t")
        if s[1].startswith("Unclassified"):
            continue
        if s[1].strip() in tax_dict:
            for i in checkm_list:
                if i[0] == s[0].strip():
                    tax_dict[s[1]].append(i)
        else:
            tax_dict[s[1]] = []
            for i in checkm_list:
                if i[0] == s[0].strip():
                    tax_dict[s[1]].append(i)

#Output the found MAGs per taxonomy and their relevant stats. While doing so, order the tuple list in the dictionaries by completion, then contamination and then size.
out = open(outfile, "w")
for tax in tax_dict:
    out.write(tax + "\n")
    out.write("MAG\tCompletion\tContamination\tSize(bp)\n")
    lst = tax_dict[tax]
    new_list = sorted(lst, reverse = True, key=lambda element: (element[3]))
    new_list = sorted(new_list, key=lambda element: (element[2]))
    new_list = sorted(new_list, reverse = True, key=lambda element: (element[1]))
    for tup in new_list:
        out.write("{}\t{}\t{}\t{}\n".format(tup[0], tup[1], tup[2], tup[3]))
        
    