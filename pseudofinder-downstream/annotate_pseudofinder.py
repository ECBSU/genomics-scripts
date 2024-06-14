#!/usr/bin/env python

"""Convert a GFF and associated FASTA file into GenBank format.

Usage:
    annotate_pseudofinder.py -i <pseudofinder intact gff file> -p <prokka gff> -o <output gff>
"""

import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--intact", help="gff with intact genes from pseudofinder",required=True)
parser.add_argument("-p", "--prokka", help="prokka gff file with gene annotation",required=True)
parser.add_argument("-o", "--output", help="output gff file",required=True)
args = parser.parse_args()

# Identify intact genes
intact = pd.read_csv(args.intact, sep='\t', comment='#', header=None)
intact_names = intact[8].tolist()

# Obtain lists with gene annotation 
with open(args.prokka, encoding='latin1') as prokka:
    id = []
    filtered_prokka = []
    fasta = []
    fasta_lines = False
    for ln in prokka:
        if fasta_lines:
            fasta.append(ln)
            continue
        if ln.startswith("##sequence-region"):
            id.append(ln)
            continue
        if any(name in ln for name in intact_names):
            filtered_prokka.append(ln)
            continue
        if len(ln.split("\t")) >= 3 and "RNA" in ln.split("\t")[2]:
            filtered_prokka.append(ln)
            continue
        if ln.strip() == "##FASTA":
            fasta_lines = True
            fasta.append(ln)
            continue

# Write file
with open(args.output, 'w',) as o:
    o.write("##gff-version 3\n")
    o.writelines(id)
    o.writelines(filtered_prokka)
    o.writelines(fasta)
