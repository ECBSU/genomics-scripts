#!/usr/bin/env python3
"""
Reads all the files in a directory finding the trimmed files, and concatenates the sequences based on the sequence names  
"""
import sys
import os

#Generator function that reads file line by line, without reading the full file into memory
def gen_line_reader(file_path):
    for line in open(file_path, "r"):
        yield line

input_dir = sys.argv[1]
output_file = sys.argv[2]
out_dict = {}

#First catalog all genes/proteins
for subdir, dirs, files in os.walk(input_dir): #iterate over all directories and subdirectories
    for file in files: #iterate over all files
        if file.endswith("_trimmed"):
            for line in gen_line_reader(os.path.join(subdir,file)):
                if line.startswith(">"):
                    sample_name = line.strip(">").rpartition("_")[0]
                else:
                    if line:
                        line = line.strip()
                        if sample_name in out_dict:
                            out_dict[sample_name].append(line)
                        else:
                            out_dict[sample_name] = [line]

outfile = open(output_file, "w")
for sample in out_dict:
    outfile.write(">" + sample + "\n")
    outfile.write("".join(out_dict[sample]) + "\n")