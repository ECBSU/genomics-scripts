#!/usr/bin/env python3
"""
Renames the protein names in the .faa prokka output to match the sample name  
"""
import sys
import os

#Generator function that reads file line by line, without reading the full file into memory
def gen_line_reader(file_path):
    for line in open(file_path, "r"):
        yield line

input_dir = sys.argv[1]
for subdir, dirs, files in os.walk(input_dir): #iterate over all directories and subdirectories
    for file in files: #iterate over all files
        #Only use .faa files
        if file.endswith(".faa"):
            #Rewrite files with samplenames instead of protein names
            temp = file + "_temp"
            out = open(os.path.join(subdir,temp), "w")
            sample = file.partition(".faa")[0]
            print(sample)
            for line in gen_line_reader(os.path.join(subdir,file)):
                if line.startswith(">"):
                    end = line.rpartition("_")[2]
                    out.write(">{}_{}".format(sample, end))
                else: 
                    out.write(line)
            #Delete the old .faa file and rename the new one 
            os.remove(os.path.join(subdir,file))
            os.rename(os.path.join(subdir,temp), os.path.join(subdir,file))