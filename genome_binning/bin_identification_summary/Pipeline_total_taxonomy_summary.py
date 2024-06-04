#!/usr/bin/env python3
"""
Reads the identification pipeline outputs and writes a summary of both MEGAN and SSU hits per sample and per bin

Usage: (Python) Pipeline_total_taxonomy_summary.py PIPELINE_OUTPUT_DIRECTORY OUTPUT.tsv
"""
import sys
import os

#Functions
def gen_line_reader(file_path):
    """
    Generator function that reads a file line by line, without reading the file into memory
    Input:
    file_path(string) = path to file to read
    Output:
    Iterator, where each item is a single line of the provided file
    """
    for line in open(file_path, "r"):
        yield line
#Main
if __name__ == "__main__":
    sample_dir = {} #Initialize dictionary for holding found data
    for root, dirs, files in os.walk(sys.argv[1]): #iterate over all directories and subdirectories
        for file in files: #iterate over all files
            if file.endswith("_SSU_blast"): #if file is SSU BLAST results
                sample = root.rpartition("/")[2].partition("_")[0] #find sample name from file path
                if sample not in sample_dir: #add sample to output dictionary
                    sample_dir[sample] = {}
                binname = file.partition("_SSU")[0] #Find bin name from file name
                if binname not in sample_dir[sample]: 
                    sample_dir[sample][binname] = {}
                SSU_hits_list = [] #initialize list to store SSU hits in output
                for line in gen_line_reader(os.path.join(root,file)): #read SSU BLAST output
                    splt = line.split("\t")  
                    if "uncultured" in splt[2] or "Uncultured" in splt[2]: #Ignore items with "uncultured" since theyre usually not informative
                        continue
                    else:
                        hit_length = float(splt[4].strip())
                        identity = float(splt[3].strip())
                        match = splt[2].strip()
                        qualityscore = hit_length*identity #Give the results a quality score (length*identity) to preferentially list better longer hits
                        SSU_hits_list.append((qualityscore,identity, hit_length, match)) 
                if SSU_hits_list != []:
                    sample_dir[sample][binname]["SSU"] = sorted(SSU_hits_list, reverse=True)[0:4] #Sort hits based on quality score to preferentially show better hits
            if file.endswith("_lca_summary"): #if file is MEGAN results
                sample = root.rpartition("/")[2].partition("_")[0] #find sample name from file path
                binname = file.partition("_lca")[0]  #Find bin name from file name
                if sample not in sample_dir: #Add sample and bin to dict if not there yet. 
                    sample_dir[sample] = {}
                if binname not in sample_dir[sample]:
                    sample_dir[sample][binname] = {}
                tax = "dom" #Initialize tax variable
                for line in gen_line_reader(os.path.join(root,file)): #Iterate over output lines, changing the tax variable depending on the found taxonomic level
                    if line.strip()=="":
                        continue
                    if line.startswith("---------"):
                        continue
                    if "Top domains" in line:
                        tax = "dom"
                    elif "Top kingdoms" in line:
                        tax = "kin"
                    elif "Top phyla" in line:
                        tax = "phy"
                    elif "Top classes" in line:
                        tax = "cla"
                    elif "Top orders" in line:
                        tax = "ord"
                    elif "Top families" in line:
                        tax = "fam"
                    elif "Top genera" in line:
                        tax = "gen"
                    elif "Top species" in line:
                        tax = "spe"
                    else:
                        if tax not in sample_dir[sample][binname]: #else is a result line, if taxonomic level has not been added yet, add it. 
                            sample_dir[sample][binname][tax] = []
                        sample_dir[sample][binname][tax].append((line.split("\t")[0],line.split("\t")[1].strip()))
                
    outfile = open(sys.argv[2], "w")  #Write output
    for sample in sorted(sample_dir):
        for binname in sorted(sample_dir[sample]):
            if "phy" not in sample_dir[sample][binname] and "SSU" not in sample_dir[sample][binname]:
                continue
            outfile.write(">>>{}\t{}\n".format(sample,binname))
            if "SSU" in sample_dir[sample][binname]:
                if sample_dir[sample][binname]["SSU"] != []:
                    outfile.write("Best_SSU_hits (name,identity,length):\n++++++++++++++++\n")
                    for tup in sample_dir[sample][binname]["SSU"]:
                        outfile.write("{}\t{}\t{}\n".format(tup[3], tup[1], tup[2]))
            else: 
                outfile.write("No SSU hits found\n")
                
            outfile.write("\nTop MEGAN results:\n++++++++++++++++\n")
            if "phy" in sample_dir[sample][binname]:
                outfile.write("Top phyla:\n-------------------\n")
                for tup in sample_dir[sample][binname]["phy"]:
                    outfile.write("{}\t{}\n".format(tup[0],tup[1]))
            if "cla" in sample_dir[sample][binname]:
                outfile.write("Top classes:\n-------------------\n")
                for tup in sample_dir[sample][binname]["cla"]:
                    outfile.write("{}\t{}\n".format(tup[0],tup[1]))
            if "ord" in sample_dir[sample][binname]:
                outfile.write("Top orders:\n-------------------\n")
                for tup in sample_dir[sample][binname]["ord"]:
                    outfile.write("{}\t{}\n".format(tup[0],tup[1]))
            if "fam" in sample_dir[sample][binname]:
                outfile.write("Top families:\n-------------------\n")
                for tup in sample_dir[sample][binname]["fam"]:
                    outfile.write("{}\t{}\n".format(tup[0],tup[1]))
            if "gen" in sample_dir[sample][binname]:
                outfile.write("Top genera:\n-------------------\n")
                for tup in sample_dir[sample][binname]["gen"]:
                    outfile.write("{}\t{}\n".format(tup[0],tup[1]))
            if "spe" in sample_dir[sample][binname]:
                outfile.write("Top species:\n-------------------\n")
                for tup in sample_dir[sample][binname]["spe"]:
                    outfile.write("{}\t{}\n".format(tup[0],tup[1]))
            outfile.write("\n")
            
                    