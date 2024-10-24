#!/usr/bin/env python3
"""
Parses the K number from an Eggnog .annotation file. Uses a database of KEGG k term entries to both
find to which (sub)category each k term belongs, and to reconstruct the KEGG hierarchy of nested
categories. Then per (sub)category checks which k terms are present based on the EggNOG output.

Usage: (python) KEGGstand_category_gene_presence.py input.emapper.annotations output_file KEGG_k_term_database

Output:
Writes a tab-delimited file that lists each KEGG category along with its subcategories. After
a tab, the number of genes found for the (sub)category is given, along with the total number of
genes known for this category. 
"""

###################
#Import statements
################### 
import sys
import os

###################
#Functions
################### 

def gen_line_reader(file_path):
    """
    Generator function that allows reading a text file line 
    by line, without reading the full file into memory
    """
    for line in open(file_path, "r"):
        yield line

def leading_dash_counter(string):
    """
    Counts the number of leading dashes in the string. Returns
    this number as a integer
    """
    dash = 0
    for char in string:
        if char == "-":
            dash += 1
        else:
            break
    return dash
    
def completion_tsv_reader(filepath):
    """
    Parser function for reading the completion values in the KEGGstand pipeline output.
    Returns a dictionary where the module is the key, and its completion fraction is the value
    """
    out_dict = {}
    for line in gen_line_reader(filepath):
        if line.startswith("#"):
            continue
        #Combine the KEGG Module code and its name into a single variable
        KEGG_name = "{} {}".format(line.split("\t")[0], line.split("\t")[1])
        completion = float(line.split("\t")[2])
        out_dict[KEGG_name] = completion
    return out_dict

def BRITE_output_parser(file_path, BRITE_list, allbrite):
    """
    Parser function for reading the BRITE/pathway gene counts given in the KEGGstand pipeline output.
    Returns a dictionary where the pathway is the key, and the value is a string: 
    "genes_found/total_genes", which will be parsed in the output function
    """
    out_dict = {}
    #Check if allbrite was given (this means every category will be returned)
    #If not, only return the categories in the BRITE_list
    if not allbrite:
        #To determine which categories are subcategories of interest, keep track of the leading dashes
        cat_leading_dash = 100
        #Initiate a variable that shows that we are processing a subcategory
        in_cat = False
        for line in gen_line_reader(file_path):
            #If there are fewer or equal leading dashes as the categories of interest,
            #then we are not in a subcategory, and the line will be ignored
            if leading_dash_counter(line) <= cat_leading_dash:
                in_cat = False
            #If the line contains a category of interest, the category will be stored. Furthermore,
            #every next line, UNTIL there are fewer or equal dashes, is a subcategory, which will
            #also be processed using the in_cat variable.
            if line.split("\t")[0].strip("-").strip() in BRITE_list:
                out_dict[line.split("\t")[0]] = line.split("\t")[1].strip()
                cat_leading_dash = leading_dash_counter(line)
                in_cat = True
            if in_cat:
                out_dict[line.split("\t")[0]] = line.split("\t")[1].strip()
    #If allbrite is given, simply output every single category
    else:
        for line in gen_line_reader(file_path):
            out_dict[line.split("\t")[0]] = line.split("\t")[1]
    return out_dict
    
def module_output(module_dict, out_file, minimum_completion, avg):
    """
    Output function to write the contents of the module dictionary to a tab-delimited text file.
    Considers the minumum completion and avg variables. Will only output modules that have at least one
    organism with the module above minumum completion. If avg is true, only considers modules where
    the average is above minimum completion
    """
    out = open(out_file, "w")
    
    #Write header line
    out.write("#Module\t")
    for org in module_dict:
        out.write(org + "\t")
    out.write("\n")
    #Write per module completions for modules above minimum completion
    #First iterate over the modules by grabbing the first organism, and iterating over the keys of its nested dictionary,
    #after which the loop is broken       
    for x in module_dict:
        for module in module_dict[x]:
            #Acquire all completion values for module
            completion_list = []
            for org in module_dict:
                if module in module_dict[org]:
                    completion_list.append(float(module_dict[org][module]))
                else:
                    completion_list.append(0)
            #If in average mode, check if the average is above minimum completion. Else check if any value is above minimum
            if avg:
                #If above required average, output to file
                if sum(completion_list)/len(completion_list) >= minimum_completion:
                    out.write(module + "\t")
                    for i in completion_list:
                        out.write(str(i) + "\t")
                    out.write("\n")
            else:
                #If any completion is above minimum, output to file
                output = False
                for i in completion_list:
                    if i >= minimum_completion:
                        output = True
                if output:
                    out.write(module + "\t")
                    for i in completion_list:
                        out.write(str(i) + "\t")
                    out.write("\n")
        break
            
def BRITE_output(brite_dict, out_file):
    """
    Output function to write the contents of the BRITE dictionary to a tab-delimited text file
    """
    out = open(out_file, "w")
    #Write header line
    out.write("#Module\t")
    for org in brite_dict:
        out.write(org + "\t")
    out.write("\n")
    #Write per category gene count. Write total genes known behind the category
    for x in brite_dict:
        for cat in brite_dict[x]:
            total_genes = brite_dict[x][cat].split("/")[1]
            out.write(cat + " ({} genes)\t".format(total_genes))
            for org in brite_dict:
                if cat in brite_dict[org]:
                    out.write(str(brite_dict[org][cat].split("/")[0]) + "\t")
                else:
                    out.write("0" + "\t")
            out.write("\n")
        break
        
            

####################################################################
#MAIN
####################################################################    
if __name__ == "__main__":
    #Acquire input
    if "-i" in sys.argv:
        indir = sys.argv[sys.argv.index("-i") + 1]
    else:
        print("No input specified, specify input directory with '-i'")
        sys.exit()
    if "-o" in sys.argv:
        out_prefix = sys.argv[sys.argv.index("-o") + 1]
    else:
        out_prefix = "KEGGstand_merged"
    #Initiate and check optional variables
    minimum_completion = 0
    all_BRITE = False
    BRITE_list = []
    avg = False
    if "--minimum" in sys.argv:
        minimum_completion = float(sys.argv[sys.argv.index("--minimum") + 1])
    if "--KEGG_cat" in sys.argv:
        KEGG_cat = sys.argv[sys.argv.index("--KEGG_cat") + 1]
        if KEGG_cat == "ALL":
            all_BRITE = True
        elif "," in KEGG_cat:
            for i in KEGG_cat.split(","):
                BRITE_list.append(i.strip())   
        else:
            BRITE_list = [KEGG_cat]
    if "--average" in sys.argv:
        avg = True
                
    #Initiate variables
    module_dict = {}
    BRITE_dict = {}
    #Iterate over the files in the input directory            
    for root, dirs, files in os.walk(indir): 
        for file in files: 
            #Obtain the completion value per module per fasta
            if file.endswith(".emapper.annotations_KEGG_completion.tsv"):
                org_name = file.partition(".emapper")[0]
                modules = completion_tsv_reader(os.path.join(root, file))
                module_dict[org_name] = modules
            #Only read the BRITE and pathway output if required
            if all_BRITE or len(BRITE_list) > 0:
                if file.endswith(".emapper.annotations_pathway_and_BRITE"):
                    org_name = file.partition(".emapper")[0]
                    BRITE_dict[org_name] = BRITE_output_parser(os.path.join(root, file), BRITE_list, all_BRITE)
    #Write output
    module_output(module_dict, out_prefix + "_merged", minimum_completion, avg)       
    
    if all_BRITE or len(BRITE_list) > 0:
        BRITE_output(BRITE_dict, out_prefix + "_merged_BRITE")
