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

def leading_space_counter(string):
    """
    Counts the number of leading spaces in the string. Returns
    this number as a integer
    """
    spaces = 0
    for char in string:
        if char == " ":
            spaces += 1
        else:
            break
    return spaces
    
def parse_KEGG_kterm_db(k_term_db):
    """
    Reads provided file and creates a list of all KEGG entries
    within. 
    """
    #Initiate variables
    K_term_dict = {}
    hierarchy_list = []
    entry = False
    read = False
    temp_list = []
    #Iterate over each line in the database
    for line in gen_line_reader(k_term_db):
        #Parse out the K term from the entry line
        if line.startswith("ENTRY"):
            entry = line.split()[1]
            read = False
            if temp_list != []:
                hierarchy_list.append(temp_list)
                temp_list = []
        #Only read once the read variable is true, which happens only in the BRITE section
        if read:
            #Ignore the lines that starts with the entry itself
            if line.strip().startswith(entry):
                continue
            #Per entry, create a list of tuples (temp_list) that is added to the "hierarchy_list". Each tuple
            #will be (number of leading spaces, KEGG category). This information will be used downstream
            #to reconstruct the KEGG hierarchy
            temp_list.append((leading_space_counter(line), line.strip()))
            #Each line in the BRITE section represents a KEGG category, and it will be stored as a key the KEGG 
            #dictionary (K_term_dict). The values are lists of all the k terms found to belong to each category
            if line.strip() not in K_term_dict:
                K_term_dict[line.strip()] = [entry]
            else:
                K_term_dict[line.strip()].append(entry)
        if line.startswith("BRITE"):
            read = True
    
    if temp_list != []:
        hierarchy_list.append(temp_list)
    return K_term_dict, hierarchy_list

def branching_dict_helper(in_dict, in_list):
    """
    Helper function to create a branching dictionary. 
    
    Note:
    I am aware there exist more elegant recursive solutions. However, I do not fully understand their logic,
    and thus chose to instead implement a more direct solution which will function exactly as I expect. At the time 
    of writing (18-10-24), the highest level of nestedness in the KEGG entries is 5. This function can deal with 
    a level of nestedness of 10, which is unlikely to ever occur. Therefore, while this solution is sub-optimal, 
    it should be quite robust. 
    """
    i = 0
    while i < len(in_list): #Highest level of nestedness is 5, so will do 7 levels to be sure
        if i == 0 and in_list[i] not in in_dict:
            in_dict[in_list[i]] = {}
        elif i == 1 and in_list[i] not in in_dict[in_list[0]]:
            in_dict[in_list[0]][in_list[i]] = {}
        elif i == 2 and in_list[i] not in in_dict[in_list[0]][in_list[1]]:
            in_dict[in_list[0]][in_list[1]][in_list[i]] = {}
        elif i == 3 and in_list[i] not in in_dict[in_list[0]][in_list[1]][in_list[2]]:
            in_dict[in_list[0]][in_list[1]][in_list[2]][in_list[i]] = {}
        elif i == 4 and in_list[i] not in in_dict[in_list[0]][in_list[1]][in_list[2]][in_list[3]]:
            in_dict[in_list[0]][in_list[1]][in_list[2]][in_list[3]][in_list[i]] = {}
        elif i == 5 and in_list[i] not in in_dict[in_list[0]][in_list[1]][in_list[2]][in_list[3]][in_list[4]]:
            in_dict[in_list[0]][in_list[1]][in_list[2]][in_list[3]][in_list[4]][in_list[i]] = {}
        elif i == 6 and in_list[i] not in in_dict[in_list[0]][in_list[1]][in_list[2]][in_list[3]][in_list[4]][in_list[5]]:
            in_dict[in_list[0]][in_list[1]][in_list[2]][in_list[3]][in_list[4]][in_list[5]][in_list[i]] = {}
        elif i == 7 and in_list[i] not in in_dict[in_list[0]][in_list[1]][in_list[2]][in_list[3]][in_list[4]][in_list[5]][in_list[6]]:
            in_dict[in_list[0]][in_list[1]][in_list[2]][in_list[3]][in_list[4]][in_list[5]][in_list[6]][in_list[i]] = {}
        elif i == 8 and in_list[i] not in in_dict[in_list[0]][in_list[1]][in_list[2]][in_list[3]][in_list[4]][in_list[5]][in_list[6]][in_list[7]]:
            in_dict[in_list[0]][in_list[1]][in_list[2]][in_list[3]][in_list[4]][in_list[5]][in_list[6]][in_list[7]][in_list[i]] = {}
        elif i == 9 and in_list[i] not in in_dict[in_list[0]][in_list[1]][in_list[2]][in_list[3]][in_list[4]][in_list[5]][in_list[6]][in_list[7]][in_list[8]]:
            in_dict[in_list[0]][in_list[1]][in_list[2]][in_list[3]][in_list[4]][in_list[5]][in_list[6]][in_list[7]][in_list[8]][in_list[i]] = {}
        elif i == 10 and in_list[i] not in in_dict[in_list[0]][in_list[1]][in_list[2]][in_list[3]][in_list[4]][in_list[5]][in_list[6]][in_list[7]][in_list[8]][in_list[9]]:
            in_dict[in_list[0]][in_list[1]][in_list[2]][in_list[3]][in_list[4]][in_list[5]][in_list[6]][in_list[7]][in_list[8]][in_list[9]][in_list[i]] = {}
        i += 1
    return in_dict
    
def reconstruct_KEGG_hierarchy(hierarchy_list):
    """
    Reconstructs the hierarchy list generated by "parse_KEGG_kterm_db" into a branching dictionary 
    correctly representing the hierarchical structures of categories and subcategories. 
    
    The output is a dictionary with multiple levels of dictionaries within dictionary. The key of each dictionary
    represents a KEGG category while the values represent the subcategories. 
    """
    hierarchy_dict = {}
    #Iterate over each item in the hierarchy list, which each represent the hierarchical categories found for a 
    #k term
    for entry in hierarchy_list:
        #Make temporary list to store potential nestedness
        temp_list = []
        for tup in entry:
            if temp_list != []:
                #If the number of leading spaces is equal or lesser than previous categories, than this is a new nested group
                #It'll get a new entry at the lowest level of nestedness
                if tup[0] <= temp_list[0][0]:
                    temp_list = []
                    if tup not in hierarchy_dict:
                        hierarchy_dict[tup] = {}
                    temp_list.append(tup)
                #If the number of leading spaces is not equal or lesser than all previous categories, then it is part of the nested
                #group of previous tuples. Check at which level of nestedness by comparing leading spaces against all tuples in the nested list.
                #A new list is made to replace temp_list, since all nested tuples below the relevant level can be discarded.
                else:
                    new_temp_list = []
                    for item in temp_list:
                        if item[0] < tup[0]:
                            new_temp_list.append(item)
                    temp_list = new_temp_list[:]
                    #Add the current category to the end of the temp_list, which is it's correct level of nestedness
                    temp_list.append(tup)
                    #Make a list of just the category names (without the number of leading spaces) and add them to the dictionary
                    hierarchy_dict = branching_dict_helper(hierarchy_dict, temp_list)
            else:
                if tup not in hierarchy_dict:
                    hierarchy_dict[tup] = {}
                temp_list.append(tup)
    return hierarchy_dict

def eggnog_parser(eggnog_path):
    """
    Parses an eggnog.annotations output file. Returns a list of all the found k terms.
    """
    out_list = []
    for line in gen_line_reader(eggnog_path):
        if line.startswith("#"):
            continue
        ko = line.split("\t")[11]
        #!!!!A comma means there is multiple ko terms. BLASTkoala appears to only save the first one.
        #This script will include both ko terms
        if ko == "-":
            ko = ""
        elif "," in ko: 
            for i in ko.split(","):
                i = i.replace("ko:", "")
                out_list.append(i)             
        else:
            ko = ko.replace("ko:", "")
            out_list.append(ko)
    return out_list   

def recursive_dict_iteration(dictionary):
    """
    Iterates over a dictionary containing multiple levels
    of nested dictionaries. 
    
    Returns a list in which the contained values are the keys and values in the dictionaries, in the 
    encountered order. This order is the same as what is found left to right if the dictionary was printed. 
    """
    out_list = []
    for key, value in dictionary.items():
        out_list.append(key)
        if isinstance(value, dict):
            for item in recursive_dict_iteration(value):
                out_list.append(item)
    return out_list
    
def output_hierarchical_gene_count(hierarchy_dict, K_term_dict, K_term_present, outfile):
    """
    Converts the hierarchical dictionary of categories into an ordered list. Per category,
    finds the associated genes from the K_term_dict, and finds which ones of these are present 
    using the K_term_present list. Outputs the number of found genes per category to the outfile.
    """
    out = open(outfile, "w")
    #Iterate over the hierarchical dictionary to find the contained items
    #in the appropriate hierarchical order
    hierarchical_tup_list = recursive_dict_iteration(hierarchy_dict)
    #Iterate over this ordered list of categories, checking the associated genes,
    #and whether these genes are present in the provided EggNOG output. Then write 
    #this per category
    for item in hierarchical_tup_list:
        leading_spaces = item[0]
        KEGG_cat = item[1]
        leading_spaces = leading_spaces - 12
        associated_kterms = K_term_dict[KEGG_cat]
        associated_kterms_present = []
        for kterm in associated_kterms:
            if kterm in K_term_present:
                associated_kterms_present.append(kterm)
        out.write("{}{}\t{}/{}\n".format((leading_spaces*"-"), item[1], len(associated_kterms_present), len(associated_kterms)))

####################################################################
#MAIN
####################################################################    
if __name__ == "__main__":
    #Acquire the input EggNOG file, output_filename, and the input database
    input_eggnog = sys.argv[1]
    out_file = sys.argv[2]
    k_term_db = sys.argv[3]
    #Parse the kterms per category into a dictionary
    K_term_dict, hierarchy_list = parse_KEGG_kterm_db(k_term_db)
    #Reconstruct the KEGG hierarchy
    hierarchy_dict = reconstruct_KEGG_hierarchy(hierarchy_list)
    #Parse the k terms found by EggNOG
    K_term_present = eggnog_parser(input_eggnog)
    #Output k terms found per category in hierarchical fashion, in full format
    output_hierarchical_gene_count(hierarchy_dict,K_term_dict, K_term_present, out_file)
    
    