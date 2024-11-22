#!/usr/bin/env python3
"""
Will query the KEGG API server to obtain all the KEGG module defintions, and writes them
to a text file

Usage:
(python) KEGG_module_db_generate.py /path/to/db
"""
##############
#Import statements
##############

from Bio.KEGG import REST
import sys
import os
import time

##############
#Functions
##############

def gen_line_reader(file_path):
    """
    Generator function that allows reading a text file line 
    by line, without reading the full file into memory
    """
    for line in open(file_path, "r"):
        yield line

def find_processed_entries(file):
    """
    Reads provided file and creates a list of all KEGG entries
    within. 
    """
    entry_list = []
    for line in gen_line_reader(file):
        if line.startswith("Module:"):
            entry_list.append(line.split()[1])
    return entry_list
    
def list_modules():
    """
    Connects to KEGG to request a list of all KEGG k terms.
    Parses out only the k term (removing the name) and stores
    it into a list.
    """
    module_list = []
    temp_list = REST.kegg_list("module", org=None)
    for item in temp_list:
        m_term = item.split()[0]
        module_list.append(m_term)
    return module_list
    
####################################################################
#MAIN
####################################################################    
if __name__ == "__main__":
    #Aquire the desired database filename
    out_file = sys.argv[1]
    #Acquire the list of k terms from KEGG
    module_list = list_modules()
    #Check which k terms are already present in the database (if it exists).
    #Remove any terms already in the database from the list, so they wont be downloaded again
    if os.path.isfile(out_file):
        for mod in find_processed_entries(out_file):
            if mod in module_list:
                module_list.remove(mod)
            else:
                print("error, found processed module that is not in list: {}".format(term))
                
    #For any undownloaded terms, download them and append them to the database
    for mod in module_list:
        print(mod)
        written = False
        while not written:
            #Here "try" is used because KEGG will sometimes refuse connection, causing an error.
            #If this happens, the script will pause for 2 minutes and try again
            tries = 0
            try:
                #obtain the KEGG entry and only write the part up to and including the BRITE hierarchy
                entry = REST.kegg_get(mod, option=None)
                defin = False
                clss = False
                with open(out_file, "a") as out:
                    for i in entry:
                        #Extract only the name and definition from the module entry
                        if i.startswith("NAME"):
                            out.write("Module: {} {}\n".format(mod, i.partition(" ")[2].strip()))
                        if i.startswith("DEFINITION"):
                            defin = True
                        if i.startswith("CLASS"):
                            clss = True
                        if i.startswith(("PATHWAY", "REACTION", "COMPOUND", "ORTHOLOGY", "REFERENCE")):
                            defin = False
                            clss = False
                        if defin:
                            out.write("Definition: {}\n".format(i.strip("DEFINITION").strip()))
                        if clss:
                            out.write("Class: {}\n".format(i.strip("CLASS").strip()))
                    written = True
                    tries = 0
            except:
                tries += 1
                print("No connection, waiting 2 minutes")
                time.sleep(120)
            if tries > 4:
                print("Connection failed 5 consecutive times. Script will exit")
                sys.exit()
        out.close()
