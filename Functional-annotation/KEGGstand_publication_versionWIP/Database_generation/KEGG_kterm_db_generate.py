#!/usr/bin/env python3
"""
Will query the KEGG API server to obtain all the KEGG k terms, and write their relevant
information into a text file. 

Usage:
(python) KEGG_kterm_db_generate.py /path/to/db
"""
###################
#Import statements
###################  

from Bio.KEGG import REST
import sys
import os
import time

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

def find_processed_entries(file):
    """
    Reads provided file and creates a list of all KEGG entries
    within. 
    """
    entry_list = []
    for line in gen_line_reader(file):
        if line.startswith("ENTRY"):
            entry_list.append(line.split()[1])
    return entry_list
    
def k_term_list():
    """
    Connects to KEGG to request a list of all KEGG k terms.
    Parses out only the k term (removing the name) and stores
    it into a list.
    """
    k_list = []
    temp_list = REST.kegg_list("ko", org=None)
    for item in temp_list:
        k_term = item.split()[0]
        k_list.append(k_term)
    return k_list
    
####################################################################
#MAIN
####################################################################    
if __name__ == "__main__":
    #Aquire the desired database filename
    out_file = sys.argv[1]
    #Acquire the list of k terms from KEGG
    k_terms = k_term_list()
    #Check which k terms are already present in the database (if it exists).
    #Remove any terms already in the database from the list, so they wont be downloaded again
    if os.path.isfile(out_file):
        for term in find_processed_entries(out_file):
            if term in k_terms:
                k_terms.remove(term)
            else:
                print("error, found processed k term that is not in list: {}".format(term))
    
    #For any undownloaded terms, download them and append them to the database
    for term in k_terms:
        written = False
        while not written:
            #Here "try" is used because KEGG will sometimes refuse connection, causing an error.
            #If this happens, the script will pause for 2 minutes and try again
            tries = 0
            try:
                #obtain the KEGG entry and only write the part up to and including the BRITE hierarchy
                entry = REST.kegg_get(term, option=None)
                read = False
                with open(out_file, "a") as out:
                    for i in entry:
                        # Many KEGG entries contain extensive external links to other resources, sequences, and literature.
                        # These are not relevant here and skipped
                        if i.startswith("ENTRY"):
                            read = True
                        if i.startswith(("DBLINKS", "GENES", "REFERENCE")):
                            read = False
                        if read:
                            out.write(i)
                written = True
                tries = 0  # Reset tries after successful execution
            except:
                tries += 1
                print("No connection, waiting 2 minutes")
                time.sleep(120)
            if tries > 4:
                print("Connection failed 5 consecutive times. Script will exit")
                sys.exit()
        out.close()

