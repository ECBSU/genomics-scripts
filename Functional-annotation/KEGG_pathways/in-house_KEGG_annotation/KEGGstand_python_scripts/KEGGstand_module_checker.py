#!/usr/bin/env python3
"""
Parses the K number from an Eggnog .annotation file. Uses a database of KEGG modules to find 
all possible combinations of genes that complete the pathway. Using the parsed K number,
checks which combination is most complete, and gives the highest completion per KEGG module.

Usage: (python) Eggnog_KEGG_completion_checker.py input.emapper.annotations (optional path_to_KEGG_module.txt)

Output:
Writes a .tsv file following the name: *input.emapper.annotations*_KEGG_completion.tsv. 
This tsv file contains a column containing the module name, the highest completion found for the module (as a fraction 0-1), 
and the k numbers of the genes that comprise this highest completion.
"""
###################
#Import statements
###################  

import sys
import os

###################
#File handling functions
###################  
def gen_line_reader(file_path):
    for line in open(file_path, "r"):
        yield line
        
def KEGG_module_reader(KEGG_module_file_path):
    """
    Reads a database of KEGG module definitions and outputs a dictionary
    where the key is "Modulenumber Modulename" and the value is the definition.
    
    Since some modules have multiple definitions, the value is given as a list.
    """
    KEGG_dict = {}
    name = False
    for line in gen_line_reader(KEGG_module_file_path):
        if line.startswith("#"):
            continue
        if not line.strip():
            continue
        if line.startswith("M00"):
            name = line.strip()
            if name not in KEGG_dict:
                KEGG_dict[name] = []
        else:
            KEGG_dict[name].append(line.strip())
    return KEGG_dict

def output_tsv(outputname, output_dict, minimum_completion):
    """
    Outputs a tsv containing per column the module name, the highest completion 
    found for the module (as a fraction 0-1), and the k numbers of the genes that 
    comprise this highest completion.
    """
    out = open(outputname, "w")
    out.write("#Entry\tName\tHighest_completion\tMost_complete_pathway\tNon-essential_genes_found\n")
    for modulename in output_dict:
        if float(output_dict[modulename][0]) >= float(minimum_completion):
            out.write("{}\t{}\t{}\t".format(modulename.partition(" ")[0], modulename.partition(" ")[2],output_dict[modulename][0]))
            out.write(",".join(output_dict[modulename][1]))
            out.write("\t")
            if len(output_dict[modulename][2]) > 0:
                out.write(",".join(output_dict[modulename][2]))
            else:
                out.write("None")
            out.write("\n")
    out.close()
    
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
    
###################
#Functions for parsing KEGG definitions
###################  

#Wrapper for KEGG definition parsing
def retrieve_all_possible_pathways(KEGG_definition):
    """
    Finds the different possible combinations of genes that complete the module in the KEGG_definition
    """
    final_list = [[]]
    #The "for reaction" structure is required since sometimes KEGG gives 2 definitions for the same module. 
    #It seems these are merging pathways required to for example generate substrates. Therefore, they will be treated as prerequisites for completion.
    for reaction in KEGG_definition:
        #Pluses and dashes are functionally the same
        if " " in reaction:
            reaction = reaction.replace(" ", "+")
        #Minuses denote non-essential genes. They will be considered here for calculating the gene combinations
        #but their absence will not count as pathway incompletion later. 
        if "-" in reaction:
            reaction = reaction.replace("-", "+")
        #To deal with nested brackets, first find the possible combination WITHIN each bracket, so these
        #can then be substituted into the appropriate brackets when parsing the complete definitions.
        #These possible combinations are put into poss_dict
        poss_dict = {}
        if "(" in reaction:
            #If there is brackets in the definition, first find the possibilities contained therein.
            poss_dict = find_bracket_possibilities(reaction)
        #Parse the different pathways
        poss_list = parse_possibilities(reaction, poss_dict)
        temp_list = []
        for option in final_list:
            for poss in poss_list:
                new_option = option[:] + poss
                temp_list.append(new_option)
        final_list = temp_list[:]   
    return final_list

#Wrapper for parsing (nested) brackets
def find_bracket_possibilities(KEGG_definition):
    """
    Finds all (nested) brackets, and iterates through each level of brackets to find the possibilties.
    Returns a completed dictionary of possibilities for the lowest level of brackets, which can be used to 
    analyze the KEGG_definition
    """
    #Finds the brackets and their level of nestedness
    brackets_list = find_bracket_contents(KEGG_definition)
    #Find the highest level of nested brackets
    highest = 0
    for element in brackets_list:
        if element[0] > highest:
            highest = element[0]
    poss_dict = {}
    #Iterate over each bracket level
    #For each bracket, make a list of lists representing each possible option and store them in the poss_dict
    for level in reversed(range(highest+1)):
        #iterate over each bracket pair in the KEGG definition
        for bracket_pair in brackets_list:
            #starting at the highest level of brackets, start searching for options
            if bracket_pair[0] == level:
                poss_dict["({})".format(bracket_pair[1])] = parse_possibilities(bracket_pair[1],poss_dict)
    return poss_dict

#Main parsing function    
def parse_possibilities(string, poss_dict):
    """
    Returns a list of lists, representing the options possible in the string
    according to the KEGG definitions. 
    
    Utilizes the poss_dict to resolve brackets
    """
    #To prevent parsing mistakes, replace brackets with empty brackets.
    #The removed contents of the brackets are present in the bracket_list,
    #in the order that they occur in the string
    bracket_list = []
    if "(" in string:
        string, bracket_list  = empty_brackets(string)
    #If the string has commas in it, it is a bracketed substring. Comma's will take priority
    #over plusses and spaces, and will thus be considered first    
    if "," in string:
        option_list = []
        for alternative in string.split(","):
            if alternative == "":
                continue
            elif "+" in alternative:
                poss_list = parse_pluses(alternative, bracket_list, poss_dict)
                #Since considered parts are separated by commas, each combination in the poss_list here 
                #denotes a branching path ALTERNATIVE to the other parts of the string being parsed here.
                #Therefore, add each possibility as a separate one to the option_list
                for poss in poss_list:
                    option_list.append(poss)     
            elif alternative == "()":
                #Again, since this is in a comma string, each possibility is an entirely different option
                bracket_possibilities = parse_bracket_possibilities(bracket_list.pop(0), poss_dict)
                for poss in bracket_possibilities:
                    option_list.append(poss)     
            #If there is only spaces, simply add the genes to existing options
            else:
                option_list.append([alternative])
    #If there is no commas in the string, then it is only comprised of additions.
    #simply add each addition, and add the options contained in each bracket
    else:
        option_list = [[]]
        for addition in string.split("+"):
            if addition == "":
                continue
            elif addition == "()":
                poss_list = parse_bracket_possibilities(bracket_list.pop(0), poss_dict)
                temp_list = []
                for poss in poss_list:
                    for option in option_list:
                        new_option = option[:]
                        for item in poss:
                            new_option.append(item)   
                        temp_list.append(new_option)
                option_list = temp_list[:]        
            #If there is only spaces, simply add the genes to existing options
            else:
                if len(option_list) == 0:
                    option_list = []
                for option in option_list:
                    option.append(addition)
    return option_list

#Secondary functions    
def parse_pluses(plus_string, bracket_list, poss_dict):
    """
    Parses a plus-containing substring found in brackets, or between 
    comma's. 
    
    Returns a list of lists where each nested list represents a possible combination
    of genes following the KEGG definition
    """
    out_list = [[]]
    for part in plus_string.split("+"):
        #If a part of the pluses are brackets, remove them and add the possibilities contained
        if part == "()":
            bracket_possibilities = parse_bracket_possibilities(bracket_list.pop(0), poss_dict)
            temp_list = []
            for poss in bracket_possibilities:
                for option in out_list:
                    new_option = option[:]
                    for item in poss:
                        new_option.append(item)   
                    temp_list.append(new_option)
            out_list = temp_list[:]        
        #If there is no brackets, simply add each part of the chain of pluses as a an addition to the existing combinations
        else:
            temp_list = []
            for option in out_list:
                new_option = option[:]
                new_option.append(part)   
                temp_list.append(new_option)
            out_list = temp_list[:]   
    return out_list

def find_bracket_pairs(string):
    """
    Returns a dictionary for each bracket pair, where the first bracket
    is the dict key, and the second is the value
    """
    match_bracket_index_dict = {}
    stack_list = []
    for index, char in enumerate(string):
        if char == "(":
            stack_list.append(index)
        elif char == ")":
            match_bracket_index_dict[stack_list.pop()] = index
    return dict(sorted(match_bracket_index_dict.items()))
    #return match_bracket_index_dict

def find_bracket_contents(string):
    """
    returns a list of tuples for each bracket pair in the string.
    Each tuple contains (bracket nestedness level, bracket contents).
    """
    stack = []
    bracket_list = []
    out_dict = {}
    for index, char in enumerate(string):
        if char == '(':
            stack.append(index)
        elif char == ')' and stack:
            start = stack.pop()
            bracket_list.append((len(stack), string[start + 1: index]))
    return bracket_list

def parse_bracket_possibilities(bracket_string, poss_dict):
    """
    If brackets are found, tries to find them in the poss_dict to give 
    the options they represent. The options are not added to an option list
    since from the brackets itself its impossible to tell if it is an addition or
    a branch.
    """
    if bracket_string in poss_dict:
        return poss_dict[bracket_string]
    else:
        print("Error, could not find {} in poss_dict".format(bracket_string))

def empty_brackets(string):
    """
    Returns the provided string, but with the brackets replaced with empty brackets. In addition, 
    returns a list containing the removed bracketed contents, so these can be retrieved later from 
    the poss_dict
    
    Only considers the lowest level of bracket present, assuming nested brackets are already included in 
    its possibilities
    """
    out_list = []
    #Make dictionary of bracket pairs
    bracket_dict = find_bracket_pairs(string)
    #Only consider the lowest level of brackets, ignoring any nested ones
    lowest_list = []
    for key in bracket_dict:
        lower_exists = False
        for key2 in bracket_dict:
            if key > key2 and bracket_dict[key] < bracket_dict[key2]:
                lower_exists = True
        if not lower_exists:
            lowest_list.append(key)
    #Store the lowest brackets in the out_list, so they can be called later
    for key in lowest_list:
        out_list.append(string[key:bracket_dict[key]+1])
    #Replace the lowest level brackets with empty brackets
    for key in reversed(lowest_list):
        string = string.replace(string[key:bracket_dict[key]+1], "()")
    return string, out_list

def non_essential_finder(KEGG_definition):
    """
    In KEGG definitions, minuses denote non-essential genes. These will be parsed here, and returned as a list.
    """
    out_list = []
    if "-" in KEGG_definition:
        #Find all minuses in the string
        i = 0
        minus_indices = []
        while i < len(KEGG_definition):
            if KEGG_definition[i] == "-":
                minus_indices.append(i)
            i+= 1
        #Find the genes after the minus, counting multiple if they are bracketed. 
        for i in minus_indices:
            string = KEGG_definition[i+1:]
            #If the minus is followed by brackets, inside are all non-essential genes
            if string.startswith("("):
                string = string.partition(")")[0].strip("(")
                splt = string.split(",")
                for comma_split in splt:
                    for plus_split in comma_split.split("+"):
                        out_list.append(plus_split)
            #Several instances of double minuses exist, such as M00814, M00840 or M00890.
            #Their meaning is unclear, but they do not seem to denote non-essential genes
            elif string.startswith("-"):
                continue
            #If it is not followed by brackets or minus, simply take the gene following the minus
            else:
                string = string.partition("+")[0].partition(",")[0].partition(" ")[0].partition("(")[0].partition(")")[0].partition("-")[0]
                if string.strip():
                    out_list.append(string)
    return out_list
    
###################
#Analysis wrapper function 
###################   

def pathway_completion_checker(KEGG_dict, kterms_list):
    """
    Uses several functions to find all potential combinations of genes that complete a KEGG module.
    Then find which of those combinations is most complete for the given genes in kterms_list.
    
    Returns a dictionary with the KEGG module name as key, and a tuple as value. The tuple
    contains the highest completion value found (as a 0-1 float),the combination of genes
    amounting to that completion as a list, and a list genes missing from the pathway that were still
    included during completion calculation because they were non-essential: (0.5, [gene1,gene2,gene3], [gene2])
    """
    out_dict = {}
    #Make a dictionary containing lists of all non-essential genes per module
    non_essential_dict = {}
    for module_name in KEGG_dict:
        for reaction in KEGG_dict[module_name]:
            #Since there can be multiple reactions, check if the dictionary entry already exists, if so append to it
            if module_name in non_essential_dict:
                for i in non_essential_finder(reaction):
                    non_essential_dict[module_name].append(i)
            else:
                non_essential_dict[module_name] = non_essential_finder(reaction)
    #Make a dictionary containing all the possible combinations of genes to complete the pathway
    for module_name in KEGG_dict:
        possible_combinations = retrieve_all_possible_pathways(KEGG_dict[module_name])
        if possible_combinations != [[]]:
            KEGG_dict[module_name] = possible_combinations
    #Iterate over the dictionary to check every pathway
    for module_name in KEGG_dict:
        #Iterate over each possible combinations, to check which is most complete
        highest_completion = -1
        best_pathway = []
        non_essential_found = []
        for combination in KEGG_dict[module_name]:
            gene_is_present = []
            #Check for each gene in the combination if it non-essential and if it is present in the tested organism
            for gene in combination:
                #If the gene is non-essential it will not be counted towards the completion of the combination
                #While not factored into the calculation, any non-essential gene associated with the module will be mentioned in the output if it was found to be present in the organism
                if gene in non_essential_dict[module_name]:
                    if gene in kterms_list and gene not in non_essential_found:
                        non_essential_found.append(gene)
                else: 
                    if gene in kterms_list:
                        gene_is_present.append(gene)
            #To accurately calculate the completion of essential genes, non-essential genes are removed from the combination
            for gene in non_essential_dict[module_name]:
                if gene in combination:
                    combination.remove(gene)
            number_of_genes = len(combination)
            completion = len(gene_is_present)/number_of_genes
            if completion > highest_completion:
                highest_completion = completion
                best_pathway = gene_is_present[:]
                out_dict[module_name] = (highest_completion, best_pathway, non_essential_found)
    return out_dict
    
####################################################################
#MAIN
####################################################################    
if __name__ == "__main__":
    #Obtain inputs
    eggnog_input = sys.argv[1]
    if len(sys.argv) > 2:
        KEGG_db = sys.argv[2]
    else:
        KEGG_db = "/bucket/HusnikU/Unit_Members/Arno/Scripts/KEGG_modules.txt"
    #Read in the KEGG module DB
    KEGG_dict = KEGG_module_reader(KEGG_db)
    #Parse out the K terms from eggnog output_name
    kterms_list = eggnog_parser(eggnog_input)
    #Check per pathway how complete it is based on the eggnog K terms
    completion_dict = pathway_completion_checker(KEGG_dict, kterms_list)
    #Output the found completion
    output_tsv(eggnog_input + "_KEGG_completion.tsv",completion_dict, 0)
    output_tsv(eggnog_input + "_KEGG_complete_modules.tsv",completion_dict, 1)
