#!/usr/bin/env python3
"""
Processes and merges all the KEGGstand output in the provided directory, and visualizes them as a heatmap

Usage: (python) KEGGstand_graph_maker.py --in_dir /directory/of/KEGGstand/output --output plot_name.extension [arguments]

Output:
Image file of the specified heatmap
"""

###################
#Import statements
################### 
import sys
import os

###################
# Functions
###################

def gen_line_reader(file_path):
    """
    Generator function that allows reading a text file line
    by line, without reading the full file into memory
    """
    for line in open(file_path, "r"):
        yield line


def completion_tsv_reader(filepath):
    """
    Parser function for reading the completion values in the KEGGstand pipeline output.
    Returns a dictionary where the module is the key, and its completion fraction is the value
    """
    out_dict = {}
    for line in gen_line_reader(filepath):
        if line.startswith("#"):
            continue
        # Combine the KEGG Module code and its name into a single variable
        KEGG_name = "{} {}".format(line.split("\t")[0], line.split("\t")[1])
        completion = float(line.split("\t")[2])
        out_dict[KEGG_name] = completion
    return out_dict

def remove_modules_below_completion(module_dict, minimum_completion, method = "avg"):
    """
    Remove any modules from the input dictionary that are below the minimum completion.
    How this is calculated depends on the specified method either considering the average (avg)
    minimum (min), or maximum (max) across samples. 
    """
    out_dict = {}
    out_modules = []
    #Check which modules need to be entered into the output dictionary
    #Use a dummy loop to enter and iterate over the "modules" dictionary first
    for dummy in module_dict:
        for module in module_dict[dummy]:
            temp_list = []
            #Iterate over each organism/fasta, adding their completion for the module in question to the temp_list
            for org in module_dict:
                #Add the organism to the out_dict
                if org not in out_dict:
                    out_dict[org] = {}
                temp_list.append(module_dict[org][module])
            #Check if the average completion (if avg) or any completion (if not avg) is higher than the minimum
            #completion
            if method == "avg":
                if sum(temp_list)/len(temp_list) >= minimum_completion:
                    out_modules.append(module)
            elif method == "min":
                if min(temp_list) >= minimum_completion:
                    out_modules.append(module)
            elif method == "max":
                if max(temp_list) >= minimum_completion:
                    out_modules.append(module)
        break
    #Make the output dictionary
    for org in out_dict:
        for module in out_modules:
            out_dict[org][module] = module_dict[org][module]
    return out_dict

def retain_only_specified_modules(module_dict, search_string_list, filter_string_list):
    """
    Outputs a new dictionary removing or retaining modules based on the strings found in the specified
    lists. Search_string_list specifies strings of modules to be retained, while filter_string_list specifies
    strings that should be removed from the module_dict
    """
    temp_dict = {}
    #Remove any modules matching string in the filter_string_list
    if len(filter_string_list) > 0:
        for org in module_dict:
            for module in module_dict[org]:
                write = True
                for string in filter_string_list:
                    if string in module:
                        write = False
                if write:    
                    if org not in temp_dict:
                        temp_dict[org] = {}
                    temp_dict[org][module] = module_dict[org][module]
    #Only use the modules with a string matching the search_string_list
    if len(search_string_list) > 0:
        out_dict = {}
        for org in temp_dict:
            for module in temp_dict[org]:
                for string in search_string_list:
                    if string in module:
                        out_dict[org][module] = module_dict[org][module]
    else:
        out_dict = temp_dict
    return out_dict
    
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
        if line.startswith("Module: "):
            name = line.strip("Module: ").strip()
            if name not in KEGG_dict:
                KEGG_dict[name] = []
        if line.startswith("Definition: "):
            KEGG_dict[name].append(line.strip("Definition: ").strip())
    return KEGG_dict    
    
def category_collapser(module_dict, module_db_path, level, method, show_module_count):
    """
    'Collapses' the modules into their overarching categories. Requires the module database. Level 
    specifies the level of category, 1 being the broadest, and 2 being the more specific categories. 
    Method species the way to consolidate the modules into a category. The value can be the average (avg), 
    the lowest value (min), or the highest value (max). Show module count specifies that the new categories should
    list the number of modules they represent. 
    """
    out_dict = {}
    #Iterate over the module database, finding which modules belong to which categories.
    #Make a list for each category, listing the completions of included modules
    for line in gen_line_reader(module_db_path):
        if line.startswith("Module: "):
            name = line.strip("Module:").strip()
        if line.startswith("Class: "):
            clss = line.partition(";")[2].split(";")
            clss = clss[level-1].strip()
            for org in module_dict:
                if name in module_dict[org]:
                    if org not in out_dict:
                        out_dict[org] = {}
                    if clss not in out_dict[org]:
                        out_dict[org][clss] = []
                    out_dict[org][clss].append(module_dict[org][name])
    #Rename the categories, adding the number of modules that they represent
    if show_module_count:
        temp_dict = {}
        for org in out_dict:
            #Add new categories
            temp_dict[org] = {}
            for clss in out_dict[org]:
                    value  = out_dict[org][clss]
                    new_clss = "{} ({} modules)".format(clss, len(value))
                    temp_dict[org][new_clss] = value
        out_dict = temp_dict             
    #Consolidate the lists into single values depending on the specified method    
    for org in out_dict:
        for clss in out_dict[org]:
            if method == "avg":
                out_dict[org][clss] = sum(out_dict[org][clss]) / len(out_dict[org][clss])
            if method == "min":
                out_dict[org][clss] = min(out_dict[org][clss])
            if method == "max":
                out_dict[org][clss] = max(out_dict[org][clss])
    return out_dict
                
def heatmap(data_dict, outname, height, width, color, outformat, cluster_samples):
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    import scipy
    
    df = pd.DataFrame(data_dict)
    ax = sns.clustermap(df, metric="euclidean", method="single", cmap= color, figsize=(float(width), float(height)), linewidths=0.5, linecolor='white', col_cluster = cluster_samples)
    ax.fig.subplots_adjust(left=0.1)
    ax.fig.subplots_adjust(bottom=0.25)
    ax.ax_cbar.set_position((0.02, 0.6, 0.03, 0.2))
    plt.savefig(outname, format = outformat)     
            

####################################################################
#MAIN
####################################################################    
if __name__ == "__main__":
    #Acquire and handle input
    #########################
    if "--in_dir" in sys.argv:
        indir = sys.argv[sys.argv.index("--in_dir") + 1]
    elif "--in_files" in sys.argv:
        files = sys.argv[sys.argv.index("--in_files") + 1].split(",")
    else:
        print("No input specified, specify input directory with '--in_dir'")
        sys.exit()
    if "--output" in sys.argv:
        graphname = sys.argv[sys.argv.index("--output") + 1]
        if "." in graphname:
            outformat = graphname.split(".")[1].strip()
            if outformat != "png" and outformat != "pdf" and outformat != "ps" and outformat != "eps" and outformat != "svg":
                print("{} is an unrecognized output format, will default to .png format".format(outformat))
        else:
            print("No extension in the outputname, will default to .png format")
    else:
        graphname = "KEGGstand_heatmap.png"
        outformat = "png"
    #Initiate and check optional variables
    minimum_completion = 0
    filter_method = "min"
    collapse_method = "avg"
    height = 10
    width = 10
    color = "coolwarm"
    search_string_list = []
    filter_string_list = []
    cat_search_string_list = []
    cat_filter_string_list = []
    collapse_level = 0
    db_path = ""
    show_module_count = False
    if "--completion" in sys.argv:
        minimum_completion = float(sys.argv[sys.argv.index("--completion") + 1])
    if "--filter_method" in sys.argv:
        filter_method = sys.argv[sys.argv.index("--filter_method") + 1]
        if filter_method != "avg" and filter_method != "min" and filter_method != "max":
            print('Error, specified method needs to be "min","max", or "avg". {} was not recognized.'.format(filter_method))
            exit()
    if "--collapse_method" in sys.argv:
        collapse_method = sys.argv[sys.argv.index("--collapse_method") + 1]
        if collapse_method != "avg" and collapse_method != "min" and collapse_method != "max":
            print('Error, specified method needs to be "min","max", or "avg". {} was not recognized.'.format(collapse_method))
            exit()
    if "--dimensions" in sys.argv:
        height,width = sys.argv[sys.argv.index("--dimensions") + 1].split(",")
    if "--module_search" in sys.argv:
        if "," in sys.argv[sys.argv.index("--module_search") + 1]:
            search_string_list = sys.argv[sys.argv.index("--module_search") + 1].split(",")
        else:
            search_string_list = [sys.argv[sys.argv.index("--module_search") + 1]]
    if "--module_filter" in sys.argv:
        if "," in sys.argv[sys.argv.index("--module_filter") + 1]:
            filter_string_list = sys.argv[sys.argv.index("--module_filter") + 1].split(",")
        else:
            filter_string_list = [sys.argv[sys.argv.index("--module_filter") + 1]]
    if "--category_search" in sys.argv:
        if "," in sys.argv[sys.argv.index("--category_search") + 1]:
            cat_search_string_list = sys.argv[sys.argv.index("--category_search") + 1].split(",")
        else:
            cat_search_string_list = [sys.argv[sys.argv.index("--category_search") + 1]]
    if "--category_filter" in sys.argv:
        if "," in sys.argv[sys.argv.index("--category_filter") + 1]:
            cat_filter_string_list = sys.argv[sys.argv.index("--category_filter") + 1].split(",")
        else:
            cat_filter_string_list = [sys.argv[sys.argv.index("--category_filter") + 1]]
    if "--color" in sys.argv:
        color = sys.argv[sys.argv.index("--color") + 1]         
    if "--db" in sys.argv:
        db_path = sys.argv[sys.argv.index("--db") + 1]   
    if "--collapse" in sys.argv:
        collapse_level = int(sys.argv[sys.argv.index("--collapse") + 1])
    if collapse_level > 2:
        print("Level of category collapse has to be lower than 3. Defaulting to 2")
        collapse_level = 2
    elif collapse_level > 0 and db_path == "":
        print("Error, collapse was specified, but no module database path was given. The module database is required for collapsing the categories")
        exit()
    if "--show_module_count" in sys.argv:
        show_module_count = True
        
    #######################################
    module_dict = {}
    #Iterate over the files in the input directory or read the specified files         
    if "--in_dir" in sys.argv:
        cluster_samples = True
        for root, dirs, files in os.walk(indir): 
            for file in files: 
                #Obtain the completion value per module per fasta
                if file.endswith(".emapper.annotations_KEGG_completion.tsv"):
                    org_name = file.partition(".emapper")[0]
                    modules = completion_tsv_reader(os.path.join(root, file))
                    module_dict[org_name] = modules
    else:
        cluster_samples = False
        for file in files: 
            #Obtain the completion value per module per fasta
            if file.endswith(".emapper.annotations_KEGG_completion.tsv"):
                org_name = file.partition(".emapper")[0].rpartition("/")[2]
                modules = completion_tsv_reader(file)
                module_dict[org_name] = modules


    #Remove any modules that dont meet the completion requirements
    if minimum_completion > 0 or filter_method != "min":
        module_dict = remove_modules_below_completion(module_dict, minimum_completion, filter_method) 
    #Remove modules based on string searches
    if search_string_list != [] or filter_string_list != []:
        module_dict = retain_only_specified_modules(module_dict, search_string_list, filter_string_list)
    #Collapse modules into categories
    if collapse_level != 0:
        module_dict = category_collapser(module_dict, db_path, collapse_level, collapse_method, show_module_count)
        #Remove categories based on string searches
        if cat_filter_string_list != [] or cat_search_string_list != []:
            module_dict = retain_only_specified_modules(module_dict, cat_search_string_list, cat_filter_string_list)
    
    heatmap(module_dict, graphname, height, width, color, outformat, cluster_samples)
