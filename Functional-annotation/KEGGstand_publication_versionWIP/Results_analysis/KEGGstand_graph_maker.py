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
import argparse

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
        if len(temp_dict) == 0:
            temp_dict = module_dict
        out_dict = {}
        for org in temp_dict:
            out_dict[org] = {}
            for string in search_string_list:
                for module in temp_dict[org]:
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
                
def heatmap(data_dict, outname, height, width, color, outformat, cluster_samples, cluster_modules):
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    import scipy
    
    df = pd.DataFrame(data_dict)
    ax = sns.clustermap(df, metric="euclidean", method="single", cmap= color, figsize=(float(width), float(height)), linewidths=0.5, linecolor='white', col_cluster = cluster_samples, row_cluster = cluster_modules)
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
    parser = argparse.ArgumentParser(description = 'Processes and merges all the KEGGstand output in the provided directory, and visualizes them as a heatmap')
    in_group = parser.add_mutually_exclusive_group(required=True)
    in_group.add_argument("-in", "--in_dir", help="Path to the directory containing KEGGstand output.", dest='indir')
    in_group.add_argument("--in_files", help="Comma-delimited list of KEGGstand output files. If you provide this list INSTEAD of --input_dir, the clustering of samples will be disabled, and instead the order of samples will reflect the order of the provided list.", dest='files')
    parser.add_argument("-out", "--output",help="output file", default="KEGGstand_heatmap.png")
    #Optional variables
    optional_grp = parser.add_argument_group("Optional parameters")
    optional_grp.add_argument("--completion",default=0,type=float, help="Desired completion values, used to avoid outputting uninformative modules. By default will only consider modules where the lowest completion found for any sample is above or equal to the specified completion.")
    optional_grp.add_argument("--filter_method", choices={"avg","max","min"}, default="min", help='Changes the calculation of the desired completion values. Options are "min", "max", and "avg". "min" is the default behaviour, requiring the lowest completion to be above or equal to the desired completion. "max" specifies the inverse, with the highest value needing to be equal or higher. "avg" specifies that the average module completion needs to be higher or equal to the specified completion.')
    optional_grp.add_argument("--collapse_method", choices={"avg","max","min"}, default="avg", help='Specifies the value to be given for the collapsed categories. Default is the average completion of the contained modules (avg). Alternatives are "min" (lowest value) and "max" (highest value)')
    optional_grp.add_argument('--height', default=10, help="Desired height of the created plot in inches")
    optional_grp.add_argument('--width', default=10, help="Desired width of the created plot in inches")
    optional_grp.add_argument("--module_search", help="Specifies a string that needs to be PRESENT in the module name. Modules WITHOUT this string are not considered. Can give multiple strings by delimiting with a comma. If spaces are included the string should be within quotation marks")
    optional_grp.add_argument("--module_filter", help="Specifies a string that needs to be ABSENT in the module name. Modules WITH this string are not considered. Can give multiple strings by delimiting with a comma. If spaces are included the string should be within quotation marks")
    optional_grp.add_argument("--category_search", help="Specifies that the modules are instead to be listed by their overarching categories (eg rather than showing the Citrate cycle and Pentose phosphate pathway completions as separate, their results are combined into a Central carbohydrate metabolism category). The number represents the level of detail for the categories. 1 is the broadest categories (eg energy metabolism), 2 is slightly more specific categories (eg methane metabolism). If this is specified, --db needs to be specified as well.")
    optional_grp.add_argument("--category_filter", help="Functions like the --module_filter argument, but instead of disregarding modules with the specified string, will disregard the categories with the specified string.")
    optional_grp.add_argument("--color", default="coolwarm", help="Desired colorscheme of the heatmap. Default is 'coolwarm'.")
    optional_grp.add_argument("--db", default="", help="Path to the module database used to generate the results. Used to find the categories to which the modules belong.")
    optional_grp.add_argument("--collapse", default=0, type=int, choices={0,1,2}, help="Specifies that the modules are instead to be listed by their overarching categories (eg rather than showing the Citrate cycle and Pentose phosphate pathway completions as separate, their results are combined into a Central carbohydrate metabolism category). The number represents the level of detail for the categories. 1 is the broadest categories (eg energy metabolism), 2 is slightly more specific categories (eg methane metabolism). If this is specified, --db needs to be specified as well.")
    optional_grp.add_argument("--show_module_count", action="store_true", help="Specifies that the number of contained modules are to be listed in the name of the categories. Only does anything if --collapse is given.")
    optional_grp.add_argument("--no_module_cluster", action="store_false", help="By providing this argument, clustering of the modules/categories will be disabled. By default, modules will then be given in order of their entry code (M00001, M00002 etc.). Alternatively, you can specify your own order by providing a list using the --module_search (or --category_search if --collapse is specified) argument. In that case, the module order will reflect the order of strings given. Note that if multiple modules match the given string(s), they will be ordered according to the entry code.")
    args = parser.parse_args()

    if "." in args.output:
            outformat = args.output.split(".")[1].strip()
            if outformat not in {"png", "pdf", "ps", "eps", "svg"}:
                print("{} is an unrecognized output format, will default to .png format".format(outformat))
    else:
        print("No extension in the outputname, will default to .png format")

    if args.collapse > 0 and args.db_path is None:
        parser.error("Error: collapse was specified, but no module database path was given. The module database is required for collapsing the categories.")
    
    search_string_list = [s.strip() for s in args.module_search.split(",")] if args.module_search else []
    filter_string_list = [s.strip() for s in args.module_filter.split(",")] if args.module_filter else []
    cat_search_string_list = [s.strip() for s in args.category_search.split(",")] if args.category_search else []
    cat_filter_string_list = [s.strip() for s in args.category_filter.split(",")] if args.category_filter else []
    #######################################
    module_dict = {}
    #Iterate over the files in the input directory or read the specified files
    if args.indir:
        cluster_samples = True
        for root, dirs, files in os.walk(args.indir): 
            for file in files: 
                #Obtain the completion value per module per fasta
                if file.endswith(".emapper.annotations_KEGG_completion.tsv"):
                    org_name = file.partition(".emapper")[0]
                    modules = completion_tsv_reader(os.path.join(root, file))
                    module_dict[org_name] = modules
    else:
        cluster_samples = False
        for file in args.files: 
            #Obtain the completion value per module per fasta
            if file.endswith(".emapper.annotations_KEGG_completion.tsv"):
                org_name = file.partition(".emapper")[0].rpartition("/")[2]
                modules = completion_tsv_reader(file)
                module_dict[org_name] = modules

    #Remove any modules that dont meet the completion requirements
    if args.completion > 0 or args.filter_method != "min":
        module_dict = remove_modules_below_completion(module_dict, args.completion, args.filter_method) 
    #Remove modules based on string searches
    if search_string_list != [] or filter_string_list != []:
        module_dict = retain_only_specified_modules(module_dict, search_string_list, filter_string_list)
    #Collapse modules into categories
    if args.collapse != 0:
        module_dict = category_collapser(module_dict, args.db, args.collapse, args.collapse_method, args.show_module_count)
        #Remove categories based on string searches
        if cat_filter_string_list != [] or cat_search_string_list != []:
            module_dict = retain_only_specified_modules(module_dict, cat_search_string_list, cat_filter_string_list)
    
    heatmap(module_dict, args.output, args.height, args.width, args.color, outformat, cluster_samples, args.no_module_cluster)
