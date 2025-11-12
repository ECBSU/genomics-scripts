#!/usr/bin/env python3
"""
Parses the K number from an Eggnog .annotation file. Uses a database of KEGG k term entries to both
find to which (sub)category each k term belongs, and to reconstruct the KEGG hierarchy of nested
categories. Then per (sub)category checks which k terms are present based on the EggNOG output.

Usage: (python) KEGGstand_BRITE_checker.py input.emapper.annotations output_file KEGG_k_term_database

Output:
Writes a tab-delimited file that lists each KEGG category along with its subcategories. After
a tab, the number of genes found for the (sub)category is given, along with the total number of
genes known for this category.
"""
from typing import List, Tuple, Set, Dict
import re
import json
from pathlib import Path
import argparse
import networkx as nx
import pandas as pd


def gen_line_reader(file_path):
    """
    Generator function that allows reading a text file line
    by line, without reading the full file into memory
    """
    for line in open(file_path, "r"):
        yield line


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



def read_kterm_json(json_path):
    with open(json_path, "r") as s:
        kterms = json.load(s)
    return kterms


def convert_dict_to_edgelist(treedict, terminal_node=None) -> List[Tuple[str, str]]:
    # Takes in the json dict of kterms and converts BRITE part {1:{2:3} to the graph edges [(1,2)(2,3)]
    edges = []
    def get_edges(treedict, parent="root"):
        name = next(iter(treedict.keys()))
        if parent is not None:
            edges.append((parent, name))

        this_child = treedict[name]

        if isinstance(this_child, str):
            edges.append((name, this_child))
        else:
            for item in this_child:
                if isinstance(item, list):
                    for el in item:
                        get_edges(el, parent=name)
                elif isinstance(item, dict):
                    get_edges(item, parent=name)
                elif isinstance(item, str):
                    get_edges(this_child, parent=name)
    get_edges(treedict)
    terminal_node_edges = (edges[-1][-1], terminal_node)
    edges.append(terminal_node_edges)
    return edges


def convert_edgelist_to_graph(kterm_edgelist, kterm_dict):
    # Main function that creates a graph from the list of adjacent nodes
    kterm_graph = nx.DiGraph()
    for i in range(0, len(kterm_edgelist)):
        kterm_id = kterm_dict[i]["ENTRY"]
        kterm_brite = kterm_dict[i]["BRITE"]
        orthology_name = "KEGG Orthology (KO) [BR:ko00001]"
        if isinstance(kterm_brite, list):
            for el in kterm_brite:
                if orthology_name in el.keys():
                    this_el_edgelist = convert_dict_to_edgelist(el, kterm_id)
                    kterm_graph.add_edges_from(this_el_edgelist)
        else:
            this_term_edgelist = convert_dict_to_edgelist(kterm_brite, kterm_id)
            kterm_graph.add_edges_from(this_term_edgelist)
    return kterm_graph


def get_terminal_nodes(g: nx.DiGraph) -> Set[str]:
    terminal_nodes = [n for n, d in g.degree if d == 1]
    if "root" in terminal_nodes:
        terminal_nodes.remove("root")
    terminal_nodes = set(terminal_nodes)
    return terminal_nodes


def calc_kterm_enrichment(kterm_graph: nx.DiGraph, present_kterms: Set[str]) -> Dict[str, Tuple[int, int, str]]:
    # For each node will find total number of kterms associated with it, and number of kterms from eggnog
    terminal_node_list = get_terminal_nodes(kterm_graph)
    enrichments = dict()
    for node in kterm_graph.nodes:
        des = list(nx.descendants(kterm_graph, node))
        this_node_total_kterms = [n for n in des if n in terminal_node_list]
        this_node_pres_kterms = set(this_node_total_kterms).intersection(present_kterms)
        this_node_pres_kterms_str = ",".join(this_node_pres_kterms)
        enrichments[node] = (len(this_node_pres_kterms), len(this_node_total_kterms), this_node_pres_kterms_str)
    return enrichments


def dedupe_hierarchy(lines: List[Tuple[str, str]]) -> List[Tuple[str, str]]:
    # keep only first instances of hierarchical levels and 0 kterms
    key_fun = lambda x: x[1]
    seen = set()
    out = []
    for line in lines:
        key = key_fun(line)
        if re.match(r"^K\d{5}", key):
            continue
        if key not in seen:
            seen.add(key)
            out.append(line)
    return out

def path_top_bottom(graph: nx.DiGraph, node: Tuple[str, str]):
    # Do the depth first search on the graph to find all the connected nodes
    # Should start from the root
    # The input is node tuple ("", node_id). The hierarchy is stored in the first value of tuple.
    # The output is a list of node tuples with hierarchy in the first value of tuple.
    all_paths = []
    g_path = [node]
    key_fun = lambda x: x[1]
    def find_sux(graph, node, g_path, depth=0):
        sux = list(graph.successors(key_fun(node)))
        if sux == []:
            all_paths.append(g_path)
        else:
            depth += 1
            for s in sux:
                new_s_name = (depth * "-", s)
                this_suc_g_path = g_path.copy()
                this_suc_g_path.append(new_s_name)
                find_sux(graph, new_s_name, this_suc_g_path, depth)

    find_sux(graph, node, g_path)

    joined_paths = []
    for path in all_paths:
        joined_paths.extend(path)
    return joined_paths


def drop_orthology_and_root_from_hier(hier_repr: List[Tuple[str, str]]):
    # Remove orthology and root to decrease number of dashes
    new_hier = []
    things_to_remove = {"KEGG Orthology (KO) [BR:ko00001]", "root"}
    for line in hier_repr:
        if line[1] in things_to_remove:
            continue
        else:
            new_line = (line[0][2:], line[1])
            new_hier.append(new_line)
    return new_hier


def create_hier_repr(kterm_graph) -> List[Tuple[str, str]]:
    # Starting from the root do depth first search
    # The output is a list of node tuples with hierarchy in the first value of tuple.
    #hier_repr = [("---", node_name)]
    full_hier = path_top_bottom(kterm_graph, ("","root"))
    dedup_hier = dedupe_hierarchy(full_hier)
    hier_repr = drop_orthology_and_root_from_hier(dedup_hier)
    return hier_repr


def add_enrichment_to_hier_list(hier_repr, enrichments) -> List[Tuple[str, str, int, int, str]]:
    #hier_enrich = ["---", node_name, 10, 20, present_kterms_str_list]
    hier_enrich = []
    for node in hier_repr:
        node_name = node[1]
        enrichment = enrichments[node_name]
        new_node = (node[0], node[1], *enrichment)
        hier_enrich.append(new_node)
    return hier_enrich


def create_kterms_graph(kterm_json) -> nx.DiGraph:
    # convert json file to the graph representation
    kterms_edges = []
    for term in kterm_json:
        kterm_edgelist = convert_dict_to_edgelist(term)
        kterms_edges.append(kterm_edgelist)
    kterms_graph = convert_edgelist_to_graph(kterms_edges, kterm_json)
    return kterms_graph

def create_hierarchical_representation(kterms_graph, eggnog_kterms) -> List[Tuple[str, str, int, int, str]]:
    # Create hierarchical representation based on the graph
    # Add information about number of kterms per hierarchical level (enrichment)
    eggnog_kterm_set = set(eggnog_kterms)
    enrichments = calc_kterm_enrichment(kterms_graph, eggnog_kterm_set)
    hier_repr = create_hier_repr(kterms_graph)
    hier_enrich = add_enrichment_to_hier_list(hier_repr, enrichments)
    return hier_enrich


def convert_to_df_and_filter_zero(hier_repr_list) -> pd.DataFrame:
    column_names = ["hierarchy_level", "kegg_orthology_name", "num_kterms_present", "num_total_kterms", "present_kterms"]
    res_df = pd.DataFrame.from_records(hier_repr_list, columns=column_names)
    zero_matches = res_df.index[res_df["num_kterms_present"] == 0]
    res_df.drop(labels=zero_matches, axis=0, inplace=True)
    return res_df


def parseargs() -> Tuple[Path, Path, Path]:
    parser = argparse.ArgumentParser(description="Estimate completion of the modules")
    parser.add_argument("-k", help="Path to the k terms json file", required=True, dest="kterm_json_path", type=Path)
    parser.add_argument("-e", help="Path to the eggnog file", required=True, dest="eggnogfile_path", type=Path)
    parser.add_argument("-o", help="Path to the output tsv table", required=True, dest="out_path", type=Path)
    args = parser.parse_args()
    return args.kterm_json_path, args.eggnogfile_path, args.out_path


def resolve_rel_path_list(path_list: List[Path]):
    res_paths = []
    for p in path_list:
        res_paths.append(p.resolve())
    return res_paths


def main():
    kterm_json_path, eggnogfile_path, out_path = parseargs()
    kterm_json_path, eggnogfile_path, out_path = resolve_rel_path_list([kterm_json_path, eggnogfile_path, out_path])

    print(f"Reading kterms from {kterm_json_path}")
    kterm_json = read_kterm_json(kterm_json_path)

    print(f"Parsing eggnog from {eggnogfile_path}")
    eggnog_kterms = eggnog_parser(eggnogfile_path)

    print("Creating summary")
    kterms_graph = create_kterms_graph(kterm_json)
    hier_repr = create_hierarchical_representation(kterms_graph, eggnog_kterms)
    #Output k terms found per category in hierarchical fashion, in full format

    print(f"Saving summary to {out_path}")
    res_df = convert_to_df_and_filter_zero(hier_repr)
    res_df.to_csv(out_path, sep="\t", index=False)


if __name__ == "__main__":
    main()
