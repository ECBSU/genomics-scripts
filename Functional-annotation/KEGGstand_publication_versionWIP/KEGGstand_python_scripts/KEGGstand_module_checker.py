from typing import List, Dict, Union, Tuple, TypedDict
import re
from pathlib import Path
from collections import defaultdict

import networkx as nx
import pandas as pd
import argparse


class Enrichment(TypedDict):
    completion: float 
    present_genes: List[str] 
    pathway: List[str]
    optional: List[str]


######################### Parsing the KEGG compressed graph representation
def tokenize(pathway_str):
    """Extract K-numbers and operators from pathway string."""
    pattern = r"K\d+|[(),+\-\s]"
    return [m.group() for m in re.finditer(pattern, pathway_str)]


def find_closing_paren(tokens, open_idx):
    """Find matching closing parenthesis for opening paren at open_idx."""
    depth = 1
    for i in range(open_idx + 1, len(tokens)):
        if tokens[i] == "(":
            depth += 1
        elif tokens[i] == ")":
            depth -= 1
            if depth == 0:
                return i
    return len(tokens) - 1


def parse_element(tokens, idx, end_idx):
    """
    Parse single element: K-number or parenthesized group.
    Returns: (edges, optional_nodes, entry_nodes, exit_nodes, next_index)
    """
    if tokens[idx] == "(":
        close_idx = find_closing_paren(tokens, idx)
        edges, optional, entry, exit = parse_alternatives(tokens, idx + 1, close_idx - 1)
        return edges, optional, entry, exit, close_idx + 1
    else:
        k_num = tokens[idx]
        return [], set(), [k_num], [k_num], idx + 1


def parse_operator_sequence(tokens, idx, end_idx):
    """
    Parse sequence connected by + (required) or - (optional) operators.
    Returns: (edges, optional_nodes, entry_nodes, exit_nodes, next_index)
    """
    all_edges = []
    all_optional = set()
    
    # Parse first element
    edges, optional, entry, exit, idx = parse_element(tokens, idx, end_idx)
    all_edges.extend(edges)
    all_optional.update(optional)
    
    sequence_entry = entry
    current_exits = exit
    skip_nodes = []  # Nodes that can skip over optional sections
    
    # Continue while we see + or - operators
    while idx <= end_idx and tokens[idx] in ("+", "-"):
        operator = tokens[idx]
        idx += 1
        
        # Parse next element
        edges, optional, next_entry, next_exit, idx = parse_element(tokens, idx, end_idx)
        all_edges.extend(edges)
        all_optional.update(optional)
        
        if operator == "-":
            # Optional section: mark nodes as optional
            all_optional.update(next_entry)
            all_optional.update(next_exit)
            
            # Track where we can skip from (start of optional chain)
            if not skip_nodes:
                skip_nodes = list(current_exits)
            
            # Connect through optional path
            for curr in current_exits:
                for nxt in next_entry:
                    all_edges.append((curr, nxt))
            
            current_exits = next_exit
            
        else:  # operator == "+"
            # Required section: connect normally
            for curr in current_exits:
                for nxt in next_entry:
                    all_edges.append((curr, nxt))
            
            # Close any skip path
            if skip_nodes:
                # Merge: skip_nodes can jump to next_exit, or go through optional path
                current_exits = list(set(skip_nodes + next_exit))
                skip_nodes = []
            else:
                current_exits = next_exit
    
    # If we ended with optional sections, enable skip path
    if skip_nodes:
        current_exits = list(set(skip_nodes + current_exits))
    
    return all_edges, all_optional, sequence_entry, current_exits, idx


def parse_space_sequence(tokens, idx, end_idx):
    """
    Parse space-separated sequence (sequential steps).
    Returns: (edges, optional_nodes, entry_nodes, exit_nodes, next_index)
    """
    all_edges = []
    all_optional = set()
    
    # Parse first operator sequence
    edges, optional, entry, exit, idx = parse_operator_sequence(tokens, idx, end_idx)
    all_edges.extend(edges)
    all_optional.update(optional)
    
    sequence_entry = entry
    current_exits = exit
    
    # Continue while we see spaces (indicating sequential steps)
    while idx <= end_idx and tokens[idx] == " ":
        idx += 1
        if idx > end_idx or tokens[idx] == ",":
            break
        
        # Parse next operator sequence
        edges, optional, next_entry, next_exit, idx = parse_operator_sequence(tokens, idx, end_idx)
        all_edges.extend(edges)
        all_optional.update(optional)
        
        # Connect previous exits to next entries
        for prev in current_exits:
            for nxt in next_entry:
                all_edges.append((prev, nxt))
        
        current_exits = next_exit
    
    return all_edges, all_optional, sequence_entry, current_exits, idx


def parse_alternatives(tokens, start_idx, end_idx):
    """
    Parse comma-separated alternatives (parallel paths).
    Returns: (edges, optional_nodes, entry_nodes, exit_nodes)
    """
    all_edges = []
    all_optional = set()
    all_entries = []
    all_exits = []
    
    idx = start_idx
    while idx <= end_idx:
        if tokens[idx] == ",":
            idx += 1
            continue
        
        # Parse one alternative (space-separated sequence)
        edges, optional, entry, exit, idx = parse_space_sequence(tokens, idx, end_idx)
        
        all_edges.extend(edges)
        all_optional.update(optional)
        all_entries.extend(entry)
        all_exits.extend(exit)
    
    return all_edges, all_optional, all_entries, all_exits


def parse_pathway(pathway_str):
    """Parse KEGG pathway string into edges and optional nodes."""
    tokens = tokenize(pathway_str)
    
    all_edges = []
    all_optional = set()
    current_exits = ["BEGIN"]
    
    idx = 0
    while idx < len(tokens):
        if tokens[idx] == " ":
            idx += 1
            continue
        
        # Parse one top-level step
        edges, optional, entry, exit, idx = parse_operator_sequence(tokens, idx, len(tokens) - 1)
        
        # Connect previous step to current step
        for prev in current_exits:
            for curr in entry:
                all_edges.append((prev, curr))
        
        all_edges.extend(edges)
        all_optional.update(optional)
        current_exits = exit
    
    # Connect to END
    for node in current_exits:
        all_edges.append((node, "END"))
    
    return all_edges, all_optional


def sanitise_pathway_str(pathway_str_input):
    pathway_str_s = re.sub("--", "", pathway_str_input)
    pathway_str_s = re.sub(r"\s-K", " K", pathway_str_s)
    pathway_str_s = re.sub(r"\s+", " ", pathway_str_s)
    pathway_str_s = pathway_str_s.strip()
    return pathway_str_s


def create_pathway_graph(pathway_str_input):
    """Create NetworkX DiGraph from KEGG pathway string."""
    pathway_str = sanitise_pathway_str(pathway_str_input)
    edges, optional_nodes = parse_pathway(pathway_str)
    
    G = nx.DiGraph()
    G.add_edges_from(edges)
    
    optional_attr = dict()
    for node in G.nodes():
        optional_attr[node] = node in optional_nodes
    nx.set_node_attributes(G, optional_attr, "is_optional")
    return G


################################# Shortest path search
def handle_ambiguous_path(graph, target_nodes):
    # When there are nodes that do not belong to the unique shortest path
    # will find the the shortest path from beginning to end that includes
    # maximum nodes from the user

    sp_list = []    
    for target in target_nodes:
        segment_1 = nx.shortest_path(graph, "BEGIN", target)
        segment_2 = nx.shortest_path(graph, target, "END")
        this_node_sp = segment_1 + segment_2
        sp_list.append(this_node_sp)
    
    optional_dict = nx.get_node_attributes(graph, "is_optional")
    path_info = defaultdict(dict)
    path_scores = []
    for i, path in enumerate(sp_list):
        path_info[i]["length"] = -len(sp_list)
        path_info[i]["present"] = 0
        for target in target_nodes:
            if target in path:
                is_optional = optional_dict[target]
                if is_optional:
                    path_info[i]["present"] += 1
                else:
                    path_info[i]["present"] += 2
        path_info[i]["score"] = path_info[i]["length"] + path_info[i]["present"]
        path_scores.append(path_info[i]["score"])		   
    best_path_id = path_scores.index(max(path_scores))
    best_path = sp_list[best_path_id]
    return best_path
            	   
        

def find_shortest_path_through(graph: nx.DiGraph, target_nodes: List[str]) -> List[str]:
    # Build path segments
    full_path = []
    extended_target_nodes = target_nodes + ["END"]
    current = "BEGIN"
    
    is_path_broken = False
    for target in extended_target_nodes:
        try:
            segment = nx.shortest_path(graph, current, target)
        except nx.NetworkXNoPath:
            is_path_broken = True
            print("Problematic assignment ", target_nodes)
            break
        if full_path != []:
            full_path.extend(segment[1:])  # Avoid repeating current node
        else:
            full_path.extend(segment)
        current = target
    
    if is_path_broken:
        full_path = handle_ambiguous_path(graph, target_nodes)
    return full_path


def process_all_kegg_modules_to_pathways(kegg_dict: Dict[str, List[str]]) -> Dict[str, nx.DiGraph]:
    kegg_pathways = dict()
    for k_id, pathway_kegg in kegg_dict.items():
        pathway_str = pathway_kegg[0]
        if "M" in pathway_str:
            continue
        pathway_graph = create_pathway_graph(pathway_str)
        kegg_pathways[k_id] = pathway_graph
    return kegg_pathways


def find_in_which_pathway(target_gene_list: List[str], kegg_pathways: Dict[str, nx.DiGraph]) -> Dict[str, List[str]]:
    target_genes = set(target_gene_list)
    pathways_with_target_genes = dict()
    for k_id, pathway_g in kegg_pathways.items():
        for gene in target_genes:  
            if pathway_g.has_node(gene):
                if k_id in pathways_with_target_genes:
                    pathways_with_target_genes[k_id] += [gene]
                else:
                    pathways_with_target_genes[k_id] = [gene]
    return pathways_with_target_genes


def list_optional_nodes(pathway_graph: nx.DiGraph, node_list: List[str]):
    optional_nodes = []
    for node_name in node_list:
        is_optional = pathway_graph.nodes[node_name]["is_optional"]
        if is_optional:
            optional_nodes.append(node_name)
    return optional_nodes


def compute_completion(pathway_graph: nx.DiGraph, target_genes: List[str]) -> Enrichment:
    shortest_path_through_nodes = find_shortest_path_through(pathway_graph, target_genes)
    full_pathway_li = shortest_path_through_nodes[1:-1]
    completion = round(len(target_genes) / len(full_pathway_li),3)
    optional_genes = list_optional_nodes(pathway_graph, shortest_path_through_nodes)
    res = {"completion": completion, "present_genes": target_genes, "pathway": full_pathway_li, "optional": optional_genes}
    return res


def sort_nodes(in_graph: nx.DiGraph, in_nodes: List[str]) -> List[str]:
    #sorted_g = list(nx.topological_sort(in_graph))
    sorted_g = list(in_graph.nodes)
    node_ids = [sorted_g.index(node) for node in in_nodes]
    sorted_nodes = [node for n_id, node in sorted(zip(node_ids, in_nodes))]
    return sorted_nodes


def compute_completion_of_all_pathways(kegg_pathways: Dict[str, nx.DiGraph], pathways_with_target_genes: Dict[str, List[str]]) -> Dict[str, Enrichment]:
    completion_res = dict()
    for k_id, pathway_g in kegg_pathways.items():
        if k_id in pathways_with_target_genes:
            target_genes = sort_nodes(pathway_g, pathways_with_target_genes[k_id])
            print(k_id)
            completion_info = compute_completion(pathway_g, target_genes)
            completion_res[k_id] = completion_info
        else:
            completion_res[k_id] = {"completion": 0.0, "present_genes": [], "pathway": [], "optional": []}
    return completion_res

#################################


def mark_chosen_path_in_graph(pathway_graph: nx.DiGraph, path_g: List[str], target_genes: List[str]):
    marked_graph = pathway_graph.copy()
    is_path_attr = dict()
    is_present_attr = dict()
    for node in pathway_graph.nodes():
        is_path_attr[node] = node in path_g
        is_present_attr[node] = node in target_genes
    nx.set_node_attributes(marked_graph, is_path_attr, "is_path")
    nx.set_node_attributes(marked_graph, is_present_attr, "is_present")
    return marked_graph


def add_attributes_to_all_graphs(pathway_graphs: Dict[str, nx.DiGraph], completion_all_pathways: Dict[str, Enrichment]) -> Dict[str, nx.DiGraph]:
    marked_graphs = dict()
    for g_name in pathway_graphs:
        graph = pathway_graphs[g_name]
        present_genes = completion_all_pathways[g_name]["present_genes"]
        if len(present_genes) == 0:
            continue
        path_g = completion_all_pathways[g_name]["pathway"]
        marked_path = mark_chosen_path_in_graph(graph, path_g, present_genes)
        marked_graphs[g_name] = marked_path
    return marked_graphs


#################################


def gen_line_reader(file_path):
    for line in open(file_path, "r"):
        yield line


def KEGG_module_reader(KEGG_module_file_path) -> Dict[str, List[str]]:
    """
    Output a dict of module_id+name: str kos
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
        if line.startswith("Module:"):
            name = line.partition("Module:")[2].strip()
            if name not in KEGG_dict:
                KEGG_dict[name] = []
        if line.startswith("Definition:"):
            KEGG_dict[name].append(line.partition("Definition:")[2].strip())
    return KEGG_dict


def eggnog_parser(eggnog_path) -> List[str]:
    """
    Output list of ids KO ids
    Parses an eggnog.annotations output file. Returns a list of all the found k terms.
    """
    out_list = []
    for line in gen_line_reader(eggnog_path):
        if line.startswith("#"):
            continue
        ko = line.split("\t")[11]
        # !!!!A comma means there is multiple ko terms. BLASTkoala appears to only save the first one.
        # This script will include both ko terms
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


def convert_completion_dict_to_df(completion_dict):
    graph_res = dict()
    for k, v in completion_dict.items():
        if v["completion"] == 0.0:
            continue
        else:
            new_k = k.split(" ")[0]
            present_genes_str = ",".join(v["present_genes"])
            pathway_str = ",".join(v["pathway"])
            optional_genes_str = ",".join(v["optional"])
            new_v_dict = {"completion": v["completion"], "present_genes": present_genes_str, "description": k,
                          "pathway": pathway_str, "optional_genes": optional_genes_str}
            graph_res[new_k] = new_v_dict
    completion_df = pd.DataFrame.from_dict(graph_res, orient="index")
    completion_df.index.names = ["module"]
    return completion_df


def parseargs() -> Tuple[Path, Path, Path]:
    parser = argparse.ArgumentParser(description="Estimate completion of the modules")
    parser.add_argument("-m", help="Path to the module file", required=True, dest="mod_path", type=Path)
    parser.add_argument("-e", help="Path to the eggnog file", required=True, dest="eggnogfile_path", type=Path)
    parser.add_argument("-o", help="Path to the output directory", required=True, dest="out_dir", type=Path)
    args = parser.parse_args()
    return args.mod_path, args.eggnogfile_path, args.out_dir


def resolve_rel_path_list(path_list: List[Path]):
    res_paths = []
    for p in path_list:
        res_paths.append(p.resolve())
    return res_paths


def main():
    mod_path, eggnogfile_path, out_dir = parseargs()
    mod_path, eggnogfile_path, out_dir = resolve_rel_path_list([mod_path, eggnogfile_path, out_dir])
    out_dir.mkdir(exist_ok=True)
    
    print(f"Parsing the module list from {mod_path}")
    kegg_modules = KEGG_module_reader(mod_path)
    kegg_pathways = process_all_kegg_modules_to_pathways(kegg_modules)

    print(f"Parsing the eggnog file {eggnogfile_path}")
    eggnog_list = eggnog_parser(eggnogfile_path)

    print("Estimating completion")
    pathways_with_target_genes = find_in_which_pathway(eggnog_list, kegg_pathways)
    completion_of_all_pathways = compute_completion_of_all_pathways(kegg_pathways, pathways_with_target_genes)

    print("Writing to a table")
    completion_df = convert_completion_dict_to_df(completion_of_all_pathways)
    

    marked_graphs = add_attributes_to_all_graphs(kegg_pathways, completion_of_all_pathways)
    graph_dir = out_dir / "graphs"
    graph_dir.mkdir(exist_ok=True)
    print(f"Saving graph information to {graph_dir}")
    for name, graph in marked_graphs.items():
        mod_id = name[:6]
        out_file_path = graph_dir / f"{mod_id}_graph.gml"
        nx.write_gml(graph, out_file_path)

    out_path = out_dir / "modules.tsv"
    print(f"Saving results to {out_path}")
    completion_df.to_csv(out_path, sep="\t", index=True)
    print("Done")


if __name__ == "__main__":
    main()
