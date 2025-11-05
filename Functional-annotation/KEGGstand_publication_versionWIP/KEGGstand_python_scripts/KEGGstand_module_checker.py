from typing import List, Dict, Union, Tuple
import networkx as nx
import pandas as pd
import argparse
from pathlib import Path

######################### Parsing the KEGG compressed graph representation

def parse_sequence(s, i=0, end_chars=None):
    """
    Parse a sequence of items from string s starting at index i, until one of end_chars is reached.
    Returns a tuple (elements_list, new_index).

    Each item in elements_list is one of:
      - a plain enzyme ID (string),
      - a list of alternatives (each alternative is itself a list of elements),
      - a tuple ("OPTIONAL", list_of_sequences), where each sequence is a list of elements.

    Rules:
      • Whitespace (spaces, tabs) and '+' both separate sequential sub‐elements.
      • A '(' kicks off a bracketed group of alternatives, parsed by parse_alternatives.
      • A ',' or ')' ends the current sequence (depending on end_chars).
      • A '-' indicates an optional subchain:
         – Skip the '-' itself,
         – If the next character is '(' → parse a bracket of alternatives → wrap those alternatives
           into one optional group,
         – Otherwise parse a bare ID and wrap it as a one‐step optional group.
    """
    if end_chars is None:
        end_chars = []
    elements = []
    length = len(s)

    while i < length:
        # 1) Skip whitespace
        if s[i].isspace():
            i += 1
            continue

        # 2) If we've reached an “end char” (',' or ')'), stop parsing this level
        if s[i] in end_chars:
            break

        # 3) A '+' just separates (complex) steps → skip
        if s[i] == '+':
            i += 1
            continue

        # 4) A '-' means “start of an optional subchain”
        if s[i] == '-':
            i += 1
            # Skip any whitespace after the '-'
            while i < length and s[i].isspace():
                i += 1

            # If next char is '(', parse the bracketed alternatives as the optional group
            if i < length and s[i] == '(':
                alts, i = parse_alternatives(s, i + 1)
                # alts is a list of alternative sequences (each itself a list of elements).
                # Wrap it into one OPTIONAL‐tuple:
                elements.append(("OPTIONAL", alts))

            else:
                # Otherwise parse a bare enzyme ID (until whitespace or special char)
                j = i
                while j < length and not s[j].isspace() and s[j] not in end_chars \
                        and s[j] not in ['+', '-', '(', ')', ',']:
                    j += 1
                token = s[i:j].strip()
                if token:
                    # Wrap that single ID into a one‐step optional group
                    elements.append(("OPTIONAL", [[token]]))
                i = j

            continue

        # 5) A '(' means “start a bracketed alternative group”
        if s[i] == '(':
            alts, i = parse_alternatives(s, i + 1)
            # alts is a list of alternative sequences; append it directly
            elements.append(alts)
            continue

        # 6) Otherwise it’s a bare enzyme ID (collect until whitespace or a special char)
        j = i
        while j < length and not s[j].isspace() and s[j] not in end_chars \
                and s[j] not in ['+', '-', '(', ')', ',']:
            j += 1
        token = s[i:j].strip()
        if token:
            elements.append(token)
        i = j

    return elements, i

def parse_alternatives(s, i):
    """
    Parse a comma‐separated list of alternatives from s starting at index i (just after '(').
    Stops when the matching ')' is found. Returns (list_of_alternatives, new_index).

    Each alternative is itself parsed by parse_sequence(...) up to the next ',' or ')'.
    """
    alternatives = []
    length = len(s)

    while i < length:
        # Parse one alternative as a sequence until we hit ',' or ')'
        alt_seq, i = parse_sequence(s, i, end_chars=[',', ')'])
        alternatives.append(alt_seq)

        if i >= length:
            break
        if s[i] == ',':
            i += 1  # skip comma and parse next alt
            continue
        if s[i] == ')':
            i += 1  # skip closing ')'
            break

    return alternatives, i

def merge_optional_elements(elements):
    """
    In a single list of parsed elements, merge consecutive ("OPTIONAL", ...) items
    into a single OPTIONAL group whose sequences are the concatenation of each pair.

    For example:
      elements = ["A", ("OPTIONAL", [["X"], ["Y"]]), ("OPTIONAL", [["Z"]]), "B"]
    → we want one OPTIONAL whose sequences = [ ["X","Z"], ["Y","Z"] ].
    """
    merged = []
    i = 0
    while i < len(elements):
        if isinstance(elements[i], tuple) and elements[i][0] == "OPTIONAL":
            combined_seqs = elements[i][1]  # this is a list of sequences (each seq is a list of elements)
            i += 1
            # As long as the next element is also an OPTIONAL, keep merging
            while i < len(elements) and isinstance(elements[i], tuple) and elements[i][0] == "OPTIONAL":
                next_seqs = elements[i][1]
                new_comb = []
                for seq1 in combined_seqs:
                    for seq2 in next_seqs:
                        new_comb.append(seq1 + seq2)
                combined_seqs = new_comb
                i += 1
            merged.append(("OPTIONAL", combined_seqs))
        else:
            merged.append(elements[i])
            i += 1
    return merged

def recursive_merge(elements):
    """
    Recursively traverse each element, merging OPTIONAL groups at every nesting level.

    If an element is:
      - a plain string → leave it,
      - a list of alternatives → recurse into each alternative (which is itself a list of elements),
      - an ("OPTIONAL", seqs) tuple → for each seq (a list of elements), recurse into that seq.

    After recursion, we also call merge_optional_elements(...) on the top‐level list.
    """
    new_elems = []
    for elem in elements:
        if isinstance(elem, str):
            new_elems.append(elem)

        elif isinstance(elem, list):
            # This is a bracketed‐alternative group: a list of alternative sequences
            merged_alts = []
            for alt in elem:   # alt is a list of elements
                processed_alt = recursive_merge(alt)
                processed_alt = merge_optional_elements(processed_alt)
                merged_alts.append(processed_alt)
            new_elems.append(merged_alts)

        else:
            # Must be ("OPTIONAL", seqs)
            tag, seqs = elem
            if tag != "OPTIONAL":
                raise ValueError(f"Unknown element type: {elem!r}")
            merged_seqs = []
            for seq in seqs:  # seq is itself a list of elements
                processed_seq = recursive_merge(seq)
                processed_seq = merge_optional_elements(processed_seq)
                merged_seqs.append(processed_seq)
            new_elems.append(("OPTIONAL", merged_seqs))

    return new_elems

def old_process_elements(elements, prev_ids, edges):
    """
    Given parsed elements (each of which may be):
       • a string (enzyme ID),
       • a list (alternatives),
       • a tuple ("OPTIONAL", list_of_sequences),
    and a list prev_ids (the “upstream” nodes to connect from),
    append (source,target) edges into `edges` and return the new leaves.

    Rules:
      1) If elem is a string "X":
         - For each p in prev_ids, append (p, "X")
         - Then current_prev = ["X"].

      2) If elem is a list (i.e. bracketed alternatives):
         - For each alternative sequence alt_seq (which itself is a list of elements),
           call process_elements(alt_seq, prev_ids, edges) to get leaves_for_that_alt.
         - Combine all those leaves into new current_prev.

      3) If elem is ("OPTIONAL", seqs):
         - We have two possibilities: “skip” or “take” the entire optional subchain.
         - Let leaves_skip = prev_ids
         - For each sequence seq in seqs:
             • call process_elements(seq, prev_ids, edges) to build edges for “taking” that chain.
             • collect the leaves from that taken path.
         - new current_prev = leaves_skip + (all leaves from taken‐paths)
    """
    current_prev = prev_ids

    for elem in elements:
        # Case 1: a plain ID
        if isinstance(elem, str):
            for p in current_prev:
                edges.append((p, elem))
            current_prev = [elem]

        # Case 2: bracketed‐alternatives
        elif isinstance(elem, list):
            all_leaves = []
            for alt_seq in elem:
                alt_leaves = process_elements(alt_seq, current_prev, edges)
                all_leaves.extend(alt_leaves)
            current_prev = all_leaves

        # Case 3: an optional‐chain marker
        else:
            tag, seqs = elem
            if tag != "OPTIONAL":
                raise ValueError(f"Unexpected element: {elem!r}")

            # 3a) If we skip the entire optional chain, leaves_skip = current_prev
            leaves_skip = list(current_prev)

            # 3b) If we take it, we must traverse each possible sub‐sequence in seqs
            leaves_taken = []
            for seq in seqs:
                taken_leaves = process_elements(seq, current_prev, edges)
                leaves_taken.extend(taken_leaves)

            # 3c) New “current_prev” is union of skip‐leaves and taken‐leaves
            current_prev = leaves_skip + leaves_taken

    return current_prev

def old_build_pathway_graph(pathway_str: str):
    """
    Top‐level function. Given a pathway_str that can contain:
      - plain IDs (e.g. "K00789"),
      - complexes joined by '+',
      - optional subchains joined by '-',
      - bracketed alternatives "(A,B,C,...)" with commas,

    this builds a directed acyclic graph (DAG) whose edges reflect:
      BEGIN → first step(s),
      step → next step(s),
      optional chains either skipped or fully taken,
      bracketed alternatives branched,
      final leaves → END.

    Returns (G, edges), where G is a networkx.DiGraph and edges is the Python list of (src, dst).
    """
    # 1) First parse into a raw “elements” list
    parsed, _ = parse_sequence(pathway_str, 0, end_chars=[])

    # 2) Recursively merge nested OPTIONAL groups
    parsed_rec = recursive_merge(parsed)

    # 3) Merge any consecutive OPTIONAL markers at this top level
    parsed_clean = merge_optional_elements(parsed_rec)

    # 4) Walk the cleaned elements, building edges
    edges = []
    leaves = process_elements(parsed_clean, ["BEGIN"], edges)

    # 5) Connect every final leaf to "END"
    for leaf in leaves:
        edges.append((leaf, "END"))

    # 6) Build the NetworkX graph
    G = nx.DiGraph()
    G.add_edges_from(edges)
    return G


def process_elements(elements, prev_ids, edges, node_optional, optional_context=False):
    """
    elements: parsed list where each item is:
        - string (enzyme ID)
        - list (bracketed alternatives; each alternative is a list of elements)
        - ("OPTIONAL", list_of_sequences) where each sequence is a list of elements
    prev_ids: list of upstream node ids
    edges: list to append (src, dst) tuples
    node_optional: dict mapping node_id -> bool (is_optional)
    optional_context: boolean, True if current call is within an optional chain

    Returns: list of leaf node ids after processing elements
    """
    current_prev = prev_ids

    for elem in elements:
        # Plain ID
        if isinstance(elem, str):
            # mark node optional if any time this element is seen under optional_context
            node_optional[elem] = node_optional.get(elem, False) or optional_context
            for p in current_prev:
                edges.append((p, elem))
            current_prev = [elem]

        # Bracketed alternatives (list)
        elif isinstance(elem, list):
            all_leaves = []
            for alt_seq in elem:
                # alt_seq processed with the same optional_context as the container
                alt_leaves = process_elements(alt_seq, current_prev, edges, node_optional, optional_context=optional_context)
                all_leaves.extend(alt_leaves)
            current_prev = all_leaves

        # OPTIONAL tuple
        else:
            tag, seqs = elem
            if tag != "OPTIONAL":
                raise ValueError(f"Unexpected element: {elem!r}")

            # skipping the optional chain: leaves_skip = current_prev
            leaves_skip = list(current_prev)

            # taking the optional chain: mark everything inside as optional_context=True
            leaves_taken = []
            for seq in seqs:
                taken_leaves = process_elements(seq, current_prev, edges, node_optional, optional_context=True)
                leaves_taken.extend(taken_leaves)

            # union skip + taken
            current_prev = leaves_skip + leaves_taken

    return current_prev


def build_pathway_graph(pathway_str):
    """
    Parse pathway_str and build graph and node optional labels.

    Returns: (G, edges, node_optional)
    - G: networkx.DiGraph with node attribute 'is_optional' set for each node
    - edges: list of (src, dst) tuples
    - node_optional: dict mapping node -> bool
    """
    # 1) parse and normalize optional groups
    parsed, _ = parse_sequence(pathway_str, 0, end_chars=[])
    parsed_rec = recursive_merge(parsed)
    parsed_clean = merge_optional_elements(parsed_rec)

    # 2) build edges while collecting node optional flags
    edges = []
    node_optional = {}

    node_optional["BEGIN"] = False
    leaves = process_elements(parsed_clean, ["BEGIN"], edges, node_optional, optional_context=False)

    # mark final leaves' optional flags if they were in optional_context already (process_elements handled it)
    # connect leaves to END
    for leaf in leaves:
        edges.append((leaf, "END"))

    node_optional["END"] = False

    # 3) create graph and set node attributes
    G = nx.DiGraph()
    G.add_edges_from(edges)

    # Ensure all nodes appear in node_optional map (if a node appears but wasn't assigned yet)
    for n in G.nodes():
        if n not in node_optional:
            node_optional[n] = False

    # Assign attribute to each node
    for n, is_opt in node_optional.items():
        if n in G:
            G.nodes[n]['is_optional'] = bool(is_opt)
    return G

################################# Shortest path search

def find_shortest_path_through(graph: nx.DiGraph, target_nodes: List[str], start='BEGIN', end='END') -> List[str]:
    # Build path segments
    full_path = []
    current = start
    for target in target_nodes:
        segment = nx.shortest_path(graph, current, target)
        if full_path:
            full_path += segment[1:]  # Avoid repeating current node
        else:
            full_path += segment
        current = target
    # Path from last target to END
    segment = nx.shortest_path(graph, current, end)
    full_path += segment[1:]
    return full_path


def process_all_kegg_modules_to_pathways(kegg_dict: Dict[str, List[str]]) -> Dict[str, nx.DiGraph]:
    kegg_pathways = dict()
    for k_id, pathway_kegg in kegg_dict.items():
        pathway_str = pathway_kegg[0]
        print(k_id, pathway_str)
        pathway_graph = build_pathway_graph(pathway_str)
        kegg_pathways[k_id] = pathway_graph
    return kegg_pathways

def find_in_which_pathway(target_gene_list: List[str], kegg_pathways: Dict[str, nx.DiGraph]) -> Dict[str, List[str]]:
    target_genes = set(target_gene_list)
    pathways_with_target_genes = dict()
    for gene in target_genes:
        for k_id, pathway_g in kegg_pathways.items():
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


def compute_completion(pathway_graph: nx.DiGraph, target_genes: List[str]) -> Dict[str, Union[str, List[str], List[str]]]:
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


def compute_completion_of_all_pathways(kegg_pathways: Dict[str, nx.DiGraph], pathways_with_target_genes: Dict[str, List[str]]) -> Dict[str, Union[str, List[str], List[str]]]:
    completion_res = dict()
    for k_id, pathway_g in kegg_pathways.items():
        if k_id in pathways_with_target_genes:
            print(k_id, pathways_with_target_genes[k_id])
            target_genes = sort_nodes(pathway_g, pathways_with_target_genes[k_id])
            completion_info = compute_completion(pathway_g, target_genes)
            completion_res[k_id] = completion_info
        else:
            completion_res[k_id] = {"completion": 0.0, "present_genes": [], "pathway": [], "optional": []}
    return completion_res



#################################


def gen_line_reader(file_path):
    for line in open(file_path, "r"):
        yield line


def KEGG_module_reader(KEGG_module_file_path) -> Dict[str, str]:
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
            print(name)
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


def parse_args() -> Tuple[Path, Path, Path]:
    parser = argparse.ArgumentParser(description="Estimate completion of the modules")
    parser.add_argument("-m", help="Path to the module file", required=True, dest="mod_path", type=Path)
    parser.add_argument("-e", help="Path to the eggnog file", required=True, dest="eggnogfile_path", type=Path)
    parser.add_argument("-o", help="Path to the output tsv table", required=True, dest="out_path", type=Path)
    args = parser.parse_args()
    return args.mod_path, args.eggnogfile_path, args.out_path


def resolve_rel_path_list(path_list: List[Path]):
    res_paths = []
    for p in path_list:
        res_paths.append(p.resolve())
    return res_paths


def main():
    mod_path, eggnogfile_path, out_path = parse_args()
    mod_path, eggnogfile_path, out_path = resolve_rel_path_list([mod_path, eggnogfile_path, out_path])
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

    print(f"Saving results to {out_path}")
    completion_df.to_csv(out_path, sep="\t", index=True)
    print("Done")


if __name__ == "__main__":
    main()
