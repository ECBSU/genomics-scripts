import re
import json
from typing import Dict, List, Union
from pathlib import Path
import argparse


class Node:
    def __init__(self, indented_line):
        self.children = []
        self.level = len(indented_line) - len(indented_line.lstrip())
        self.text = indented_line.strip()

    def add_children(self, nodes):
        childlevel = nodes[0].level
        while nodes:
            node = nodes.pop(0)
            if node.level == childlevel: # add node as a child
                self.children.append(node)
            elif node.level > childlevel: # add nodes as grandchildren of the last child
                nodes.insert(0,node)
                self.children[-1].add_children(nodes)
            elif node.level <= self.level: # this node is a sibling, no more children
                nodes.insert(0,node)
                return

    def as_dict(self):
        if len(self.children) > 1:
            return {self.text: [node.as_dict() for node in self.children]}
        elif len(self.children) == 1:
            return {self.text: self.children[0].as_dict()}
        else:
            return self.text


def read_k_term_lines(path):
    lines = []
    with open(path, "r") as s:
        lines = s.readlines()
    return lines


k_term = {
    "ENTRY": "",
    "SYMBOL": "",
    "NAME": "",
    "PATHWAY": [],
    "REACTION": [],
    "BRITE": {}
}


def convert_k_terms_lines_to_list_of_dict(lines) -> List[Dict[str, Union[str, dict]]]:
    entries = []
    read_this_part = slice(12, None)
    num_lines = len(lines)
    for i in range(0, num_lines):
        line = lines[i]
        if i == num_lines-1:
            next_line = "STOP"
        else:
            next_line = lines[i+1]

        if line.startswith("ENTRY"):
            pathway_started = False
            reaction_started = False
            brite_started = False

            entry_id = re.search(r"K\d{5}", line).group(0)
            this_entry_id = entry_id
            this_entry = k_term.copy()
            this_entry["ENTRY"] = this_entry_id
            entries.append(this_entry)
        else:
            if line.startswith("SYMBOL"):
                sym_info = line[read_this_part].strip()
                this_entry["SYMBOL"] = sym_info
            elif line.startswith("NAME"):
                name_info = line[read_this_part].strip()
                this_entry["NAME"] = name_info

            elif line.startswith("PATHWAY"):
                pathway_started = True
                this_entry_pathway = [line[read_this_part].strip()]
            elif pathway_started:
                pathway_info = line[read_this_part].strip()
                this_entry_pathway.append(pathway_info)
                if re.match(r"^\S", next_line):
                    pathway_started = False
                    this_entry["PATHWAY"] = this_entry_pathway

            elif line.startswith("REACTION"):
                reaction_started = True
                this_entry_reaction = [line[read_this_part].strip()]
            elif reaction_started:
                reaction_info = line[read_this_part].strip()
                this_entry_reaction.append(reaction_info)
                if re.match(r"^\S", next_line):
                    reaction_started = False
                    this_entry["REACTION"] = this_entry_reaction

            elif line.startswith("BRITE"):
                root = Node("root")
                brite_started = True
                brite_lines = [line[read_this_part]]
            elif brite_started:
                line_info_with_spaces = line[read_this_part]
                brite_lines.append(line_info_with_spaces)
                if re.match(r"^\S", next_line):
                    brite_started = False
                    root.add_children([Node(bl) for bl in brite_lines])
                    this_entry["BRITE"] = root.as_dict()["root"]
    return entries


def save_to_json(out_path, entries):
    with open(out_path, "w") as s:
        json.dump(entries, s, sort_keys=False, indent=4)


def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-k", help="Path to k terms file", type=Path, dest="k_terms_file_path")
    parser.add_argument("-o", help="Path to output json", type=Path, dest="out_path")
    args = parser.parse_args()
    return args.k_terms_file_path, args.out_path


def resolve_rel_path_list(path_list: List[Path]):
    res_paths = []
    for p in path_list:
        res_paths.append(p.resolve())
    return res_paths


def main():
    k_terms_file_path, out_path = parseargs()
    k_terms_file_path, out_path = resolve_rel_path_list([k_terms_file_path, out_path])
    k_term_lines = read_k_term_lines(k_terms_file_path)
    k_term_list_of_dict = convert_k_terms_lines_to_list_of_dict(k_term_lines)
    save_to_json(out_path, k_term_list_of_dict)


if __name__ == "__main__":
    main()


# example for parsing
# indented_text = \
# """
# apple
#     colours
#         red
#         yellow
#         green
#     type
#         granny smith
#     price
#         0.10
# 123
#     456
#         789
# """
#
# root = Node('root')
# root.add_children([Node(line) for line in indented_text.splitlines() if line.strip()])
# d = root.as_dict()['root']
# print(d)
