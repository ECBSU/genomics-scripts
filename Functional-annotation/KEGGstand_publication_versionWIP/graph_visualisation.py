from pathlib import Path
from typing import List, Tuple
import argparse

import networkx as nx


def extract_gml_files(input_dir: Path):
    gml_files = []
    listing = input_dir.iterdir()
    for filename in listing:
        if filename.suffix == ".gml":
            gml_files.append(filename)
    return gml_files


def transform_gml_to_graphs(gml_files: List[Path]):
    graphs = []
    for path in gml_files:
        g_name = path.name
        G = nx.read_gml(path)
        graphs.append((g_name, G))
    return graphs


def get_pydot_attrs(node):
    attrs = node.obj_dict["attributes"]
    return_attr = dict(
        is_optional = attrs["is_optional"] == "1", 
        is_path = attrs["is_path"] == "1",
        is_present = attrs["is_present"] == "1"  
    )
    return return_attr

c = dict(
    path_fill = "#00ff00",
    path_border = "#00ff00",
    optional_fill = "#bdffbd",
    optional_pborder = "#00ff00",
    default_fill = "#ffffff",
    default_border = "#d2d2d2"
)

def process_graph_pydot(graph):
    pdot = nx.drawing.nx_pydot.to_pydot(graph)
    termini = ["BEGIN", "END"]
    nodes_on_path = [] + termini

    for node in pdot.get_nodes():
        node_name = node.get_name()
        node_attrs = get_pydot_attrs(node)
        node.set_style("filled")
        # border colour
        if node_attrs["is_path"]:
            node.set_color(c["path_border"])
        elif node_name in termini:
            node.set_color(c["path_border"])
        else:
            node.set_color(c["default_border"])

        # fill colour
        if node_attrs["is_present"]:
            node.set_fillcolor(c["path_fill"])
            if node_attrs["is_optional"]:
                node.set_fillcolor(c["optional_fill"])
        else:
            node.set_fillcolor(c["default_fill"])

        node.set_shape("box")
        if node_attrs["is_path"]:
            nodes_on_path.append(node_name)
    
    nodes_on_path = set(nodes_on_path)

    for edge in pdot.get_edges():
        source, target = edge.obj_dict["points"]
        if source in nodes_on_path and target in nodes_on_path:
            edge.set_style("filled")
            edge.set_color(c["path_border"])
        else:
            edge.set_color(c["default_border"])
    return pdot


def write_pydot_graph(pdot, out_path, file_format):
    writers = {"png": pdot.write_png, "svg": pdot.write_svg}
    writer = writers[file_format]
    writer(out_path)


def parseargs() -> Tuple[str, Path, Path]:
    parser = argparse.ArgumentParser(description="Estimate completion of the modules")
    parser.add_argument("-f", help="File format png or svg", required=True, dest="file_format")
    parser.add_argument("-g", help="Path to the directory with graph files", required=True, dest="mod_path", type=Path)
    parser.add_argument("-o", help="Path to the output directory", required=True, dest="out_dir", type=Path)
    args = parser.parse_args()
    return args.file_format, args.graph_dir, args.out_dir


def resolve_rel_path_list(path_list: List[Path]):
    res_paths = []
    for p in path_list:
        res_paths.append(p.resolve())
    return res_paths


def main():
    file_format, graph_dir, out_dir = parseargs()
    graph_dir, out_dir = resolve_rel_path_list([graph_dir, out_dir])

    out_dir.mkdir(exist_ok=True)

    gml_files = extract_gml_files(graph_dir)
    graph_info = transform_gml_to_graphs(gml_files)
    
    num_graphs = len(graph_info)
    i = 0
    for name, graph in graph_info:
        print(f"Plotting graph {i}/{num_graphs} {name} ", end="\r")
        out_path = out_dir / f"{name}.{file_format}"
        pydot_graph = process_graph_pydot(graph)
        write_pydot_graph(pydot_graph, out_path, file_format)
        i += 1

if __name__ == "__main__":
    main()
