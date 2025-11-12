#!/usr/bin/env python3
"""
Will query the KEGG API server to obtain all the KEGG module defintions, and writes them
to a text file

Usage:
(python) KEGG_db_generate.py -m [MODULE or KTERM] -o /path/to/db.txt
"""
import time
import re
import argparse
from pathlib import Path
from copy import deepcopy
from _io import TextIOWrapper
from enum import Enum
from typing import List, Dict, Tuple, Generator, Union

from Bio.KEGG import REST


def gen_line_reader(file_path: Path) -> Generator[str]:
    """
    Generator function that allows reading a text file line
    by line, without reading the full file into memory
    """
    for line in open(file_path, "r"):
        yield line


def find_processed_entries(file_path, search_prefix: str) -> List[str]:
    """
    Reads provided file and creates a list of all KEGG entries
    within.
    """
    entry_list = []
    for line in gen_line_reader(file_path):
        if line.startswith(search_prefix):
            entry_list.append(line.split()[1])
    return entry_list


def get_kegg_id_list(kegg_list_type: str) -> List[str]:
    """
    Connects to KEGG to request a list of all KEGG k terms.
    Parses out only the k term (removing the name) and stores
    it into a list.
    """
    id_list = []
    response: TextIOWrapper = REST.kegg_list(kegg_list_type, org=None)
    response_message = response.buffer.msg
    if response_message != "OK":
        raise ConnectionError(f"Problem with connecting to kegg, got response: {response_message}")
    for item in response:
        term_id = item.split()[0] # M00001\tGlycolysis (Embden-Meyerhof pathway), glucose => pyruvate\n -> M00001
        id_list.append(term_id)
    return id_list


class _KeggAccessMode:
    def __init__(self, kegg_list, search_prefix):
        self.kegg_list = kegg_list
        self.search_prefix = search_prefix

class KeggAccModes(Enum):
    MODULE = _KeggAccessMode(kegg_list="module", search_prefix="Module:")
    KTERM = _KeggAccessMode(kegg_list="ko", search_prefix="ENTRY:")


def check_if_item_already_present(file_path: Path, id_list: List[str], search_prefix: str) -> List[str]:
    """
    Check which k terms are already present in the database (if it exists).
    """
    # Remove any terms already in the database from the list, so they will not be downloaded again

    filtered_id_list = deepcopy(id_list)
    if file_path.is_file():
        for _id in find_processed_entries(file_path, search_prefix):
            if _id in filtered_id_list:
                filtered_id_list.remove(_id)
            else:
                print(f"error, found processed module that is not in list: {_id}")
    return filtered_id_list


def parse_entry(entry: TextIOWrapper) -> Dict[str, str]:
    """
    Parse entry into a dictionary -> {"NAME": "string contents"}
    example of entry https://rest.kegg.jp/get/M00993
    """
    entry_text = entry.buffer.read().decode("UTF8")
    entry_lines = entry_text.split("\n")

    entry_dict = dict()
    this_chapter = ""
    prev_chapter = ""
    for line in entry_lines:
        line_parts = line.split()
        if line_parts == [] or line.startswith("///"):
            break
        line_title = line_parts[0]
        if line.startswith(" "):
            entry_dict[prev_chapter] += line + "\n"
        else:
            line_text = re.sub(rf"^{line_title}\s+", "", line + "\n")
            this_chapter = line_title
            prev_chapter = line_title
            entry_dict[this_chapter] = line_text
    return entry_dict


def entry_dict_to_text(entry_dict: Dict[str, str]):
    txt = ""
    for ch_name, ch_text in entry_dict.items():
        chapter_str = ch_name.ljust(12) + ch_text
        txt += chapter_str
    return txt


def process_module_entry_info(mod_id: str, entry_dict: Dict[str, str]) -> str:
    name = f"Module: {mod_id} {entry_dict['NAME']}"  # 'Module: M00993 Dimethylsulfoniopropionate (DMSP) degradation, cleavage pathway, DMSP => propionyl-CoA\n'
    defin = f"Definition: {entry_dict['DEFINITION']}"  # 'Definition: K28072 K19745\n'
    cls = f"Class: {entry_dict['CLASS']}"  # 'Class: Pathway modules; Energy metabolism; Sulfur metabolism\n'
    res = name + defin + cls
    return res


def process_kterm_entry_info(entry_dict: Dict[str, str]) -> str:
    entry_dict_proc = entry_dict.copy()
    for chapter in ("DBLINKS", "GENES", "REFERENCE"):
        if chapter in entry_dict_proc:
            del entry_dict_proc[chapter]
    return entry_dict_to_text(entry_dict_proc)


def download_entry(_id: str) -> Union[TextIOWrapper, str]:
    for tries in range(1, 5):
        try:
            response: TextIOWrapper = REST.kegg_get(_id, option=None)
            response_message = response.buffer.msg
            if response_message == "OK":
                return response
            else:
                print(response_message)
                return "BAD"
        except Exception as e:
            t_delay = tries * 5
            print("Exception:", e)
            print(f"No connection, retry in {t_delay} seconds")
            time.sleep(t_delay)
    print(f"Connection failed {tries+1} consecutive times. Script will exit")
    return "BAD"


def process_entry(_id: str, entry_dict: Dict[str, str], mode_name: str) -> str:
    if mode_name == "MODULE":
        entry_str = process_module_entry_info(_id, entry_dict)
    elif mode_name == "KTERM":
        entry_str = process_kterm_entry_info(entry_dict)
    else:
        raise ValueError("Unknown mode")
    return entry_str


def write_entries(entries_list: List[str], out_file: Path):
    entries_str = "".join(entries_list)
    with open(out_file, "a") as s:
        s.write(entries_str)
    return


def download_and_process_and_save_entries(id_list: List[str], out_file: Path, mode_name: str):
    # Will download Kterms and save to disk every 100, and then wait 5 sec
    # If the server is refusing will wait 5 seconds and retry 5 times
    print("Downloading entries")
    proc_entries_list = []
    num_ids = len(id_list)
    for i, _id in enumerate(id_list):
        print(f"Entry: {i}/{num_ids}", end="\r")
        entry = download_entry(_id)
        if entry == "BAD":
            write_entries(proc_entries_list, out_file)
            raise ConnectionError("Problems with connection to the server")
        entry_dict = parse_entry(entry)
        processed_entry = process_entry(_id, entry_dict, mode_name)
        proc_entries_list.append(processed_entry)

        if len(proc_entries_list) == 100 or i == num_ids - 1:
            write_entries(proc_entries_list, out_file)
            proc_entries_list = []
            time.sleep(5)
    return


def parse_cmd_args() -> Tuple[Path, str]:
    parser = argparse.ArgumentParser(description="Downloads all KEGG modules utilizing the KEGG API server")
    parser.add_argument("-o", help="Path to output file", required=True, dest="out_file", type=Path)
    parser.add_argument("-m", help="Mode: which db to download: KTERM or MODULE", required=True, dest="mode", choices=["KTERM", "MODULE"])
    cmd_args = parser.parse_args()
    return cmd_args.out_file, cmd_args.mode


def main():
    out_file, mode_name = parse_cmd_args()
    print("Started")
    mode = KeggAccModes[mode_name]

    id_list = get_kegg_id_list(mode.value.kegg_list) # Acquire the list of k terms from KEGG
    filtered_id_list = check_if_item_already_present(out_file, id_list, mode.value.search_prefix)
    download_and_process_and_save_entries(filtered_id_list, out_file, mode.name)
    print("Finished")


if __name__ == "__main__":
    main()
