#! /bin/python3
from Bio import SeqIO
import argparse

def split_fasta(input_file, output_prefix, entries_per_file=100):
    records = list(SeqIO.parse(input_file, "fasta"))

    # Split the records into chunks of size entries_per_file
    chunks = [records[i:i + entries_per_file] for i in range(0, len(records), entries_per_file)]

    # Write each chunk to a separate output file
    for i, chunk in enumerate(chunks):
        output_file = f"{output_prefix}_{i + 1}.fasta"
        SeqIO.write(chunk, output_file, "fasta")
        print(f"File {output_file} created with {len(chunk)} entries.")
        print(f"You can use --array=0-{len(chunk) - 1}]")

if __name__ == "__main__":
    parser = argparse.ArgumentParser (
            prog='Split fasta',
            description='Split fasta file into smaller chunks')

    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output_dir')
    parser.add_argument('-n', '--entries',nargs='?',default='1000', type=int)
    args = parser.parse_args()
    
    split_fasta(args.input, args.output_dir, args.entries)
