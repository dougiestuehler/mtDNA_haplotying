import pandas as pd
import sys

def calculate_percentages(gff_file):
    # Read the GFF3 file into a pandas DataFrame
    df = pd.read_csv(gff_file, sep="\t", comment='#', header=None, 
                     names=["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"])

    # Calculate the length of each feature
    df["length"] = df["end"] - df["start"] + 1

    # Calculate total length of the sequence
    total_length = df[df["type"] == "region"]["end"].iloc[0]

    # Sum lengths of different feature types
    protein_coding_length = df[df["type"] == "CDS"]["length"].sum()
    rrna_length = df[df["type"] == "rRNA"]["length"].sum()
    trna_length = df[df["type"] == "tRNA"]["length"].sum()

    # Non-coding length is total length minus the sum of the other lengths
    non_coding_length = total_length - (protein_coding_length + rrna_length + trna_length)

    # Calculate percentages
    protein_coding_percentage = (protein_coding_length / total_length) * 100
    rrna_percentage = (rrna_length / total_length) * 100
    trna_percentage = (trna_length / total_length) * 100
    non_coding_percentage = (non_coding_length / total_length) * 100

    print(f"Protein-coding: {protein_coding_percentage:.2f}%")
    print(f"rRNA: {rrna_percentage:.2f}%")
    print(f"tRNA: {trna_percentage:.2f}%")
    print(f"Non-coding: {non_coding_percentage:.2f}%")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <input_gff3_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    calculate_percentages(input_file)
