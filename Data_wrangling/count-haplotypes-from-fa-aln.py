import argparse
from Bio import SeqIO

def remove_gaps_and_count_haplotypes(fasta_file):
    # Read the sequences from the FASTA file
    sequences = list(SeqIO.parse(fasta_file, "fasta"))

    # Convert sequences to a list of strings
    sequence_strs = [str(record.seq) for record in sequences]

    # Transpose the alignment to columns
    transposed = list(zip(*sequence_strs))

    # Identify columns with gaps
    gap_columns = {i for i, column in enumerate(transposed) if '-' in column or 'N' in column}

    # Remove gap columns
    filtered_sequences = [''.join(base for i, base in enumerate(seq) if i not in gap_columns) for seq in sequence_strs]

    # Create a set to store unique sequences (haplotypes)
    unique_sequences = set(filtered_sequences)

    # The number of unique sequences is the number of haplotypes
    num_haplotypes = len(unique_sequences)
    return num_haplotypes

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Count the number of haplotypes in a FASTA file, ignoring columns with gaps.")
    parser.add_argument("fasta_file", help="Path to the input FASTA file")
    args = parser.parse_args()

    # Count the number of haplotypes
    num_haplotypes = remove_gaps_and_count_haplotypes(args.fasta_file)

    # Print the result
    print(f"Number of haplotypes: {num_haplotypes}")

if __name__ == "__main__":
    main()
