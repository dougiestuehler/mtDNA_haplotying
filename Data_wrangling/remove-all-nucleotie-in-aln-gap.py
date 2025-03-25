def read_fasta(filename):
    with open(filename, 'r') as file:
        sequences = []
        sequence = ''
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if sequence:
                    sequences.append(sequence)
                    sequence = ''
                sequences.append(line)
            else:
                sequence += line
        if sequence:
            sequences.append(sequence)
    return sequences

def write_fasta(filename, headers, sequences):
    with open(filename, 'w') as file:
        for header, sequence in zip(headers, sequences):
            file.write(header + '\n')
            for i in range(0, len(sequence), 80):  # Wrap sequence every 80 characters
                file.write(sequence[i:i+80] + '\n')

def remove_invalid_positions(sequences):
    sequence_length = len(sequences[0])
    invalid_positions = set()
    
    for seq in sequences:
        invalid_positions.update({i for i, char in enumerate(seq) if char == '-' or char == 'N'})
    
    invalid_positions = sorted(invalid_positions)
    
    new_sequences = []
    for seq in sequences:
        new_seq = ''.join([char for i, char in enumerate(seq) if i not in invalid_positions])
        new_sequences.append(new_seq)
    
    return new_sequences

def process_fasta(input_file, output_file):
    fasta_data = read_fasta(input_file)
    headers = [line for line in fasta_data if line.startswith('>')]
    sequences = [line for line in fasta_data if not line.startswith('>')]
    
    if not all(len(seq) == len(sequences[0]) for seq in sequences):
        raise ValueError("All sequences are not the same length")

    new_sequences = remove_invalid_positions(sequences)
    write_fasta(output_file, headers, new_sequences)

# Example usage
input_fasta = '/mnt/c/Bioinformatics/Heck_USDA_Pathology/Dcitri_mtDNA_haplotyping/Examining_the_effect_of_kmer_size/ska2/Combined_ACP/ska-align-k15-combined-out-p90-nofilter-justNs.aln'
output_fasta = '/mnt/c/Bioinformatics/Heck_USDA_Pathology/Dcitri_mtDNA_haplotyping/Examining_the_effect_of_kmer_size/ska2/Combined_ACP/ska-align-k15-combined-out-p90-nofilter-justNs-no-gap-nucl.fasta'
process_fasta(input_fasta, output_fasta)
