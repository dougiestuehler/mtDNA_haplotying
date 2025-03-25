import gzip
import sys
import os

def open_file(filename):
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')
    else:
        return open(filename, 'r')

def count_samples_with_alt_allele(gvcf_file, output_file):
    with open_file(gvcf_file) as file, open(output_file, 'w') as out:
        num_samples = 0  # Initialize num_samples

        for line in file:
            # Skip header lines
            if line.startswith('##'):
                continue
            elif line.startswith('#'):
                # Parse header line to get sample IDs
                sample_ids = line.strip().split('\t')[9:]
                num_samples = len(sample_ids)
                continue
            
            # Split the line into columns
            columns = line.strip().split('\t')
            
            # Extract SNP location (chromosome and position)
            chrom = columns[0]
            pos = columns[1]
            snp_location = f'"{pos}"'  # Surround position with quotes
            
            # Extract genotype information for each sample
            genotypes = columns[9:]
            
            # Count the number of samples with an alternate allele
            alt_count = 0
            for genotype in genotypes:
                # Genotype format is GT:...
                gt = genotype.split(':')[0]
                if '1' in gt:  # '1' represents the presence of an alternate allele
                    alt_count += 1
            
            # Calculate percentage of samples with alternate allele
            alt_percentage = (alt_count / num_samples) * 100 if num_samples > 0 else 0
            
            # Write SNP location, count, and percentage of samples with alternate allele to output file
            out.write(f"{chrom}\t{snp_location}\t{alt_count}\t{alt_percentage:.2f}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_gvcf_file> <output_file>")
        sys.exit(1)

    gvcf_file = sys.argv[1]
    output_file = sys.argv[2]
    
    if not os.path.isfile(gvcf_file):
        print(f"Error: File '{gvcf_file}' does not exist.")
        sys.exit(1)
    
    count_samples_with_alt_allele(gvcf_file, output_file)
