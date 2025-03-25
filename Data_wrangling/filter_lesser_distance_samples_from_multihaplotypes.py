import csv
import argparse

def load_file2(file2_path):
    """Load file2 into a dictionary where key is sample_name and value is average_distance."""
    file2_data = {}
    with open(file2_path, 'r') as file2:
        reader = csv.reader(file2, delimiter='\t')
        for row in reader:
            sample_name, avg_distance = row[0], float(row[1])
            file2_data[sample_name] = avg_distance
    return file2_data

def process_file1(file1_path, file2_data):
    """Process file1 and return the sample with the greatest distance for each haplogroup."""
    samples_to_keep = set()
    samples_to_remove = set()
    
    with open(file1_path, 'r') as file1:
        for line in file1:
            # Each line in file1 is a comma-separated list of sample names
            samples = line.strip().split(',')
            
            # Find the sample with the greatest distance in the haplogroup
            max_sample = max(samples, key=lambda s: file2_data.get(s, -float('inf')))
            samples_to_keep.add(max_sample)
            
            # Add the lesser distance samples to the remove list
            samples_to_remove.update(set(samples) - {max_sample})
    
    return samples_to_keep, samples_to_remove

def filter_file2(file2_path, samples_to_remove, output_path):
    """Filter file2 and write the results without the lesser distance samples from file1."""
    with open(file2_path, 'r') as file2, open(output_path, 'w', newline='') as output_file:
        reader = csv.reader(file2, delimiter='\t')
        writer = csv.writer(output_file, delimiter='\t')
        
        for row in reader:
            sample_name = row[0]
            # Write the row to output if it's not in the remove list
            if sample_name not in samples_to_remove:
                writer.writerow(row)

def main():
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(description='Filter file2 based on file1 haplogroups.')
    parser.add_argument('file1', help='Path to file1 (haplogroup sample names)')
    parser.add_argument('file2', help='Path to file2 (sample names and distances)')
    parser.add_argument('output', help='Path to the output file')
    
    args = parser.parse_args()
    
    # Load file2 and file1, and filter file2 based on the greatest distance samples
    file2_data = load_file2(args.file2)
    _, samples_to_remove = process_file1(args.file1, file2_data)
    filter_file2(args.file2, samples_to_remove, args.output)
    
    print(f"Filtered results have been saved to {args.output}")

if __name__ == "__main__":
    main()
