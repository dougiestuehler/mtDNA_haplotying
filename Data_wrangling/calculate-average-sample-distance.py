import sys

def calculate_averages(filename):
    # Initialize a dictionary to store sums and counts for each unique string
    value_sums = {}
    value_counts = {}

    # Open the file and process each row
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith("Sample1\tSample2\tDistance\tMismatches"):
                continue  # Skip header lines
            row = line.split()
            for i in [0, 1]:  # Check both columns 1 and 2
                key = row[i]
                value = float(row[2])
                if key in value_sums:
                    value_sums[key] += value
                    value_counts[key] += 1
                else:
                    value_sums[key] = value
                    value_counts[key] = 1

    # Calculate averages
    averages = {key: value_sums[key] / value_counts[key] for key in value_sums}
    return averages

def main():
    if len(sys.argv) != 2:
        print("Usage: python script.py <input_file>")
        sys.exit(1)
    
    filename = sys.argv[1]
    averages = calculate_averages(filename)
    
    # Sort the averages by value
    sorted_averages = sorted(averages.items(), key=lambda item: item[1], reverse=True)
    
    # Print the results
    for key, avg in sorted_averages:
        print(f"{key}: {avg:.2f}")

if __name__ == "__main__":
    main()
