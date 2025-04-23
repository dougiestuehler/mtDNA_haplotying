#!/usr/bin/env python3

import argparse
import subprocess
import dendropy
import pandas as pd
from pathlib import Path
import glob
import os
import logging
import matplotlib.pyplot as plt
import numpy as np

# Setup logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Optimize ska k-mer length and min sample fraction using all *.fasta files with FastTree.")
    parser.add_argument(
        "-k", "--kmer", required=True, help="Comma-separated k-mer lengths (e.g., 13,15,17,21)"
    )
    parser.add_argument(
        "-m", "--minfrac", required=True, help="Comma-separated min sample fractions (e.g., 0.1,0.2,0.3)"
    )
    return parser.parse_args()

def get_fasta_files():
    """Get all *.fasta files in current directory."""
    fasta_files = glob.glob("*.fa")
    if not fasta_files:
        raise FileNotFoundError("No *.fasta files found in current directory")
    return [Path(f) for f in fasta_files]

def validate_alignment(aln_file):
    """Check if alignment file is non-empty and contains variable sites."""
    if not aln_file.is_file() or aln_file.stat().st_size == 0:
        return False, "Alignment file is empty", 0, 0.0, 0.0
    
    sequences = []
    headers = []
    with open(aln_file) as f:
        seq = ""
        header = ""
        for line in f:
            if line.startswith(">"):
                if seq:
                    sequences.append(seq)
                    headers.append(header)
                header = line.strip()[1:]
                seq = ""
            else:
                seq += line.strip()
        if seq:
            sequences.append(seq)
            headers.append(header)
    
    if not sequences:
        return False, "No sequences in alignment", 0, 0.0, 0.0
    
    seq_count = len(sequences)
    if seq_count < 2:
        return False, "Too few sequences", seq_count, 0.0, 0.0
    
    seq_len = len(sequences[0])
    if not all(len(s) == seq_len for s in sequences):
        return False, "Sequences have unequal lengths", seq_count, 0.0, 0.0
    
    variable_sites = 0
    gap_count = 0
    n_count = 0
    total_bases = seq_count * seq_len
    for i in range(seq_len):
        column = [s[i].upper() for s in sequences]
        bases = {c for c in column if c not in {"-", "N"}}
        if len(bases) > 1:
            variable_sites += 1
        gap_count += column.count("-")
        n_count += column.count("N")
    
    gap_prop = gap_count / total_bases if total_bases > 0 else 0.0
    n_prop = n_count / total_bases if total_bases > 0 else 0.0
    
    if variable_sites == 0:
        return False, "No variable sites in alignment", seq_count, gap_prop, n_prop
    
    return True, f"{seq_count} sequences, {variable_sites} variable sites, {gap_prop:.2%} gaps, {n_prop:.2%} Ns", seq_count, gap_prop, n_prop

def compute_haplotype_count(aln_file):
    """Compute the number of haplotypes, ignoring positions with gaps or Ns."""
    sequences = []
    with open(aln_file) as f:
        seq = ""
        for line in f:
            if line.startswith(">"):
                if seq:
                    sequences.append(seq)
                seq = ""
            else:
                seq += line.strip()
        if seq:
            sequences.append(seq)
    
    if not sequences or len(sequences) < 1:
        return 0
    
    seq_len = len(sequences[0])
    # Create a filtered sequence by ignoring positions with gaps or Ns
    filtered_seqs = []
    for seq in sequences:
        filtered_seqs.append(list(seq.upper()))
    
    # Identify positions to ignore (those with gaps or Ns)
    positions_to_ignore = set()
    for i in range(seq_len):
        column = [seq[i] for seq in sequences]
        if '-' in column or 'N' in column:
            positions_to_ignore.add(i)
    
    # Create haplotype strings by excluding ignored positions
    haplotype_strings = []
    for seq in filtered_seqs:
        haplotype = ''.join(seq[i] for i in range(seq_len) if i not in positions_to_ignore)
        haplotype_strings.append(haplotype)
    
    # Count unique haplotypes
    unique_haplotypes = len(set(haplotype_strings))
    return unique_haplotypes

def run_ska_build(k, fasta_files, outdir):
    """Run ska build for given k-mer length with multiple FASTA files."""
    skf_file = outdir / f"Combined-k{k}.skf"
    cmd = ["ska", "build", "-k", str(k), "-o", str(skf_file)] + [str(f) for f in fasta_files]
    subprocess.run(cmd, check=True, capture_output=True)
    return skf_file

def run_ska_align(skf_file, m, outdir):
    """Run ska align for given min sample fraction."""
    m_str = str(m).replace(".", "p")
    aln_file = outdir / f"Combined-k{skf_file.stem.split('-k')[1]}-m{m_str}.aln"
    cmd = [
        "ska", "align", str(skf_file), "-m", str(m),
        "--filter", "no-filter", "--ambig-mask", "-o", str(aln_file)
    ]
    subprocess.run(cmd, check=True, capture_output=True)
    return aln_file

def run_fasttree(aln_file, outdir):
    """Run FastTree to build phylogenetic tree with bootstrap."""
    tree_file = outdir / f"{aln_file.stem}.tree"
    cmd = ["FastTree", "-gtr", "-nt", "-quote", "-quiet", "-boot", "100"]
    cmd.append(str(aln_file))
    try:
        with open(tree_file, "w") as f:
            result = subprocess.run(cmd, check=True, stdout=f, stderr=subprocess.PIPE, text=True)
        return tree_file, ""
    except subprocess.CalledProcessError as e:
        return None, e.stderr

def compute_bootstrap_support(tree_file):
    """Compute average bootstrap support from FastTree output."""
    try:
        tree = dendropy.Tree.get(path=tree_file, schema="newick")
        bootstrap_values = []
        for node in tree.internal_nodes():
            if node.label and node.label.replace(".", "").isdigit():
                bootstrap_values.append(float(node.label))
        return sum(bootstrap_values) / len(bootstrap_values) if bootstrap_values else 0
    except Exception as e:
        return 0, str(e)

def plot_heatmap(results_df, outdir):
    """Plot a heatmap of k-mer length vs. min sample fraction with bootstrap scores and haplotype counts."""
    # Pivot the data for the heatmap
    all_k = sorted(results_df["k"].unique())
    all_m = sorted(results_df["m"].unique())
    pivot_score = results_df.pivot(index="m", columns="k", values="score").reindex(index=all_m, columns=all_k, fill_value=0)
    pivot_haplo = results_df.pivot(index="m", columns="k", values="haplotype_count").reindex(index=all_m, columns=all_k, fill_value=0)
    
    plt.figure(figsize=(12, 8))
    plt.imshow(pivot_score, cmap="YlOrRd", aspect="auto", interpolation="nearest")
    plt.colorbar(label="Average Bootstrap Support")
    plt.xticks(ticks=range(len(all_k)), labels=all_k, rotation=45)
    plt.yticks(ticks=range(len(all_m)), labels=[f"{m:.1f}" for m in all_m])
    plt.xlabel("k-mer Length")
    plt.ylabel("Min Sample Fraction")
    plt.title("Heatmap of Bootstrap Support by k-mer Length and Min Sample Fraction")
    
    # Annotate each cell with haplotype count
    for i in range(len(all_m)):
        for j in range(len(all_k)):
            haplo_count = pivot_haplo.iloc[i, j]
            if not np.isnan(haplo_count):  # Only annotate if data exists
                plt.text(j, i, f"{int(haplo_count)}", ha="center", va="center", color="black", fontsize=10)
    
    plt.tight_layout()
    plt.savefig(outdir / "bootstrap_heatmap.png")
    plt.close()

def main():
    args = parse_args()
    
    # Parse k-mer lengths and min fractions
    k_values = [int(k) for k in args.kmer.split(",")]
    m_values = [float(m) for m in args.minfrac.split(",")]
    
    # Get FASTA files
    try:
        fasta_files = get_fasta_files()
        logging.info(f"Found {len(fasta_files)} *.fasta files")
    except FileNotFoundError as e:
        logging.error(str(e))
        return
    
    # Create output directory
    outdir = Path("output")
    outdir.mkdir(exist_ok=True)
    
    results = []
    
    # Iterate over parameter grid
    for k in k_values:
        try:
            # Run ska build with all FASTA files
            skf_file = run_ska_build(k, fasta_files, outdir)
        except subprocess.CalledProcessError as e:
            logging.error(f"ska build failed for k={k}: {e.stderr.decode()}")
            continue
        
        for m in m_values:
            try:
                # Run ska align
                aln_file = run_ska_align(skf_file, m, outdir)
                
                # Validate alignment
                is_valid, message, seq_count, gap_prop, n_prop = validate_alignment(aln_file)
                logging.info(f"Alignment for k={k}, m={m}: {message}")
                if not is_valid:
                    logging.warning(f"Skipping k={k}, m={m}: {message}")
                    continue
                
                # Compute haplotype count
                haplotype_count = compute_haplotype_count(aln_file)
                logging.info(f"Haplotype count for k={k}, m={m}: {haplotype_count}")
                
                # Run FastTree
                tree_file, error = run_fasttree(aln_file, outdir)
                if tree_file is None:
                    logging.error(f"FastTree failed for k={k}, m={m}: {error}")
                    continue
                
                # Score the tree
                score = compute_bootstrap_support(tree_file)
                
                results.append({
                    "k": k,
                    "m": m,
                    "score": score,
                    "tree_file": str(tree_file),
                    "aln_file": str(aln_file),
                    "seq_count": seq_count,
                    "gap_prop": gap_prop,
                    "n_prop": n_prop,
                    "haplotype_count": haplotype_count
                })
                logging.info(f"Processed k={k}, m={m}, Score={score:.4f}")
                
            except Exception as e:
                logging.error(f"Error for k={k}, m={m}: {str(e)}")
                continue
    
    # Check if any results were obtained
    if not results:
        logging.error("No valid results obtained. Check input data or parameters (e.g., lower m or increase k).")
        print("\nNo valid parameter combinations produced usable trees.")
        print("Suggestions:")
        print("- Lower -m (e.g., 0.1, 0.2) to retain more k-mers.")
        print("- Increase -k (e.g., 13, 15, 17) for specificity.")
        print("- Verify *.fasta files contain one sequence each with valid nucleotides.")
        print("- Check output/Combined-k11-m0p9.aln for formatting issues.")
        return
    
    # Convert results to DataFrame
    df = pd.DataFrame(results)
    
    # Sort all results by score in descending order
    df = df.sort_values(by="score", ascending=False)
    
    # Print all results in order of score
    print("\nAll parameter combinations (Metric: Average Bootstrap Support)")
    print("-" * 50)
    for _, row in df.iterrows():
        print(f"k={int(row['k'])}, m={row['m']:.2f}, Score={row['score']:.4f}, Haplotypes={row['haplotype_count']}, Tree={row['tree_file']}")
    
    # Save all results to CSV
    df.to_csv("results.csv", index=False)
    print("\nAll results saved to results.csv")
    
    # Plot heatmap
    plot_heatmap(df, outdir)
    print(f"Heatmap saved as output/bootstrap_heatmap.png")
    print(f"Intermediate files saved in {outdir}")

if __name__ == "__main__":
    main()