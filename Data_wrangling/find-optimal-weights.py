#!/usr/bin/env python3

import argparse
import subprocess
import dendropy
import pandas as pd
from pathlib import Path
import glob
import logging
import numpy as np
from scipy.optimize import minimize

# Setup logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Find optimal weights to maximize score for k=11, m=0.4.")
    parser.add_argument(
        "-k", "--kylo", required=True, help="Comma-separated k-mer lengths (e.g., 15,11,11)"
    )
    parser.add_argument(
        "-m", "--minifys", required=True, help="Comma-separated min sample fractions (e.g., 0.4,0.9)"
    )
    parser.add_argument(
        "--phylo", action="store_true", help="Run in phylo mode; default is network mode"
    )
    parser.add_argument(
        "--fastmode", action="store_true", help="Run in fast network mode (skip bootstrap and FastTree)"
    )
    args = parser.parse_args()

    if args.phylo and args.fastmode:
        parser.error("Cannot use --phylo and --fastmode together")

    return args

def get_fasta_files():
    """Get all *.fasta files in current directory."""
    fasta_files = glob.glob("*.fasta")
    if not fasta_files:
        raise FileNotFoundError("No *.fa files found in current directory")
    return [Path(f) for f in fasta_files]

def validate_alignment(aln_file):
    """Check if alignment file is valid and compute metrics, with PopArt-style segregating sites."""
    if not aln_file.is_file() or aln_file.stat().st_size == 0:
        return False, "Alignment file is empty or missing", 0, 0.0, 0.0, 0, 0
    
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
        return False, "No sequences in alignment", 0, 0.0, 0.0, 0, 0
    
    seq_count = len(sequences)
    if seq_count < 2:
        return False, "Too few sequences", seq_count, 0.0, 0.0, 0, 0
    
    seq_len = len(sequences[0])
    if not all(len(s) == seq_len for s in sequences):
        return False, "Sequences have unequal lengths", seq_count, 0.0, 0.0, 0, 0
    
    segregating_site = 0
    gap_count = 0
    n_count = 0
    total_bases = seq_count * seq_len
    valid_bases = {'A', 'T', 'C', 'G'}
    
    for i in range(seq_len):
        column = [s[i].upper() for s in sequences]
        if any(c not in valid_bases for c in column):
            gap_count += column.count("-")
            n_count += sum(1 for c in column if c not in valid_bases and c != "-")
            continue
        bases = set(column)
        if len(bases) > 1:
            segregating_site += 1
    
    gap_prop = gap_count / total_bases if total_bases > 0 else 0.0
    n_prop = n_count / total_bases if total_bases > 0 else 0.0
    total_nucleotide_site = seq_len  # Use alignment length
    
    if segregating_site == 0:
        return False, "No segregating sites", seq_count, gap_prop, n_prop, total_nucleotide_site, segregating_site
    
    return True, f"{seq_count} sequences, {segregating_site} sites", seq_count, gap_prop, n_prop, total_nucleotide_site, segregating_site

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
    filtered_seqs = [list(seq.upper()) for seq in sequences]
    
    positions_to_ignore = set()
    for i in range(seq_len):
        column = [seq[i] for seq in sequences]
        if '-' in column or 'N' in column:
            positions_to_ignore.add(i)
    
    haplotype_strings = [''.join(seq[i] for i in range(seq_len) if i not in positions_to_ignore) for seq in filtered_seqs]
    unique_haplotypes = len(set(haplotype_strings))
    return unique_haplotypes

def run_ska_build(k, fasta_files, outdir):
    """Run ska build for given k-mer length if .skf file doesn't exist."""
    skf_file = outdir / f"Combined-k{k}.skf"
    if skf_file.is_file() and skf_file.stat().st_size > 0:
        logging.info(f"Using existing {skf_file}")
        return skf_file
    cmd = ["ska", "build", "-k", str(k), "-o", str(skf_file)] + [str(f) for f in fasta_files]
    subprocess.run(cmd, check=True, capture_output=True)
    return skf_file

def run_ska_align(skf_file, m, outdir):
    """Run ska align for given min sample fraction if .aln file doesn't exist."""
    m_str = str(m).replace(".", "p")
    aln_file = outdir / f"Combined-k{skf_file.stem.split('-k')[1]}-m{m_str}.aln"
    if aln_file.is_file() and aln_file.stat().st_size > 0:
        logging.info(f"Using existing {aln_file}")
        return aln_file
    cmd = [
        "ska", "align", str(skf_file), "-m", str(m),
        "--filter", "no-filter", "--ambig-mask", "-o", str(aln_file)
    ]
    subprocess.run(cmd, check=True, capture_output=True)
    return aln_file

def run_phylogenetic_analysis(aln_file, k, m, outdir):
    """Run FastTree if .tree file doesn't exist."""
    m_str = str(m).replace(".", "p")
    tree_file = outdir / f"Combined-k{k}-m{m_str}.tree"
    if tree_file.is_file() and tree_file.stat().st_size > 0:
        logging.info(f"Using existing {tree_file} for k={k}, m={m}")
    else:
        logging.info(f"Generating {tree_file} for k={k}, m={m}")
        cmd = ["FastTree", "-gtr", "-nt", "-quote", "-quiet", "-boot", "100", str(aln_file)]
        try:
            with open(tree_file, "w") as f:
                subprocess.run(cmd, check=True, stdout=f, stderr=subprocess.PIPE, text=True)
        except subprocess.CalledProcessError as e:
            logging.error(f"FastTree failed for {tree_file}: {e.stderr.decode()}")
            return None, 0, e.stderr

    try:
        tree = dendropy.Tree.get(path=tree_file, schema="newick")
        bootstrap_values = [float(node.label) for node in tree.internal_nodes() if node.label and node.label.replace(".", "").isdigit()]
        bootstrap_score = sum(bootstrap_values) / len(bootstrap_values) if bootstrap_values else 0
    except Exception as e:
        logging.error(f"Failed to parse {tree_file}: {str(e)}")
        return tree_file, 0, str(e)

    return tree_file, bootstrap_score, ""

def compute_composite_score(df, w_bootstrap, w_haplotypes, w_nucleotide_sites, w_segregating_sites, w_gap_prop):
    """Compute a composite score based on normalized metrics."""
    def normalize(series, tolerance=0.0):
        min_val = series.min()
        max_val = series.max()
        if max_val > min_val:
            norm = (series - min_val) / (max_val - min_val)
            if tolerance > 0 and series.name in ["segregating_sites", "haplotype_count"]:
                norm = np.where(series >= max_val * (1 - tolerance), 1.0, norm)
            return norm
        return series / series
    
    score = 0
    if w_bootstrap != 0:
        norm_bootstrap = normalize(df["bootstrap_score"])
        score += w_bootstrap * norm_bootstrap
    if w_haplotypes != 0:
        norm_haplotypes = normalize(df["haplotype_count"], tolerance=0.02)  # 2% tolerance for haplotypes
        score += w_haplotypes * norm_haplotypes
    if w_nucleotide_sites != 0:
        norm_nucleotide_sites = normalize(df["total_nucleotide_sites"])
        score += w_nucleotide_sites * norm_nucleotide_sites
    if w_segregating_sites != 0:
        norm_segregating_sites = normalize(df["segregating_sites"], tolerance=0.05)  # 5% tolerance
        score += w_segregating_sites * norm_segregating_sites
    if w_gap_prop != 0:
        norm_gap_prop = normalize(df["gap_prop"])
        score += w_gap_prop * norm_gap_prop
    
    score = np.clip(score, 0, 1)
    return score

def optimize_weights(df, optimal_k=11, optimal_m=0.4, fastmode=False):
    """Optimize weights to maximize score for k=11, m=0.4 relative to others."""
    def objective(weights):
        w_b, w_h, w_n, w_s, w_g = weights
        df["temp_score"] = compute_composite_score(df, w_b, w_h, w_n, w_s, w_g)
        
        optimal_score = df[(df["k"] == optimal_k) & (df["m"] == optimal_m)]["temp_score"].iloc[0] if not df[(df["k"] == optimal_k) & (df["m"] == optimal_m)].empty else 0
        other_scores = df[(df["k"] != optimal_k) | (df["m"] != optimal_m)]["temp_score"]
        max_other_score = other_scores.max() if not other_scores.empty else 0
        optimal_row = df[(df["k"] == optimal_k) & (df["m"] == optimal_m)].iloc[0] if not df[(df["k"] == optimal_k) & (df["m"] == optimal_m)].empty else None
        
        if optimal_score == 0:
            return 1e10  # Large penalty if optimal pair missing
        
        # Stronger penalty for others scoring higher
        ranking_penalty = 0
        if max_other_score > optimal_score:
            ranking_penalty = 10 * (max_other_score - optimal_score)
        
        # Small gap proportion penalty
        gap_penalty = 0.05 * optimal_row["gap_prop"] if optimal_row is not None else 0
        
        # Log progress
        logging.debug(f"Objective: optimal_score={optimal_score:.4f}, max_other_score={max_other_score:.4f}, "
                      f"ranking_penalty={ranking_penalty:.4f}, gap_penalty={gap_penalty:.4f}, "
                      f"weights=[{w_b:.3f}, {w_h:.3f}, {w_n:.3f}, {w_s:.3f}, {w_g:.3f}]")
        
        return -(optimal_score - max_other_score) + ranking_penalty + gap_penalty
    
    # Initial weights and bounds
    initial_weights = [0.85, 0.05, 0.1, 0.0, -0.05]  # wb=0.45, wh=0.05, wn=0.1, ws=0.0, wg=-0.0
    bounds = [
        (0.80, 0.90),  # wb: 0.45 ± 0.05
        (0.00, 0.10),  # wh: 0.05 ± 0.05
        (0.05, 0.15),  # wn: 0.1 ± 0.05
        (0.0, 0.05), # ws: 0.0 ± 0.05
        (-0.15, 0.0)   # wg: -0.0 to 0.0 to penalize gaps
    ]
    constraints = [
        {"type": "eq", "fun": lambda w: w[0] + w[1] + w[2] + w[3] - 1},  # wb + wh + wn + ws = 1
    ]
    
    result = minimize(
        objective,
        initial_weights,
        method="SLSQP",
        bounds=bounds,
        constraints=constraints,
        options={"disp": True, "maxiter": 2000, "ftol": 1e-8}
    )
    
    if result.success:
        optimal_weights = result.x
        logging.info(f"Optimization weights sum (wb+wh+wn+ws): {sum(optimal_weights[:4]):.3f}, wg={optimal_weights[4]:.3f}")
    else:
        logging.warning(f"Optimization failed: {result.message}")
        logging.info(f"Final weights: {result.x}, Sum (wb+wh+wn+ws): {sum(result.x[:4]):.3f}, wg={result.x[4]:.3f}")
        optimal_weights = initial_weights  # Fallback to initial weights
    
    df["temp_score"] = compute_composite_score(df, *optimal_weights)
    optimal_score = df[(df["k"] == optimal_k) & (df["m"] == optimal_m)]["temp_score"].iloc[0] if not df[(df["k"] == optimal_k) & (df["m"] == optimal_m)].empty else 0
    other_scores = df[(df["k"] != optimal_k) | (df["m"] != optimal_m)]["temp_score"]
    avg_other_score = other_scores.mean() if not other_scores.empty else 0
    max_other_score = other_scores.max() if not other_scores.empty else 0
    max_score_row = df.loc[df["temp_score"].idxmax()]
    optimal_row = df[(df["k"] == optimal_k) & (df["m"] == optimal_m)].iloc[0] if not df[(df["k"] == optimal_k) & (df["m"] == optimal_m)].empty else None
    
    # Compute normalized metrics and contributions
    def normalize(series, tolerance=0.0):
        min_val = series.min()
        max_val = series.max()
        if max_val > min_val:
            norm = (series - min_val) / (max_val - min_val)
            if tolerance > 0 and series.name in ["segregating_sites", "haplotype_count"]:
                norm = np.where(series >= max_val * (1 - tolerance), 1.0, norm)
            return norm
        return series / series
    norm_metrics = {
        "norm_bootstrap": normalize(df["bootstrap_score"]),
        "norm_haplotypes": normalize(df["haplotype_count"], tolerance=0.02),
        "norm_nucleotide_sites": normalize(df["total_nucleotide_sites"]),
        "norm_segregating_sites": normalize(df["segregating_sites"], tolerance=0.05),
        "norm_gap_prop": normalize(df["gap_prop"])
    }
    optimal_idx = df[(df["k"] == optimal_k) & (df["m"] == optimal_m)].index[0] if not df[(df["k"] == optimal_k) & (df["m"] == optimal_m)].empty else -1
    max_idx = df["temp_score"].idxmax()
    
    # Compute score contributions
    optimal_contributions = {
        "bootstrap": optimal_weights[0] * norm_metrics["norm_bootstrap"][optimal_idx] if optimal_idx != -1 else 0,
        "haplotypes": optimal_weights[1] * norm_metrics["norm_haplotypes"][optimal_idx] if optimal_idx != -1 else 0,
        "nucleotide_sites": optimal_weights[2] * norm_metrics["norm_nucleotide_sites"][optimal_idx] if optimal_idx != -1 else 0,
        "segregating_sites": optimal_weights[3] * norm_metrics["norm_segregating_sites"][optimal_idx] if optimal_idx != -1 else 0,
        "gap_prop": optimal_weights[4] * norm_metrics["norm_gap_prop"][optimal_idx] if optimal_idx != -1 else 0
    }
    max_contributions = {
        "bootstrap": optimal_weights[0] * norm_metrics["norm_bootstrap"][max_idx],
        "haplotypes": optimal_weights[1] * norm_metrics["norm_haplotypes"][max_idx],
        "nucleotide_sites": optimal_weights[2] * norm_metrics["norm_nucleotide_sites"][max_idx],
        "segregating_sites": optimal_weights[3] * norm_metrics["norm_segregating_sites"][max_idx],
        "gap_prop": optimal_weights[4] * norm_metrics["norm_gap_prop"][max_idx]
    }
    
    # Add normalized metrics to DataFrame
    df["norm_bootstrap"] = norm_metrics["norm_bootstrap"]
    df["norm_haplotypes"] = norm_metrics["norm_haplotypes"]
    df["norm_nucleotide_sites"] = norm_metrics["norm_nucleotide_sites"]
    df["norm_segregating_sites"] = norm_metrics["norm_segregating_sites"]
    df["norm_gap_prop"] = norm_metrics["norm_gap_prop"]
    df["contrib_bootstrap"] = df["norm_bootstrap"] * optimal_weights[0]
    df["contrib_haplotypes"] = df["norm_haplotypes"] * optimal_weights[1]
    df["contrib_nucleotide_sites"] = df["norm_nucleotide_sites"] * optimal_weights[2]
    df["contrib_segregating_sites"] = df["norm_segregating_sites"] * optimal_weights[3]
    df["contrib_gap_prop"] = df["norm_gap_prop"] * optimal_weights[4]
    
    # Log metric ranges
    logging.info(f"Metric ranges: bootstrap_score=[{df['bootstrap_score'].min():.4f}, {df['bootstrap_score'].max():.4f}], "
                 f"haplotypes=[{df['haplotype_count'].min()}, {df['haplotype_count'].max()}], "
                 f"nucleotide_sites=[{df['total_nucleotide_sites'].min()}, {df['total_nucleotide_sites'].max()}], "
                 f"segregating_sites=[{df['segregating_sites'].min()}, {df['segregating_sites'].max()}], "
                 f"gap_prop=[{df['gap_prop'].min():.4f}, {df['gap_prop'].max():.4f}]")
    
    logging.info(f"Optimization result: wb={optimal_weights[0]:.3f}, wh={optimal_weights[1]:.3f}, wn={optimal_weights[2]:.3f}, ws={optimal_weights[3]:.3f}, wg={optimal_weights[4]:.3f}")
    logging.info(f"Optimal score (k={optimal_k}, m={optimal_m}): {optimal_score:.4f}, Avg other score: {avg_other_score:.4f}, Max other score: {max_other_score:.4f}")
    if optimal_row is not None:
        logging.info(f"Optimal pair metrics: k={optimal_k}, m={optimal_m:.2f}, score={optimal_score:.4f}, "
                     f"haplotypes={optimal_row['haplotype_count']}, "
                     f"nucleotide_sites={optimal_row['total_nucleotide_sites']}, "
                     f"segregating_sites={optimal_row['segregating_sites']}, "
                     f"gap_prop={optimal_row['gap_prop']:.4f}, "
                     f"bootstrap_score={optimal_row['bootstrap_score']:.4f}, "
                     f"norm_bootstrap={norm_metrics['norm_bootstrap'][optimal_idx]:.4f}, "
                     f"norm_haplotypes={norm_metrics['norm_haplotypes'][optimal_idx]:.4f}, "
                     f"norm_nucleotide_sites={norm_metrics['norm_nucleotide_sites'][optimal_idx]:.4f}, "
                     f"norm_segregating_sites={norm_metrics['norm_segregating_sites'][optimal_idx]:.4f}, "
                     f"norm_gap_prop={norm_metrics['norm_gap_prop'][optimal_idx]:.4f}, "
                     f"contrib_bootstrap={optimal_contributions['bootstrap']:.4f}, "
                     f"contrib_haplotypes={optimal_contributions['haplotypes']:.4f}, "
                     f"contrib_nucleotide_sites={optimal_contributions['nucleotide_sites']:.4f}, "
                     f"contrib_segregating_sites={optimal_contributions['segregating_sites']:.4f}, "
                     f"contrib_gap_prop={optimal_contributions['gap_prop']:.4f}")
    logging.info(f"Top-scoring pair: k={int(max_score_row['k'])}, m={max_score_row['m']:.2f}, score={max_score_row['temp_score']:.4f}, "
                 f"haplotypes={max_score_row['haplotype_count']}, "
                 f"nucleotide_sites={max_score_row['total_nucleotide_sites']}, "
                 f"segregating_sites={max_score_row['segregating_sites']}, "
                 f"gap_prop={max_score_row['gap_prop']:.4f}, "
                 f"bootstrap_score={max_score_row['bootstrap_score']:.4f}, "
                 f"norm_bootstrap={norm_metrics['norm_bootstrap'][max_idx]:.4f}, "
                 f"norm_haplotypes={norm_metrics['norm_haplotypes'][max_idx]:.4f}, "
                 f"norm_nucleotide_sites={norm_metrics['norm_nucleotide_sites'][max_idx]:.4f}, "
                 f"norm_segregating_sites={norm_metrics['norm_segregating_sites'][max_idx]:.4f}, "
                 f"norm_gap_prop={norm_metrics['norm_gap_prop'][max_idx]:.4f}, "
                 f"contrib_bootstrap={max_contributions['bootstrap']:.4f}, "
                 f"contrib_haplotypes={max_contributions['haplotypes']:.4f}, "
                 f"contrib_nucleotide_sites={max_contributions['nucleotide_sites']:.4f}, "
                 f"contrib_segregating_sites={max_contributions['segregating_sites']:.4f}, "
                 f"contrib_gap_prop={max_contributions['gap_prop']:.4f}")
    
    # Log top 5 scoring pairs
    top_5 = df.nlargest(5, "temp_score")[["k", "m", "temp_score", "haplotype_count", "total_nucleotide_sites", "segregating_sites", "gap_prop", "bootstrap_score"]]
    logging.info("Top 5 scoring pairs:\n" + top_5.to_string(index=False))
    
    # Warn if optimal pair has suboptimal metrics
    if optimal_row is not None:
        if optimal_row["segregating_sites"] < df["segregating_sites"].max() * 0.9:
            logging.warning(f"Optimal pair (k={optimal_k}, m={optimal_m}) has low segregating sites: {optimal_row['segregating_sites']} vs max {df['segregating_sites'].max()}")
        if optimal_row["gap_prop"] > df["gap_prop"].min() + 0.01:
            logging.warning(f"Optimal pair (k={optimal_k}, m={optimal_m}) has high gap proportion: {optimal_row['gap_prop']:.4f} vs min {df['gap_prop'].min():.4f}")
        if max_score_row["total_nucleotide_sites"] > 14500:
            logging.warning(f"Top scorer (k={int(max_score_row['k'])}, m={max_score_row['m']:.2f}) has high nucleotide sites: {max_score_row['total_nucleotide_sites']}")
    
    return optimal_weights

def main():
    args = parse_args()
    
    k_values = [int(k) for k in args.kylo.split(",")]
    m_values = [float(m) for m in args.minifys.split(",")]
    
    logging.info(f"Input parameters: k={k_values}, m={m_values}")
    
    if 11 not in k_values or 0.4 not in m_values:
        logging.error("k=11 and m=0.4 must be included in parameters")
        print("Error: Must include k=11 and m=0.4 in -k and -m arguments")
        return
    
    try:
        fasta_files = get_fasta_files()
        logging.info(f"Found {len(fasta_files)} *.fasta files")
    except FileNotFoundError as e:
        logging.error(str(e))
        return
    
    outdir = Path("output_repeats_masked")
    outdir.mkdir(exist_ok=True)
    
    results = []
    mode = "fastmode" if args.fastmode else ("phylo" if args.phylo else "network")
    
    for k in k_values:
        try:
            skf_file = run_ska_build(k, fasta_files, outdir)
        except subprocess.CalledProcessError as e:
            logging.error(f"ska build failed for k={k}: {e.stderr.decode()}")
            continue
        
        for m in m_values:
            try:
                aln_file = run_ska_align(skf_file, m, outdir)
                is_valid, message, seq_count, gap_prop, n_prop, total_nucleotide_site, segregating_site = validate_alignment(aln_file)
                logging.info(f"Alignment for k={k}, m={m}: {message}")
                if not is_valid:
                    logging.warning(f"Skipping k={k}, m={m}: {message}")
                    continue
                
                haplotype_count = compute_haplotype_count(aln_file)
                logging.info(f"Haplotype count for k={k}, m={m}: {haplotype_count}")
                
                if args.fastmode:
                    tree_file = ""
                    bootstrap_score = 0
                    error = ""
                else:
                    tree_file, bootstrap_score, error = run_phylogenetic_analysis(aln_file, k, m, outdir)
                    if tree_file is None:
                        logging.error(f"Phylogenetic analysis failed for k={k}, m={m}: {error}")
                        continue
                    logging.info(f"Bootstrap score for k={k}, m={m}: {bootstrap_score:.4f}")
                
                results.append({
                    "k": k,
                    "m": m,
                    "bootstrap_score": bootstrap_score,
                    "tree_file": str(tree_file),
                    "aln_file": str(aln_file),
                    "seq_count": seq_count,
                    "gap_prop": gap_prop,
                    "n_prop": n_prop,
                    "haplotype_count": haplotype_count,
                    "total_nucleotide_sites": total_nucleotide_site,
                    "segregating_sites": segregating_site
                })
                logging.info(f"Processed k={k}, m={m}, Bootstrap Score={bootstrap_score:.4f}, Segregating Sites={segregating_site}, Gap Prop={gap_prop:.4f}")
                
            except Exception as e:
                logging.error(f"Error for k={k}, m={m}: {str(e)}")
                continue
    
    if not results:
        logging.error("No valid results obtained. Check input data or parameters.")
        print("\nNo valid parameter combinations produced usable results.")
        print("Suggestions:")
        print("- Ensure k=11, m=0.4 is included in parameters.")
        print("- Lower -m (e.g., 0.4, 0.9) or increase -k (e.g., 15, 11, 11).")
        print("- Verify *.fasta files contain one sequence each.")
        return
    
    df = pd.DataFrame(results)
    
    if not ((df["k"] == 11) & (df["m"] == 0.4)).any():
        logging.error("No results for k=11, m=0.4. Ensure these parameters produce a valid alignment.")
        print("\nError: No valid alignment for k=11, m=0.4. Try lowering m or adjusting k.")
        return
    
    optimal_weights = optimize_weights(df, fastmode=args.fastmode)
    w_b, w_h, w_n, w_s, w_g = optimal_weights
    
    df["score"] = compute_composite_score(df, w_b, w_h, w_n, w_s, w_g)
    optimal_score = df[(df["k"] == 11) & (df["m"] == 0.4)]["score"].iloc[0] if not df[(df["k"] == 11) & (df["m"] == 0.4)].empty else 0
    max_score_row = df.loc[df["score"].idxmax()]
    optimal_row = df[(df["k"] == 11) & (df["m"] == 0.4)].iloc[0] if not df[(df["k"] == 11) & (df["m"] == 0.4)].empty else None
    
    df.to_csv(f"weight_optimization_results_{mode}.csv", index=False)
    print(f"\nOptimized Weights for {mode.capitalize()} Mode:")
    print(f"wb (bootstrap) = {w_b:.3f}")
    print(f"wh (haplotypes) = {w_h:.3f}")
    print(f"wn (nucleotide sites) = {w_n:.3f}")
    print(f"ws (segregating sites) = {w_s:.3f}")
    print(f"wg (gap proportion) = {w_g:.3f}")
    print(f"Composite Score for k=11, m=0.4: {optimal_score:.4f}")
    if optimal_row is not None:
        print(f"Metrics for k=11, m=0.4: haplotypes={optimal_row['haplotype_count']}, "
              f"nucleotide_sites={optimal_row['total_nucleotide_sites']}, segregating_sites={optimal_row['segregating_sites']}, "
              f"gap_prop={optimal_row['gap_prop']:.4f}")
    if max_score_row["k"] != 11 or max_score_row["m"] != 0.4:
        print(f"Warning: k=11, m=0.4 is not the top scorer. Top pair: k={int(max_score_row['k'])}, m={max_score_row['m']:.2f}, score={max_score_row['score']:.4f}")
        print(f"Metrics for k={int(max_score_row['k'])}, m={max_score_row['m']:.2f}: haplotypes={max_score_row['haplotype_count']}, "
              f"nucleotide_sites={max_score_row['total_nucleotide_sites']}, segregating_sites={max_score_row['segregating_sites']}, "
              f"gap_prop={max_score_row['gap_prop']:.4f}")
    print(f"Full results saved to weight_optimization_results_{mode}.csv")
    print(f"Use these weights in new_optimize-ska-mtDNA.py with --wb {w_b:.3f} --wh {w_h:.3f} --wn {w_n:.3f} --ws {w_s:.3f} --wg {w_g:.3f}")

if __name__ == "__main__":
    main()