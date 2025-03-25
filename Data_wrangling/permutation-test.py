import numpy as np
import random

# Observed haplotype diversity values
h_old = 1.000
h_new = 0.985
n_samples = 17

# Observed difference
observed_diff = h_old - h_new

# Simulate haplotype distributions
# Since you only have the summarized values, we'll simulate haplotype frequencies 
# consistent with the given diversity values.

# Function to calculate Nei's haplotype diversity
def nei_diversity(freqs):
    n = sum(freqs)
    pi_sq = sum((f/n)**2 for f in freqs)
    return (n/(n-1)) * (1 - pi_sq)

# Simulate haplotype frequencies close to observed diversities
def simulate_haplotypes(target_diversity, n_individuals, n_haplotypes):
    for _ in range(10000):
        freqs = np.random.multinomial(n_individuals, np.ones(n_haplotypes)/n_haplotypes)
        h = nei_diversity(freqs)
        if np.abs(h - target_diversity) < 0.001:
            return freqs
    raise ValueError("Could not simulate suitable haplotypes")

# Simulate initial haplotype frequencies for each method
freqs_old = simulate_haplotypes(h_old, n_samples, 17)
freqs_new = simulate_haplotypes(h_new, n_samples, 15)

# Combine data
combined = np.concatenate([np.repeat(range(len(freqs_old)), freqs_old),
                           np.repeat(range(len(freqs_new)), freqs_new)])

# Permutation test
n_permutations = 10000
perm_diffs = []
for _ in range(n_permutations):
    np.random.shuffle(combined)
    perm_old = combined[:n_samples]
    perm_new = combined[n_samples:]

    h_perm_old = nei_diversity(np.bincount(perm_old))
    h_perm_new = nei_diversity(np.bincount(perm_new))
    
    perm_diffs.append(h_perm_old - h_perm_new)

# Calculate p-value
perm_diffs = np.array(perm_diffs)
p_value = np.mean(np.abs(perm_diffs) >= np.abs(observed_diff))

# Results
print(f"Observed difference in haplotype diversity: {observed_diff:.4f}")
print(f"P-value from permutation test: {p_value:.4f}")

if p_value < 0.05:
    print("Significant difference: the new method differs significantly from the old method.")
else:
    print("No significant difference: the new method performs similarly to the old method.")
