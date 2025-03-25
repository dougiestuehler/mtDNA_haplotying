import numpy as np
from scipy.stats import norm

def test_stats(method_h, method_pi, method_S, method_Hd, 
               ska_h, ska_pi, ska_S, ska_Hd, n, L):
    results = {}
    
    # h: Number of Haplotypes
    if method_h is not None and ska_h is not None:
        p1_h = method_h / n
        p2_h = ska_h / n
        p_pooled_h = (method_h + ska_h) / (2 * n)
        se_h = np.sqrt(p_pooled_h * (1 - p_pooled_h) * (2 / n))
        diff_h = p1_h - p2_h
        z_h = diff_h / se_h
        p_h = 2 * (1 - norm.cdf(abs(z_h)))
        results['h'] = (diff_h, se_h, p_h)
    else:
        results['h'] = (None, None, None)
    
    # pi: Nucleotide Diversity
    if method_pi is not None and ska_pi is not None and L is not None:
        k = n * (n - 1) / 2
        var1_pi = method_pi / (k * L)
        var2_pi = ska_pi / (k * L)
        se_pi = np.sqrt(var1_pi + var2_pi)
        diff_pi = method_pi - ska_pi
        z_pi = diff_pi / se_pi
        p_pi = 2 * (1 - norm.cdf(abs(z_pi)))
        results['pi'] = (diff_pi, se_pi, p_pi)
    else:
        results['pi'] = (None, None, None)
    
    # S: Segregating Sites
    if method_S is not None and ska_S is not None and L is not None:
        p1_S = method_S / L
        p2_S = ska_S / L
        p_pooled_S = (method_S + ska_S) / (2 * L)
        se_S = np.sqrt(p_pooled_S * (1 - p_pooled_S) * (2 / L))
        diff_S = p1_S - p2_S
        z_S = diff_S / se_S
        p_S = 2 * (1 - norm.cdf(abs(z_S)))
        results['S'] = (diff_S, se_S, p_S)
    else:
        results['S'] = (None, None, None)
    
    # Hd: Haplotype Diversity
    if method_Hd is not None and ska_Hd is not None:
        var1_Hd = 0 if method_Hd == 1 else method_Hd * (1 - method_Hd) / n
        var2_Hd = ska_Hd * (1 - ska_Hd) / n
        se_Hd = np.sqrt(var1_Hd + var2_Hd)
        diff_Hd = method_Hd - ska_Hd
        z_Hd = diff_Hd / se_Hd
        p_Hd = 2 * (1 - norm.cdf(abs(z_Hd)))
        results['Hd'] = (diff_Hd, se_Hd, p_Hd)
    else:
        results['Hd'] = (None, None, None)
    
    return results

# Datasets
ovis = {'WGA': {'h': 17, 'pi': 0.016, 'S': 1387, 'Hd': 1.0}, 'MLST': {'h': None, 'pi': None, 'S': None, 'Hd': None}, 'SKA': {'h': 15, 'pi': 0.0113, 'S': 753, 'Hd': 0.985}, 'n': 17, 'L': 16463}
f_intonsa = {'WGA': {'h': None, 'pi': None, 'S': None, 'Hd': None}, 'MLST': {'h': 124, 'pi': 0.00456, 'S': 474, 'Hd': 0.997}, 'SKA': {'h': 123, 'pi': 0.00347, 'S': 434, 'Hd': 0.996}, 'n': 149, 'L': 15220}
d_citri = {'WGA': {'h': None, 'pi': None, 'S': None, 'Hd': None}, 'MLST': {'h': 17, 'pi': 0.00099, 'S': None, 'Hd': 0.938}, 'SKA': {'h': 21, 'pi': 0.00098, 'S': 74, 'Hd': 0.974}, 'n': 31, 'L': 15038}

# Run tests
ovis_wga = test_stats(ovis['WGA']['h'], ovis['WGA']['pi'], ovis['WGA']['S'], ovis['WGA']['Hd'], ovis['SKA']['h'], ovis['SKA']['pi'], ovis['SKA']['S'], ovis['SKA']['Hd'], ovis['n'], ovis['L'])
ovis_mlst = test_stats(ovis['MLST']['h'], ovis['MLST']['pi'], ovis['MLST']['S'], ovis['MLST']['Hd'], ovis['SKA']['h'], ovis['SKA']['pi'], ovis['SKA']['S'], ovis['SKA']['Hd'], ovis['n'], ovis['L'])
f_intonsa_wga = test_stats(f_intonsa['WGA']['h'], f_intonsa['WGA']['pi'], f_intonsa['WGA']['S'], f_intonsa['WGA']['Hd'], f_intonsa['SKA']['h'], f_intonsa['SKA']['pi'], f_intonsa['SKA']['S'], f_intonsa['SKA']['Hd'], f_intonsa['n'], f_intonsa['L'])
f_intonsa_mlst = test_stats(f_intonsa['MLST']['h'], f_intonsa['MLST']['pi'], f_intonsa['MLST']['S'], f_intonsa['MLST']['Hd'], f_intonsa['SKA']['h'], f_intonsa['SKA']['pi'], f_intonsa['SKA']['S'], f_intonsa['SKA']['Hd'], f_intonsa['n'], f_intonsa['L'])
d_citri_wga = test_stats(d_citri['WGA']['h'], d_citri['WGA']['pi'], d_citri['WGA']['S'], d_citri['WGA']['Hd'], d_citri['SKA']['h'], d_citri['SKA']['pi'], d_citri['SKA']['S'], d_citri['SKA']['Hd'], d_citri['n'], d_citri['L'])
d_citri_mlst = test_stats(d_citri['MLST']['h'], d_citri['MLST']['pi'], d_citri['MLST']['S'], d_citri['MLST']['Hd'], d_citri['SKA']['h'], d_citri['SKA']['pi'], d_citri['SKA']['S'], d_citri['SKA']['Hd'], d_citri['n'], d_citri['L'])

# Print results
print("Ovis spp. Results:")
print("WGA vs. SKA:")
print(f"h: Difference = {ovis_wga['h'][0]:.4f} ± {ovis_wga['h'][1]:.4f}, p = {ovis_wga['h'][2]:.4f}")
print(f"π: Difference = {ovis_wga['pi'][0]:.4f} ± {ovis_wga['pi'][1]:.4f}, p = {ovis_wga['pi'][2]:.4e}")
print(f"S: Difference = {ovis_wga['S'][0]:.4f} ± {ovis_wga['S'][1]:.4f}, p = {ovis_wga['S'][2]:.4e}")
print(f"Hd: Difference = {ovis_wga['Hd'][0]:.4f} ± {ovis_wga['Hd'][1]:.4f}, p = {ovis_wga['Hd'][2]:.4f}")
print("MLST vs. SKA:")
for stat in ['h', 'pi', 'S', 'Hd']:
    diff, se, p = ovis_mlst[stat]
    print(f"{stat}: Difference = {diff if diff is None else f'{diff:.4f}'} ± {se if se is None else f'{se:.4f}'}, p = {p if p is None else f'{p:.4f}'}")

print("\nF. intonsa Results:")
print("WGA vs. SKA:")
for stat in ['h', 'pi', 'S', 'Hd']:
    diff, se, p = f_intonsa_wga[stat]
    print(f"{stat}: Difference = {diff if diff is None else f'{diff:.4f}'} ± {se if se is None else f'{se:.4f}'}, p = {p if p is None else f'{p:.4f}'}")
print("MLST vs. SKA:")
print(f"h: Difference = {f_intonsa_mlst['h'][0]:.4f} ± {f_intonsa_mlst['h'][1]:.4f}, p = {f_intonsa_mlst['h'][2]:.4f}")
print(f"π: Difference = {f_intonsa_mlst['pi'][0]:.4f} ± {f_intonsa_mlst['pi'][1]:.4f}, p = {f_intonsa_mlst['pi'][2]:.4e}")
print(f"S: Difference = {f_intonsa_mlst['S'][0]:.4f} ± {f_intonsa_mlst['S'][1]:.4f}, p = {f_intonsa_mlst['S'][2]:.4f}")
print(f"Hd: Difference = {f_intonsa_mlst['Hd'][0]:.4f} ± {f_intonsa_mlst['Hd'][1]:.4f}, p = {f_intonsa_mlst['Hd'][2]:.4f}")

print("\nD. citri Genbank Results:")
print("WGA vs. SKA:")
for stat in ['h', 'pi', 'S', 'Hd']:
    diff, se, p = d_citri_wga[stat]
    print(f"{stat}: Difference = {diff if diff is None else f'{diff:.4f}'} ± {se if se is None else f'{se:.4f}'}, p = {p if p is None else f'{p:.4f}'}")
print("MLST vs. SKA:")
print(f"h: Difference = {d_citri_mlst['h'][0]:.4f} ± {d_citri_mlst['h'][1]:.4f}, p = {d_citri_mlst['h'][2]:.4f}")
print(f"π: Difference = {d_citri_mlst['pi'][0]:.4f} ± {d_citri_mlst['pi'][1]:.4f}, p = {d_citri_mlst['pi'][2]:.4f}")
print(f"S: Difference = {d_citri_mlst['S'][0] if d_citri_mlst['S'][0] is None else f'{d_citri_mlst['S'][0]:.4f}'} ± {d_citri_mlst['S'][1] if d_citri_mlst['S'][1] is None else f'{d_citri_mlst['S'][1]:.4f}'}, p = {d_citri_mlst['S'][2] if d_citri_mlst['S'][2] is None else f'{d_citri_mlst['S'][2]:.4f}'}")
print(f"Hd: Difference = {d_citri_mlst['Hd'][0]:.4f} ± {d_citri_mlst['Hd'][1]:.4f}, p = {d_citri_mlst['Hd'][2]:.4f}")