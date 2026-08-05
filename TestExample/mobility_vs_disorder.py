"""
Parameter sweep for intra-chain mobility.

For both methods (steady-state and diffusivity), sweeps:
    beta in {0.05, 0.10}
    sigma_alpha_static over a range starting from 0.07 eV, the value set as starting_sigma

For each (method, beta, sigma_alpha_static) it:
    1. Writes an input file (based on the parent DEFAULT_INPUT template).
    2. Runs polymer_intra_chain_mobility.py in the parent directory.
    3. Reads mobility_list.dat and averages the mobility over all chains.

Outputs two PNG plots inside Test_case:
    mobility_steady_state.png   - avg mobility vs sigma_alpha_static (two beta lines)
    mobility_diffusivity.png    - avg mobility vs sigma_alpha_static (two beta lines)

IMPORTANT: Run this after executing f2py compilation script and loading Intel API environment.
"""
starting_sigma= 0.07
upto_sigma= 0.10
number_of_points= 4
beta_list= [ 0.1]

import os
import sys
import subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')            # no display needed, just save files
import matplotlib.pyplot as plt

# ============================================================
# PATHS
# ============================================================

script_dir  = os.path.dirname(os.path.abspath(__file__))      # TestExample/
parent_dir  = os.path.dirname(script_dir)                     # contains main program
main_script = os.path.join(parent_dir, 'polymer_intra_chain_mobility.py')
template_in = os.path.join(parent_dir, 'DEFAULT_INPUT')
mobility_list = os.path.join(parent_dir, 'mobility_list.dat')

# ============================================================
# SWEEP PARAMETERS
# ============================================================

methods = ['mobility_steady_state', 'mobility_diffusivity']
beta_values = beta_list
sigma_values = np.array( np.round( np.linspace( starting_sigma, upto_sigma, number_of_points), 4 ) )
# ============================================================
# INPUT FILE HANDLING
# ============================================================
# Base input template written from scratch (not read from DEFAULT_INPUT).
# Lines 3 (property), 7 (beta), and 8 (sigma_alpha_static) are overwritten
# for each run; everything else stays fixed.

template = [
    "# Main Input",   # 0
    "",               # 1
    "4",              # 2  n_chain
    "mobility_diffusivity ipr ll",  # 3  property
    "600",            # 4  n_site
    "pbc",            # 5  obc_or_pbc
    "0.0",            # 6  alpha
    "0.10",           # 7  beta
    "0.07",           # 8  sigma_alpha_static
    "0.0",            # 9  sigma_beta_static
    "marcus",         # 10 rate_equation_type
    "0.01",           # 11 sigma_alpha_dynamic
    "0.01",           # 12 sigma_beta_dynamic
    "0.1",            # 13 lambda
    "300.0",          # 14 T
    "100.0",          # 15 field
    "10.0",           # 16 sru_length
    "h",              # 17 carrier
    ""
]


def run_case(property_str, beta, sigma_alpha, tag):
    """
    Write the input file with the given property, beta, and sigma_alpha_static,
    run the main program, and return the average mobility over all chains.
    n_chain and n_site depend on the method.
    """
    lines = list(template)

    # Method-dependent n_chain and n_site
    if property_str == 'mobility_steady_state':
        lines[2] = '10'      # n_chain
        lines[4] = '200'     # n_site
    elif property_str == 'mobility_diffusivity':
        lines[2] = '25'     # n_chain
        lines[4] = '1000'    # n_site

    lines[3] = property_str          # property
    lines[7] = str(beta)             # beta
    lines[8] = str(sigma_alpha)      # sigma_alpha_static

    input_path = os.path.join(script_dir, 'input.txt')  # single reused file
    with open(input_path, 'w') as f:
        f.write('\n'.join(lines) + '\n')

    # Run the main program in the parent directory (where the compiled
    # Fortran modules and Output_Folder live). Use the same interpreter.
    print(f"    Running: {property_str}, beta={beta}, sigma_alpha={sigma_alpha}")
    subprocess.run([sys.executable, main_script, input_path],
                   cwd=parent_dir, check=True)

    # Read mobility_list.dat and average over chains
    mob = np.loadtxt(mobility_list)
    avg = float(np.mean(mob))
    return avg


# ============================================================
# RUN SWEEP
# ============================================================

# results[method][beta] = list of avg mobility, one per sigma value
results = {method: {beta: [] for beta in beta_values} for method in methods}

for method in methods:
    print(f"\n=== Method: {method} ===")
    for beta in beta_values:
        for sigma in sigma_values:
            tag = f"{method}_b{beta}_s{sigma}"
            avg = run_case(method, beta, sigma, tag)
            results[method][beta].append(avg)
            print(f"      avg mobility = {avg:.6e} cm^2/(V.s)")


# ============================================================
# SAVE NUMERIC RESULTS
# ============================================================

results_file = os.path.join(script_dir, 'mobility_values.dat')
with open(results_file, 'w') as f:
    f.write("# method  beta  sigma_alpha_static  avg_mobility(cm^2/Vs)\n")
    for method in methods:
        for beta in beta_values:
            for sigma, avg in zip(sigma_values, results[method][beta]):
                f.write(f"{method}  {beta}  {sigma}  {avg:.10e}\n")
print(f"\nNumeric results saved to {results_file}")


# ============================================================
# PLOTS
# ============================================================

def make_plot(method, filename, title):
    plt.figure(figsize=(7, 5))
    marker_cycle = ['o-', 's-', '^-', 'D-', 'v-', 'p-', '*-', 'x-', '+-', 'h-']
    for i, beta in enumerate(beta_values):
        marker = marker_cycle[i % len(marker_cycle)]
        plt.plot(sigma_values, results[method][beta], marker,
                 label=f'beta = {beta}', markersize=7, linewidth=1.5)
    plt.xlabel(r'$\sigma_{\alpha(\mathrm{static})}$ /eV')
    plt.ylabel('Mobility /cm$^2$/V$\cdot$s')
    plt.title(title)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    out_path = os.path.join(script_dir, filename)
    plt.savefig(out_path, dpi=200)
    plt.close()
    print(f"Saved plot: {out_path}")


make_plot('mobility_steady_state', 'mobility_steady_state.png',
          'Steady-State Mobility')
make_plot('mobility_diffusivity', 'mobility_diffusivity.png',
          'Mobility by Diffusivity')

print("\nSweep complete.")
