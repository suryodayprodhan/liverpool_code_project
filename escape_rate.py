
"""
Inter-chain escape rate k_esc, eq 16 of Troisi & Burke,
J. Chem. Phys. 162, 244706 (2025), evaluated as a double sum over eigenstates:

    k_AB = (2pi/hbar) V_inter^2 (1/Q_A) Sum_{i in A} exp(-beta E_i)
           x Sum_{j in B} |C_bj|^2 f(E_i - E_j, lambda_i + lambda_j)

    Q_A = Sum_{i in A} exp(-beta E_i),   lambda_j = lambda_0 x IPR_j,
    f(dE, L) = (4 pi L kB T)^(-1/2) exp[-(dE - L)^2 / (4 L kB T)]   (Marcus)

PDOS(E) = Sum_j |C_bj|^2 delta(E - E_j) is projected on a contact site b
of chain B, drawn at random for each pair. Normalisation follows from
Sum_j |C_bj|^2 = 1.
k_esc is the mean of k_AB over all ordered chain pairs A != B.
Eigenstates are read from the static_disordered_hamiltonian output;
energies are negated for holes.
"""

import os
import numpy as np

#################################################################
hbar = 6.582119569e-16   # eV.s   (reduced Planck constant)
kB = 8.617333262e-5      # eV/K   (Boltzmann constant)
#################################################################

##########################READ EIGENVALUES EIGENVECTORS####################################################################################
def read_eigenstates(n_chain, n_site, carrier, folder_path):
    energies = np.empty((n_chain, n_site))
    coefficients = np.empty((n_chain, n_site, n_site))

    for k, ichain in enumerate(range(1, n_chain + 1)):
        energies[k] = np.loadtxt(os.path.join(folder_path, f"localized_state_energy_{ichain}.out"))
        coefficients[k] = np.loadtxt(os.path.join(folder_path, f"localized_state_coefficient_{ichain}.out")).T

    if carrier == 'h':
        energies = -energies

    return energies, coefficients

########################################################################################################################

######################IPR###########################
def calculate_ipr(coefficients):
    c2 = coefficients**2
    return np.sum(c2**2, axis=-2) / np.sum(c2, axis=-2)**2
##############################################################

######################CALCULATE Q###########################
def calculate_partition_function(energies, T):
    beta = 1.0 / (kB * T)
    return np.sum(np.exp(-beta * energies), axis=-1)
##############################################################

##################ESCAPE RATE################################################################
def k_esc_pair(E_A, E_B, C_A, C_B, b, lambda_0, T, V_inter):
    kT = kB * T
    lam = lambda_0 * (calculate_ipr(C_A)[:, None] + calculate_ipr(C_B)[None, :])
    dE = E_A[:, None] - E_B[None, :]
    f = np.exp(-(dE - lam)**2 / (4.0 * lam * kT)) / np.sqrt(4.0 * np.pi * lam * kT)
    p = np.exp(-E_A / kT) / calculate_partition_function(E_A, T)
    w = C_B[b, :]**2
    return (2.0 * np.pi / hbar) * V_inter**2 * (p @ f @ w)
##############################################################################################

##############################EXECUTE###################################################
def main(n_chain, n_site, V_inter, T, carrier, output_folder, lambda_0,
         seed=0, verbose=True):
    energies, coefficients = read_eigenstates(n_chain, n_site, carrier, output_folder)
    rng = np.random.default_rng(seed)
    b = rng.integers(n_site, size=(n_chain, n_chain))

    k = np.full((n_chain, n_chain), np.nan)
    for A in range(n_chain):
        for B in range(n_chain):
            if A != B:
                k[A, B] = k_esc_pair(energies[A], energies[B], coefficients[A], coefficients[B],
                                     b[A, B], lambda_0, T, V_inter)

    k_esc = np.nanmean(k)

    if verbose:
        k_jack = np.array([np.nanmean(np.delete(np.delete(k, m, 0), m, 1)) for m in range(n_chain)])
        k_err = np.sqrt((n_chain - 1) / n_chain * np.sum((k_jack - k_jack.mean())**2))
        print(f"Escape rate (Marcus, lambda_0={lambda_0} eV, {n_chain*(n_chain-1)} pairs): "
              f"{k_esc:.4e} +/- {k_err:.4e} s^-1")
        print(f"  tau_esc = 1/k_esc = {1.0/k_esc:.4e} s "
              f"= {1.0/k_esc * 1e12:.2f} ps")

    return k_esc

#########################################################################################

