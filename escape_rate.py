"""
Escape rate of a charge carrier from a semiconducting polymer chain.

Implements the escape rate expressions from Troisi & Burke,
J. Chem. Phys. 162, 244706 (2025) (doi: 10.1063/5.0273118).

Two methods are provided:

1. Numerical PDOS (paper eq 17, lambda->0 limit):
   k_esc = (2pi/hbar) <|V_ab|^2> <PDOS>_T

   where <PDOS>_T is the Boltzmann-weighted thermal average of PDOS(E),
   computed numerically from the eigenvalue spectra of all chains.

2. Gaussian analytical (paper eq 19):
   k_esc = (2pi/hbar) |V|^2 (1/(2 sigma sqrt(pi))) exp(-beta^2 sigma^2 / 4)

   where sigma is the width of the PDOS, fitted from the eigenvalue spectra.

The eigenvalue spectra are read from the localized_state_energy files
already produced by the static_disordered_hamiltonian Fortran module.

Usage:
    k_esc, index = escape_rate.main(n_chain, n_site, V_inter, T, carrier,
                                     output_folder)
"""

import numpy as np
from scipy.stats import gaussian_kde


# ============================================================
# PHYSICAL CONSTANTS
# ============================================================

hbar = 6.582119569e-16   # eV.s   (reduced Planck constant)
kB = 8.617333262e-5      # eV/K   (Boltzmann constant)


# ============================================================
# READ EIGENVALUE SPECTRA
# ============================================================

def read_all_eigenvalues(n_chain, n_site, carrier, folder_path):
    """
    Read eigenvalues from all chains and return a single array.
    For holes, energies are negated (paper: "we have used the negative
    of the energy as the carriers are holes").

    Returns
    -------
    energies : (n_chain * n_site,)   all eigenvalues (eV)
    """
    all_E = []
    for ichain in range(1, n_chain + 1):
        filename = f"{folder_path}localized_state_energy_{ichain}.out"
        E = np.loadtxt(filename)
        assert len(E) == n_site, (
            f"Chain {ichain}: expected {n_site} eigenvalues, got {len(E)}")
        all_E.append(E)

    energies = np.concatenate(all_E)

    # For holes, negate energies so Boltzmann weighting exp(-beta E)
    # correctly favours the top of the spectrum
    if carrier == 'h':
        energies = -energies

    return energies


# ============================================================
# PDOS CONSTRUCTION
# ============================================================

def build_pdos(energies, n_site, bandwidth=None):
    """
    Build the partial density of states (PDOS) as a KDE-smoothed
    function from the eigenvalue spectra.

    PDOS(E) = DOS(E) / N_monomers   (paper eq 11)

    Since the KDE is built from all eigenvalues across all chains,
    it estimates DOS(E) normalized to integrate to 1. To get PDOS
    we divide by n_site (number of monomers per chain).

    Parameters
    ----------
    energies  : array    all eigenvalues (already sign-corrected for carrier)
    n_site    : int      monomers per chain
    bandwidth : float    KDE bandwidth in eV (None = automatic Scott's rule)

    Returns
    -------
    pdos_func : callable   PDOS(E) in units of eV^-1
    """
    if bandwidth is not None:
        kde = gaussian_kde(energies, bw_method=bandwidth / energies.std())
    else:
        kde = gaussian_kde(energies)

    # kde(E) integrates to 1 over all eigenvalues, so it's the normalized
    # DOS. PDOS = DOS / N = kde / n_site  (per monomer).
    def pdos_func(E):
        return kde(E) / n_site

    return pdos_func


# ============================================================
# ESCAPE RATE: NUMERICAL PDOS  (paper eq 17, lambda->0)
# ============================================================

def kesc_numerical(energies, n_site, V_sq, T):
    """
    Escape rate from eq 17 (vanishing reorganization energy limit):

    k_esc = (2pi/hbar) <|V|^2> <PDOS>_T

    where <PDOS>_T = integral[ exp(-beta E) PDOS(E)^2 dE ]
                   / integral[ exp(-beta E) PDOS(E) dE ]

    Parameters
    ----------
    energies : array   all eigenvalues (sign-corrected)
    n_site   : int     monomers per chain
    V_sq     : float   <|V_ab|^2> in eV^2
    T        : float   temperature in K

    Returns
    -------
    k_esc : float   escape rate in s^-1
    pdos_T : float  thermal average <PDOS>_T in eV^-1
    """
    beta = 1.0 / (kB * T)
    pdos_func = build_pdos(energies, n_site)

    # Integration grid: cover the energy range with margin
    E_min = energies.min() - 0.5
    E_max = energies.max() + 0.5
    E_grid = np.linspace(E_min, E_max, 5000)
    dE = E_grid[1] - E_grid[0]

    pdos_vals = pdos_func(E_grid)

    # Shift energies for numerical stability of exp(-beta E)
    E_shift = E_grid.min()
    boltz = np.exp(-beta * (E_grid - E_shift))

    numerator = np.sum(boltz * pdos_vals**2) * dE
    denominator = np.sum(boltz * pdos_vals) * dE

    pdos_T = numerator / denominator
    prefactor = 2.0 * np.pi / hbar
    k_esc = prefactor * V_sq * pdos_T

    return k_esc, pdos_T


# ============================================================
# ESCAPE RATE: GAUSSIAN ANALYTICAL  (paper eq 19)
# ============================================================

def kesc_gaussian(sigma, V_sq, T):
    """
    Escape rate from eq 19 (Gaussian PDOS, lambda->0):

    k_esc = (2pi/hbar) |V|^2 (1/(2 sigma sqrt(pi))) exp(-beta^2 sigma^2 / 4)

    Parameters
    ----------
    sigma : float   PDOS width (eV)
    V_sq  : float   <|V_ab|^2> in eV^2
    T     : float   temperature in K

    Returns
    -------
    k_esc : float   escape rate in s^-1
    """
    beta = 1.0 / (kB * T)
    prefactor = 2.0 * np.pi / hbar
    k_esc = prefactor * V_sq * (1.0 / (2.0 * sigma * np.sqrt(np.pi))) \
            * np.exp(-beta**2 * sigma**2 / 4.0)
    return k_esc


# ============================================================
# MAIN
# ============================================================

def main(n_chain, n_site, V_inter, T, carrier, output_folder,
         method='numerical', verbose=True):
    """
    Compute the escape rate from a polymer chain.

    Parameters
    ----------
    n_chain       : int     number of chain realizations
    n_site        : int     monomers per chain
    V_inter       : float   sqrt(<|V_ab|^2>), inter-chain coupling (eV)
    T             : float   temperature (K)
    carrier       : str     'e' or 'h'
    output_folder : str     path to Fortran output files (with trailing /)
    method        : str     'numerical' (eq 17) or 'gaussian' (eq 19)
    verbose       : bool    print results

    Returns
    -------
    k_esc : float   escape rate (s^-1)
    index : int     0 on success, 1 on error
    """
    try:
        V_sq = V_inter**2         # <|V_ab|^2> in eV^2

        # Read all eigenvalue spectra
        energies = read_all_eigenvalues(n_chain, n_site, carrier, output_folder)

        if method == 'numerical':
            k_esc, pdos_T = kesc_numerical(energies, n_site, V_sq, T)
            if verbose:
                print(f"Escape rate (numerical PDOS): {k_esc:.4e} s^-1")
                print(f"  <PDOS>_T = {pdos_T:.6f} eV^-1")
                print(f"  tau_esc = 1/k_esc = {1.0/k_esc:.4e} s "
                      f"= {1.0/k_esc * 1e12:.2f} ps")

        elif method == 'gaussian':
            sigma = np.std(energies)
            k_esc = kesc_gaussian(sigma, V_sq, T)
            if verbose:
                print(f"Escape rate (Gaussian, sigma={sigma:.4f} eV): "
                      f"{k_esc:.4e} s^-1")
                print(f"  tau_esc = 1/k_esc = {1.0/k_esc:.4e} s "
                      f"= {1.0/k_esc * 1e12:.2f} ps")

        else:
            raise ValueError(f"Unknown method '{method}'. "
                             f"Use 'numerical' or 'gaussian'.")

        # Also compute the other method for comparison
        if verbose:
            sigma = np.std(energies)
            k_gauss = kesc_gaussian(sigma, V_sq, T)
            k_num, _ = kesc_numerical(energies, n_site, V_sq, T)
            print(f"  (comparison: numerical {k_num:.4e}, "
                  f"Gaussian {k_gauss:.4e} s^-1)")

        return k_esc, 0

    except Exception as e:
        print(f"Error in escape_rate.main: {e}", flush=True)
        return 0.0, 1
