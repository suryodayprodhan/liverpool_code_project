"""
3D charge diffusivity from the complete intra-chain + inter-chain model.

Orchestrates the full pipeline:
  1. diffusivity.py  -> contour diffusivity D_contour (per chain, averaged)
  2. escape_rate.py  -> inter-chain escape rate k_esc (from eigenvalue spectra)
  3. WLC integral    -> 3D diffusivity D_3D (from D_contour, k_esc, and P)

Implements the worm-like chain (WLC) model of Carbone & Troisi,
J. Phys. Chem. Lett. 2014, 5, 2637-2641, combined with the escape rate
expression of Troisi & Burke, J. Chem. Phys. 162, 244706 (2025).

The charge performs a 1D random walk along the chain contour with
diffusion coefficient D_contour. The chain is folded in 3D as a WLC with
persistence length P. With escape rate k_esc = 1/tau2 the charge hops to
a new (uncorrelated) chain. The 3D diffusivity is:

    D_3D = (1/(6 tau2)) Integral_0^inf <R^2(t)> (1/tau2) exp(-t/tau2) dt

where <R^2(t)> is the WLC mean square 3D displacement (Carbone eq 8).

Limits:
    Rigid rod  (P -> inf):  D_3D -> D_contour / 3
    Gaussian   (P -> 0):    D_3D -> (P/3) sqrt(D_contour / tau2)

Usage from the central executing program:

    import diffusivity_3d
    D_3D, mu_3D, index = diffusivity_3d.main(
        n_chain, n_site, sru_length, carrier, T,
        V_inter, P, output_folder)

NOTE: The WLC model assumes ZERO external field. The rates feeding
D_contour should be computed at zero (or negligibly small) field.
"""

import numpy as np
from scipy.integrate import quad
from scipy.special import erfcx

import diffusivity
import escape_rate


# ============================================================
# PHYSICAL CONSTANTS
# ============================================================

kB = 8.617333262e-5      # eV/K   (Boltzmann constant)
q = 1.0                  # elementary charge (eV per volt units)

# Numerical parameters for the WLC integral
QUAD_LIMIT = 200
REL_TOL = 1.0e-10


# ============================================================
# WLC MEAN SQUARE 3D DISPLACEMENT  (Carbone eq 8)
# ============================================================

def msd_3d_wlc(t, D_contour, P):
    """
    Mean square 3D displacement of a charge diffusing along a WLC contour.

    <R^2(t)> = 2 P^2 [ (2/sqrt(pi)) z - 1 + erfcx(z) ]
    z = sqrt(D_contour t) / P

    erfcx(z) = exp(z^2)(1 - erf(z)) avoids overflow at large z.
    """
    z = np.sqrt(D_contour * np.asarray(t, dtype=float)) / P
    return 2.0 * P**2 * ((2.0 / np.sqrt(np.pi)) * z - 1.0 + erfcx(z))


# ============================================================
# LIMITING CASES  (Carbone eqs 14, 15)
# ============================================================

def d3d_rigid_rod(D_contour):
    """Rigid-rod limit: D_3D = D_contour / 3. tau2 drops out."""
    return D_contour / 3.0


def d3d_gaussian(D_contour, P, k_esc):
    """Gaussian-chain limit: D_3D = (P/3) sqrt(D_contour k_esc)."""
    return (P / 3.0) * np.sqrt(D_contour * k_esc)


# ============================================================
# 3D DIFFUSIVITY  (Carbone eq 13 with eq 8 as integrand)
# ============================================================

def calculate_d3d(D_contour, P, k_esc):
    """
    3D diffusion coefficient from the residence-time-averaged WLC MSD.

    With substitution u = t/tau2:
        D_3D = (1/(6 tau2)) Integral_0^inf <R^2(u tau2)> exp(-u) du

    Parameters
    ----------
    D_contour : float   contour diffusivity (Angstrom^2/s)
    P         : float   persistence length (Angstrom)
    k_esc     : float   interchain escape rate (s^-1)

    Returns
    -------
    D_3D : float   3D diffusivity (Angstrom^2/s)
    """
    if D_contour <= 0.0:
        raise ValueError("D_contour must be positive.")
    if P <= 0.0:
        raise ValueError("Persistence length P must be positive.")
    if k_esc <= 0.0:
        raise ValueError("Escape rate k_esc must be positive.")

    tau2 = 1.0 / k_esc

    def integrand(u):
        return msd_3d_wlc(u * tau2, D_contour, P) * np.exp(-u)

    integral, _ = quad(integrand, 0.0, np.inf,
                       limit=QUAD_LIMIT, epsrel=REL_TOL)

    return integral / (6.0 * tau2)


# ============================================================
# COMPUTE D_CONTOUR FROM ALL CHAINS  (calls diffusivity.py)
# ============================================================

def compute_average_d_contour(n_chain, n_site, sru_length, carrier, T,
                               output_folder ):
    """
    Run the master-equation diffusivity calculation for each chain and
    return the average contour diffusivity in Angstrom^2/s.

    Calls diffusivity.py's sub-functions directly to extract D_contour
    before the cm^2 conversion.
    """
    verbose=True
    
    D_list = []
    
    for ichain in range(1, n_chain + 1):
        if verbose:
            print(f"  Chain {ichain}/{n_chain}: master equation...", flush=True)

        eigenvalues, eigenvectors, rates, co_ords = diffusivity.read_and_sort(
            ichain, n_site, carrier, output_folder)

        nc = diffusivity.find_initial_state(eigenvalues, n_site, carrier)

        msd_arr, dt = diffusivity.run_master_equation(
            eigenvectors, rates, nc, n_site, sru_length,
            n_steps=diffusivity.N_STEPS)

        mobility, D_cm2 = diffusivity.calculate_msd(
            msd_arr, dt, T, n_steps=diffusivity.N_STEPS)

        D_contour = D_cm2 * 1e16     # convert cm^2/s -> Angstrom^2/s
        D_list.append(D_contour)

        if verbose:
            print(f"    D_contour = {D_contour:.4e} Å²/s")

    D_avg = np.mean(D_list)
    D_std = np.std(D_list)

    if verbose:
        print(f"  Average D_contour = {D_avg:.4e} ± {D_std:.4e} Å²/s")

    return D_avg


# ============================================================
# MAIN
# ============================================================

def main(n_chain, n_site, sru_length, carrier, T, 
         output_folder, V_inter, P):
    """
    Compute the 3D diffusivity and mobility from the full pipeline.

    Parameters
    ----------
    n_chain       : int     number of chain realizations
    n_site        : int     monomers per chain
    sru_length    : float   site length (Angstrom)
    carrier       : str     'e' or 'h'
    T             : float   temperature (K)
    V_inter       : float   sqrt(<|V_ab|^2>), inter-chain coupling (eV)
    P             : float   persistence length (Angstrom)
    output_folder : str     path to Fortran output files (with trailing /)
    verbose       : bool    print progress and results

    Returns
    -------
    D_3D_cm2 : float   3D diffusivity in cm^2/s
    mu_3D    : float   3D mobility in cm^2/(V.s)
    index    : int     0 on success, 1 on error
    """
    verbose=True

    try:
        if verbose:
            print("=" * 60)
            print("3D DIFFUSIVITY CALCULATION")
            print("=" * 60)
            print(f"  n_chain={n_chain}, n_site={n_site}, "
                  f"sru_length={sru_length} Å")
            print(f"  carrier={carrier}, T={T} K")
            print(f"  V_inter={V_inter} eV, P={P} Å")

        # ---- Step 1: contour diffusivity from master equation ----
        if verbose:
            print(f"\nStep 1: Contour diffusivity (averaging {n_chain} chains)")

        D_contour = compute_average_d_contour(
            n_chain, n_site, sru_length, carrier, T, output_folder)

        # ---- Step 2: escape rate from eigenvalue spectra ----
        if verbose:
            print(f"\nStep 2: Escape rate from eigenvalue spectra")

        k_esc, esc_index = escape_rate.main(
            n_chain, n_site, V_inter, T, carrier, output_folder,
            method='numerical', verbose=verbose)

        if esc_index != 0:
            raise RuntimeError("escape_rate.main failed.")

        tau2 = 1.0 / k_esc

        # ---- Step 3: 3D diffusivity via WLC integral ----
        if verbose:
            print(f"\nStep 3: WLC 3D diffusivity integral")

        D_3D = calculate_d3d(D_contour, P, k_esc)          # Angstrom^2/s
        D_3D_cm2 = D_3D * 1e-16                            # cm^2/s

        # Mobility via Einstein relation
        mu_3D = D_3D * (q / (kB * T)) * 1e-16              # cm^2/(V.s)

        # Limiting values for context
        D_rr = d3d_rigid_rod(D_contour) * 1e-16
        D_gc = d3d_gaussian(D_contour, P, k_esc) * 1e-16

        # Regime indicator
        z_star = np.sqrt(D_contour * tau2) / P

        if verbose:
            print(f"\n{'=' * 60}")
            print(f"RESULTS")
            print(f"{'=' * 60}")
            print(f"  D_contour     = {D_contour * 1e-16:.6e} cm²/s")
            print(f"  k_esc         = {k_esc:.4e} s⁻¹ "
                  f"(tau2 = {tau2:.4e} s = {tau2*1e12:.2f} ps)")
            print(f"  D_3D          = {D_3D_cm2:.6e} cm²/s")
            print(f"  mu_3D         = {mu_3D:.6e} cm²/(V·s)")
            print(f"  (rigid-rod limit: {D_rr:.6e} cm²/s)")
            print(f"  (Gaussian limit:  {D_gc:.6e} cm²/s)")
            print(f"  z* = sqrt(D_contour tau2)/P = {z_star:.4f}")
            print(f"  (z* << 1: rod-like; z* >> 1: coil-like)")
            print(f"{'=' * 60}")

        return D_3D_cm2, mu_3D, 0

    except Exception as e:
        print(f"Error in diffusivity_3d.main: {e}", flush=True)
        return 0.0, 0.0, 1
