"""
Intra-chain charge mobility from diffusivity.

Reads precomputed eigenstates, eigenvectors, hopping rates, and coordinates
(from the static_disordered_hamiltonian and hopping_rate modules), solves the
master equation for the time-independent (static disorder) case, computes the
mean squared displacement, and extracts diffusivity and mobility.

Designed to be called from the central executing program as:

    mobility[ichain-1], mobility_index = diffusivity.main(
        ichain, n_site, sru_length, carrier, T, output_folder)

Master equation:
    dP_i(t)/dt = Sum_j [ P_j(t) k(j->i) ] - Sum_j [ P_i(t) k(i->j) ]

Diffusivity:
    D = <x^2(t)> / (2t)   in the diffusive regime

Mobility (Einstein relation):
    mu = q D / (kB T)


"""

import numpy as np
import scipy.sparse as sp


# ============================================================
# PHYSICAL CONSTANTS
# ============================================================

kB = 8.617333262e-5      # eV/K   (Boltzmann constant)
q = 1.0                  # elementary charge (eV per volt units)

# Numerical parameters
N_STEPS = 10000          # number of time steps for the master equation
DT_FRACTION = 0.4        # timestep as fraction of inverse fastest total rate
FIT_FRACTION = 3 / 4     # fit diffusivity from this fraction of the run onward

# Sparse option: 0.0 = exact dense (default). If > 0, rates below
# RATE_THRESHOLD * max_rate are set to zero and a sparse matvec is used.
# This is an APPROXIMATION - validate against a dense run before trusting it.
RATE_THRESHOLD = 0.0


# ============================================================
# READ AND PREPARE DATA
# ============================================================

def read_and_sort(ichain, n_site, carrier, folder_path):
    """
    Read eigenvalues, eigenvectors, hopping rates, and coordinates, then sort
    all of them by ascending coordinate. For holes, first un-reverse the rate
    matrix (Fortran stored hole rates with reversed indices).

    Returns
    -------
    eigenvalues  : (n_site,)          sorted by coordinate
    eigenvectors : (n_site, n_site)   rows = eigenstates, sorted by coordinate
    rates        : (n_site, n_site)   rates[i,j] = rate FROM state i TO state j
    co_ords      : (n_site,)          ascending
    """
    eigenvalues = np.loadtxt(f"{folder_path}localized_state_energy_{ichain}.out")
    eigenvectors = np.loadtxt(f"{folder_path}localized_state_coefficient_{ichain}.out")
    rates = np.loadtxt(f"{folder_path}localized_state_hopping_rate_{ichain}.out").reshape(n_site, n_site)
    ll_data = np.loadtxt(f"{folder_path}localized_state_ll_{ichain}.out")
    co_ords = ll_data[:, 0]

    # Un-reverse hole rate matrix to match eigenstate ordering
    if carrier == 'h':
        rates = rates[::-1, ::-1]

    # Sort everything by ascending coordinate
    sort_idx = np.argsort(co_ords)
    co_ords = co_ords[sort_idx]
    eigenvalues = eigenvalues[sort_idx]
    eigenvectors = eigenvectors[sort_idx]
    rates = rates[np.ix_(sort_idx, sort_idx)]

    return eigenvalues, eigenvectors, rates, co_ords


def find_initial_state(eigenvalues, n_site, carrier):
    """
    Find nc, the initial eigenstate near the chain midpoint.
    Electron: lowest energy in the window (bottom of spectrum).
    Hole:     highest energy in the window (top of spectrum).
    Window is n_site//2 +/- n_site//20 in the coordinate-sorted array.
    """
    mid = n_site // 2
    half_range = n_site // 20
    lo = max(0, mid - half_range)
    hi = min(n_site, mid + half_range)

    if carrier == 'e':
        nc = lo + np.argmin(eigenvalues[lo:hi])
    elif carrier == 'h':
        nc = lo + np.argmax(eigenvalues[lo:hi])
    else:
        raise ValueError(f"Unknown carrier type: {carrier}. Use 'e' or 'h'.")

    return nc


# ============================================================
# MASTER EQUATION  (returns MSD time series directly)
# ============================================================

def run_master_equation(eigenvectors, rates, nc, n_site, sru_length,
                        n_steps=N_STEPS):
    """
    Solve the master equation for time-independent rates by explicit Euler
    integration, returning the mean squared displacement at each time step.

    The site occupancy is never stored. Instead we precompute the per-eigenstate
    weight  w(alpha) = sum_site x_site^2 |C(site,alpha)|^2, so that
    <x^2(t)> = w . P(t) is a single dot product each step.

    Returns
    -------
    msd : (n_steps,)   mean squared displacement (Angstrom^2) at each step
    dt  : float        time step (s)
    """
    # Master equation needs k[i,j] = rate FROM state j TO state i, so transpose.
    k = np.ascontiguousarray(rates.T)
    k_out = np.sum(k, axis=0)          # total outgoing rate from each state

    # Time step from fastest total outgoing rate
    max_out = np.max(k_out)
    if max_out <= 0:
        raise ValueError("All hopping rates are zero; cannot set time step.")
    dt = DT_FRACTION / max_out

    # Optional sparse rate matrix (approximation; see RATE_THRESHOLD note)
    use_sparse = RATE_THRESHOLD > 0.0
    if use_sparse:
        cutoff = RATE_THRESHOLD * np.max(k)
        k_sparse = k.copy()
        k_sparse[k_sparse < cutoff] = 0.0
        k_sparse = sp.csr_matrix(k_sparse)

    # Per-eigenstate x^2 weight.
    # eigenvectors[alpha, site]; x_site centered so chain middle is at 0.
    x = (np.arange(n_site) - n_site // 2) * sru_length     # Angstrom
    w = (eigenvectors ** 2) @ (x ** 2)                     # length n_site

    # Initial condition
    P = np.zeros(n_site)
    P[nc] = 1.0

    msd = np.zeros(n_steps)
    msd[0] = w @ P

    for t in range(1, n_steps):
        if use_sparse:
            P = P + dt * (k_sparse @ P - k_out * P)
        else:
            P = P + dt * (k @ P - k_out * P)

        # Clamp negatives from numerical error and renormalize
        if np.any(P < 0):
            P = np.maximum(P, 0.0)
        P /= np.sum(P)

        msd[t] = w @ P

    return msd, dt


# ============================================================
# DIFFUSIVITY -> MOBILITY  (fit the MSD time series)
# ============================================================

def calculate_msd(msd, dt, T, n_steps=N_STEPS):
    """
    Extract diffusivity as the slope of <x^2> vs 2t in the diffusive regime,
    then mobility via the Einstein relation.

    Returns
    -------
    mobility : float   in cm^2/(V.s)
    D_cm2    : float   diffusivity in cm^2/s
    """
    leave = int(n_steps * FIT_FRACTION)
    times = np.arange(n_steps) * dt
    slope, _ = np.polyfit(times[leave:], msd[leave:], 1)

    D = slope * 0.5                               # Angstrom^2 / s
    D_cm2 = D * 1e-16                             # cm^2 / s

    # Einstein relation: mu = q D / (kB T)
    mu = D * (q / (kB * T))                       # Angstrom^2 / (V.s)
    mobility = mu * 1e-16                         # cm^2 / (V.s)

    return mobility, D_cm2


# ============================================================
# MAIN (called by polymer_intra_chain_mobility.py)
# ============================================================

def main(ichain, n_site, sru_length, carrier, T, output_folder):
    """
    Compute intra-chain mobility from diffusivity for one chain.

    Returns
    -------
    mobility : float   in cm^2/(V.s)
    index    : int     0 on success, 1 on error
    """
    try:
        eigenvalues, eigenvectors, rates, co_ords = read_and_sort(
            ichain, n_site, carrier, output_folder)

        nc = find_initial_state(eigenvalues, n_site, carrier)

        msd, dt = run_master_equation(
            eigenvectors, rates, nc, n_site, sru_length, n_steps=N_STEPS)

        mobility, D_cm2 = calculate_msd(msd, dt, T, n_steps=N_STEPS)

        print(f"Diffusivity {D_cm2:.6e} cm^2/s")
        print(f"Mobility {mobility:.10f} cm^2/(V.s)")

        return mobility, 0

    except Exception as e:
        print(f"Error in diffusivity.main: {e}", flush=True)
        return 0.0, 1