#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""2D Ising-model Monte-Carlo engine (Metropolis algorithm), Numba-accelerated.

Simulates the square-lattice Ising model with the Metropolis algorithm in
natural units (k_B = 1, J = 1). For a range of temperatures it measures the
internal energy U, the specific heat C and the absolute magnetisation |M| per
spin, for lattice sizes N = 10, 50, 100, with free and periodic boundaries and
a number of measurement sweeps N_s = 1e3 / 1e4.

A few notes on how the measurement is done, since these matter for getting
curves that actually follow the exact Onsager result:

  * Each run starts from the fully ordered state (all spins +1). Below T_c this
    avoids the system freezing into two domains with |M| ~ 0.5 (a metastable
    trap that single-spin Metropolis falls into easily with free boundaries);
    above T_c the ordered start relaxes to disorder within the equilibration
    phase anyway.
  * Each temperature is first equilibrated for N_EQUIL sweeps (one sweep = N*N
    attempted flips), and only then are U, E^2 and |M| averaged over the
    following N_s sweeps. So |M| is a proper thermal average <|M|>, not a single
    end-of-run snapshot.
  * The five possible Boltzmann factors are precomputed per temperature, so the
    inner loop never calls exp().

Data is written to ``Data/`` as ``2D_{U,C,M}_{N}N_{Ns}NS_{free,per}.npy``.
``plot_2d.py`` in the project root reads those files and plots them against the exact
theory.

Run from the project root:  ``python ising_2d.py``  (takes a few minutes).
"""
import time
import numpy as np
from numba import njit


# Temperature grid. Kept the same in plot_2d.py so the data and theory line up.
# A bit denser than a plain coarse grid so the specific-heat peak near
# T_c ~ 2.269 is actually resolved.
T_GRID = np.linspace(0.2, 4.0, 40)

N_EQUIL = 2000          # equilibration sweeps before measuring


@njit(cache=True)
def _neighbour_sum(S, N, i, j, periodic):
    """Sum of the (up to four) nearest-neighbour spins of site (i, j)."""
    total = 0
    # up / down
    if i > 0:
        total += S[i - 1, j]
    elif periodic:
        total += S[N - 1, j]
    if i < N - 1:
        total += S[i + 1, j]
    elif periodic:
        total += S[0, j]
    # left / right
    if j > 0:
        total += S[i, j - 1]
    elif periodic:
        total += S[i, N - 1]
    if j < N - 1:
        total += S[i, j + 1]
    elif periodic:
        total += S[i, 0]
    return total


@njit(cache=True)
def _total_energy(S, N, periodic):
    """Total energy E = -sum_<ij> S_i S_j, counting each bond once."""
    E = 0
    for i in range(N):
        for j in range(N):
            # right and down bonds only, so every bond is counted once
            if j < N - 1:
                E -= S[i, j] * S[i, j + 1]
            elif periodic:
                E -= S[i, j] * S[i, 0]
            if i < N - 1:
                E -= S[i, j] * S[i + 1, j]
            elif periodic:
                E -= S[i, j] * S[0, j]
    return E


@njit(cache=True)
def _run(N, T, n_measure, periodic, seed):
    """One Metropolis run at temperature T. Returns U/N^2, C/N^2, <|M|>/N^2."""
    np.random.seed(seed)
    S = np.ones((N, N), dtype=np.int64)          # ordered start, all spins +1

    # Precompute the only two positive energy changes' acceptance probabilities
    # (dE can be -8, -4, 0, +4, +8; only +4 and +8 need a Boltzmann factor).
    boltz4 = np.exp(-4.0 / T)
    boltz8 = np.exp(-8.0 / T)

    n_sites = N * N

    def sweep(S, E, M):
        for _ in range(n_sites):
            i = np.random.randint(N)
            j = np.random.randint(N)
            nb = _neighbour_sum(S, N, i, j, periodic)
            dE = 2 * S[i, j] * nb
            accept = False
            if dE <= 0:
                accept = True
            elif dE == 4:
                accept = np.random.random() < boltz4
            else:  # dE == 8
                accept = np.random.random() < boltz8
            if accept:
                M -= 2 * S[i, j]
                S[i, j] = -S[i, j]
                E += dE
        return E, M

    # equilibration
    E = _total_energy(S, N, periodic)
    M = 0
    for i in range(N):
        for j in range(N):
            M += S[i, j]
    for _ in range(N_EQUIL):
        E, M = sweep(S, E, M)

    # measurement
    sumE = 0.0
    sumE2 = 0.0
    sumM = 0.0
    for _ in range(n_measure):
        E, M = sweep(S, E, M)
        sumE += E
        sumE2 += E * E
        sumM += abs(M)

    meanE = sumE / n_measure
    meanE2 = sumE2 / n_measure
    U = meanE / n_sites
    C = (meanE2 - meanE * meanE) / (T * T) / n_sites
    Mabs = sumM / n_measure / n_sites
    return U, C, Mabs


def run_curve(N, n_measure, periodic, seed=1556):
    """Run the whole temperature sweep for one (size, N_s, boundary) setting."""
    U = np.zeros(len(T_GRID))
    C = np.zeros(len(T_GRID))
    M = np.zeros(len(T_GRID))
    for k, T in enumerate(T_GRID):
        # a different seed per temperature keeps the runs independent
        U[k], C[k], M[k] = _run(N, T, n_measure, periodic, seed + k)
    return U, C, M


def gen_data():
    sizes = [10, 50, 100]
    samples = [1000, 10000]
    boundaries = [("per", True), ("free", False)]
    for N in sizes:
        for n_measure in samples:
            for name, periodic in boundaries:
                t0 = time.time()
                U, C, M = run_curve(N, n_measure, periodic)
                np.save(f"Data/2D_U_{N}N_{n_measure}NS_{name}.npy", U)
                np.save(f"Data/2D_C_{N}N_{n_measure}NS_{name}.npy", C)
                np.save(f"Data/2D_M_{N}N_{n_measure}NS_{name}.npy", M)
                print(f"N={N:3d}  Ns={n_measure:5d}  {name:4s}  "
                      f"({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    gen_data()
