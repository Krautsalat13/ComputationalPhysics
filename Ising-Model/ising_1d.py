#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""1D Ising-model Monte-Carlo (Metropolis), Numba-accelerated.

Simulates an open chain of spins with the Metropolis algorithm in natural units
(k_B = 1, J = 1) and measures the internal energy U and specific heat C per spin
as functions of temperature, for chain lengths N = 10, 100, 1000 and N_s = 1e3 /
1e4 measurement sweeps. In one dimension there is no phase transition; the
simulation is compared against the exact results

    U/N = -(N-1)/N * tanh(1/T),
    C/N =  (N-1)/N * (1 / (T * cosh(1/T)))**2.

Each temperature starts from the ordered state, is equilibrated, and the
observables are then averaged over the following sweeps (one sweep = N attempted
flips). Figures are written to ``Plots/``.

Run from the project root:  ``python ising_1d.py``
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")   # non-interactive backend
import matplotlib.pyplot as plt
from numba import njit

# plain mathtext labels, no LaTeX needed
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = "18"


def Color(color_name):
    """RWTH Aachen corporate colour palette, returned as a 0-1 RGB tuple."""
    colors = {
        "RWTH": (0, 84, 159),
        "Schwarz": (0, 0, 0),
        "Petrol": (0, 97, 101),
        "Türkis": (0, 152, 161),
        "Grün": (87, 171, 39),
        "Maigrün": (189, 205, 0),
        "Gelb": (255, 237, 0),
        "Orange": (246, 168, 0),
        "Rot": (204, 7, 30),
        "Magenta": (227, 0, 102),
        "Bordeaux": (161, 16, 53),
        "Violett": (97, 33, 88),
        "Lila": (122, 111, 172),
        "RWTHlight": (142, 186, 229),
        "Bordeauxlight": (205, 139, 135),
    }
    if color_name not in colors:
        raise ValueError(f"Color '{color_name}' not found.")
    r, g, b = colors[color_name]
    return (r / 255, g / 255, b / 255)


SIM_COLORS = [Color(c) for c in ("RWTH", "Bordeaux")]

T_GRID = np.linspace(0.2, 4.0, 40)
N_EQUIL = 2000          # equilibration sweeps before measuring


@njit(cache=True)
def _run(N, T, n_measure, seed):
    """One Metropolis run for a chain of N spins. Returns U/N, C/N."""
    np.random.seed(seed)
    S = np.ones(N, dtype=np.int64)          # ordered start

    # Boltzmann factors for the only positive energy changes (dE = 2 at the
    # chain ends, dE = 4 in the interior).
    b2 = np.exp(-2.0 / T)
    b4 = np.exp(-4.0 / T)

    E = 0
    for i in range(N - 1):
        E -= S[i] * S[i + 1]

    def sweep(S, E):
        for _ in range(N):
            i = np.random.randint(N)
            nb = 0
            if i > 0:
                nb += S[i - 1]
            if i < N - 1:
                nb += S[i + 1]
            dE = 2 * S[i] * nb
            accept = False
            if dE <= 0:
                accept = True
            elif dE == 2:
                accept = np.random.random() < b2
            else:  # dE == 4
                accept = np.random.random() < b4
            if accept:
                S[i] = -S[i]
                E += dE
        return E

    for _ in range(N_EQUIL):
        E = sweep(S, E)

    sumE = 0.0
    sumE2 = 0.0
    for _ in range(n_measure):
        E = sweep(S, E)
        sumE += E
        sumE2 += E * E
    meanE = sumE / n_measure
    meanE2 = sumE2 / n_measure
    return meanE / N, (meanE2 - meanE * meanE) / (T * T) / N


def run_curve(N, n_measure, seed=1069):
    U = np.empty(len(T_GRID))
    C = np.empty(len(T_GRID))
    for k, T in enumerate(T_GRID):
        U[k], C[k] = _run(N, T, n_measure, seed + k)
    return U, C


def U_theory(N, T):
    return -(N - 1) / N * np.tanh(1 / T)


def C_theory(N, T):
    return (N - 1) / N * (1 / (T * np.cosh(1 / T))) ** 2


SAMPLES = [1000, 10000]

for N in [10, 100, 1000]:
    curves = {Ns: run_curve(N, Ns) for Ns in SAMPLES}

    # internal energy
    plt.figure(figsize=(10, 6))
    plt.plot(T_GRID, U_theory(N, T_GRID), "k--", lw=2, label="theory")
    for c, Ns in zip(SIM_COLORS, SAMPLES):
        plt.plot(T_GRID, curves[Ns][0], color=c, marker="o", ms=4, lw=1.3,
                 label=f"$N_s = {Ns}$")
    plt.xlabel("temperature $T$")
    plt.ylabel("$U/N$")
    plt.title(f"1D internal energy per spin, N = {N}")
    plt.grid(ls="--", alpha=0.5)
    plt.legend()
    plt.savefig(f"Plots/1D_U{N}.png", dpi=150, bbox_inches="tight")
    plt.close()

    # specific heat
    plt.figure(figsize=(10, 6))
    plt.plot(T_GRID, C_theory(N, T_GRID), "k--", lw=2, label="theory")
    for c, Ns in zip(SIM_COLORS, SAMPLES):
        plt.plot(T_GRID, curves[Ns][1], color=c, marker="o", ms=4, lw=1.3,
                 label=f"$N_s = {Ns}$")
    plt.xlabel("temperature $T$")
    plt.ylabel("$C/N$")
    plt.title(f"1D specific heat per spin, N = {N}")
    plt.grid(ls="--", alpha=0.5)
    plt.legend()
    plt.savefig(f"Plots/1D_C{N}.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"N={N} done")
