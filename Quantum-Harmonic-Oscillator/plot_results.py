#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Plot stored quantum-harmonic-oscillator probability snapshots.

``QHO.py`` runs the (expensive) time evolution and stores the probability
density |Phi(x,t)|^2 at times t = 0, 2, 4, 6, 8, 10 for several parameter
sets under ``Data/`` (file name pattern ``tami_t{t}_params{Omega}{sigma}{x0}``).

This script simply loads those snapshots and draws the probability density
P(x,t) = |Phi|^2 * Delta against position for a chosen parameter set, so the
figures can be regenerated in a second without rerunning the full solver.

Output: PNG figures to ``Plots/``.  Run from the project folder:
``python plot_results.py``
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")   # non-interactive backend
import matplotlib.pyplot as plt

# plain mathtext labels, no LaTeX needed
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = "16"


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


# one colour per snapshot time, RWTH and Bordeaux first
SIM_COLORS = [Color(c) for c in
              ("RWTH", "Bordeaux", "Petrol", "Orange", "Grün", "Violett")]

# Spatial grid used by QHO.py
L = 1201
Delta = 0.025
x = np.linspace(-15, 15, L)
times = [0, 2, 4, 6, 8, 10]


def plot_evolution(params, title, outfile, xlim=(-6, 6)):
    """Plot |Phi|^2 * Delta at each stored time for one parameter set."""
    plt.figure(figsize=(10, 6))
    for c, t in zip(SIM_COLORS, times):
        P = np.load(f"Data/tami_t{t}_params{params}.npy")
        plt.plot(x, P * Delta, lw=1.6, color=c, label=f"$t = {t}$")
    plt.xlabel("position $x$")
    plt.ylabel(r"probability $|\Phi|^2\,\Delta$")
    plt.title(title)
    plt.xlim(*xlim)
    plt.grid(ls="--", alpha=0.5)
    plt.legend(ncol=2)
    plt.savefig(outfile, dpi=150, bbox_inches="tight")
    plt.close()


# Displaced packet in the well (Omega=1, sigma=1, x0=1): it sloshes back and
# forth like a classical particle while keeping its shape.
plot_evolution("111",
               r"Harmonic oscillator, $\Omega=1,\ \sigma=1,\ x_0=1$",
               "Plots/qho_probability_111.png")

# Wider initial packet (Omega=1, sigma=2, x0=0): a near-stationary breathing
# state centred in the well.
plot_evolution("120",
               r"Harmonic oscillator, $\Omega=1,\ \sigma=2,\ x_0=0$",
               "Plots/qho_probability_120.png")

print("wrote Plots/qho_probability_111.png and Plots/qho_probability_120.png")
