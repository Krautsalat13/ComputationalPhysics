#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Plot stored quantum-tunnelling wavefunctions from the barrier simulation.

``SE_Potential_Barrier.py`` runs the (expensive) time-dependent Schrodinger
solver and stores the complex wavefunction Phi(x, t) under ``Data/`` at a set
of times, for two systems:
  * free    - free particle (no potential), and
  * barrier - square potential barrier at 50 <= x <= 50.5.

This script loads those snapshots and plots the probability density
P(x,t) = |Phi|^2 * Delta, so the figures can be regenerated instantly without
rerunning the solver. It also plots the transmitted probability versus time
(stored in ``Data/Psum_with.npy`` / ``Psum_without.npy``), which quantifies how
much of the packet gets past the barrier.

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


# colours for simulation curves, RWTH and Bordeaux first
SIM_COLORS = [Color(c) for c in ("RWTH", "Bordeaux", "Petrol", "Orange", "Grün")]

# Spatial grid and time step used by SE_Potential_Barrier.py
L = 1001
Delta = 0.1
tau = 0.001
x = np.linspace(0, 100, L)
barrier = (50.0, 50.5)          # location of the square barrier


def density(system, time):
    """Load Phi(x, time) for the given system and return P = |Phi|^2 * Delta."""
    phi = np.load(f"Data/TDSE_{system}_times{time}.npy")
    return np.abs(phi) ** 2 * Delta


def plot_snapshots(system, title, outfile, times=(0, 15, 30, 45), show_barrier=False):
    plt.figure(figsize=(10, 6))
    for c, t in zip(SIM_COLORS, times):
        plt.plot(x, density(system, t), lw=1.6, color=c, label=f"$t = {t}$")
    if show_barrier:
        plt.axvspan(*barrier, color="tab:green", alpha=0.4, label="barrier")
    plt.xlabel("position $x$")
    plt.ylabel(r"probability $|\Phi|^2\,\Delta$")
    plt.title(title)
    plt.xlim(0, 100)
    plt.grid(ls="--", alpha=0.5)
    plt.legend()
    plt.savefig(outfile, dpi=150, bbox_inches="tight")
    plt.close()


# Free particle: the packet travels to the right and spreads.
plot_snapshots("free", "Free wave packet (no barrier)",
               "Plots/free_particle.png")

# With the barrier: the packet partly reflects and partly tunnels through.
plot_snapshots("barrier", "Wave packet meeting a square barrier",
               "Plots/tunnelling.png", show_barrier=True)

# Transmitted probability (fraction of the packet beyond the barrier) vs time.
t = np.arange(len(np.load("Data/Psum_with.npy"))) * tau
plt.figure(figsize=(10, 6))
plt.plot(t, np.load("Data/Psum_without.npy"), lw=1.8, color=Color("RWTH"), label="no barrier")
plt.plot(t, np.load("Data/Psum_with.npy"), lw=1.8, color=Color("Bordeaux"), label="with barrier")
plt.xlabel("time $t$")
plt.ylabel("probability beyond $x > 50$")
plt.title("Transmitted probability vs time")
plt.grid(ls="--", alpha=0.5)
plt.legend()
plt.savefig("Plots/transmission.png", dpi=150, bbox_inches="tight")
plt.close()

print("wrote Plots/free_particle.png, Plots/tunnelling.png, Plots/transmission.png")
