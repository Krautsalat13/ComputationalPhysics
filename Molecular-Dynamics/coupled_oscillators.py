#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Velocity-Verlet simulation of a chain of coupled harmonic oscillators.

A line of N equal masses is connected by identical springs (unit mass and
stiffness, free ends). Oscillator n feels the force from its neighbours,
    F_0     = -(x_0 - x_1),
    F_{N-1} = -(x_{N-1} - x_{N-2}),
    F_n     = -(2 x_n - x_{n-1} - x_{n+1})   for the interior masses,
and the whole chain is advanced with the symplectic Velocity-Verlet scheme.

Two initial conditions are shown:
  * a single displaced mass in the middle of the chain, and
  * a standing-wave (normal mode) x_n(0) = sin(pi*j*(n+1)/(N+1)).

The script writes to ``Plots/``:
  * ``coupled_positions.png`` - displacement of selected masses vs time, and
  * ``coupled_energy.png``    - total energy vs time (a check that the
                                symplectic integrator conserves energy).

Run from the project folder:  ``python coupled_oscillators.py``
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


def forces(x):
    """Spring forces on every mass of a free-ended chain (vectorised)."""
    f = np.empty_like(x)
    f[1:-1] = -(2 * x[1:-1] - x[:-2] - x[2:])
    f[0] = -(x[0] - x[1])
    f[-1] = -(x[-1] - x[-2])
    return f


def initial_condition(N, kind):
    """Return the initial displacements for the requested configuration.

    kind="middle" displaces the central mass; kind="mode<j>" sets the j-th
    standing wave, e.g. "mode1" or "mode8".
    """
    x = np.zeros(N)
    if kind == "middle":
        x[N // 2 - 1] = 1.0
    elif kind.startswith("mode"):
        j = int(kind[4:])
        n = np.arange(N) + 1
        x[:] = np.sin(np.pi * j * n / (N + 1))
    else:
        raise ValueError(f"unknown initial condition: {kind}")
    return x


def simulate(N, dt, steps, kind):
    """Velocity-Verlet evolution of the chain; returns t, X (steps+1, N), E."""
    x = initial_condition(N, kind)
    v = np.zeros(N)

    X = np.empty((steps + 1, N))
    E = np.empty(steps + 1)
    X[0] = x

    a = forces(x)
    for i in range(steps):
        x = x + v * dt + 0.5 * a * dt ** 2
        a_new = forces(x)
        v = v + 0.5 * (a + a_new) * dt
        a = a_new
        X[i + 1] = x
        # total energy: kinetic + spring potential energy of the bonds
        E[i + 1] = 0.5 * np.sum(v ** 2) + 0.5 * np.sum(np.diff(x) ** 2)
    E[0] = 0.5 * np.sum(forces(initial_condition(N, kind)) * 0)  # KE(0)=0
    E[0] = 0.5 * np.sum(np.diff(X[0]) ** 2)
    t = np.arange(steps + 1) * dt
    return t, X, E


N = 16
dt = 0.01
steps = 100000      # t_max = 1000

# A single standing wave (normal mode j=1) stays coherent; a localized kick
# spreads energy through the chain.
t, X, E = simulate(N, dt, steps, "mode1")

# --- Figure 1: positions of a few representative masses ------------------------
# Pick masses with different amplitudes in the j=1 mode (x_1 and x_N coincide by
# symmetry, so showing both would just draw the same curve twice).
shown = [(0, Color("RWTH"), "-"),
         (3, Color("Bordeaux"), "--"),
         (7, Color("Petrol"), "-.")]
plt.figure(figsize=(10, 6))
for n, c, ls in shown:
    plt.plot(t[: 2000], X[: 2000, n], color=c, ls=ls, lw=1.6,
             label=f"mass $x_{{{n + 1}}}$")
plt.xlabel("time $t$")
plt.ylabel("displacement")
plt.title(f"Coupled chain, normal mode $j=1$ ($N={N}$, $\\Delta t={dt}$)")
plt.grid(ls="--", alpha=0.5)
plt.legend()
plt.savefig("Plots/coupled_positions.png", dpi=150, bbox_inches="tight")
plt.close()

# --- Figure 2: energy conservation over a long run ----------------------------
plt.figure(figsize=(10, 6))
plt.plot(t, E, lw=1.3, color=Color("RWTH"))
plt.xlabel("time $t$")
plt.ylabel("total energy $E(t)$")
plt.title(f"Energy conservation of the chain ($N={N}$, $\\Delta t={dt}$)")
plt.ylim(0, max(E) * 1.5)
plt.grid(ls="--", alpha=0.5)
plt.savefig("Plots/coupled_energy.png", dpi=150, bbox_inches="tight")
plt.close()

print(f"energy drift over t=0..{t[-1]:.0f}: "
      f"{(E.max() - E.min()) / E.mean():.2e} (relative)")
print("wrote Plots/coupled_positions.png and Plots/coupled_energy.png")
