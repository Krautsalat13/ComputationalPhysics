#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Compare Euler, Euler-Cromer and Velocity-Verlet on the harmonic oscillator.

The 1D harmonic oscillator (unit mass and spring constant) obeys
    d2x/dt2 = a(x) = -x,
with initial conditions x(0) = 0, v(0) = 1, so the exact solution is
x(t) = sin(t), v(t) = cos(t) and the total energy E = v^2/2 + x^2/2 = 1/2 is
conserved for all time.

Three explicit integrators are compared at a fixed time step:
  * Euler            - not symplectic; the energy grows without bound.
  * Euler-Cromer     - symplectic; the energy stays bounded (oscillates).
  * Velocity-Verlet  - symplectic and second order; the energy is conserved
                       to high accuracy.

The script writes two figures to ``Plots/``:
  * ``integrators_trajectory.png`` - x(t) for each method vs the exact sine.
  * ``integrators_energy.png``     - total energy E(t) for each method.

Run from the project folder:  ``python oscillator_integrators.py``
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


def acceleration(x):
    """Harmonic-oscillator acceleration a(x) = -x (unit mass and stiffness)."""
    return -x


def integrate(method, dt, t_max, x0=0.0, v0=1.0):
    """Integrate the oscillator with the chosen method.

    Parameters
    ----------
    method : {"euler", "euler-cromer", "velocity-verlet"}
    dt     : time step
    t_max  : final time
    x0, v0 : initial position and velocity

    Returns
    -------
    t, x, v : equal-length arrays of time, position and velocity.
    """
    n = int(round(t_max / dt)) + 1
    t = np.linspace(0.0, (n - 1) * dt, n)
    x = np.empty(n)
    v = np.empty(n)
    x[0], v[0] = x0, v0

    for i in range(n - 1):
        if method == "euler":
            # Position and velocity both use the OLD state -> not symplectic.
            x[i + 1] = x[i] + v[i] * dt
            v[i + 1] = v[i] + acceleration(x[i]) * dt
        elif method == "euler-cromer":
            # Update velocity first, then use the NEW velocity for position.
            v[i + 1] = v[i] + acceleration(x[i]) * dt
            x[i + 1] = x[i] + v[i + 1] * dt
        elif method == "velocity-verlet":
            # Half-kick / drift / half-kick; second order and symplectic.
            a_i = acceleration(x[i])
            x[i + 1] = x[i] + v[i] * dt + 0.5 * a_i * dt ** 2
            v[i + 1] = v[i] + 0.5 * (a_i + acceleration(x[i + 1])) * dt
        else:
            raise ValueError(f"unknown method: {method}")
    return t, x, v


def energy(x, v):
    """Total energy E = v^2/2 + x^2/2 of the unit harmonic oscillator."""
    return 0.5 * v ** 2 + 0.5 * x ** 2


METHODS = ["euler", "euler-cromer", "velocity-verlet"]
LABELS = {"euler": "Euler",
          "euler-cromer": "Euler-Cromer",
          "velocity-verlet": "Velocity-Verlet"}

# Colour + line style per method, so the symplectic methods (which track each
# other very closely) can still be told apart where they overlap.
STYLE = {"euler":           (Color("RWTH"),     "-"),
         "euler-cromer":    (Color("Bordeaux"), "--"),
         "velocity-verlet": (Color("Petrol"),   "-.")}

dt = 0.1
t_max = 40.0

results = {m: integrate(m, dt, t_max) for m in METHODS}

# --- Figure 1: trajectory vs the exact solution --------------------------------
plt.figure(figsize=(10, 6))
t_ref = np.linspace(0, t_max, 2000)
plt.plot(t_ref, np.sin(t_ref), color="black", ls=":", lw=1.5, label="exact  $\\sin(t)$")
for m in METHODS:
    t, x, _ = results[m]
    c, ls = STYLE[m]
    plt.plot(t, x, color=c, ls=ls, lw=1.8, label=LABELS[m])
plt.xlabel("time $t$")
plt.ylabel("position $x(t)$")
plt.title(f"Harmonic oscillator trajectory ($\\Delta t = {dt}$)")
plt.ylim(-2.6, 2.6)
plt.grid(ls="--", alpha=0.5)
plt.legend()
plt.savefig("Plots/integrators_trajectory.png", dpi=150, bbox_inches="tight")
plt.close()

# --- Figure 2: energy drift ----------------------------------------------------
# Euler runs away to ~30 while the symplectic methods sit on E = 1/2, so a single
# y-axis hides one or the other. Left panel: full range (Euler blowing up); right
# panel: zoomed to E = 1/2 so Euler-Cromer and Velocity-Verlet are distinguishable.
fig, (ax_full, ax_zoom) = plt.subplots(1, 2, figsize=(13, 5.5))
for ax in (ax_full, ax_zoom):
    ax.axhline(0.5, color="black", ls=":", lw=1.5, label="exact  $E = 1/2$")
    for m in METHODS:
        t, x, v = results[m]
        c, ls = STYLE[m]
        ax.plot(t, energy(x, v), color=c, ls=ls, lw=1.8, label=LABELS[m])
    ax.set_xlabel("time $t$")
    ax.grid(ls="--", alpha=0.5)
ax_full.set_ylabel("total energy $E(t)$")
ax_full.set_title("full range")
ax_zoom.set_title("zoom near $E = 1/2$")
ax_zoom.set_ylim(0.45, 0.62)
ax_full.legend()
fig.suptitle(f"Energy conservation by integrator ($\\Delta t = {dt}$)")
fig.tight_layout()
fig.savefig("Plots/integrators_energy.png", dpi=150, bbox_inches="tight")
plt.close(fig)

print("wrote Plots/integrators_trajectory.png and Plots/integrators_energy.png")
