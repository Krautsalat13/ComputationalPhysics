#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Time-dependent Schrodinger equation: quantum tunnelling through a barrier.

Propagates a Gaussian wave packet (units m = hbar = 1) governed by the 1D
time-dependent Schrodinger equation using the second-order product formula.
The spatial domain 0 <= x <= 100 is discretised on L = 1001 points; one time
step is split into kinetic sub-steps acting on neighbouring pairs of grid
points (K1 on even pairs, K2 on odd pairs) and a diagonal potential step (V),
applied symmetrically so the scheme is second-order accurate and unitary.

Two systems are selected with the ``task`` flag:
  * task = 0 - free particle (no potential), and
  * task = 1 - square potential barrier V = 2 on 50 <= x <= 50.5, higher than
               the packet's mean energy, so part of the packet tunnels through
               and part is reflected.

Snapshots of the wavefunction are written to ``Data/`` and probability-density
figures to ``Plots/``.  Run from the project folder:
``python SE_Potential_Barrier.py``  (expensive: m = 50000 time steps).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")   # non-interactive backend
import matplotlib.pyplot as plt
from numba import njit, vectorize, float64

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


pi      = np.pi
rand    = np.random

# --- physical / numerical parameters ---
sigma = 3           # width of the initial Gaussian wave packet
q = 1               # central wavenumber (sets the packet's momentum)
x0 = 20             # initial centre of the packet
task = 1            # 0 = free particle, 1 = square potential barrier
Delta = 0.1         # spatial grid spacing
L = 1001            # number of grid points (domain 0..100)
tau = 0.001         # time step
m = 50000           # number of time steps
# cos/sin entries of the 2x2 kinetic propagator for one pair of grid points
c = np.cos(tau/(4*Delta**2))
s = 1j*np.sin(tau/(4*Delta**2))

@njit
def Phi0(x):
    return (2*pi*sigma**2)**(-1/4) *np.exp(1j*q*(x-x0))*np.exp(-(x-x0)**2/(4*sigma**2))
    
x = np.linspace(0,100,L) 

@njit
def V(x):
    if task == 0:
        return 0
    if task == 1:
        if 50<= x<= 50.5:
            return 2
        else:
            return 0

Vvec= np.vectorize(V)    
Phi = np.array(Phi0(x)).astype(np.complex128)

M = np.array([[c,s],[s,c]], dtype=np.complex128) 
#Vvec = np.heaviside(x-50, 1)*2 - np.heaviside(x-50.5, 0)*2
Vx = np.exp(-1j*tau*(Delta**(-2)+Vvec(x)))
#Vx = np.exp(-1j*tau*(Delta**(-2)+Vvec))
@njit
def dt(phi):
    """Advance the wavefunction by one time step tau (second-order product formula).

    The kinetic operator is split into two parts acting on even (K1) and odd
    (K2) neighbouring pairs of grid points; the diagonal potential step (V) is
    sandwiched symmetrically between them: K1/2 K2/2 V K2/2 K1/2. The symmetric
    ordering makes the step second-order accurate in tau and exactly unitary.
    """
    for k in range(0,len(phi)-1,2):             #K1: kinetic half-step, even pairs
        phi[k:k+2] = np.dot(M,phi[k:k+2])
    for k in range(0,len(phi)-1,2):             #K2: kinetic half-step, odd pairs
        phi[k+1:k+3] = np.dot(M,phi[k+1:k+3])
    phi = Vx*phi                                #V: potential step (diagonal phase)
    for k in range(0,len(phi)-1,2):             #K2: kinetic half-step, odd pairs
        phi[k+1:k+3] = np.dot(M,phi[k+1:k+3])
    for k in range(0,len(phi)-1,2):             #K1: kinetic half-step, even pairs
        phi[k:k+2] = np.dot(M,phi[k:k+2])
    return phi

# Time-stepping loop. At ten evenly spaced times we store the wavefunction and
# save a snapshot of the probability density; "above" tracks the peak
# transmitted probability (the part of the packet beyond the barrier, x > 50).
above = 0
Psum = []                                       # transmitted probability vs time
Ptot = []                                       # total probability vs time (unitarity check)
for i in range(m):
    P = np.abs(Phi)**2 *Delta
    Ptot +=[np.sum(P)]
    if above < np.sum(P[506:]):
        above = np.sum(P[506:])
    Psum +=[np.sum(P[506:])]
    if (i)%(m//10) == 0:
        plt.figure()
        plt.plot(x,P, color=Color("RWTH"))
        plt.title("t = "+str(i*tau))
        plt.xlabel("x")
        plt.ylabel("P(x,t)")
        plt.xlim(0,100)
        plt.locator_params(nbins=8)
        plt.grid()
        if task ==1:
            plt.axvspan(50, 50.5, alpha=0.5, color="tab:green")
        plt.ylim(0,0.008)
        np.save(f"Data/TDSE_task{task}_times{int(i*tau)}.npy", Phi)
        plt.savefig(f"Plots/TDSE_task{task}_times{int(i*tau)}.png", dpi=150, bbox_inches="tight")
        plt.close()
        print(i)
    Phi = dt(Phi)
 
P = np.abs(Phi)**2 *Delta
Psum +=[np.sum(P[506:])]
Ptot +=[np.sum(P)]

# Save the transmitted- and total-probability time series. The "with"/"without"
# suffix refers to the presence of the barrier (task 1 vs task 0); plot_results.py
# reads these back to plot the transmitted probability versus time.
suffix = "with" if task == 1 else "without"
np.save(f"Data/Psum_{suffix}.npy", np.array(Psum))
np.save(f"Data/Ptot_{suffix}.npy", np.array(Ptot))

# Final probability-density snapshot.
plt.figure()
plt.plot(x, P, color=Color("RWTH"))
plt.xlabel("x")
plt.ylabel("P(x,t)")
plt.title(f"final state, t = {m*tau}")
plt.xlim(0, 100)
if task == 1:
    plt.axvspan(50, 50.5, alpha=0.5, color="tab:green")
plt.ylim(0, 0.015)
plt.grid()
plt.savefig(f"Plots/TDSE_task{task}_final.png", dpi=150, bbox_inches="tight")
plt.close()