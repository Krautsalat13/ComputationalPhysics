#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Solve the 1D diffusion equation with a second-order product-formula scheme.

The density N(x,t) obeys dN/dt = D d^2N/dx^2. The spatial Laplacian is split
into two block-diagonal operators A and B, and the time evolution over one
step tau is approximated by the symmetric (second-order) product formula
    exp(alpha*tau*(A+B)) ~ exp(alpha*tau*A/2) exp(alpha*tau*B) exp(alpha*tau*A/2),
each factor acting on 2x2 blocks in closed form. The script evolves a delta
initial condition, plots the density profile N(x,t), and checks that the
variance grows linearly in time with the expected slope 2D/Delta^2.

Outputs: variance arrays to ``Data/``, figures to ``Plots/``.
Run from the project folder:  ``python product_formula.py``
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")   # non-interactive backend
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from numba import njit

# plain mathtext labels, no LaTeX needed
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"]   = "27"


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

def setup(title, xlabel, ylabel):
    plt.figure(figsize=(16, 9))
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xticks()
    plt.yticks()
    plt.locator_params(nbins=6)
    plt.grid(ls="--")

# Define Parameters
L       = 1001
Delta   = 0.1
tau     = 0.001

D       = 1
alpha   = -D/Delta**2

# Calculate the mean of x**p/Delta**2
@njit  
def x_mean_p(p, phi, i0):
    i_arr   = np.arange(1, L+1)
    return np.sum((i_arr-i0)**int(p)/np.sum(phi) *phi)

# Calculate variance/Delta**2
@njit
def var(phi, i0):
    return x_mean_p(2, phi, i0) - x_mean_p(1, phi, i0)**2

# Define the exponential Block matrices
eA2 = np.array([[1+np.exp(alpha*tau), 1-np.exp(alpha*tau)],[1-np.exp(alpha*tau), 1+np.exp(alpha*tau)]])/2
eB  = np.array([[1+np.exp(2*alpha*tau), 1-np.exp(2*alpha*tau)],[1-np.exp(2*alpha*tau), 1+np.exp(2*alpha*tau)]])/2

# Calculate the matrix-vector multiplication eA2*Phi_in
@njit
def Phi_A2(Phi_in):
    Phi = np.zeros(L)
    for i in range(L//2+1):
        if i-L//2 == 0:
            Phi[L-1] = np.exp(alpha*tau/2)*Phi_in[L-1]
        else:
            temp            = Phi_in[2*i:2*i+2]
            Phi[2*i:2*i+2]  = np.dot(eA2, temp)
    return Phi

# Calculate the matrix-vector multiplication eB*Phi_in
@njit
def Phi_B(Phi_in):
    Phi = np.zeros(L)
    for i in range(L//2+1):
        if i == 0:
            Phi[0] = np.exp(alpha*tau)*Phi_in[0]
        else:
            temp = Phi_in[2*i-1:2*i+1]
            Phi[2*i-1:2*i+1] = np.dot(eB, temp)
    return Phi

# do the three matrix-vector multiplications subsequently 3 times and iterate m times 
# to obtain the soluztion at time t = m*tau
# After each steps calculate the variance
@njit
def solve_product(m, i0):
    # Define initian condition 
    Phi0        = np.zeros(L)
    Phi0[i0-1]  = 1
    
    var_arr = np.zeros(m+1)
    temp    = Phi0
    for i in range(m):
        temp = Phi_A2(temp)
        temp = Phi_B(temp)
        temp = Phi_A2(temp)
        var_arr[i+1] = var(temp, i0)
    return var_arr


# Plot the variance growth for one initial condition (label: "centre"/"boundary").
def plot_variance(i0, m, label):
    t       = np.linspace(0, m*tau, m+1)
    var_arr = solve_product(m, i0)

    slope = (var_arr[-1]-var_arr[0])/(t[-1]-t[0])

    np.save(f"Data/variance_{label}.npy", var_arr)
    setup(" ", "time t", r"$\Delta^{-2}$ var($x(t)$)")
    plt.plot(t, var_arr, lw=5, color=Color("RWTH"), label="simulation")
    plt.plot(t, slope*t, lw=3, ls="--", color="black", label=f"slope = {slope:.2f}")
    plt.legend()
    plt.savefig(f"Plots/variance_{label}.png", dpi=150, bbox_inches="tight")
    plt.close()

# Plot the density profile N(x, t) at a few times for one initial condition.
def plot_profile(m, i0, label):
    x           = np.arange(L)+1
    Phi0        = np.zeros(L)
    Phi0[i0-1]  = 1

    setup(" ", "position $x$", r"$N(x, t)$")
    plt.ylim(0,1)
    temp = Phi0
    ci = 0
    for i in range(m+1):
        if (i)%(m//2) == 0:
            plt.plot(x, temp, color=SIM_COLORS[ci],
                     label=rf"t = {i*tau:.2f}, $\sum \Phi$ = {np.sum(temp):.3f}")
            ci += 1
        temp = Phi_A2(temp)
        temp = Phi_B(temp)
        temp = Phi_A2(temp)

    if label == "centre":
        plt.xlim(486, 516)      # zoom on the spreading peak centred at x=501
    else:
        plt.xlim(1, 14)         # boundary start: show the decaying profile
    # drop the lowest x tick so its label does not collide with the y-axis
    plt.gca().xaxis.set_major_locator(MaxNLocator(nbins=6, prune="lower"))

    plt.legend()
    plt.savefig(f"Plots/profile_{label}.png", dpi=150, bbox_inches="tight")
    plt.close()

# Define initial condition
i01     = int(L+1)//2
i02     = 1

# set time: t = m*tau
m = 10000

# time for the other plots
m2 = 60

# Generate plots. i01 is the centre of the chain, i02 is the boundary.
plot_variance(i01, m, "centre")
plot_variance(i02, m, "boundary")

plot_profile(m2, i01, "centre")
plot_profile(m2, i02, "boundary")