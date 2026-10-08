#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Simulate 1D light propagation through a glass plate with the Yee (FDTD) scheme.

Solves the 1D Maxwell equations on a staggered Yee grid, where the electric
field E_z and magnetic field H_y are offset by half a cell in both space and
time. A Gaussian-modulated sinusoidal source launches a wave that partially
reflects and transmits at a dielectric slab (refractive index n_d = 1.46),
with absorbing layers (conductivity sigma, magnetic loss sigma*) at the
boundaries. The script produces snapshots of E_z(t,x) for a thin and a thick
glass plate and illustrates the Courant stability limit by contrasting a
stable time step (tau = 0.9*Delta) with an unstable one (tau = 1.05*Delta).

Output: PNG snapshots to ``Plots/``.  Run from the project folder:
``python EM_Waves.py``  (Numba-accelerated; the full set of runs takes a
minute or two on first run while Numba compiles.)
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")   # non-interactive backend
import matplotlib.pyplot as plt

#numba for faster code compilation
from numba import njit, vectorize, float64

# plain mathtext labels, no LaTeX needed
plt.rcParams["font.family"] = "serif"


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


#parameters
lamb = 1
grid = 50
Delta = lamb/grid

tau1 = 0.9*Delta
tau2 = 1.05*Delta

X = 100*lamb
L = 5000
f = 1/lamb
w = 2*pi*f
m = 10000

#left and right insulator boundary
leftinsu = 6*lamb
rightinsu = (L*Delta -6*lamb)


#function to calculate sigma and sigma*
@vectorize([float64(float64)])
def sigma(x):
    if leftinsu < x< rightinsu:
        return 0
    else:
        return 1
    
#function to calculate epsilon
#takes in position as well as the end of the right glass
@vectorize([float64(float64,float64)])
def epsilon(x,rightglass):
    n_d = 1.46
    leftglass = (L*Delta/2)
    if leftglass <= x < rightglass:
        return n_d**2
    else:
        return 1

#function defining the mu    
@vectorize([float64(float64)])    
def mu(x):
    return 1

#position of source
x_s = 20*lamb
i_s = x_s/Delta
#source function, takes in time t = n\tau
@njit
def J_s(t):
    return np.sin(w*t)*np.exp(-((t-30)/10)**2)

#positions        
x = np.linspace(0,X,L+1)


#the four coefficient introduced in the lecture
@njit
def A_l12(l,tau):
    x = (l+1/2)*Delta
    num = 1-(sigma(x)*tau/(2*mu(x)))
    den = 1+(sigma(x)*tau/(2*mu(x)))
    return num/den
@njit
def C_l(l,tau,rightglass):
    x = l*Delta
    num = 1-(sigma(x)*tau/(2*epsilon(x,rightglass)))
    den = 1+(sigma(x)*tau/(2*epsilon(x,rightglass)))
    return num/den
    
@njit
def B_l12(l,tau):
    x = (l+1/2)*Delta
    num = (tau/mu(x))
    den = 1+(sigma(x)*tau/(2*mu(x)))
    return num/den

@njit
def D_l(l,tau,rightglass):
    x = l*Delta
    num = (tau/epsilon(x,rightglass))
    den = 1+(sigma(x)*tau/(2*epsilon(x,rightglass)))
    return num/den


#update function for H
@njit
def H_n1_l12(l,n, tau, E,H):    
    Bl12 = B_l12(l,tau)
    Al12 = A_l12(l,tau)
    return  Bl12[:-1]*(E[1:]-E[:-1])/Delta + Al12[:-1]*H

#update function for E
@njit
def E_n12_l(l,n, tau, E,H,rightglass):
    Dl = D_l(l,tau,rightglass)
    Cl = C_l(l,tau,rightglass)

    E[1:-1] = Dl[1:-1]*(H[1:]-H[:-1])/Delta + Cl[1:-1]*E[1:-1]
    E[int(i_s)] -=  Dl[int(i_s)]*J_s(n*tau1)
    return E
    
#function that calculates the fields for a given n_max, tau
#glass takes in arguements 1 for "thin" and 0 for "thick"
def maxwell(n_max,tau, glass):
    
    #defines the boundaries of the glass depending on the chosen glass
    leftglass = (L*Delta/2)
    if glass == 1:
        rightglass = (L*Delta/2 + 2*lamb)
    else:
        rightglass = L*Delta
    
    #defining the grid and the initial field values
    l = np.arange(0,L+1)
    H = np.zeros(L)
    E = np.zeros(L+1)
    
    #array only needed if the reflection coefficients are calculated
    E_incoming = []
    E_reflected = []
    #iteration over time
    for n in range(n_max):
        #updating the fields
        E = E_n12_l(l,n,tau,E,H,rightglass)
        H = H_n1_l12(l,n,tau,E,H)
        #finds maxima in specific time window (only for "thick" glas)
        #needed to calc R
        if glass == 0 and n_max >5000:
            if 1700 < n < 2000:
                E_incoming += [np.max((E[1000:2000])**2)]
            elif 4700 < n <4950:
                E_reflected += [np.max((E[1000:2000])**2)]
    #prints the Reflection coefficient            
    if glass == 0 and n_max >5000:
        R = np.mean(E_reflected)/np.mean(E_incoming)
        print("Reflection Coefficient: "+ str(R))
    return E, H
        
#function needed to set certain plotting arguments
def setup(title, xlabel, ylabel):
    plt.figure(figsize=(16, 9))
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(linestyle="--")
    plt.tight_layout()
    
#function generating the plots for given values n_max, tau, glass   
def maxwell_plots(n_max, tau, glass, set_legend):
    E, H    = maxwell(n_max, tau, glass)
    x       = np.linspace(0, X, L+1)

    leftglass = (L*Delta/2)
    
    if set_legend:
        plt.rcParams["font.size"]   = "25"
    else:
        plt.rcParams["font.size"]   = "42"
    
    if tau == tau1:
        savetau = 0.90
    else:
        savetau = 1.05
    
    if glass == 1.:
        rightglass = (L*Delta/2 + 2*lamb)
    else:
        rightglass = X


    setup(rf"$n={n_max}$, $\tau={savetau} \Delta$", "Position $x$", "$E_z(t,x)$")

    lw = 3
    plt.plot(x, E, c=Color("RWTH"), lw = lw)
    plt.scatter(x_s, 0, marker="o", color="red", label="source", s=100)
    
    plt.vlines(leftglass, -0.2, 0.2, color="tab:green", lw=lw)
    plt.vlines(rightglass, -0.2, 0.2, color="tab:green", lw=lw)
    plt.vlines(leftinsu, -0.2, 0.2, color="grey", lw=lw)
    plt.vlines(rightinsu, -0.2, 0.2, color="grey", lw=lw)
    
    plt.axvspan(leftglass, rightglass, alpha=0.5, color="tab:green", label="glass")
    plt.axvspan(0, leftinsu, alpha=0.5, color="grey", label="insulator")
    plt.axvspan(rightinsu, X, alpha=0.5, color="grey")
    
    plt.locator_params(axis='y', nbins=5)
    plt.xlim(0, 100)
    
    if tau == tau1:
        plt.ylim(-0.015, 0.015)
    
    if set_legend:
        plt.legend(bbox_to_anchor=(1.01, 0.65), fancybox=True, shadow=True, fontsize=25)
        
    plt.tight_layout()
    if glass == 1.:
        plt.savefig(f"Plots/maxwell_tau{savetau}_thinglass_nmax{n_max}.png", dpi=150, bbox_inches="tight")
    else:
        plt.savefig(f"Plots/maxwell_tau{savetau}_thickglass_nmax{n_max}.png", dpi=150, bbox_inches="tight")
    plt.close()


#values at which we are interested to generate plot
n_max0 = 0
n_max1 = 2500
n_max2 = 3500
n_max3 = 3510
n_max4 = 4500
n_max5 = 20000


#generation of plots
maxwell_plots(n_max0, tau1, 0., False)
maxwell_plots(n_max1, tau1, 0., False)
maxwell_plots(n_max2, tau1, 0., False)
maxwell_plots(n_max3, tau1, 0., False)
maxwell_plots(n_max4, tau1, 0., False)
maxwell_plots(n_max5, tau1, 0., False)

maxwell_plots(500, tau2, 0., False)

maxwell_plots(500, tau2, 1., False)

maxwell_plots(n_max0, tau1, 1., True)
maxwell_plots(n_max1, tau1, 1., False)
maxwell_plots(n_max2, tau1, 1., False)
maxwell_plots(n_max3, tau1, 1., False)
maxwell_plots(n_max4, tau1, 1., False)
maxwell_plots(n_max5, tau1, 1., False)