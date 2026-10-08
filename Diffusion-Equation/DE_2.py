#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Random-walk (Monte-Carlo) check of the diffusion equation.

Simulates n = 10000 independent particles performing unbiased random walks on
a lattice (fixed seed 1069 for reproducibility) and measures how the spatial
variance of the particle cloud grows with time. The result is compared with
the macroscopic diffusion-equation prediction var(x) = 2*D/Delta^2 * t,
demonstrating that the microscopic random walk reproduces the continuum law.

Output: figure to ``Plots/``.  Run from the project folder:  ``python DE_2.py``
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")   # non-interactive backend
import matplotlib.pyplot as plt

# plain mathtext labels, no LaTeX needed
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = "27"


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
#set the seed
rand.seed(1069)
#the constants given by the exercise
L = 1001
n = 10000
Delta = 0.1
D = 5  
#inital conditions
i0  = (L+1)//2
#array consisting the current position of each particle
N = np.ones(n)*i0
tau = 0.001
#time array
t = np.arange(11)


#a function returning the direction each particle will be moved in
def walk(n):
    return rand.choice([-1,1],n)

#function that does one iteration of t = t+1
def onet(N):
    F = N.copy()
    for i in range(int(1/tau)):
        F = F+walk(len(F))
    return F

#array below will contain the position of the particles at each time
N_total = [N]
#loops over all times
for T in t:
    N_total +=[onet(N_total[T])]

#function calculating the averages of x^p    
def x_mean_p(p, N):
    i_arr   = np.arange(1, L+1)
    return np.sum((i_arr-i0)**int(p) *N)/np.sum(N)

# Calculate the variance/Delta^
def var(N):
    return x_mean_p(2, N) - x_mean_p(1, N)

#array consisting the variances at each time t
Var = []
for i in range(11):
    #histogram is used to calculate the distribution of the particles at each position
    Var +=[var(np.histogram(N_total[i], bins = 1001, range=(0,1001))[0])]

#plots
plt.figure(figsize =(16 , 9))
plt.plot(t,Var, lw=5, marker="o", markersize=7, markerfacecolor="white", color=Color("RWTH"), label= "result random walk simulation")
#theoretical prediction
plt.plot (t , 2*D/(Delta **2) * t , lw =3 , color = "black" , ls = "--", label= "theory diffusion equation" )
plt.title("Random Walk")
plt.xlabel(r"time $t$")
plt.grid()
plt.legend()
plt.ylabel(r"Var$(x)$")
plt.savefig("Plots/task2.png", dpi=150, bbox_inches="tight")
plt.close()