#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 18 16:23:33 2023

@author: tamilarasan
"""

import numpy as np
import matplotlib.pyplot as plt
from numba import njit

# For fancy plots :)
plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "times new roman"
plt.rcParams["font.size"]   = "27"

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


# Do the plots
def plot_var(i0, m):
    t       = np.linspace(0, m*tau, m+1)
    var_arr = solve_product(m, i0)

    slope = (var_arr[-1]-var_arr[0])/(t[-1]-t[0])

    np.save(f"var_{i0}.npy", var_arr)
    setup(" ", "time t", r"$\Delta^{-2}$ var($x(t)$)")
    plt.plot(t, var_arr, lw=5)
    plt.plot(t, slope*t, lw=3, ls="--", color="black", label=f"slope = {slope:.2f}")
    plt.legend()
    plt.savefig(f"plots/var_{i0}.pdf")

# Similar function to solve_product(m, i0)
# Here we do plots of the array Phi vs. x 
def plot_N(m, i0):
    x           = np.arange(L)+1
    Phi0        = np.zeros(L)
    Phi0[i0-1]  = 1
    
    setup(" ", "position $x$", r"$N(x, t)$")
    plt.ylim(0,1)
    temp = Phi0
    for i in range(m+1):
        if (i)%(m//2) == 0:
            plt.plot(x, temp, label=rf"t = {i*tau:.2f}, $\sum \Phi$ = {np.sum(temp):.3f}")
        temp = Phi_A2(temp)
        temp = Phi_B(temp)
        temp = Phi_A2(temp)
    
    if i0==501:
        plt.xlim(400, 600)
    else:
        plt.xlim(1, 20)
        
    plt.legend()
    plt.savefig(f"plots/Phi_{i0}.pdf")
    plt.close()

# Define initial condition
i01     = int(L+1)//2
i02     = 1

# set time: t = m*tau
m = 10000

# time for the other plots
m2 = 60

# Generate plots
plot_var(i01, m)
plot_var(i02, m)

plot_N(m2, i01)
plot_N(m2, i02)