#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 13:20:53 2023

@author: tamilarasan
"""
import numpy as np
import matplotlib.pyplot as plt

#numba for faster code compilation
from numba import njit, vectorize, float64

#for fancier plots
import os 
os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin'

plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "times new roman"


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
    plt.plot(x, E, c="tab:blue", lw = lw)
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
        plt.savefig(f"Plots/maxwell_tau{savetau}_thinglass_nmax{n_max}.pdf")
    else:
        plt.savefig(f"Plots/maxwell_tau{savetau}_thickglass_nmax{n_max}.pdf")
    plt.show()


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