#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Time-dependent Schrodinger equation for the 1D quantum harmonic oscillator.

Evolves a Gaussian wave packet in a harmonic potential V(x) = Omega^2 x^2 / 2
(units m = hbar = 1) using a matrix decomposition together with the
second-order product formula for the kinetic part. For several parameter sets
(frequency Omega, initial width sigma, initial centre x0) it tracks the mean
position <x(t)> and the variance Var(x(t)), and compares them with the exact
analytic expectation values. Runs are distributed over CPU cores with
``multiprocessing``; this is the expensive step.

The mean/variance arrays are written to ``Data/`` and the figures to
``Plots/``. To just plot already-stored snapshots, use
``plot_results.py`` instead.

Run from the project folder:  ``python QHO.py``
"""
#import important modules
import numpy as np
import matplotlib
matplotlib.use("Agg")   # non-interactive backend
import matplotlib.pyplot as plt
from numba import njit, vectorize, float64

from multiprocessing import Pool

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


SIM_COLORS = [Color(c) for c in
              ("RWTH", "Bordeaux", "Petrol", "Orange", "Grün", "Violett")]

pi      = np.pi


#parameters for the discretizations
Delta = 0.025
L = 1201
tau = 0.00025
m = 40000


#function that solves the system Input is an array consisting of [Omega, sigma, x0]
def solve(Input): 
    Omega, sigma, x0 = Input
    
    #Matrix M as explained in simulation model
    c = np.cos(tau/(4*Delta**2))
    s = 1j*np.sin(tau/(4*Delta**2))
    M = np.array([[c,s],[s,c]], dtype=np.complex128) 
    
    #discretization of space
    x = np.linspace(-15,15,L)
    
    #function that calculates initial wavepacket
    @njit
    def Phi0(x):
        return (pi*sigma**2)**(-1/4)*np.exp(-(x-x0)**2/(2*sigma**2))
    
    #the initial wavepacket
    Phi = np.array(Phi0(x)).astype(np.complex128) 
    
    #function for the the harmonic potential
    #function is vectorized to input array
    @njit
    def V(x):
        return 1/2*Omega**2*x**2
    Vvec= np.vectorize(V)    
    Vx = np.exp(-1j*tau*(Delta**(-2)+Vvec(x)))
    
    #function that calculates one time step tau evolution
    #input is a phi(t) output is phi(t+tau)
    @njit
    def dt(phi):
        for k in range(0,len(phi)-1,2):             #K1
            phi[k:k+2] = np.dot(M,phi[k:k+2])
        for k in range(0,len(phi)-1,2):             #K2
            phi[k+1:k+3] = np.dot(M,phi[k+1:k+3])
        phi = Vx*phi                                #V part
        for k in range(0,len(phi)-1,2):             #K2
            phi[k+1:k+3] = np.dot(M,phi[k+1:k+3])
        for k in range(0,len(phi)-1,2):             #K1
            phi[k:k+2] = np.dot(M,phi[k:k+2])
        return phi
    
    #array that will consist of mean values of x and x^2
    X = np.zeros(m+1)
    Xsq = np.zeros(m+1)
    #time at which the plots are plot
    tprint = [0,2,4,6,8,10]
    #array in which the probabilities at tprint will be saved
    Pt = []
    #doing m iterations
    for i in range(m):
        P = np.abs(Phi)**2
        if round(i*tau,5) in tprint:
            Pt+= [P]
        X[i] = np.sum(x*P*Delta)
        Xsq[i] = np.sum(x**2*P*Delta)
        Phi = dt(Phi)
    #after last iteration the last values (t=0) have to be saved seperately
    P = np.abs(Phi)**2
    Pt += [P]
    X[m] = np.sum(x*P*Delta)
    Xsq[m] = np.sum(x**2*P*Delta)
    np.save(f"Data/mean_params{Omega}{sigma}{x0}.npy", X)
    np.save(f"Data/var_params{Omega}{sigma}{x0}.npy", Xsq-X**2)
    return X,Xsq,Pt


    
#function that plots the values obtained
def plot(result,arg):
    Omega, sigma, x0 = arg
    legend = [[1,1,0],[1,1,1],[1,2,0],[2,1,1],[2,2,2]]         #for which plot the universal legend is attached
    X,Xsq,P = result                  #solutions from function solve(arg)
    t = np.linspace(0,m*tau,m+1)    #timearray
    x = np.linspace(-15,15,L)       #position discretization
    
    #theoretical expectation of x and x^2
    xth = x0*np.cos(Omega*t)
    xsqth = 1/(2*Omega**2*sigma**2)*(np.sin(Omega*t))**2+ 1/2*(sigma**2+2*x0**2)*(np.cos(Omega*t))**2
    tprint = [0,2,4,6,8,10] #at which time the probabilities are plot
    
    #plotting the averages (both theory and simulation)
    plt.grid()
    plt.plot(t, X, label =r" $\langle x(t)_{sim} \rangle$", color =Color("RWTH"), lw = 2.5)
    plt.plot(t, xth, label = r"$\langle x(t)_{theo} \rangle$", color = "black", ls = "-.", lw = 1.5)
    plt.plot(t, Xsq-X**2,label =r"$Var(x_{sim})$", color =Color("Bordeaux"), lw = 2.5)
    plt.plot(t, xsqth-xth**2, label = r"$Var(x_{theo})$", color = "black",ls = "--", lw = 1.5)
    plt.xlabel("t")
    plt.ylabel("Average and Variance")
    plt.title(r"$\Omega$ = "+str(Omega)+r", $\sigma$ = "+str(sigma)+r", $x_0$ = "+str(x0))
    if arg in legend:
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(f"Plots/Omega{Omega}_sigma{sigma}_x0{x0}_Averages.png", dpi=150, bbox_inches="tight")
    plt.close()
    
    #plotting the differences between theory and simulation
    plt.grid()
    plt.plot(t, X-xth, label = r"$\langle x(t)_{sim} \rangle - \langle x(t)_{theo} \rangle$", color = Color("RWTH"), lw = 2.5)
    plt.plot(t, Xsq-X**2 - (xsqth-xth**2),label =r"$Var(x_{sim}) -Var(x_{theo})$", color =Color("Bordeaux"), lw = 2.5)
    plt.xlabel("t")
    plt.ylabel("Difference of Average and Variance")
    plt.title(r"$\Omega$ = "+str(Omega)+r", $\sigma$ = "+str(sigma)+r", $x_0$ = "+str(x0))
    if arg in legend:
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(f"Plots/Omega{Omega}_sigma{sigma}_x0{x0}_Averages_expect.png", dpi=150, bbox_inches="tight")
    plt.close()
    
    #plotting the probability density from -5 to 5
    for i in range(6):
        plt.plot(x,P[i]*Delta, color=SIM_COLORS[i], label = "t = "+str(tprint[i]))
    plt.title(r"$\Omega$ = "+str(Omega)+r", $\sigma$ = "+str(sigma)+r", $x_0$ = "+str(x0))
    plt.xlabel("x")
    plt.xlim(-5,5)
    plt.grid()
    plt.ylabel(r"Probability $|\Phi|^2\cdot \Delta$")
    if arg in legend:
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(f"Plots/Omega{Omega}_sigma{sigma}_x0{x0}_Probabilities.png", dpi=150, bbox_inches="tight")
    plt.close()
    
    #plotting the probability density from -15 to 15
    for i in range(6):
        plt.plot(x,P[i]*Delta, color=SIM_COLORS[i], label = "t = "+str(tprint[i]))
    plt.title(r"$\Omega$ = "+str(Omega)+r", $\sigma$ = "+str(sigma)+r", $x_0$ = "+str(x0))
    plt.xlabel("x")
    plt.xlim(-15,15)
    plt.grid()
    plt.ylabel(r"Probability $|\Phi|^2\cdot \Delta$")
    if arg in legend:
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(f"Plots/Omega{Omega}_sigma{sigma}_x0{x0}_Probabilities_full.png", dpi=150, bbox_inches="tight")
    plt.close()
    return 0

#multiprocessing to calculate all the initial values at the same time
if __name__ == '__main__':
    pool = Pool()
    args = [[1,1,0],[1,1,1],[1,2,0],[2,1,1],[2,2,2]]
    Sol = pool.map(solve, args)

    for i in range(5):
        plot(Sol[i],args[i])
