#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 15:46:58 2023

@author: tamilarasan
"""
#import important modules
import numpy as np
import matplotlib.pyplot as plt
from numba import njit, vectorize, float64

from multiprocessing import Pool

#some plotting arguments for better visualization
import os 
os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin'

plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "times new roman"
plt.rcParams["font.size"] = "18"

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
    np.save(f"DataNew/tami_mean_params{Omega}{sigma}{x0}.npy", X)
    np.save(f"DataNew/tami_var_params{Omega}{sigma}{x0}.npy", Xsq-X**2)
    return X,Xsq,Pt


    
#function that plots the values obtained
def plot(task,arg):
    Omega, sigma, x0 = arg
    legend = [[1,1,0],[1,1,1],[1,2,0],[2,1,1],[2,2,2]]         #for which plot the universal legend is attached
    X,Xsq,P = task                  #solutions from function solve(arg)
    t = np.linspace(0,m*tau,m+1)    #timearray
    x = np.linspace(-15,15,L)       #position discretization
    
    #theoretical expectation of x and x^2
    xth = x0*np.cos(Omega*t)
    xsqth = 1/(2*Omega**2*sigma**2)*(np.sin(Omega*t))**2+ 1/2*(sigma**2+2*x0**2)*(np.cos(Omega*t))**2
    tprint = [0,2,4,6,8,10] #at which time the probabilities are plot
    
    #plotting the averages (both theory and simulation)
    plt.grid()
    plt.plot(t, X, label =r" $\langle x(t)_{sim} \rangle$", color ="tab:cyan", lw = 2.5)
    plt.plot(t, xth, label = r"$\langle x(t)_{theo} \rangle$", color = "black", ls = "-.", lw = 1.5)
    plt.plot(t, Xsq-X**2,label =r"$Var(x_{sim})$", color ="tab:orange", lw = 2.5)
    plt.plot(t, xsqth-xth**2, label = r"$Var(x_{theo})$", color = "black",ls = "--", lw = 1.5)
    plt.xlabel("t")
    plt.ylabel("Average and Variance")
    plt.title(r"$\Omega$ = "+str(Omega)+r", $\sigma$ = "+str(sigma)+r", $x_0$ = "+str(x0))
    if arg in legend:
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(f"plot/Omega{Omega}_sigma{sigma}_x0{x0}_Averages.pdf",bbox_inches='tight')
    plt.show()
    
    #plotting the differences between theory and simulation
    plt.grid()
    plt.plot(t, X-xth, label = r"$\langle x(t)_{sim} \rangle - \langle x(t)_{theo} \rangle$", color = "tab:cyan", lw = 2.5)
    plt.plot(t, Xsq-X**2 - (xsqth-xth**2),label =r"$Var(x_{sim}) -Var(x_{theo})$", color ="tab:orange", lw = 2.5)
    plt.xlabel("t")
    plt.ylabel("Difference of Average and Variance")
    plt.title(r"$\Omega$ = "+str(Omega)+r", $\sigma$ = "+str(sigma)+r", $x_0$ = "+str(x0))
    if arg in legend:
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(f"plot/Omega{Omega}_sigma{sigma}_x0{x0}_Averages_expect.pdf",bbox_inches='tight')
    plt.show()
    
    #plotting the probability density from -5 to 5
    for i in range(6):
        plt.plot(x,P[i]*Delta, label = "t = "+str(tprint[i]))
    plt.title(r"$\Omega$ = "+str(Omega)+r", $\sigma$ = "+str(sigma)+r", $x_0$ = "+str(x0))
    plt.xlabel("x")
    plt.xlim(-5,5)
    plt.grid()
    plt.ylabel(r"Probability $|\Phi|^2\cdot \Delta$")
    if arg in legend:
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(f"plot/Omega{Omega}_sigma{sigma}_x0{x0}_Probabilities.pdf",bbox_inches='tight')
    plt.show()
    
    #plotting the probability density from -15 to 15
    for i in range(6):
        plt.plot(x,P[i]*Delta, label = "t = "+str(tprint[i]))
    plt.title(r"$\Omega$ = "+str(Omega)+r", $\sigma$ = "+str(sigma)+r", $x_0$ = "+str(x0))
    plt.xlabel("x")
    plt.xlim(-15,15)
    plt.grid()
    plt.ylabel(r"Probability $|\Phi|^2\cdot \Delta$")
    if arg in legend:
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig(f"plot/Omega{Omega}_sigma{sigma}_x0{x0}_Probabilities_full.pdf",bbox_inches='tight')
    plt.show()
    return 0

#multiprocessing to calculate all the initial values at the same time
if __name__ == '__main__':
    pool = Pool()
    args = [[1,1,0],[1,1,1],[1,2,0],[2,1,1],[2,2,2]]
    Sol = pool.map(solve, args)

    for i in range(5):
        plot(Sol[i],args[i])
