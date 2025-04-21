#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 13:20:53 2023

@author: tamilarasan
"""
import numpy as np
import matplotlib.pyplot as plt

from numba import njit, vectorize, float64

import os 
os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin'

plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "times new roman"
plt.rcParams["font.size"] = "18"

plt.rcdefaults()
pi      = np.pi  

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




leftcond = 6*lamb
rightcond = (L*Delta -6*lamb)

glass = "thick"

leftglass = (L*Delta/2)
if glass == "thin":
    rightglass = (L*Delta/2 + 2*lamb)
else:
    rightglass = L*Delta

@vectorize([float64(float64)])
def sigma(x):
    if leftcond < x< rightcond:
        return 0
    else:
        return 1

@vectorize([float64(float64)])
def epsilon(x):
    n_d = 1.46
    if leftglass <= x < rightglass:
        return n_d**2
    else:
        return 1
 
@vectorize([float64(float64)])    
def mu(x):
    return 1


x_s = 20*lamb
i_s = x_s/Delta
@njit
def J_s(t):
    return np.sin(w*t)*np.exp(-((t-30)/10)**2)
        
x = np.linspace(0,X,1000)

@njit
def A_l12(l,tau):
    x = (l+1/2)*Delta
    num = 1-(sigma(x)*tau/(2*mu(x)))
    den = 1+(sigma(x)*tau/(2*mu(x)))
    return num/den
@njit
def C_l(l,tau):
    x = l*Delta
    num = 1-(sigma(x)*tau/(2*epsilon(x)))
    den = 1+(sigma(x)*tau/(2*epsilon(x)))
    return num/den
    
@njit
def B_l12(l,tau):
    x = (l+1/2)*Delta
    num = (tau/mu(x))
    den = 1+(sigma(x)*tau/(2*mu(x)))
    return num/den

@njit
def D_l(l,tau):
    x = l*Delta
    num = (tau/epsilon(x))
    den = 1+(sigma(x)*tau/(2*epsilon(x)))
    return num/den



@njit
def H_n1_l12(l,n, tau, E,H):    
    Bl12 = B_l12(l,tau)
    Al12 = A_l12(l,tau)
    return  Bl12[:-1]*(E[1:]-E[:-1])/Delta + Al12[:-1]*H

@njit
def E_n12_l(l,n, tau, E,H):
    Dl = D_l(l,tau)
    Cl = C_l(l,tau)

    E[1:-1] = Dl[1:-1]*(H[1:]-H[:-1])/Delta + Cl[1:-1]*E[1:-1]
    E[int(i_s)] -=  Dl[int(i_s)]*J_s(n*tau1)
    return E
    

def maxwell(n_max,tau):
    l = np.arange(0,L+1)
    H = np.zeros(L)
    E = np.zeros(L+1)
    k = 0
    E_incoming = []
    E_reflected = []
    for n in range(n_max):
        E = E_n12_l(l,n,tau,E,H)
        H = H_n1_l12(l,n,tau,E,H)
        if glass == "thick":
            if 1700 < n < 2000:
                E_incoming += [np.max((E[1000:2000])**2)]
            elif 4700 < n <4950:
                E_reflected += [np.max((E[1000:2000])**2)]
        if (n)%1000 ==0:
            plt.grid()
            plt.plot(l*Delta,E, linewidth = 1)
            #plt.title(r"$t = {:.2f}$".format((n+1/2)*tau))
            plt.title(r"$t = {:.2f}$".format(n))
            plt.axvspan(leftglass ,rightglass, color="tab:green", alpha=.5)
            plt.axvspan(0 ,leftcond, color="grey", alpha=.5)
            plt.axvspan(rightcond ,L*Delta, color="grey", alpha=.5)
            plt.ylim(-0.0178, 0.0178)
            plt.xlim(0,100)
            plt.scatter(x_s, 0, marker="o", color="red")
            #plt.legend()
            plt.savefig('{}.png'.format(n), transparent=True, dpi=300)
            plt.show()
            
    if glass == "thick":
        R = np.mean(E_reflected)/np.mean(E_incoming)
        print("Reflection Coefficient: "+ str(R))
    return E, H
        

