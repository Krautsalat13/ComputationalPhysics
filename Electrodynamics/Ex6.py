#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 13:20:53 2023

@author: tamilarasan
"""
import numpy as np
import matplotlib.pyplot as plt
import csv

import os 
os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin'

plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "times new roman"
plt.rcParams["font.size"] = "18"

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


l_max = L

def sigma(x):
    def sig(x):
        if 6*lamb < x< (L*Delta -6*lamb):
            return 0
        else:
            return 1
    return np.vectorize(sig)(x)    


def epsilon(x, n_d = 1.46):  
    def eps(x, n_d):
        if (L*Delta/2) <= x < (L*Delta/2 + 2*lamb):
            return n_d**2
        else:
            return 1
    return np.vectorize(eps)(x,n_d)
    
def mu(x):
    return 1

x_s = 20*lamb
i_s = x_s/Delta

def J_s(t):
    return np.sin(w*t)*np.exp(-((t-30)/10)**2)
        
x = np.linspace(0,X,1000)

def A_l12(l,tau):
    x = (l+1/2)*Delta
    num = 1-(sigma(x)*tau/(2*mu(x)))
    den = 1+(sigma(x)*tau/(2*mu(x)))
    return num/den

def C_l(l,tau):
    x = l*Delta
    num = 1-(sigma(x)*tau/(2*epsilon(x)))
    den = 1+(sigma(x)*tau/(2*epsilon(x)))
    return num/den
    

def B_l12(l,tau):
    x = (l+1/2)*Delta
    num = (tau/mu(x))
    den = 1+(sigma(x)*tau/(2*mu(x)))
    return num/den

def D_l(l,tau):
    x = l*Delta
    num = (tau/epsilon(x))
    den = 1+(sigma(x)*tau/(2*epsilon(x)))
    return num/den


def H_n1_l12(l,n, tau, E,H):
    return  B_l12(l,tau)[:-1]*(E[1:]-E[:-1])/Delta + A_l12(l,tau)[:-1]*H

def E_n12_l(l,n, tau, E,H):
    Eold = D_l(l,tau)[1:-1]*(H[1:]-H[:-1])/Delta + C_l(l,tau)[1:-1]*E[1:-1]
    Eold[int(i_s)] -=  D_l(x_s,tau)*J_s(n*tau)
    return Eold
    
def maxwell(n_max):
    H = np.zeros(l_max)
    E = np.zeros(l_max+1)
    l = np.arange(0,L+1)
    k = 0
    for n in range(n_max):
        E[1:-1] = E_n12_l(l,n,tau1,E,H)
        H = H_n1_l12(l,n,tau1,E,H)
        
        if k ==10:
            plt.plot(l*Delta,E, linewidth = 1)
            plt.plot(l*Delta, sigma(l*Delta)-0.5, color="tab:blue")
            plt.plot(l*Delta, epsilon(l*Delta)-1.1, color="tab:green")
            plt.ylim(-0.015, 0.015)
            plt.scatter(x_s, 0, marker="o", color="red")
            plt.show()
            k = 0
        k+=1
    return E,H
        
  
maxwell(10000)           