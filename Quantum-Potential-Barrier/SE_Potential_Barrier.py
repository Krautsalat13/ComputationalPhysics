#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 10:41:51 2023

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

pi      = np.pi  
rand    = np.random

sigma = 3
q = 1
x0 = 20
task = 1    #0 for first, 1 for second
Delta = 0.1
L = 1001
tau = 0.001
m = 50000
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

above = 0
Psum = []
Ptot = []
for i in range(m):
    P = np.abs(Phi)**2 *Delta
    Ptot +=[np.sum(P)]
    if above < np.sum(P[506:]):
        above = np.sum(P[506:])
    Psum +=[np.sum(P[506:])]
    if (i)%(m//10) == 0:
        plt.plot(x,P)
        plt.title("t = "+str(i*tau))
        plt.xlabel("x")
        plt.ylabel("P(x,t)")
        plt.xlim(0,100)
        plt.locator_params(nbins=8)
        plt.grid()
        if task ==1:
            plt.axvspan(50, 50.5, alpha=0.5, color="tab:green")
        plt.ylim(0,0.015)
        plt.ylim(0,0.008)
        np.save(f"Data/TDSE_task{task}_times{int(i*tau)}.npy", Phi)
        plt.show()
        print(i)
    Phi = dt(Phi)
 
P = np.abs(Phi)**2 *Delta       
Psum +=[np.sum(P[506:])]
Ptot +=[np.sum(P)]      
plt.plot(x,P)
plt.xlim(0,100)
#plt.axvspan(50, 50.5, alpha=0.5, color="tab:green")
plt.ylim(0,0.015)
plt.show()