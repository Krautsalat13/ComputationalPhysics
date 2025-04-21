#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 13:49:29 2023

@author: tamilarasan
"""


import numpy as np
import matplotlib.pyplot as plt
import time as time

import os 
os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin'

plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "times new roman"
plt.rcParams["font.size"] = "18"
start = time.time()
pi      = np.pi  
rand    = np.random
rand.seed(1069)


N = np.array([10,100,1000])
T = np.linspace(0.2,4,20)
N_samples = np.array([1000,10000])
N_samples = N_samples[0]


part = 100

def S1d(n,M):
    return rand.choice([-1,1],(M,n))

def E1d(S):
    return -np.sum(S[:,:-1]*S[:,1:], axis = 1)

def flip(S,n):
    i = (rand.random(size = N_samples)*n).astype(int)
    S[np.arange(len(S)), i] *= -1
    return S

def U(n,T):    
    Sold = S1d(n,N_samples)
    beta = 1/T
    for i in range(N_samples):
        Eold = E1d(Sold)
        Snew = np.copy(Sold)
        Snew = flip(Snew,n)
        Enew = E1d(Snew)
        deltaE = Enew-Eold
        q = np.exp(-beta*deltaE)
        r = rand.random(N_samples)
        Sold = np.where(q[:,np.newaxis]>r[:,np.newaxis], Snew,Sold)
        
    return np.sum(E1d(Sold))/N_samples
U_vec = np.vectorize(U)

U10 = U_vec(N[0],T)
plt.plot(T,U10/10)
plt.plot(T, -(N[0]-1)/N[0] *np.tanh(1/T))
"""
U10 = U_vec(N[0],T)
U100 = U_vec(N[1],T)
U1000 = U_vec(N[2],T)

plt.plot(T,U10/10)
plt.plot(T, -(N[0]-1)/N[0] *np.tanh(1/T))

plt.plot(T,U100/100)
plt.plot(T, -(N[1]-1)/N[1] *np.tanh(1/T))

plt.plot(T,U1000/1000)
plt.plot(T, -(N[2]-1)/N[2] *np.tanh(1/T))
"""


end = time.time()
print(end-start)
    