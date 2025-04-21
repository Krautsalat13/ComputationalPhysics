#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 14 16:35:58 2023

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
rand.seed(1556)


N = 50
T = np.linspace(0.2,4,20)
T = 0.2
N_samples = 10000
N_wait = N_samples/10

def S2d(n):
    return rand.choice([-1,1],(n,n))
def E2d(S,n):
    E = 0
    for i in range(n-1):
        for j in range(n-1):
            E+= S[i,j]*(S[i+1,j] + S[i,j+1])
    last = n-1
    E += S[last,last]*(S[last-1,last]+S[last,last-1])
    return -E  
def E2d(S,n):
    E = 0
    for i in range(n-1):
        for j in range(n-1):
            E+= S[i,j]*(S[i+1,j] + S[i,j+1])
    last = n-1
    E += S[last,last]*(S[last-1,last]+S[last,last-1])
    return -E  
    


def flip(S,N):
    i,j = int(rand.random()*N),int(rand.random()*N)
    S[i,j] *= -1
    return S

def energycrit(Sold,N,beta):
    Eold = E2d(Sold,N)
    Snew = np.copy(Sold)
    Snew = flip(Snew,N)
    Enew = E2d(Snew,N)
    deltaE = Enew-Eold
    q = np.exp(-beta*deltaE)
    r = rand.random()
    return np.where(q>r, Snew,Sold)

Sold = S2d(N)
plt.imshow(Sold)
plt.show()
beta = 1/T
E = np.array([])
for i in range(int(N_wait + N_samples)):
    Sold = energycrit(Sold,N,beta)
    #plt.imshow(Sold)
    #plt.show()
    E =np.append(E,E2d(Sold,N))
plt.imshow(Sold) 
    


end = time.time()
print(end-start)
    