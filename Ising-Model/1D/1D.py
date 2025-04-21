#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 13:49:29 2023

@author: tamilarasan
"""


import numpy as np
import matplotlib.pyplot as plt
import time as time
import csv

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
T = np.linspace(0.2,4.2,21)
T[-1] = 4
N_samples = np.array([1000,10000])



def S1d(n):
    return rand.choice([-1,1],n)

def E1d(S):
    return -np.sum(S[:-1]*S[1:])

def energycrit(Sold,beta):
    n = Sold.shape[0]
    i = int(rand.random()*n)
    if (i ==0 or i == n-1):
        if (i == 0):
            deltaE = 2*Sold[0]* Sold[1]
        else:
            deltaE = 2*Sold[n-2]* Sold[n-1]
    else:
        deltaE = 2*Sold[i]*(Sold[i-1]+Sold[i+1])
    q = np.exp(-beta*deltaE)
    r = rand.random()
    if (q>r):
        Sold[i] *=-1
        return Sold, deltaE
    else:
        return Sold, 0      


def E(Sold,T, N_samp):    
    beta = 1/T
    U = 0
    C = 0
    N_wait = N_samp
    n = Sold.shape[0]
    #n = 1
    e = E1d(Sold)
    for i in range(N_samp*n):
            a = energycrit(Sold,beta)
            Sold = a[0]
            e += a[1]
            U += e/(N_samp*n)
            C += e**2/(N_samp*n)  
    return Sold, [U,C]


"""
S10 = S1d(N[0])
ET10 = []

#E_vec = np.vectorize(E, excluded= ["Sold","N_samp"])
for i in range(len(T)):
    et10 = E(S10,T[-(i+1)],N_samples[0])
    ET10 += [et10[1]]
    S10 = et10[0]
"""   
 

S10 = S1d(N[0])
S100 = S1d(N[1])
S1000 = S1d(N[2])
E10_1000 = []
E100_1000 = []
E1000_1000 = []
S102 = S1d(N[0])
S1002 = S1d(N[1])
S10002 = S1d(N[2])
E10_10000 = []
E100_10000 = []
E1000_10000 = []
i = 0
#E_vec = np.vectorize(E, excluded= ["Sold","N_samp"])
for i in range(len(T)):
    print(T[-(i+1)])
    et10 = E(S10,T[-(i+1)],N_samples[0])
    et100 = E(S100,T[-(i+1)],N_samples[0])
    et1000 = E(S1000,T[-(i+1)],N_samples[0])
    E10_1000 += [et10[1]]
    E100_1000 += [et100[1]]
    E1000_1000 += [et1000[1]]
    S10 = et10[0]
    S100 = et100[0]
    S1000 = et1000[0]

    
k = 0   
for k in range(len(T)):
    print(T[-(k+1)])
    et10 = E(S102,T[-(k+1)],N_samples[1])
    et100 = E(S1002,T[-(k+1)],N_samples[1])
    et1000 = E(S10002,T[-(k+1)],N_samples[1])
    E10_10000 += [et10[1]]
    E100_10000 += [et100[1]]
    E1000_10000 += [et1000[1]]
    S102 = et10[0]
    S1002 = et100[0]
    S10002 = et1000[0]


#plt.plot(T,(np.array(ET10).T[0])[::-1]/10)
#plt.plot(T, -(N[0]-1)/N[0]*np.tanh(1/T), label = "Theory")



E10_1000 = np.array(E10_1000).T
E100_1000 = np.array(E100_1000).T
E1000_1000 = np.array(E1000_1000).T
E10_10000 = np.array(E10_10000).T
E100_10000 = np.array(E100_10000).T
E1000_10000 = np.array(E1000_10000).T


#np.savetxt('U10.csv', [E10_1000[0],E10_1000[1],E10_10000[0],E10_10000[1]], delimiter=';') 
#np.savetxt('U100.csv', [E100_1000[0],E100_1000[1],E100_10000[0],E100_10000[1]], delimiter=';') 
#np.savetxt('U1000.csv', [E1000_1000[0],E1000_1000[1],E1000_10000[0],E1000_10000[1]], delimiter=';')



end = time.time()
print(end-start)    