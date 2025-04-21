#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 18 16:23:33 2023

@author: tamilarasan
"""
import numpy as np
import matplotlib.pyplot as plt
import csv
from numba import njit, vectorize, float64
import os 
os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin'

plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "times new roman"
plt.rcParams["font.size"] = "18"

pi      = np.pi  
rand    = np.random

rand.seed(1069)

N = 1
L = 1001
D = 1
Delta = 0.1
m = 10000
start  = (L+1)//2+1
alpha = -D*Delta**(-2)
tau  = 0.001
t = 10


X = [[1,-1],[-1,1]]


eA2 = 1/2 * np.array([[1+np.exp(alpha*t/(tau*m)),1-np.exp(alpha*t/(tau*m))],[1-np.exp(alpha*t/(tau*m)),1+np.exp(alpha*t/(tau*m))]])

eB = 1/2 * np.array([[1+np.exp(2*alpha*t/(tau*m)),1-np.exp(2*alpha*t/(tau*m))],[1-np.exp(2*alpha*t/(tau*m)),1+np.exp(2*alpha*t/(tau*m))]])

c = np.exp(tau*alpha/2)
Phi = np.zeros(L)
Phi[start] = N
x = np.arange(len(Phi))
plt.plot(x,Phi)
#@njit
def dt(Phi):
    for i in range(m):
        for k in range(0,len(Phi)-1,2):
            if (k ==0):
                Phi[-1] = c*Phi[-1]
            Phi[k:k+2] = np.dot(eA2,Phi[k:k+2])
        for k in range(0,len(Phi)-1,2):
            if (k ==0):
                Phi[0] = c*Phi[0]
            Phi[k+1:k+3] = np.dot(eB,Phi[k+1:k+3])
        for k in range(0,len(Phi)-1,2):
            Phi[k:k+2] = np.dot(eA2,Phi[k:k+2])
            if (k ==0):
                Phi[-1] = c*Phi[-1]
    plt.plot()
    return Phi

dt(Phi)        

