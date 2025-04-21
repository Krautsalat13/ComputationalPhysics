#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 18 16:23:33 2023

@author: fynn@tamilarasan
"""
import numpy as np
import matplotlib.pyplot as plt
import csv
from numba import njit, vectorize, float64
import os 
os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin'

plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "times new roman"
plt.rcParams["font.size"] = "27"


pi      = np.pi  
rand    = np.random
#set the seed
rand.seed(1069)
#the constants given by the exercise
L = 1001
n = 10000
Delta = 0.1
D = 5  
#inital conditions
i0  = (L+1)//2
#array consisting the current position of each particle
N = np.ones(n)*i0
tau = 0.001
#time array
t = np.arange(11)


#a function returning the direction each particle will be moved in
def walk(n):
    return rand.choice([-1,1],n)

#function that does one iteration of t = t+1
def onet(N):
    F = N.copy()
    for i in range(int(1/tau)):
        F = F+walk(len(F))
    return F

#array below will contain the position of the particles at each time
N_total = [N]
#loops over all times
for T in t:
    N_total +=[onet(N_total[T])]

#function calculating the averages of x^p    
def x_mean_p(p, N):
    i_arr   = np.arange(1, L+1)
    return np.sum((i_arr-i0)**int(p) *N)/np.sum(N)

# Calculate the variance/Delta²
def var(N):
    return x_mean_p(2, N) - x_mean_p(1, N)

#array consisting the variances at each time t
Var = []
for i in range(11):
    #histogram is used to calculate the distribution of the particles at each position
    Var +=[var(np.histogram(N_total[i], bins = 1001, range=(0,1001))[0])]

#plots
plt.figure(figsize =(16 , 9))
plt.plot(t,Var, lw=5, marker="o", markersize=7, markerfacecolor="white", label= "result random walk simulation")
#theoretical prediction
plt.plot (t , 2*D/(Delta **2) * t , lw =3 , color = "black" , ls = "--", label= "theory diffusion equation" )
plt.title("Random Walk")
plt.xlabel(r"time $t$")
plt.grid()
plt.legend()
plt.ylabel(r"Var$(x)$")
plt.savefig("task2.pdf")