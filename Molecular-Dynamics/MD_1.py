#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 24 10:13:53 2023

@author: tamilarasan
"""
import numpy as np
import matplotlib.pyplot as plt
import csv
import sys
sys.setrecursionlimit(20000)
import os 
os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin'

plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "times new roman"
plt.rcParams["font.size"] = "18"


def a(x):
    return -x

r0 = 0
v0 = 1
a0 = a(r0)
Deltat = np.array([0.1,0.01,0.001])


def solve(r,v,i,deltat):
    r+= [r[i]+ v[i]*deltat]
    v+= [v[i]+a(r[i])*deltat]
    if i == 10000:
        return r
    return solve(r,v,i+1,deltat)

r = [r0]
v = [v0]
rend = []
t = []
for j in range(3):
    rend += [solve(r,v,0,Deltat[j])]    
    t += [Deltat[j]*np.arange(len(rend[j]))]
    plt.plot(rend[0])

#plt.plot(t[0],np.sin(t[0]))

