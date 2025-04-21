#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 24 11:35:43 2023

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

def solvea(r,v,i,deltat):
    
    v+= [v[i]+a(r[i])*deltat]
    r+= [r[i]+ v[i+1]*deltat]
    
    if i == 10000:
        return r,v
    return solvea(r,v,i+1,deltat)

r = [r0]
v = [v0]

renda, venda = np.array(solvea(r,v,0,0.01))

#plt.plot(venda[:1000])
#plt.plot(renda[:1000])
#plt.plot(venda[:1000]**2/2 + renda[:1000]**2/2)
#plt.show()

def solveb(r,v,i,deltat):
    
    r+= [r[i]+ v[i]*deltat]
    v+= [v[i]+a(r[i+1])*deltat]
    
    
    if i == 10000:
        return r,v
    return solveb(r,v,i+1,deltat)

r = [r0]
v = [v0]

rendb, vendb = np.array(solveb(r,v,0,0.01))
t = np.linspace(0,100,10000)
#plt.plot(vendb[:1000])
#plt.plot(rendb[:1000])
plt.plot(t,vendb[:10000]**2/2 + rendb[:10000]**2/2-0.5 - 0.01**2/2*())
#plt.plot(t,-1/4*0.01*np.sin(2*t))
plt.show()