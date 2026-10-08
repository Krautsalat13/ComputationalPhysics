#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Plot 2D Ising observables from the CSV results in ``Data/2D/``.

Reads the magnetisation / energy / heat-capacity CSVs written by ``2D.py`` and
plots them against temperature. ``IM.py`` (in the project root) is the main 2D
plotting script; this one plots the sizes separately.

Run from the project root:  ``python 2D/plot2.py``
"""

import numpy as np
import matplotlib.pyplot as plt
import time
import csv

import os 
#os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin'

plt.rcParams["text.usetex"] = False
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = "18"


T = np.linspace(0.2,4.2,21)
N = np.array([10,50,100])
N = N[2]

a, b, c= np.loadtxt('Data/2D/2d100_1000.csv', delimiter=';') #MUC

U1000_1000 = np.array(b[::-1])
M1000_1000 = np.array(a[::-1])
C1000_1000 = np.array(c[::-1])

def M2d(T):
    Tc = 2/(np.log(1+np.sqrt(2)))
    if T < Tc:
        return (1-(np.sinh(2*1/T))**(-4))**(1/8)
    else:
        return 0
    
M2dvec = np.vectorize(M2d)

plt.plot(T,M1000_1000/N**2)    
plt.plot(T,M2dvec(T))  


a10, b10, c10= np.loadtxt('Data/2D/2d10_1000.csv', delimiter=';') #MUC
a, b, c= np.loadtxt('Data/2D/2d10_10000.csv', delimiter=';') #MUC