#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 18 13:56:14 2023

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
rand    = np.random

rand.seed(1069)

lamb = 1
grid = 50
Delta = lamb/grid

tau1 = 0.9*Delta
tau2 = 1.05*Delta

def t(n ,tau = tau1):
    return (n+1/2)*tau1