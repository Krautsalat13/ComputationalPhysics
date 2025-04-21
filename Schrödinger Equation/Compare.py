#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  1 15:38:07 2023

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

f5 = np.load("/Users/tamilarasan/Downloads/Data/TDSE_task2_times5.npy")
f10 = np.load("/Users/tamilarasan/Downloads/Data/TDSE_task2_times10.npy")
f25 = np.load("/Users/tamilarasan/Downloads/Data/TDSE_task2_times25.npy")
f40 = np.load("/Users/tamilarasan/Downloads/Data/TDSE_task2_times40.npy")
f45 = np.load("/Users/tamilarasan/Downloads/Data/TDSE_task2_times45.npy")
f50 = np.load("/Users/tamilarasan/Downloads/Data/TDSE_task2_times50.npy")



t5 = np.load("Data/TDSE_task1_times5.npy")
t10 = np.load("Data/TDSE_task1_times10.npy")
t25 = np.load("Data/TDSE_task1_times25.npy")
t40 = np.load("Data/TDSE_task1_times40.npy")
t45 = np.load("Data/TDSE_task1_times45.npy")
#t50 = np.load("Data/TDSE_task1_times50.npy")

def P(t):
    return np.absolute(t)**2

plt.plot(P(t5)-P(f5),label = "t = 5")
plt.plot(P(t10)-P(f10),label = "t = 10")
plt.plot(P(t25)-P(f25),label = "t = 25")
#plt.plot(P(t40)-P(f40),label = "t = 40")
#plt.plot(P(t45)-P(f45),label = "t = 45")
plt.ylabel("Tami-Fynn")
plt.legend()

"""
f5 = np.load("/Users/tamilarasan/Downloads/Data/TDSE_task1_times5.npy")
f10 = np.load("/Users/tamilarasan/Downloads/Data/TDSE_task1_times10.npy")
f40 = np.load("/Users/tamilarasan/Downloads/Data/TDSE_task1_times40.npy")
f45 = np.load("/Users/tamilarasan/Downloads/Data/TDSE_task1_times45.npy")




t5 = np.load("Data/TDSE_task0_times5.npy")
t10 = np.load("Data/TDSE_task0_times10.npy")
t40 = np.load("Data/TDSE_task0_times40.npy")
t45 = np.load("Data/TDSE_task0_times45.npy")

def P(t):
    return np.abs(t)**2

plt.plot(P(t5)-P(f5),label = "t = 5")
plt.plot(P(t10)-P(f10),label = "t = 10")
plt.plot(P(t40)-P(f40),label = "t = 40")
plt.plot(P(t45)-P(f45),label = "t = 45")
plt.ylabel("Tami-Fynn")
plt.legend()
"""

"""
twithout = np.load("Psum_without.npy")
twith = np.load("Psum_with.npy")

fwithout = np.load("/Users/tamilarasan/Downloads/P_max_ohneV.npy")
fwith = np.load("/Users/tamilarasan/Downloads/P_max_mitV.npy")

plt.plot(twith-fwith)
plt.plot(twithout-fwithout)
"""

"""
twith = np.load("Ptot_with.npy")
fwith = np.load("/Users/tamilarasan/Downloads/P_tot.npy")

twithout = np.load("Ptot_without.npy")
fwithout = np.load("/Users/tamilarasan/Downloads/P_tot_ohneV.npy")

#plt.plot(twith-fwith)
plt.plot(twithout-fwithout)
"""