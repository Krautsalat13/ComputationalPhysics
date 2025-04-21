#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 14 18:57:39 2023

@author: tamilarasan
"""

import numpy as np
import matplotlib.pyplot as plt
import time
import csv
#plt.rcdefaults()

import os 
#os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin'

plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "times new roman"
plt.rcParams["font.size"] = "18"


T = np.linspace(0.2,4.2,21)
N = np.array([10,100,1000])


a, b, c, d= np.loadtxt('U10.csv', delimiter=';')

E10_1000 = np.array([a[::-1],b[::-1]])
E10_10000 = np.array([c[::-1],d[::-1]])


a, b, c, d = np.loadtxt('U100.csv', delimiter=';')

E100_1000 = np.array([a[::-1],b[::-1]])
E100_10000 = np.array([c[::-1],d[::-1]])


a, b, c, d= np.loadtxt('U1000.csv', delimiter=';')

E1000_1000 = np.array([a[::-1],b[::-1]])
E1000_10000 = np.array([c[::-1],d[::-1]])



#10
plt.plot(T, E10_1000[0]/10, label = r"10 Spins $N_{S} = 1000$")
plt.plot(T, E10_10000[0]/10, label = r"10 Spins $N_{S} = 10000$")

plt.plot(T, -(N[0]-1)/N[0] *np.tanh(1/T), label = "Theory")

plt.title("Average Energy")
plt.xlabel("Temperature")
plt.ylabel(r"U/N")
#plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig("U10.pdf", bbox_inches='tight')

plt.show()


plt.plot(T, 1/T**2*(E10_1000[1]- E10_1000[0]**2)/10, label = r"10 Spins $N_{S} = 1000$")
plt.plot(T, 1/T**2*(E10_10000[1]- E10_10000[0]**2)/10, label = r"10 Spins $N_{S} = 10000$")

plt.plot(T, (N[0]-1)/N[0] *1/(T*np.cosh(1/T))**2, label = "Theory")

plt.title("Specific Heat")
plt.xlabel("Temperature")
plt.ylabel(r"C/N")
#plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig("C10.pdf", bbox_inches='tight')

plt.show()

#100
plt.plot(T, E100_1000[0]/100, label = r"100 Spins $N_{S} = 1000$")
plt.plot(T, E100_10000[0]/100, label = r"100 Spins $N_{S} = 10000$")

plt.plot(T, -(N[1]-1)/N[1] *np.tanh(1/T), label = "Theory")

plt.title("Average Energy")
plt.xlabel("Temperature")
plt.ylabel(r"U/N")
#plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig("U100.pdf", bbox_inches='tight')

plt.show()
plt.plot(T, 1/T**2*(E100_1000[1]- E100_1000[0]**2)/100, label = r"100 Spins $N_{S} = 1000$")
plt.plot(T, 1/T**2*(E100_10000[1]- E100_10000[0]**2)/100, label = r"100 Spins $N_{S} = 10000$")

plt.plot(T, (N[1]-1)/N[1] *1/(T*np.cosh(1/T))**2, label = "Theory")

plt.title("Specific Heat")
plt.xlabel("Temperature")
plt.ylabel(r"C/N")
#plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig("C100.pdf", bbox_inches='tight')

plt.show()


#1000
plt.plot(T, E1000_1000[0]/1000, label = r"1000 Spins $N_{S} = 1000$")
plt.plot(T, E1000_10000[0]/1000, label = r"1000 Spins $N_{S} = 10000$")

plt.plot(T, -(N[2]-1)/N[2] *np.tanh(1/T), label = "Theory")

plt.title("Average Energy")
plt.xlabel("Temperature")
plt.ylabel(r"U/N")
#plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig("U1000.pdf", bbox_inches='tight')

plt.show()

plt.plot(T, 1/T**2*(E1000_1000[1]- E1000_1000[0]**2)/1000, label = r"10 Spins $N_{S} = 1000$")
plt.plot(T, 1/T**2*(E1000_10000[1]- E1000_10000[0]**2)/1000, label = r"10 Spins $N_{S} = 10000$")

plt.plot(T, (N[2]-1)/N[2] *1/(T*np.cosh(1/T))**2, label = "Theory")

plt.title("Specific Heat")
plt.xlabel("Temperature")
plt.ylabel(r"C/N")
#plt.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig("C1000.pdf", bbox_inches='tight')

plt.show()

