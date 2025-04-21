#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 24 13:30:53 2023

@author: tamilarasan
"""
import numpy as np
import matplotlib.pyplot as plt

import os 
os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin'
#some plotting arguments to make the plot prettier
plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "times new roman"
plt.rcParams["font.size"] = "40"
plt.rcParams['lines.linewidth'] = 4
plt.rcParams["font.weight"] = "bold"


Deltat = np.array([0.1,0.01]) # the delta t we use 

#function solving the equation for given delta t = dt
def Exercise3(dt):
    #function defining the acceleration
    def a(x):
        return -x
    #length of array depends on the delta t
    if dt == 0.1:
        M = 101
    elif dt == 0.01:
        M = 1001
    #initializing the arrays for position x and velocity v
    x = np.zeros(M)
    v = np.zeros(M)
    #initial condition
    v[0] =1
    
    #velocity verlet algorithm, xdt = x(t+dt) and vdt = v(t+dt)
    def Velocity_Verlet(x,v):
        xdt= x+ v*dt+ 1/2 *a(x)*dt**2
        vdt= v+1/2*(a(x) + a(xdt))*dt
        return xdt,vdt
    
    #iterating over all the indexes in the array
    for i in range(M-1):
        x[i+1],v[i+1] = Velocity_Verlet(x[i],v[i])
    
    
    #calculate the energy
    Et = v**2/2 + x**2/2
    #defining the time array
    t = np.linspace(0,(M-1)*dt,M)
    
    #below only code needed for the plots
    plt.figure(figsize= [16,9])
    plt.title(r"Velocity Verlet algorithm - $\Delta t = $ "+str(dt))
    plt.plot(t,x, label = r"Position - $x(t)$")
    plt.xlabel("time t", loc="right")
    plt.ylabel(r"$x(t)$, $v(t)$ and $E(t)$")
    plt.plot(t,v,label = r"Velocity - $v(t)$", color = "tab:green")
    plt.plot(t,np.sin(t), linestyle='dashed', label = r"Position - $x_{theory}(t)$", color ="black")
    plt.plot(t,Et, label = r"Energy - $E(t)$", color = "tab:red")
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.14),fancybox=True, shadow=True, ncol=2)
    plt.grid()
    plt.savefig("Velocity-Verlet_x_v_"+str(dt)+".pdf",bbox_inches='tight')
    plt.show()
    plt.clf()
    plt.figure(figsize= [16,9])
    plt.plot(t,(Et - 1/2), label = r"$E(t)$ - $E_{theory}$", color = "tab:red")
    plt.xlabel("time t", loc="right")
    plt.ylabel(r"$E(t) - E_{theory}$")
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.08),fancybox=True, shadow=True, ncol=3)
    plt.ticklabel_format(axis='y', scilimits=[-3, 3])
    plt.grid()
    plt.savefig("Velocity-Verlet_E_"+str(dt)+".pdf",bbox_inches='tight')
    plt.show()
    plt.clf()
    return 0    

#executing the task for both delta ts
Exercise3(Deltat[0])
Exercise3(Deltat[1])

