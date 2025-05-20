#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 24 13:54:05 2023

@author: tamilarasan
"""
import numpy as np
import matplotlib.pyplot as plt
import csv


import os 
os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin'

plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "times new roman"
plt.rcParams["font.size"] = "40"
plt.rcParams['lines.linewidth'] = 4
plt.rcParams["font.weight"] = "bold"


pi      = np.pi  

M = [100,1000]      #number of steps
N = [4,16,128]      #number of oscillators
dt = [0.1, 0.01]    #delta t


def Exercise4(N,dt,M,config,plot, z ="selected"):
    #depending on the initial condition the arrays are initialized differently
    def initial(config = "a"):
        x = np.zeros((N,M))
        v = np.zeros((N,M))
        if config == "a":           #first condition
            x[N//2 -1,0] = 1 
        elif config == "b":         #second condtion with j = 1
            j = 1
            x[:,0] = np.sin(pi*j*(np.arange(N)+1)/(N+1)) 
        elif config =="c":          #second condtion with j = N/2
            j = N/2
            x[:,0] = np.sin(pi*N/2*(np.arange(N)+1)/(N+1))
        return x,v
            
    x,v = initial(str(config))  #initialize arrays for configuration
    
    #function that yields the forces acting oscillator depending on oscillator
    def Fn(n,i):
        if n == 0:
            return -(x[0,i]-x[1,i])
        elif n == (N-1):
            return -(x[N-1,i]-x[N-2,i])
        else:
            return -(2*x[n,i] - x[n-1,i] - x[n+1,i])
    
    #algorithm that calculates position and velocity
    def Velocity_Verlet(x,v,i,N):
        #iterate over all the oscillators 
        for k in range(N):
            x[k,i+1]= x[k,i]+dt*(v[k,i]+dt/2*Fn(k,i))   #x[k,i+1] = x_k (t+dt)
        for k in range(N):
            v[k,i+1]= v[k,i]+dt/2*(Fn(k,i) + Fn(k,i+1)) #v[v,i+1] = v_k (t+dt)
        return x,v
    
    #call function iteratively
    for i in range(M-1):
        x,v = Velocity_Verlet(x,v,i,N)
    #energy of system
    def E(x, v):
        KE = 1/2*np.sum(v**2, axis=0)
        PE = 1/2* np.sum(np.diff(x, axis = 0)**2, axis = 0)
        return KE+PE
    
    #calculates velocity and position for the arguments given to the function
    X,V = Velocity_Verlet(x,v,0,N)
    #time array
    t = np.linspace(0,dt*(M-1),M)
    
    #code that chooses which functions are plotted 
    if plot =="x":
        plt.figure(figsize= [16,9])
        plt.grid()
        if z == "all":
            for k in range(N):
                plt.plot(t,X[k])
            
        elif N ==4 and config != "c":
            plt.plot(t,X[0], label = r"$x_{1}(t)$")
            plt.plot(t,X[1], label = r"$x_{2}(t)$")
            plt.plot(t,X[2], label = r"$x_{3}(t)$", ls = "--")
            plt.plot(t,X[3], label = r"$x_{4}(t)$",  ls = "--")
        elif config == "a":
            if N!=4:
                plt.plot(t, X[0], label =r"$x_{1}(t)$")
                plt.plot(t, X[N//2 -2], label =r"$x_{N/2-1}(t)$", ls = "-.")
                plt.plot(t, X[N//2 -1], label =r"$x_{N/2}(t)$") #oscillator which is extendend
                plt.plot(t, X[N//2 ], label =r"$x_{N/2+1}(t)$", ls = "--")
                plt.plot(t, X[-1], label =r"$x_{N}(t)$")
                
        elif config == "b":
            if N!=4:
                plt.plot(t, X[0], label =r"$x_{1}(t)$", ls = "-")
                plt.plot(t, X[1], label =r"$x_{2}(t)$", ls = "-")
                plt.plot(t, X[N//2 -2], label =r"$x_{N/2-1}(t)$", ls = "-")
                plt.plot(t, X[N//2 -1], label =r"$x_{N/2}(t)$")
                plt.plot(t, X[N//2 ], label =r"$x_{N/2+1}(t)$", ls = "--")
                plt.plot(t, X[N//2 +1], label =r"$x_{N/2+2}(t)$", ls = "--")
                plt.plot(t, X[-2], label =r"$x_{N-1}(t)$",ls = "--")
                plt.plot(t, X[-1], label =r"$x_{N}(t)$",ls = "--")
        elif config == "c":
            if N ==4:
                plt.plot(t,X[0], label = r"$x_{1}(t)$")
                plt.plot(t,X[1], label = r"$x_{2}(t)$")
                plt.plot(t,X[2], label = r"$x_{3}(t)$")
                plt.plot(t,X[3], label = r"$x_{4}(t)$")
            elif N!=4:
                plt.plot(t, X[0], label =r"$x_{1}(t)$")
                plt.plot(t, X[1], label =r"$x_{2}(t)$")
                plt.plot(t, X[N//2 -2], label =r"$x_{N/2-1}(t)$")
                plt.plot(t, X[N//2 -1], label =r"$x_{N/2}(t)$") 
                plt.plot(t, X[N//2 ], label =r"$x_{N/2+1}(t)$")
                plt.plot(t, X[N//2 +1], label =r"$x_{N/2+2}(t)$")
                plt.plot(t, X[-2], label =r"$x_{N-1}(t)$")
                plt.plot(t, X[-1], label =r"$x_{N}(t)$")
        plt.title(r"$N =$"+str(N)+", $\Delta t = $ "+str(dt))
        plt.xlabel("time t", loc="right")

        plt.ylabel(r"$x(t)$")
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.tight_layout()
        plt.savefig("images/"+str(N)+"_"+str(dt)+"_"+config+".pdf",bbox_inches='tight')
        plt.show()
    elif plot =="E":
        plt.figure(figsize= [16,9])
        plt.grid()
        plt.plot(t, E(X,V))
        plt.xlabel("time t")
        plt.ylabel(r"$E(t)$")
        plt.title(r" $N = $ "+str(N)+", $\Delta t = $ "+str(dt), y=1.1)
        plt.tight_layout()
        plt.savefig("images/"+str(N)+"_E_"+str(dt)+"_"+config+".pdf",bbox_inches='tight')
        plt.show()
 
#vectorize function to call function with arrays        
Exercise4vec = np.vectorize(Exercise4, excluded= "dt,M,N")


#calclulate the function for each parameter
for n in N:
    for it in range(len(dt)):
        #execute code with M+1 because we do M steps
        Exercise4vec(n,dt[it],M[it]+1, ["a","b", "c"],plot = "E", z = "all")
        Exercise4vec(n,dt[it],M[it]+1, ["a","b", "c"],plot = "x") 
        Exercise4vec(n,dt[it],M[it]+1, ["a","b", "c"],plot = "E")
        
