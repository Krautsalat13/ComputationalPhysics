#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 14 16:35:58 2023

@author: tamilarasan
"""

import numpy as np
import matplotlib.pyplot as plt
import time as time
import scipy as sp

import os 
os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin'

plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "times new roman"
plt.rcParams["font.size"] = "18"
start = time.time()
pi      = np.pi  
rand    = np.random

rand.seed(1556)


N = np.array([10,50,100])
T = np.linspace(0.2,4.2,21)
T[-1] = 4

N_samples = np.array([1000,10000])

def S2d(n):
    return rand.choice([-1,1],(n,n))

def E2d(S,n):
    E = 0
    for i in range(int(n-1)):
        for j in range(int(n-1)):
            E+= S[i,j]*(S[i+1,j] + S[i,j+1])
    last = n-1
    E += S[last,last]*(S[last-1,last]+S[last,last-1])
    return -E  



def energycrit(S,N,beta, z = "periodic"):
    N = S.shape[0]
    n = N
    i,j = int(rand.random()*N),int(rand.random()*N)
    if (z == "free"):
        first_term = S[i-1,j] if i > 0 else 0
        second_term = S[i+1,j] if i < N-1 else 0 
        third_term = S[i,j-1] if j>0 else 0 
        fourth_term = S[i,j+1] if j < N-1 else 0  
    else:
        first_term = S[i-1,j]
        second_term = S[i+1,j] if i < N-1 else S[0,j] 
        third_term = S[i,j-1]
        fourth_term = S[i,j+1] if j < N-1 else S[i,0] 
        #deltaE = 2*S[i,j]*(S[i-1,j]+S[i+1,j]+S[i,j-1]+S[i,j+1])
        
    deltaE = 2*S[i,j]*(first_term+second_term+third_term+fourth_term)
    q = np.exp(-beta*deltaE)
    r = rand.random()
    if (q>r):
        S[i,j] *=-1
        return S, deltaE, 2*S[i,j]
    else:
        return S, 0, 0

def E(Sold, T, N_samples):
    #plt.imshow(Sold)
    #plt.show()
    N = Sold.shape[0]
    beta = 1/T
    e = E2d(Sold,N)
    U = 0
    C = 0
    M = 0
    Mold = np.sum(Sold)
    for i in range(N_samples):
        for j in range(N**2):
            a = energycrit(Sold,N,beta)
            Sold = a[0]
            e += a[1]
            U += e/(N_samples*N*N)
            C += e**2/(N_samples*N*N)
            Mold = Mold +a[2]
            M += Mold
    return Sold,U,C,M/(N_samples*N*N)


S10_1000 = S2d(N[0])
#S10_10000 = S2d(N[1])
#S100 = S2d(N[1])
#S1000 = S2d(N[2])
M10_1000 = []
U10_1000 = []
C10_1000 = []
"""
M10_10000 = []
U10_10000 = []
C10_10000 = []
"""
for i in range(len(T)):
    print(T[-(i+1)])
    Et1000 = E(S10_1000,T[-(i+1)], N_samples[0])
    S10_1000 = Et1000[0]
    U10_1000 += [Et1000[1]]
    C10_1000 += [Et1000[2]]
    M10_1000 += [np.abs(Et1000[3])]
    """
    Et10000 = E(S10_10000,T[-(i+1)], N_samples[1])
    S10_10000 = Et10000[0]
    U10_10000 += [Et10000[1]]
    C10_10000 += [Et10000[2]]
    M10_10000 += [np.abs(Et10000[3])]
    """

    
np.savetxt('2d10free_1000.csv', [np.array(M10_1000),np.array(U10_1000),np.array(C10_1000)], delimiter=';')

#np.savetxt('2d100free_10000.csv', [np.array(M10_10000),np.array(U10_10000),np.array(C10_10000)], delimiter=';')

def M2d(T):
    Tc = 2/(np.log(1+np.sqrt(2)))
    if T < Tc:
        return (1-(np.sinh(2*1/T))**(-4))**(1/8)
    else:
        return 0
    
    
M2dvec = np.vectorize(M2d)

#plt.plot(T,np.array(M10_1000)[::-1]/N[0]**2)    
#plt.plot(T,M2dvec(T))    

plt.show()



T = np.linspace(0.2,4.0,20)
def k(T):
    return 1/np.sinh(2/T)**2

def integrand(k,x):
    return 1/np.sqrt(1 - 4*k/(1+k)**2 * np.sin(x)**2)


def U_(T):
    I = sp.integrate.quad(lambda x: integrand(k(T), x), 0, np.pi/2)[0]
    U = -1/np.tanh(2/T)*( 1 + 2/np.pi *(2*np.tanh(2/T)**2 -1 )*I )
    return U

U = np.vectorize(U_)

plt.plot(T,np.array(U10_1000)[::-1][:-1]/N[0]**2)    
plt.plot(T,U(T))   

end = time.time()
print(end-start)
    