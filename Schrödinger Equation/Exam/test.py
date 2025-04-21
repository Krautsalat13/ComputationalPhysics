#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 16:47:22 2023

@author: tamilarasan
"""


from multiprocessing import Pool
import time
import math
import numpy as np
import matplotlib.pyplot as plt

N = 5000000

def cube(x):
    return math.sqrt(x)
"""
if __name__ == '__main__':
    pool = Pool()
    args = [10, 11]
    results = pool.map(cube, args)
"""   
    
a = [0,2,4,6,8,10]

for i in range(10):
    plt.plot(a,a)
