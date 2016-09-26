# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 17:15:47 2016

@author: danielko
"""

import numpy as np
import matplotlib.pyplot as p

def Euler(a, b, n, S, I, R):
    h = (b-a)/n
    m = int(n)
    t = np.array([])
    sp = np.array([])
    ip = np.array([])
    rp = np.array([])
    sp = np.append(sp, S)
    ip = np.append(ip, I)
    rp = np.append(rp, R)
    t = np.append(t, a)
    for i in range(1, m + 1):
        t = np.append(t, a + i*h)
        sp = np.append(sp, sp[i-1] - (2.0*sp[i-1]*ip[i-1]*h))
        rp = np.append(rp, rp[i-1] + (0.5*ip[i-1]*h))
        ip = np.append(ip, ip[i-1] + (2.0*sp[i-1]*ip[i-1] - 0.5*ip[i-1]*h))
    p.plot(t, ip, label = "infected")
    p.plot(t, sp, label = "susceptible")
    p.plot(t, rp, label = "recovered")
    p.xlabel("Time in Days")
    p.ylabel("Population Proportion")
    p.title("SIR Model through Euler's Method")
    p.legend(loc = 7)
    p.show()
    return t

Euler(0, 15, 100.0, 0.99, 0.001, 0)