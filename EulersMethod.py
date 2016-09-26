# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 12:40:21 2016

@author: danielko
"""

import math as m
import matplotlib.pyplot as p

##EULER's METHOD

def EuMethod(f, a, b, n, alpha):
    #f is the function, a & b the bounds, and alpha the initial condition
    h = (b-a)/n
    tlist = []
    for i in range(0, n):
        tlist.append(a+(i*h))
    wlist = []
    wlist.append(alpha) #so the w0 is the initial condition
    for j in range(0, n-1):
        wvalue = wlist[j]+h*f(tlist[j],wlist[j])
        wlist.append(wvalue)
    return wlist



'''
a = 0
b = ending time in days
n = time period in which you want values specified
alpha = infectious rate
beta = recovery rate 
h = steps in time
p = total population
tot = totalpopulation = s+i+r
'''


def EuSIR(a, b, n, S, I, R):
    h = (b-a)/float(n)
    y = int(n)
    t = []  #time
    Sus = []  #Susceptible population
    Inf = []  #Infectious population
    Rec = []  #Recovery population
    t.append(a)
    Sus.append(S)
    Inf.append(I)
    Rec.append(R)
    for i in range(1, y+1):
        tvalue = a + i*h
        Svalue = Sus[i-1] - (2.0*Sus[i-1]*Inf[i-1]*h)
        Rvalue = Rec[i-1] + (0.5*Inf[i-1]*h)
        Ivalue = Inf[i-1] + (2.0*Sus[i-1]*Inf[i-1] - 0.5*Inf[i-1]*h)
        t.append(tvalue)
        Sus.append(Svalue)
        Inf.append(Ivalue)
        Rec.append(Rvalue)
    return t, Sus, Inf, Rec

T,_,_,_ = EuSIR(0, 15.0, 100.0, 0.99, 0.001, 0)
_,S,_,_ = EuSIR(0, 15.0, 100.0, 0.99, 0.001, 0)
_,_,I,_ = EuSIR(0, 15.0, 100.0, 0.99, 0.001, 0)
_,_,_,R = EuSIR(0, 15.0, 100.0, 0.99, 0.001, 0)


p.plot(T, S, label = "Susceptible")
p.plot(T, I, label = "Infectious")
p.plot(T, R, label = "Recovered")
p.xlabel("Time in Days")
p.ylabel("Population Proportion")
p.title("SIR Model through Euler's Method")
p.legend(loc = 7)
p.show()

        
    