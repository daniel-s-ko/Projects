# -*- coding: utf-8 -*-
"""
Created on Wed May  4 23:48:31 2016

@author: danielko
"""

## RUNGE KUTTA METHOD FOR SYSTEMS

# take in a list as number of initial conditions
# incond starts at alpha1, so alpha0 = 0

# for problems 1 and 2
# only 2 initial conditions
# f is a LIST of functions
def RKS2(f, a, b, m, h, incond):
    n = (b-a)/h
    end = int(n)
    t = a
    w1 = [0]
    w2 = [0]
    #values for w start at w1, reiterate by one spot, so w0= 0
    for i in range(1, m+1):
        w1.append(incond[i])
        w2.append(incond[i])
    
    
    k1 = [0]
    k2 = [0]
    k3 = [0]
    k4 = [0]
    for c in range(1, end+1):   
        for j in range(1, m+1):
            k1value = h*f[j-1](t, w1[1], w2[2])
            k1.append(k1value)
        for k in range(1, m+1):
            k2value = h*f[k-1](t+(h/2.0),w1[1]+k1[1]/2.0, w2[2]+k1[2]/2.0)
            k2.append(k2value)
        for l in range(1, m+1):
            k3value = h*f[l-1](t+(h/2.0),w1[1]+k2[1]/2.0, w2[2]+k2[2]/2.0)
            k3.append(k3value)
        for p in range(1, m+1):
            k4value = h*f[p-1](t+(h),w1[1]+k3[1], w2[2]+k3[2])
            k4.append(k4value)
        for q in range(1, m+1):
            w1[q] = w1[q-1]+(1/6.0)*(k1[1]+2.0*k2[1]+2.0*k3[1]+k4[1])
            w2[q] = w2[q-1]+(1/6.0)*(k1[2]+2.0*k2[2]+2.0*k3[2]+k4[2])
        t = a + c*h
    return t, w1, w2
   
    
## EXAMPLE PROBLEM
f = [lambda i, j, t: -4.0*i+3.0*j-6.0, lambda k, l, t: -2.4*k+1.6*l+3.6]    
initialconditions = [0, 0, 0, 0, 0]

print RKS2(f, 0, 1, 2, 0.1, initialconditions)