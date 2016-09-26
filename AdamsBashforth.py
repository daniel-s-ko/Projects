# -*- coding: utf-8 -*-
"""
Created on Wed May  4 21:05:14 2016

@author: danielko
"""
import math as m

## HOMEWORK PROBLEM 1
## ADAMS - BASHFORTH

def AdBash(f, a, b, h, z):
    #z is alpha, the initial condition
    #f is a function of two variables, w and t
    N = (b-a)/h
    end = int(N)
    tnode = [a]
    wnode = [z]
    for i in range(1, 4):
        Kone = h*f(tnode[i-1],wnode[i-1])
        Ktwo = h*f(tnode[i-1]+(h/2.0), wnode[i-1]+(Kone/2.0))
        Kthree = h*f(tnode[i-1]+(h/2.0), wnode[i-1]+(Ktwo/2.0))
        Kfour = h*f(tnode[i-1]+(h), wnode[i-1]+(Kthree))
        
        wvalue = wnode[i-1] + (Kone+(2.0*Ktwo)+(2.0*Kthree)+Kfour)/6.0
        wnode.append(wvalue)
        tvalue = a + i*h
        tnode.append(tvalue)
    t2 = tnode
    w2 = wnode

    for i in range(4, end+1):
        t = a+(i*h)
        wpred = w2[i-1]+h*(55.0*f(t2[i-1],w2[i-1])-59.0*f(t2[i-2],w2[i-2])+37.0*f(t2[i-3],w2[i-3])-9.0*f(t2[i-4],w2[i-4]))/24.0
        wreal = w2[i-1]+h*(9.0*f(t,wpred)+19.0*f(t2[i-1],w2[i-1])-5*f(t2[i-2],w2[i-2])+f(t2[i-3],w2[i-3]))/24.0
        wnode.append(wreal)
        tnode.append(t)
      
    return wnode, tnode

'''
#PROBLEM 1
#a
f = lambda t, w: t*(m.e**(3.0*t))-2.0*w
print AdBash(f, 0, 1, 0.2, 0)

#c
f = lambda t, w: 1 + w/t
print AdBash(f, 1, 2, 0.2, 2)

#d
f = lambda t, w: m.cos(2.0*t)+m.sin(2.0*t)
print AdBash(f, 0, 1, 0.2, 1)

#PROBLEM 2



#b
f = lambda t, w: (w**2)/(1+t)
print AdBash(f, 1, 2, 0.1, -1.0/m.log(2))
'''
#c 
f = lambda t, w: ((w**2)+w)/t
print AdBash(f, 1, 3, 0.2, -2.0)