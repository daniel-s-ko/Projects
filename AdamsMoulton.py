# -*- coding: utf-8 -*-
"""
Created on Mon May  9 21:39:44 2016

@author: danielko
"""

import math as m

## ADAMS MOULTON

def AdMoul4(f, a, b, h, z):
    N = (b-a)/h
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
        wreal = w2[i-1] + (h/720)*(251*f(t2[i+1],w2[i+1])+646*f(t2[i],w2[i])-264*f(t2[i-1],w2[i-1])+106*f(t2[i-2],w2[i-2])-19*f(t2[i-3],w2[i-3]))
        wnode.append(wreal)
        tnode.append(t)
    return wnode, tnode

def AdMoul3(f, a, b, h, z):
    N = (b-a)/h
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
        wreal = w2[i-1] + (h/24.0)*(9*f(t2[i+1],w2[i+1])+19*f(t2[i],w2[i])-5*f(t2[i-1],w2[i-1])+f(t2[i-2],w2[i-1]))
        wnode.append(wreal)
        tnode.append(t)
    return wnode, tnode

def AdMoul2(f, a, b, h, z):
    N = (b-a)/h
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
        tnode.append(t)
        wreal = w2[i-1] + (h/12)*(5*f(t2[i+1],w2[i+1])+8*f(t2[i],w2[i])-f(t2[i-1],w2[i-1]))
        wnode.append(wreal)
    return wnode, tnode
    
## PROBLEM 4
#a
f = lambda t, w: t*(m.e**(3.0*t))-2.0*w
print AdMoul2(f, 0, 1, 0.2, 0)
print AdMoul3(f, 0, 1, 0.2, 0)
print AdMoul4(f, 0, 1, 0.2, 0)

'''
#c
f = lambda t, w: 1 + w/t
print AdMoul2(f, 1, 2, 0.2, 2)
print AdMoul3(f, 1, 2, 0.2, 2)
print AdMoul4(f, 1, 2, 0.2, 2)


#d
f = lambda t, w: m.cos(2.0*t)+m.sin(2.0*t)
print AdMoul2(f, 0, 1, 0.2, 1)
print AdMoul3(f, 0, 1, 0.2, 1)
print AdMoul4(f, 0, 1, 0.2, 1)

'''