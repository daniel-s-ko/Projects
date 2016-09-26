# -*- coding: utf-8 -*-
"""
Created on Mon May 23 18:42:23 2016

@author: danielko
"""

import math as m
import numpy as np
import matplotlib.pyplot as p

## TRI DIAGONAL SOLVER
def GaussianTri(A, B):
    n = len(A)
    #FORWARD ELIMINATION
    for r in range(0, n):
        for i in range(r+1, n):
            multiplier = A[i][r]/A[r][r]
            for j in range(r+1, r+2):
                A[i][j] = A[i][j]-(multiplier*A[r][j])
            B[i] = B[i]-(multiplier*B[r])
    
    #MAKING THE X VECTOR WITH 0's
    X = np.zeros((n,1))
    
    #BACKWARDS SUBSTITUTION
    X[n-1] = B[n-1]/A[n-1][n-1]
    for j in range(n-1, -1, -1):
        X[j]=B[j]
        for k in range(j+1, n):
            X[j]=X[j]-(A[j][k]*X[k])
        X[j]=X[j]/A[j][j]    
    return X
    




def GaussTrin(n, A):
    X =[0]*n
    for i in range(0 , n-1):
        for j in range(i+1, n+1):
            mult = (A[j][i])/(A[i][i])
            A[j] = A[j] - mult*A[i]
    X[n-1] = (A[n-1][n])/(A[n-1][n-1])
    for c in range(n-2, 0, -1):
        ## need to do the summation in the formula
        sumval = 0
        for k in range(c+1, n):
            sumval += A[c][k]*X[k]
        X[c] = (A[c][n] - sumval)/A[c][c]
    return X

'''
## B = vector of answers, non augmented matrix
def GaussTri(n, A, B):   ## n = no. of equations, A =  matrix using np.matrix([[]])
    X = [0]*n # answer vector
    for r in range(0, n):
        for i in range(r+1, n):
            multiplier = A[i,r]/A[r,r]
            for j in range(r+1, r+2):
                A[i,j] = A[i,j]-(multiplier*A[r, j])
            X[i] = X[i]-(multiplier*X[r])
    
 '''
  