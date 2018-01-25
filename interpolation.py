# -*- coding: utf-8 -*-
" Sometimes we need to fit a function for perdiction or other usage from several given data. In order to achieve the objective, we use technique called interpolation or fitting. This pyfile can help us calculate the function required."

import numpy as np
import sympy as sp

""" Data Input: We input data in an array which has the form of (x,y)

data = np.array([[1,-1],
                 [2,-1],
                 [3,1]
                 ])
"""

" Polynomial Interpolation f(x) = c0 + c1x + c2x^2 + ... + cnx^n "


# In fact, this is a problem solving system of linear equation. Therefore, we have to meet the case that there are multiple solutions or probably no solution. To simplify our analysis, we only consider the case that the number of equation is equal to the number of variable. For example, if our data contains three (x,y), we can just solve f(x)=c0+c1x+c2x^2.

def Polynomial(data):
    dim = data.shape[0]
    A = np.zeros((dim,dim))
    for c in range(dim):
        g = lambda x:x**c
        A[:,c] = g(data[:,0])
    soln = np.linalg.solve(A,data[:,1])
    #for index in range(dim):
    #    print("c%d = %r"%(index,soln[index]))
    return soln

def Polynomial_predict(data,x):
    coefficient = Polynomial(data)
    result = 0
    for index in range(len(coefficient)):
        g = lambda x: x**index
        result = result + g(x)*coefficient[index]
    return result
    
# Polynomial_predict(data,-1)
data = np.array([[1,0],
                 [2,4],
                 [3,8]])
Polynomial(data)
Polynomial_predict(data,1.5)
" Lagrangian Interpolation "

def li_calc(data,i):
    dim = data.shape[0]
    x = sp.Symbol('x')
    li = 1
    for index in range(dim):
        if index != i:
            li = li*(x-data[index,0])/(data[i,0]-data[index,0])
        else:
            li = li
    return li

def Lagrangian(data):
    dim = data.shape[0]
    L = 0
    for index in range(dim):
        L = L +li_calc(data,index)*data[index,1]
    return L.expand()

def Lagrangian_predict(data,y):
    result = Lagrangian(data)
    x = sp.Symbol('x')
    return result.evalf(subs = {x:y})

data = np.array([[1,-1],
                 [2,-1],
                 [3,1]])
Lagrangian(data)
Lagrangian_predict(data,1.5)
    
" Newton Interpolation "

def Divided_difference(data):
    dim = data.shape[0]
    D = np.zeros((dim,dim+1))
    D[:,0] = data[:,0]
    D[:,1] = data[:,1]
    diff = 1
    for column in range(2,dim+1):
        for row in range(diff,dim):
            D[row,column] = (D[row,column-1]-D[row-1,column-1])/(D[row,0]-D[row-diff,0])
        diff = diff+1
    return D
    
# Divided_difference(data) 

def Newton(data):
    Coefficient = Divided_difference(data)
    x = sp.Symbol('x')
    Newton = 0
    X = 1
    dim = data.shape[0]
    for index in range(dim):
        Newton = Newton + Coefficient[index,index+1]*X
        X = X*(x-Coefficient[index,0])
    return Newton.expand()

"""
data = np.array([[0,2],
                [1,-3],
                [2,-6],
                [3,11]])
Newton(data)
"""

def Newton_predict(data,y):
    result = Newton(data)
    x = sp.Symbol('x')
    return result.evalf(subs = {x:y})

Divided_difference(data)
Newton(data)
Newton_predict(data,1.5)