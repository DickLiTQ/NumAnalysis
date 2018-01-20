# -*- coding: utf-8 -*-

"""
This file is designed to find numerical solution the non-linear equation through recursive method and Newton method.

Bibliograhpy: Numerical Method(计算方法) Damei Li(李大美), Wuhan University.
"""
import numpy as np

""" Recursive Method 
 f(x) = 0 <=> x = g(x) 

Notation:
g = lambda x: x**2-3*x**(1/2)+4 # Here input your recursive function
x0 = 1 # Your initial point
n = 10 # Your iteration times
residual = 0.001 # Your tolerance 
"""
def Recursive_Method(g,x0,n,residual):
    x = np.zeros(n+1)
    x[0] = x0
    for index in range(1,n+1):
        x[index] = g(x[index-1])
        difference = x[index] - x[index-1]
        if np.abs(x[index]-x[index-1])<residual:
            print("Iteration %r: x=%r, difference=%r <= Residual=%r" %(index,x[index],difference,residual))
            break
        else:
            print("Iteration %r: x=%r, difference=%r" %(index,x[index],difference))
    print("Terminal: The final result is %r" % x[index])
    return x[index]

""" Test
Recursive_Method(g,x0,n,residual)
g = lambda x: 0.5*(10-x**3)**0.5
Recursive_Method(g,1.5,10,10e-6)
"""

""" 
    Newton Method
    Input your function f(x) and f'(x) in x

Notation:
f = lambda x: x**2+3*x+1 # Your function f(x)
f1 = lambda x: 2*x+3 # First order deviation of f(x)
x0 = 1 # The initial point
n = 100 # Iteration time
residual = 10e-6 # Tolerance
"""
 
def Newton_Recursive(f,f1,x0,n,residual):
    x = np.zeros(n+1)
    x[0] = x0
    for index in range(1,n+1):
        x[index] = x[index-1]-f(x[index-1])/f1(x[index-1])
        difference = x[index] - x[index-1]
        if np.abs(x[index]-x[index-1])<residual:
            print("Iteration %r: x=%r, difference=%r <= Residual=%r" %(index,x[index],difference,residual))
            break
        else:
            print("Iteration %r: x=%r, difference=%r" %(index,x[index],difference))
    print("Terminal: The final result is %r" % x[index])
    return x[index]

""" Test 
Newton_Recursive(f,f1,x0,n,residual)
f = lambda x: x**3+4*x**2-10
f1 = lambda x: 3*x**2+8*x
Newton_Recursive(f,f1,1.5,10,10e-6)
"""