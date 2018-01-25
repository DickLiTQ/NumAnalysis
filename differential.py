# -*- coding: utf-8 -*-

" This pyfile is used to give a numerical solution of ordinary differential equation "

import numpy as np
import matplotlib.pyplot as plt

""" Input: 
Interval: [a,b]
Step: h
Differential: f(x,y)
Initial point: y(a) = y0
"""

def EulerMethod(a,b,h,f,y0):
    if b>a:
        step = int((b-a)/h+1)
        x = np.linspace(a,b,step)
        y = np.zeros_like(x)
        y[0]=y0
        for index in range(1,step):
            y[index] = y[index-1] + h*f(x[index-1],y[index-1])
        return y
    else:
        return "Error"
    
def I_EulerMethod(a,b,h,f,y0):
    if b>a:
        step = int((b-a)/h+1)
        x = np.linspace(a,b,step)
        y = np.zeros_like(x)
        ye = np.zeros_like(x)
        y[0] = y0
        ye[0] = y0
        for index in range(1,step):
            ye[index] = ye[index-1] + h*f(x[index-1],ye[index-1])
            y[index] = y[index-1] + 0.5*h*(f(x[index-1],y[index-1])+f(x[index],ye[index]))
        return y
    else:
        return "Error"
    
def Draw(a,b,h,f,y0,ft):
    step = int((b-a)/h+1)
    x = np.linspace(a,b,step)
    fig, ax = plt.subplots(dpi=120)
    ye = EulerMethod(a,b,h,f,y0)
    yie = I_EulerMethod(a,b,h,f,y0)
    yt = ft(x)
    ax.plot(x,ye,marker='o',color='r',alpha=0.5,label='Euler Method')
    ax.plot(x,yie,marker='x',color='b',alpha=0.5,label='Improved Euler Method')
    ax.plot(x,yt,marker='v',color='g',alpha=0.5,label='Explicit Function')
    ax.legend()
    fig.show()
    return 0
    
"""
a = 0
b = 1
h = 0.1
f = lambda x,y: y-2*x/y
y0 = 1
ft = lambda x: np.sqrt(1+2*x)

EulerMethod(a,b,h,f,y0)
I_EulerMethod(a,b,h,f,y0)
Draw(a,b,h,f,y0,ft)
"""