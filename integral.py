# -*- coding: utf-8 -*-
" We solve the numerical integral in this pyfile. Newton-Cotes integral rule, Composite Newton-Cotes integral rule are mentioned " 

import numpy as np

""" Input data Instance: 
data = np.array([[0,3],
                 [1,5],
                 [2,6],
                 [3,8]])    

"""

def trapezium(data,a,b):
    " Integral from ath data to bth data "
    #dim = b-a
    x = data[a-1,0]
    y = data[b-1,0]
    fx = data[a-1,1]
    fy = data[b-1,1]
    return (y-x)/2*(fx+fy)

# trapezium(data,1,2)

def sinpson(data,a,b):
    " Integral from ath data to bth data "
    dim = b-a
    if dim%2 != 0:
        return "Error"
    else:
        x = data[a-1,0]
        #xy = data[(b+a)/2-1,0]
        y = data[b-1,0]
        fx = data[a-1,0]
        fy = data[b-1,0]
        fxy = data[(b+a)/2-1,0]
        return (y-x)/6*(fx+4*fxy+fy)

"""
data = np.array([[0,1],
                 [0.125,0.997398],
                 [0.25,0.989688],
                 [0.375,0.976727],
                 [0.5,0.958851],
                 [0.625,0.936156],
                 [0.75,0.908858],
                 [0.875,0.877193],
                 [1,0.841471]])
"""

def com_trapezium(data):
    dim = data.shape[0]
    a = data[0,0]
    b = data[dim-1,0]
    h = (b-a)/(dim-1)
    result = data[0,1]+data[dim-1,1]
    for index in range(1,dim-1):
        result = result + 2*data[index,1]
    return result*h/2

def com_sinpson(data):
    dim = data.shape[0]
    a = data[0,0]
    b = data[dim-1,0]
    h = (b-a)/(dim-1)*2
    result = data[0,1]+data[dim-1,1]
    #times = (dim-1)/2
    for index in range(1,dim-1):
        if index%2 != 0:
            result = result + 4*data[index,1]
        else:
            result = result + 2*data[index,1]
    return result*h/6
        

#com_trapezium(data)
#com_sinpson(data)