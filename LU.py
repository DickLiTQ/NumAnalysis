# -*- coding: utf-8 -*-

" In this pyfile we will give the solution of linear equation system in details. That is, not only the final result will be given, but also the L, U matrix and the determination of A will be listed here for your convinience. "

import numpy as np

""" Input: Ax=b 
A = np.array([[1,3,5],
              [2,5,2],
              [9,3,4]
              ]) # Here input your coefficient matrix
b = np.array([10,24,31]) # Here input the vector
"""

def check(A):
    row = A.shape[0]
    col = A.shape[1]
    if row!=col:
        print("Input error: A is not a square matrix")
        return 0
    else:
        if np.linalg.det(A)==0:
            print("The determination of A is equal to zero")
            return 0
        else:
            if row == 1:
                print("The dimension of matrix is 1")
                return 0
            else:
                return row
def Decomposition(A):
    if check(A) == 0:
        print("Error")
    else:
        print("det(A)=%r"%np.linalg.det(A))
        dim = check(A)
        L = np.eye(dim)    
        U = np.zeros_like(A)
        U[0,:]=A[0,:]
        L[1:,0]=A[1:,0]/U[0,0]
        for r in range(1,dim):
            for l in range(1,r):
                L[r,l]=1/U[l,l]*(A[r,l]-L[r,:l]@U[:l,l])
            for u in range(r,dim):
                U[r,u]=A[r,u]-L[r,:r]@U[:r,u]
    print("L=\n",L,"\n","U=\n",U)
    return L,U

#Decomposition(A)

def Solve(A,b):
    L,U = Decomposition(A)
    y = np.linalg.solve(L,b)
    x = np.linalg.solve(U,y)
    print("y=\n",y,"\n","x=\n",x)
    return y,x

#Solve(A,b)
    

    

