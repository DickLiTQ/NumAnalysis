## 数值分析/NumAnalysis

我们希望将《计算方法》课程中的理论内容通过Python语言进行实现，并将其运用到实际问题之中。
作为一个对计算机稍微有一定了解的学习者，深知将理论转化为实际往往是存在困难的。通常发展了一门新的理论，并将其转化为技术，这期间需要很长的时间来逐步完善并提高效率。

更重要的是，我们认为这门课程中的一些方法是形象生动而有趣的。可能在更高级的模块中集成了比这些先进得多的算法，但是这并不重要，我们希望用这本书中的方法来解决一些适合于我们的需求的问题。


This repository is used to create a project on lesson "Numerical Analysis" and practically use the method we have learned in problem solving.
As a learner who knows only a little about computer and programming, I know that there must exist difficulties in the transformation from theory to practice. It takes long time when the theory generally improves efficiency as technical tools.

From our perspective, some methods in this lesson are really interesting and fantastic. Maybe some much advanced technique has been used in popular modules, however, it is not important while we hope that the method in our textbook can solve the problem we faced and adapted to our demand.

## 目录/Table of Contents
  * 非线性方程数值解/Numerical Solution to Nonlinear Equation
	* [迭代法/Recursive Method](#recursive-method)
	* [牛顿法/Newton Method](#newton-method)
	* [双点快速截弦法/](#双点快速截弦法/)
  * 线性方程组LU分解/System of Linear Equation——LU Decomposition
	* [杜利特尔分解/Doolittle Decomposition](#doolittle-decomposition)
	* [性态分析/Property Analysis](#property-analysis)
  * 线性方程组迭代法/Recursive way to solve the System of Linear Equation
    * [Jacobi迭代/Jacobi Recursion](#jacobi-recursion)
	* [Gauss-Seidel迭代/Gauss-Seidel Recursion](#gauss-seidel-recursion)

## 非线性方程数值解/Numerical Solution to Nonlinear Equation

### 迭代法
#### Recursive Method 
基本原理
![](http://latex.codecogs.com/gif.latex?f(x)=0\\Rightarrow~x=g(x)~\\Rightarrow~x_{n+1}=g(x_n))
在这个部分中我们将利用循环完成迭代的操作，需要手动输入迭代公式![](http://latex.codecogs.com/gif.latex?g(x))，具体代码如下：
``` python
import numpy as np
import sympy as sp
g = lambda x: x**2-3*x**(1/2)+4 # Here input your recursive function
x0 = 1 # Your initial point
n = 10 # Your iteration times
residual = 0.001 # Your tolerance 

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

Recursive_Method(g,x0,n,residual)
```
例如我们要使用![](http://latex.codecogs.com/gif.latex?x=\frac{1}{2}(10-x^3)^{\frac{1}{2}})求解![](http://latex.codecogs.com/gif.latex?x^3+4x^2-10=0)在![](http://latex.codecogs.com/gif.latex?[0,1])上的根，迭代10次或两次误差小于![](http://latex.codecogs.com/gif.latex?10^{-5})，则使用以下的代码
``` python
g = lambda x: 0.5*(10-x**3)**0.5
Recursive_Method(g,1.5,10,1e-5)
```
得到的结果为：
```
Iteration 1: x=1.2869537676233751, difference=-0.2130462323766249
Iteration 2: x=1.4025408035395783, difference=0.11558703591620323
Iteration 3: x=1.3454583740232942, difference=-0.057082429516284172
Iteration 4: x=1.3751702528160383, difference=0.029711878792744173
Iteration 5: x=1.3600941927617329, difference=-0.015076060054305396
Iteration 6: x=1.3678469675921328, difference=0.0077527748303998223
Iteration 7: x=1.3638870038840212, difference=-0.0039599637081115802
Iteration 8: x=1.3659167333900399, difference=0.0020297295060187626
Iteration 9: x=1.3648782171936771, difference=-0.0010385161963628597
Iteration 10: x=1.3654100611699569, difference=0.00053184397627981106
Terminal: The final result is 1.3654100611699569
```

**注意**
1. 本部分尚未加入是否发散的判断，因此要根据迭代结果进行发散的判断
2. 请自行计算迭代式![](http://latex.codecogs.com/gif.latex?x_{n+1}=g(x_n))，在输入过程中注意指数与根号的输入规范

### 牛顿法
#### Newton Method
具有特殊迭代格式的迭代法，对于![](http://latex.codecogs.com/gif.latex?f(x)=0)的求解，采用格式：

![](http://latex.codecogs.com/gif.latex?g(x)=x-\dfrac{f(x)}{f'(x)})

能保证更快的收敛速度。这里，我们为了得到更为精确的结果，避免数值微分带来的误差和符号计算可能碰到的不确定性，要求手动输入函数![](http://latex.codecogs.com/gif.latex?f(x),f'(x))，并同理使用[迭代法](#recursive-method)中的方法：
``` python
f = lambda x: x**2+3*x+1 # Your function f(x)
f1 = lambda x: 2*x+3 # First order deviation of f(x)
x0 = 1 # The initial point
n = 100 # Iteration time
residual = 1e-5 # Tolerance

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

Newton_Recursive(f,f1,x0,n,residual)
```
作为对比，我们依旧解决迭代法中10次尚未达到精度的例子：
>求解![](http://latex.codecogs.com/gif.latex?x^3+4x^2-10=0)在![](http://latex.codecogs.com/gif.latex?[0,1])上的根，迭代10次或两次误差小于![](http://latex.codecogs.com/gif.latex?10^{-5})。
此时，在牛顿法中，![](http://latex.codecogs.com/gif.latex?f(x)=x^3+4x^2-10,f'(x)=3x^2+8x)。对应代码为：
``` python
f = lambda x: x**3+4*x**2-10
f1 = lambda x: 3*x**2+8*x
Newton_Recursive(f,f1,1.5,10,1e-5)
```
仅通过三次迭代就得到收敛结果
```
Iteration 1: x=1.3733333333333333, difference=-0.12666666666666671
Iteration 2: x=1.3652620148746266, difference=-0.0080713184587066777
Iteration 3: x=1.3652300139161466, difference=-3.2000958479994068e-05
Iteration 4: x=1.3652300134140969, difference=-5.0204973511824846e-10 <= Residual=1e-05
Terminal: The final result is 1.3652300134140969
```
### 双点快速截弦法/

## 线性方程组LU分解/System of Linear Equation——LU Decomposition
本部分求解的线性方程组为行列式不为0的实方阵，用数学符号来表示，为：

![](http://latex.codecogs.com/gif.latex?Ax=b)

其中![](http://latex.codecogs.com/gif.latex?A\in\mathbb{R}^{n~\times~n},b\in\mathbb{R}^n,det(A)\neq0)
### 杜利特尔分解
#### Doolittle Decomposition
我们考虑将系数矩阵![](http://latex.codecogs.com/gif.latex?A\in\mathbb{R}^{n~\times~n})分解为一个同维度的下三角矩阵![](http://latex.codecogs.com/gif.latex?L)和一个同维度的上三角矩阵![](http://latex.codecogs.com/gif.latex?U)的乘积，即

![LU Decomposition](https://raw.githubusercontent.com/DickLiTQ/NumAnalysis/4ad6e37ba2565a0f628ed0965fefe81593c84547/LUDecomposition.gif)

之所以这样做，是因为对于![](http://latex.codecogs.com/gif.latex?Ly=b,~Ux=y)的求解是相对容易的，可以先通过求解![](http://latex.codecogs.com/gif.latex?y)再来求解![](http://latex.codecogs.com/gif.latex?x)。

首先我们先对系数矩阵A进行判断并返回维度供后续使用
``` python
import numpy as np
A = np.array([[1,3,5,9],
              [2,5,2,5],
              [9,3,4,1],
              [1,10,0,19]]) # Here input your coefficient matrix
b = np.array([10,24,31,42]) # Here input the vector

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
```
获得检验矩阵并获得维度信息后，再根据矩阵乘法定义进行LU分解：
``` python
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
```
例如我们对下方矩阵A进行分解：
``` python
A = np.array([[1,3,5],
              [2,5,2],
              [9,3,4]
              ])
Decomposition(A)
```
得到结果为
``` python
det(A)=-151.0
L=
 [[  1.   0.   0.]
 [  2.   1.   0.]
 [  9.  24.   1.]] 
 U=
 [[  1   3   5]
 [  0  -1  -8]
 [  0   0 151]]
 ```
 再分别求解：![](http://latex.codecogs.com/gif.latex?Ly=b,~Ux=y)
 ``` python
 def Solve(A,b):
    L,U = Decomposition(A)
    y = np.linalg.solve(L,b)
    x = np.linalg.solve(U,y)
    print("y=\n",y,"\n","x=\n",x)
    return y,x
 ```
### 性态分析
#### Property Analysis

## 线性方程组迭代法/Recursive way to solve the System of Linear Equation
目标依然是求解线性方程组![](http://latex.codecogs.com/gif.latex?Ax=b,~A\in\mathbb{R}^{n~\times~n},~b\in\mathbb{R}^n)
此处的迭代法类似于[迭代法](#recursive-method)，但迭代对象换成了矩阵和向量，依旧选择合适的 ![](http://latex.codecogs.com/gif.latex?g(\cdot)) 使得

![](http://latex.codecogs.com/gif.latex?X_{n+1}=g(X_n),~X\in\mathbb{R}^n)

以上迭代思想导出不同收敛速度的Jacobi迭代法和Gauss-Seidel迭代法
### Jacobi迭代
#### Jacobi Recursion
### Gauss-Seidel迭代
#### Gauss-Seidel Recursion

<!--	## 插值与拟合/Interpolation and Fitting
	### Lagrange插值/Lagrangian Interpolation
	### Newton插值/Newton Interpolation
	### 最小二乘拟合/Least Square
	## 数值积分/Numerical Integral
	### 插值型数值积分
	### Newton-Cotes公式/Newton-Cotes Method
	#### 梯形公式
	#### Sinpson公式
	### 复化求积公式
	#### 复化梯形公式
	#### 复化Sinpson公式
	### Gauss求积公式
	## 常微分方程数值解/Numerical Solution in Ordinary Differential Equation
	### Euler法/Euler Method
	### 改进的Euler法/Improved Euler Method
-->

