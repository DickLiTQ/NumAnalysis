## 数值分析/NumAnalysis

我们希望将《计算方法》课程中的理论内容通过Python语言进行实现，并将其运用到实际问题之中。
作为一个对计算机稍微有一定了解的学习者，深知将理论转化为实际往往是存在困难的。通常发展了一门新的理论，并将其转化为技术，这期间需要很长的时间来逐步完善并提高效率。

更重要的是，我们认为这门课程中的一些方法是形象生动而有趣的。可能在更高级的模块中集成了比这些先进得多的算法，但是这并不重要，我们希望用这本书中的方法来解决一些适合于我们的需求的问题。


This repository is used to create a project on lesson "Numerical Analysis" and practically use the method we have learned in problem solving.
As a learner who knows only a little about computer and programming, I know that there must exist difficulties in the transformation from theory to practice. It takes long time when the theory generally improves efficiency as technical tools.

From our perspective, some methods in this lesson are really interesting and fantastic. Maybe some much advanced technique has been used in popular modules, however, it is not important while we hope that the method in our textbook can solve the problem we faced and adapted to our demand.

## 目录/Table of Contents
* [非线性方程数值解/Numerical Solution to Nonlinear Equation](##非线性方程数值解/Numerical-Solution-to-Nonlinear-Equation)
	* [迭代法/Recursive Method](###迭代法/Recursive-Method)
	* [牛顿法/Newton Method](###牛顿法/Newton-Method)
	* [双点快速截弦法/](###双点快速截弦法/)


## 非线性方程数值解/Numerical Solution to Nonlinear Equation
### 迭代法/Recursive Method
基本原理
![](http://latex.codecogs.com/gif.latex?f(x)=0\\Rightarrow~x=g(x)~\\Rightarrow~x_{n+1}=g(x_n))
在这个部分中我们将利用循环完成迭代的操作，需要手动输入迭代公式![](http://latex.codecogs.com/gif.latex?g(x))，具体代码如下：
""" python
import numpy as np
import sympy as sp
g = lambda x: x**2-3*x**(1/2)+4 # Here input your recursive function
x0 = 1 # Your initial point
n = 10 # Your iteration times
residual = 0.001 # Your tolerance 

def Recursive_Method(g,x0,n,residual):
    x = np.zeros(n)
    x[0] = x0
    for index in range(1,n):
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
"""
### 牛顿法/Newton Method
### 双点快速截弦法/

<!--
	## 线性方程组LU分解/System of Linear Equation——LU Decomposition
	### 杜利特尔分解/Doolittle Decomposition
	### 形态分析/Property Analysis
	## 线性方程组迭代法/Recursive way to solve the System of Linear Equation
	### Jacobi迭代/Jacobi Recursion
	### Gauss-Seidel迭代/Gauss-Seidel Recursion
	## 插值与拟合/Interpolation and Fitting
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

