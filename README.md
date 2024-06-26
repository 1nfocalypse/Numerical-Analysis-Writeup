<p align="center">
  <a href="https://github.com/1nfocalypse/Numerical-Analysis-Writeup">
	<img alt="Numerical Analysis" src="https://i.imgur.com/nt6Lu1t.png"/>
  </a>
</p>
<p align="center">
  <a href="https://choosealicense.com/licenses/cc-by-sa-4.0/">
  	<img alt="License: CC-BY-SA-4.0" src="https://img.shields.io/github/license/1nfocalypse/Numerical-Analysis-Writeup"/>
  </a>
</p>
<h2 align="center">Numerical Analysis</h3>
<h3 align="center">
  A Writeup on Introductory Numerical Methods and Techniques
</h2>
<p align="center">
  By <a href="https://github.com/1nfocalypse">1nfocalypse</a>
</p>


## Contents

- [Summary](#summary)
- [Floating Point Numbers](#floating-point-numbers)
- [Nonlinear Equations](#nonlinear-equations)
- [Interpolation](#interpolation)
- [Splines](#splines)
- [Numerical Derivatives](#numerical-derivatives)
- [Numerical Integration](#numerical-integration)
- [Linear Systems](#linear-systems)
- [Initial Value Problems](#initial-value-problems)
- [Boundary Value Problems](#boundary-value-problems)
- [Partial Differential Equations](#partial-differential-equations)
- [Conclusion](#conclusion)
- [Thanks and Acknowledgements](#thanks-and-acknowledgements)

## Summary

Numerical analysis is a field of study regarding approximation of solutions to problems via algorithms. Via this method, we are able to obtain approximated solutions to complex problems 
quickly with computers, rather than attempting to analytically solve them. As such, it is incredibly relevant to several fields of Computer Science research and industry. Examples of this 
include Computational Fluid Dynamics, Finance, Physics calculations in Engineering, such as Finite Element Analysis, and many more.

This document is intended to serve as a reference guide along with as a resource for introductory numerical analysis. While it may be useful to utilize it to supplement for self study, it 
is not intended to replace an appropriate curriculum. This document is also not meant to be a comprehensive guide to numerical analysis, and as such, not all methods or types of problems 
will be discussed, i.e. certain optimization techniques, the multivariate case, etc. There will, however, be recommended reading delineated at the bottom of this guide that will be much 
more comprehensive. This guide assumes prior knowledge of Calculus I, II, Multivariable Calculus, Linear Algebra, and Ordinary Differential Equations. Brief refreshers may be provided, but 
it is strongly advised to have background in this prior to beginning a study of numerical analysis.

Most informationation is sourced from [Numerical Mathematics and Computing: 7ED](https://www.amazon.com/Numerical-Mathematics-Computing-published-Learning/dp/B00E28JH28/), 
Wikipedia, and personal notes from a university mathematics course. Some of the methods here are or will be implemented in the EnkiNum python library, found [here](https://github.com/1nfocalypse/EnkiNum)
. This is an ongoing project, and as such, may not be complete at the time of reading.

## Floating Point Numbers
Obviously, computers utilize binary as a means of storing information. A consequence of this is that it becomes necessary to store decimals in base two, which is usually done in 
accordance with [IEEE standard 754](#https://en.wikipedia.org/wiki/IEEE_754). However, given the limited space in which a decimal is able to be stored, this leads to problems with
exactness, as while the format could approach exactness if afforded enough space, it is limited by practicality. Some of you may have encountered this, such as with the number 0.2 not 
being able to be exactly stored in the single or double precision float. However, approximate solutions are often suitable for practical usage, and as such, this field provides a large
amount of value.

### Accuracy, Precision, and Error
Accuracy is the term used to describe the amount of precision, i.e. how many decimal places, are accurate and to be trusted. This refers to the level of precision that a method is capable
of producing. Error is defined in two manners, which are absolute and relative errors. Absolute error is defined as $E_{a} = |Exact - Measured|$, and relative error is defined as
$E_{r} = \frac{|Exact - Measured|}{|Exact|}$. In both cases, exact is the best known value or the analytical solution.

### Rounding and Chopping
Rounding is the familiar process of approximating a solution via changing a digit to obtain a certain level of precision. For example, 9.4 would be rounded to 9, and 9.6 would be rounded 
to 10. However, there are alternative means of rounding, such as round to even or round to odd, in which the number is rounded to the closest even or odd number, respectively. This would
cause 9.4 to be rounded to 9 with round to odd, and 10 to round to even. This is done to reduce potential stastical bias, among other reasons. Chopping is simply the complete removal of 
numbers after a certain point, i.e. 7.5472 chopped to two decimal places would be 7.54, whereas typically rounded it would be 7.55.

### Subtractive Cancellation 
Subtractive cancellation is a phenomenon in which when two nearly identicial numbers are subtracted, a significant loss of precision occurs due to chopping. This is unique to floating
point representations, and can cause significant error in calculations. To combat this, calculations are usually rearranged in order to prevent this subtraction operation from occuring. 

## Nonlinear Equations
We will begin a discussion of algorithms and methods by talking about finding solutions to nonlinear equations. We will discuss two root finding algorithms, which are the Bisection method
and Newton's Method, also known as the Newton-Raphson Method. These two methods are both approaches to approximating solutions, or roots, of nonlinear equations, and are particularly 
useful in exemplifying differences in speed from different algorithms that arrive at similar ends to further drive home the diversity of this field.

### Bisection Method
The Bisection Method, or alternatively, the Binary Search method, is a logarithmic time algorithm used to locate roots of nonlinear functions. It operates by taking a function, a lower x
value a, and an upper x value b. Necessarily, there must exist a root between these two bounds. A simple test for this is to evaluate if the function changes sign between them, i.e.
f(a) > 0 and f(b) < 0 or f(a) < 0 and f(b) > 0. Once this has been established, the following algorithm is conducted.
```
- for n in n_max // set a number of iterations
- c = (a + b) / 2 // calculate the midpoint
- if f(c) == 0 return c // if the midpoint is 0, a root has been found.
- sign = f(a) * f(c) // see if a root still exists between a and c
- If sign > 0 // determine which side of midpoint root is on
- a = c
- else
- b = c
```
As a note, it's highly recommended to check for a tolerance instead of absolute equality on step (3) of the algorithm in order to prevent excessive time spent finding a solution.
The Bisection Method is particularly slow, executing in logarithmic time, however, quite robust.

### Newton's Method
Newton's Method utilizes the derivative of a function in order to create tangent lines. These tangent lines intersect with y = 0, which then restarts the iteration. The algorithm requires
an initial guess $x_{0}$, a function, and the derivative of said function. Failure conditions for the function are cases where the initial guess lies on a minimum or maximum value in the
original function, or when the original function is not differentiable. It also slows down considerably when the original function is expensive to accurately differentiate.
The pseudocode for the algorithm is as follows:
```
- if f(x) == 0 return x // if guess is a root, return
- while f(x) > 0.001 or f(x) < -0.001 // set tolerance bounds aka convergence critera
- x = x - (f(x) / f'(x)) // establish new x
```
This pseudocode performs no checks for validity, but edge cases aside, will function. This algorithm performs best when the initial guess is close to the actual root. As such, it is not 
uncommon to find hybrid implementations featuring the bisection method to provide a very rough approximation quickly, and then handing that approximation to a Newton's method algorithm to
find a much more precise approximation. Newton's method converges in quadratic time, which is markedly superior to the linear convergence of the Bisection method. 

## Interpolation
Interpolation is a technique used to find a polynomial that intersects a given set of points, with at most degree $n$, where $n$ is the number of given points. This is given by the 
existence and uniqueness theorem for interpolating polynomials, which states that there necessarily exists a unique polynomial of at most degree $n$. Interpolation has many uses in that
it provides a polynomial that fits data, however, it also has a very unique problem in the sense that it is subject to something called the Runge Phenomenon, which causes intense
oscillations in the polynomial when a high number of nodes are interpolated, generating a great deal of noise. This is fixed via something called splines, which will be touched on later.
This is not to say interpolation does not find use; it is still heavily employed in graphics, data analysis, physics, engineering, signal processing, and more, and as such, will be 
discussed below. 

### Lagrange Form
The Lagrange Form of an interpolating polynomial is a relatively straightforward means of interpolation. However, it does not adapt well, with the addition of any new points requiring
a complete recalculation of the polynomial and allowing for larger degrees of error. However, it adequately serves for low-degree tasks. Lagrangian Interpolation is described as follows:

The Lagrangian Interpolating Polynomial is comprised of Lagrange basis polynomials $L_{n}$ and the y values of the data points as 
$P(x)=L_{0}(x) * y_{0} +L_{1}(x)*y_{1}+…+L_{n}(x)*y_{n}$.

Lagrange basis polynomials are defined as $L_i(x) = \displaystyle\prod_{j=0, j \neq i}^{n}\frac{x-x_{i}}{x_{i}-x_{j}}$. This gives them the unique quality of the Kronecker delta, meaning that in 
the case that $i=j$, the value is 1, and is otherwise 0, meaning that only one term survives when a point is evaulated, guaranteeing interpolation. Pseudocode for Lagrange interpolation
is as follows, creating a list of the Lagrangian Basis Polynomials as coefficients to be used later.
```
- Requires x values, y values
- n = length(x_values)
- coefficients = initialize_list_of_zeros(n)
- for i from 0 to n-1:
-   term = y_values\[i]
-   for j from 0 to n-1:
-       if j != i:
-       term = term * 1 / (x_values[i] - x_values[j])
-   for k from 0 to n-1:
-       if k != i:
-       coefficients[k] = coefficients[k] + term * (-x_values[k]) / (x_values[i] - x_values[k])
```
This provides a list of coefficients that can later be used to calculate arbitrary values along the function. However, these will need to be recalculated if more points are added that need
to be interpolated, and is also subject to the aforementioned Runge Phenomenon.


## Splines
Splines are an alternative means of interpolating data via piecewise functions. Splines are useful because they avoid the Runge Phenomenon that occurs in interpolation, as well as can
maintain certain characteristics, such as differentiability. Splines find heavy use in contemporary computing, seeing heavy employment in computer graphics, engineering, CAD software, 
Finite Element Analysis, Data Visualization, game development, and more.

### Natural Cubic Splines
Natural Cubic Splines are both optimal and unique given a set of points. They are optimal in the sense that they minimize the curve required to get between two points, unique in the sense
that they are thus minimized, and cubic in the sense that the piecewise polynomials are cubic in nature. To further simplify, Natural Cubic Splines are a set of cubic polynomials between
a set of points that interpolate these points while minimizing the amount of bend between them, which in turn makes this interpolation unique. Natural Cubic Splines, as a result of strict
boundary conditions and their cubic nature, must have $C^{2}$ continuity. The boundary conditions for Natural Cubic Splines are such that there is zero curvature at the given points, such
that it remains twice differentiable. The equation for a Natural Cubic Spline is defined as:

$S_{i}(x) = a_{i} + b_{i}(x-x_{i}) + c_{i}(x-x_{i})^{2} + d_{i}(x-x_{i})^3$

This results in a tridiagonal banded matrix, which can be solved for the coefficients $a_{i}, b_{i}, c_{i}, d_{i}$.

Pseudocode for solving for the coefficients of a spline is as follows:
```
- Requires data points (x,y), and real arrays x,y,z such that x vals are in x, y vals in y, and z is empty.
- n = size(data) - 1
- allocate arrays d, e, f, g w/ size n-1
- for i in range(n-1)
-   d[i] = x[i + 1] - x[i]
-   e[i] = (y[i+1]-y[i])/d[i]
- f[1] = 2(a[0] + a[1])
- g[1] = 6(e[1] - e[0])
- for i = 2 in range(n-1)
-   f[i] = 2(d[i] + d[i-1]) - ((d[i-1] * d[i-1]) / f[i-1])
-   g[i] = 6(e[i] -e[i-1]) - ((d[i-1] * g[i-1]) / f[i-1])
- z[n] = 0
- for i in range(n-1, 1, -1)
-   z[i] = (g[i] - (x[i] * z[i+1]))/f[i]
- z[0] = 0
```

Once this algorithm has been run, you are able to evaluate the value of the spline at actual locations. The following pseudocode details how this is accomplished.
```
- Requires n data points, real arrays x,y,z, h, which is step size, and x, which is the position to be evaluated
- int tmploc
- for i in range(n-1, 0, ,-1)
-   if x - t[i] > 0
-       tmploc = i
-       break;
- i = tmploc
- h = t[tmploc + 1] - t[tmploc]
- tmp = (z[i] / 2) + (((x - t[i]) * (z[i+1] - z[i]))/(6 * h))
- tmp = -(h/6)(z[i+1] + 2 * z[i]) + ((y[i+1] - y[i])/h) + ((x-t[i]) * tmp)
- return y[i] + ((x - t[i]) * tmp)
```

At this point, you will have constructed a natural cubic spline for a function as well as a means to locate specific points on the spline. This allows for interpolation to occur without
the aforementioned Runge Phenomenon, as well as ensures $C^{2}$ differentiability.

### B Splines
B splines, or Basis Splines, also known as Bell Splines from their bell shape in higher orders, are another form of spline, utilized for similar purposes as NCSs. However, B splines of higher 
order are able to offer more control over differentiability than Natural Cubics, along with some finer control of curve type. These properties have made them popular in production and research
where higher orders or finer control are needed, and optimality of curve is not. B splines derive their name from the fact that they form a basis of spline functions. Drawing back from linear
algebra, this denotes that all possible spline functions can be built from linear combinations of B splines, and for each function, a unique linear combination of B splines exists to form it.
B splines are defined recursively by the Cox-de Boor function, which is elaborated below along with the recursive definition for higher order splines.

A basis spline $B^{0}\_{i}$ is defined as $\\{ 1 : t_{i} \leq x < t_{i + 1} || 0 : otherwise\\}$.
From this, we find the recursive definition:


$B^{k}\_{i}(x) = (\frac{x-t_{i}}{t_{i + k}-t_{i}})B^{k-1}\_{i}(x) + (\frac{t_{i+k+1} - x}{t_{i + k + 1} - t_{i + 1}})B^{k-1}\_{i+1}(x)$ for $k \geq 1$.

### Bezier Curves


## Numerical Derivatives
Numerical differentiation is highly pertinent, as calculus is often imperative in engineering tasks and simulations. However, given the discrete nature of computers, we find ourselves with
limited precision, meaning that the analytic path towards calculus implementation is not possible. Some may be thinking back to their introductory calculus course, remembering the basic
definition of a derivative as $f(x) = \frac{f(x + h) - f(x)}{h}$. However, in the context of this writeup, the problem becomes immediately evident with this definition. For a sufficiently
small value $h$, we find ourselves quickly speeding down the road to errors and instability. How then, are we expected to find derivatives of continuous functions? The rest of this section
discusses methods to approximate derivatives within the bounds of computers.

### Symmetric Differences


### Truncated Taylor Series


### Richardson Extrapolation


## Numerical Integration


### Trapezoid Method


### Romberg Method


### Simpson's Rules


### Gaussian Quadrature


## Linear Systems
Linear systems frequently occur in natural life, engineering, computer graphics, flight, and many other facets of modern life. As such, the ability to solve non-trivial linear systems 
with computers is incredibly important in many aspects of contemporary life. While seemingly innocuous enough, linear systems are able to pose a distinct challenge to computers due to
their ability to contain small values and the need for subtraction as a method of solving these systems. From this, we necessarily observe that traditional, naive Gaussian elimination,
as you likely learned in your linear algebra course, is not sufficient. There are other methods available, however, that do not encounter the same pitfalls that you can imagine Naive
Gaussian Elimination does, which are further discussed below, along with some notes on various matrix forms, operations, and characteristics. 


### Scaled Partial Pivoting


### Determinants


### Computing Inverses


### Tridiagonal and Banded Matrices


### LU Decomposition


### LDLT Decoposition


### Cholesky Decomposition


### Eigenvalues


### Eigenvectors


### Single-Value Decomposition


### Power Method


### Iterative Methods


## Initial Value Problems


### Truncated Taylor Series


### Runge-Kutta Methods


### Adams-Bashforth-Moulton Methods


## Boundary Value Problems


### Shooting Method


### Finite Differences


## Partial Differential Equations


### Parabolic Problems


### Hyperbolic Problems


### Elliptic Problems


## Conclusion
Numerical analysis is a very widely utilized and deep field, with active research occurring globally today. This writeup is far from all-encompassing, however, that was not its purpose.
I hope that it has inspired you to further look into the subject, along given you more tools for approaching problems and solving them. As mentioned, please check out my python library,
[EnkiNum](https://github.com/1nfocalypse/EnkiNum), as it implements these techniques in an approachable way such that the reader can understand application. If you have any concerns or 
questions, please open an issue, pull request, or otherwise contact me, and I will be happy to resolve them. 

## Thanks and Acknowledgements
Thank you for taking the time to read this writeup! This has been a one-person project thus far, and as such, there may be errors. If you find one, please submit a pull request. If it's
merged, your name will be added below! Another thanks to another one of my professors, who unfortunately shall go unnamed due to the uniqueness of his name posing an OPSEC risk. This would
have never happened without his efforts.

### Merged PR Contributors
- Could be you!
