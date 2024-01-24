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

- for n in n_max // set a number of iterations
- c = (a + b) / 2 // calculate the midpoint
- if f(c) == 0 return c // if the midpoint is 0, a root has been found.
- sign = f(a) * f(c) // see if a root still exists between a and c
- If sign > 0 // determine which side of midpoint root is on
- a = c
- else
- b = c

As a note, it's highly recommended to check for a tolerance instead of absolute equality on step (3) of the algorithm in order to prevent excessive time spent finding a solution.
The Bisection Method is particularly slow, executing in logarithmic time, however, quite robust.

### Newton's Method
Newton's Method utilizes the derivative of a function in order to create tangent lines. These tangent lines intersect with y = 0, 


## Interpolation


### Lagrange Form


### Divided Differences


## Splines


### Natural Cubic Splines


### B Splines


### Bezier Curves


## Numerical Derivatives


### Truncated Taylor Series


### Richardson Extrapolation


## Numerical Integration


### Trapezoid Method


### Romberg Method


### Simpson's Rules


### Gaussian Quadrature


## Linear Systems


### Naive Gaussian Elimination


### Scaled Partial Pivoting


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


## Thanks and Acknowledgements
Thank you for taking the time to read this writeup! This has been a one-person project thus far, and as such, there may be errors. If you find one, please submit a pull request. If it's
merged, your name will be added below! Another thanks to another one of my professors, who unfortunately shall go unnamed due to the uniqueness of his name posing an OPSEC risk. This would
have never happened without his efforts.

### Merged PR Contributors
- Could be you!
