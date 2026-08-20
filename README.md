# Numerical Methods in Python

A collection of numerical methods and algorithms implemented in Python for solving problems in numerical analysis, linear algebra, interpolation, numerical integration, ordinary differential equations, and partial differential equations.

The repository is primarily intended as a learning and reference resource, with implementations designed to illustrate the underlying numerical algorithms rather than replace established scientific-computing libraries.

## Overview

This project contains implementations of fundamental numerical methods, including:

* Nonlinear equation solving
* Linear systems and matrix factorization
* Numerical interpolation and approximation
* Cubic splines
* Numerical integration and quadrature
* Ordinary differential equations (ODEs)
* Partial differential equations (PDEs)
* Iterative methods for linear systems
* Numerical linear algebra

The repository also includes a Jupyter Notebook covering several numerical linear algebra algorithms with examples and comparisons against NumPy implementations.

## Contents

### 1. Nonlinear Equations

Methods for finding numerical solutions to nonlinear equations and systems:

* Newton–Raphson method for nonlinear equations
* Newton's method for systems of nonlinear equations
* Modified Newton–Raphson methods

**Relevant files:**

* `newton_raphson_eq_no_lineals.py`
* `metode_newton_eq_no_lineals.py`

### 2. Least Squares and Interpolation

Methods for approximating and interpolating data:

* Polynomial least-squares approximation
* Polynomial interpolation
* Error analysis and comparison between approximation methods

**Relevant file:**

* `minims_quadrats_i_interpolacio.py`

### 3. Cubic Splines

Implementations and examples involving spline interpolation:

* Cubic spline construction
* Derivative approximation using finite differences
* Cubic Hermite splines
* `C2` cubic splines

**Relevant file:**

* `splines.py`

### 4. Numerical Quadrature

Numerical integration methods, including:

* Gaussian quadrature
* Composite trapezoidal rule
* Composite Simpson's rule
* Composite Gaussian quadrature

**Relevant files:**

* `quadraturaGauss.py`
* `quadratures.py`
* `quadratures_trapeziCompost_SimpsonComp_GaussComp.py`

### 5. Ordinary Differential Equations

Numerical methods for initial-value and boundary-value problems:

* Explicit Euler method
* Backward Euler method
* Fourth-order Runge–Kutta (RK4)
* Euler's method for systems of ODEs
* Heun's method
* Shooting method

**Relevant files:**

* `EDOs_RungeKutta4_Euler.py`
* `EulerEdos.py`
* `Edo_Heune.py`
* `shootingMethod.py`

For example, the ODE implementations include explicit Euler, RK4, and an implicit backward Euler method, with error comparisons against analytical or reference solutions.

### 6. Partial Differential Equations

Finite-difference approaches for PDE problems, including:

* One-dimensional diffusion equation
* Explicit Euler time integration
* Dirichlet boundary conditions
* Neumann boundary conditions
* Explicit and implicit discretizations
* Finite-difference methods

**Relevant files:**

* `EDPs_Neumann_Dirichlet.py`
* `EDPs_euler.py`
* `FD1D_diffusionEq.py`

### 7. Numerical Linear Algebra

The repository includes a dedicated Jupyter Notebook covering fundamental numerical linear algebra algorithms.

Topics include:

#### Triangular Systems

* Forward substitution
* Backward substitution
* Forward substitution for unit lower-triangular matrices

#### Gaussian Elimination and LU Factorization

* Gaussian elimination
* LU factorization
* Doolittle factorization
* LU factorization with partial pivoting
* Partial pivoting
* Matrix inversion using LU factorization

#### Matrix Factorizations

* Cholesky factorization
* `A = LDLᵀ` factorization
* QR factorization
* Modified Gram–Schmidt
* Householder transformations

#### Vector and Matrix Norms

* 1-norm
* 2-norm
* Infinity norm
* p-norm
* Matrix 1-norm
* Matrix infinity norm
* Frobenius norm

#### Iterative Methods

* Jacobi method
* Gauss–Seidel method
* Successive Over-Relaxation (SOR)
* Gradient method
* Conjugate gradient method

**Notebook:**

* `algebra_lineal_numerica.ipynb`

The notebook implements these algorithms directly and includes numerical examples and comparisons with NumPy routines.

## Requirements

The code is written in Python and uses standard scientific-computing packages such as:

* [NumPy](https://numpy.org/)
* [SciPy](https://scipy.org/)
* [Matplotlib](https://matplotlib.org/)
* [Jupyter](https://jupyter.org/) for the notebook
