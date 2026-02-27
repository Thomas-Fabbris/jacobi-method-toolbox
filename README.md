# Jacobi Method Toolbox

A MATLAB® toolbox that provides a robust implementation of the **Jacobi iterative method** for solving large-scale systems of linear equations. The toolbox supports both standard MATLAB arrays and **distributed arrays**, enabling parallel computation via MATLAB Parallel Computing Toolbox.

---

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
  - [Syntax](#syntax)
  - [Input Arguments](#input-arguments)
  - [Output Arguments](#output-arguments)
  - [Convergence Flags](#convergence-flags)
- [Examples](#examples)
  - [Basic Usage](#basic-usage)
  - [Distributed Array Usage](#distributed-array-usage)
- [Testing](#testing)
- [License](#license)

---

## Overview

The **Jacobi Method** is a classical iterative algorithm for solving a system of linear equations of the form:

```
A·x = b
```

where `A` is an *n×n* square coefficient matrix and `b` is a column vector of length *n*. The method converges when `A` is **strictly row diagonally dominant**.

This toolbox exposes a `jacobi` function with an interface consistent with MATLAB's built-in iterative solvers (e.g., `pcg`, `gmres`), including support for tolerance, maximum iterations, initial guess, and detailed output flags.

---

## Features

- ✅ Solves dense or sparse linear systems `A·x = b` using the Jacobi iterative method
- ✅ Supports **single** and **double** precision floating-point inputs
- ✅ Full support for **distributed arrays** (MATLAB Parallel Computing Toolbox)
- ✅ Configurable tolerance, maximum iterations, and initial guess
- ✅ Returns convergence flag, relative residual, iteration count, and residual history
- ✅ Interface consistent with MATLAB built-in iterative solvers

---

## Requirements

| Requirement | Version |
|---|---|
| MATLAB® | R2025a or later |
| Parallel Computing Toolbox *(optional)* | For distributed array support |

---

## Installation

1. Clone or download this repository:
   ```bash
   git clone https://github.com/Thomas-Fabbris/jacobi-method-toolbox.git
   ```
2. Open MATLAB and navigate to the `Code/jacobi` directory.
3. Add the toolbox folder to the MATLAB path:
   ```matlab
   addpath('toolbox')
   ```
4. Alternatively, open `jacobi.prj` in MATLAB to install the toolbox directly from the project.

---

## Usage

### Syntax

```matlab
x = jacobi(A, b)
x = jacobi(A, b, tol)
x = jacobi(A, b, tol, maxit)
x = jacobi(A, b, tol, maxit, x0)
[x, flag] = jacobi(___)
[x, flag, relres] = jacobi(___)
[x, flag, relres, iter] = jacobi(___)
[x, flag, relres, iter, resvec] = jacobi(___)
```

### Input Arguments

| Argument | Description | Default |
|---|---|---|
| `A` | Square coefficient matrix (floating-point, real) | *(required)* |
| `b` | Right-hand side column vector (floating-point, real) | *(required)* |
| `tol` | Convergence tolerance | `1e-6` (double), `1e-3` (single) |
| `maxit` | Maximum number of iterations | `min(n, 100)` |
| `x0` | Initial guess (column vector of length *n*) | zeros vector |

Pass `[]` for any optional argument to use its default value.

### Output Arguments

| Argument | Description |
|---|---|
| `x` | Approximate solution vector |
| `flag` | Convergence flag (see table below) |
| `relres` | Relative residual: `norm(b - A*x) / norm(b)` |
| `iter` | Iteration number at which `x` was computed |
| `resvec` | Vector of residual norms at each iteration |

### Convergence Flags

| Flag | Meaning |
|---|---|
| `0` | `jacobi` converged to the desired tolerance |
| `1` | `jacobi` did not converge — maximum iterations reached |
| `2` | `jacobi` stopped — a zero was found on the diagonal of `A` |

> **Note:** When `flag` is omitted from the output, `jacobi` prints a diagnostic message to the command window.

---

## Examples

### Basic Usage

Solve a system of linear equations using a diagonally dominant matrix:

```matlab
A = gallery('poisson', 50);
b = full(sum(A, 2));
tol   = 1e-9;
maxit = 500;

[x, flag, relres, iter] = jacobi(A, b, tol, maxit);

fprintf('Converged: flag=%d, relres=%.2e, iterations=%d\n', flag, relres, iter);
```

### Distributed Array Usage

Leverage MATLAB Parallel Computing Toolbox for large-scale problems:

```matlab
A  = gallery('poisson', 50);
b  = full(sum(A, 2));
dA = distributed(A);
db = distributed(b);

tol   = 1e-9;
maxit = 500;

dX = jacobi(dA, db, tol, maxit);
x  = gather(dX);
```

---

## Testing

Unit tests are located in `Code/jacobi/tests/testJacobi.m`. To run them in MATLAB:

```matlab
results = runtests('Code/jacobi/tests/testJacobi.m');
disp(results);
```

The test suite covers input validation, convergence, single precision, zero right-hand side, and more.

---

## License

This project is licensed under the terms described in the [LICENSE](LICENSE) file.

> **Author:** Thomas Fabbris  
> **Version:** 1.0.0  
> **Release date:** 05-Sep-2025
