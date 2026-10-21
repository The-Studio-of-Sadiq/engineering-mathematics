# Notation and Conventions

[← Table of Contents](table-of-contents.md)

These conventions hold throughout the book unless a chapter states otherwise.

## Numbers and sets

| Symbol | Meaning |
|---|---|
| $\mathbb{N}, \mathbb{Z}, \mathbb{Q}, \mathbb{R}, \mathbb{C}$ | naturals, integers, rationals, reals, complex numbers |
| $j$ | imaginary unit, $j^2=-1$ (engineering convention; mathematics texts use $i$) |
| $z^*$ or $\bar z$ | complex conjugate |
| $\lvert z\rvert,\ \angle z$ | magnitude and argument (phase) of $z$ |
| $\operatorname{Re}\{z\},\ \operatorname{Im}\{z\}$ | real and imaginary parts |

## Scalars, vectors, matrices

| Symbol | Meaning |
|---|---|
| $a,\ x,\ \lambda$ | scalars (italic) |
| $\mathbf{x},\ \mathbf{E}$ | vectors (bold lowercase or field quantities) |
| $\hat{\mathbf{x}},\ \hat{\mathbf{n}}$ | unit vectors |
| $\mathbf{A},\ \mathbf{L}$ | matrices (bold uppercase) |
| $\mathbf{A}^\mathsf{T},\ \mathbf{A}^\mathsf{H}$ | transpose, conjugate (Hermitian) transpose |
| $\mathbf{I}$ | identity matrix |
| $\langle \mathbf{u},\mathbf{v}\rangle$ | inner product |
| $\lVert\mathbf{x}\rVert$ | norm (Euclidean unless subscripted) |

## Coordinate systems

| System | Coordinates | Notes |
|---|---|---|
| Cartesian | $(x, y, z)$ | |
| Cylindrical | $(r, \phi, z)$ | $r$ = distance from $z$-axis; keeps $\rho$ free for charge density |
| Spherical | $(R, \theta, \phi)$ | $\theta$ = polar angle from $+z$, $\phi$ = azimuth (physics/Cheng convention) |

## Calculus and operators

| Symbol | Meaning |
|---|---|
| $\dfrac{dy}{dt},\ \dot y,\ y'$ | derivative (dot = time derivative) |
| $\dfrac{\partial f}{\partial x},\ f_x$ | partial derivative |
| $\nabla f,\ \nabla\cdot\mathbf{F},\ \nabla\times\mathbf{F},\ \nabla^2 f$ | gradient, divergence, curl, Laplacian |
| $\oint_C$ | integral around a closed contour/curve |

## Signals and transforms

| Symbol | Meaning |
|---|---|
| $t$, $n$ | continuous time, discrete-time index |
| $u(t)$, $u[n]$ | unit step |
| $\delta(t)$, $\delta[n]$ | Dirac impulse, Kronecker (unit) impulse |
| $\omega = 2\pi f$ | angular frequency (rad/s), $f$ in Hz |
| $X(s) = \mathcal{L}\{x(t)\}$ | Laplace transform, $s=\sigma+j\omega$ |
| $X(\omega) = \mathcal{F}\{x(t)\}$ | Fourier transform (angular-frequency convention, see Ch. 16) |
| $X(z) = \mathcal{Z}\{x[n]\}$ | Z-transform |
| $x * h$ | convolution |

## Probability

| Symbol | Meaning |
|---|---|
| $P(A)$, $P(A\mid B)$ | probability, conditional probability |
| $X,\ Y$ | random variables (uppercase); $x,\ y$ their values |
| $E[X],\ \mu$ | expectation |
| $\operatorname{Var}(X),\ \sigma^2$ | variance |
| $X\sim\mathcal{N}(\mu,\sigma^2)$ | $X$ is Gaussian |
| $R_{XX}(\tau),\ S_{XX}(\omega)$ | autocorrelation, power spectral density |

## Units and dimensions

SI units throughout. Dimensions are written in brackets: $[V] = \mathrm{M\,L^2\,T^{-3}\,I^{-1}}$. Every equation in this book should be dimensionally consistent — checking this is the cheapest error detector you have (Ch. 1, Ch. 25).

## Callouts

The book uses Obsidian callouts (they degrade to blockquotes on GitHub):

> [!definition] Definition
> A precise statement of a concept.

> [!theorem] Theorem
> A result you can rely on.

> [!example] Worked Example
> A solved problem.

> [!tip] Engineering Note
> Where the method shows up in EE practice.

> [!warning] Common Pitfall
> A mistake worth avoiding.
