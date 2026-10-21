# Chapter 91 — Inverse Problems: Reconstructing the Invisible

[← Previous](90-spectral-theory.md) · [Track J](README.md) · [Table of Contents](../../00-front/table-of-contents.md) · [Next →](92-mathematical-physics.md)

> **Level 3 · Track J — Mathematics at the Frontier**
> **Builds on:** [Ch. 9](../../part-3-linear-algebra/09-vector-spaces-and-linear-transformations.md), [Ch. 10](../../part-3-linear-algebra/10-eigenvalues-and-eigenvectors.md), [Ch. 13](../../part-4-differential-equations/13-partial-differential-equations.md), [Ch. 61](../e-abstract-math/61-functional-analysis.md), [Ch. 90](90-spectral-theory.md) · **Leads to:** [Ch. 92](92-mathematical-physics.md), [Ch. 97](97-measure-theory.md) · **Related:** [Ch. 24](../../part-7-numerical-computational/24-computational-mathematics.md), [Ch. 35](../b-probability-information/35-bayesian-inference.md), [Ch. 76](../g-graphics-simulation/76-signal-processing-audio-visual.md), [Ch. 89](89-pde-beyond-the-basics.md)
> **Reading time:** ~30 min

## 91.1 The hook

On 24 April 1990 the Space Shuttle *Discovery* carried the Hubble Space Telescope into orbit. Its 2.4-metre primary mirror was the most precisely polished large mirror ever made. Within weeks the first images came back blurred. Stars that should have been sharp points were surrounded by a wide, faint halo. The investigation found that the mirror had been ground to the wrong shape: near its edge it was too flat by about $2.2$ micrometres, roughly a fiftieth of the width of a human hair. A testing instrument used during polishing had been assembled with one lens about a millimetre out of place.

The error was tiny, but it was also *perfectly known*. Engineers could compute exactly how the faulty mirror spread the light of a point star, the **point-spread function**. Every Hubble image was therefore the true sky convolved with a known kernel. In principle, undoing the blur is just solving a linear system: the blurred image $\mathbf{g}$ equals $\mathbf{K}\mathbf{f}$, with $\mathbf{K}$ known, so $\mathbf{f}=\mathbf{K}^{-1}\mathbf{g}$. Astronomers spent three years doing exactly this kind of computation, until a servicing mission in December 1993 fitted corrective optics. They learned the hard way that the direct inverse is useless. Reconstructions came out as fields of violent speckle, with every pixel of sensor noise magnified into garbage.

Here is the same experience in one dimension. Take a signal with a box, a ramp and a bump; blur it with a narrow Gaussian; add noise of 1%, smaller than the line width on a plot. The blurring matrix is perfectly invertible. Its exact inverse returns a signal whose error is **273,000 times** the size of the signal itself.

And this is not an astronomy problem. It is the generic situation whenever you measure effects and want causes:

- **CT and MRI**: X-ray attenuations along lines, or Fourier samples, in; a tissue map out.
- **Seismic imaging**: echoes recorded at the surface in; rock layers kilometres down out.
- **Electrical impedance tomography**: voltages on the skin in; the conductivity of the lungs out.
- **Parameter estimation**: a step response in; the time constants of a circuit or a thermal model out ([§25.8](../../part-7-numerical-computational/25-mathematical-modeling.md#258-parameter-identification)).

> **If you know the forward process exactly, why can't you just run it backwards? And if the exact inverse is worthless, what are you allowed to compute instead, and how do you know it isn't fiction?**

## 91.2 The big idea

**Forward processes forget.** Blurring averages neighbouring values together. Diffusion smears heat. A CT projection sums along a line. Physics is mostly smoothing, and smoothing destroys fine detail: two very different signals, one of them with a high-frequency wiggle added, produce almost identical measurements. The forward map is stable because it forgets. That same forgetting makes the inverse unstable: to recover the wiggle, the inverse must amplify the tiny trace it left in the data, and it amplifies the noise sitting in the same place just as much.

In the language of [Ch. 13](../../part-4-differential-equations/13-partial-differential-equations.md#132-boundary-and-initial-conditions), the forward problem is well posed and the inverse problem is **ill posed**: its solution does not depend continuously on the data. Hadamard thought such problems were physically meaningless. A century later, they are most of medical imaging.

**The singular value decomposition shows exactly where the forgetting happens.** Write the forward matrix as $\mathbf{K}=\sum_k\sigma_k\mathbf{u}_k\mathbf{v}_k^\mathsf{T}$ ([§10.7](../../part-3-linear-algebra/10-eigenvalues-and-eigenvectors.md#107-singular-value-decomposition)). Each input pattern $\mathbf{v}_k$ is mapped to output pattern $\mathbf{u}_k$, shrunk by $\sigma_k$. For a blur the $\mathbf{v}_k$ are oscillations of increasing frequency and the $\sigma_k$ fall off a cliff. The exact inverse divides by $\sigma_k$, and in the components where $\sigma_k$ is below the noise level it is dividing noise by almost nothing.

![Picard plot on a log scale: the singular values fall from 1 to below 10 to the minus 8; the data coefficients follow them down until about index 60 and then flatten at a noise floor; the naive solution coefficients, their ratio, grow from index 60 onward to about 10 to the 6](../../images/ch91-picard.svg)

*The Picard plot. Where the data coefficients (orange) still track the singular values (blue), the data carry information about $\mathbf{f}$. Beyond about $k=60$ they hit the noise floor and carry only noise, and dividing by $\sigma_k$ (red) inflates that noise by up to a million.*

**The cure is to want less.** Since some components of $\mathbf{f}$ cannot be recovered, stop trying. Replace "find $\mathbf{f}$ with $\mathbf{K}\mathbf{f}=\mathbf{g}$" by "find a *reasonable* $\mathbf{f}$ that *nearly* fits the data". Reasonable means extra information you supply: the signal is not huge, or is smooth, or has few edges. That extra information is called **regularization**, and the art is in how much of it to add. An analogy: reading a smudged handwritten word. You cannot recover the ink exactly. You recover it by combining what is legible with what words usually look like, and if you trust your expectations too much you read the word you expected instead of the one that was written.

## 91.3 The mathematics

Throughout, the data are $\mathbf{g}=\mathbf{K}\mathbf{f}+\mathbf{e}$ with noise $\lVert\mathbf{e}\rVert\le\delta$, and $\mathbf{K}=\sum_k\sigma_k\mathbf{u}_k\mathbf{v}_k^\mathsf{T}$ with $\sigma_1\ge\sigma_2\ge\cdots>0$. The same formulas hold for compact operators between Hilbert spaces, where $\sigma_k\to0$ is forced ([Ch. 61](../e-abstract-math/61-functional-analysis.md)) and the inverse is unbounded rather than merely large.

### Why the naive inverse explodes

> [!theorem] Noise amplification
> The naive solution is $\mathbf{K}^{-1}\mathbf{g}=\mathbf{f}+\sum_k\dfrac{\mathbf{u}_k^\mathsf{T}\mathbf{e}}{\sigma_k}\mathbf{v}_k$. Its worst-case error over all noise with $\lVert\mathbf{e}\rVert\le\delta$ is $\delta/\sigma_{\min}$. If the noise is white with variance $\eta^2$ per component, the expected squared error is $\eta^2\sum_k\sigma_k^{-2}$.

*Proof.* Since $\{\mathbf{u}_k\}$ is orthonormal, $\mathbf{K}^{-1}=\sum_k\sigma_k^{-1}\mathbf{v}_k\mathbf{u}_k^\mathsf{T}$, and applying it to $\mathbf{K}\mathbf{f}+\mathbf{e}$ gives the formula. The error has coordinates $\mathbf{u}_k^\mathsf{T}\mathbf{e}/\sigma_k$ in the orthonormal basis $\{\mathbf{v}_k\}$, so its squared norm is $\sum_k(\mathbf{u}_k^\mathsf{T}\mathbf{e})^2/\sigma_k^2\le\sigma_{\min}^{-2}\sum_k(\mathbf{u}_k^\mathsf{T}\mathbf{e})^2=\delta^2/\sigma_{\min}^2$, with equality for $\mathbf{e}=\delta\mathbf{u}_{\min}$. For white noise, $E[(\mathbf{u}_k^\mathsf{T}\mathbf{e})^2]=\eta^2$ for each $k$. ∎

In the hook the condition number of $\mathbf{K}$ is $1.8\times10^8$, so 1% noise can become an error of over a million times the signal; the observed $2.7\times10^5$ is the typical case rather than the worst. Making the grid finer does not help. It adds more components with tiny $\sigma_k$, and in the continuum limit $\sigma_{\min}\to0$.

**How bad is bad?** The decay of $\sigma_k$ grades the difficulty. If $\sigma_k\sim k^{-p}$ the problem is **mildly ill posed**: differentiation is $p=1$ (the Volterra operator of [§61.4](../e-abstract-math/61-functional-analysis.md#614-the-algorithm) has $\sigma_k=2/((2k-1)\pi)$), and the CT scanner's Radon transform behaves like $p=\tfrac12$ in two dimensions. If $\sigma_k$ decays exponentially, the problem is **severely ill posed**: a Gaussian blur damps frequency $\omega$ by $e^{-s^2\omega^2/2}$, and running the heat equation backwards for time $T$ requires undoing $e^{-k^2\pi^2T}$. In a severely ill-posed problem, halving the noise buys only a few more recoverable components. That is why CT images are sharp and why nobody reconstructs the temperature history of a rod from its present state.

> [!tip] The Picard condition
> A solution exists in the continuum only if $\sum_k(\mathbf{u}_k^\mathsf{T}\mathbf{g})^2/\sigma_k^2<\infty$: the data coefficients must decay faster than the singular values. Exact data always satisfy it. Noisy data never do, because noise coefficients stop decaying. The Picard plot shows the index where the condition breaks, and that index is roughly how many components any method can hope to recover.

### Tikhonov regularization

Andrey Tikhonov's 1963 proposal: instead of solving $\mathbf{K}\mathbf{f}=\mathbf{g}$, minimize

$$
J_\lambda(\mathbf{f})=\lVert\mathbf{K}\mathbf{f}-\mathbf{g}\rVert^2+\lambda\lVert\mathbf{f}\rVert^2,\qquad\lambda>0 .
$$

The first term asks for data fit, the second penalizes size, and $\lambda$ sets the exchange rate. Setting the gradient to zero gives the regularized normal equations $(\mathbf{K}^\mathsf{T}\mathbf{K}+\lambda\mathbf{I})\mathbf{f}_\lambda=\mathbf{K}^\mathsf{T}\mathbf{g}$, which always have a unique solution because $\mathbf{K}^\mathsf{T}\mathbf{K}+\lambda\mathbf{I}$ is positive definite. It is ordinary least squares ([§9.11](../../part-3-linear-algebra/09-vector-spaces-and-linear-transformations.md#911-least-squares-approximation)) for the stacked system $\begin{bmatrix}\mathbf{K}\\ \sqrt\lambda\,\mathbf{I}\end{bmatrix}\mathbf{f}\approx\begin{bmatrix}\mathbf{g}\\ \mathbf{0}\end{bmatrix}$. In the SVD basis, $\mathbf{K}^\mathsf{T}\mathbf{K}+\lambda\mathbf{I}$ is diagonal with entries $\sigma_k^2+\lambda$, so

$$
\mathbf{f}_\lambda=\sum_k\underbrace{\frac{\sigma_k^2}{\sigma_k^2+\lambda}}_{\text{filter factor}}\;\frac{\mathbf{u}_k^\mathsf{T}\mathbf{g}}{\sigma_k}\,\mathbf{v}_k .
$$

Where $\sigma_k^2\gg\lambda$ the filter factor is nearly $1$ and the component is kept. Where $\sigma_k^2\ll\lambda$ it is nearly $0$ and the component is suppressed. Tikhonov is a smooth version of truncating the SVD, which [Ch. 61](../e-abstract-math/61-functional-analysis.md) tried, with the cutoff at $\sigma_k\approx\sqrt\lambda$.

![Three panels: the true signal and its slightly noisy blurred measurement, almost indistinguishable; the naive inverse, a violent oscillation of amplitude 2 times 10 to the 5; and Tikhonov reconstructions at lambda = 1e-6, which is still noisy, lambda = 5e-3, which tracks the box, tent and bump with small ripples, and lambda = 0.3, which is smooth but too low and too rounded](../../images/ch91-deblurring.svg)

*Deblurring a signal. The measurement (left) looks like the truth, yet the exact inverse (centre) is pure amplified noise. Tikhonov (right) with $\lambda$ too small is still noise-dominated (relative error $116\%$); with $\lambda=5\times10^{-3}$ the error is $9.8\%$; with $\lambda=0.3$ the noise is gone but so are the corners and the peak heights (error $27\%$).*

### How good can it get?

The regularized solution has two errors with opposite appetites for $\lambda$, and the theorem below makes their trade exact. It needs one assumption about the truth, because without some smoothness of $\mathbf{f}$ no method can do better than arbitrarily badly.

> [!theorem] Tikhonov error bound
> Suppose $\mathbf{f}=\mathbf{K}^\mathsf{T}\mathbf{w}$ for some $\mathbf{w}$ (a **source condition**: $\mathbf{f}$ is smoother than a generic vector, in the sense the forward map defines). Then
> $$
> \lVert\mathbf{f}_\lambda-\mathbf{f}\rVert\le\frac{\delta}{2\sqrt\lambda}+\frac{\sqrt\lambda}{2}\lVert\mathbf{w}\rVert ,
> $$
> and choosing $\lambda=\delta/\lVert\mathbf{w}\rVert$ gives $\lVert\mathbf{f}_\lambda-\mathbf{f}\rVert\le\sqrt{\delta\,\lVert\mathbf{w}\rVert}$.

*Proof.* Insert $\mathbf{g}=\mathbf{K}\mathbf{f}+\mathbf{e}$ into the filtered formula and use $\mathbf{u}_k^\mathsf{T}\mathbf{K}\mathbf{f}=\sigma_k\mathbf{v}_k^\mathsf{T}\mathbf{f}$:

$$
\mathbf{f}_\lambda-\mathbf{f}=\sum_k\frac{\sigma_k}{\sigma_k^2+\lambda}(\mathbf{u}_k^\mathsf{T}\mathbf{e})\,\mathbf{v}_k\;-\;\sum_k\frac{\lambda}{\sigma_k^2+\lambda}(\mathbf{v}_k^\mathsf{T}\mathbf{f})\,\mathbf{v}_k .
$$

*Noise term.* By the AM–GM inequality, $\sigma^2+\lambda\ge2\sigma\sqrt\lambda$, so every coefficient $\sigma_k/(\sigma_k^2+\lambda)\le1/(2\sqrt\lambda)$. The first sum therefore has norm at most $\lVert\mathbf{e}\rVert/(2\sqrt\lambda)\le\delta/(2\sqrt\lambda)$.

*Bias term.* The source condition gives $\mathbf{v}_k^\mathsf{T}\mathbf{f}=\mathbf{v}_k^\mathsf{T}\mathbf{K}^\mathsf{T}\mathbf{w}=\sigma_k\,\mathbf{u}_k^\mathsf{T}\mathbf{w}$. So the second sum has coefficients $\lambda\sigma_k(\mathbf{u}_k^\mathsf{T}\mathbf{w})/(\sigma_k^2+\lambda)$, and by the same inequality $\lambda\sigma/(\sigma^2+\lambda)\le\sqrt\lambda/2$. Its norm is at most $\tfrac{\sqrt\lambda}{2}\lVert\mathbf{w}\rVert$.

*Balance.* The sum $a/\sqrt\lambda+b\sqrt\lambda$ with $a=\delta/2$, $b=\lVert\mathbf{w}\rVert/2$ is minimized at $\lambda=a/b=\delta/\lVert\mathbf{w}\rVert$, where it equals $2\sqrt{ab}=\sqrt{\delta\lVert\mathbf{w}\rVert}$. ∎

Read the result carefully, because it is the central fact of the subject. With exact data ($\delta\to0$) and $\lambda$ chosen proportionally to $\delta$, the reconstruction converges to the truth, so the method is a genuine **regularization**. But the error goes like $\sqrt\delta$, not $\delta$: to gain one more correct digit you need two more digits of data accuracy. Stronger smoothness assumptions improve the exponent up to $\delta^{2/3}$ and no further for Tikhonov (its **saturation**). Ill-posedness is never fully cured; it is only negotiated.

### Choosing $\lambda$ from the data

The optimal $\lambda$ in the theorem depends on $\lVert\mathbf{w}\rVert$, which describes the unknown truth. In practice $\lambda$ must come from what you have.

- **The discrepancy principle** (Morozov, 1966). Do not fit the data better than the noise deserves: pick $\lambda$ so that $\lVert\mathbf{K}\mathbf{f}_\lambda-\mathbf{g}\rVert=\tau\delta$ with $\tau\ge1$. It needs the noise level $\delta$, and with it, it provably achieves the optimal rate $O(\sqrt\delta)$ under the source condition. In the hook's example it picks $\lambda\approx5.6\times10^{-3}$ and an error of $9.8\%$, essentially the best possible.
- **The L-curve** (Hansen, 1992). Plot the solution size against the residual on log–log axes for all $\lambda$. Small $\lambda$ gives large solutions with small residuals; large $\lambda$, the reverse. The curve has a corner where neither can be reduced without a large increase in the other, and the corner is the choice. It needs no knowledge of $\delta$.
- **Generalized cross-validation** chooses $\lambda$ so that the fit predicts each data point well when that point is left out. It also ignores $\delta$.

![Left: the L-curve on log-log axes, nearly vertical for small lambda and nearly horizontal for large lambda, with the L-curve corner, the discrepancy-principle point and the best possible point marked. Right: relative error against lambda, a valley with its minimum near 4e-3; the discrepancy choice sits at the bottom, the L-curve corner partway up the noisy side](../../images/ch91-l-curve.svg)

*Three ways of choosing $\lambda$ on the same problem. The discrepancy principle, which knows $\delta$, lands at the bottom of the error valley. The L-curve corner, which doesn't, picks $\lambda\approx5\times10^{-5}$, about 100 times too small, and more than doubles the error to $23\%$.*

Is there a rule that never needs $\delta$ and always works? No, and the reason is a theorem. **Bakushinskii's veto** (1984): for an ill-posed problem, any regularization method whose parameter is chosen from the data $\mathbf{g}$ alone, without knowledge of $\delta$, fails to converge for some true solutions as $\delta\to0$. The noise level is information that cannot be manufactured from the data. Heuristic rules like the L-curve and GCV still work very well in most practical cases, and understanding exactly when is an active line of research.

### The Bayesian reading

Model the unknown as random, with a Gaussian prior $\mathbf{f}\sim\mathcal{N}(\mathbf{0},\tau^2\mathbf{I})$, and the noise as $\mathbf{e}\sim\mathcal{N}(\mathbf{0},\eta^2\mathbf{I})$. By Bayes' rule ([Ch. 35](../b-probability-information/35-bayesian-inference.md)) the posterior density is proportional to

$$
\exp\Big(-\frac{\lVert\mathbf{K}\mathbf{f}-\mathbf{g}\rVert^2}{2\eta^2}-\frac{\lVert\mathbf{f}\rVert^2}{2\tau^2}\Big),
$$

and its maximizer is exactly the Tikhonov solution with $\lambda=\eta^2/\tau^2$. The regularization parameter is a ratio of noise variance to signal variance, and the penalty $\lVert\mathbf{f}\rVert^2$ is a belief about what signals look like. Changing the penalty changes the belief: $\lVert\mathbf{D}\mathbf{f}\rVert^2$ with a difference matrix $\mathbf{D}$ says "smooth", the total variation $\lVert\mathbf{D}\mathbf{f}\rVert_1$ says "piecewise constant, with sharp edges", and $\lVert\mathbf{f}\rVert_1$ in a wavelet basis says "sparse", the idea behind compressed sensing ([Ch. 76](../g-graphics-simulation/76-signal-processing-audio-visual.md)).

The Bayesian view gives something Tikhonov alone does not: an uncertainty. The posterior covariance $(\mathbf{K}^\mathsf{T}\mathbf{K}/\eta^2+\mathbf{I}/\tau^2)^{-1}$ says which features of the reconstruction are determined by the data and which are only the prior talking. But it also creates the tension [Ch. 61](../e-abstract-math/61-functional-analysis.md) warned about. In the Bayesian view $\lambda$ is not tuned; it is fixed by $\eta$ and $\tau$. And the prior $\mathcal{N}(\mathbf{0},\tau^2\mathbf{I})$ says each pixel is independent, so as the grid is refined, a "typical" signal under the prior becomes white noise, nothing like the box and bump we were reconstructing.

> [!mystery] A prior that survives refinement
> A good prior should describe the same beliefs whether the signal is stored on 200 points or 200 million. That requires a probability measure on a space of *functions*, and on an infinite-dimensional space there is no uniform "Lebesgue" measure to write densities against. Gaussian measures with carefully decaying covariances work. In 2004 Lassas and Siltanen showed that the natural discretized total-variation prior, the edge-preserving favourite of image processing, does not: as the grid is refined, its reconstructions converge either to something blurry or to something meaningless, depending on how the prior is scaled. Which priors are "discretization invariant", and why, is a question about measures on function spaces: [Ch. 97](97-measure-theory.md).

### Nonlinear inverse problems and parameter estimation

Most real inverse problems are nonlinear: the data $\mathbf{g}=F(\mathbf{p})+\mathbf{e}$ depend nonlinearly on parameters $\mathbf{p}$, such as the conductivity $k(x)$ of the rod in [Ch. 89](89-pde-beyond-the-basics.md) or the time constants of a model. The workhorse is to linearize and regularize at each step. Gauss–Newton ([§24.7](../../part-7-numerical-computational/24-computational-mathematics.md#247-parameter-estimation)) solves $\mathbf{J}\,\Delta\mathbf{p}\approx\mathbf{g}-F(\mathbf{p})$ with the Jacobian $\mathbf{J}$, and the Levenberg–Marquardt modification adds $\mu\lVert\Delta\mathbf{p}\rVert^2$: it is Tikhonov regularization applied to every linearized step. A new difficulty appears, **identifiability**: different parameters can produce *identical* data, not just similar data, and then no amount of regularization can tell them apart.

The model problem of the field is **Calderón's problem** (1980): recover the conductivity $\gamma(\mathbf{x})$ inside a body from all possible pairs of boundary voltages and boundary currents, the Dirichlet-to-Neumann map. Uniqueness was proved by Sylvester and Uhlmann in 1987 for smooth conductivities in three or more dimensions, and by Astala and Päivärinta in 2006 for any bounded conductivity in the plane. Stability is another matter: Mandache proved in 2001 that it can be at best *logarithmic*, so that improving the data accuracy from $\delta$ to $\delta^2$ may only double the reconstruction accuracy. Impedance tomography is exponentially ill posed, however clever the algorithm.

> [!mystery] How rough can a conductivity be and still be heard?
> In three dimensions uniqueness for Calderón's problem is known for Lipschitz conductivities (Caro and Rogers, 2016) and somewhat beyond, but not for merely bounded ones, which is what real tissue boundaries and composite materials look like. Whether two different bounded conductivities in 3-D can produce exactly the same boundary measurements is an open problem. So is a satisfying answer to the practical question: given logarithmic stability, how much prior information buys how much resolution in a real impedance scanner?

## 91.4 The algorithm

For a moderate problem, compute the SVD once, at $O(n^3)$ cost. After that, $\mathbf{f}_\lambda$ for any $\lambda$ costs only $O(n^2)$, so scanning a hundred values of $\lambda$ is cheap. That is what the experiment below does. Without an SVD, the stable way to compute one Tikhonov solution is the stacked least-squares problem, solved by QR rather than by forming $\mathbf{K}^\mathsf{T}\mathbf{K}$, which would square the condition number ([§24.3](../../part-7-numerical-computational/24-computational-mathematics.md#243-numerical-linear-algebra)). The discrepancy principle then becomes a one-dimensional root-finding problem in $\log\lambda$, because the residual $\lVert\mathbf{K}\mathbf{f}_\lambda-\mathbf{g}\rVert$ increases monotonically with $\lambda$.

```python
import numpy as np

def tikhonov(K, g, lam):
    """argmin ||K f - g||^2 + lam ||f||^2, via a stacked least-squares problem (QR inside)."""
    n = K.shape[1]
    A = np.vstack([K, np.sqrt(lam) * np.eye(n)])
    b = np.concatenate([g, np.zeros(n)])
    return np.linalg.lstsq(A, b, rcond=None)[0]

def discrepancy(K, g, delta, tau=1.0, lo=1e-12, hi=1e2, steps=40):
    """Largest lam with ||K f_lam - g|| <= tau * delta, by bisection in log(lam)."""
    for _ in range(steps):
        mid = np.sqrt(lo * hi)
        if np.linalg.norm(K @ tikhonov(K, g, mid) - g) > tau * delta:
            hi = mid                      # residual too big: regularizing too much
        else:
            lo = mid
    return lo
```

On the hook's problem, `discrepancy` returns $\lambda=5.6\times10^{-3}$ and a reconstruction with relative error $9.8\%$, matching the SVD scan.

**Large problems.** A CT image of $512\times512$ pixels has $n\approx2.6\times10^5$ unknowns, and $\mathbf{K}$ is never stored as a dense matrix: it is a routine that computes projections. Then the method of choice is an iterative least-squares solver such as CGLS or LSQR, which needs only products with $\mathbf{K}$ and $\mathbf{K}^\mathsf{T}$. These methods regularize by themselves. Early iterations pick up the large-$\sigma_k$ components; later ones reach the small-$\sigma_k$, noise-dominated components. So the error first falls and then rises again, which is called **semi-convergence**, and the iteration count plays the role of $1/\lambda$. Stopping by the discrepancy principle is standard.

**Failure modes.**

- **The inverse crime.** If you generate synthetic data with the same discretized $\mathbf{K}$ that you then invert, the model is perfect and the reconstruction looks far better than it will on real data. Honest tests generate data on a finer grid, or with a different method, than the one used for inversion. The experiment below commits this crime on purpose; one of its knobs removes it.
- **Wrong noise model.** Astronomical images have Poisson photon noise, not Gaussian noise. The least-squares misfit is then the wrong measure, and the Richardson–Lucy iteration, which maximizes the Poisson likelihood and keeps intensities positive, was the tool astronomers actually used on Hubble's flawed images.
- **Wrong prior.** The penalty $\lVert\mathbf{f}\rVert^2$ rounds every edge, as the $\lambda=0.3$ curve shows. Edge-preserving penalties such as total variation are better, but they are non-smooth and need convex-optimization algorithms ([Ch. 47](../d-optimization/47-convexity.md)).
- **Boundary effects.** Near the ends of the interval, the blur averages over values the model assumes are zero. If the true scene continues beyond the edge, the reconstruction invents artefacts there.

## 91.5 The experiment

**See the filter.** The blue curve is Tikhonov's filter factor $\sigma^2/(\sigma^2+\lambda)$ with $\lambda=0.01$, plotted against the singular value $\sigma$; the red step is truncated SVD with cutoff $\sigma=\sqrt\lambda=0.1$. The green curve, $\sigma/(\sigma^2+\lambda)$, is how much Tikhonov amplifies noise in a component with singular value $\sigma$, and the dashed line is the bound $1/(2\sqrt\lambda)=5$ from §91.3. Change `0.01` to `0.0001` everywhere (and the cutoff `0.1` to `0.01`, the bound `5` to `50`) and watch the filter sharpen while the noise amplification grows tenfold.

```desmos-graph
left=0; right=1; bottom=-0.1; top=6
---
y=\frac{x^{2}}{x^{2}+0.01}|blue|0<x<1
y=1|red|0.1<x<1
y=0|red|0<x<0.1
y=\frac{x}{x^{2}+0.01}|green|0<x<1
y=5|grey|dashed
```

*Tikhonov's filter (blue) keeps components with $\sigma\gg\sqrt\lambda$ and smoothly suppresses the rest, unlike the hard cut of truncated SVD (red). The amplification of noise (green) peaks exactly at $\sigma=\sqrt\lambda$, at the height the error bound uses.*

**Deblur it yourself.** This script rebuilds the hook's example, prints the error for a range of $\lambda$, and lets the discrepancy principle choose.

```python
import numpy as np

n, width, noise_level = 200, 0.01, 0.01   # knobs: grid, blur width, relative noise
rng = np.random.default_rng(0)
x = (np.arange(n) + 0.5) / n
K = np.exp(-(x[:, None] - x[None, :])**2 / (2 * width**2)) / (np.sqrt(2 * np.pi) * width * n)
f = ((x > 0.15) & (x < 0.35)) + 0.8 * np.maximum(0, 1 - np.abs(x - 0.6) / 0.1) \
    + 0.6 * np.exp(-((x - 0.85) / 0.04)**2)          # box, tent, bump
g_clean = K @ f
e = noise_level * np.linalg.norm(g_clean) / np.sqrt(n) * rng.standard_normal(n)
g = g_clean + e
delta = np.linalg.norm(e)

U, s, Vt = np.linalg.svd(K)
print(f"condition number {s[0] / s[-1]:.2e}")
rel = lambda fr: np.linalg.norm(fr - f) / np.linalg.norm(f)
print(f"naive inverse: relative error {rel(np.linalg.solve(K, g)):.2e}")

def tik(lam):                                       # Tikhonov through the SVD filter
    return Vt.T @ (s / (s**2 + lam) * (U.T @ g))

print("   lambda    residual/delta   rel. error")
for lam in [1e-8, 1e-6, 1e-4, 1e-3, 5e-3, 1e-2, 1e-1, 1]:
    fl = tik(lam)
    print(f"{lam:9.0e}   {np.linalg.norm(K @ fl - g) / delta:10.2f}   {rel(fl):10.3f}")

lams = np.logspace(-10, 1, 400)                     # discrepancy principle
res = np.array([np.linalg.norm(K @ tik(l) - g) for l in lams])
lam_d = lams[res <= delta].max()
print(f"discrepancy principle picks lambda = {lam_d:.1e}, error {rel(tik(lam_d)):.3f}")
```

With the default knobs the condition number is $1.83\times10^8$ and the naive error is $2.73\times10^5$. The table's error column runs $12.8$, $1.16$, $0.178$, $0.106$, $0.098$, $0.100$, $0.149$, $0.522$ as $\lambda$ grows: a valley with a flat floor between about $10^{-3}$ and $10^{-2}$. The discrepancy principle picks $\lambda=5.2\times10^{-3}$, with error $0.098$. Knobs to turn:

- Set `noise_level` to $10^{-3}$ and then $10^{-4}$. By how much does the best error improve each time? Compare with the $\sqrt\delta$ rate of the theorem: is this signal smooth enough to satisfy the source condition? (The box has jumps.)
- Double `width`. The problem becomes far more severely ill posed. How many singular values are now above the noise floor?
- Commit no crime: build the data with a finer kernel matrix (say $3n$ points, then average each group of three), and invert with the $n$-point `K`. How much worse is the best achievable error?
- Replace the penalty $\lVert\mathbf{f}\rVert^2$ by $\lVert\mathbf{D}\mathbf{f}\rVert^2$ with `D = np.diff(np.eye(n), axis=0)`, solving the normal equations $(\mathbf{K}^\mathsf{T}\mathbf{K}+\lambda\mathbf{D}^\mathsf{T}\mathbf{D})\mathbf{f}=\mathbf{K}^\mathsf{T}\mathbf{g}$ directly. Which features improve, and which get worse?

## 91.6 The challenge

1. **Intuition.** A thermometer at the midpoint of a rod records the temperature for one hour. Explain without equations why estimating the rod's *initial* temperature profile from this record is severely ill posed, while estimating a single unknown heat-loss coefficient from it is well posed. What is different about the number of unknowns compared with the number of components the data can carry?
2. **Derivation.** For truncated SVD with cutoff index $r$, $\mathbf{f}_r=\sum_{k\le r}\sigma_k^{-1}(\mathbf{u}_k^\mathsf{T}\mathbf{g})\mathbf{v}_k$, derive the analogue of the Tikhonov error bound under the source condition $\mathbf{f}=\mathbf{K}^\mathsf{T}\mathbf{w}$, and show that the best choice of $r$ again gives an error of order $\sqrt{\delta\lVert\mathbf{w}\rVert}$. Then show that under the stronger condition $\mathbf{f}=(\mathbf{K}^\mathsf{T}\mathbf{K})\mathbf{z}$, truncated SVD reaches $O(\delta^{2/3})$ and keeps improving under still stronger conditions, while Tikhonov stalls at $\delta^{2/3}$. Find the step of the Tikhonov argument that causes the stall.
3. **Implementation.** Implement total-variation deblurring for the hook's signal: minimize $\lVert\mathbf{K}\mathbf{f}-\mathbf{g}\rVert^2+\lambda\sum_i\sqrt{(f_{i+1}-f_i)^2+\epsilon^2}$ by gradient descent or by lagged diffusivity (re-solve a weighted least-squares problem at each step). Compare its reconstruction of the box with Tikhonov's at the discrepancy-principle $\lambda$, and report both errors.
4. **Open investigation.** Train the simplest possible learned reconstructor: a linear map $\mathbf{f}\approx\mathbf{B}\mathbf{g}$ fitted by least squares on a thousand random training signals made of boxes and bumps. Compare it with Tikhonov on new signals of the same kind, then on signals of a different kind (a sine wave, a single spike). Then search for a small perturbation of the data, $\lVert\Delta\mathbf{g}\rVert\le\delta$, that maximizes the change in $\mathbf{B}\mathbf{g}$. Is the learned map more or less stable than Tikhonov, and does it trade stability for accuracy in any predictable way?

## 91.7 Loose threads

An inverse problem is a forward problem that forgot something. The SVD shows exactly which components were forgotten, regularization fills them in with assumptions, and the Tikhonov bound shows the price: errors that shrink like $\sqrt\delta$ at best, and a parameter whose choice needs information, the noise level, that the data cannot supply. The Bayesian reading turns the assumptions into explicit probability and buys honest error bars. Calderón's problem and the priors that survive refinement are still open. One more clue comes from the present day.

> [!mystery] Reconstructions that hallucinate
> Since about 2017, neural networks trained on large image collections have produced MRI and CT reconstructions sharper than any regularization formula, from fewer measurements. In 2020 Antun, Renna, Poon, Adcock and Hansen showed that several such networks are unstable in a way Tikhonov never is. A tiny, carefully chosen perturbation of the measurements can make a lesion vanish or a false structure appear, and adding *more* measurements can make the result worse. Their follow-up work showed that, for some inverse problems, stable and accurate neural networks exist but no training algorithm can be guaranteed to find them. When a learned reconstructor can be trusted, and how to certify it on a particular scan, is an open problem with patients waiting on the answer.

Further pointers:

- **The physics behind the operators.** Most forward maps in this chapter, from blur to radiation transport to wave scattering, come from variational principles and conservation laws, and the symmetries of the physics decide which parameters can be identified: [Ch. 92](92-mathematical-physics.md).
- **Uncertainty at scale.** Turning a posterior over a million unknowns into error bars needs sampling methods that work in high dimension, such as the preconditioned Crank–Nicolson MCMC, whose design comes straight from the infinite-dimensional view: [Ch. 38](../b-probability-information/38-monte-carlo.md).
- **Geometry from images.** Recovering depth from photographs, the inverse problem left open in [Ch. 69](../g-graphics-simulation/69-mathematics-of-3d-worlds.md), is nonlinear, and its ill-posedness shows up as near-degenerate camera configurations.

## Further reading

- P. C. Hansen, *Discrete Inverse Problems: Insight and Algorithms*, SIAM, 2010. The Picard plot, the L-curve and Tikhonov regularization in practice, with MATLAB code throughout.
- J. L. Mueller and S. Siltanen, *Linear and Nonlinear Inverse Problems with Practical Applications*, SIAM, 2012. From deblurring and CT to electrical impedance tomography, with real data.
- H. W. Engl, M. Hanke and A. Neubauer, *Regularization of Inverse Problems*, Kluwer, 1996. The rigorous theory: source conditions, convergence rates, saturation and parameter choice.
- J. Kaipio and E. Somersalo, *Statistical and Computational Inverse Problems*, Springer, 2005. The Bayesian approach, including priors, error modelling and the inverse crime.
- A. M. Stuart, "Inverse problems: a Bayesian perspective", *Acta Numerica* 19 (2010), 451–559. Bayesian inversion formulated directly on function spaces.
- V. Antun, F. Renna, C. Poon, B. Adcock and A. C. Hansen, "On instabilities of deep learning in image reconstruction and the potential costs of AI", *PNAS* 117 (2020), 30088–30095. The paper behind the last mystery.
