# Chapter 3 — Trigonometry and Periodic Functions

[← Previous](02-complex-numbers.md) · [Table of Contents](../00-front/table-of-contents.md) · [Next →](../part-2-calculus/04-differential-calculus.md)

> **Prerequisites:** [Ch. 1](01-mathematical-language-and-algebra.md) (functions, graphs, basic trigonometric ratios), [Ch. 2](02-complex-numbers.md) (Euler's formula, phasors)
> **Learning outcomes:** After this chapter you can
> - measure angles in radians, degrees and revolutions, convert between them, and explain why calculus insists on radians;
> - read amplitude, offset, frequency and phase directly from a sinusoidal expression or graph, and convert any $a\cos\omega t+b\sin\omega t$ into amplitude–phase form;
> - compute phase differences, time delays and leads/lags between sinusoids, and add same-frequency sinusoids exactly;
> - decide whether a sum of sinusoids (continuous or discrete-time) is periodic and find its fundamental period;
> - use the compound-angle, power-reduction and product-to-sum identities to move between products, powers and harmonics;
> - describe periodic waveforms (square, triangle, sawtooth, rectified sine) by formula, symmetry, mean and RMS, and explain how they are built from harmonics.

Almost every quantity an electrical engineer measures either *is* a sinusoid or is analysed as a sum of sinusoids. The mains supply is a 50 Hz or 60 Hz sinusoid; a radio carrier is a sinusoid whose amplitude, frequency or phase is nudged to carry information; the steady-state response of any linear circuit to a sinusoid is another sinusoid of the same frequency; and a distorted current drawn by a rectifier is described by the list of sinusoidal *harmonics* it contains. Trigonometry is the language of all of this.

This chapter builds the working fluency that later chapters assume. It is not a chapter on triangle-solving: the emphasis is on the circular functions as *functions of time*, on the three numbers (amplitude, frequency, phase) that describe a sinusoid, on the identities that convert products and powers of sinusoids into sums of sinusoids, and on the idea that a periodic waveform is a superposition of harmonically related sinusoids. The full theory of that last idea — Fourier series, convergence and spectra — is developed in [Ch. 15](../part-5-transforms/15-fourier-series.md); here we build the intuition and the algebra that make Ch. 15 easy to read.

Throughout, $t$ is time in seconds, $\omega$ is angular frequency in rad/s, $f$ is frequency in Hz, and $j$ is the imaginary unit.

## 3.1 Angular measurement

### Why angle measure needs care

An angle measures *rotation*: how far a ray has turned about a point. In electrical engineering, rotation appears literally (the shaft of a motor or generator, the rotor position read by an encoder, the bearing of an antenna) and figuratively (the phase of a sinusoid, which is the angle of a rotating phasor). Getting the *unit* of angle right is the single most common source of numerical error in trigonometric calculations, because the same number — say $1.5$ — is a modest angle in radians, a tiny angle in degrees and an enormous one in revolutions, and a calculator or software library will silently accept any of them.

### Units of angle

> [!definition] Radian
> Let a circle of radius $r$ be centred at the vertex of an angle, and let $s$ be the length of arc the angle cuts off. The angle in **radians** is
> $$
> \theta = \frac{s}{r}.
> $$
> Because it is a ratio of two lengths, the radian is dimensionless: $[\theta]=1$. One radian is the angle that cuts off an arc equal to the radius.

A full turn cuts off the whole circumference $2\pi r$, so one revolution is $2\pi$ rad. All other units are defined by dividing the full turn into a conventional number of parts:

| Unit | Full turn | Right angle | Typical use |
|---|---|---|---|
| radian (rad) | $2\pi\approx 6.2832$ | $\pi/2$ | all of calculus, $\omega t$, software libraries |
| degree ($^\circ$) | $360$ | $90$ | phase on datasheets and phasor diagrams |
| revolution (rev, turn) | $1$ | $1/4$ | shaft speed (rpm), encoder counts |
| gradian (gon) | $400$ | $100$ | surveying; a trap on some calculators |
| milliradian (mrad) | $\approx 6283$ | $\approx 1571$ | small pointing errors, optics |

The conversion rules follow from equating full turns:

$$
\theta_{\text{deg}} = \frac{180}{\pi}\,\theta_{\text{rad}},\qquad
\theta_{\text{rad}} = \frac{\pi}{180}\,\theta_{\text{deg}},\qquad
\theta_{\text{rev}} = \frac{\theta_{\text{rad}}}{2\pi}.
$$

A useful anchor is $1\ \text{rad}\approx 57.296^\circ$ and $1^\circ\approx 0.017453\ \text{rad}$.

> [!example] Worked Example 3.1 — Unit conversions
> (a) Convert $2.5$ rad to degrees. (b) Convert $37^\circ$ to radians. (c) A shaft turns through $1350^\circ$. How many revolutions is that, and what angle in $[0^\circ,360^\circ)$ does it end at?
>
> **(a)** $\theta = 2.5\times 180/\pi = 143.24^\circ$.
>
> **(b)** $\theta = 37\times\pi/180 = 0.6458$ rad.
>
> **(c)** $1350/360 = 3.75$ rev. Three whole turns ($1080^\circ$) bring the shaft back to where it started, leaving $1350-1080=270^\circ$. The final position is $270^\circ$ (equivalently $3\pi/2$ rad, or $-90^\circ$).
>
> **Check:** $143.24^\circ/57.296 = 2.500$; $0.6458\times 57.296=37.00^\circ$.

### Arc length, sector area and angular speed

The radian definition gives three formulas with no conversion factors — which is precisely why radians are the natural unit:

$$
s = r\theta,\qquad A_{\text{sector}} = \tfrac12 r^2\theta,\qquad v = r\,\Omega,
$$

where $\Omega = d\theta/dt$ is the angular speed in rad/s and $v$ the linear speed of a point at radius $r$. The sector area follows by proportion: the sector is the fraction $\theta/(2\pi)$ of the disc of area $\pi r^2$, so $A=\tfrac{\theta}{2\pi}\pi r^2=\tfrac12 r^2\theta$. If $\theta$ were in degrees each formula would acquire a factor $\pi/180$.

Shaft speeds are usually quoted in revolutions per minute (rpm). Since one revolution is $2\pi$ rad and one minute is $60$ s,

$$
\Omega\ [\text{rad/s}] = \frac{2\pi}{60}\,N\ [\text{rpm}] \approx 0.10472\,N.
$$

> [!example] Worked Example 3.2 — Encoder resolution and tip speed
> An incremental encoder on a motor shaft gives 2048 counts per revolution. The shaft drives a pulley of radius $r=40$ mm and runs at $N=3000$ rpm.
>
> (a) What angle does one count represent, in degrees and in milliradians? (b) What arc length of belt passes per count? (c) What is the belt speed?
>
> **(a)** One count is $1/2048$ rev, i.e. $360/2048 = 0.17578^\circ$, or $2\pi/2048 = 3.0680\times10^{-3}$ rad $=3.068$ mrad.
>
> **(b)** $s=r\theta = 0.040\times 3.0680\times10^{-3} = 1.227\times10^{-4}$ m $=0.1227$ mm.
>
> **(c)** $\Omega = 2\pi\times 3000/60 = 314.16$ rad/s, so $v=r\Omega = 0.040\times314.16 = 12.57$ m/s.
>
> **Check:** At 3000 rpm the shaft makes 50 rev/s; the circumference is $2\pi(0.04)=0.2513$ m; $50\times0.2513=12.57$ m/s. ✓

### Electrical and mechanical degrees

In a rotating machine with $p$ pole *pairs*, the magnetic field pattern repeats $p$ times around the air gap. One mechanical revolution therefore carries the generated voltage through $p$ complete electrical cycles. Engineers distinguish

$$
\theta_e = p\,\theta_m,\qquad \omega_e = p\,\Omega_m,
$$

where the subscript $e$ denotes *electrical* angle/frequency and $m$ *mechanical*. A four-pole machine ($p=2$) running at 1500 rpm has $\Omega_m = 157.08$ rad/s and produces $\omega_e = 314.16$ rad/s, i.e. $f=50$ Hz. Phase angles of voltages and currents are always electrical angles.

> [!tip] Engineering Note — Which "degree" is on the datasheet?
> Motor-control datasheets quote commutation angles in electrical degrees; encoder datasheets quote resolution in mechanical counts. Mixing them up in a field-oriented controller for a $p$-pole-pair motor produces a commutation error multiplied by $p$, which shows up as a motor that vibrates, runs hot, or will not start.

### The circular functions on the unit circle

The six trigonometric functions are defined for *any* real angle, not only acute angles of right triangles, by placing the angle in standard position (vertex at the origin, initial ray along the positive $x$-axis, positive angles anticlockwise) and reading off the point where the terminal ray meets the unit circle.

> [!definition] Circular functions
> Let the terminal ray of the angle $\theta$ (in radians) meet the unit circle $x^2+y^2=1$ at $P=(x,y)$. Then
> $$
> \cos\theta = x,\qquad \sin\theta = y,\qquad \tan\theta=\frac{y}{x}\ (x\neq0),
> $$
> $$
> \sec\theta=\frac1x,\qquad \csc\theta=\frac1y,\qquad \cot\theta=\frac{x}{y}.
> $$
> For a point at distance $r$ from the origin, $(x,y)=(r\cos\theta,\ r\sin\theta)$.

Several facts follow at once from the picture:

- **Range.** $-1\le\cos\theta\le1$ and $-1\le\sin\theta\le1$, because $P$ lies on the unit circle.
- **Pythagorean identity.** $\cos^2\theta+\sin^2\theta=1$, because $x^2+y^2=1$.
- **Periodicity.** Adding $2\pi$ to $\theta$ returns $P$ to the same point, so $\sin(\theta+2\pi)=\sin\theta$ and $\cos(\theta+2\pi)=\cos\theta$. The tangent repeats after only $\pi$, because $\theta+\pi$ sends $P$ to $-P$ and $(-y)/(-x)=y/x$.
- **Symmetry.** Reflecting in the $x$-axis replaces $\theta$ by $-\theta$ and $y$ by $-y$: $\cos(-\theta)=\cos\theta$ (cosine is *even*), $\sin(-\theta)=-\sin\theta$ (sine is *odd*).
- **Complementarity.** Reflecting in the line $y=x$ swaps $x$ and $y$ and sends $\theta$ to $\pi/2-\theta$: $\sin(\pi/2-\theta)=\cos\theta$, $\cos(\pi/2-\theta)=\sin\theta$.

The signs of $\sin$ and $\cos$ in each quadrant follow from the signs of $y$ and $x$:

| Quadrant | $\theta$ range | $\sin\theta$ | $\cos\theta$ | $\tan\theta$ |
|---|---|---|---|---|
| I | $0<\theta<\pi/2$ | $+$ | $+$ | $+$ |
| II | $\pi/2<\theta<\pi$ | $+$ | $-$ | $-$ |
| III | $\pi<\theta<3\pi/2$ | $-$ | $-$ | $+$ |
| IV | $3\pi/2<\theta<2\pi$ | $-$ | $+$ | $-$ |

### Special angles and reference angles

The values at multiples of $30^\circ$ and $45^\circ$ come from the equilateral triangle (halved) and the isosceles right triangle, and are worth knowing without a calculator:

| $\theta$ | $0$ | $\pi/6$ ($30^\circ$) | $\pi/4$ ($45^\circ$) | $\pi/3$ ($60^\circ$) | $\pi/2$ ($90^\circ$) |
|---|---|---|---|---|---|
| $\sin\theta$ | $0$ | $1/2$ | $\sqrt2/2$ | $\sqrt3/2$ | $1$ |
| $\cos\theta$ | $1$ | $\sqrt3/2$ | $\sqrt2/2$ | $1/2$ | $0$ |
| $\tan\theta$ | $0$ | $1/\sqrt3$ | $1$ | $\sqrt3$ | undefined |

For any other quadrant, the **reference angle** $\theta_r\in[0,\pi/2]$ is the acute angle between the terminal ray and the $x$-axis. Then $\lvert\sin\theta\rvert=\sin\theta_r$, $\lvert\cos\theta\rvert=\cos\theta_r$, and the sign comes from the quadrant table. For example, $\theta=240^\circ$ lies in quadrant III with $\theta_r=60^\circ$, so $\sin240^\circ=-\sqrt3/2$ and $\cos240^\circ=-1/2$. These two numbers reappear constantly in three-phase work, where phases are $0^\circ$, $-120^\circ$ and $-240^\circ$ (equivalently $+120^\circ$).

### Inverse functions and the quadrant problem

Because $\sin$, $\cos$ and $\tan$ are periodic they are not one-to-one, so their inverses are defined by restricting the domain to a **principal branch**:

> [!definition] Principal inverse trigonometric functions
> $$
> \begin{aligned}
> \arcsin x &\in[-\pi/2,\ \pi/2], & -1\le x\le 1,\\
> \arccos x &\in[0,\ \pi], & -1\le x\le 1,\\
> \arctan x &\in(-\pi/2,\ \pi/2), & x\in\mathbb{R}.
> \end{aligned}
> $$
> Each returns the *one* angle in its range with the given sine, cosine or tangent. Every other solution differs by the symmetries above: $\sin\theta=s$ has solutions $\theta=\arcsin s+2\pi k$ and $\theta=\pi-\arcsin s+2\pi k$; $\cos\theta=c$ has $\theta=\pm\arccos c+2\pi k$; $\tan\theta=t$ has $\theta=\arctan t+\pi k$.

The engineering problem is usually: *given the coordinates $(x,y)$ of a point (or the real and imaginary parts of a phasor), find its angle.* The ratio $y/x$ loses the signs of $x$ and $y$ separately, so $\arctan(y/x)$ cannot distinguish $(1,1)$ from $(-1,-1)$. The two-argument arctangent fixes this.

> [!definition] Two-argument arctangent
> For $(x,y)\neq(0,0)$, $\operatorname{atan2}(y,x)$ is the unique $\theta\in(-\pi,\pi]$ with $\cos\theta = x/\sqrt{x^2+y^2}$ and $\sin\theta=y/\sqrt{x^2+y^2}$. In terms of $\arctan$:
> $$
> \operatorname{atan2}(y,x)=\begin{cases}
> \arctan(y/x), & x>0,\\
> \arctan(y/x)+\pi, & x<0,\ y\ge0,\\
> \arctan(y/x)-\pi, & x<0,\ y<0,\\
> +\pi/2, & x=0,\ y>0,\\
> -\pi/2, & x=0,\ y<0.
> \end{cases}
> $$

> [!example] Worked Example 3.3 — Solving trigonometric equations on a full period
> (a) Find all $\theta\in[0,2\pi)$ with $\sin\theta=-0.5$. (b) Find the angle of the point $(-4,-3)$. (c) Solve $\cos(2\theta)=0.3$ for $\theta\in[0,2\pi)$.
>
> **(a)** $\arcsin(-0.5)=-\pi/6$. The two families are $\theta=-\pi/6+2\pi k$ and $\theta=\pi-(-\pi/6)+2\pi k = 7\pi/6+2\pi k$. In $[0,2\pi)$: $\theta=7\pi/6$ ($210^\circ$) and $\theta=11\pi/6$ ($330^\circ$). Both are in quadrants III and IV, where sine is negative. ✓
>
> **(b)** $\arctan(-3/-4)=\arctan(0.75)=36.87^\circ$ — a first-quadrant answer for a third-quadrant point. Since $x<0,\ y<0$: $\operatorname{atan2}(-3,-4)=36.87^\circ-180^\circ=-143.13^\circ$ (or $216.87^\circ$). Check: $5\cos(-143.13^\circ)=-4.00$, $5\sin(-143.13^\circ)=-3.00$. ✓
>
> **(c)** Let $u=2\theta$; as $\theta$ ranges over $[0,2\pi)$, $u$ ranges over $[0,4\pi)$ — *two* full turns, so expect four solutions. $\arccos0.3=1.2661$ rad. Solutions $u=\pm1.2661+2\pi k$ in $[0,4\pi)$: $u=1.2661,\ 5.0171,\ 7.5493,\ 11.3003$. Dividing by 2: $\theta = 0.6331,\ 2.5086,\ 3.7746,\ 5.6502$ rad.
>
> **Check:** $\cos(2\times2.5086)=\cos(5.0171)=0.300$. ✓

> [!warning] Common Pitfall — Calculator and library angle modes
> `sin(30)` is $-0.988$ in radian mode and $0.5$ in degree mode. Every mainstream programming library (C `math.h`, Python `math`, NumPy, MATLAB `sin`) works in radians; MATLAB's `sind`, `cosd` and spreadsheet `RADIANS()` exist precisely because people forget. A second trap is `atan` vs `atan2`: using `atan(y/x)` to compute the phase of a phasor gives an answer off by exactly $180^\circ$ whenever the real part is negative.

### Why calculus insists on radians

The decisive reason for radians is the limit

$$
\lim_{\theta\to0}\frac{\sin\theta}{\theta}=1\qquad(\theta\ \text{in radians}).
$$

**Proof sketch.** For $0<\theta<\pi/2$ draw the unit circle, the point $P=(\cos\theta,\sin\theta)$ and the point $Q=(1,\tan\theta)$ on the tangent line $x=1$. The triangle $OAP$ (with $A=(1,0)$), the sector $OAP$ and the triangle $OAQ$ are nested, so their areas satisfy

$$
\tfrac12\sin\theta\ \le\ \tfrac12\theta\ \le\ \tfrac12\tan\theta.
$$

The middle term uses the radian sector-area formula $\tfrac12r^2\theta$. Dividing by $\tfrac12\sin\theta>0$ and taking reciprocals gives $\cos\theta\le\dfrac{\sin\theta}{\theta}\le1$. As $\theta\to0^+$, $\cos\theta\to1$, and the squeeze gives the limit; evenness of $\sin\theta/\theta$ covers $\theta\to0^-$. $\square$

This limit is exactly what makes $\dfrac{d}{d\theta}\sin\theta=\cos\theta$ (derived in [Ch. 4](../part-2-calculus/04-differential-calculus.md)). If $\theta$ were measured in degrees, every derivative would carry a factor $\pi/180$: $\frac{d}{d\theta^\circ}\sin\theta^\circ=\frac{\pi}{180}\cos\theta^\circ$. The exponential form $e^{j\theta}=\cos\theta+j\sin\theta$ of [Ch. 2](02-complex-numbers.md) likewise holds only with $\theta$ in radians.

### Small-angle approximations

The limit above says $\sin\theta\approx\theta$ for small $\theta$. Combined with the Taylor series (developed in [Ch. 4](../part-2-calculus/04-differential-calculus.md))

$$
\sin\theta=\theta-\frac{\theta^3}{6}+\cdots,\qquad \cos\theta=1-\frac{\theta^2}{2}+\cdots,\qquad \tan\theta=\theta+\frac{\theta^3}{3}+\cdots,
$$

we have the working approximations $\sin\theta\approx\tan\theta\approx\theta$ and $\cos\theta\approx1-\theta^2/2$, with relative error of $\sin\theta\approx\theta$ about $\theta^2/6$.

> [!example] Worked Example 3.4 — How small is "small"?
> For what angle is the relative error of $\sin\theta\approx\theta$ equal to $1\%$? Evaluate the approximation at $\theta=0.1$ rad and at $10^\circ$.
>
> **Solution.** Relative error $\approx\theta^2/6=0.01$ gives $\theta\approx\sqrt{0.06}=0.245$ rad $\approx14^\circ$.
>
> At $\theta=0.1$: $\sin0.1=0.099833$, error $0.1-0.099833=1.67\times10^{-4}$, relative error $0.167\%$, matching $\theta^2/6=0.00167$. ✓
>
> At $10^\circ=0.17453$ rad: $\sin10^\circ=0.17365$; relative error $(0.17453-0.17365)/0.17365=0.51\%$; prediction $0.17453^2/6=0.51\%$. ✓
>
> So the approximation is good to about 1 % up to roughly $14^\circ$, which is why it is used freely for small phase errors, small pointing errors and small-signal analysis of phase-locked loops.

```desmos-graph
left=-3.2; right=3.2; bottom=-2; top=2
---
y=\sin(x)|blue
y=x|red|dashed
y=x-\frac{x^{3}}{6}|green|dotted
```
*$\sin x$ (blue) with its linear approximation $x$ (red) and cubic approximation $x-x^3/6$ (green): both hug the sine near the origin and peel away for larger angles.*

## 3.2 Sinusoidal functions

### Motivation

A sinusoid is the only periodic waveform whose shape survives every linear, time-invariant operation: scaling it, delaying it, differentiating it, integrating it, or adding it to another sinusoid of the same frequency all produce a sinusoid of the same frequency. Resistors scale, capacitors and inductors differentiate or integrate, transmission lines delay — so in a linear circuit driven by a sinusoid, *every* voltage and current is a sinusoid of the source frequency, differing only in amplitude and phase. That is the mathematical reason why AC analysis, phasors ([Ch. 2](02-complex-numbers.md)) and frequency response ([Ch. 14](../part-5-transforms/14-laplace-transforms.md), [Ch. 16](../part-5-transforms/16-fourier-transforms.md)) all reduce to the bookkeeping of two numbers per signal.

### The general sinusoid

> [!definition] Sinusoid
> A **sinusoid** (sinusoidal function of time) is any function of the form
> $$
> x(t)=A\cos(\omega t+\phi)+C,
> $$
> with
> - **amplitude** (peak value) $A>0$, in the units of $x$;
> - **angular frequency** $\omega>0$, in rad/s;
> - **phase** (phase angle, initial phase) $\phi$, in rad or degrees, conventionally reduced to $(-\pi,\pi]$;
> - **offset** (DC component, mean value) $C$, in the units of $x$.
>
> The quantity $\Theta(t)=\omega t+\phi$ is the **instantaneous phase** (total phase angle) at time $t$; $\phi=\Theta(0)$.

The electrical-engineering convention is to use the **cosine** as the reference sinusoid, because $\cos(\omega t+\phi)=\operatorname{Re}\{e^{j(\omega t+\phi)}\}$ connects directly to the phasor $Ae^{j\phi}$. Many physics and mathematics texts use $\sin$ as the reference; both describe the same family of functions, and you must convert before comparing phases (§3.3).

The four parameters act on the graph of $\cos t$ as elementary transformations (Ch. 1):

| Parameter | Effect on the graph of $\cos t$ |
|---|---|
| $A$ | vertical stretch: peaks at $C\pm A$ |
| $C$ | vertical shift: the waveform oscillates about the line $x=C$ |
| $\omega$ | horizontal compression by the factor $\omega$: one cycle takes $2\pi/\omega$ seconds |
| $\phi$ | horizontal shift by $-\phi/\omega$: the peak that $\cos t$ has at $t=0$ moves to $t=-\phi/\omega$ |

Note that the horizontal shift is $\phi/\omega$, **not** $\phi$: writing $\cos(\omega t+\phi)=\cos\big(\omega(t+\phi/\omega)\big)$ shows that the graph is shifted left by $\phi/\omega$ seconds. A positive phase moves the waveform *earlier* in time.

```desmos-graph
left=-1; right=7; bottom=-2.5; top=4.5
---
y=\cos(x)|grey|dashed
y=1+2\cos(2x-\frac{\pi}{3})|blue
y=1|red|dotted
(\frac{\pi}{6},3)|label:first peak
```
*The reference $\cos t$ (grey) and $1+2\cos(2t-\pi/3)$ (blue): amplitude 2, offset 1 (red line), angular frequency 2 so the period halves to $\pi$, and a phase of $-\pi/3$ that delays the first peak to $t=\pi/6$.*

### Elementary properties

For $x(t)=A\cos(\omega t+\phi)+C$:

- **Extremes.** Maximum $C+A$ where $\omega t+\phi=2\pi k$; minimum $C-A$ where $\omega t+\phi=\pi+2\pi k$.
- **Peak-to-peak value.** $x_{pp}=2A$.
- **Crossings of the mean.** $x=C$ where $\omega t+\phi=\pi/2+\pi k$, i.e. twice per cycle, a quarter-cycle after each extreme.
- **Period.** $T=2\pi/\omega$ (§3.4).
- **Boundedness and smoothness.** $x$ is bounded and infinitely differentiable; all of its derivatives are sinusoids of the same frequency.

### Converting between sine and cosine and normalising the phase

Because the sine curve is the cosine curve delayed by a quarter cycle, every sinusoid can be written in the cosine reference with a positive amplitude. The conversion rules, read off the unit circle, are:

$$
\begin{aligned}
\sin\theta &= \cos(\theta-\tfrac{\pi}{2}), & -\sin\theta &= \cos(\theta+\tfrac{\pi}{2}),\\
-\cos\theta &= \cos(\theta\pm\pi), & \cos\theta &= \sin(\theta+\tfrac{\pi}{2}).
\end{aligned}
$$

In degrees: $\sin\theta=\cos(\theta-90^\circ)$, $-\sin\theta=\cos(\theta+90^\circ)$, $-\cos\theta=\cos(\theta\pm180^\circ)$. A negative amplitude is absorbed into the phase as $\pm180^\circ$, chosen so that the resulting phase lands in $(-180^\circ,180^\circ]$.

> [!definition] Standard (canonical) form
> A sinusoid is in **standard form** when written $A\cos(\omega t+\phi)+C$ with $A>0$, $\omega>0$ and $-\pi<\phi\le\pi$. In standard form the parameters are unique, so two sinusoids are equal if and only if their standard-form parameters agree.

A negative frequency can likewise be removed, using evenness of cosine: $\cos(-\omega t+\phi)=\cos(\omega t-\phi)$.

> [!example] Worked Example 3.5 — Reading parameters and putting a sinusoid into standard form
> Write each signal in standard form and state amplitude, angular frequency, frequency, phase and offset.
>
> (a) $v(t)=-170\sin(377t+20^\circ)$ V. (b) $i(t)=5-3\cos(2000\pi t-60^\circ)$ mA. (c) $x(t)=4\sin(-50t+0.3)$.
>
> **(a)** $-\sin\theta=\cos(\theta+90^\circ)$, so $v=170\cos(377t+110^\circ)$ V. Amplitude 170 V, $\omega=377$ rad/s, $f=377/2\pi=60.0$ Hz, phase $+110^\circ$, offset 0. (This is a 120 V RMS, 60 Hz supply: $170/\sqrt2=120.2$ V.)
>
> **(b)** $-\cos\theta=\cos(\theta+180^\circ)$: $-3\cos(2000\pi t-60^\circ)=3\cos(2000\pi t+120^\circ)$. So $i=5+3\cos(2000\pi t+120^\circ)$ mA: amplitude 3 mA, $\omega=2000\pi=6283$ rad/s, $f=1$ kHz, phase $120^\circ$, offset 5 mA. The current never reverses: it swings between 2 mA and 8 mA.
>
> **(c)** Sine is odd: $4\sin(-50t+0.3)=-4\sin(50t-0.3)=4\cos(50t-0.3+\pi/2)=4\cos(50t+1.2708)$. Amplitude 4, $\omega=50$ rad/s, $f=7.958$ Hz, phase $1.2708$ rad $=72.8^\circ$.
>
> **Check (c) at $t=0.01$:** original $4\sin(-0.5+0.3)=4\sin(-0.2)=-0.7947$; standard form $4\cos(0.5+1.2708)=4\cos(1.7708)=-0.7947$. ✓

### Amplitude–phase form: combining $a\cos\omega t+b\sin\omega t$

Solutions of differential equations ([Ch. 11](../part-4-differential-equations/11-ordinary-differential-equations.md), [Ch. 12](../part-4-differential-equations/12-linear-differential-equations.md)) and Fourier series ([Ch. 15](../part-5-transforms/15-fourier-series.md)) naturally come out as $a\cos\omega t+b\sin\omega t$. For interpretation we want a single sinusoid.

> [!theorem] Amplitude–phase form
> For real $a,b$ not both zero,
> $$
> a\cos\omega t+b\sin\omega t = R\cos(\omega t-\alpha),\qquad R=\sqrt{a^2+b^2},\quad \alpha=\operatorname{atan2}(b,a).
> $$
> Equivalently $a\cos\omega t+b\sin\omega t=R\cos(\omega t+\phi)$ with $\phi=-\alpha=\operatorname{atan2}(-b,a)$. Conversely, $R\cos(\omega t+\phi)=a\cos\omega t+b\sin\omega t$ with $a=R\cos\phi$, $b=-R\sin\phi$.

**Derivation 1 (compound angle).** Expand the right-hand side with $\cos(X-Y)=\cos X\cos Y+\sin X\sin Y$ (proved in §3.6):

$$
R\cos(\omega t-\alpha)=R\cos\alpha\,\cos\omega t+R\sin\alpha\,\sin\omega t.
$$

This equals $a\cos\omega t+b\sin\omega t$ for all $t$ if and only if $R\cos\alpha=a$ and $R\sin\alpha=b$ (compare values at $t=0$ and at $\omega t=\pi/2$). These say that $(a,b)$ is the point with polar coordinates $(R,\alpha)$, so $R=\sqrt{a^2+b^2}$ and $\alpha=\operatorname{atan2}(b,a)$.

**Derivation 2 (phasor).** Since $\cos\omega t=\operatorname{Re}\{e^{j\omega t}\}$ and $\sin\omega t=\operatorname{Re}\{-je^{j\omega t}\}$,

$$
a\cos\omega t+b\sin\omega t=\operatorname{Re}\{(a-jb)e^{j\omega t}\}.
$$

The phasor is $a-jb=Re^{-j\alpha}$, so the sum is $\operatorname{Re}\{Re^{j(\omega t-\alpha)}\}=R\cos(\omega t-\alpha)$. The phasor route makes the sign of the phase automatic: **the sine coefficient enters the phasor with a minus sign.**

> [!example] Worked Example 3.6 — Amplitude–phase conversions
> (a) $3\cos\omega t+4\sin\omega t$. (b) $5\cos\omega t-12\sin\omega t$. (c) $-2\sqrt3\cos\omega t-2\sin\omega t$.
>
> **(a)** $R=\sqrt{9+16}=5$, $\alpha=\operatorname{atan2}(4,3)=53.13^\circ$: $\ 3\cos\omega t+4\sin\omega t=5\cos(\omega t-53.13^\circ)$.
>
> **(b)** $R=\sqrt{25+144}=13$, $\alpha=\operatorname{atan2}(-12,5)=-67.38^\circ$: $\ 5\cos\omega t-12\sin\omega t=13\cos(\omega t+67.38^\circ)$.
>
> **(c)** $R=\sqrt{12+4}=4$. Here $a<0$, $b<0$, so the point is in quadrant III: $\alpha=\operatorname{atan2}(-2,-2\sqrt3)=-150^\circ$. Hence $-2\sqrt3\cos\omega t-2\sin\omega t=4\cos(\omega t+150^\circ)$.
>
> **Check (c):** $4\cos(\omega t+150^\circ)=4\cos150^\circ\cos\omega t-4\sin150^\circ\sin\omega t=4(-0.8660)\cos\omega t-4(0.5)\sin\omega t=-3.464\cos\omega t-2\sin\omega t$. ✓ A naive $\arctan(b/a)=\arctan(1/\sqrt3)=30^\circ$ would have given $4\cos(\omega t-30^\circ)$ — the negative of the correct answer.

### Derivatives and integrals of sinusoids

Using $\frac{d}{dt}\cos(\omega t+\phi)=-\omega\sin(\omega t+\phi)$ and the conversion rule $-\sin\theta=\cos(\theta+\pi/2)$:

$$
\frac{d}{dt}\big[A\cos(\omega t+\phi)\big]=\omega A\cos\!\big(\omega t+\phi+\tfrac{\pi}{2}\big),\qquad
\int A\cos(\omega t+\phi)\,dt=\frac{A}{\omega}\cos\!\big(\omega t+\phi-\tfrac{\pi}{2}\big)+K.
$$

**Differentiation multiplies the amplitude by $\omega$ and advances the phase by $90^\circ$; integration divides by $\omega$ and retards the phase by $90^\circ$.** In phasor language these are multiplication and division by $j\omega$, since $j=e^{j\pi/2}$. Differentiating twice gives back the original sinusoid times $-\omega^2$, so every sinusoid of angular frequency $\omega$ (without offset) satisfies

$$
\ddot x+\omega^2x=0,
$$

the equation of simple harmonic motion (§3.5).

> [!example] Worked Example 3.7 — Capacitor current from capacitor voltage
> A $C=1\ \mu$F capacitor has $v(t)=10\cos(1000t+30^\circ)$ V. Find $i(t)=C\,dv/dt$ and state which leads.
>
> **Solution.** $dv/dt = 1000\times10\cos(1000t+30^\circ+90^\circ)=10^4\cos(1000t+120^\circ)$ V/s. So $i=10^{-6}\times10^4\cos(1000t+120^\circ)=10\cos(1000t+120^\circ)$ mA.
>
> The current phase ($120^\circ$) exceeds the voltage phase ($30^\circ$) by $90^\circ$: **current leads voltage by $90^\circ$** in a capacitor. The amplitude ratio $V/I=10/0.01=1000\ \Omega=1/(\omega C)$ is the capacitive reactance. ✓

### Mean, RMS and related measures

The **mean** (average) of a periodic function over one period $T$ is $\bar x=\frac1T\int_0^Tx\,dt$, and the **root-mean-square** value is $X_{\text{rms}}=\sqrt{\frac1T\int_0^Tx^2\,dt}$. RMS matters because the average power delivered to a resistor is $\overline{v^2}/R=V_{\text{rms}}^2/R$: the RMS value is the DC value that would deliver the same power. The general theory of mean and RMS values belongs to [Ch. 5](../part-2-calculus/05-integral-calculus.md); here we need the sinusoidal case.

> [!theorem] Mean and RMS of a sinusoid
> Over any whole number of periods,
> $$
> \overline{\cos(\omega t+\phi)}=0,\qquad \overline{\cos^2(\omega t+\phi)}=\tfrac12.
> $$
> Hence for $x=A\cos(\omega t+\phi)+C$:
> $$
> \bar x = C,\qquad X_{\text{rms}}=\sqrt{C^2+\tfrac12A^2}.
> $$
> For a pure sinusoid ($C=0$), $X_{\text{rms}}=A/\sqrt2\approx0.7071A$.

**Proof.** The integral of $\cos(\omega t+\phi)$ over one period is $\frac1\omega[\sin(\omega t+\phi)]_0^T=\frac1\omega[\sin(2\pi+\phi)-\sin\phi]=0$. For the square, use the power-reduction identity $\cos^2\theta=\frac12(1+\cos2\theta)$ (§3.6): the constant $\frac12$ averages to $\frac12$ and $\cos2\theta$ averages to $0$. Then

$$
\overline{x^2}=\overline{C^2+2AC\cos(\cdot)+A^2\cos^2(\cdot)}=C^2+0+\tfrac12A^2.\qquad\square
$$

There is a neat way to see $\overline{\cos^2}=\tfrac12$ without any identity: $\cos^2$ and $\sin^2$ have the same graph shifted by a quarter period, so they have the same average over a full period; and they add to 1. Each must average $\tfrac12$.

Two further ratios describe waveform shape and are printed on meter datasheets:

$$
\text{crest factor}=\frac{X_{\text{peak}}}{X_{\text{rms}}},\qquad \text{form factor}=\frac{X_{\text{rms}}}{\overline{\lvert x\rvert}}.
$$

For a pure sinusoid, the rectified average is $\overline{\lvert A\cos\omega t\rvert}=2A/\pi$ (computed in §3.8), so crest factor $=\sqrt2\approx1.414$ and form factor $=\dfrac{A/\sqrt2}{2A/\pi}=\dfrac{\pi}{2\sqrt2}\approx1.1107$.

> [!example] Worked Example 3.8 — Mains voltage and a biased sinusoid
> (a) A 230 V (RMS), 50 Hz supply: what is the peak voltage, the peak-to-peak voltage, and the instantaneous voltage 3 ms after a positive peak?
> (b) A sensor output is $v(t)=2+3\sin(100t)$ V. Find its mean, its RMS value and the power it would deliver to a 50 $\Omega$ resistor.
>
> **(a)** $V_p=\sqrt2\times230=325.3$ V; $V_{pp}=650.5$ V. With the peak at $t=0$, $v=325.3\cos(2\pi\cdot50\,t)$. At $t=3$ ms: $\omega t=2\pi(50)(0.003)=0.9425$ rad $=54^\circ$, so $v=325.3\cos54^\circ=191.2$ V.
>
> **(b)** Mean $=2$ V. RMS $=\sqrt{2^2+3^2/2}=\sqrt{8.5}=2.915$ V. Power $=V_{\text{rms}}^2/R=8.5/50=0.170$ W. Note that the DC and AC powers simply add: $2^2/50+(3/\sqrt2)^2/50=0.08+0.09=0.17$ W — a first glimpse of the orthogonality result of §3.7.

> [!example] Worked Example 3.9 — Recovering a sinusoid from two observations
> An oscilloscope shows a sinusoid with a maximum of $7$ V at $t=2$ ms and the next minimum of $-1$ V at $t=7$ ms. Write it in standard form.
>
> **Offset and amplitude.** $C=\frac{7+(-1)}2=3$ V, $A=\frac{7-(-1)}2=4$ V.
>
> **Frequency.** A maximum and the next minimum are half a period apart: $T/2=5$ ms, so $T=10$ ms, $f=100$ Hz, $\omega=200\pi=628.3$ rad/s.
>
> **Phase.** A maximum occurs where $\omega t+\phi=0$ (mod $2\pi$): $\phi=-\omega(0.002)=-0.4\pi$ rad $=-72^\circ$.
>
> **Result.** $v(t)=3+4\cos(200\pi t-0.4\pi)$ V.
>
> **Check:** at $t=7$ ms, $200\pi(0.007)-0.4\pi=1.0\pi$, $\cos\pi=-1$, $v=3-4=-1$ V. ✓

> [!warning] Common Pitfall — Peak, peak-to-peak or RMS?
> "10 V AC" is ambiguous. Function generators often display peak-to-peak (into a matched load!), multimeters display RMS, and the formula $A\cos(\omega t+\phi)$ uses peak. Convert explicitly: $V_{pp}=2V_p=2\sqrt2\,V_{\text{rms}}\approx2.83\,V_{\text{rms}}$. Also note that the $1/\sqrt2$ rule is for *sinusoids only*: an "average-responding, RMS-calibrated" meter multiplies the rectified average by 1.1107 and is wrong for every non-sinusoidal waveform.

## 3.3 Phase and phase difference

### Motivation

Amplitude tells you *how big* a sinusoid is; phase tells you *when* it peaks. In a single isolated signal the phase is arbitrary — it depends only on where you chose $t=0$. What carries physical meaning is the **phase difference** between two sinusoids of the same frequency: the angle by which the current in a motor lags its voltage (which fixes the power factor), the phase shift through an amplifier (which decides whether a feedback loop oscillates), the phase between the three lines of a power system, or the phase between signals at two antennas (which steers a beam). This section makes "leads", "lags" and "phase difference" precise and shows how to compute them without sign errors.

### Phase of a single sinusoid

For $x(t)=A\cos(\omega t+\phi)$ the phase $\phi$ is the instantaneous phase at $t=0$. Moving the time origin by $t_0$ (that is, writing $t=t'+t_0$) changes the phase to $\phi+\omega t_0$ but leaves the waveform itself untouched. Hence:

- the phase of one sinusoid is meaningful only relative to a stated time reference (or a reference signal);
- phases are defined modulo $2\pi$ ($360^\circ$), and are normally reported in the principal interval $(-180^\circ,180^\circ]$ — this reduction is called **phase wrapping**.

### Phase difference, lead and lag

> [!definition] Phase difference, lead, lag
> Let $x_1(t)=A_1\cos(\omega t+\phi_1)$ and $x_2(t)=A_2\cos(\omega t+\phi_2)$ have the **same** angular frequency $\omega$ and be written in standard form ($A_1,A_2>0$, both cosines). The **phase difference** of $x_1$ relative to $x_2$ is
> $$
> \Delta\phi=\phi_1-\phi_2,\qquad\text{reduced to }(-\pi,\pi].
> $$
> If $0<\Delta\phi<\pi$, $x_1$ **leads** $x_2$ by $\Delta\phi$; if $-\pi<\Delta\phi<0$, $x_1$ **lags** $x_2$ by $\lvert\Delta\phi\rvert$. The signals are **in phase** if $\Delta\phi=0$, in **quadrature** if $\Delta\phi=\pm\pi/2$, and in **antiphase** (phase opposition) if $\Delta\phi=\pi$.

"Leads" means *happens earlier*: if $x_1$ leads $x_2$ by $\Delta\phi$, each feature of $x_1$ (peak, upward zero crossing) occurs $\Delta\phi/\omega$ seconds before the corresponding feature of $x_2$. Indeed

$$
x_1(t)=A_1\cos\big(\omega t+\phi_2+\Delta\phi\big)=\frac{A_1}{A_2}\,x_2\!\left(t+\frac{\Delta\phi}{\omega}\right),
$$

so $x_1$ is a scaled copy of $x_2$ advanced by $\Delta\phi/\omega$.

The restriction to $(-\pi,\pi]$ matters: a lead of $270^\circ$ is the same physical relationship as a lag of $90^\circ$, and engineers always quote the smaller one. A phase difference of exactly $180^\circ$ is neither a lead nor a lag.

> [!theorem] Phase difference as time delay
> For two sinusoids of the same frequency $f$ (period $T$), a phase difference $\Delta\phi$ corresponds to a time shift
> $$
> \Delta t=\frac{\Delta\phi}{\omega}=\frac{\Delta\phi}{2\pi f}=\frac{\Delta\phi_{\text{deg}}}{360^\circ}\,T .
> $$
> Conversely a pure time delay $\tau$ produces a phase lag $\omega\tau$, which grows linearly with frequency.

```desmos-graph
left=-0.5; right=7; bottom=-1.6; top=1.6
---
y=\cos(x)|blue
y=0.7\cos(x-\frac{\pi}{3})|red
x=0|blue|dashed
x=\frac{\pi}{3}|red|dashed
```
*$\cos t$ (blue) leads $0.7\cos(t-\pi/3)$ (red) by $60^\circ$: every red peak arrives $\pi/3$ rad of $\omega t$ after the corresponding blue peak (dashed lines).*

> [!example] Worked Example 3.10 — Comparing phases written in different forms
> $v(t)=10\sin(\omega t+30^\circ)$ V and $i(t)=-2\cos(\omega t-20^\circ)$ A. Find the phase of $i$ relative to $v$.
>
> **Step 1 — convert both to standard cosine form.**
> $v=10\cos(\omega t+30^\circ-90^\circ)=10\cos(\omega t-60^\circ)$.
> $i=2\cos(\omega t-20^\circ+180^\circ)=2\cos(\omega t+160^\circ)$.
>
> **Step 2 — subtract.** $\Delta\phi=\phi_i-\phi_v=160^\circ-(-60^\circ)=220^\circ$.
>
> **Step 3 — wrap.** $220^\circ-360^\circ=-140^\circ$. So **$i$ lags $v$ by $140^\circ$** (equivalently, $v$ leads $i$ by $140^\circ$).
>
> **Check.** Comparing the original expressions directly ($30^\circ$ vs $-20^\circ$, "$v$ leads by $50^\circ$") is wrong on two counts: one is a sine and the other a cosine, and $i$ has a negative amplitude. Numerically, $v$ has a positive peak at $\omega t=60^\circ$; $i$ has a positive peak at $\omega t=-160^\circ\equiv200^\circ$, which is $140^\circ$ later. ✓

> [!example] Worked Example 3.11 — From phase to time and back
> (a) Two 50 Hz voltages differ in phase by $36^\circ$. What is the time between their peaks?
> (b) A coaxial cable delays a signal by $5$ ns per metre. What phase lag does 2 m of cable introduce at 10 MHz, at 100 MHz and at 150 MHz?
> (c) An RC low-pass filter with $R=1$ k$\Omega$, $C=100$ nF has output phase $-\arctan(\omega RC)$ relative to its input. At 1 kHz, what are the phase lag and the corresponding time delay?
>
> **(a)** $T=20$ ms; $\Delta t=(36/360)\times20=2$ ms.
>
> **(b)** $\tau=10$ ns. Lag $=360^\circ f\tau$: at 10 MHz, $36^\circ$; at 100 MHz, $360^\circ$ (indistinguishable from $0^\circ$!); at 150 MHz, $540^\circ\equiv180^\circ$.
>
> **(c)** $\omega RC=2\pi(1000)(1000)(10^{-7})=0.6283$. Lag $=\arctan0.6283=32.14^\circ$. Delay $=\frac{32.14}{360}\times1\text{ ms}=89.3\ \mu$s.
>
> **Interpretation.** Part (b) shows that a measured phase cannot distinguish delays that differ by whole periods; part (c) shows that the phase lag of a filter, unlike that of a cable, is *not* proportional to frequency, so different frequencies are delayed by different amounts.

> [!warning] Common Pitfall — Phase ambiguity
> A phase measurement determines a time delay only modulo one period. Phase meters, network analysers and FFT routines report a wrapped phase in $(-180^\circ,180^\circ]$; to recover total delay you must *unwrap* — track the phase continuously as frequency is swept and add $360^\circ$ at each jump. NumPy's `unwrap` and MATLAB's `unwrap` do exactly this, assuming that consecutive samples differ by less than $180^\circ$.

### Adding sinusoids of the same frequency

The sum of two sinusoids of the same frequency is a sinusoid of that frequency. The cleanest proof uses phasors ([Ch. 2](02-complex-numbers.md)):

$$
A_1\cos(\omega t+\phi_1)+A_2\cos(\omega t+\phi_2)=\operatorname{Re}\big\{(A_1e^{j\phi_1}+A_2e^{j\phi_2})e^{j\omega t}\big\}=R\cos(\omega t+\psi),
$$

where $Re^{j\psi}=A_1e^{j\phi_1}+A_2e^{j\phi_2}$. Writing out real and imaginary parts gives explicit formulas:

> [!theorem] Sum of two same-frequency sinusoids
> $$
> R^2=A_1^2+A_2^2+2A_1A_2\cos(\phi_1-\phi_2),\qquad
> \psi=\operatorname{atan2}\big(A_1\sin\phi_1+A_2\sin\phi_2,\ A_1\cos\phi_1+A_2\cos\phi_2\big).
> $$

The magnitude formula is the law of cosines applied to the phasor triangle. It shows that the resultant amplitude depends only on the *phase difference*: it ranges from $A_1+A_2$ (in phase, **constructive** interference) down to $\lvert A_1-A_2\rvert$ (antiphase, **destructive** interference), and equals $\sqrt{A_1^2+A_2^2}$ in quadrature.

> [!example] Worked Example 3.12 — Adding two currents at a node
> Two branch currents enter a node: $i_1=10\cos(\omega t+30^\circ)$ A and $i_2=6\cos(\omega t-45^\circ)$ A. Find $i_1+i_2$.
>
> **Rectangular components.**
> $10\angle30^\circ=8.660+j5.000$; $\ 6\angle-45^\circ=4.243-j4.243$.
> Sum: $12.903+j0.757$.
>
> **Polar.** $R=\sqrt{12.903^2+0.757^2}=12.925$; $\psi=\operatorname{atan2}(0.757,12.903)=3.36^\circ$.
>
> **Result.** $i_1+i_2=12.93\cos(\omega t+3.36^\circ)$ A.
>
> **Check with the law of cosines:** $R^2=100+36+120\cos75^\circ=136+31.06=167.06$, $R=12.925$. ✓ Also at $\omega t=0.3$ rad the original sum evaluates to $12.103$ and the result to $12.103$. ✓ Note that the resultant (12.9 A) is less than the arithmetic sum (16 A) — an ammeter reading RMS values on each branch does *not* add.

### Three-phase sets

A **balanced three-phase set** consists of three sinusoids of equal amplitude and frequency displaced by $120^\circ$:

$$
v_a=V\cos\omega t,\qquad v_b=V\cos(\omega t-120^\circ),\qquad v_c=V\cos(\omega t+120^\circ).
$$

> [!theorem] Balanced sets sum to zero
> $v_a+v_b+v_c=0$ for all $t$. More generally, for any integer $N\ge2$, $\displaystyle\sum_{k=0}^{N-1}\cos\!\Big(\omega t+\phi-\frac{2\pi k}{N}\Big)=0$.

**Proof.** The phasor sum is $V\sum_{k=0}^{N-1}e^{-j2\pi k/N}$, a geometric series with ratio $r=e^{-j2\pi/N}\neq1$: it equals $V\frac{1-r^N}{1-r}=0$ because $r^N=e^{-j2\pi}=1$. (These are the $N$-th roots of unity of Ch. 2, which sum to zero.) $\square$

This is why a balanced three-phase load needs no neutral return current.

> [!example] Worked Example 3.13 — Line-to-line voltage
> With $v_a$ and $v_b$ as above, find the line-to-line voltage $v_{ab}=v_a-v_b$.
>
> **Phasors.** $V\angle0^\circ-V\angle-120^\circ=V\big[1-(-\tfrac12-j\tfrac{\sqrt3}{2})\big]=V\big(\tfrac32+j\tfrac{\sqrt3}2\big)$.
>
> Magnitude $V\sqrt{9/4+3/4}=\sqrt3\,V$; angle $\arctan\big(\frac{\sqrt3/2}{3/2}\big)=30^\circ$.
>
> **Result.** $v_{ab}=\sqrt3\,V\cos(\omega t+30^\circ)$: the line voltage is $\sqrt3$ times the phase voltage and leads $v_a$ by $30^\circ$. For a 230 V phase supply, the line voltage is $\sqrt3\times230=398\approx400$ V.

### Measuring phase: zero crossings and Lissajous figures

Two practical measurement methods follow directly from the definitions.

**Zero-crossing timing.** Measure the time $\Delta t$ between corresponding upward zero crossings of the two signals and the period $T$; then $\lvert\Delta\phi\rvert=360^\circ\,\Delta t/T$. The signal whose crossing comes first leads.

**Lissajous (X–Y) method.** Apply $x=X\sin\omega t$ to the horizontal input of an oscilloscope and $y=Y\sin(\omega t+\phi)$ to the vertical. Eliminating $t$ gives an ellipse:

$$
\frac{x^2}{X^2}-\frac{2xy\cos\phi}{XY}+\frac{y^2}{Y^2}=\sin^2\phi .
$$

**Derivation.** Write $s=\sin\omega t=x/X$ and expand $y/Y=s\cos\phi+\cos\omega t\,\sin\phi$, so $\cos\omega t\,\sin\phi=y/Y-(x/X)\cos\phi$. Squaring and using $\cos^2\omega t=1-x^2/X^2$:

$$
\Big(1-\frac{x^2}{X^2}\Big)\sin^2\phi=\frac{y^2}{Y^2}-\frac{2xy}{XY}\cos\phi+\frac{x^2}{X^2}\cos^2\phi,
$$

which rearranges to the ellipse equation. Setting $x=0$ gives $y=\pm Y\lvert\sin\phi\rvert$. Hence

$$
\lvert\sin\phi\rvert=\frac{y_0}{Y},
$$

where $y_0$ is the height at which the ellipse crosses the vertical axis and $Y$ is its maximum height. The ellipse collapses to a line of positive slope at $\phi=0$, to a line of negative slope at $\phi=180^\circ$, and becomes an axis-aligned ellipse at $\phi=\pm90^\circ$. The sign of $\phi$ (lead or lag) cannot be read from a static picture; the tilt distinguishes $\lvert\phi\rvert<90^\circ$ (major axis in quadrants I/III) from $\lvert\phi\rvert>90^\circ$ (quadrants II/IV).

```desmos-graph
left=-1.6; right=1.6; bottom=-1.6; top=1.6
---
\frac{x^{2}}{1}-2xy\cos(0.6435)+y^{2}=\sin(0.6435)^{2}|blue
y=0.6|red|dashed
y=1|grey|dotted
```
*Lissajous ellipse for equal amplitudes and $\phi=36.87^\circ$: it crosses the vertical axis at $y_0=0.6$ (red), its maximum height is $Y=1$ (grey), and $\sin\phi=0.6$.*

> [!example] Worked Example 3.14 — Reading a Lissajous figure
> An X–Y display shows an ellipse tilted into quadrants I and III, with maximum height 5.0 divisions and vertical-axis intercepts at $\pm3.0$ divisions. Find the magnitude of the phase difference.
>
> **Solution.** $\lvert\sin\phi\rvert=3.0/5.0=0.6$, so $\lvert\phi\rvert=\arcsin0.6=36.87^\circ$ or $180^\circ-36.87^\circ=143.13^\circ$. The tilt into quadrants I/III means $\lvert\phi\rvert<90^\circ$: $\lvert\phi\rvert=36.9^\circ$.
>
> **Check:** at $x=0$ we have $\sin\omega t=0$, $\cos\omega t=\pm1$, so $y=\pm Y\sin\phi=\pm5(0.6)=\pm3$. ✓

> [!tip] Engineering Note — Phase in feedback and power systems
> Two phase numbers dominate EE practice. In power systems the **power factor** $\cos\Delta\phi$ (with $\Delta\phi$ the voltage–current phase difference) sets how much of the volt-amperes delivered is real power: $P=V_{\text{rms}}I_{\text{rms}}\cos\Delta\phi$ — a consequence of the product-to-sum identity in §3.6. In feedback amplifiers the **phase margin** is $180^\circ$ plus the loop phase at the frequency where loop gain is 1; when the loop phase reaches $-180^\circ$ negative feedback has become positive feedback and the circuit oscillates. Both are phase *differences* of same-frequency sinusoids.

## 3.4 Period, frequency, and angular frequency

### Motivation

"How often does it repeat?" is the first question asked of any signal. The answer comes in three interchangeable currencies — period in seconds, frequency in hertz, angular frequency in radians per second — and each is natural in a different place: the oscilloscope shows periods, the spectrum analyser and the datasheet show hertz, and every formula derived from calculus or phasors uses rad/s. A factor of $2\pi$ lost between them is the most common numerical error in circuit and signal calculations. Beyond single sinusoids, we need to know when a *combination* of sinusoids repeats, because that decides whether Fourier series (Ch. 15) or Fourier transforms (Ch. 16) is the right tool, and how sampled (discrete-time) sinusoids differ from continuous ones.

### Periodic functions

> [!definition] Periodic function, fundamental period
> A function $x(t)$ defined for all real $t$ is **periodic** with period $T>0$ if
> $$
> x(t+T)=x(t)\qquad\text{for all }t .
> $$
> If $T$ is a period then so is $kT$ for every positive integer $k$. The smallest positive period, when one exists, is the **fundamental period** $T_0$; its reciprocal $f_0=1/T_0$ is the **fundamental frequency** and $\omega_0=2\pi f_0$ the **fundamental angular frequency**.

A constant function satisfies $x(t+T)=x(t)$ for every $T>0$ and so has no smallest period; by convention it is periodic with any period but has no fundamental period. (Its "frequency" is zero, which is why a constant is called the DC component.)

### Period, frequency and angular frequency of a sinusoid

For $x(t)=A\cos(\omega t+\phi)$, the cosine repeats when its argument advances by $2\pi$. Solving $\omega(t+T)+\phi=\omega t+\phi+2\pi$ gives

$$
T=\frac{2\pi}{\omega},\qquad f=\frac1T=\frac{\omega}{2\pi},\qquad \omega=2\pi f=\frac{2\pi}{T}.
$$

| Quantity | Symbol | Unit | Meaning |
|---|---|---|---|
| period | $T$ | s | time for one cycle |
| frequency | $f$ | Hz $=$ s$^{-1}$ (cycles per second) | cycles per unit time |
| angular frequency | $\omega$ | rad/s | rate of change of phase, $\omega=d\Theta/dt$ |

Dimensionally $[f]=[\omega]=\mathrm{T^{-1}}$; the radian is dimensionless, so the unit "rad/s" is formally s$^{-1}$ as well. This is exactly why a factor of $2\pi$ is never caught by a dimensional check — you must track it by keeping the units "Hz" and "rad/s" distinct in your working. The product $\omega t$ is dimensionless, as the argument of any function like $\cos$ must be.

> [!example] Worked Example 3.15 — Conversions in practice
> (a) North American mains has $\omega=377$ rad/s: find $f$ and $T$. (b) A Wi-Fi carrier is at $2.4$ GHz: find $T$ and $\omega$. (c) A scope shows 4 complete cycles across 10 divisions at 50 $\mu$s/div: find $f$ and $\omega$.
>
> **(a)** $f=377/2\pi=60.0$ Hz, $T=1/60=16.67$ ms.
>
> **(b)** $T=1/(2.4\times10^9)=416.7$ ps, $\omega=2\pi\times2.4\times10^9=1.508\times10^{10}$ rad/s.
>
> **(c)** 10 divisions $=500\ \mu$s contains 4 periods, so $T=125\ \mu$s, $f=8$ kHz, $\omega=5.027\times10^4$ rad/s.

### Time scaling and related functions

If $x(t)$ has fundamental period $T_0$, then $x(at)$ with $a>0$ has fundamental period $T_0/a$: compressing time by $a$ multiplies all frequencies by $a$. Several functions built from sinusoids have periods *shorter* than the sinusoid they come from:

| Function | Fundamental period | Reason |
|---|---|---|
| $\cos(\omega t+\phi)$ | $2\pi/\omega$ | definition |
| $\tan\omega t$ | $\pi/\omega$ | $\tan(\theta+\pi)=\tan\theta$ |
| $\lvert\cos\omega t\rvert$ | $\pi/\omega$ | full-wave rectification removes the sign |
| $\cos^2\omega t$ | $\pi/\omega$ | $\cos^2\theta=\tfrac12+\tfrac12\cos2\theta$ |
| $\cos\omega t\,\sin\omega t$ | $\pi/\omega$ | $=\tfrac12\sin2\omega t$ |

The full-wave-rectifier entry is the reason the ripple on a 50 Hz full-wave rectifier's output is at 100 Hz.

### When is a sum of sinusoids periodic?

> [!theorem] Periodicity of a sum
> Let $x_1$ and $x_2$ be periodic with fundamental periods $T_1$ and $T_2$.
> 1. If $T_1/T_2=p/q$ is **rational**, with $p,q$ coprime positive integers, then $x_1+x_2$ is periodic with period $T=qT_1=pT_2$ (the least common multiple of $T_1$ and $T_2$). Equivalently, the common frequency is the greatest common divisor of $f_1$ and $f_2$.
> 2. For sums of sinusoids, if $T_1/T_2$ is **irrational**, $x_1+x_2$ is **not** periodic.
>
> In case 1, $T$ is the fundamental period of the sum unless cancellations occur (for example $x_2=-x_1$); for a sum of sinusoids with nonzero amplitudes and distinct frequencies the fundamental frequency is exactly $\gcd(f_1,f_2)$.

**Proof of 1.** $qT_1$ is a whole number of periods of $x_1$, and $pT_2=qT_1$ is a whole number of periods of $x_2$, so both terms, hence their sum, repeat after $T$.

**Proof of 2 (for $\cos\omega_1t+\cos\omega_2t$).** The sum equals its maximum value 2 at $t=0$. It equals 2 at any other time $t$ only if both cosines equal 1, i.e. $\omega_1t=2\pi m$ and $\omega_2t=2\pi n$ for integers $m,n$ with $t\ne0$. Dividing, $\omega_1/\omega_2=m/n$ would be rational — contradiction. So the value 2 is attained only at $t=0$. A periodic function attains each of its values infinitely often (at $t$, $t+T$, $t+2T,\dots$), so the sum is not periodic. $\square$

Such sums are called **almost periodic** (or quasi-periodic): they come arbitrarily close to repeating but never do exactly. In practice, any measured frequency ratio is rational to within instrument precision, but the theoretical distinction matters: Fourier *series* need exact periodicity, whereas the Fourier *transform* (Ch. 16) handles any combination of frequencies.

> [!example] Worked Example 3.16 — Fundamental periods of combinations
> Find the fundamental period, or show none exists.
>
> (a) $x(t)=\sin3t+\cos4.5t$. (b) $v(t)=\cos(100\pi t)+0.2\cos(300\pi t)+0.1\cos(500\pi t)$. (c) $y(t)=\cos t+\cos\sqrt2\,t$. (d) $p(t)=\cos^2(100\pi t)$.
>
> **(a)** $\omega_1=3$, $\omega_2=4.5$; ratio $3/4.5=2/3$, rational. The common angular frequency is $\gcd(3,4.5)=1.5$ rad/s ($3=2\times1.5$, $4.5=3\times1.5$). So $T_0=2\pi/1.5=4\pi/3\approx4.189$ s. Check: in $4\pi/3$ s the first term completes $3\times\frac{4\pi}{3}/2\pi=2$ cycles and the second $4.5\times\frac{4\pi}{3}/2\pi=3$ cycles. ✓
>
> **(b)** Frequencies $50,150,250$ Hz; $\gcd=50$ Hz, so $T_0=20$ ms. This is a 50 Hz waveform with 3rd and 5th harmonic distortion (§3.5, §3.7).
>
> **(c)** $1/\sqrt2$ is irrational: **not periodic**.
>
> **(d)** $\cos^2(100\pi t)=\tfrac12+\tfrac12\cos(200\pi t)$: fundamental period $2\pi/(200\pi)=10$ ms — half the 20 ms period of $\cos(100\pi t)$. This is the instantaneous power in a resistor fed from 50 Hz: it pulses at 100 Hz.

### Spatial periodicity: wavelength and wavenumber

Signals on transmission lines, in waveguides and in free space vary in both time and position. A sinusoidal wave travelling in the $+z$ direction is

$$
u(z,t)=A\cos(\omega t-kz+\phi).
$$

At a fixed position it is a sinusoid in time with angular frequency $\omega$; at a fixed instant it is a sinusoid in space with **wavenumber** (spatial angular frequency) $k$ in rad/m and spatial period, the **wavelength**,

$$
\lambda=\frac{2\pi}{k}.
$$

A point of constant phase satisfies $\omega t-kz=\text{const}$, so it moves at the **phase velocity** $v_p=dz/dt=\omega/k$. Combining,

$$
v_p=\frac{\omega}{k}=f\lambda,\qquad \lambda=\frac{v_p}{f}.
$$

The sign convention is worth memorising: $\omega t-kz$ travels towards $+z$; $\omega t+kz$ travels towards $-z$.

> [!example] Worked Example 3.17 — Wavelengths and electrical length
> (a) Find the free-space wavelength ($v_p=c=3.00\times10^8$ m/s) at 50 Hz and at 2.4 GHz. (b) A PCB trace has $v_p=0.5c$. How long is a trace that is a quarter-wavelength at 2.4 GHz, and what phase shift does a 10 mm trace introduce at that frequency?
>
> **(a)** At 50 Hz: $\lambda=3\times10^8/50=6000$ km. At 2.4 GHz: $\lambda=0.125$ m $=12.5$ cm.
>
> **(b)** On the trace $\lambda=0.5\times12.5=6.25$ cm, so $\lambda/4=15.6$ mm. A 10 mm trace has $kz=2\pi(10/62.5)=1.005$ rad $=57.6^\circ$ of phase lag.
>
> **Interpretation.** A circuit is "electrically small" — phase differences across it negligible, lumped-element analysis valid — when its dimensions are a small fraction of $\lambda$. Mains wiring is tiny compared with 6000 km; a 10 mm PCB trace at 2.4 GHz is not small at all.

### Discrete-time sinusoids

When a sinusoid $\cos(\omega t+\phi)$ is sampled every $T_s$ seconds, the samples are

$$
x[n]=\cos(\omega nT_s+\phi)=\cos(\Omega n+\phi),\qquad \Omega=\omega T_s=2\pi\frac{f}{f_s},
$$

where $f_s=1/T_s$ is the sampling rate and $\Omega$ is the **digital (normalised) angular frequency** in rad/sample. Discrete-time sinusoids behave differently from continuous ones in two important ways.

> [!theorem] Periodicity of a discrete-time sinusoid
> $x[n]=\cos(\Omega n+\phi)$ is periodic in $n$ if and only if $\Omega/(2\pi)$ is rational. If $\Omega/(2\pi)=m/N$ in lowest terms ($N>0$), the fundamental period is $N$ samples.

**Proof.** $x[n+N]=x[n]$ for all $n$ requires $\Omega N=2\pi m$ for some integer $m$, i.e. $\Omega/2\pi=m/N$. The smallest such positive $N$ is the denominator in lowest terms. $\square$

> [!theorem] Frequency ambiguity (aliasing)
> $\cos\big((\Omega+2\pi k)n+\phi\big)=\cos(\Omega n+\phi)$ for every integer $k$, and $\cos(-\Omega n-\phi)=\cos(\Omega n+\phi)$. Hence discrete-time frequencies are only distinguishable within $0\le\Omega\le\pi$; in terms of the original analogue frequency, tones at $f$ and $\lvert f-kf_s\rvert$ produce identical samples.

The highest distinguishable frequency, $\Omega=\pi$ (i.e. $f=f_s/2$), is the **Nyquist frequency**; sampling theory proper is in Ch. 16–17.

> [!example] Worked Example 3.18 — Periods and aliases of sampled sinusoids
> (a) Find the fundamental periods of $\cos(\pi n/6)$, $\cos(3\pi n/8)$ and $\cos(0.5n)$.
> (b) A 1 kHz tone and a 7 kHz tone are each sampled at $f_s=8$ kHz. Compare the samples.
>
> **(a)** $\frac{\pi/6}{2\pi}=\frac1{12}$: period 12 samples. $\frac{3\pi/8}{2\pi}=\frac3{16}$: period 16 samples (it completes 3 cycles of the underlying cosine in 16 samples). $\frac{0.5}{2\pi}=\frac1{4\pi}$ is irrational: **not periodic**, even though the continuous-time $\cos(0.5t)$ is periodic with period $4\pi$.
>
> **(b)** $\Omega_1=2\pi(1/8)=\pi/4$; $\Omega_2=2\pi(7/8)=7\pi/4=2\pi-\pi/4$. Since $\cos(7\pi n/4)=\cos(2\pi n-\pi n/4)=\cos(\pi n/4)$, the two sample sequences are **identical**: the 7 kHz tone aliases to 1 kHz ($\lvert7-8\rvert=1$ kHz). Both have period 8 samples.
>
> **Check:** computing $\cos(2\pi\cdot7000\,n/8000)$ and $\cos(2\pi\cdot1000\,n/8000)$ for $n=0,\dots,11$ gives identical values to machine precision. ✓

> [!warning] Common Pitfall — Hz versus rad/s in formulas
> Reactance is $X_C=1/(\omega C)=1/(2\pi fC)$, not $1/(fC)$; a resonant frequency $\omega_0=1/\sqrt{LC}$ is in rad/s, so $f_0=1/(2\pi\sqrt{LC})$; filter "corner frequencies" on datasheets are usually in Hz while transfer functions $H(j\omega)$ use rad/s. Because the radian is dimensionless, no unit check will catch a missing $2\pi$ — always label which one a number is.

### Instantaneous frequency

When the phase of a signal is not a linear function of time — as in frequency modulation or a radar "chirp" — we generalise frequency to the rate of change of phase.

> [!definition] Instantaneous angular frequency
> For $x(t)=A\cos\Theta(t)$, the **instantaneous angular frequency** is $\omega_i(t)=\dfrac{d\Theta}{dt}$ and the instantaneous frequency is $f_i(t)=\dfrac{1}{2\pi}\dfrac{d\Theta}{dt}$.

For a pure sinusoid $\Theta=\omega t+\phi$ and $\omega_i=\omega$, as it should. For the linear chirp $\Theta(t)=2\pi\big(f_0t+\tfrac12\beta t^2\big)$, $f_i(t)=f_0+\beta t$: the frequency sweeps linearly at $\beta$ Hz/s. With $f_0=1$ kHz and $\beta=10$ kHz/s, at $t=0.5$ s the instantaneous frequency is $1+10(0.5)=6$ kHz. Note that the "frequency" of a chirp is *not* the coefficient of $t$ inside the cosine at a given instant: writing $\cos(2\pi f(t)\,t)$ with $f(t)=f_0+\beta t$ would give $f_i=f_0+2\beta t$, twice the intended sweep rate.

```desmos-graph
left=0; right=6; bottom=-1.5; top=1.5
---
y=\cos(2\pi(0.5x+0.25x^{2}))|blue
```
*A linear chirp $\cos\big(2\pi(0.5t+0.25t^2)\big)$: the instantaneous frequency $0.5+0.5t$ Hz rises steadily, so the oscillations crowd together as $t$ increases.*

## 3.5 Harmonic functions

### Two meanings of "harmonic"

The word *harmonic* carries two related meanings in engineering mathematics, and this section is about both:

1. A **harmonic function of time** (a *harmonic oscillation*) is a sinusoid $A\cos(\omega t+\phi)$ — the motion of an ideal mass on a spring, the voltage across an ideal LC tank, the projection of a uniformly rotating vector. It is characterised as the solution of the simple-harmonic equation $\ddot x+\omega^2x=0$.
2. The **harmonics** of a periodic signal with fundamental frequency $f_0$ are the sinusoids at the integer multiples $nf_0$. A distorted 50 Hz current "contains a 5th harmonic" if it has a 250 Hz sinusoidal component.

(A third, unrelated, meaning — a function of several variables satisfying Laplace's equation $\nabla^2u=0$ — appears in [Ch. 13](../part-4-differential-equations/13-partial-differential-equations.md) and [Ch. 18](../part-5-transforms/18-complex-analysis.md). The shared word comes from the fact that Laplace's equation in a disc is solved by sums of $r^n\cos n\phi$ — harmonics again.)

### Simple harmonic motion

> [!definition] Simple harmonic oscillator
> A quantity $x(t)$ undergoes **simple harmonic motion** (SHM) with natural angular frequency $\omega_0>0$ if
> $$
> \ddot x+\omega_0^2x=0 .
> $$
> The acceleration is proportional to the displacement and directed back towards zero.

> [!theorem] Solutions of the simple-harmonic equation
> Every solution of $\ddot x+\omega_0^2x=0$ can be written in each of the equivalent forms
> $$
> x(t)=c_1\cos\omega_0t+c_2\sin\omega_0t=A\cos(\omega_0t+\phi)=\operatorname{Re}\{\mathbf{X}e^{j\omega_0t}\},
> $$
> with $\mathbf{X}=Ae^{j\phi}=c_1-jc_2$. Given initial values $x(0)=x_0$ and $\dot x(0)=v_0$,
> $$
> c_1=x_0,\quad c_2=\frac{v_0}{\omega_0},\qquad A=\sqrt{x_0^2+\Big(\frac{v_0}{\omega_0}\Big)^2},\quad \phi=\operatorname{atan2}\!\Big(-\frac{v_0}{\omega_0},\ x_0\Big).
> $$

**Derivation.** *That these are solutions* follows from §3.2: differentiating $A\cos(\omega_0t+\phi)$ twice multiplies it by $-\omega_0^2$. *That they are all the solutions* follows from an energy argument. Suppose $x$ is any solution and let $y=x-(c_1\cos\omega_0t+c_2\sin\omega_0t)$ with $c_1=x(0)$, $c_2=\dot x(0)/\omega_0$. Then $y$ solves the same equation with $y(0)=\dot y(0)=0$. Multiply $\ddot y+\omega_0^2y=0$ by $\dot y$:

$$
\frac{d}{dt}\Big[\tfrac12\dot y^2+\tfrac12\omega_0^2y^2\Big]=\dot y\,\ddot y+\omega_0^2y\,\dot y=0 .
$$

So $E=\tfrac12\dot y^2+\tfrac12\omega_0^2y^2$ is constant; it is zero at $t=0$, hence zero for all $t$, so $y\equiv0$. The initial-value formulas come from $x(0)=c_1$, $\dot x(0)=\omega_0c_2$ and the amplitude–phase conversion of §3.2. The general theory of such equations is in [Ch. 11](../part-4-differential-equations/11-ordinary-differential-equations.md)–[12](../part-4-differential-equations/12-linear-differential-equations.md). $\square$

The quantity $E$ in the proof is (up to a constant factor) the physical energy, which in SHM sloshes between two forms — kinetic and potential, or electric and magnetic — while its total stays fixed.

> [!example] Worked Example 3.19 — Fitting initial conditions
> Solve $\ddot x+9x=0$ with $x(0)=2$, $\dot x(0)=-6$ and write the answer in amplitude–phase form.
>
> **Solution.** $\omega_0=3$. $c_1=2$, $c_2=-6/3=-2$, so $x=2\cos3t-2\sin3t$.
> Amplitude $A=\sqrt{4+4}=2\sqrt2=2.828$; phase $\phi=\operatorname{atan2}(2,2)=45^\circ$.
> $x(t)=2\sqrt2\cos(3t+45^\circ)$.
>
> **Check:** $x(0)=2\sqrt2\cos45^\circ=2$ ✓; $\dot x(0)=-3(2\sqrt2)\sin45^\circ=-6$ ✓.

### The LC tank: an electrical harmonic oscillator

For an ideal inductor $L$ connected across an ideal capacitor $C$, the capacitor voltage $v$ and the current $i$ (taken as the current leaving the capacitor's positive plate and flowing through $L$) satisfy $i=-C\,\dot v$ and $v=L\,\dot i$. Eliminating $i$:

$$
\ddot v+\frac{1}{LC}v=0,\qquad \omega_0=\frac1{\sqrt{LC}},\qquad f_0=\frac{1}{2\pi\sqrt{LC}} .
$$

The stored energy $W=\tfrac12Cv^2+\tfrac12Li^2$ is constant. With $v=V\cos(\omega_0t+\phi)$, the capacitor energy is $\tfrac12CV^2\cos^2(\cdot)=\tfrac14CV^2[1+\cos(2\omega_0t+2\phi)]$: **each stored energy oscillates at twice the natural frequency**, the two oscillating parts cancelling exactly.

> [!example] Worked Example 3.20 — Free oscillation of an LC tank
> $L=10$ mH, $C=100$ nF. At $t=0$ the capacitor voltage is $5$ V and the current through the inductor is $i(0)=20$ mA (so $C\dot v(0)=-20$ mA). Find $v(t)$, the peak voltage, and verify energy conservation.
>
> **Natural frequency.** $\omega_0=1/\sqrt{10^{-2}\times10^{-7}}=1/\sqrt{10^{-9}}=31\,623$ rad/s; $f_0=5.033$ kHz.
>
> **Initial slope.** $\dot v(0)=-0.020/10^{-7}=-2.0\times10^5$ V/s, so $v_0/\omega_0=-6.325$ V.
>
> **Amplitude and phase.** $V=\sqrt{5^2+6.325^2}=\sqrt{65}=8.062$ V; $\phi=\operatorname{atan2}(6.325,5)=51.67^\circ$.
>
> **Result.** $v(t)=8.062\cos(31\,623\,t+51.67^\circ)$ V.
>
> **Energy check.** Initially $\tfrac12(10^{-7})(25)+\tfrac12(10^{-2})(0.02)^2=1.25+2.00=3.25\ \mu$J. At a voltage peak the current is zero and $W=\tfrac12(10^{-7})(65)=3.25\ \mu$J. ✓

### The rotating phasor picture

Euler's formula gives a geometric picture of every harmonic oscillation. The complex number $\mathbf{X}e^{j\omega t}=Ae^{j(\omega t+\phi)}$ moves anticlockwise around a circle of radius $A$, completing one revolution every $T=2\pi/\omega$ seconds. Its projection on the real axis is $A\cos(\omega t+\phi)$ and on the imaginary axis $A\sin(\omega t+\phi)$:

- **amplitude** = radius of the circle,
- **angular frequency** = angular speed of rotation,
- **phase** = angle of the vector at $t=0$,
- **phase difference** between two same-frequency signals = the fixed angle between their rotating vectors.

Since both vectors rotate at the same rate, the angle between them never changes; this is why a phasor diagram can be drawn as a static picture "frozen" at $t=0$.

A real cosine can also be written as the sum of *two* counter-rotating vectors:

$$
A\cos(\omega t+\phi)=\frac A2e^{j(\omega t+\phi)}+\frac A2e^{-j(\omega t+\phi)} .
$$

The second term rotates clockwise — it has *negative frequency* $-\omega$. The imaginary parts cancel and the real parts add. This two-sided view is the basis of the complex Fourier series ([Ch. 15](../part-5-transforms/15-fourier-series.md)) and the Fourier transform ([Ch. 16](../part-5-transforms/16-fourier-transforms.md)), in which every real sinusoid appears as a pair of spectral lines at $\pm\omega$ with conjugate complex amplitudes $\tfrac A2e^{\pm j\phi}$.

### Damped harmonic oscillations

Real oscillators lose energy. Adding a resistance to the LC tank, or friction to the spring, gives the damped equation

$$
\ddot x+2\alpha\dot x+\omega_0^2x=0,\qquad \alpha>0,
$$

whose solution for light damping ($\alpha<\omega_0$) is derived in [Ch. 12](../part-4-differential-equations/12-linear-differential-equations.md):

> [!theorem] Underdamped oscillation
> For $0<\alpha<\omega_0$ every solution has the form
> $$
> x(t)=Ae^{-\alpha t}\cos(\omega_dt+\phi),\qquad \omega_d=\sqrt{\omega_0^2-\alpha^2} .
> $$
> It is a sinusoid of angular frequency $\omega_d$ (the **damped natural frequency**) whose amplitude decays inside the **envelope** $\pm Ae^{-\alpha t}$. The **quality factor** is $Q=\omega_0/(2\alpha)$.

The verification is a direct (if tedious) differentiation; more elegantly, $x=\operatorname{Re}\{\mathbf{X}e^{st}\}$ with $s=-\alpha+j\omega_d$, and $s$ is a root of $s^2+2\alpha s+\omega_0^2=0$.

Such a signal is **not periodic** — it never repeats exactly — but its zero crossings are evenly spaced, half a *pseudo-period* $T_d=2\pi/\omega_d$ apart, and the ratio of successive peaks is the constant $e^{-\alpha T_d}$. The natural logarithm of that ratio, $\delta=\alpha T_d$, is the **logarithmic decrement**, a standard way to measure damping from an oscilloscope trace.

```desmos-graph
left=0; right=10; bottom=-1.3; top=1.3
---
y=e^{-0.3x}\cos(4x)|blue
y=e^{-0.3x}|red|dashed
y=-e^{-0.3x}|red|dashed
```
*An underdamped oscillation $e^{-0.3t}\cos4t$ (blue) confined within its exponential envelope $\pm e^{-0.3t}$ (red).*

> [!example] Worked Example 3.21 — Ringing of a series RLC circuit
> A series RLC circuit with $R=20\ \Omega$, $L=10$ mH, $C=100$ nF rings after a disturbance. For a series circuit $\alpha=R/(2L)$ and $\omega_0=1/\sqrt{LC}$. Find $\omega_d$, $Q$, the logarithmic decrement, and the fraction of the initial amplitude left after 10 cycles.
>
> **Parameters.** $\alpha=20/(0.02)=1000$ s$^{-1}$; $\omega_0=31\,623$ rad/s; $\omega_d=\sqrt{31\,623^2-1000^2}=31\,607$ rad/s (only 0.05 % below $\omega_0$).
>
> **Quality factor.** $Q=31\,623/2000=15.8$.
>
> **Decrement.** $T_d=2\pi/\omega_d=198.8\ \mu$s; $\delta=\alpha T_d=0.1988$, so each peak is $e^{-0.1988}=0.820$ of the previous one.
>
> **After 10 cycles.** $e^{-10\delta}=e^{-1.988}=0.137$: about 14 % of the initial amplitude remains.
>
> **Check:** a useful rule of thumb is that the amplitude falls by $e^{-\pi}\approx4.3\%$ in $Q$ cycles; here $Q\approx16$ cycles gives $e^{-16\times0.1988}=e^{-3.14}$. ✓

### Harmonics of a fundamental

> [!definition] Harmonics
> Given a fundamental frequency $f_0$ (angular frequency $\omega_0$), the **$n$-th harmonic** is a sinusoid at $nf_0$, $n=1,2,3,\dots$:
> $$
> h_n(t)=A_n\cos(n\omega_0t+\phi_n).
> $$
> The 1st harmonic is the fundamental itself; $n$ is the **harmonic number** (harmonic order). Harmonics with even $n$ are **even harmonics**, those with odd $n$ **odd harmonics**. A component at a frequency that is not an integer multiple of $f_0$ is an **interharmonic**; one at $f_0/m$ is a **subharmonic**.

Musicians call the $n$-th harmonic the $(n-1)$-th *overtone*; engineers almost always count harmonics, so the "3rd harmonic" of 50 Hz is 150 Hz.

> [!theorem] A sum of harmonics is periodic with the fundamental period
> Any finite sum
> $$
> x(t)=A_0+\sum_{n=1}^{N}A_n\cos(n\omega_0t+\phi_n)
> $$
> is periodic with period $T_0=2\pi/\omega_0$. If $A_1\neq0$ it is the fundamental period; in general the fundamental period is $T_0/g$, where $g$ is the greatest common divisor of the harmonic numbers $n$ with $A_n\neq0$.

**Proof.** In time $T_0$ the $n$-th harmonic advances its phase by $n\omega_0T_0=2\pi n$, a whole number of cycles, so every term repeats. The statement about $g$ follows from §3.4: the frequencies present are multiples of $gf_0$ and their greatest common divisor is exactly $gf_0$. $\square$

So, for example, $\cos2\omega_0t+\cos4\omega_0t$ repeats every $T_0/2$: it "contains only even harmonics of $f_0$", but is better described as the 1st and 2nd harmonics of $2f_0$.

> [!example] Worked Example 3.22 — Identifying harmonics
> A power-quality analyser on a 60 Hz feeder reports components at 60, 180, 300, 420, 540 and 210 Hz. Classify each.
>
> **Solution.** $180/60=3$, $300/60=5$, $420/60=7$, $540/60=9$: the 3rd, 5th, 7th and 9th harmonics (all odd — typical of loads with half-wave symmetric current, §3.8). $210/60=3.5$ is not an integer: an **interharmonic**, often produced by variable-frequency drives or arc furnaces, and not part of the 60 Hz harmonic series.

### Harmonic content versus waveform shape

The amplitudes $A_n$ determine the *power* in each harmonic (§3.7) but the phases $\phi_n$ strongly affect the *shape*. Adding a third harmonic of amplitude $\tfrac13$ to a fundamental gives a flattened, nearly square top if the third harmonic is in the right phase, and a peaked waveform if it is inverted:

```desmos-graph
left=0; right=6.3; bottom=-1.6; top=1.6
---
y=\sin(x)+\frac{1}{3}\sin(3x)|blue
y=\sin(x)-\frac{1}{3}\sin(3x)|red
y=\sin(x)|grey|dashed
```
*Fundamental plus a one-third-amplitude third harmonic: in phase (blue) the peaks flatten towards a square wave; inverted (red) the peaks sharpen. Both have exactly the same RMS value.*

> [!tip] Engineering Note — Harmonics in power systems and RF
> In power systems, nonlinear loads (rectifiers, switch-mode supplies, LED drivers) draw currents rich in odd harmonics; standards such as IEC 61000-3-2 limit each harmonic current by order. In RF systems, a power amplifier driven at $f_0$ radiates harmonics at $2f_0, 3f_0,\dots$, which regulations require to be filtered. In both cases the engineer works with the *list* of harmonic amplitudes, which is why the rest of this chapter develops the algebra of harmonics.

## 3.6 Trigonometric expansions

### Motivation

Linear systems are analysed one frequency at a time, so the natural currency of signal analysis is the *sum* of sinusoids. But physics keeps handing us *products* and *powers*: instantaneous power is the product $v(t)\,i(t)$; a mixer multiplies two signals; a diode or transistor responds to $v^2$, $v^3,\dots$; a modulator multiplies a carrier by a message. Trigonometric identities are the dictionary that translates products and powers back into sums, and the frequencies that appear in those sums are exactly the frequencies a spectrum analyser will show. This section collects the identities, derives them all from one source — Euler's formula — and shows how to apply them systematically.

### The master method: Euler's formula

From $e^{j\theta}=\cos\theta+j\sin\theta$ ([Ch. 2](02-complex-numbers.md)) and its conjugate,

$$
\cos\theta=\frac{e^{j\theta}+e^{-j\theta}}{2},\qquad \sin\theta=\frac{e^{j\theta}-e^{-j\theta}}{2j}.
$$

Every identity in this section can be proved by (i) replacing each $\cos$ and $\sin$ by these exponentials, (ii) multiplying out using $e^{ja}e^{jb}=e^{j(a+b)}$, and (iii) regrouping conjugate pairs back into cosines and sines. Products of exponentials *add* their exponents, which is the whole reason products of sinusoids become sinusoids at sum and difference frequencies.

### Compound-angle formulas

> [!theorem] Compound-angle (addition) formulas
> $$
> \begin{aligned}
> \cos(A\pm B)&=\cos A\cos B\mp\sin A\sin B,\\
> \sin(A\pm B)&=\sin A\cos B\pm\cos A\sin B,\\
> \tan(A\pm B)&=\frac{\tan A\pm\tan B}{1\mp\tan A\tan B}.
> \end{aligned}
> $$

**Proof.** $e^{j(A+B)}=e^{jA}e^{jB}$. Expand both sides:

$$
\cos(A+B)+j\sin(A+B)=(\cos A+j\sin A)(\cos B+j\sin B)=(\cos A\cos B-\sin A\sin B)+j(\sin A\cos B+\cos A\sin B).
$$

Equating real and imaginary parts gives the $+$ formulas; replacing $B$ by $-B$ and using evenness of $\cos$ and oddness of $\sin$ gives the $-$ formulas. Dividing the sine formula by the cosine formula, then numerator and denominator by $\cos A\cos B$, gives the tangent formula. $\square$

A purely geometric proof (rotating a unit vector by $A$ and then by $B$) is also standard; the rotation view is worth keeping in mind because it says that the compound-angle formula *is* the rule for multiplying by the rotation matrix $\begin{pmatrix}\cos B&-\sin B\\ \sin B&\cos B\end{pmatrix}$ of [Ch. 9](../part-3-linear-algebra/09-vector-spaces-and-linear-transformations.md).

> [!example] Worked Example 3.23 — Exact values from compound angles
> Find exact values of $\sin75^\circ$, $\cos15^\circ$ and $\tan75^\circ$.
>
> **Solution.** $75^\circ=45^\circ+30^\circ$:
> $$
> \sin75^\circ=\sin45^\circ\cos30^\circ+\cos45^\circ\sin30^\circ=\frac{\sqrt2}{2}\cdot\frac{\sqrt3}{2}+\frac{\sqrt2}{2}\cdot\frac12=\frac{\sqrt6+\sqrt2}{4}.
> $$
> $\cos15^\circ=\cos(90^\circ-75^\circ)=\sin75^\circ=\dfrac{\sqrt6+\sqrt2}{4}$.
> $$
> \tan75^\circ=\frac{\tan45^\circ+\tan30^\circ}{1-\tan45^\circ\tan30^\circ}=\frac{1+1/\sqrt3}{1-1/\sqrt3}=\frac{\sqrt3+1}{\sqrt3-1}=2+\sqrt3 .
> $$
> **Check:** $(\sqrt6+\sqrt2)/4=0.965926=\sin75^\circ$; $2+\sqrt3=3.73205=\tan75^\circ$. ✓

### Double-angle, half-angle and power-reduction formulas

Setting $B=A$ in the compound-angle formulas:

$$
\sin2A=2\sin A\cos A,\qquad \cos2A=\cos^2A-\sin^2A=2\cos^2A-1=1-2\sin^2A,\qquad \tan2A=\frac{2\tan A}{1-\tan^2A}.
$$

Solving the two alternative forms of $\cos2A$ for the squares gives the **power-reduction** formulas, which convert a square (a "power") into a constant plus a second harmonic:

$$
\cos^2A=\frac{1+\cos2A}{2},\qquad \sin^2A=\frac{1-\cos2A}{2},\qquad \sin A\cos A=\frac{\sin2A}{2}.
$$

Replacing $A$ by $A/2$ gives the **half-angle** formulas $\cos\frac A2=\pm\sqrt{\frac{1+\cos A}{2}}$, $\sin\frac A2=\pm\sqrt{\frac{1-\cos A}{2}}$, with the sign fixed by the quadrant of $A/2$, and

$$
\tan\frac A2=\frac{\sin A}{1+\cos A}=\frac{1-\cos A}{\sin A}.
$$

The **tangent half-angle substitution** $u=\tan(A/2)$ expresses everything rationally,

$$
\sin A=\frac{2u}{1+u^2},\qquad \cos A=\frac{1-u^2}{1+u^2},
$$

which turns any rational function of $\sin A$ and $\cos A$ into a rational function of $u$ — a standard integration technique ([Ch. 5](../part-2-calculus/05-integral-calculus.md)) and, in signal processing, the origin of the bilinear transform ([Ch. 17](../part-5-transforms/17-z-transforms.md)).

Applying the double-angle formula twice gives triple angles:

$$
\cos3A=4\cos^3A-3\cos A,\qquad \sin3A=3\sin A-4\sin^3A,
$$

and inverting these gives the cubic power-reduction formulas

$$
\cos^3A=\frac{3\cos A+\cos3A}{4},\qquad \sin^3A=\frac{3\sin A-\sin3A}{4}.
$$

**These are the identities that explain harmonic distortion.** If a device has a cubic term in its transfer characteristic, a pure input $\cos\omega t$ produces $\cos^3\omega t$, which contains a component at $3\omega$ (§3.7).

> [!example] Worked Example 3.24 — Reducing a fourth power
> Express $\sin^4\theta$ as a sum of cosines of multiple angles, and deduce its average value.
>
> **Solution.** Square the power-reduction formula:
> $$
> \sin^4\theta=\Big(\frac{1-\cos2\theta}{2}\Big)^2=\frac{1-2\cos2\theta+\cos^22\theta}{4}.
> $$
> Reduce $\cos^22\theta=\frac{1+\cos4\theta}{2}$:
> $$
> \sin^4\theta=\frac14-\frac12\cos2\theta+\frac{1+\cos4\theta}{8}=\frac38-\frac12\cos2\theta+\frac18\cos4\theta .
> $$
> The cosine terms average to zero over a period, so $\overline{\sin^4\theta}=3/8$.
>
> **Check:** at $\theta=\pi/2$: LHS $=1$; RHS $=\frac38+\frac12+\frac18=1$ ✓. At $\theta=\pi/4$: LHS $=\frac14$; RHS $=\frac38-0-\frac18=\frac14$ ✓. A numerical comparison at 50 points agrees to $10^{-16}$.

### Powers of cosine in general

> [!theorem] Harmonic expansion of $\cos^n\theta$
> For a positive integer $n$,
> $$
> \cos^n\theta=\frac1{2^n}\sum_{k=0}^{n}\binom nk\cos\big((n-2k)\theta\big).
> $$
> Grouping the terms $k$ and $n-k$ (which give equal cosines), $\cos^n\theta$ contains only the harmonics $n,\,n-2,\,n-4,\dots$ — all odd if $n$ is odd, and all even (including a constant) if $n$ is even. The constant term for even $n$ is $\binom{n}{n/2}/2^n$.

**Proof.** By the binomial theorem,

$$
\cos^n\theta=\frac{(e^{j\theta}+e^{-j\theta})^n}{2^n}=\frac1{2^n}\sum_{k=0}^n\binom nk e^{j(n-k)\theta}e^{-jk\theta}=\frac1{2^n}\sum_{k=0}^n\binom nk e^{j(n-2k)\theta}.
$$

The left side is real, so it equals the real part of the right, which replaces each exponential by a cosine. $\square$

The same argument with $(e^{j\theta}-e^{-j\theta})^n/(2j)^n$ handles $\sin^n\theta$: odd powers of sine give sines of odd multiples, even powers give cosines of even multiples.

### Multiple angles and Chebyshev polynomials

Going the other way, De Moivre's theorem $(\cos\theta+j\sin\theta)^n=\cos n\theta+j\sin n\theta$ ([Ch. 2](02-complex-numbers.md)) expresses $\cos n\theta$ as a polynomial in $\cos\theta$.

> [!theorem] Chebyshev polynomials
> For each $n\ge0$ there is a polynomial $T_n$ of degree $n$ with
> $$
> \cos n\theta=T_n(\cos\theta).
> $$
> They satisfy $T_0(x)=1$, $T_1(x)=x$ and the recurrence
> $$
> T_{n+1}(x)=2x\,T_n(x)-T_{n-1}(x),
> $$
> giving $T_2=2x^2-1$, $T_3=4x^3-3x$, $T_4=8x^4-8x^2+1$, $T_5=16x^5-20x^3+5x$.

**Proof of the recurrence.** Add the compound-angle formulas $\cos(n+1)\theta=\cos n\theta\cos\theta-\sin n\theta\sin\theta$ and $\cos(n-1)\theta=\cos n\theta\cos\theta+\sin n\theta\sin\theta$:

$$
\cos(n+1)\theta+\cos(n-1)\theta=2\cos\theta\cos n\theta .
$$

With $x=\cos\theta$ this is $T_{n+1}(x)+T_{n-1}(x)=2xT_n(x)$. Since $T_0$ and $T_1$ are polynomials, induction shows every $T_n$ is a polynomial of degree $n$ with leading coefficient $2^{n-1}$ ($n\ge1$). $\square$

Because $\lvert T_n(x)\rvert=\lvert\cos n\theta\rvert\le1$ on $-1\le x\le1$ while growing rapidly outside, Chebyshev polynomials are the basis of Chebyshev filters (equal-ripple passband) and of near-optimal polynomial approximation ([Ch. 22](../part-7-numerical-computational/22-numerical-methods.md)).

> [!example] Worked Example 3.25 — Two directions: $\cos5\theta$ and $\cos^5\theta$
> (a) Express $\cos5\theta$ as a polynomial in $c=\cos\theta$. (b) Express $\cos^5\theta$ as a sum of harmonics.
>
> **(a)** From the recurrence: $T_4=2x(4x^3-3x)-(2x^2-1)=8x^4-8x^2+1$, then $T_5=2x(8x^4-8x^2+1)-(4x^3-3x)=16x^5-20x^3+5x$. So $\cos5\theta=16c^5-20c^3+5c$.
>
> **(b)** With $n=5$: $\binom50,\dots,\binom55=1,5,10,10,5,1$ multiply $\cos5\theta,\cos3\theta,\cos\theta,\cos(-\theta),\cos(-3\theta),\cos(-5\theta)$. Pairing:
> $$
> \cos^5\theta=\frac{1}{32}\big(2\cos5\theta+10\cos3\theta+20\cos\theta\big)=\frac{10\cos\theta+5\cos3\theta+\cos5\theta}{16}.
> $$
> **Consistency check.** Substitute (b)'s $\cos5\theta=16\cos^5\theta-5\cos3\theta-10\cos\theta$ and $\cos3\theta=4c^3-3c$: $16c^5-20c^3+15c-10c=16c^5-20c^3+5c$, agreeing with (a). ✓ At $\theta=0$ both sides of (b) equal $16/16=1$. ✓

### Product-to-sum and sum-to-product formulas

Adding and subtracting the compound-angle formulas in pairs gives:

> [!theorem] Product-to-sum formulas
> $$
> \begin{aligned}
> \cos A\cos B&=\tfrac12\big[\cos(A-B)+\cos(A+B)\big],\\
> \sin A\sin B&=\tfrac12\big[\cos(A-B)-\cos(A+B)\big],\\
> \sin A\cos B&=\tfrac12\big[\sin(A+B)+\sin(A-B)\big].
> \end{aligned}
> $$

For example, the first comes from adding $\cos(A-B)=\cos A\cos B+\sin A\sin B$ and $\cos(A+B)=\cos A\cos B-\sin A\sin B$. Applied to signals, $A=\omega_1t$ and $B=\omega_2t$: **the product of two sinusoids is the sum of two sinusoids at the sum and difference frequencies, each with half the product amplitude.** This single fact underlies mixing, modulation, demodulation and the computation of average power.

Setting $A=(P+Q)/2$, $B=(P-Q)/2$ reverses the formulas:

> [!theorem] Sum-to-product formulas
> $$
> \begin{aligned}
> \cos P+\cos Q&=2\cos\tfrac{P+Q}2\cos\tfrac{P-Q}2, & \cos P-\cos Q&=-2\sin\tfrac{P+Q}2\sin\tfrac{P-Q}2,\\
> \sin P+\sin Q&=2\sin\tfrac{P+Q}2\cos\tfrac{P-Q}2, & \sin P-\sin Q&=2\cos\tfrac{P+Q}2\sin\tfrac{P-Q}2.
> \end{aligned}
> $$

Sum-to-product is the tool for equations (a product is zero iff a factor is) and for the beat phenomenon (§3.7).

> [!example] Worked Example 3.26 — Instantaneous and average power
> A load has $v(t)=325\cos\omega t$ V and $i(t)=10\cos(\omega t-30^\circ)$ A. Find the instantaneous power $p(t)=v\,i$, its average, and the frequency at which it pulsates.
>
> **Solution.** Product-to-sum with $A=\omega t$, $B=\omega t-30^\circ$:
> $$
> p(t)=3250\cdot\tfrac12\big[\cos30^\circ+\cos(2\omega t-30^\circ)\big]=1407+1625\cos(2\omega t-30^\circ)\ \text{W}.
> $$
> The average (real) power is $P=1407$ W and $p(t)$ pulsates at **twice** the supply frequency with amplitude 1625 W. Because $1625>1407$, $p$ is briefly negative each half-cycle: energy flows back from the load's reactive part to the source.
>
> **General result.** For $v=V_m\cos(\omega t+\phi_v)$, $i=I_m\cos(\omega t+\phi_i)$:
> $$
> P=\overline{p}=\tfrac12V_mI_m\cos(\phi_v-\phi_i)=V_{\text{rms}}I_{\text{rms}}\cos(\phi_v-\phi_i),
> $$
> the power-factor formula quoted in §3.3. Check: $\frac{325}{\sqrt2}\cdot\frac{10}{\sqrt2}\cos30^\circ=1625\times0.8660=1407$ W. ✓

> [!example] Worked Example 3.27 — Solving an equation with sum-to-product
> Solve $\sin3x+\sin x=\cos x$ on $[0,\pi]$.
>
> **Solution.** $\sin3x+\sin x=2\sin2x\cos x$. The equation becomes $2\sin2x\cos x-\cos x=0$, i.e. $\cos x\,(2\sin2x-1)=0$.
>
> - $\cos x=0$: $x=\pi/2$.
> - $\sin2x=\tfrac12$: $2x=\pi/6$ or $5\pi/6$ (within $[0,2\pi]$), so $x=\pi/12$ or $5\pi/12$.
>
> Solutions: $x=\pi/12,\ 5\pi/12,\ \pi/2$.
>
> **Check** $x=\pi/12=15^\circ$: $\sin45^\circ+\sin15^\circ=0.7071+0.2588=0.9659=\cos15^\circ$. ✓

### Sums of sinusoids with phases in arithmetic progression

A sum that appears in antenna arrays, multi-tap filters and the discrete Fourier transform is

$$
S=\sum_{k=0}^{N-1}e^{jk\delta}.
$$

> [!theorem] Geometric phasor sum
> For $\delta$ not a multiple of $2\pi$,
> $$
> \sum_{k=0}^{N-1}e^{jk\delta}=e^{j(N-1)\delta/2}\,\frac{\sin(N\delta/2)}{\sin(\delta/2)},
> $$
> and consequently
> $$
> \sum_{k=0}^{N-1}\cos(\theta+k\delta)=\frac{\sin(N\delta/2)}{\sin(\delta/2)}\cos\!\Big(\theta+\frac{(N-1)\delta}{2}\Big).
> $$

**Proof.** The geometric-series formula gives $S=\dfrac{e^{jN\delta}-1}{e^{j\delta}-1}$. Factor $e^{jN\delta/2}$ from the numerator and $e^{j\delta/2}$ from the denominator:

$$
S=\frac{e^{jN\delta/2}\big(e^{jN\delta/2}-e^{-jN\delta/2}\big)}{e^{j\delta/2}\big(e^{j\delta/2}-e^{-j\delta/2}\big)}=e^{j(N-1)\delta/2}\frac{2j\sin(N\delta/2)}{2j\sin(\delta/2)} .
$$

For the cosine sum, multiply by $e^{j\theta}$ and take real parts. $\square$

The magnitude $\lvert S\rvert=\lvert\sin(N\delta/2)/\sin(\delta/2)\rvert$ equals $N$ at $\delta=0$ (all phasors aligned) and vanishes when $N\delta/2$ is a nonzero multiple of $\pi$ (phasors spread evenly round a full circle, as in the three-phase theorem of §3.3).

```desmos-graph
left=-6.3; right=6.3; bottom=-1; top=9
---
y=\left|\frac{\sin(4x)}{\sin(x/2)}\right||blue
y=8|red|dashed
```
*The array factor $\lvert\sin(N\delta/2)/\sin(\delta/2)\rvert$ for $N=8$ against the inter-element phase $\delta$: a main lobe of height $N=8$ at $\delta=0$, nulls at multiples of $\delta=2\pi/8$, and smaller side lobes between them.*

> [!example] Worked Example 3.28 — Eight-element array
> Eight equal antenna elements are fed with a progressive phase $\delta$ between neighbours. Find the magnitude of the combined phasor (relative to one element) for $\delta=0^\circ$, $10^\circ$ and $45^\circ$.
>
> **Solution.** $\lvert S\rvert=\lvert\sin(4\delta)/\sin(\delta/2)\rvert$.
> - $\delta=0$: limit $=N=8$ (full coherent addition).
> - $\delta=10^\circ$: $\sin40^\circ/\sin5^\circ=0.6428/0.08716=7.375$.
> - $\delta=45^\circ$: $\sin180^\circ=0$: complete cancellation — the eight phasors are spaced $45^\circ$ apart round the full circle.
>
> **Check** for $\delta=10^\circ$: summing $\cos(10^\circ k)$ and $\sin(10^\circ k)$ for $k=0..7$ directly gives $6.041$ and $4.230$, magnitude $\sqrt{6.041^2+4.230^2}=7.375$. ✓

### Series expansions

The circular functions also have power-series expansions, valid for every real (indeed complex) $\theta$:

$$
\sin\theta=\sum_{m=0}^{\infty}\frac{(-1)^m\theta^{2m+1}}{(2m+1)!},\qquad
\cos\theta=\sum_{m=0}^{\infty}\frac{(-1)^m\theta^{2m}}{(2m)!}.
$$

They are derived as Maclaurin series in [Ch. 4](../part-2-calculus/04-differential-calculus.md); substituting them into $e^{j\theta}=\sum(j\theta)^k/k!$ and separating even and odd powers is one standard proof of Euler's formula. For computation they converge quickly only for small $\lvert\theta\rvert$, so libraries first reduce the argument into $[-\pi/4,\pi/4]$ using periodicity and the identities above, then apply a short polynomial — a direct engineering use of this section.

> [!warning] Common Pitfall — Squaring away the sign
> Half-angle formulas involve a square root, and solving $\sin^2x=k$ or $\cos2x=c$ yields *all* sign choices; check which roots are genuine. Conversely, identities like $\sqrt{1-\cos^2\theta}=\sin\theta$ are only true where $\sin\theta\ge0$. In an algorithm that computes $\sin\phi$ from $\cos\phi$, the lost sign is a lost lead/lag distinction — exactly the error that a power meter makes if it computes reactive power as $\sqrt{S^2-P^2}$ instead of $V_{\text{rms}}I_{\text{rms}}\sin\Delta\phi$.

## 3.7 Harmonic relationships

### Motivation

Section 3.5 introduced harmonics as sinusoids at integer multiples of a fundamental; §3.6 supplied the identities that turn products and powers into sums. This section puts the two together and answers the questions engineers actually ask about collections of sinusoids:

- How do the powers of several sinusoids combine? (Orthogonality.)
- Where do harmonics come from, and how big are they? (Nonlinearity.)
- What new frequencies appear when two signals pass through the same nonlinear device or are multiplied? (Intermodulation, mixing, modulation.)
- What does the sum of two *nearly* equal frequencies look like? (Beats.)
- How is "how distorted" quantified? (Total harmonic distortion.)
- How do harmonics behave in three-phase systems? (Phase sequence of harmonics.)

### Orthogonality of harmonics

> [!theorem] Orthogonality relations
> Let $\omega_0=2\pi/T_0$ and let $m,n$ be non-negative integers. Over one period,
> $$
> \begin{aligned}
> \int_0^{T_0}\cos m\omega_0t\,\cos n\omega_0t\,dt&=\begin{cases}0,&m\neq n,\\ T_0/2,&m=n\geq1,\\ T_0,&m=n=0,\end{cases}\\[4pt]
> \int_0^{T_0}\sin m\omega_0t\,\sin n\omega_0t\,dt&=\begin{cases}0,&m\neq n,\\ T_0/2,&m=n\geq1,\end{cases}\\[4pt]
> \int_0^{T_0}\sin m\omega_0t\,\cos n\omega_0t\,dt&=0\quad\text{for all }m,n.
> \end{aligned}
> $$
> The same results hold over any interval of length $T_0$.

**Proof.** By product-to-sum, $\cos m\omega_0t\cos n\omega_0t=\frac12[\cos(m-n)\omega_0t+\cos(m+n)\omega_0t]$. A cosine $\cos k\omega_0t$ with integer $k\neq0$ completes exactly $\lvert k\rvert$ cycles in $T_0$, so its integral is zero; for $k=0$ it is the constant 1 with integral $T_0$. If $m\neq n$ both $m-n$ and $m+n$ are nonzero and the integral vanishes. If $m=n\geq1$, only $\cos0=1$ survives, giving $\frac12T_0$; if $m=n=0$ the integrand is 1. The sine–sine case is identical with $\frac12[\cos(m-n)\omega_0t-\cos(m+n)\omega_0t]$. For sine–cosine, $\frac12[\sin(m+n)\omega_0t+\sin(m-n)\omega_0t]$ integrates to zero in every case (including $\sin0=0$). $\square$

In the language of [Ch. 9](../part-3-linear-algebra/09-vector-spaces-and-linear-transformations.md), the functions $1,\cos\omega_0t,\sin\omega_0t,\cos2\omega_0t,\dots$ are mutually **orthogonal** with respect to the inner product $\langle x,y\rangle=\int_0^{T_0}x(t)y(t)\,dt$. That is the structural fact behind Fourier series: it is what lets us *extract* each harmonic's coefficient from a waveform by a single integral (§3.8, Ch. 15).

In complex form the relations are even more compact:

$$
\frac1{T_0}\int_0^{T_0}e^{jm\omega_0t}\,e^{-jn\omega_0t}\,dt=\begin{cases}1,&m=n,\\0,&m\neq n.\end{cases}
$$

### Power of a sum of harmonics

> [!theorem] Mean-square value of a harmonic sum
> If $x(t)=A_0+\sum_{n=1}^NA_n\cos(n\omega_0t+\phi_n)$, then
> $$
> \overline{x^2}=\frac1{T_0}\int_0^{T_0}x^2\,dt=A_0^2+\sum_{n=1}^N\frac{A_n^2}{2},\qquad X_{\text{rms}}=\sqrt{A_0^2+\sum_{n=1}^N\Big(\frac{A_n}{\sqrt2}\Big)^2}.
> $$
> The mean-square value — the average power in a 1 $\Omega$ resistor — is the **sum of the powers of the individual harmonics**, independent of their phases.

**Proof.** Squaring $x$ produces squared terms and cross terms. Each cross term $2A_mA_n\cos(m\omega_0t+\phi_m)\cos(n\omega_0t+\phi_n)$ with $m\neq n$ is, after expanding the phases with the compound-angle formula, a combination of the orthogonal products above, so it averages to zero; the cross terms with the constant average to zero because each harmonic does. Each squared harmonic averages to $A_n^2/2$ (§3.2). $\square$

This is the finite version of **Parseval's theorem** (Ch. 15). Two warnings: the result needs *distinct* frequencies (two components at the same frequency must first be combined with phasors, §3.3, because their cross term does *not* vanish); and it concerns power, not peak value — the peak of $x$ depends strongly on the phases (recall the two waveforms in §3.5 with identical RMS but different peaks).

> [!example] Worked Example 3.29 — RMS of a multi-harmonic waveform
> Find the RMS value of $v(t)=10+8\cos\omega t+6\cos(3\omega t+40^\circ)$ V and the average power it delivers to $R=5\ \Omega$.
>
> **Solution.** $V_{\text{rms}}^2=10^2+\frac{8^2}{2}+\frac{6^2}{2}=100+32+18=150$ V$^2$, so $V_{\text{rms}}=12.25$ V. Power $=150/5=30$ W, of which 20 W is DC, 6.4 W fundamental, 3.6 W third harmonic. The $40^\circ$ phase plays no role.
>
> **Check:** numerical averaging of $v^2$ over a period on 200 000 points gives $V_{\text{rms}}=12.2474$ V. ✓

### Nonlinearity generates harmonics

A linear system cannot create new frequencies. A nonlinear one — a diode, a transistor driven hard, a saturating transformer core, an amplifier clipping at its supply rails — can and does. Many weakly nonlinear devices are well modelled by a polynomial (the first few terms of a Taylor expansion about the operating point, Ch. 4):

$$
y=a_1x+a_2x^2+a_3x^3 .
$$

Drive it with $x=A\cos\omega t$ and reduce the powers with §3.6:

$$
x^2=\frac{A^2}{2}+\frac{A^2}{2}\cos2\omega t,\qquad x^3=\frac{3A^3}{4}\cos\omega t+\frac{A^3}{4}\cos3\omega t.
$$

> [!theorem] Harmonics from a cubic nonlinearity
> For $y=a_1x+a_2x^2+a_3x^3$ and $x=A\cos\omega t$,
> $$
> y=\underbrace{\frac{a_2A^2}{2}}_{\text{DC}}+\underbrace{\Big(a_1A+\frac{3a_3A^3}{4}\Big)}_{\text{fundamental}}\cos\omega t+\underbrace{\frac{a_2A^2}{2}}_{\text{2nd}}\cos2\omega t+\underbrace{\frac{a_3A^3}{4}}_{\text{3rd}}\cos3\omega t .
> $$

Four features generalise to higher-order polynomials:

1. A term $x^k$ produces harmonics up to the $k$-th, of the same parity as $k$ (the $\cos^n$ theorem of §3.6).
2. **Even-order terms** produce even harmonics and a **DC shift** — the principle of rectification and of envelope detection.
3. **Odd-order terms** produce odd harmonics and also change the fundamental: if $a_3<0$ (a compressive, saturating device) the fundamental gain falls as $A$ grows — **gain compression**.
4. The $k$-th harmonic grows as $A^k$: halving the drive level (a 6 dB reduction) lowers the fundamental by 6 dB, the 2nd harmonic by 12 dB and the 3rd by 18 dB, so *relative to the fundamental* the 2nd harmonic improves by 6 dB and the 3rd by 12 dB. Plotting harmonic levels against drive level on log axes gives straight lines of slope $k$ — the basis of intercept-point specifications.

A device with odd symmetry, $y(-x)=-y(x)$ (a symmetric clipper, a push–pull amplifier), has $a_2=a_4=\dots=0$ and generates **only odd harmonics**; this is the waveform-symmetry result of §3.8 seen from the device side.

> [!example] Worked Example 3.30 — Harmonic distortion of a weakly nonlinear amplifier
> An amplifier has $y=x+0.1x^2+0.05x^3$ (volts) and is driven by $x=2\cos\omega t$ V. Find each output component and the second- and third-harmonic distortion relative to the output fundamental.
>
> **Solution.** With $A=2$, $a_1=1$, $a_2=0.1$, $a_3=0.05$:
> - DC: $a_2A^2/2=0.1(4)/2=0.20$ V,
> - fundamental: $a_1A+\frac34a_3A^3=2+\frac34(0.05)(8)=2.30$ V,
> - 2nd harmonic: $a_2A^2/2=0.20$ V,
> - 3rd harmonic: $a_3A^3/4=0.05(8)/4=0.10$ V.
>
> $\mathrm{HD}_2=0.20/2.30=8.7\%$, $\mathrm{HD}_3=0.10/2.30=4.3\%$.
>
> **Check:** at $\omega t=0$, $x=2$ and $y=2+0.4+0.4=2.8$; the harmonic sum gives $0.2+2.3+0.2+0.1=2.8$. ✓ Numerical projection onto $\cos n\omega t$ also returns $0.2,\ 2.3,\ 0.2,\ 0.1$.

### Two tones in a nonlinearity: intermodulation

Now drive the same nonlinearity with two tones, $x=A\cos\omega_1t+A\cos\omega_2t$. The square term gives, by product-to-sum,

$$
x^2=A^2+\frac{A^2}{2}\cos2\omega_1t+\frac{A^2}{2}\cos2\omega_2t+A^2\cos(\omega_1+\omega_2)t+A^2\cos(\omega_1-\omega_2)t .
$$

Besides the harmonics $2\omega_1$, $2\omega_2$ there are new **intermodulation (IM) products** at $\omega_1\pm\omega_2$. The cube term, expanded the same way, produces components at

$$
\omega_1,\ \omega_2,\ 3\omega_1,\ 3\omega_2,\ 2\omega_1\pm\omega_2,\ 2\omega_2\pm\omega_1 .
$$

In general an order-$k$ nonlinearity produces every frequency $\lvert p\omega_1+q\omega_2\rvert$ with integers $p,q$ satisfying $\lvert p\rvert+\lvert q\rvert\le k$ and of the same parity as $k$. The **order** of a product is $\lvert p\rvert+\lvert q\rvert$.

> [!example] Worked Example 3.31 — Third-order intermodulation
> Two equal carriers at $f_1=100$ MHz and $f_2=101$ MHz pass through an amplifier with second- and third-order nonlinearity. List the output frequencies and identify those that fall near the original band.
>
> **Second order** ($\lvert p\rvert+\lvert q\rvert=2$): $2f_1=200$, $2f_2=202$, $f_1+f_2=201$, $f_2-f_1=1$ MHz, plus DC.
>
> **Third order** ($\lvert p\rvert+\lvert q\rvert=3$): $3f_1=300$, $3f_2=303$, $2f_1+f_2=301$, $2f_2+f_1=302$, $2f_1-f_2=99$, $2f_2-f_1=102$ MHz, plus contributions at $f_1$ and $f_2$ themselves.
>
> **Interpretation.** Every second-order product and every harmonic lies far outside a 99–102 MHz band and can be filtered out. The **third-order IM products at 99 MHz and 102 MHz** sit right next to the wanted signals, cannot be filtered, and land on neighbouring channels. This is why RF amplifiers and mixers are specified by their third-order intercept point rather than by their harmonic distortion.

### Mixing and modulation

A **mixer** multiplies two signals. With product-to-sum,

$$
\cos\omega_1t\cdot\cos\omega_2t=\tfrac12\cos(\omega_1-\omega_2)t+\tfrac12\cos(\omega_1+\omega_2)t,
$$

so multiplication *shifts* frequencies: a signal at $\omega_1$ mixed with a local oscillator at $\omega_2$ reappears at the **difference** (intermediate) frequency $\lvert\omega_1-\omega_2\rvert$ and the **sum** frequency $\omega_1+\omega_2$. A filter then keeps the one that is wanted. This is the superheterodyne principle behind essentially every radio receiver.

**Amplitude modulation** is mixing of a carrier with $(1+\text{message})$:

$$
s(t)=A_c\big[1+m\cos\omega_mt\big]\cos\omega_ct=A_c\cos\omega_ct+\frac{mA_c}{2}\cos(\omega_c-\omega_m)t+\frac{mA_c}{2}\cos(\omega_c+\omega_m)t .
$$

The spectrum is a carrier with two **sidebands** at $\omega_c\pm\omega_m$; $0<m\le1$ is the **modulation index**.

> [!example] Worked Example 3.32 — Spectrum and power of an AM signal
> A 1 MHz carrier of amplitude 10 V is amplitude-modulated by a 5 kHz tone with $m=0.5$. Find the spectral components and the fraction of the total power in the sidebands.
>
> **Components.** Carrier: 1000 kHz, 10 V. Sidebands: 995 kHz and 1005 kHz, each $mA_c/2=2.5$ V.
>
> **Power** (into 1 $\Omega$, using orthogonality): carrier $10^2/2=50$ W; each sideband $2.5^2/2=3.125$ W; total $56.25$ W. Sideband fraction $6.25/56.25=11.1\%$.
>
> **General formula.** Sideband fraction $=\dfrac{m^2/2}{1+m^2/2}$; at $m=0.5$ this is $0.125/1.125=11.1\%$ ✓, and even at $m=1$ only one third of the power carries information — the reason for suppressed-carrier and single-sideband systems.

### Beats

The sum of two equal-amplitude sinusoids at *close* frequencies is best understood through sum-to-product:

$$
\cos\omega_1t+\cos\omega_2t=2\cos\Big(\frac{\omega_1-\omega_2}{2}t\Big)\cos\Big(\frac{\omega_1+\omega_2}{2}t\Big).
$$

When $\lvert\omega_1-\omega_2\rvert\ll\omega_1$, this is a fast oscillation at the **average** frequency $\bar\omega=(\omega_1+\omega_2)/2$, whose amplitude is modulated by the slow envelope $\lvert2\cos(\Delta\omega\,t/2)\rvert$, with $\Delta\omega=\omega_1-\omega_2$. The envelope magnitude peaks twice per cycle of $\cos(\Delta\omega t/2)$, so the loudness (or the reading of an envelope detector) pulses at the **beat frequency**

$$
f_{\text{beat}}=\lvert f_1-f_2\rvert .
$$

Physically: the two phasors rotate at slightly different speeds; once per beat period they line up (constructive, amplitude $2A$) and once they oppose (destructive, amplitude 0).

```desmos-graph
left=0; right=4; bottom=-2.3; top=2.3
---
y=\cos(2\pi\cdot 10x)+\cos(2\pi\cdot 11x)|blue
y=2\cos(\pi x)|red|dashed
y=-2\cos(\pi x)|red|dashed
```
*Beats: $\cos(2\pi\,10t)+\cos(2\pi\,11t)$ (blue) inside the envelope $\pm2\cos(\pi t)$ (red). The envelope magnitude peaks once per second, the beat frequency $\lvert11-10\rvert=1$ Hz.*

> [!example] Worked Example 3.33 — Tuning by beats, and the difference with AM
> (a) A reference tone at 440 Hz and a slightly sharp tone at 444 Hz sound together. What is heard?
> (b) Compare $\cos\omega_1t+\cos\omega_2t$ with the AM signal of Example 3.32. Why do they look similar?
>
> **(a)** A 442 Hz tone whose loudness pulses $\lvert444-440\rvert=4$ times per second. Adjusting until the beats slow to zero brings the tones into unison — the classical way to tune an instrument, and the principle of a *zero-beat* frequency comparison.
>
> **(b)** By sum-to-product, the beat signal is a product $2\cos(\Delta\omega t/2)\cos\bar\omega t$: a carrier at $\bar\omega$ multiplied by a sinusoid at $\Delta\omega/2$. That is AM with *suppressed* carrier (double-sideband suppressed carrier, DSB-SC): two sidebands at $\bar\omega\pm\Delta\omega/2$ and no component at $\bar\omega$ itself. The "envelope" crosses zero and the fast oscillation flips sign there, unlike true AM with $m<1$, whose envelope $1+m\cos\omega_mt$ stays positive.

### Total harmonic distortion

> [!definition] Total harmonic distortion
> For a periodic signal with fundamental RMS amplitude $X_1$ and harmonic RMS amplitudes $X_2,X_3,\dots$,
> $$
> \mathrm{THD}=\frac{\sqrt{X_2^2+X_3^2+X_4^2+\cdots}}{X_1}.
> $$
> It is often quoted in percent or in dB ($20\log_{10}\mathrm{THD}$). By orthogonality the total RMS value is $X_{\text{rms}}=X_1\sqrt{1+\mathrm{THD}^2}$ when there is no DC component.

(Some instruments normalise by the total RMS instead of the fundamental — "THD-R" or THD+N — so state which definition is used.)

> [!example] Worked Example 3.34 — THD of a rectifier load current
> A single-phase rectifier draws a 50 Hz current with RMS harmonic content $I_1=10$ A, $I_3=3$ A, $I_5=2$ A, $I_7=1$ A. Find the THD, the total RMS current, and the RMS current a meter responding only to the fundamental would read.
>
> **Solution.** $\mathrm{THD}=\sqrt{9+4+1}/10=\sqrt{14}/10=37.4\%$.
> $I_{\text{rms}}=\sqrt{100+9+4+1}=\sqrt{114}=10.68$ A $=10\sqrt{1+0.374^2}$ ✓.
> A fundamental-only meter reads 10 A, under-reading the heating current by 6.4 %.
>
> **Engineering consequence.** The harmonic currents carry no average power with a sinusoidal supply voltage (by orthogonality, $V_1I_n\cos(\cdot)$ averages to zero for $n\neq1$), yet they heat cables and transformers. The *true* power factor therefore drops below $\cos\Delta\phi_1$ by the factor $1/\sqrt{1+\mathrm{THD}^2}$.

### Harmonics in three-phase systems

Take a balanced set of identical, distorted phase currents, each a copy of the previous one delayed by one-third of a period: $i_b(t)=i_a(t-T_0/3)$, $i_c(t)=i_a(t-2T_0/3)$. A delay of $T_0/3$ shifts the $n$-th harmonic by $n\omega_0T_0/3=n\cdot120^\circ$. So the $n$-th harmonics of the three phases have phases $0^\circ$, $-120^\circ n$ and $-240^\circ n$, which reduce as follows:

| Harmonic $n$ | Phases (mod $360^\circ$) | Sequence | Sum of the three phases |
|---|---|---|---|
| $1, 4, 7, \dots$ ($n=3k+1$) | $0^\circ,-120^\circ,+120^\circ$ | positive ($abc$) | 0 |
| $2, 5, 8, \dots$ ($n=3k+2$) | $0^\circ,+120^\circ,-120^\circ$ | negative ($acb$) | 0 |
| $3, 6, 9, \dots$ ($n=3k$) | $0^\circ,0^\circ,0^\circ$ | zero (in phase) | $3\times$ one phase |

> [!theorem] Triplen harmonics
> In a balanced three-phase set, the harmonics of order $3,6,9,\dots$ (**triplens**) are in phase in all three lines. They do not cancel in the neutral: the neutral current at each triplen frequency is three times the per-phase value. All other harmonics form balanced sets (positive or negative sequence) and sum to zero.

> [!example] Worked Example 3.35 — Neutral current from single-phase loads
> A four-wire 50 Hz system feeds three identical banks of switch-mode supplies, one per phase. Each phase current has $I_1=20$ A, $I_3=12$ A, $I_5=6$ A (RMS). Find the RMS neutral current.
>
> **Solution.** Fundamentals and 5th harmonics form balanced sets and cancel in the neutral. The 3rd harmonics add in phase: $I_{N}=3\times12=36$ A at 150 Hz.
>
> **Comparison.** Each phase carries $\sqrt{400+144+36}=24.1$ A RMS; the neutral carries 36 A — **more than any phase conductor**, although the load is "balanced". This is why neutral conductors in buildings full of computers are oversized, and why transformers feeding such loads use a delta winding, in which triplen currents circulate instead of reaching the supply.

> [!tip] Engineering Note — Thinking in harmonics
> A good habit when meeting any nonlinear or time-varying operation is to ask "what does it do to a single sinusoid?" and then "what does it do to two?". Squaring doubles frequencies and creates DC; cubing triples them and creates $2f_1-f_2$; multiplying shifts; rectifying (§3.8) creates all even harmonics. Once you know the frequencies a device produces, filter design, channel planning and EMC compliance become bookkeeping.

## 3.8 Periodic function representation

### Motivation

Real periodic signals are rarely pure sinusoids: a clock is a square wave, a switching converter's inductor current is a triangle, a scanning circuit produces a sawtooth, a rectifier delivers humps of sine. We need to (i) write such waveforms down unambiguously, (ii) compute their mean and RMS values, (iii) recognise their symmetries, and (iv) express them as sums of harmonics so that the sinusoidal machinery of the previous sections — phasors, filtering, power by orthogonality — can be applied to them. Items (i)–(iii) are completed here; item (iv) is introduced here and developed fully as Fourier series in [Ch. 15](../part-5-transforms/15-fourier-series.md).

### Defining a periodic function from one period

A periodic function is completely determined by its values on any interval of length $T$ (one period). The standard way to specify it is to give a formula on, say, $0\le t<T$ (or $-T/2\le t<T/2$) and then declare the **periodic extension**

$$
x(t+kT)=x(t),\qquad k\in\mathbb{Z}.
$$

To evaluate the extension at an arbitrary $t$, reduce $t$ into the base interval with the floor function: $t\mapsto t-T\lfloor t/T\rfloor\in[0,T)$, which is what `mod(t, T)` computes for $T>0$.

The common engineering waveforms, with amplitude $A$, period $T$ and $\omega_0=2\pi/T$:

| Waveform | One period | Closed form for all $t$ |
|---|---|---|
| square wave ($\pm A$) | $A$ on $(0,T/2)$, $-A$ on $(T/2,T)$ | $A\operatorname{sgn}(\sin\omega_0t)$ |
| sawtooth ($-A$ to $A$) | $2At/T$ on $(-T/2,T/2)$ | $2A\big(\tfrac tT-\lfloor\tfrac tT+\tfrac12\rfloor\big)$ |
| triangle ($\pm A$) | linear between $-A$ and $A$ | $\tfrac{2A}{\pi}\arcsin(\sin\omega_0t)$ |
| full-wave rectified sine | $A\lvert\sin\omega_0t\rvert$ | $A\lvert\sin\omega_0t\rvert$ (period $T/2$) |
| half-wave rectified sine | $A\sin\omega_0t$ on $(0,T/2)$, $0$ on $(T/2,T)$ | $\tfrac A2\big(\sin\omega_0t+\lvert\sin\omega_0t\rvert\big)$ |
| pulse train, duty cycle $D$ | $A$ on $(0,DT)$, $0$ on $(DT,T)$ | $A\cdot\mathbf{1}\big[\operatorname{mod}(t,T)<DT\big]$ |

The triangle-wave formula works because $\arcsin(\sin\theta)$ equals $\theta$ on $[-\pi/2,\pi/2]$ and folds back linearly outside it — an application of the principal branch of §3.1.

```desmos-graph
left=-7; right=13; bottom=-1.5; top=1.5
---
y=\operatorname{sign}(\sin(x))|blue
y=\frac{2}{\pi}\arcsin(\sin(x))|red
y=2(\frac{x}{2\pi}-\operatorname{floor}(\frac{x}{2\pi}+\frac{1}{2}))|green
```
*Three standard waveforms of period $2\pi$ and amplitude 1: square (blue), triangle (red) and sawtooth (green), all generated from closed-form expressions.*

```desmos-graph
left=-1; right=13; bottom=-0.3; top=1.3
---
y=\left|\sin(x)\right||blue
y=\max(\sin(x),0)|red|dashed
y=\frac{2}{\pi}|blue|dotted
y=\frac{1}{\pi}|red|dotted
```
*Full-wave (blue) and half-wave (red, dashed) rectified sines with their mean values $2/\pi$ and $1/\pi$ (dotted). The full-wave output repeats every $\pi$ — twice the input frequency.*

> [!example] Worked Example 3.36 — Evaluating a periodic extension
> A sawtooth generator produces a waveform of amplitude 5 V and period 2 ms that rises linearly from $-5$ V to $+5$ V. Find its value at $t=3.7$ ms and at $t=-0.3$ ms.
>
> **Solution.** $x(t)=10\big(t/T-\lfloor t/T+\tfrac12\rfloor\big)$ V with $T=2$ ms.
> - $t=3.7$ ms: $t/T=1.85$, $\lfloor2.35\rfloor=2$, so $x=10(1.85-2)=-1.5$ V.
> - $t=-0.3$ ms: $t/T=-0.15$, $\lfloor0.35\rfloor=0$, so $x=10(-0.15)=-1.5$ V.
>
> **Check:** the two times differ by $4$ ms $=2T$, so the values must agree. ✓ Note that $\lfloor\cdot\rfloor$ rounds *down* for negative arguments ($\lfloor-0.65\rfloor=-1$), which is why the floor form handles negative $t$ correctly where truncation would not.

### Mean and RMS of standard waveforms

Using $\bar x=\frac1T\int_0^Tx\,dt$ and $X_{\text{rms}}^2=\frac1T\int_0^Tx^2\,dt$ (Ch. 5):

| Waveform (amplitude $A$) | Mean | RMS | Crest factor |
|---|---|---|---|
| sine | $0$ | $A/\sqrt2$ | $\sqrt2=1.414$ |
| square ($\pm A$) | $0$ | $A$ | $1$ |
| triangle or sawtooth ($\pm A$) | $0$ | $A/\sqrt3$ | $\sqrt3=1.732$ |
| full-wave rectified sine | $2A/\pi$ | $A/\sqrt2$ | $\sqrt2$ |
| half-wave rectified sine | $A/\pi$ | $A/2$ | $2$ |
| pulse train, $0$ to $A$, duty $D$ | $DA$ | $A\sqrt D$ | $1/\sqrt D$ |

**Derivations.** *Full-wave mean:* over one hump, $\frac1\pi\int_0^\pi A\sin\theta\,d\theta=\frac A\pi[-\cos\theta]_0^\pi=\frac{2A}{\pi}$. *Full-wave RMS:* $\lvert\sin\rvert^2=\sin^2$, so the RMS is that of the sine. *Half-wave:* the waveform equals the sine on half the period and zero on the other half, so its mean is half the full-wave mean, $A/\pi$, and its mean square is half of $A^2/2$, giving RMS $A/2$. *Sawtooth:* on $(-T/2,T/2)$, $x=2At/T$, so

$$
X_{\text{rms}}^2=\frac1T\int_{-T/2}^{T/2}\frac{4A^2t^2}{T^2}\,dt=\frac{4A^2}{T^3}\cdot\frac{2}{3}\Big(\frac T2\Big)^3=\frac{A^2}{3}.
$$

The triangle wave consists of the same linear ramps rearranged in time, so it has the same RMS. *Pulse train:* $x^2=A^2$ for a fraction $D$ of the time, so $X_{\text{rms}}^2=DA^2$.

> [!example] Worked Example 3.37 — Rectifier and PWM levels
> (a) A half-wave rectifier fed from 230 V RMS mains (peak 325 V) drives a resistive load. What DC voltmeter and true-RMS voltmeter readings are expected across the load?
> (b) A PWM output switches between 0 and 5 V at 25 % duty cycle. Find its mean and RMS, and the power delivered to a 10 $\Omega$ heater.
>
> **(a)** Mean $=325/\pi=103.5$ V. RMS $=325/2=162.5$ V.
>
> **(b)** Mean $=0.25\times5=1.25$ V. RMS $=5\sqrt{0.25}=2.5$ V. Power $=2.5^2/10=0.625$ W — which is $D\times5^2/10$, as it must be, and *not* $1.25^2/10=0.156$ W. Heating depends on the RMS value, not the mean.

### Symmetry

Symmetries of a waveform constrain which harmonics it can contain, often halving the work of analysis and giving an instant sanity check.

> [!definition] Even, odd and half-wave symmetry
> A periodic function $x$ with period $T$ is
> - **even** if $x(-t)=x(t)$ (mirror symmetry about $t=0$);
> - **odd** if $x(-t)=-x(t)$ (point symmetry about the origin);
> - **half-wave symmetric** if $x(t+T/2)=-x(t)$ (the second half-cycle is the first one inverted).

Every function splits uniquely into even and odd parts,

$$
x(t)=x_e(t)+x_o(t),\qquad x_e(t)=\frac{x(t)+x(-t)}2,\qquad x_o(t)=\frac{x(t)-x(-t)}2 .
$$

> [!theorem] Symmetry and harmonic content
> For a periodic function written as a sum of harmonics $a_0+\sum_{n\ge1}(a_n\cos n\omega_0t+b_n\sin n\omega_0t)$:
> 1. If $x$ is **even**, all $b_n=0$: only cosines (and possibly DC) appear.
> 2. If $x$ is **odd**, $a_0=0$ and all $a_n=0$: only sines appear.
> 3. If $x$ has **half-wave symmetry**, $a_0=0$ and all **even** harmonics vanish: only odd harmonics appear.
>
> Shifting the time origin changes even/odd symmetry but not half-wave symmetry, so property 3 is the more robust one.

**Proof of 3.** If $x(t+T/2)=-x(t)$, then the mean over a period is minus itself, so $a_0=0$. For the $n$-th harmonic coefficient (formula below), split $\int_0^T$ into $\int_0^{T/2}+\int_{T/2}^T$ and substitute $t=\tau+T/2$ in the second:

$$
\int_{T/2}^{T}x(t)\cos n\omega_0t\,dt=\int_0^{T/2}x(\tau+\tfrac T2)\cos(n\omega_0\tau+n\pi)\,d\tau=(-1)^{n+1}\int_0^{T/2}x(\tau)\cos n\omega_0\tau\,d\tau .
$$

So $a_n\propto\big(1+(-1)^{n+1}\big)\int_0^{T/2}x\cos n\omega_0t\,dt$, which is zero for even $n$; the same holds for $b_n$. Properties 1 and 2 follow because $\cos$ is even, $\sin$ is odd, and the integral over a symmetric interval of an odd function is zero. $\square$

> [!example] Worked Example 3.38 — Predicting harmonic content from symmetry
> Without computing any integrals, state which harmonics can appear in: (a) the square wave $A\operatorname{sgn}(\sin\omega_0t)$; (b) the triangle wave; (c) the sawtooth; (d) the full-wave rectified sine $A\lvert\sin\omega_0t\rvert$ (period $T$ of the *original* sine); (e) the half-wave rectified sine.
>
> - **(a)** Odd and half-wave symmetric: only sines of odd harmonics.
> - **(b)** $\frac{2A}{\pi}\arcsin(\sin\omega_0t)$ is odd and half-wave symmetric: only sines of odd harmonics.
> - **(c)** Odd but *not* half-wave symmetric: sines of all harmonics.
> - **(d)** Even, and it repeats every $T/2$, so in terms of $\omega_0$ it contains a DC term and only **even** harmonics $2\omega_0,4\omega_0,\dots$, all cosines. (Its own fundamental is $2\omega_0$.)
> - **(e)** Neither even nor odd and not half-wave symmetric. But $\tfrac A2(\sin\omega_0t+\lvert\sin\omega_0t\rvert)$ shows it is half a sine plus half of (d): DC, the fundamental $\sin\omega_0t$, and even cosine harmonics — and **no odd harmonics other than the fundamental**.
>
> These predictions are confirmed by the expansions in the summary table below and by numerical projection.

### Representing a periodic function by harmonics

The orthogonality relations of §3.7 suggest how a periodic function might be built from harmonics. *Suppose* that

$$
x(t)=a_0+\sum_{n=1}^{\infty}\big(a_n\cos n\omega_0t+b_n\sin n\omega_0t\big).
$$

Multiply both sides by $\cos m\omega_0t$ and integrate over a period. By orthogonality every term on the right vanishes except the $n=m$ cosine term, which contributes $a_mT/2$. Doing the same with $\sin m\omega_0t$ and with 1 gives:

> [!theorem] Harmonic (Fourier) coefficients
> If $x$ has a harmonic representation as above, its coefficients are
> $$
> a_0=\frac1T\int_0^Tx(t)\,dt,\qquad a_n=\frac2T\int_0^Tx(t)\cos n\omega_0t\,dt,\qquad b_n=\frac2T\int_0^Tx(t)\sin n\omega_0t\,dt .
> $$
> In amplitude–phase form (§3.2) the $n$-th harmonic is $A_n\cos(n\omega_0t+\phi_n)$ with $A_n=\sqrt{a_n^2+b_n^2}$ and $\phi_n=\operatorname{atan2}(-b_n,a_n)$. In complex form, $x(t)=\sum_{n=-\infty}^{\infty}c_ne^{jn\omega_0t}$ with $c_n=\frac1T\int_0^Tx(t)e^{-jn\omega_0t}\,dt$, $c_0=a_0$ and $c_{\pm n}=\tfrac12(a_n\mp jb_n)$ for $n\ge1$.

Each coefficient is a *correlation*: $a_n$ measures how much $x$ "looks like" $\cos n\omega_0t$. The questions of *which* functions have such representations, in what sense the infinite series converges, and what happens at discontinuities (the Gibbs phenomenon) are the subject of [Ch. 15](../part-5-transforms/15-fourier-series.md). For all the waveforms in this section the series converges to $x(t)$ at every point of continuity, and to the midpoint of the jump at a discontinuity.

> [!example] Worked Example 3.39 — Harmonics of a square wave
> Find the harmonic representation of the square wave $x=A$ on $(0,T/2)$, $x=-A$ on $(T/2,T)$, its THD, and check Parseval's relation.
>
> **Symmetry first.** Odd and half-wave symmetric, so $a_0=a_n=0$ and only odd-$n$ sines appear.
>
> **Coefficients.** Directly from the definition:
> $$
> b_n=\frac2T\Big[\int_0^{T/2}A\sin n\omega_0t\,dt-\int_{T/2}^TA\sin n\omega_0t\,dt\Big]=\frac{2A}{T}\cdot\frac{1}{n\omega_0}\Big[(1-\cos n\pi)-(\cos n\pi-\cos2n\pi)\Big]=\frac{2A}{n\pi}\big(1-(-1)^n\big).
> $$
> So $b_n=\dfrac{4A}{n\pi}$ for odd $n$ and $0$ for even $n$:
> $$
> x(t)=\frac{4A}{\pi}\Big(\sin\omega_0t+\frac13\sin3\omega_0t+\frac15\sin5\omega_0t+\cdots\Big).
> $$
>
> **THD.** The fundamental has RMS $\frac{4A}{\pi\sqrt2}=0.9003A$; the total RMS is $A$. So $\mathrm{THD}=\sqrt{A^2-0.8106A^2}/(0.9003A)=\sqrt{\pi^2/8-1}=48.3\%$.
>
> **Parseval check.** $\sum_{n\text{ odd}}\frac12\big(\frac{4A}{n\pi}\big)^2=\frac{8A^2}{\pi^2}\sum_{n\text{ odd}}\frac1{n^2}=\frac{8A^2}{\pi^2}\cdot\frac{\pi^2}{8}=A^2=X_{\text{rms}}^2$ ✓, using the classical sum $1+\frac19+\frac1{25}+\cdots=\pi^2/8$ (itself provable from this very example in Ch. 15). Numerical projection gives $b_1=1.2732A=4A/\pi$, $b_3=0.4244A$, $b_5=0.2546A$. ✓

```desmos-graph
left=-0.5; right=6.8; bottom=-1.5; top=1.5
---
y=\operatorname{sign}(\sin(x))|grey|dashed
y=\frac{4}{\pi}\sin(x)|red|dotted
y=\sum_{n=1}^{3}\frac{4}{\pi(2n-1)}\sin((2n-1)x)|green
y=\sum_{n=1}^{10}\frac{4}{\pi(2n-1)}\sin((2n-1)x)|blue
```
*Building a square wave (grey) from harmonics: the fundamental alone (red), the sum through the 5th harmonic (green), and through the 19th (blue). The ripple narrows toward the jumps but its overshoot stays near 9 % of the jump — the Gibbs phenomenon of Ch. 15.*

> [!example] Worked Example 3.40 — Harmonics of a full-wave rectified sine
> Find the DC and the first two nonzero harmonics of $x(t)=A\lvert\sin\omega_0t\rvert$, and interpret for a full-wave rectifier fed at 50 Hz.
>
> **Symmetry.** Even with period $T/2$: only DC and cosines at $2\omega_0,4\omega_0,\dots$ (Example 3.38d).
>
> **DC.** $a_0=2A/\pi$ (computed above).
>
> **Harmonics.** By evenness and the $T/2$ periodicity, $a_n=\frac{2}{T}\int_0^T A\lvert\sin\omega_0t\rvert\cos n\omega_0t\,dt=\frac{2A}{\pi}\int_0^{\pi}\sin\theta\cos n\theta\,d\theta$ for even $n$ (with $\theta=\omega_0t$, using that the two half-periods contribute equally). Product-to-sum: $\sin\theta\cos n\theta=\tfrac12[\sin(n+1)\theta-\sin(n-1)\theta]$, and $\int_0^\pi\sin k\theta\,d\theta=\frac{1-\cos k\pi}{k}=\frac2k$ for odd $k$. With $n$ even, $n\pm1$ are odd:
> $$
> a_n=\frac{2A}{\pi}\cdot\frac12\Big[\frac{2}{n+1}-\frac{2}{n-1}\Big]=\frac{2A}{\pi}\cdot\frac{-2}{n^2-1}=-\frac{4A}{\pi(n^2-1)}.
> $$
> So $a_2=-\frac{4A}{3\pi}=-0.4244A$ and $a_4=-\frac{4A}{15\pi}=-0.0849A$:
> $$
> A\lvert\sin\omega_0t\rvert=\frac{2A}{\pi}-\frac{4A}{3\pi}\cos2\omega_0t-\frac{4A}{15\pi}\cos4\omega_0t-\cdots
> $$
>
> **Interpretation.** From a 50 Hz supply the rectifier output is DC ($0.637A$) plus a dominant 100 Hz ripple of amplitude $0.424A$ and a 200 Hz component one-fifth as large. A smoothing filter is designed against the 100 Hz term.
>
> **Check:** numerical projection gives $a_2=-0.42441A$, $a_4=-0.08488A$, $a_6=-0.03638A$ $(=-4A/35\pi)$. ✓

### Effect of time shift, scaling and smoothness on the harmonics

Three general rules let you reuse known expansions:

- **Time shift.** If $x(t)$ has harmonics $A_n\cos(n\omega_0t+\phi_n)$, then $x(t-t_0)$ has harmonics $A_n\cos(n\omega_0t+\phi_n-n\omega_0t_0)$: amplitudes are unchanged; the $n$-th phase shifts by $-n\omega_0t_0$ — proportional to $n$, the "linear phase" of a pure delay seen in §3.3. Example: the square wave shifted by $T/4$ becomes even, and its expansion becomes $\frac{4A}{\pi}\big(\cos\omega_0t-\frac13\cos3\omega_0t+\frac15\cos5\omega_0t-\cdots\big)$.
- **Amplitude scaling and offset.** $\alpha x(t)+\beta$ multiplies every harmonic by $\alpha$ and adds $\beta$ to the DC term. A 0-to-$A$ square wave is $\frac A2+\frac12(\pm A\text{ square wave})$.
- **Smoothness and decay.** The smoother the waveform, the faster the harmonic amplitudes decay with $n$. A waveform with jumps (square, sawtooth) has $A_n\sim1/n$; a continuous waveform with corners (triangle, rectified sine) has $A_n\sim1/n^2$; each further continuous derivative buys another factor of $1/n$. Precise statements are in Ch. 15.

| Waveform | Harmonic expansion | Decay |
|---|---|---|
| square $\pm A$ | $\dfrac{4A}{\pi}\displaystyle\sum_{n\text{ odd}}\frac{\sin n\omega_0t}{n}$ | $1/n$ |
| sawtooth $-A$ to $A$ | $\dfrac{2A}{\pi}\displaystyle\sum_{n\ge1}\frac{(-1)^{n+1}\sin n\omega_0t}{n}$ | $1/n$ |
| triangle $\pm A$ (odd, peak at $T/4$) | $\dfrac{8A}{\pi^2}\displaystyle\sum_{n\text{ odd}}\frac{(-1)^{(n-1)/2}\sin n\omega_0t}{n^2}$ | $1/n^2$ |
| full-wave $A\lvert\sin\omega_0t\rvert$ | $\dfrac{2A}{\pi}-\dfrac{4A}{\pi}\displaystyle\sum_{k\ge1}\frac{\cos2k\omega_0t}{4k^2-1}$ | $1/n^2$ |
| half-wave rectified | $\dfrac{A}{\pi}+\dfrac A2\sin\omega_0t-\dfrac{2A}{\pi}\displaystyle\sum_{k\ge1}\frac{\cos2k\omega_0t}{4k^2-1}$ | $1/n^2$ |

> [!tip] Engineering Note — Rise time, bandwidth and EMI
> The $1/n$ decay of a square wave means that a digital clock edge contains significant energy at many harmonics of the clock frequency; a 100 MHz clock with fast edges produces radiated emissions at 300, 500, 700 MHz and beyond. Slowing the edges (rounding the corners) converts the spectrum to $1/n^2$ decay above a frequency of roughly $1/(\pi t_r)$, where $t_r$ is the rise time — a standard EMC design rule. Conversely, an oscilloscope or channel that passes only the first few harmonics turns a square wave into the rippled shape in the partial-sum graph above.

> [!warning] Common Pitfall — Which fundamental?
> When a waveform has a hidden shorter period (such as $\lvert\sin\omega_0t\rvert$, whose period is $T/2$), state whether harmonic numbers refer to the original $\omega_0$ or to the waveform's own fundamental $2\omega_0$. "The 2nd harmonic of the rectifier output" is ambiguous: it is 200 Hz if counted from the 100 Hz ripple, and 100 Hz if counted from the 50 Hz supply. Power-system practice counts from the supply frequency.

## Summary

| Topic | Key formulas |
|---|---|
| Angle units | $\theta_{\text{rad}}=s/r$; $\pi\text{ rad}=180^\circ$; $\Omega=2\pi N/60$ (rpm); $\theta_e=p\,\theta_m$ |
| Arc, sector | $s=r\theta$; $A=\tfrac12r^2\theta$; $v=r\Omega$ |
| Fundamental limit | $\lim_{\theta\to0}\sin\theta/\theta=1$ (radians); $\sin\theta\approx\theta$, $\cos\theta\approx1-\theta^2/2$ |
| Quadrant-safe angle | $\operatorname{atan2}(y,x)\in(-\pi,\pi]$ |
| General sinusoid | $x=A\cos(\omega t+\phi)+C$; shift $=-\phi/\omega$ |
| Sine/cosine conversion | $\sin\theta=\cos(\theta-90^\circ)$; $-\cos\theta=\cos(\theta\pm180^\circ)$ |
| Amplitude–phase form | $a\cos\omega t+b\sin\omega t=R\cos(\omega t-\alpha)$, $R=\sqrt{a^2+b^2}$, $\alpha=\operatorname{atan2}(b,a)$ |
| Derivative / integral | $\frac{d}{dt}$: $\times\omega$, $+90^\circ$; $\int dt$: $\div\omega$, $-90^\circ$ |
| Mean, RMS of sinusoid | $\bar x=C$; $X_{\text{rms}}=\sqrt{C^2+A^2/2}$; pure sine $A/\sqrt2$ |
| Phase difference | $\Delta\phi=\phi_1-\phi_2\in(-180^\circ,180^\circ]$; $\Delta t=\Delta\phi/\omega$ |
| Same-frequency sum | $R^2=A_1^2+A_2^2+2A_1A_2\cos(\phi_1-\phi_2)$ |
| Three-phase | $\sum_{k=0}^{2}\cos(\omega t-120^\circ k)=0$; $V_{LL}=\sqrt3V_{ph}\angle30^\circ$ |
| Lissajous | $\lvert\sin\phi\rvert=y_0/Y$ |
| Period, frequency | $T=2\pi/\omega=1/f$; $\omega=2\pi f$ |
| Sum periodicity | periodic iff $f_1/f_2$ rational; $f_0=\gcd(f_1,f_2)$ |
| Waves | $\cos(\omega t-kz)$; $\lambda=2\pi/k=v_p/f$; $v_p=\omega/k$ |
| Discrete sinusoid | $\Omega=\omega T_s$; periodic iff $\Omega/2\pi=m/N$; $\Omega\equiv\Omega+2\pi$ |
| Instantaneous frequency | $\omega_i=d\Theta/dt$ |
| SHM | $\ddot x+\omega_0^2x=0\Rightarrow x=A\cos(\omega_0t+\phi)$; $A=\sqrt{x_0^2+(v_0/\omega_0)^2}$ |
| LC tank | $\omega_0=1/\sqrt{LC}$ |
| Damped oscillation | $Ae^{-\alpha t}\cos(\omega_dt+\phi)$, $\omega_d=\sqrt{\omega_0^2-\alpha^2}$, $Q=\omega_0/2\alpha$ |
| Compound angle | $\cos(A\pm B)=\cos A\cos B\mp\sin A\sin B$; $\sin(A\pm B)=\sin A\cos B\pm\cos A\sin B$ |
| Power reduction | $\cos^2A=\frac{1+\cos2A}2$; $\sin^2A=\frac{1-\cos2A}2$; $\cos^3A=\frac{3\cos A+\cos3A}4$ |
| Powers of cosine | $\cos^n\theta=2^{-n}\sum_k\binom nk\cos(n-2k)\theta$ |
| Chebyshev | $\cos n\theta=T_n(\cos\theta)$, $T_{n+1}=2xT_n-T_{n-1}$ |
| Product-to-sum | $\cos A\cos B=\frac12[\cos(A-B)+\cos(A+B)]$ |
| Sum-to-product | $\cos P+\cos Q=2\cos\frac{P+Q}2\cos\frac{P-Q}2$ |
| Phasor-array sum | $\big\lvert\sum_{k=0}^{N-1}e^{jk\delta}\big\rvert=\lvert\sin(N\delta/2)/\sin(\delta/2)\rvert$ |
| Average power | $P=V_{\text{rms}}I_{\text{rms}}\cos(\phi_v-\phi_i)$ |
| Orthogonality | $\frac{2}{T}\int_0^T\cos m\omega_0t\cos n\omega_0t\,dt=\delta_{mn}$ ($m,n\ge1$) |
| Power of harmonic sum | $X_{\text{rms}}^2=A_0^2+\sum A_n^2/2$ |
| Cubic nonlinearity | $x^3\to\frac{3A^3}{4}\cos\omega t+\frac{A^3}{4}\cos3\omega t$ |
| IM products | $\lvert p f_1+q f_2\rvert$, order $\lvert p\rvert+\lvert q\rvert$; IM3 at $2f_1-f_2$, $2f_2-f_1$ |
| Beats | $f_{\text{beat}}=\lvert f_1-f_2\rvert$ |
| AM | $A_c(1+m\cos\omega_mt)\cos\omega_ct$: sidebands $mA_c/2$ at $\omega_c\pm\omega_m$ |
| THD | $\sqrt{\sum_{n\ge2}X_n^2}/X_1$ |
| Triplens | harmonics $3k$ are zero-sequence: neutral current $=3\times$ per phase |
| Standard waveform RMS | square $A$; triangle/sawtooth $A/\sqrt3$; half-wave $A/2$; pulse $A\sqrt D$ |
| Symmetry | even → cosines; odd → sines; half-wave → odd harmonics only |
| Harmonic coefficients | $a_0=\bar x$; $a_n,b_n=\frac2T\int_0^Tx\,\{\cos,\sin\}\,n\omega_0t\,dt$ |
| Square wave | $\frac{4A}{\pi}\sum_{n\text{ odd}}\frac{\sin n\omega_0t}{n}$; THD $48.3\%$ |

## Further reading

- E. Kreyszig, *Advanced Engineering Mathematics*, 10th ed., Wiley, 2011 — Appendix on trigonometric identities and Ch. 11 on Fourier series, the natural continuation of §3.7–3.8.
- K. A. Stroud and D. J. Booth, *Engineering Mathematics*, 8th ed., Palgrave/Bloomsbury, 2020 — patient, example-driven treatment of trigonometry, sinusoids and phasors.
- A. V. Oppenheim and A. S. Willsky, *Signals and Systems*, 2nd ed., Prentice Hall, 1997 — Ch. 1 on continuous- and discrete-time sinusoids, periodicity and aliasing; Ch. 3 on harmonic representation.
- C. K. Alexander and M. N. O. Sadiku, *Fundamentals of Electric Circuits*, 7th ed., McGraw-Hill, 2021 — sinusoids, phase, average power, three-phase systems and harmonics from the circuit side.
