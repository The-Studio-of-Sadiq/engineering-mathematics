# Chapter 7 — Vector Calculus

[← Previous](../part-2-calculus/06-multivariable-calculus.md) · [Table of Contents](../00-front/table-of-contents.md) · [Next →](../part-3-linear-algebra/08-matrices-and-linear-systems.md)

> **Prerequisites:** [Ch. 4](../part-2-calculus/04-differential-calculus.md), [Ch. 5](../part-2-calculus/05-integral-calculus.md), [Ch. 6](../part-2-calculus/06-multivariable-calculus.md)
> **Learning outcomes:** After this chapter you can
> - manipulate vectors in Cartesian, cylindrical and spherical bases, including dot, cross and triple products;
> - differentiate vector-valued functions and describe curves, velocities and field lines;
> - compute gradient, divergence, curl and Laplacian in all three coordinate systems and interpret them physically;
> - set up and evaluate line, surface and volume integrals, choosing parametrisations and orientations correctly;
> - state, apply and sketch the proofs of the theorems of Green, Stokes and Gauss (divergence theorem);
> - convert Maxwell's equations between integral and differential form, and use them to solve symmetric field problems.

Chapter 6 gave us functions of several variables and a first look at the gradient. This chapter adds the missing ingredient for field theory: **direction**. An electric field, a current density, a heat flux and the velocity of a fluid each attach a *vector* to every point of space. Engineers need to ask three kinds of question about such fields:

1. **Local questions** — at a point, how fast does a potential change, how much does a field spread out, how much does it swirl? The answers are the gradient, divergence and curl.
2. **Global questions** — how much work does a field do along a wire, how much of it flows through a surface, how much charge sits in a volume? The answers are line, surface and volume integrals.
3. **The bridge** — Green's, Stokes' and the divergence theorem say that a global integral over a boundary equals the integral of a local derivative over the inside. They are the multivariable versions of the fundamental theorem of calculus.

Put the three together and Maxwell's equations fall out in two equivalent forms — the integral form you measure, and the differential form you solve. That is the payoff in Section 7.12. Throughout, the methods are general; electromagnetics is the running motivation.

**Conventions.** Cartesian unit vectors are $\hat{\mathbf{x}},\hat{\mathbf{y}},\hat{\mathbf{z}}$. Cylindrical coordinates are $(r,\phi,z)$ with unit vectors $\hat{\mathbf{r}},\hat{\boldsymbol{\phi}},\hat{\mathbf{z}}$; spherical coordinates are $(R,\theta,\phi)$ with $\theta$ the polar angle from $+z$ and unit vectors $\hat{\mathbf{R}},\hat{\boldsymbol{\theta}},\hat{\boldsymbol{\phi}}$ (as in [Ch. 6](../part-2-calculus/06-multivariable-calculus.md), §6.9). The position vector of a point is written $\mathbf{R}=x\hat{\mathbf{x}}+y\hat{\mathbf{y}}+z\hat{\mathbf{z}}$, with length $R=\lVert\mathbf{R}\rVert$, so that $\mathbf{R}=R\hat{\mathbf{R}}$ in spherical coordinates. A curve is written $\mathbf{r}(t)$ when it is parametrised by a variable $t$.

## 7.1 Scalars and vectors

### Why the distinction matters

A **scalar** is a quantity fully described by one number (with a unit): temperature, voltage, charge, energy, resistance. A **vector** needs a magnitude *and* a direction: force, velocity, electric field, current density, magnetic flux density. Mixing them up is not a matter of taste — adding a scalar to a vector is meaningless, and the direction of a field is often the whole point (it decides, for example, whether a force accelerates or decelerates an electron, or whether flux enters or leaves a surface).

> [!definition] Vector
> A **vector** in three-dimensional space is a quantity with magnitude and direction, which adds by the parallelogram rule. Relative to a right-handed orthonormal basis $\hat{\mathbf{x}},\hat{\mathbf{y}},\hat{\mathbf{z}}$ it is written
> $$
> \mathbf{A}=A_x\hat{\mathbf{x}}+A_y\hat{\mathbf{y}}+A_z\hat{\mathbf{z}}=(A_x,A_y,A_z),
> $$
> where the **components** $A_x,A_y,A_z$ are scalars. Its **magnitude** is $\lVert\mathbf{A}\rVert=A=\sqrt{A_x^2+A_y^2+A_z^2}$ and, if $A\neq0$, its **unit vector** is $\hat{\mathbf{a}}_A=\mathbf{A}/A$.

A vector is a geometric object; its components depend on the basis you choose. Rotating the axes changes the three numbers but not the arrow. This is why physical laws are written in vector form ($\mathbf{F}=q\mathbf{E}$) — they then hold in every coordinate system.

**Algebra.** Addition and scalar multiplication act component by component:

$$
\mathbf{A}+\mathbf{B}=(A_x+B_x,\ A_y+B_y,\ A_z+B_z),\qquad c\mathbf{A}=(cA_x,\ cA_y,\ cA_z).
$$

These obey the familiar laws: $\mathbf{A}+\mathbf{B}=\mathbf{B}+\mathbf{A}$, $(\mathbf{A}+\mathbf{B})+\mathbf{C}=\mathbf{A}+(\mathbf{B}+\mathbf{C})$, $c(\mathbf{A}+\mathbf{B})=c\mathbf{A}+c\mathbf{B}$. The zero vector $\mathbf{0}$ has no direction. Geometrically, $\mathbf{A}+\mathbf{B}$ is the diagonal of the parallelogram spanned by $\mathbf{A}$ and $\mathbf{B}$, or equivalently the result of placing the tail of $\mathbf{B}$ at the head of $\mathbf{A}$.

### Position and separation vectors

The **position vector** $\mathbf{R}=(x,y,z)$ points from the origin to the point $P(x,y,z)$. The vector from a *source point* $P'$ at $\mathbf{R}'$ to a *field point* $P$ at $\mathbf{R}$ is the **separation vector**

$$
\mathbf{R}-\mathbf{R}'=(x-x',\ y-y',\ z-z'),\qquad \lVert\mathbf{R}-\mathbf{R}'\rVert=\sqrt{(x-x')^2+(y-y')^2+(z-z')^2}.
$$

Every inverse-square law in field theory is written in terms of it. Coulomb's law for the force on charge $q$ at $\mathbf{R}$ due to charge $q'$ at $\mathbf{R}'$ is

$$
\mathbf{F}=\frac{qq'}{4\pi\varepsilon_0}\,\frac{\mathbf{R}-\mathbf{R}'}{\lVert\mathbf{R}-\mathbf{R}'\rVert^3}.
$$

The cube in the denominator is not a typo: one power of the distance normalises the numerator to a unit vector, and the other two give the inverse-square law.

> [!example] Worked Example — separation vector and unit vector
> A point charge sits at $P'(1,2,0)$ m and we observe at $P(3,-1,6)$ m. Find the separation vector, its length and its unit vector.
>
> **Step 1.** $\mathbf{R}-\mathbf{R}'=(3-1,\ -1-2,\ 6-0)=(2,-3,6)$ m.
>
> **Step 2.** $\lVert\mathbf{R}-\mathbf{R}'\rVert=\sqrt{4+9+36}=\sqrt{49}=7$ m.
>
> **Step 3.** $\hat{\mathbf{a}}=(2,-3,6)/7=(0.2857,\ -0.4286,\ 0.8571)$. Check: $0.2857^2+0.4286^2+0.8571^2=0.0816+0.1837+0.7347=1.000$.

> [!example] Worked Example — superposition of Coulomb forces
> Charges $q_1=+2\ \mathrm{nC}$ at $(0,0,0)$ and $q_2=-1\ \mathrm{nC}$ at $(2,0,0)$ m act on a test charge $q=+1\ \mathrm{nC}$ at $(0,1,0)$ m. Find the net force. Use $1/(4\pi\varepsilon_0)=8.988\times10^{9}\ \mathrm{N\,m^2/C^2}$.
>
> **Step 1 — force from $q_1$.** Separation $(0,1,0)$, length $1$. 
> $\mathbf{F}_1=8.988\times10^{9}\cdot(2\times10^{-9})(1\times10^{-9})\,(0,1,0)/1^3=(0,\ 17.98,\ 0)\ \mathrm{nN}$.
>
> **Step 2 — force from $q_2$.** Separation $(0-2,\ 1-0,\ 0)=(-2,1,0)$, length $\sqrt5$, cube $5\sqrt5=11.180$. 
> $\mathbf{F}_2=8.988\times10^{9}\cdot(-1\times10^{-9})(1\times10^{-9})\,(-2,1,0)/11.180=(1.608,\ -0.804,\ 0)\ \mathrm{nN}$.
>
> **Step 3 — superpose.** $\mathbf{F}=\mathbf{F}_1+\mathbf{F}_2=(1.608,\ 17.17,\ 0)\ \mathrm{nN}$, magnitude $17.25\ \mathrm{nN}$, at $\arctan(17.17/1.608)=84.6^\circ$ from $+x$.
>
> The negative $q_2$ *attracts* the test charge towards $(2,0,0)$, which is why $F_{2x}>0$ and $F_{2y}<0$. Forces add as vectors, never as magnitudes.

### Curvilinear bases

In cylindrical and spherical coordinates we also use orthonormal bases, but the basis vectors **change direction from point to point**. At a point with azimuth $\phi$,

$$
\hat{\mathbf{r}}=\cos\phi\,\hat{\mathbf{x}}+\sin\phi\,\hat{\mathbf{y}},\qquad
\hat{\boldsymbol{\phi}}=-\sin\phi\,\hat{\mathbf{x}}+\cos\phi\,\hat{\mathbf{y}},\qquad
\hat{\mathbf{z}}=\hat{\mathbf{z}},
$$

and in spherical coordinates

$$
\begin{aligned}
\hat{\mathbf{R}}&=\sin\theta\cos\phi\,\hat{\mathbf{x}}+\sin\theta\sin\phi\,\hat{\mathbf{y}}+\cos\theta\,\hat{\mathbf{z}},\\
\hat{\boldsymbol{\theta}}&=\cos\theta\cos\phi\,\hat{\mathbf{x}}+\cos\theta\sin\phi\,\hat{\mathbf{y}}-\sin\theta\,\hat{\mathbf{z}},\\
\hat{\boldsymbol{\phi}}&=-\sin\phi\,\hat{\mathbf{x}}+\cos\phi\,\hat{\mathbf{y}}.
\end{aligned}
$$

Each unit vector points in the direction in which its own coordinate increases. Both triads are right-handed: $\hat{\mathbf{r}}\times\hat{\boldsymbol{\phi}}=\hat{\mathbf{z}}$ and $\hat{\mathbf{R}}\times\hat{\boldsymbol{\theta}}=\hat{\boldsymbol{\phi}}$.

Because each triad is orthonormal, the component of $\mathbf{A}$ along a basis vector is a dot product (Section 7.2), and converting components is a rotation:

$$
\begin{bmatrix}A_r\\A_\phi\\A_z\end{bmatrix}=
\begin{bmatrix}\cos\phi&\sin\phi&0\\-\sin\phi&\cos\phi&0\\0&0&1\end{bmatrix}
\begin{bmatrix}A_x\\A_y\\A_z\end{bmatrix},
$$

$$
\begin{bmatrix}A_R\\A_\theta\\A_\phi\end{bmatrix}=
\begin{bmatrix}\sin\theta\cos\phi&\sin\theta\sin\phi&\cos\theta\\\cos\theta\cos\phi&\cos\theta\sin\phi&-\sin\theta\\-\sin\phi&\cos\phi&0\end{bmatrix}
\begin{bmatrix}A_x\\A_y\\A_z\end{bmatrix}.
$$

The matrices are orthogonal (rows are orthonormal), so the inverse conversion uses the transpose ([Ch. 9](../part-3-linear-algebra/09-vector-spaces-and-linear-transformations.md)).

> [!example] Worked Example — Cartesian to cylindrical components
> Express $\mathbf{A}=y\,\hat{\mathbf{x}}-x\,\hat{\mathbf{y}}+z\,\hat{\mathbf{z}}$ in cylindrical coordinates, and evaluate it at $(x,y,z)=(1,\sqrt3,2)$.
>
> **Step 1.** Substitute $x=r\cos\phi$, $y=r\sin\phi$: $A_x=r\sin\phi$, $A_y=-r\cos\phi$, $A_z=z$.
>
> **Step 2.** Rotate:
> $A_r=\cos\phi\,(r\sin\phi)+\sin\phi\,(-r\cos\phi)=0$,
> $A_\phi=-\sin\phi\,(r\sin\phi)+\cos\phi\,(-r\cos\phi)=-r$,
> $A_z=z$.
> So $\mathbf{A}=-r\,\hat{\boldsymbol{\phi}}+z\,\hat{\mathbf{z}}$ — a clockwise swirl (seen from $+z$) plus a vertical part. The cylindrical form makes the geometry obvious, which the Cartesian form hid.
>
> **Step 3.** At $(1,\sqrt3,2)$: $r=2$, $\phi=60^\circ$, so $\mathbf{A}=-2\hat{\boldsymbol{\phi}}+2\hat{\mathbf{z}}$. Check in Cartesian: $\mathbf{A}=(\sqrt3,-1,2)$, $\lVert\mathbf{A}\rVert=\sqrt{3+1+4}=\sqrt8$; cylindrical: $\sqrt{4+4}=\sqrt8$. ✓

> [!warning] Common Pitfall
> You may **not** add the cylindrical or spherical components of vectors located at *different points*. $\hat{\mathbf{r}}$ at $\phi=0$ is $\hat{\mathbf{x}}$, but at $\phi=\pi$ it is $-\hat{\mathbf{x}}$. Two vectors "$1\,\hat{\mathbf{r}}$" at those two points sum to zero, not to $2\hat{\mathbf{r}}$. Convert to Cartesian components before superposing contributions from many source points — this is exactly why field integrals such as the field of a ring of charge are done in Cartesian components.

> [!tip] Engineering Note
> Field solvers and antenna-pattern software report fields in whichever basis matches the symmetry: a coaxial cable's $\mathbf{E}$ is purely $\hat{\mathbf{r}}$ and its $\mathbf{H}$ purely $\hat{\boldsymbol{\phi}}$; an antenna's far field is purely $\hat{\boldsymbol{\theta}}$ and $\hat{\boldsymbol{\phi}}$. Choosing the basis that matches the symmetry reduces a three-component problem to a one-component one.

## 7.2 Vector operations

There are two ways to multiply two vectors, and each answers a different geometric question. The dot product measures *how much two vectors agree*; the cross product measures *how much they span*.

### The dot (scalar) product

> [!definition] Dot product
> $$
> \mathbf{A}\cdot\mathbf{B}=AB\cos\theta_{AB}=A_xB_x+A_yB_y+A_zB_z,
> $$
> where $\theta_{AB}\in[0,\pi]$ is the angle between the vectors.

**Why the two forms agree.** Apply the law of cosines to the triangle with sides $\mathbf{A}$, $\mathbf{B}$ and $\mathbf{A}-\mathbf{B}$: $\lVert\mathbf{A}-\mathbf{B}\rVert^2=A^2+B^2-2AB\cos\theta_{AB}$. Expanding the left side in components gives $A^2+B^2-2(A_xB_x+A_yB_y+A_zB_z)$. Comparing the two proves the identity.

Properties: commutative, distributive, $\mathbf{A}\cdot\mathbf{A}=A^2$, and $\mathbf{A}\cdot\mathbf{B}=0$ iff the vectors are perpendicular (or one is zero). For orthonormal bases $\hat{\mathbf{x}}\cdot\hat{\mathbf{x}}=1$, $\hat{\mathbf{x}}\cdot\hat{\mathbf{y}}=0$, and so on.

**Projection.** The scalar component of $\mathbf{A}$ along a unit vector $\hat{\mathbf{u}}$ is $\mathbf{A}\cdot\hat{\mathbf{u}}$, and the vector projection is $(\mathbf{A}\cdot\hat{\mathbf{u}})\hat{\mathbf{u}}$. The rest, $\mathbf{A}-(\mathbf{A}\cdot\hat{\mathbf{u}})\hat{\mathbf{u}}$, is perpendicular to $\hat{\mathbf{u}}$. This decomposition into normal and tangential parts is used constantly at material boundaries.

Physically, the dot product appears whenever only the *aligned* part of a vector matters: work $dW=\mathbf{F}\cdot d\mathbf{l}$, flux $d\Phi=\mathbf{B}\cdot d\mathbf{S}$, power $p=\mathbf{F}\cdot\mathbf{v}$.

### The cross (vector) product

> [!definition] Cross product
> $$
> \mathbf{A}\times\mathbf{B}=AB\sin\theta_{AB}\,\hat{\mathbf{n}}
> =\begin{vmatrix}\hat{\mathbf{x}}&\hat{\mathbf{y}}&\hat{\mathbf{z}}\\A_x&A_y&A_z\\B_x&B_y&B_z\end{vmatrix}
> =(A_yB_z-A_zB_y,\ A_zB_x-A_xB_z,\ A_xB_y-A_yB_x),
> $$
> where $\hat{\mathbf{n}}$ is the unit normal to the plane of $\mathbf{A}$ and $\mathbf{B}$ chosen by the **right-hand rule** (fingers from $\mathbf{A}$ to $\mathbf{B}$, thumb along $\hat{\mathbf{n}}$).

Key facts:

- $\lVert\mathbf{A}\times\mathbf{B}\rVert=AB\sin\theta_{AB}$ is the **area of the parallelogram** spanned by $\mathbf{A}$ and $\mathbf{B}$.
- $\mathbf{A}\times\mathbf{B}$ is perpendicular to both $\mathbf{A}$ and $\mathbf{B}$.
- **Anticommutative:** $\mathbf{B}\times\mathbf{A}=-\mathbf{A}\times\mathbf{B}$, hence $\mathbf{A}\times\mathbf{A}=\mathbf{0}$; parallel vectors have zero cross product.
- Distributive over addition, but **not associative**: $(\mathbf{A}\times\mathbf{B})\times\mathbf{C}\ne\mathbf{A}\times(\mathbf{B}\times\mathbf{C})$ in general.
- Cyclic basis rule: $\hat{\mathbf{x}}\times\hat{\mathbf{y}}=\hat{\mathbf{z}}$, $\hat{\mathbf{y}}\times\hat{\mathbf{z}}=\hat{\mathbf{x}}$, $\hat{\mathbf{z}}\times\hat{\mathbf{x}}=\hat{\mathbf{y}}$; the same cyclic rule holds for $(\hat{\mathbf{r}},\hat{\boldsymbol{\phi}},\hat{\mathbf{z}})$ and $(\hat{\mathbf{R}},\hat{\boldsymbol{\theta}},\hat{\boldsymbol{\phi}})$ at a given point.

Physically, the cross product appears whenever a "turning" or "sideways" effect occurs: torque $\boldsymbol{\tau}=\mathbf{r}\times\mathbf{F}$, the magnetic force $\mathbf{F}=q\mathbf{v}\times\mathbf{B}$, the Biot–Savart element $d\mathbf{B}\propto I\,d\mathbf{l}\times\hat{\mathbf{R}}/R^2$, the Poynting vector $\mathbf{S}=\mathbf{E}\times\mathbf{H}$, and the normal to a surface $\mathbf{r}_u\times\mathbf{r}_v$ (Section 7.7).

### Triple products

> [!theorem] Scalar and vector triple products
> **Scalar triple product.**
> $$
> \mathbf{A}\cdot(\mathbf{B}\times\mathbf{C})=\begin{vmatrix}A_x&A_y&A_z\\B_x&B_y&B_z\\C_x&C_y&C_z\end{vmatrix}
> =\mathbf{B}\cdot(\mathbf{C}\times\mathbf{A})=\mathbf{C}\cdot(\mathbf{A}\times\mathbf{B}).
> $$
> Its absolute value is the volume of the parallelepiped with edges $\mathbf{A},\mathbf{B},\mathbf{C}$; it vanishes iff the three vectors are coplanar.
>
> **Vector triple product ("BAC–CAB").**
> $$
> \mathbf{A}\times(\mathbf{B}\times\mathbf{C})=\mathbf{B}(\mathbf{A}\cdot\mathbf{C})-\mathbf{C}(\mathbf{A}\cdot\mathbf{B}).
> $$

*Proof sketch of BAC–CAB.* Both sides are linear in each of $\mathbf{A},\mathbf{B},\mathbf{C}$, so it suffices to check basis vectors. The left side is perpendicular to $\mathbf{B}\times\mathbf{C}$, hence lies in the plane of $\mathbf{B}$ and $\mathbf{C}$; e.g. $\hat{\mathbf{x}}\times(\hat{\mathbf{x}}\times\hat{\mathbf{y}})=\hat{\mathbf{x}}\times\hat{\mathbf{z}}=-\hat{\mathbf{y}}$, while the right side gives $\hat{\mathbf{x}}(0)-\hat{\mathbf{y}}(1)=-\hat{\mathbf{y}}$. The remaining cases are checked the same way. The identity reappears later with $\nabla$ in place of $\mathbf{A}$ and $\mathbf{B}$: $\nabla\times(\nabla\times\mathbf{F})=\nabla(\nabla\cdot\mathbf{F})-\nabla^2\mathbf{F}$ (Section 7.5), the step that turns Maxwell's equations into the wave equation.

> [!example] Worked Example — angle, projection and area
> Let $\mathbf{A}=(1,2,2)$ and $\mathbf{B}=(3,0,4)$.
>
> **Angle.** $\mathbf{A}\cdot\mathbf{B}=3+0+8=11$, $A=3$, $B=5$, so $\cos\theta=11/15$ and $\theta=42.83^\circ$.
>
> **Projection of $\mathbf{A}$ on $\mathbf{B}$.** $\mathbf{A}\cdot\hat{\mathbf{B}}=11/5=2.2$; vector projection $2.2\,(0.6,0,0.8)=(1.32,0,1.76)$. The perpendicular remainder is $(1,2,2)-(1.32,0,1.76)=(-0.32,2,0.24)$; check $(-0.32)(3)+0+(0.24)(4)=0$. ✓
>
> **Cross product.** $\mathbf{A}\times\mathbf{B}=(2\cdot4-2\cdot0,\ 2\cdot3-1\cdot4,\ 1\cdot0-2\cdot3)=(8,2,-6)$.
> Its length $\sqrt{64+4+36}=\sqrt{104}=10.198$ should equal $AB\sin\theta=15\sin42.83^\circ=15\times0.6799=10.198$. ✓
> So the triangle with sides $\mathbf{A}$ and $\mathbf{B}$ has area $\tfrac12\sqrt{104}=5.10$.

> [!example] Worked Example — volume and coplanarity
> Are $\mathbf{A}=(1,0,1)$, $\mathbf{B}=(2,1,0)$, $\mathbf{C}=(0,1,-2)$ coplanar?
> $$
> \begin{vmatrix}1&0&1\\2&1&0\\0&1&-2\end{vmatrix}=1(1\cdot(-2)-0\cdot1)-0+1(2\cdot1-1\cdot0)=-2+2=0 .
> $$
> Yes: $\mathbf{C}=\mathbf{B}-2\mathbf{A}$. Replacing $\mathbf{C}$ by $(0,1,2)$ gives determinant $1(2-0)+1(2)=4$, a parallelepiped of volume $4$.

> [!example] Worked Example — Lorentz force on an electron
> An electron ($q=-1.602\times10^{-19}$ C) moves with $\mathbf{v}=(2\times10^{5},\ 0,\ 1\times10^{5})$ m/s through $\mathbf{E}=(0,\ 500,\ 0)$ V/m and $\mathbf{B}=(0,\ 0,\ 0.01)$ T. Find $\mathbf{F}=q(\mathbf{E}+\mathbf{v}\times\mathbf{B})$.
>
> **Step 1.** $\mathbf{v}\times\mathbf{B}=(v_yB_z-v_zB_y,\ v_zB_x-v_xB_z,\ v_xB_y-v_yB_x)=(0,\ -2\times10^{3},\ 0)$ V/m. (Units: $\mathrm{m/s\cdot T}=\mathrm{V/m}$.)
>
> **Step 2.** $\mathbf{E}+\mathbf{v}\times\mathbf{B}=(0,\ 500-2000,\ 0)=(0,-1500,0)$ V/m.
>
> **Step 3.** $\mathbf{F}=-1.602\times10^{-19}\cdot(0,-1500,0)=(0,\ 2.403\times10^{-16},\ 0)$ N.
>
> The $v_z$ component, parallel to $\mathbf{B}$, contributes nothing to the magnetic force — a cross product only sees the perpendicular part. Choosing $\mathbf{E}=(0,-v_xB_z,0)$ would cancel the force entirely: this is the velocity-selector (Wien filter) condition $E=vB$.

> [!warning] Common Pitfall
> The cross product is defined only in three dimensions and depends on handedness. If you set up a left-handed axis system (for instance by drawing $\hat{\mathbf{y}}$ the wrong way on a page), every cross product — and every curl — changes sign. Always check $\hat{\mathbf{x}}\times\hat{\mathbf{y}}=\hat{\mathbf{z}}$ for your drawing.

> [!tip] Engineering Note
> The dot product is the continuous cousin of the inner product used in signal processing ([Ch. 9](../part-3-linear-algebra/09-vector-spaces-and-linear-transformations.md)): correlating two signals is "dotting" them. For complex phasor fields the time-average power uses a conjugate, $\langle\mathbf{S}\rangle=\tfrac12\operatorname{Re}\{\mathbf{E}\times\mathbf{H}^*\}$; forgetting the conjugate is a classic error.

## 7.3 Vector-valued functions

### Definition and derivative

A trajectory, a wire, the tip of a rotating phasor, the path of an electron in a magnetron — each is a point moving through space as a parameter changes. We describe it with a vector-valued function.

> [!definition] Vector-valued function
> A **vector-valued function** of one real variable assigns a vector to each $t$ in an interval:
> $$
> \mathbf{r}(t)=x(t)\hat{\mathbf{x}}+y(t)\hat{\mathbf{y}}+z(t)\hat{\mathbf{z}} .
> $$
> Its **derivative** is
> $$
> \mathbf{r}'(t)=\lim_{h\to0}\frac{\mathbf{r}(t+h)-\mathbf{r}(t)}{h}=x'(t)\hat{\mathbf{x}}+y'(t)\hat{\mathbf{y}}+z'(t)\hat{\mathbf{z}},
> $$
> and its integral is likewise taken component by component. The curve is **smooth** if $\mathbf{r}'$ is continuous and never zero.

The componentwise rule works *only because* $\hat{\mathbf{x}},\hat{\mathbf{y}},\hat{\mathbf{z}}$ are constant. Geometrically, $\mathbf{r}(t+h)-\mathbf{r}(t)$ is a chord of the curve; dividing by $h$ and letting $h\to0$ turns it into a **tangent vector**. If $t$ is time, $\mathbf{r}'=\mathbf{v}$ is the velocity, $\lVert\mathbf{v}\rVert$ the speed, and $\mathbf{r}''=\mathbf{a}$ the acceleration.

> [!theorem] Differentiation rules
> For differentiable vector functions $\mathbf{A}(t),\mathbf{B}(t)$ and scalar $f(t)$:
> $$
> \begin{aligned}
> (\mathbf{A}+\mathbf{B})'&=\mathbf{A}'+\mathbf{B}', & (f\mathbf{A})'&=f'\mathbf{A}+f\mathbf{A}',\\
> (\mathbf{A}\cdot\mathbf{B})'&=\mathbf{A}'\cdot\mathbf{B}+\mathbf{A}\cdot\mathbf{B}', & (\mathbf{A}\times\mathbf{B})'&=\mathbf{A}'\times\mathbf{B}+\mathbf{A}\times\mathbf{B}',\\
> \frac{d}{dt}\mathbf{A}(u(t))&=\mathbf{A}'(u)\,u'(t). &&
> \end{aligned}
> $$
> In the cross-product rule the **order matters**.

Each rule follows by writing the product in components and applying the scalar product rule. A useful corollary: if $\lVert\mathbf{A}(t)\rVert$ is constant then $\mathbf{A}\cdot\mathbf{A}=\text{const}$, so $2\mathbf{A}\cdot\mathbf{A}'=0$ — **a vector of constant length is perpendicular to its own derivative.** A unit vector can only rotate, never stretch.

### Arc length, unit tangent and curvature

The distance travelled along the curve from $t=a$ to $t=b$ is the **arc length**

$$
s=\int_a^b\lVert\mathbf{r}'(t)\rVert\,dt=\int_a^b\sqrt{x'(t)^2+y'(t)^2+z'(t)^2}\,dt ,
$$

because a small step $d\mathbf{r}=\mathbf{r}'\,dt$ has length $ds=\lVert\mathbf{r}'\rVert\,dt$. The **unit tangent** is $\hat{\mathbf{T}}=\mathbf{r}'/\lVert\mathbf{r}'\rVert$, and the **curvature** $\kappa=\lVert d\hat{\mathbf{T}}/ds\rVert$ measures how fast the direction turns per metre. A convenient formula is

$$
\kappa=\frac{\lVert\mathbf{r}'\times\mathbf{r}''\rVert}{\lVert\mathbf{r}'\rVert^3}.
$$

A straight line has $\kappa=0$; a circle of radius $a$ has $\kappa=1/a$. The arc length is independent of how the curve is parametrised: substituting $t=g(u)$ with $g'>0$ multiplies $\lVert\mathbf{r}'\rVert$ by $g'(u)$ and $dt$ by $1/g'(u)$ in reverse — the same cancellation that will make line integrals parametrisation-independent in Section 7.6.

> [!example] Worked Example — length of a solenoid winding
> A solenoid has radius $a=5$ mm and pitch $p=1$ mm per turn, with 100 turns. Find the wire length and the curvature of the wire.
>
> **Step 1 — parametrise.** One turn is $2\pi$ of angle and advances $p$ along $z$: $\mathbf{r}(t)=(a\cos t,\ a\sin t,\ bt)$ with $b=p/(2\pi)=0.15915$ mm.
>
> **Step 2 — speed.** $\mathbf{r}'=(-a\sin t,\ a\cos t,\ b)$, $\lVert\mathbf{r}'\rVert=\sqrt{a^2+b^2}=\sqrt{25+0.02533}=5.00253$ mm (constant).
>
> **Step 3 — length.** $s=\int_0^{200\pi}\sqrt{a^2+b^2}\,dt=200\pi\times5.00253=3143.2$ mm $=3.1432$ m. A flat-circle estimate $100\times2\pi a=3.1416$ m is only 0.05 % low, because the pitch is small compared with the circumference.
>
> **Step 4 — curvature.** $\mathbf{r}''=(-a\cos t,-a\sin t,0)$ and $\mathbf{r}'\times\mathbf{r}''=(ab\sin t,\ -ab\cos t,\ a^2)$ with length $a\sqrt{a^2+b^2}$. So $\kappa=a/(a^2+b^2)=0.19980\ \mathrm{mm}^{-1}$, a radius of curvature $5.005$ mm — very slightly larger than the coil radius.

> [!example] Worked Example — arc length that needs care
> Find the length of one arch of the cycloid $\mathbf{r}(t)=(t-\sin t,\ 1-\cos t)$, $0\le t\le2\pi$ — the path of a point on the rim of a rolling unit wheel.
>
> **Step 1.** $\mathbf{r}'=(1-\cos t,\ \sin t)$, so $\lVert\mathbf{r}'\rVert^2=1-2\cos t+\cos^2t+\sin^2t=2-2\cos t=4\sin^2(t/2)$.
>
> **Step 2.** On $[0,2\pi]$, $\sin(t/2)\ge0$, so $\lVert\mathbf{r}'\rVert=2\sin(t/2)$ — no absolute-value trouble.
>
> **Step 3.** $s=\int_0^{2\pi}2\sin(t/2)\,dt=\big[-4\cos(t/2)\big]_0^{2\pi}=4+4=8$.
>
> Note that $\mathbf{r}'(0)=\mathbf{0}$: the cycloid has a cusp there and is not smooth at the ends, but the arc-length integral is still perfectly finite.

```desmos-graph
left=-0.5; right=6.8; bottom=-0.5; top=2.6
---
(t-\sin(t),1-\cos(t))|blue|0<t<6.2832
(t,1)|grey|dashed|0<t<6.2832
(\pi,2)|label:top of arch
```

*One arch of the cycloid traced by a point on a rolling unit wheel (axle path dashed). Its length is exactly 8, four times the wheel's diameter.*

### Moving bases: velocity in polar and cylindrical coordinates

If a point moves with cylindrical coordinates $(r(t),\phi(t),z(t))$, its position is $\mathbf{r}=r\hat{\mathbf{r}}+z\hat{\mathbf{z}}$ — but $\hat{\mathbf{r}}$ depends on $\phi$, which depends on $t$. Differentiating the expressions of Section 7.1 with the chain rule,

$$
\frac{d\hat{\mathbf{r}}}{dt}=\dot\phi\,(-\sin\phi\,\hat{\mathbf{x}}+\cos\phi\,\hat{\mathbf{y}})=\dot\phi\,\hat{\boldsymbol{\phi}},\qquad
\frac{d\hat{\boldsymbol{\phi}}}{dt}=-\dot\phi\,\hat{\mathbf{r}} .
$$

(Both derivatives are perpendicular to the vector being differentiated, as the constant-length corollary demands.) Hence

$$
\begin{aligned}
\mathbf{v}&=\dot r\,\hat{\mathbf{r}}+r\dot\phi\,\hat{\boldsymbol{\phi}}+\dot z\,\hat{\mathbf{z}},\\
\mathbf{a}&=(\ddot r-r\dot\phi^2)\,\hat{\mathbf{r}}+(r\ddot\phi+2\dot r\dot\phi)\,\hat{\boldsymbol{\phi}}+\ddot z\,\hat{\mathbf{z}} .
\end{aligned}
$$

The term $-r\dot\phi^2\hat{\mathbf{r}}$ is the centripetal acceleration and $2\dot r\dot\phi\,\hat{\boldsymbol{\phi}}$ the Coriolis term. Neither appears if you forget that the basis vectors move. The speed formula $\lVert\mathbf{v}\rVert^2=\dot r^2+r^2\dot\phi^2+\dot z^2$ is the time-derivative version of the cylindrical line element $ds^2=dr^2+r^2d\phi^2+dz^2$.

> [!example] Worked Example — an electron in a uniform magnetic field
> A particle of charge $q$ and mass $m$ moves in $\mathbf{B}=B_0\hat{\mathbf{z}}$ with no electric field. Newton's law is $m\,\dot{\mathbf{v}}=q\,\mathbf{v}\times\mathbf{B}$. Show that the motion is a helix and find its radius and frequency for an electron with $v_\perp=10^6$ m/s in $B_0=0.01$ T.
>
> **Step 1 — energy is conserved.** $\frac{d}{dt}\big(\tfrac12m\,\mathbf{v}\cdot\mathbf{v}\big)=m\,\dot{\mathbf{v}}\cdot\mathbf{v}=q(\mathbf{v}\times\mathbf{B})\cdot\mathbf{v}=0$ by the scalar triple product (two equal vectors). A magnetic field never changes the speed.
>
> **Step 2 — components.** $\mathbf{v}\times\mathbf{B}=(v_yB_0,\ -v_xB_0,\ 0)$. With $\omega_c=qB_0/m$: $\dot v_x=\omega_cv_y$, $\dot v_y=-\omega_cv_x$, $\dot v_z=0$.
>
> **Step 3 — solve.** Try $v_x=v_\perp\cos\omega_ct$, $v_y=-v_\perp\sin\omega_ct$: then $\dot v_x=-v_\perp\omega_c\sin\omega_ct=\omega_cv_y$ ✓ and $\dot v_y=-v_\perp\omega_c\cos\omega_ct=-\omega_cv_x$ ✓. Integrating, the orbit is a circle of radius $r_c=v_\perp/\lvert\omega_c\rvert$ while $z$ grows uniformly — a helix like the solenoid winding above.
>
> **Step 4 — numbers.** $\lvert\omega_c\rvert=1.602\times10^{-19}\times0.01/9.109\times10^{-31}=1.759\times10^{9}$ rad/s, i.e. $f_c=\lvert\omega_c\rvert/2\pi=279.9$ MHz. Radius $r_c=10^6/1.759\times10^9=0.569$ mm.
>
> **Step 5 — check with the polar formula.** For uniform circular motion $\dot r=\ddot r=0$, $\ddot\phi=0$, so $\mathbf{a}=-r\dot\phi^2\hat{\mathbf{r}}$: magnitude $r_c\omega_c^2=v_\perp\lvert\omega_c\rvert=1.759\times10^{15}\ \mathrm{m/s^2}$, which equals $\lvert q\rvert v_\perp B_0/m$ ✓.
>
> The cyclotron frequency $\omega_c=qB/m$ is independent of speed — the fact that makes cyclotrons and magnetrons work.

> [!warning] Common Pitfall
> "Differentiate component by component" is valid only for a **fixed** basis. Writing $\frac{d}{dt}(r\hat{\mathbf{r}})=\dot r\hat{\mathbf{r}}$ loses the $r\dot\phi\,\hat{\boldsymbol{\phi}}$ term — the whole of the velocity for circular motion.

> [!tip] Engineering Note
> A phasor $V_m e^{j\omega t}$ is a vector-valued function in the complex plane: $\mathbf{r}(t)=(V_m\cos\omega t,\ V_m\sin\omega t)$. Its derivative $j\omega V_me^{j\omega t}$ is the tangent vector — rotated $90^\circ$ ahead and scaled by $\omega$ — which is why differentiation becomes multiplication by $j\omega$ in AC analysis ([Ch. 2](../part-1-foundations/02-complex-numbers.md)). The polarisation of a plane wave is likewise the curve traced by the tip of $\mathbf{E}(t)$ at a fixed point: a line, circle or ellipse.

## 7.4 Vector fields

### Definition and examples

> [!definition] Scalar and vector fields
> A **scalar field** assigns a number $f(\mathbf{R})$ to each point of a region of space. A **vector field** assigns a vector
> $$
> \mathbf{F}(x,y,z)=F_x(x,y,z)\,\hat{\mathbf{x}}+F_y(x,y,z)\,\hat{\mathbf{y}}+F_z(x,y,z)\,\hat{\mathbf{z}}
> $$
> to each point. A field is **static** if it does not depend on time and **time-varying** otherwise; it is **uniform** if it is the same at every point.

Potential $V$, temperature $T$ and charge density $\rho$ are scalar fields; $\mathbf{E}$, $\mathbf{B}$, current density $\mathbf{J}$ and heat flux $\mathbf{q}$ are vector fields. Picture a vector field as an arrow attached to every point, with length showing magnitude. The components are ordinary scalar functions of several variables, so everything in [Ch. 6](../part-2-calculus/06-multivariable-calculus.md) — continuity, partial derivatives, the chain rule — applies to each component separately.

A short catalogue of fields that recur throughout the chapter:

| Field | Formula | Picture |
|---|---|---|
| Uniform | $\mathbf{F}=E_0\hat{\mathbf{x}}$ | parallel arrows, all equal (parallel-plate capacitor) |
| Radial (3-D) | $\mathbf{F}=\mathbf{R}=x\hat{\mathbf{x}}+y\hat{\mathbf{y}}+z\hat{\mathbf{z}}$ | arrows point outward, growing with distance |
| Point charge | $\mathbf{E}=\dfrac{q}{4\pi\varepsilon_0}\dfrac{\hat{\mathbf{R}}}{R^2}$ | outward, decaying as $1/R^2$ |
| Line charge | $\mathbf{E}=\dfrac{\rho_\ell}{2\pi\varepsilon_0 r}\hat{\mathbf{r}}$ | outward from the $z$-axis, decaying as $1/r$ |
| Rigid rotation | $\mathbf{F}=-y\hat{\mathbf{x}}+x\hat{\mathbf{y}}=r\hat{\boldsymbol{\phi}}$ | circles, speed growing with $r$ |
| Line current | $\mathbf{B}=\dfrac{\mu_0I}{2\pi r}\hat{\boldsymbol{\phi}}$ | circles, strength decaying as $1/r$ |

The last two both circulate around the $z$-axis, but (as we shall see) one has curl and the other does not. Pictures alone can mislead; the differential operators of Section 7.5 make the distinction precise.

### Field lines

> [!definition] Field line (streamline, flux line)
> A **field line** of $\mathbf{F}$ is a curve that is tangent to $\mathbf{F}$ at every point. If $\mathbf{r}(s)$ is a field line then $d\mathbf{r}/ds$ is parallel to $\mathbf{F}(\mathbf{r})$, i.e.
> $$
> \frac{dx}{F_x}=\frac{dy}{F_y}=\frac{dz}{F_z}.
> $$

In two dimensions this is the first-order ODE $dy/dx=F_y/F_x$, solved by the methods of [Ch. 11](../part-4-differential-equations/11-ordinary-differential-equations.md). Field lines never cross where $\mathbf{F}\neq\mathbf{0}$ (the tangent direction would be ambiguous); they can start or stop only where $\mathbf{F}$ vanishes or is singular — at charges, for instance. Their density is drawn proportional to the field magnitude, which is why "lines of force" became Faraday's mental model of the field.

> [!example] Worked Example — field lines of three planar fields
> Find the field lines of (a) $\mathbf{F}=-y\hat{\mathbf{x}}+x\hat{\mathbf{y}}$, (b) $\mathbf{G}=x\hat{\mathbf{x}}-y\hat{\mathbf{y}}$, (c) $\mathbf{H}=\hat{\mathbf{x}}+2x\hat{\mathbf{y}}$.
>
> **(a)** $dy/dx=x/(-y)$, so $y\,dy=-x\,dx$ and $x^2+y^2=C$: **circles** about the origin, traversed counter-clockwise (at $(1,0)$ the field is $+\hat{\mathbf{y}}$).
>
> **(b)** $dy/dx=-y/x$, so $dy/y=-dx/x$, $\ln\lvert y\rvert=-\ln\lvert x\rvert+c$ and $xy=C$: **rectangular hyperbolas**. Flow comes in along the $y$-axis and leaves along the $x$-axis — the field near a saddle point of a potential, e.g. the centre of a quadrupole lens.
>
> **(c)** $dy/dx=2x$, so $y=x^2+C$: a family of **parabolas** shifted vertically.

```desmos-graph
left=-3; right=3; bottom=-3; top=3
---
x^2+y^2=1|blue
x^2+y^2=4|blue
xy=1|red
xy=-1|red
xy=3|red|dashed
xy=-3|red|dashed
```

*Field lines of the rotation field $-y\hat{\mathbf{x}}+x\hat{\mathbf{y}}$ (blue circles) and of the saddle field $x\hat{\mathbf{x}}-y\hat{\mathbf{y}}$ (red hyperbolas $xy=C$).*

> [!example] Worked Example — field lines and equipotentials of a dipole
> A small electric dipole on the $z$-axis has, in spherical coordinates,
> $$
> \mathbf{E}=\frac{p}{4\pi\varepsilon_0R^3}\left(2\cos\theta\,\hat{\mathbf{R}}+\sin\theta\,\hat{\boldsymbol{\theta}}\right),\qquad V=\frac{p\cos\theta}{4\pi\varepsilon_0R^2}.
> $$
> Find the field lines.
>
> **Step 1 — the field-line condition in spherical coordinates.** A displacement along a curve in the plane $\phi=\text{const}$ is $dR\,\hat{\mathbf{R}}+R\,d\theta\,\hat{\boldsymbol{\theta}}$ (the scale factors of [Ch. 6](../part-2-calculus/06-multivariable-calculus.md)). Parallel to $\mathbf{E}$ means
> $$
> \frac{dR}{E_R}=\frac{R\,d\theta}{E_\theta}\quad\Longrightarrow\quad\frac{dR}{2\cos\theta}=\frac{R\,d\theta}{\sin\theta}.
> $$
>
> **Step 2 — separate.** $\dfrac{dR}{R}=\dfrac{2\cos\theta}{\sin\theta}d\theta$, so $\ln R=2\ln\sin\theta+c$ and
> $$
> R=C\sin^2\theta .
> $$
>
> **Step 3 — check orthogonality.** Equipotentials are $R^2=K\cos\theta$. On them $2R\,dR=-K\sin\theta\,d\theta$, so $dR/(R\,d\theta)=-\tfrac12\tan\theta$; on field lines $dR/(R\,d\theta)=2\cot\theta$. The product of these two slopes is $-1$, so the families cross at right angles, as $\mathbf{E}=-\nabla V$ requires.
>
> In the Desmos plot below the dipole axis is drawn vertically; with Desmos's polar angle $\vartheta$ measured from the horizontal, $\sin\theta=\cos\vartheta$ and $\cos\theta=\sin\vartheta$.

```desmos-graph
left=-2.5; right=2.5; bottom=-2; top=2
---
r=0.5\cos(\theta)^2|blue
r=\cos(\theta)^2|blue
r=2\cos(\theta)^2|blue
r=4\cos(\theta)^2|blue
r=\sqrt{0.8\sin(\theta)}|red|dashed
r=\sqrt{2\sin(\theta)}|red|dashed
r=\sqrt{-0.8\sin(\theta)}|red|dashed
r=\sqrt{-2\sin(\theta)}|red|dashed
```

*Dipole with its axis vertical: field lines $R=C\sin^2\theta$ (blue) leave the upper (positive) end and return to the lower end, crossing the equipotentials $R^2=K\cos\theta$ (red, dashed) at right angles.*

> [!example] Worked Example — the magnetic field of a wire in Cartesian form
> Convert $\mathbf{B}=\dfrac{\mu_0I}{2\pi r}\hat{\boldsymbol{\phi}}$ to Cartesian components and evaluate it at $(3,4,0)$ cm for $I=10$ A.
>
> **Step 1.** $\hat{\boldsymbol{\phi}}=(-\sin\phi,\cos\phi,0)=(-y/r,\ x/r,\ 0)$, so
> $$
> \mathbf{B}=\frac{\mu_0I}{2\pi}\,\frac{-y\,\hat{\mathbf{x}}+x\,\hat{\mathbf{y}}}{x^2+y^2}.
> $$
>
> **Step 2.** At $(0.03,0.04)$ m, $r=0.05$ m and $\mu_0I/(2\pi r)=(4\pi\times10^{-7}\times10)/(2\pi\times0.05)=40\ \mu\mathrm{T}$. Direction $\hat{\boldsymbol{\phi}}=(-0.8,0.6,0)$, so $\mathbf{B}=(-32,\ 24,\ 0)\ \mu\mathrm{T}$.
>
> The field lines are again circles, exactly like the rotation field $-y\hat{\mathbf{x}}+x\hat{\mathbf{y}}$, but the magnitude falls off as $1/r$ instead of growing as $r$. That difference will decide whether the curl is zero (Section 7.5). For comparison, the Earth's field is about $50\ \mu\mathrm{T}$ — which is why a 10 A cable a few centimetres from a compass upsets it.

> [!warning] Common Pitfall
> A field-line diagram shows direction and (by density) magnitude, but **not** the path a charged particle follows. A particle's acceleration, not its velocity, is along $\mathbf{E}$; except when it starts from rest in a straight-line field, its trajectory leaves the field line.

> [!tip] Engineering Note
> Finite-element field solvers (for PCB crosstalk, motor design, antenna feeds) compute a field at mesh nodes and display it as arrow plots, streamlines and colour maps of magnitude. Reading those plots is reading a vector field: where lines crowd, the field (and the risk of dielectric breakdown) is high; where they circulate, there is curl, hence current or changing flux.

## 7.5 Gradient, divergence, and curl

### The del operator

All three local derivatives of a field are built from one symbolic vector, the **del** (or nabla) operator

$$
\nabla=\hat{\mathbf{x}}\frac{\partial}{\partial x}+\hat{\mathbf{y}}\frac{\partial}{\partial y}+\hat{\mathbf{z}}\frac{\partial}{\partial z}.
$$

$\nabla$ is not a vector in its own right — it has no magnitude until it acts on something — but in Cartesian coordinates it can be combined with a field using the three products we already know:

| Operation | Acts on | Produces | Cartesian form | Measures |
|---|---|---|---|---|
| gradient $\nabla f$ | scalar | vector | $(f_x,\ f_y,\ f_z)$ | rate and direction of steepest increase |
| divergence $\nabla\cdot\mathbf{F}$ | vector | scalar | $\partial_xF_x+\partial_yF_y+\partial_zF_z$ | net outflow per unit volume |
| curl $\nabla\times\mathbf{F}$ | vector | vector | determinant below | circulation per unit area (rotation) |
| Laplacian $\nabla^2 f=\nabla\cdot\nabla f$ | scalar | scalar | $f_{xx}+f_{yy}+f_{zz}$ | how $f$ differs from its local average |

Each operator has a dimension of $\mathrm{m^{-1}}$ (or $\mathrm{m^{-2}}$ for the Laplacian). If $V$ is in volts, $\nabla V$ is in V/m; if $\mathbf{D}$ is in $\mathrm{C/m^2}$, $\nabla\cdot\mathbf{D}$ is in $\mathrm{C/m^3}$, a charge density — the first hint of Gauss's law.

The best way to understand the three operators is through their **coordinate-free definitions** as limits of integrals over shrinking regions. These definitions explain what the operators mean, tell you how to compute them in any coordinate system, and are exactly the seeds from which the integral theorems of Sections 7.9–7.11 grow.

### Gradient

The gradient was introduced in [Ch. 6](../part-2-calculus/06-multivariable-calculus.md), §6.4: $df=\nabla f\cdot d\mathbf{r}$, $\nabla f$ points uphill, is normal to level surfaces, and its length is the maximum rate of increase. In integral language, for any path from $\mathbf{a}$ to $\mathbf{b}$,

$$
f(\mathbf{b})-f(\mathbf{a})=\int_{\mathbf{a}}^{\mathbf{b}}\nabla f\cdot d\mathbf{r},
$$

which is the fundamental theorem for line integrals (Section 7.6).

**In curvilinear coordinates.** A step $du_i$ in coordinate $u_i$ moves a distance $h_i\,du_i$, where the $h_i$ are the **scale factors** of Ch. 6: $(1,r,1)$ for cylindrical $(r,\phi,z)$ and $(1,R,R\sin\theta)$ for spherical $(R,\theta,\phi)$. The component of $\nabla f$ along $\hat{\mathbf{u}}_i$ is the rate of change of $f$ *per metre* in that direction, $\frac{1}{h_i}\frac{\partial f}{\partial u_i}$. Hence

$$
\begin{aligned}
\nabla f&=\frac{\partial f}{\partial r}\hat{\mathbf{r}}+\frac1r\frac{\partial f}{\partial\phi}\hat{\boldsymbol{\phi}}+\frac{\partial f}{\partial z}\hat{\mathbf{z}} &&\text{(cylindrical)},\\
\nabla f&=\frac{\partial f}{\partial R}\hat{\mathbf{R}}+\frac1R\frac{\partial f}{\partial\theta}\hat{\boldsymbol{\theta}}+\frac{1}{R\sin\theta}\frac{\partial f}{\partial\phi}\hat{\boldsymbol{\phi}} &&\text{(spherical)}.
\end{aligned}
$$

### Divergence

> [!definition] Divergence
> The **divergence** of $\mathbf{F}$ at a point $P$ is the net outward flux per unit volume of a small closed surface $S$ enclosing $P$:
> $$
> \nabla\cdot\mathbf{F}=\lim_{\Delta V\to0}\frac{1}{\Delta V}\oint_S\mathbf{F}\cdot d\mathbf{S}.
> $$
> In Cartesian coordinates
> $$
> \nabla\cdot\mathbf{F}=\frac{\partial F_x}{\partial x}+\frac{\partial F_y}{\partial y}+\frac{\partial F_z}{\partial z}.
> $$

**Derivation of the Cartesian formula.** Take a small box with one corner at $(x,y,z)$ and sides $\Delta x,\Delta y,\Delta z$. The flux out through the face at $x+\Delta x$ is approximately $F_x(x+\Delta x,y,z)\,\Delta y\,\Delta z$; through the opposite face at $x$ the outward normal is $-\hat{\mathbf{x}}$, giving $-F_x(x,y,z)\,\Delta y\,\Delta z$. Together

$$
\big[F_x(x+\Delta x,y,z)-F_x(x,y,z)\big]\Delta y\,\Delta z\approx\frac{\partial F_x}{\partial x}\,\Delta x\,\Delta y\,\Delta z .
$$

The other two pairs of faces contribute $\partial F_y/\partial y\,\Delta V$ and $\partial F_z/\partial z\,\Delta V$. Dividing by $\Delta V$ and letting the box shrink gives the formula. The idea: only the *change* of $F_x$ across the box in the $x$ direction produces net outflow; a uniform field carries in as much as it carries out.

**Interpretation.** $\nabla\cdot\mathbf{F}>0$ at a **source** (field lines start there, net outflow), $<0$ at a **sink**, and $=0$ where the field is **solenoidal** (whatever flows in flows out). For an electric field, sources are positive charges; for the velocity of an incompressible fluid, the divergence is zero everywhere; for current density, $\nabla\cdot\mathbf{J}=-\partial\rho/\partial t$ (charge conservation, Section 7.11).

**In curvilinear coordinates.** Repeat the box argument with a curvilinear box whose faces have areas $h_2h_3\,du_2\,du_3$, and so on. The faces now have different areas on the two sides, so the area factors go *inside* the derivatives:

$$
\nabla\cdot\mathbf{F}=\frac{1}{h_1h_2h_3}\left[\frac{\partial}{\partial u_1}(h_2h_3F_1)+\frac{\partial}{\partial u_2}(h_3h_1F_2)+\frac{\partial}{\partial u_3}(h_1h_2F_3)\right].
$$

In particular

$$
\begin{aligned}
\nabla\cdot\mathbf{F}&=\frac1r\frac{\partial}{\partial r}(rF_r)+\frac1r\frac{\partial F_\phi}{\partial\phi}+\frac{\partial F_z}{\partial z} &&\text{(cylindrical)},\\
\nabla\cdot\mathbf{F}&=\frac{1}{R^2}\frac{\partial}{\partial R}(R^2F_R)+\frac{1}{R\sin\theta}\frac{\partial}{\partial\theta}(\sin\theta\,F_\theta)+\frac{1}{R\sin\theta}\frac{\partial F_\phi}{\partial\phi} &&\text{(spherical)}.
\end{aligned}
$$

The factor $R^2$ in the spherical radial term is the area of a sphere growing as $R^2$: a field $\propto1/R^2$ spreads over exactly that growing area and so has zero divergence away from the origin.

### Curl

> [!definition] Curl
> The component of the **curl** of $\mathbf{F}$ along a unit vector $\hat{\mathbf{n}}$ at $P$ is the circulation per unit area around a small loop $C$ through $P$ in the plane perpendicular to $\hat{\mathbf{n}}$, traversed counter-clockwise when viewed from the tip of $\hat{\mathbf{n}}$ (right-hand rule):
> $$
> (\nabla\times\mathbf{F})\cdot\hat{\mathbf{n}}=\lim_{\Delta S\to0}\frac{1}{\Delta S}\oint_C\mathbf{F}\cdot d\mathbf{r}.
> $$
> In Cartesian coordinates
> $$
> \nabla\times\mathbf{F}=\begin{vmatrix}\hat{\mathbf{x}}&\hat{\mathbf{y}}&\hat{\mathbf{z}}\\\partial_x&\partial_y&\partial_z\\F_x&F_y&F_z\end{vmatrix}
> =\left(\frac{\partial F_z}{\partial y}-\frac{\partial F_y}{\partial z}\right)\hat{\mathbf{x}}+\left(\frac{\partial F_x}{\partial z}-\frac{\partial F_z}{\partial x}\right)\hat{\mathbf{y}}+\left(\frac{\partial F_y}{\partial x}-\frac{\partial F_x}{\partial y}\right)\hat{\mathbf{z}} .
> $$

**Derivation of the $z$-component.** Take a small rectangle in the plane $z=\text{const}$ with corner $(x,y)$ and sides $\Delta x,\Delta y$, traversed counter-clockwise seen from $+z$. The four sides contribute approximately

$$
F_x(x,y)\Delta x+F_y(x+\Delta x,y)\Delta y-F_x(x,y+\Delta y)\Delta x-F_y(x,y)\Delta y .
$$

Pairing the opposite sides,

$$
\big[F_y(x+\Delta x,y)-F_y(x,y)\big]\Delta y-\big[F_x(x,y+\Delta y)-F_x(x,y)\big]\Delta x\approx\left(\frac{\partial F_y}{\partial x}-\frac{\partial F_x}{\partial y}\right)\Delta x\,\Delta y .
$$

Divide by the area $\Delta x\,\Delta y$ and take the limit. The $x$- and $y$-components follow by cyclic permutation $x\to y\to z\to x$.

**Interpretation.** Put a tiny paddle wheel in a flowing fluid with velocity $\mathbf{v}$. It spins if the flow on one side is faster than on the other; its axis aligns with $\nabla\times\mathbf{v}$ and its angular velocity is $\tfrac12\lVert\nabla\times\mathbf{v}\rVert$. Curl is about **shear across the flow**, not about the field lines being curved: a field can circulate without curl (the wire field, below) or have curl with straight field lines (a shear flow $\mathbf{F}=y\hat{\mathbf{x}}$ has curl $-\hat{\mathbf{z}}$).

**In curvilinear coordinates.** The loop argument with sides $h_i\,du_i$ gives

$$
\nabla\times\mathbf{F}=\frac{1}{h_1h_2h_3}\begin{vmatrix}h_1\hat{\mathbf{u}}_1&h_2\hat{\mathbf{u}}_2&h_3\hat{\mathbf{u}}_3\\\partial_{u_1}&\partial_{u_2}&\partial_{u_3}\\h_1F_1&h_2F_2&h_3F_3\end{vmatrix},
$$

which written out is

$$
\begin{aligned}
\nabla\times\mathbf{F}&=\left(\frac1r\frac{\partial F_z}{\partial\phi}-\frac{\partial F_\phi}{\partial z}\right)\hat{\mathbf{r}}+\left(\frac{\partial F_r}{\partial z}-\frac{\partial F_z}{\partial r}\right)\hat{\boldsymbol{\phi}}+\frac1r\left(\frac{\partial(rF_\phi)}{\partial r}-\frac{\partial F_r}{\partial\phi}\right)\hat{\mathbf{z}}\quad\text{(cylindrical)},\\[4pt]
\nabla\times\mathbf{F}&=\frac{1}{R\sin\theta}\left(\frac{\partial(\sin\theta\,F_\phi)}{\partial\theta}-\frac{\partial F_\theta}{\partial\phi}\right)\hat{\mathbf{R}}+\frac1R\left(\frac{1}{\sin\theta}\frac{\partial F_R}{\partial\phi}-\frac{\partial(RF_\phi)}{\partial R}\right)\hat{\boldsymbol{\theta}}+\frac1R\left(\frac{\partial(RF_\theta)}{\partial R}-\frac{\partial F_R}{\partial\theta}\right)\hat{\boldsymbol{\phi}}\quad\text{(spherical)}.
\end{aligned}
$$

### The Laplacian

The divergence of a gradient is the **Laplacian**:

$$
\nabla^2f=\nabla\cdot\nabla f=\frac{\partial^2f}{\partial x^2}+\frac{\partial^2f}{\partial y^2}+\frac{\partial^2f}{\partial z^2}.
$$

Combining the curvilinear gradient and divergence,

$$
\begin{aligned}
\nabla^2f&=\frac1r\frac{\partial}{\partial r}\left(r\frac{\partial f}{\partial r}\right)+\frac1{r^2}\frac{\partial^2f}{\partial\phi^2}+\frac{\partial^2f}{\partial z^2} &&\text{(cylindrical)},\\
\nabla^2f&=\frac1{R^2}\frac{\partial}{\partial R}\left(R^2\frac{\partial f}{\partial R}\right)+\frac{1}{R^2\sin\theta}\frac{\partial}{\partial\theta}\left(\sin\theta\frac{\partial f}{\partial\theta}\right)+\frac{1}{R^2\sin^2\theta}\frac{\partial^2f}{\partial\phi^2} &&\text{(spherical)}.
\end{aligned}
$$

**Interpretation.** In one dimension, $f''(x)\approx[f(x+h)+f(x-h)-2f(x)]/h^2$: the second derivative compares $f$ with the average of its neighbours. In three dimensions, $\nabla^2f$ at a point is proportional to (average of $f$ on a small sphere) $-$ ($f$ at the centre). A function with $\nabla^2f=0$ — a **harmonic** function — equals the average of its surroundings, so it can have no local maximum or minimum inside a region. That is why the electrostatic potential in a charge-free region has its extremes on the electrodes.

For a vector field, the **vector Laplacian** in Cartesian coordinates is the Laplacian of each component: $\nabla^2\mathbf{F}=(\nabla^2F_x,\nabla^2F_y,\nabla^2F_z)$. In curvilinear coordinates it is *not* the Laplacian of each component (the basis vectors vary); use the identity $\nabla^2\mathbf{F}=\nabla(\nabla\cdot\mathbf{F})-\nabla\times(\nabla\times\mathbf{F})$ instead.

### Identities

> [!theorem] Vector-calculus identities
> For twice continuously differentiable $f,g$ and $\mathbf{A},\mathbf{B}$:
> 1. $\nabla\times(\nabla f)=\mathbf{0}$ — *the curl of a gradient is zero.*
> 2. $\nabla\cdot(\nabla\times\mathbf{A})=0$ — *the divergence of a curl is zero.*
> 3. $\nabla(fg)=f\nabla g+g\nabla f$.
> 4. $\nabla\cdot(f\mathbf{A})=f\,\nabla\cdot\mathbf{A}+\mathbf{A}\cdot\nabla f$.
> 5. $\nabla\times(f\mathbf{A})=f\,\nabla\times\mathbf{A}+\nabla f\times\mathbf{A}$.
> 6. $\nabla\cdot(\mathbf{A}\times\mathbf{B})=\mathbf{B}\cdot(\nabla\times\mathbf{A})-\mathbf{A}\cdot(\nabla\times\mathbf{B})$.
> 7. $\nabla\times(\nabla\times\mathbf{A})=\nabla(\nabla\cdot\mathbf{A})-\nabla^2\mathbf{A}$.

*Proof of 1 and 2.* For 1, the $z$-component is $\partial_x\partial_yf-\partial_y\partial_xf=0$ by equality of mixed partials (Clairaut, Ch. 6); the other components are the same. For 2,

$$
\nabla\cdot(\nabla\times\mathbf{A})=\partial_x(\partial_yA_z-\partial_zA_y)+\partial_y(\partial_zA_x-\partial_xA_z)+\partial_z(\partial_xA_y-\partial_yA_x)=0,
$$

since the six terms cancel in pairs. Formally both are "$\mathbf{a}\times\mathbf{a}=\mathbf{0}$" and "$\mathbf{a}\cdot(\mathbf{a}\times\mathbf{b})=0$" with $\mathbf{a}=\nabla$. Identities 4–6 are product rules: each term has $\nabla$ acting on one factor while the other is held fixed. Identity 6 is the one behind Poynting's theorem; identity 7 is BAC–CAB with $\nabla$ and is the step to the wave equation.

These two "zero" identities have converses that are the foundation of potential theory:

> [!theorem] Potentials
> On a **simply connected** region (every closed loop can be shrunk to a point without leaving it):
> - if $\nabla\times\mathbf{F}=\mathbf{0}$, there is a scalar potential $f$ with $\mathbf{F}=\nabla f$ (the field is **irrotational** or **conservative**);
> - on a region where every closed surface encloses only points of the region (e.g. all of space), if $\nabla\cdot\mathbf{F}=0$, there is a vector potential $\mathbf{A}$ with $\mathbf{F}=\nabla\times\mathbf{A}$ (the field is **solenoidal**).
>
> **Helmholtz theorem.** A field that is smooth and vanishes suitably fast at infinity is uniquely determined by its divergence and its curl, and can be split as $\mathbf{F}=-\nabla V+\nabla\times\mathbf{A}$ — an irrotational part plus a solenoidal part.

The Helmholtz theorem explains the shape of Maxwell's equations: to fix $\mathbf{E}$ and $\mathbf{B}$ you need their divergences (the charges and the absence of magnetic charge) and their curls (changing fields and currents) — exactly four equations.

> [!example] Worked Example — divergence and curl in Cartesian coordinates
> For $\mathbf{F}=x^2y\,\hat{\mathbf{x}}+2yz\,\hat{\mathbf{y}}+xz^2\,\hat{\mathbf{z}}$ compute $\nabla\cdot\mathbf{F}$ and $\nabla\times\mathbf{F}$ at $(1,2,3)$, and confirm identity 2.
>
> **Divergence.** $\partial_x(x^2y)+\partial_y(2yz)+\partial_z(xz^2)=2xy+2z+2xz$. At $(1,2,3)$: $4+6+6=16$.
>
> **Curl.**
> $(\partial_yF_z-\partial_zF_y,\ \partial_zF_x-\partial_xF_z,\ \partial_xF_y-\partial_yF_x)=(0-2y,\ 0-z^2,\ 0-x^2)=(-2y,-z^2,-x^2)$. At $(1,2,3)$: $(-4,-9,-1)$.
>
> **Identity check.** $\nabla\cdot(-2y,-z^2,-x^2)=\partial_x(-2y)+\partial_y(-z^2)+\partial_z(-x^2)=0$. ✓
>
> So near $(1,2,3)$ the field is a strong source (net outflow 16 units per unit volume) and its local rotation axis points along $(-4,-9,-1)/\sqrt{98}$.

> [!example] Worked Example — rotation field versus wire field
> Compare $\mathbf{F}_1=r\,\hat{\boldsymbol{\phi}}$ (rigid rotation) and $\mathbf{F}_2=\dfrac1r\hat{\boldsymbol{\phi}}$ (the shape of the field of a line current).
>
> Both have only a $\phi$-component that depends only on $r$, so $\nabla\cdot\mathbf{F}=\frac1r\partial_\phi F_\phi=0$ for both: neither has sources.
>
> For the curl, only the $\hat{\mathbf{z}}$ term survives: $(\nabla\times\mathbf{F})_z=\frac1r\frac{\partial}{\partial r}(rF_\phi)$.
> - $\mathbf{F}_1$: $\frac1r\frac{d}{dr}(r^2)=2$, so $\nabla\times\mathbf{F}_1=2\hat{\mathbf{z}}$ — uniform rotation everywhere (a paddle wheel spins at angular rate 1 = the rotation rate of the fluid, consistent with $\tfrac12\lVert\nabla\times\mathbf{v}\rVert$).
> - $\mathbf{F}_2$: $\frac1r\frac{d}{dr}(1)=0$ for $r>0$ — **irrotational**, even though its field lines are circles. The faster flow near the axis exactly compensates for the shorter path: a paddle wheel is carried around the axis without spinning.
>
> **Cartesian cross-check for $\mathbf{F}_2$.** $\mathbf{F}_2=(-y,x)/(x^2+y^2)$. $\partial_xF_y=\dfrac{(x^2+y^2)-2x^2}{(x^2+y^2)^2}=\dfrac{y^2-x^2}{(x^2+y^2)^2}$ and $\partial_yF_x=\dfrac{-(x^2+y^2)+2y^2}{(x^2+y^2)^2}=\dfrac{y^2-x^2}{(x^2+y^2)^2}$. Their difference is $0$. ✓
>
> Yet in Section 7.6 we will find that the circulation of $\mathbf{F}_2$ around the unit circle is $2\pi$, not zero. The resolution is that $\mathbf{F}_2$ is undefined on the $z$-axis, so the region is not simply connected — all its curl is concentrated on the axis, where the current flows.

```desmos-graph
left=0; right=3; bottom=0; top=3
---
y=x|blue
y=\frac{1}{x}|red
```

*Field magnitude against distance $r$ from the axis: rigid rotation $\lVert\mathbf{F}_1\rVert=r$ (blue) has uniform curl; the wire-type field $\lVert\mathbf{F}_2\rVert=1/r$ (red) is curl-free for $r>0$, although both have circular field lines.*

> [!example] Worked Example — the inverse-square field and $\nabla^2(1/R)$
> Show that $\nabla\cdot(\hat{\mathbf{R}}/R^2)=0$ and $\nabla^2(1/R)=0$ for $R\ne0$, first in spherical and then in Cartesian coordinates.
>
> **Spherical.** $\nabla\cdot(\hat{\mathbf{R}}/R^2)=\dfrac1{R^2}\dfrac{\partial}{\partial R}\left(R^2\cdot\dfrac1{R^2}\right)=\dfrac1{R^2}\dfrac{\partial}{\partial R}(1)=0$. Since $\nabla(1/R)=-\hat{\mathbf{R}}/R^2$, also $\nabla^2(1/R)=0$.
>
> **Cartesian.** With $R^2=x^2+y^2+z^2$, $\partial_x(1/R)=-x/R^3$ and
> $$
> \frac{\partial^2}{\partial x^2}\frac1R=-\frac1{R^3}+\frac{3x^2}{R^5}=\frac{3x^2-R^2}{R^5}.
> $$
> Adding the $y$ and $z$ terms: $\nabla^2(1/R)=\dfrac{3(x^2+y^2+z^2)-3R^2}{R^5}=0$. ✓
>
> The spherical calculation took one line, the Cartesian one five. The flux of $\hat{\mathbf{R}}/R^2$ through a sphere is nevertheless $4\pi$ (Section 7.7): all the divergence sits at $R=0$. Section 7.11 makes this precise with the Dirac delta, $\nabla^2(1/R)=-4\pi\delta(\mathbf{R})$.

> [!example] Worked Example — gradient and divergence in spherical coordinates
> (a) Show that the dipole potential $V=\dfrac{p\cos\theta}{4\pi\varepsilon_0R^2}$ gives the field quoted in Section 7.4. (b) Inside a sphere of uniform charge density $\rho_0$, $\mathbf{E}=\dfrac{\rho_0R}{3\varepsilon_0}\hat{\mathbf{R}}$; find $\nabla\cdot\mathbf{E}$ and $\nabla\times\mathbf{E}$.
>
> **(a)** $E_R=-\partial V/\partial R=\dfrac{2p\cos\theta}{4\pi\varepsilon_0R^3}$, $E_\theta=-\dfrac1R\dfrac{\partial V}{\partial\theta}=\dfrac{p\sin\theta}{4\pi\varepsilon_0R^3}$, $E_\phi=0$. ✓
>
> **(b)** $\nabla\cdot\mathbf{E}=\dfrac1{R^2}\dfrac{\partial}{\partial R}\left(R^2\cdot\dfrac{\rho_0R}{3\varepsilon_0}\right)=\dfrac1{R^2}\cdot\dfrac{3\rho_0R^2}{3\varepsilon_0}=\dfrac{\rho_0}{\varepsilon_0}$ — exactly the charge density divided by $\varepsilon_0$ (Gauss's law in differential form, Section 7.12). The curl has components built from $\partial_\theta E_R$ and $\partial_\phi E_R$ and from $E_\theta,E_\phi$, all zero, so $\nabla\times\mathbf{E}=\mathbf{0}$: the field is conservative, as every electrostatic field must be.
>
> **Units.** $\rho_0/\varepsilon_0$ is $(\mathrm{C/m^3})/(\mathrm{F/m})=\mathrm{V/m^2}$, which is $[\mathbf{E}]/\mathrm{m}$ ✓.

> [!warning] Common Pitfall
> Applying the *Cartesian* formulas to curvilinear components. For $\mathbf{F}=\hat{\mathbf{R}}$ (all $F_R=1$) the naive "$\partial F_R/\partial R$" gives zero divergence, but the correct result is $\frac1{R^2}\partial_R(R^2)=2/R$: unit arrows pointing outward everywhere *do* spread apart. Likewise $\nabla^2\mathbf{F}$ is not $(\nabla^2F_r,\nabla^2F_\phi,\nabla^2F_z)$.

> [!warning] Common Pitfall
> $\nabla\cdot\mathbf{F}$ is a scalar and $\nabla\times\mathbf{F}$ a vector; $\nabla\cdot f$ (divergence of a scalar) and $\nabla\times f$ are meaningless, and $\mathbf{F}\cdot\nabla$ is an *operator*, not the divergence: $(\mathbf{F}\cdot\nabla)g=F_x\partial_xg+F_y\partial_yg+F_z\partial_zg$ is the directional derivative of $g$ along $\mathbf{F}$, scaled by $\lVert\mathbf{F}\rVert$.

> [!tip] Engineering Note
> Finite-difference field solvers replace $\nabla^2V$ by the "five-point stencil" $[V_{i+1,j}+V_{i-1,j}+V_{i,j+1}+V_{i,j-1}-4V_{i,j}]/h^2$. Setting it to zero says each node equals the average of its four neighbours — the discrete version of the mean-value property — and iterating that averaging (Jacobi or Gauss–Seidel relaxation, [Ch. 22](../part-7-numerical-computational/22-numerical-methods.md)) solves Laplace's equation for the potential between electrodes.

## 7.6 Line integrals

### Motivation

The voltage between two points is the work per unit charge done *against* the electric field along a path: $V_{ab}=-\int_b^a\mathbf{E}\cdot d\mathbf{l}$. The EMF around a loop is $\oint\mathbf{E}\cdot d\mathbf{l}$. The total charge on a curved wire with non-uniform line density is $\int\rho_\ell\,ds$. All are integrals along a curve — **line integrals** — and they come in two kinds: integrals of a scalar with respect to arc length, and integrals of the tangential component of a vector field.

### Scalar line integrals

> [!definition] Line integral of a scalar field
> Let $C$ be a smooth curve parametrised by $\mathbf{r}(t)$, $a\le t\le b$. The line integral of $f$ along $C$ with respect to arc length is
> $$
> \int_Cf\,ds=\lim_{\max\Delta s_k\to0}\sum_kf(\mathbf{r}_k)\,\Delta s_k=\int_a^bf(\mathbf{r}(t))\,\lVert\mathbf{r}'(t)\rVert\,dt .
> $$

The value does not depend on the parametrisation *or on the direction* of travel ($ds$ is a length, always positive). With $f=1$ it gives the arc length; with $f=\rho_\ell$ the total charge; with $f$ a linear resistance per metre, the total resistance.

> [!example] Worked Example — charge on a curved wire
> A quarter-circle wire of radius $a$ in the first quadrant carries line charge density $\rho_\ell=\rho_0\,xy/a^2$. Find its total charge.
>
> **Step 1 — parametrise.** $\mathbf{r}(\phi)=(a\cos\phi,\ a\sin\phi)$, $0\le\phi\le\pi/2$; $\lVert\mathbf{r}'\rVert=a$, so $ds=a\,d\phi$.
>
> **Step 2 — integrand.** $\rho_\ell=\rho_0\cos\phi\sin\phi$.
>
> **Step 3 — integrate.** $Q=\int_0^{\pi/2}\rho_0\cos\phi\sin\phi\,a\,d\phi=\rho_0a\left[\tfrac12\sin^2\phi\right]_0^{\pi/2}=\tfrac12\rho_0a$.
>
> **Units.** $\rho_0$ in C/m times $a$ in m gives coulombs ✓.

### Line integrals of vector fields

> [!definition] Line integral of a vector field (work, circulation)
> For a vector field $\mathbf{F}$ and an **oriented** curve $C$ parametrised by $\mathbf{r}(t)$, $a\le t\le b$ (direction of increasing $t$),
> $$
> \int_C\mathbf{F}\cdot d\mathbf{r}=\int_C\mathbf{F}\cdot\hat{\mathbf{T}}\,ds=\int_a^b\mathbf{F}(\mathbf{r}(t))\cdot\mathbf{r}'(t)\,dt=\int_C(F_x\,dx+F_y\,dy+F_z\,dz).
> $$
> If $C$ is closed the integral is called the **circulation** and written $\oint_C\mathbf{F}\cdot d\mathbf{r}$.

Only the tangential component $\mathbf{F}\cdot\hat{\mathbf{T}}$ contributes. Reversing the orientation of $C$ reverses $\hat{\mathbf{T}}$ and so **changes the sign**: $\int_{-C}\mathbf{F}\cdot d\mathbf{r}=-\int_C\mathbf{F}\cdot d\mathbf{r}$. Any other parametrisation with the same orientation gives the same value: if $t=g(u)$ with $g'>0$, then $\mathbf{r}'(g(u))g'(u)\,du=\mathbf{r}'(t)\,dt$. For a curve made of smooth pieces (a polygon, say) add the integrals over the pieces.

**The procedure** — almost every line-integral computation follows these four steps:

1. Parametrise each piece of $C$ as $\mathbf{r}(t)$, with $t$ increasing in the direction of travel.
2. Compute $d\mathbf{r}=\mathbf{r}'(t)\,dt$.
3. Substitute $\mathbf{r}(t)$ into $\mathbf{F}$ and form the scalar $\mathbf{F}(\mathbf{r}(t))\cdot\mathbf{r}'(t)$.
4. Integrate in $t$ from the start value to the end value.

> [!example] Worked Example — path dependence
> Compute $\int_C\mathbf{F}\cdot d\mathbf{r}$ from $(0,0)$ to $(1,1)$ for (a) $\mathbf{F}=y\hat{\mathbf{x}}+x\hat{\mathbf{y}}$ and (b) $\mathbf{G}=-y\hat{\mathbf{x}}+x\hat{\mathbf{y}}$, along $C_1$: the line $y=x$, and $C_2$: the parabola $y=x^2$.
>
> **$C_1$:** $\mathbf{r}=(t,t)$, $\mathbf{r}'=(1,1)$, $0\le t\le1$.
> (a) $\mathbf{F}\cdot\mathbf{r}'=t+t=2t$, integral $1$.
> (b) $\mathbf{G}\cdot\mathbf{r}'=-t+t=0$, integral $0$.
>
> **$C_2$:** $\mathbf{r}=(t,t^2)$, $\mathbf{r}'=(1,2t)$.
> (a) $\mathbf{F}=(t^2,t)$, $\mathbf{F}\cdot\mathbf{r}'=t^2+2t^2=3t^2$, integral $1$.
> (b) $\mathbf{G}=(-t^2,t)$, $\mathbf{G}\cdot\mathbf{r}'=-t^2+2t^2=t^2$, integral $\tfrac13$.
>
> $\mathbf{F}$ gives the same answer on both paths; $\mathbf{G}$ does not. The reason: $\mathbf{F}=\nabla(xy)$ is a gradient, while $\nabla\times\mathbf{G}=2\hat{\mathbf{z}}\ne\mathbf{0}$. Going out along $C_1$ and back along $C_2$ gives a circulation $0-\tfrac13=-\tfrac13$ for $\mathbf{G}$: by Green's theorem (Section 7.9) this is $-2\times$ (area between the curves $=\tfrac16$), the minus sign because that loop runs clockwise.

```desmos-graph
left=-0.2; right=1.3; bottom=-0.2; top=1.3
---
y=x|blue|0<x<1
y=x^2|red|0<x<1
(0,0)|label:start
(1,1)|label:end
```

*Two paths from $(0,0)$ to $(1,1)$: the line $C_1$ (blue) and the parabola $C_2$ (red). A gradient field gives the same work along both; the rotation field $\mathbf{G}$ does not.*

### The fundamental theorem for line integrals and conservative fields

> [!theorem] Fundamental theorem for line integrals
> If $f$ is continuously differentiable on a region containing a smooth curve $C$ from $\mathbf{a}$ to $\mathbf{b}$, then
> $$
> \int_C\nabla f\cdot d\mathbf{r}=f(\mathbf{b})-f(\mathbf{a}).
> $$

*Proof.* By the chain rule, $\frac{d}{dt}f(\mathbf{r}(t))=\nabla f\cdot\mathbf{r}'(t)$, so the integral is $\int_a^b\frac{d}{dt}f(\mathbf{r}(t))\,dt=f(\mathbf{r}(b))-f(\mathbf{r}(a))$ by the ordinary fundamental theorem of calculus. ∎

> [!theorem] Conservative fields — equivalent statements
> For a continuous field $\mathbf{F}$ on an open connected region $D$, the following are equivalent:
> 1. $\mathbf{F}=\nabla f$ for some scalar $f$ on $D$;
> 2. $\int_C\mathbf{F}\cdot d\mathbf{r}$ depends only on the endpoints of $C$ (path independence);
> 3. $\oint_C\mathbf{F}\cdot d\mathbf{r}=0$ for every closed curve $C$ in $D$.
>
> Each implies $\nabla\times\mathbf{F}=\mathbf{0}$; conversely $\nabla\times\mathbf{F}=\mathbf{0}$ implies 1–3 **if $D$ is simply connected**.

*Sketch.* 1⇒2 is the fundamental theorem. 2⇔3: two paths with the same ends form a closed loop, and the difference of their integrals is the circulation. 2⇒1: fix $\mathbf{a}$ and define $f(\mathbf{x})=\int_{\mathbf{a}}^{\mathbf{x}}\mathbf{F}\cdot d\mathbf{r}$ (well defined by path independence); moving the endpoint a little in the $x$ direction along a straight segment shows $\partial f/\partial x=F_x$, and likewise for $y,z$. The curl statement is Stokes' theorem (Section 7.10).

**Finding a potential.** Given an irrotational $\mathbf{F}$ on a simply connected region, integrate one component at a time:

1. $f=\int F_x\,dx+g(y,z)$, where $g$ is an unknown "constant" of integration;
2. differentiate with respect to $y$, set equal to $F_y$, and solve for $\partial g/\partial y$;
3. repeat for $z$; the leftover constant is arbitrary (only potential *differences* are physical).

> [!example] Worked Example — finding and using a potential
> Show that $\mathbf{F}=(2xy+z)\hat{\mathbf{x}}+x^2\hat{\mathbf{y}}+(x+3z^2)\hat{\mathbf{z}}$ is conservative, find $f$ with $\mathbf{F}=\nabla f$, and evaluate the work from $(0,0,0)$ to $(1,2,1)$.
>
> **Step 1 — curl test.**
> $\partial_yF_z-\partial_zF_y=0-0$, $\partial_zF_x-\partial_xF_z=1-1$, $\partial_xF_y-\partial_yF_x=2x-2x$. All zero, and the domain ($\mathbb{R}^3$) is simply connected, so $\mathbf{F}$ is conservative.
>
> **Step 2 — integrate $F_x$.** $f=\int(2xy+z)\,dx=x^2y+xz+g(y,z)$.
>
> **Step 3 — match $F_y$.** $\partial_yf=x^2+\partial_yg=x^2$, so $g=h(z)$.
>
> **Step 4 — match $F_z$.** $\partial_zf=x+h'(z)=x+3z^2$, so $h=z^3+c$. Thus $f=x^2y+xz+z^3$.
>
> **Step 5 — work.** $f(1,2,1)-f(0,0,0)=2+1+1=4$.
>
> **Check along the straight line** $\mathbf{r}=(t,2t,t)$: $\mathbf{F}=(4t^2+t,\ t^2,\ t+3t^2)$, $\mathbf{r}'=(1,2,1)$, $\mathbf{F}\cdot\mathbf{r}'=9t^2+2t$, $\int_0^1(9t^2+2t)\,dt=3+1=4$. ✓

> [!example] Worked Example — circulation and the non-simply-connected trap
> Compute the circulation counter-clockwise around the circle $x^2+y^2=a^2$ of (a) the rotation field $\mathbf{G}=-y\hat{\mathbf{x}}+x\hat{\mathbf{y}}$ and (b) the magnetic field of a long wire on the $z$-axis, $\mathbf{B}=\dfrac{\mu_0I}{2\pi r}\hat{\boldsymbol{\phi}}$.
>
> **Parametrise.** $\mathbf{r}=(a\cos\phi,a\sin\phi,0)$, $\mathbf{r}'=a(-\sin\phi,\cos\phi,0)=a\hat{\boldsymbol{\phi}}$, so $d\mathbf{r}=a\,d\phi\,\hat{\boldsymbol{\phi}}$.
>
> **(a)** $\mathbf{G}=a\hat{\boldsymbol{\phi}}$ on the circle, so $\oint\mathbf{G}\cdot d\mathbf{r}=\int_0^{2\pi}a\cdot a\,d\phi=2\pi a^2$. This equals $(\nabla\times\mathbf{G})_z\times\text{area}=2\times\pi a^2$, as the curl definition predicts.
>
> **(b)** $\oint\mathbf{B}\cdot d\mathbf{r}=\int_0^{2\pi}\dfrac{\mu_0I}{2\pi a}\,a\,d\phi=\mu_0I$ — **independent of the radius**. This is Ampère's circuital law.
>
> So $\mathbf{B}$ has zero curl wherever it is defined ($r>0$, Section 7.5), yet its circulation around the wire is not zero. There is no contradiction: the region $r>0$ has a hole (the wire), loops around it cannot be shrunk to a point, and "curl-free ⇒ conservative" fails. A loop that does **not** encircle the wire has zero circulation. The "potential" $\mu_0I\phi/(2\pi)$ exists only locally; it increases by $\mu_0I$ on every trip around the wire, like a multi-storey car park ramp.

> [!example] Worked Example — voltage as a line integral
> A point charge $q$ at the origin produces $\mathbf{E}=\dfrac{q}{4\pi\varepsilon_0R^2}\hat{\mathbf{R}}$. Find the voltage $V_{ab}=V_a-V_b=-\int_b^a\mathbf{E}\cdot d\mathbf{l}$ between $a$ at radius $R_a$ and $b$ at radius $R_b$, along any path.
>
> **Step 1.** In spherical coordinates $d\mathbf{l}=dR\,\hat{\mathbf{R}}+R\,d\theta\,\hat{\boldsymbol{\theta}}+R\sin\theta\,d\phi\,\hat{\boldsymbol{\phi}}$, so $\mathbf{E}\cdot d\mathbf{l}=\dfrac{q}{4\pi\varepsilon_0}\dfrac{dR}{R^2}$ — only the radial motion counts.
>
> **Step 2.** $V_{ab}=-\dfrac{q}{4\pi\varepsilon_0}\int_{R_b}^{R_a}\dfrac{dR}{R^2}=\dfrac{q}{4\pi\varepsilon_0}\left(\dfrac1{R_a}-\dfrac1{R_b}\right)$.
>
> **Step 3 — numbers.** $q=1$ nC, $R_a=1$ cm, $R_b=10$ cm: $V_{ab}=8.988\times10^{9}\times10^{-9}\times(100-10)=809$ V.
>
> Letting $R_b\to\infty$ gives the familiar $V=q/(4\pi\varepsilon_0R)$. The path never mattered: $\mathbf{E}=-\nabla V$ is a gradient.

> [!warning] Common Pitfall
> Writing $d\mathbf{r}=-dR\,\hat{\mathbf{R}}$ "because the path goes inward" *and* also integrating from the larger to the smaller radius. The parametrisation and the limits together encode the direction — put it in once. If you integrate from $R_b$ to $R_a$, use $d\mathbf{l}=dR\,\hat{\mathbf{R}}$ and let the limits do the work.

> [!tip] Engineering Note
> Kirchhoff's voltage law, $\sum v_k=0$ around a loop, is the statement $\oint\mathbf{E}\cdot d\mathbf{l}=0$. It holds when the loop links no changing magnetic flux. When it does (a loop around a transformer core, or a large ground loop near a switching converter), $\oint\mathbf{E}\cdot d\mathbf{l}=-d\Phi/dt\ne0$ and "voltage between two points" becomes path-dependent: two voltmeters with differently routed leads read different values. This is the Faraday-law effect behind much measurement interference (Section 7.12).

## 7.7 Surface integrals

### Parametrised surfaces

A surface in space is described by two parameters, just as a curve is described by one:

$$
\mathbf{r}(u,v)=x(u,v)\hat{\mathbf{x}}+y(u,v)\hat{\mathbf{y}}+z(u,v)\hat{\mathbf{z}},\qquad(u,v)\in D .
$$

Holding $v$ fixed and varying $u$ traces a curve on the surface with tangent $\mathbf{r}_u=\partial\mathbf{r}/\partial u$; similarly $\mathbf{r}_v$. A small parameter rectangle $du\times dv$ maps to a small parallelogram with edges $\mathbf{r}_u\,du$ and $\mathbf{r}_v\,dv$, whose **vector area** is

$$
d\mathbf{S}=(\mathbf{r}_u\times\mathbf{r}_v)\,du\,dv,\qquad dS=\lVert\mathbf{r}_u\times\mathbf{r}_v\rVert\,du\,dv .
$$

The vector $\mathbf{N}=\mathbf{r}_u\times\mathbf{r}_v$ is normal to the surface; the surface is smooth where $\mathbf{N}\ne\mathbf{0}$. The standard surfaces of field theory:

| Surface | Parameters | $\mathbf{r}_u\times\mathbf{r}_v$ | $dS$ |
|---|---|---|---|
| graph $z=g(x,y)$ | $x,y$ | $(-g_x,\ -g_y,\ 1)$ (upward) | $\sqrt{1+g_x^2+g_y^2}\,dx\,dy$ |
| plane $z=\text{const}$ | $x,y$ | $\hat{\mathbf{z}}$ | $dx\,dy$ |
| cylinder $r=a$ | $\phi,z$ | $a\,\hat{\mathbf{r}}$ (outward) | $a\,d\phi\,dz$ |
| sphere $R=a$ | $\theta,\phi$ | $a^2\sin\theta\,\hat{\mathbf{R}}$ (outward) | $a^2\sin\theta\,d\theta\,d\phi$ |
| disk/annulus $z=\text{const}$ | $r,\phi$ | $r\,\hat{\mathbf{z}}$ | $r\,dr\,d\phi$ |

*Derivation for the sphere.* $\mathbf{r}=a(\sin\theta\cos\phi,\ \sin\theta\sin\phi,\ \cos\theta)$. Then $\mathbf{r}_\theta=a\hat{\boldsymbol{\theta}}$ and $\mathbf{r}_\phi=a\sin\theta\,\hat{\boldsymbol{\phi}}$ (lengths $h_\theta=a$, $h_\phi=a\sin\theta$), and $\hat{\boldsymbol{\theta}}\times\hat{\boldsymbol{\phi}}=\hat{\mathbf{R}}$, so $\mathbf{r}_\theta\times\mathbf{r}_\phi=a^2\sin\theta\,\hat{\mathbf{R}}$.

### Scalar surface integrals and flux

> [!definition] Surface integrals
> For a smooth surface $S$ parametrised by $\mathbf{r}(u,v)$ over $D$:
> - **scalar surface integral:** $\displaystyle\iint_Sf\,dS=\iint_Df(\mathbf{r}(u,v))\,\lVert\mathbf{r}_u\times\mathbf{r}_v\rVert\,du\,dv$;
> - **flux** of $\mathbf{F}$ through $S$ with unit normal $\hat{\mathbf{n}}$: $\displaystyle\iint_S\mathbf{F}\cdot d\mathbf{S}=\iint_S\mathbf{F}\cdot\hat{\mathbf{n}}\,dS=\iint_D\mathbf{F}(\mathbf{r}(u,v))\cdot(\mathbf{r}_u\times\mathbf{r}_v)\,du\,dv$,
>
> the last form being valid when $\mathbf{r}_u\times\mathbf{r}_v$ points along the chosen $\hat{\mathbf{n}}$ (otherwise change the sign). For a closed surface the flux is written $\oint_S\mathbf{F}\cdot d\mathbf{S}$ and $\hat{\mathbf{n}}$ is taken **outward**.

Flux counts how much of the field passes *through* the surface: only the normal component $\mathbf{F}\cdot\hat{\mathbf{n}}$ matters; a field tangent to the surface contributes nothing. If $\mathbf{F}$ is the velocity of a fluid, the flux is the volume flow rate ($\mathrm{m^3/s}$); if it is the current density $\mathbf{J}$, the flux is the current $I=\iint\mathbf{J}\cdot d\mathbf{S}$; if it is $\mathbf{B}$, the flux $\Phi$ is in webers.

**Orientation.** A flux needs a choice of side. An **orientable** surface admits a continuous choice of unit normal; a closed surface is oriented outward by convention; an open surface bounded by a curve is oriented consistently with the curve's direction by the right-hand rule (Section 7.10). The Möbius strip is the standard example of a non-orientable surface — you cannot define the flux through it.

> [!example] Worked Example — areas: spherical cap and paraboloid
> **(a)** Find the area of the cap $0\le\theta\le\alpha$ of a sphere of radius $a$.
> $$
> A=\int_0^{2\pi}\int_0^\alpha a^2\sin\theta\,d\theta\,d\phi=2\pi a^2(1-\cos\alpha).
> $$
> For $\alpha=\pi$ this is $4\pi a^2$ ✓; for $\alpha=\pi/2$, the hemisphere $2\pi a^2$ ✓. The cap's **solid angle** is $A/a^2=2\pi(1-\cos\alpha)$ steradians, the quantity used for antenna beam widths.
>
> **(b)** Find the area of the paraboloid $z=x^2+y^2$ below $z=1$.
> Here $g_x=2x$, $g_y=2y$, so $dS=\sqrt{1+4x^2+4y^2}\,dx\,dy=\sqrt{1+4r^2}\,r\,dr\,d\phi$ over the unit disk:
> $$
> A=2\pi\int_0^1\sqrt{1+4r^2}\,r\,dr=2\pi\left[\frac{(1+4r^2)^{3/2}}{12}\right]_0^1=\frac{\pi}{6}\left(5\sqrt5-1\right)=5.330 .
> $$
> Sanity check: larger than the unit disk ($\pi=3.14$) it sits above, as a curved surface must be.

> [!example] Worked Example — outward flux through a sphere
> Find the outward flux of $\mathbf{F}=x\hat{\mathbf{x}}+y\hat{\mathbf{y}}+z\hat{\mathbf{z}}=\mathbf{R}$ through the sphere $R=a$.
>
> On the sphere $\mathbf{F}=a\hat{\mathbf{R}}$ and $d\mathbf{S}=a^2\sin\theta\,d\theta\,d\phi\,\hat{\mathbf{R}}$, so $\mathbf{F}\cdot d\mathbf{S}=a^3\sin\theta\,d\theta\,d\phi$ and
> $$
> \oint_S\mathbf{F}\cdot d\mathbf{S}=a^3\int_0^{2\pi}\int_0^\pi\sin\theta\,d\theta\,d\phi=4\pi a^3 .
> $$
> Since $\nabla\cdot\mathbf{F}=3$ and the ball has volume $\tfrac43\pi a^3$, the flux is $3\times$ volume — a preview of the divergence theorem.
>
> By contrast, for the inverse-square field $\hat{\mathbf{R}}/R^2$ the flux is $\frac{1}{a^2}\cdot4\pi a^2=4\pi$, **independent of $a$**. Every sphere around the origin catches all the field lines.

> [!example] Worked Example — flux through a curved open surface
> Find the upward flux of $\mathbf{F}=z\,\hat{\mathbf{z}}$ through the paraboloid cap $S$: $z=1-x^2-y^2$, $z\ge0$.
>
> **Step 1 — normal.** With $g=1-x^2-y^2$, $\mathbf{N}=(-g_x,-g_y,1)=(2x,2y,1)$, which points upward ✓.
>
> **Step 2 — integrand.** $\mathbf{F}\cdot\mathbf{N}=z\cdot1=1-x^2-y^2=1-r^2$.
>
> **Step 3 — integrate over the projection** (the unit disk):
> $$
> \iint_S\mathbf{F}\cdot d\mathbf{S}=\int_0^{2\pi}\int_0^1(1-r^2)\,r\,dr\,d\phi=2\pi\left(\tfrac12-\tfrac14\right)=\frac{\pi}{2}.
> $$
> We will recover this from the divergence theorem in Section 7.11.

> [!example] Worked Example — flux of a point charge through one face of a cube
> A charge $q$ sits at the centre of a cube of side $2d$. Find the flux of $\mathbf{E}$ through the top face, both by symmetry and directly.
>
> **By symmetry.** The total outward flux through any closed surface around $q$ is $q/\varepsilon_0$ (Gauss's law, Section 7.12); the six faces are equivalent, so each receives $q/(6\varepsilon_0)$.
>
> **Directly.** On the top face $z=d$, $\hat{\mathbf{n}}=\hat{\mathbf{z}}$, and $\mathbf{E}\cdot\hat{\mathbf{z}}=\dfrac{q}{4\pi\varepsilon_0}\dfrac{d}{(x^2+y^2+d^2)^{3/2}}$. So
> $$
> \Phi_E=\frac{q}{4\pi\varepsilon_0}\int_{-d}^{d}\int_{-d}^{d}\frac{d\,dx\,dy}{(x^2+y^2+d^2)^{3/2}}=\frac{q}{4\pi\varepsilon_0}\cdot\Omega ,
> $$
> where $\Omega$ is the solid angle the face subtends. Substituting $x=du$, $y=dv$ removes $d$, and the standard result for a rectangle $\int_{-a}^{a}\int_{-b}^{b}\frac{du\,dv}{(u^2+v^2+1)^{3/2}}=4\arctan\dfrac{ab}{\sqrt{a^2+b^2+1}}$ with $a=b=1$ gives
> $$
> \Omega=4\arctan\frac{1}{\sqrt3}=4\cdot\frac{\pi}{6}=\frac{2\pi}{3}.
> $$
> Hence $\Phi_E=\dfrac{q}{4\pi\varepsilon_0}\cdot\dfrac{2\pi}{3}=\dfrac{q}{6\varepsilon_0}$ ✓ (and $2\pi/3=2.0944$ was confirmed by numerical quadrature).

> [!example] Worked Example — current through a conductor cross-section
> A round wire of radius $a$ carries current density $\mathbf{J}=J_0(1-r^2/a^2)\hat{\mathbf{z}}$ (peaked on the axis). Find the total current, and the radius inside which half of it flows.
>
> **Total.** $I=\iint\mathbf{J}\cdot\hat{\mathbf{z}}\,dS=\int_0^{2\pi}\int_0^aJ_0\left(1-\frac{r^2}{a^2}\right)r\,dr\,d\phi=2\pi J_0\left(\frac{a^2}{2}-\frac{a^2}{4}\right)=\frac{\pi a^2J_0}{2}$ — half of what a uniform $J_0$ would give.
>
> **Half-current radius.** $I(\rho)=2\pi J_0\left(\frac{\rho^2}{2}-\frac{\rho^4}{4a^2}\right)=\frac{\pi a^2J_0}{4}$ gives, with $s=\rho^2/a^2$, $s-\tfrac{s^2}{2}=\tfrac14$, i.e. $s^2-2s+\tfrac12=0$, $s=1-1/\sqrt2=0.2929$, so $\rho=0.541a$.

> [!warning] Common Pitfall
> Forgetting the area factor. On a sphere $dS=a^2\sin\theta\,d\theta\,d\phi$, not $d\theta\,d\phi$; on a cylinder $dS=a\,d\phi\,dz$. The quickest check is to integrate $1$ and see whether you get the known area. A second pitfall is projecting onto the $xy$-plane without the factor $\sqrt{1+g_x^2+g_y^2}$ (for scalar integrals) — for *flux* the vector $(-g_x,-g_y,1)$ already contains the right factor.

> [!tip] Engineering Note
> Radiated power is a flux: an antenna's total radiated power is $P=\oint\langle\mathbf{S}\rangle\cdot d\mathbf{S}$ over a large sphere, and its **directivity** compares the power per unit solid angle in the main beam with the average $P/(4\pi)$. The $\sin\theta$ in $dS$ is why a beam near the pole of the coordinate system carries less power than one of the same angular width near the equator.

## 7.8 Volume integrals

### Definition and the three volume elements

Volume integrals were set up in [Ch. 6](../part-2-calculus/06-multivariable-calculus.md) as triple integrals. In field theory they answer questions such as: how much charge is in this region, how much energy is stored in this field, where is the centre of this charge distribution?

> [!definition] Volume integral
> For a scalar field $f$ over a region $V$,
> $$
> \iiint_Vf\,dV=\lim_{\max\Delta V_k\to0}\sum_kf(\mathbf{R}_k)\,\Delta V_k ,
> $$
> evaluated as an iterated integral with
> $$
> dV=dx\,dy\,dz=r\,dr\,d\phi\,dz=R^2\sin\theta\,dR\,d\theta\,d\phi .
> $$
> A **vector** volume integral $\iiint_V\mathbf{F}\,dV$ is computed component by component **in Cartesian components**.

The volume elements are $h_1h_2h_3\,du_1\,du_2\,du_3$: the product of the three edge lengths of an infinitesimal coordinate box. The choice of coordinates is dictated by the shape of the region more than by the integrand: spheres and cones in spherical coordinates, cylinders and pipes in cylindrical, boxes in Cartesian.

| Quantity | Volume integral |
|---|---|
| total charge | $Q=\iiint\rho\,dV$ |
| mass, from density $\rho_m$ | $M=\iiint\rho_m\,dV$ |
| electrostatic energy | $W_e=\tfrac12\iiint\varepsilon\lVert\mathbf{E}\rVert^2dV$ |
| magnetic energy | $W_m=\tfrac12\iiint\mu\lVert\mathbf{H}\rVert^2dV$ |
| ohmic power | $P=\iiint\mathbf{J}\cdot\mathbf{E}\,dV=\iiint\sigma\lVert\mathbf{E}\rVert^2dV$ |
| dipole moment | $\mathbf{p}=\iiint\mathbf{R}\,\rho\,dV$ |
| potential of a distribution | $V(\mathbf{R})=\dfrac1{4\pi\varepsilon_0}\iiint\dfrac{\rho(\mathbf{R}')}{\lVert\mathbf{R}-\mathbf{R}'\rVert}dV'$ |

> [!example] Worked Example — charge in a non-uniform sphere
> A sphere of radius $a$ has charge density $\rho=\rho_0(1-R/a)$. Find the total charge.
>
> **Step 1.** The integrand depends only on $R$, so the angular integrals give $\int_0^{2\pi}\int_0^\pi\sin\theta\,d\theta\,d\phi=4\pi$.
>
> **Step 2.** $Q=4\pi\rho_0\int_0^a\left(1-\dfrac Ra\right)R^2\,dR=4\pi\rho_0\left(\dfrac{a^3}{3}-\dfrac{a^3}{4}\right)=\dfrac{\pi\rho_0a^3}{3}$.
>
> **Check.** A uniform density $\rho_0$ would give $\tfrac43\pi\rho_0a^3$; the tapered sphere holds a quarter of that. Plausible: most of a ball's volume is near its surface, where this density is small.

> [!example] Worked Example — charge in a cylinder (separable integrand)
> A cylinder $0\le r\le a$, $0\le z\le L$ holds $\rho=\rho_0\,\dfrac ra\,e^{-z/L}$. Find $Q$.
>
> The integrand is a product of functions of $r$, $\phi$ and $z$ separately and the limits are constants, so the triple integral is a product of three single integrals:
> $$
> Q=\rho_0\underbrace{\int_0^{2\pi}d\phi}_{2\pi}\ \underbrace{\int_0^a\frac ra\,r\,dr}_{a^2/3}\ \underbrace{\int_0^Le^{-z/L}dz}_{L(1-e^{-1})}=\frac{2\pi}{3}\left(1-e^{-1}\right)\rho_0a^2L=1.324\,\rho_0a^2L .
> $$
> Do not forget the $r$ in $dV$: it turns $\int r\,dr$ into $\int r^2\,dr$ here.

> [!example] Worked Example — energy stored in the field of a charged sphere
> A sphere of radius $a$ carries total charge $Q$ uniformly through its volume. Its field (derived in Section 7.12) is $E=\dfrac{QR}{4\pi\varepsilon_0a^3}$ inside and $E=\dfrac{Q}{4\pi\varepsilon_0R^2}$ outside. Find the stored energy $W=\tfrac{\varepsilon_0}{2}\iiint E^2\,dV$ over all space.
>
> **Inside.**
> $$
> W_{\text{in}}=\frac{\varepsilon_0}{2}\left(\frac{Q}{4\pi\varepsilon_0a^3}\right)^2\int_0^aR^2\cdot4\pi R^2\,dR=\frac{\varepsilon_0}{2}\cdot\frac{Q^2}{16\pi^2\varepsilon_0^2a^6}\cdot\frac{4\pi a^5}{5}=\frac{Q^2}{40\pi\varepsilon_0a}.
> $$
> **Outside.**
> $$
> W_{\text{out}}=\frac{\varepsilon_0}{2}\cdot\frac{Q^2}{16\pi^2\varepsilon_0^2}\int_a^\infty\frac{4\pi R^2}{R^4}\,dR=\frac{Q^2}{8\pi\varepsilon_0a}.
> $$
> **Total.** $W=\dfrac{Q^2}{8\pi\varepsilon_0a}\left(\dfrac15+1\right)=\dfrac{3Q^2}{20\pi\varepsilon_0a}$.
>
> **Numbers.** $Q=1$ nC, $a=1$ cm: $W=3\times10^{-18}/(20\pi\times8.854\times10^{-12}\times0.01)=5.39\times10^{-7}$ J. Five-sixths of the energy is stored *outside* the charge — the field, not the charge, holds the energy.

> [!example] Worked Example — a vector volume integral: dipole moment
> A sphere of radius $a$ has $\rho=\rho_0\cos\theta$ (positive in the upper half, negative in the lower). Find the total charge and the dipole moment $\mathbf{p}=\iiint\mathbf{R}\rho\,dV$.
>
> **Charge.** $Q=\rho_0\int_0^a R^2dR\int_0^\pi\cos\theta\sin\theta\,d\theta\int_0^{2\pi}d\phi=0$, since $\int_0^\pi\cos\theta\sin\theta\,d\theta=0$.
>
> **Dipole moment.** $\mathbf{R}=R\hat{\mathbf{R}}$, but $\hat{\mathbf{R}}$ varies, so use Cartesian components: $\mathbf{R}=(R\sin\theta\cos\phi,\ R\sin\theta\sin\phi,\ R\cos\theta)$. The $x$- and $y$-components contain $\int_0^{2\pi}\cos\phi\,d\phi=0$ and $\int_0^{2\pi}\sin\phi\,d\phi=0$. The $z$-component is
> $$
> p_z=\rho_0\int_0^aR^3\,dR\int_0^\pi\cos^2\theta\sin\theta\,d\theta\int_0^{2\pi}d\phi=\rho_0\cdot\frac{a^4}{4}\cdot\frac23\cdot2\pi=\frac{\pi\rho_0a^4}{3}.
> $$
> So $\mathbf{p}=\tfrac13\pi\rho_0a^4\,\hat{\mathbf{z}}$, pointing from the negative half to the positive half, as expected. Far away, this neutral sphere looks like the dipole of Section 7.4.

> [!warning] Common Pitfall
> Integrating the spherical *components* of a vector field over a volume, e.g. writing $\iiint\hat{\mathbf{R}}\,dV=\hat{\mathbf{R}}\iiint dV$. $\hat{\mathbf{R}}$ is not a constant and cannot leave the integral; over a full sphere $\iiint\hat{\mathbf{R}}\,dV=\mathbf{0}$ by symmetry. Convert to Cartesian components first.

> [!tip] Engineering Note
> Capacitance and inductance can be computed from energy: $C=2W_e/V^2$ and $L=2W_m/I^2$. Field solvers often extract $L$ and $C$ of interconnects this way, by integrating $\tfrac12\varepsilon E^2$ and $\tfrac12\mu H^2$ over the mesh — a large, numerical volume integral ([Ch. 24](../part-7-numerical-computational/24-computational-mathematics.md)).

## 7.9 Green's theorem

### Statement

Green's theorem is the two-dimensional prototype of all the integral theorems. It relates a line integral around a closed plane curve to a double integral over the region inside.

> [!theorem] Green's theorem
> Let $D$ be a bounded region of the plane whose boundary $C$ consists of finitely many piecewise-smooth simple closed curves, oriented **positively** (as you walk along $C$, $D$ is on your left — counter-clockwise for an outer boundary, clockwise for the boundary of a hole). If $P(x,y)$ and $Q(x,y)$ have continuous first partial derivatives on an open set containing $D$, then
> $$
> \oint_C(P\,dx+Q\,dy)=\iint_D\left(\frac{\partial Q}{\partial x}-\frac{\partial P}{\partial y}\right)dA .
> $$

With $\mathbf{F}=P\hat{\mathbf{x}}+Q\hat{\mathbf{y}}$, the left side is the circulation $\oint\mathbf{F}\cdot d\mathbf{r}$ and the integrand on the right is $(\nabla\times\mathbf{F})\cdot\hat{\mathbf{z}}$. In words: **the circulation around the boundary equals the total curl inside.**

### Proof sketch

*Step 1 — the $P$ term on a "vertically simple" region.* Suppose $D=\{a\le x\le b,\ g_1(x)\le y\le g_2(x)\}$. Then

$$
-\iint_D\frac{\partial P}{\partial y}\,dA=-\int_a^b\Big[P(x,g_2(x))-P(x,g_1(x))\Big]dx .
$$

Now compute $\oint_CP\,dx$ directly. On the bottom curve $y=g_1(x)$, traversed left to right, it gives $\int_a^bP(x,g_1(x))\,dx$. On the top curve, traversed right to left, it gives $-\int_a^bP(x,g_2(x))\,dx$. On any vertical sides $dx=0$. The two expressions agree.

*Step 2 — the $Q$ term* on a "horizontally simple" region is identical with the roles of $x$ and $y$ exchanged (and a sign flip that produces $+\partial Q/\partial x$).

*Step 3 — general regions.* Cut $D$ into pieces that are both vertically and horizontally simple. Apply the theorem to each piece and add. The integrals along each internal cut appear twice, in opposite directions, and cancel, leaving only the outer boundary. ∎

The cancellation in Step 3 is the idea behind every integral theorem: **interior contributions cancel in pairs and only the boundary survives.** It is also how the coordinate-free definition of curl (circulation per unit area of a tiny loop) adds up to the theorem: tile $D$ with tiny squares, each contributing $(\nabla\times\mathbf{F})_z\,\Delta A$, and the internal edges cancel.

### The normal (flux) form

Apply Green's theorem to $P=-F_y$, $Q=F_x$. The outward unit normal to a counter-clockwise curve is $\hat{\mathbf{n}}=(dy,-dx)/ds$, so $\mathbf{F}\cdot\hat{\mathbf{n}}\,ds=F_x\,dy-F_y\,dx$. Hence

$$
\oint_C\mathbf{F}\cdot\hat{\mathbf{n}}\,ds=\iint_D\left(\frac{\partial F_x}{\partial x}+\frac{\partial F_y}{\partial y}\right)dA=\iint_D\nabla\cdot\mathbf{F}\,dA :
$$

**the outward flux through the boundary equals the total divergence inside** — the two-dimensional divergence theorem.

### Area by a line integral

Choosing $P,Q$ with $\partial_xQ-\partial_yP=1$ turns the right side into the area of $D$. Three convenient choices:

$$
A=\oint_Cx\,dy=-\oint_Cy\,dx=\frac12\oint_C(x\,dy-y\,dx).
$$

For a polygon with vertices $(x_1,y_1),\dots,(x_n,y_n)$ listed counter-clockwise, applying the last formula edge by edge gives the **shoelace formula**

$$
A=\frac12\sum_{k=1}^{n}\left(x_ky_{k+1}-x_{k+1}y_k\right),\qquad(x_{n+1},y_{n+1})=(x_1,y_1).
$$

This is how CAD and PCB-layout tools compute the area of an arbitrary copper pour, and how a mechanical planimeter measures area by tracing the boundary.

> [!example] Worked Example — verifying Green's theorem on a square
> Let $P=x-y$, $Q=x^2$ on the unit square $[0,1]\times[0,1]$.
>
> **Double integral.** $\partial_xQ-\partial_yP=2x-(-1)=2x+1$, and $\int_0^1\int_0^1(2x+1)\,dx\,dy=1+1=2$.
>
> **Line integral**, counter-clockwise from the origin:
> - bottom ($y=0$, $x:0\to1$, $dy=0$): $\int_0^1x\,dx=\tfrac12$;
> - right ($x=1$, $y:0\to1$, $dx=0$): $\int_0^1Q\,dy=\int_0^11\,dy=1$;
> - top ($y=1$, $x:1\to0$, $dy=0$): $\int_1^0(x-1)\,dx=\tfrac12$;
> - left ($x=0$, $y:1\to0$, $dx=0$): $Q=0$, contribution $0$.
>
> Total $\tfrac12+1+\tfrac12+0=2$ ✓.

> [!example] Worked Example — area of an ellipse and of a polygon
> **Ellipse.** $x=a\cos t$, $y=b\sin t$, $0\le t\le2\pi$ (counter-clockwise). Then $x\,dy-y\,dx=(ab\cos^2t+ab\sin^2t)\,dt=ab\,dt$, so $A=\tfrac12\int_0^{2\pi}ab\,dt=\pi ab$.
>
> **Polygon.** A copper pad has vertices $(0,0),(4,0),(4,3),(1,5)$ mm, in counter-clockwise order. Shoelace terms $x_ky_{k+1}-x_{k+1}y_k$:
> $(0)(0)-(4)(0)=0$; $(4)(3)-(4)(0)=12$; $(4)(5)-(1)(3)=17$; $(1)(0)-(0)(5)=0$.
> $A=\tfrac12(0+12+17+0)=14.5\ \mathrm{mm}^2$.
> **Check by splitting** along the diagonal from $(0,0)$ to $(4,3)$: triangle $(0,0),(4,0),(4,3)$ has area $\tfrac12\cdot4\cdot3=6$, and triangle $(0,0),(4,3),(1,5)$ has area $\tfrac12\lvert4\cdot5-3\cdot1\rvert=8.5$ (half the cross product of its edge vectors, Section 7.2); total $14.5$ ✓. Listing the vertices clockwise would give $-14.5$: the sign of the shoelace sum reports the orientation.

```desmos-graph
left=-1; right=5; bottom=-1; top=6
---
y=0|blue|0<x<4
x=4|blue|0<y<3
y=3-\frac{2}{3}(x-4)|blue|1<x<4
y=5x|blue|0<x<1
y=\frac{3}{4}x|grey|dashed|0<x<4
```

*The pad of the shoelace example. The dashed diagonal splits it into triangles of area 6 and 8.5, confirming $A=14.5$.*

> [!example] Worked Example — a region with a hole
> Let $\mathbf{F}=\dfrac{-y\,\hat{\mathbf{x}}+x\,\hat{\mathbf{y}}}{x^2+y^2}$ (the wire field of Section 7.5, without constants). Evaluate $\oint\mathbf{F}\cdot d\mathbf{r}$ counter-clockwise around the ellipse $C$: $x^2/4+y^2=1$.
>
> **Why the direct approach is unpleasant.** On the ellipse $x^2+y^2$ is not constant, and the resulting integral $\int_0^{2\pi}\frac{2\,dt}{4\cos^2t+\sin^2t}$ needs a substitution.
>
> **Green's theorem with a hole.** $\nabla\times\mathbf{F}=\mathbf{0}$ everywhere except the origin, where $\mathbf{F}$ is undefined, so we cannot apply Green's theorem to the inside of the ellipse. Instead cut out a small circle $C_\epsilon$ of radius $\epsilon$ and apply it to the region $D$ between $C_\epsilon$ and $C$. Its positively oriented boundary is $C$ counter-clockwise plus $C_\epsilon$ clockwise:
> $$
> \oint_C\mathbf{F}\cdot d\mathbf{r}-\oint_{C_\epsilon,\ \circlearrowleft}\mathbf{F}\cdot d\mathbf{r}=\iint_D0\,dA=0 .
> $$
> So the integral around the ellipse equals the integral around the small circle, which is $\int_0^{2\pi}\frac{\epsilon^2}{\epsilon^2}\,d\phi=2\pi$.
>
> **Numerical check.** $\int_0^{2\pi}\frac{2\,dt}{4\cos^2t+\sin^2t}=2\cdot\frac{2\pi}{\sqrt{4\cdot1}}=2\pi$ using $\int_0^{2\pi}\frac{dt}{A\cos^2t+B\sin^2t}=\frac{2\pi}{\sqrt{AB}}$ ✓.
>
> Any closed curve winding once around the origin gives $2\pi$; any curve not enclosing it gives $0$. The integral counts how many times the curve **winds** around the singularity — the same idea as the residue theorem of complex analysis ([Ch. 18](../part-5-transforms/18-complex-analysis.md)).

```desmos-graph
left=-3; right=3; bottom=-2; top=2
---
\frac{x^2}{4}+y^2=1|blue
x^2+y^2=0.09|red
(0,0)|label:singularity
```

*The ellipse (blue) and a small circle (red) around the singular point: the region between them is curl-free, so the two circulations are equal ($2\pi$).*

> [!warning] Common Pitfall
> Applying Green's theorem when $P$ or $Q$ is singular inside $C$. The hypotheses demand continuous partial derivatives on the **whole** region. Ignoring the singularity of the wire field would "prove" that $\oint\mathbf{F}\cdot d\mathbf{r}=0$, contradicting the direct calculation.

> [!tip] Engineering Note
> In two-dimensional field problems (long transmission lines, the cross-section of a motor), Green's theorem is the tool that turns a circulation around a conductor's cross-section into the current it carries, and a flux through a closed contour into the charge per unit length inside. Boundary-element methods for 2-D electrostatics are built on its identities.

## 7.10 Stokes' theorem

### Statement

Green's theorem lives in a plane. Stokes' theorem lifts it to any curved surface in space: the circulation around the rim of a surface equals the flux of the curl through the surface.

> [!theorem] Stokes' theorem
> Let $S$ be a piecewise-smooth oriented surface bounded by a piecewise-smooth simple closed curve $C$, with the orientation of $C$ **induced** by that of $S$: if your right thumb points along the normal $\hat{\mathbf{n}}$, your fingers curl in the direction of $C$ (equivalently, walking along $C$ with your head along $\hat{\mathbf{n}}$, the surface is on your left). If $\mathbf{F}$ has continuous first partial derivatives on an open region containing $S$, then
> $$
> \oint_C\mathbf{F}\cdot d\mathbf{r}=\iint_S(\nabla\times\mathbf{F})\cdot d\mathbf{S}.
> $$

For a flat surface in the $xy$-plane with $\hat{\mathbf{n}}=\hat{\mathbf{z}}$ this is exactly Green's theorem.

### Proof sketches

*Sketch 1 — tiling.* Cover $S$ with a fine mesh of small, nearly flat patches $\Delta S_k$ with normals $\hat{\mathbf{n}}_k$. By the definition of curl (Section 7.5), the circulation around the edge of patch $k$ is approximately $(\nabla\times\mathbf{F})\cdot\hat{\mathbf{n}}_k\,\Delta S_k$. Add over all patches. Each internal edge is shared by two neighbouring patches and is traversed in opposite directions, so its contributions cancel; only the edges on the rim $C$ survive. In the limit the sum of the right-hand sides becomes $\iint_S(\nabla\times\mathbf{F})\cdot d\mathbf{S}$. ∎

*Sketch 2 — reduction to Green.* For a graph $z=g(x,y)$ over a region $D$, the parametrisation $(x,y,g(x,y))$ turns $\oint_C\mathbf{F}\cdot d\mathbf{r}$ into a plane line integral around the boundary of $D$, with $dz=g_x\,dx+g_y\,dy$:

$$
\oint_{\partial D}\big[(F_x+F_zg_x)\,dx+(F_y+F_zg_y)\,dy\big].
$$

Green's theorem converts this into a double integral over $D$. Expanding with the chain rule (where $F_x$ means $F_x(x,y,g(x,y))$, so $\partial_y$ acting on it produces $\partial_yF_x+\partial_zF_x\,g_y$), the terms with $g_{xy}$ cancel by Clairaut and what remains is precisely $(\nabla\times\mathbf{F})\cdot(-g_x,-g_y,1)$ integrated over $D$ — the flux of the curl through the graph. ∎

### Consequences

> [!theorem] Corollaries of Stokes' theorem
> 1. **Surface independence.** $\iint_S(\nabla\times\mathbf{F})\cdot d\mathbf{S}$ depends only on the boundary curve $C$, not on which surface spanning $C$ you choose. Pick the easiest one (usually a flat disk).
> 2. **Closed surfaces.** For a closed surface (no boundary), $\oint_S(\nabla\times\mathbf{F})\cdot d\mathbf{S}=0$ — consistent with $\nabla\cdot(\nabla\times\mathbf{F})=0$ and the divergence theorem.
> 3. **Curl-free implies conservative** on a simply connected region: every loop bounds a surface inside the region, through which the flux of $\nabla\times\mathbf{F}=\mathbf{0}$ is zero, so every circulation vanishes (Section 7.6).
> 4. **The curl definition is recovered.** Applying the theorem to a tiny flat disk of area $\Delta S$ and normal $\hat{\mathbf{n}}$ and dividing by $\Delta S$ gives back $(\nabla\times\mathbf{F})\cdot\hat{\mathbf{n}}=\lim\frac{1}{\Delta S}\oint\mathbf{F}\cdot d\mathbf{r}$.

Corollary 1 has a vivid electromagnetic meaning: the current "through a loop" is well defined only because $\nabla\cdot\mathbf{J}=0$ for steady currents — whichever surface you stretch across the loop, the same current crosses it. When charge accumulates (a charging capacitor), two surfaces spanning the same loop can carry different conduction currents, which is the puzzle Maxwell solved with the displacement current (Section 7.12).

> [!example] Worked Example — verifying Stokes' theorem on a hemisphere
> Let $\mathbf{F}=z\,\hat{\mathbf{x}}+x\,\hat{\mathbf{y}}+y\,\hat{\mathbf{z}}$ and $S$ the upper unit hemisphere $x^2+y^2+z^2=1$, $z\ge0$, with outward (upward) normal. The boundary is the unit circle in the $xy$-plane, counter-clockwise seen from above.
>
> **Line integral.** $\mathbf{r}=(\cos t,\sin t,0)$, $\mathbf{r}'=(-\sin t,\cos t,0)$, and on $C$, $\mathbf{F}=(0,\cos t,\sin t)$. So $\mathbf{F}\cdot\mathbf{r}'=\cos^2t$ and $\oint=\int_0^{2\pi}\cos^2t\,dt=\pi$.
>
> **Curl.** $\nabla\times\mathbf{F}=(\partial_yy-\partial_zx,\ \partial_zz-\partial_xy,\ \partial_xx-\partial_yz)=(1,1,1)$.
>
> **Flux of the curl through the hemisphere.** $d\mathbf{S}=\sin\theta\,d\theta\,d\phi\,\hat{\mathbf{R}}$ and $(1,1,1)\cdot\hat{\mathbf{R}}=\sin\theta\cos\phi+\sin\theta\sin\phi+\cos\theta$. The $\cos\phi$ and $\sin\phi$ terms integrate to zero over $0\le\phi\le2\pi$, leaving
> $$
> \int_0^{2\pi}\int_0^{\pi/2}\cos\theta\sin\theta\,d\theta\,d\phi=2\pi\cdot\tfrac12=\pi\ \checkmark
> $$
>
> **Using surface independence instead.** The flat unit disk also spans $C$, with normal $\hat{\mathbf{z}}$: $\iint(1,1,1)\cdot\hat{\mathbf{z}}\,dA=\pi$. Same answer in one line.

> [!example] Worked Example — a line integral made easy
> Evaluate $\oint_C\mathbf{F}\cdot d\mathbf{r}$ for $\mathbf{F}=z^2\,\hat{\mathbf{x}}+y^2\,\hat{\mathbf{y}}+x\,\hat{\mathbf{z}}$, where $C$ is the triangle $(1,0,0)\to(0,1,0)\to(0,0,1)\to(1,0,0)$.
>
> **Step 1 — curl.** $\nabla\times\mathbf{F}=(\partial_yx-\partial_zy^2,\ \partial_zz^2-\partial_xx,\ \partial_xy^2-\partial_yz^2)=(0,\ 2z-1,\ 0)$.
>
> **Step 2 — surface and orientation.** Take the flat triangle in the plane $x+y+z=1$, i.e. $z=1-x-y$ over the projected triangle $T$: $x,y\ge0$, $x+y\le1$. The edge vectors $(0,1,0)-(1,0,0)=(-1,1,0)$ and $(0,0,1)-(1,0,0)=(-1,0,1)$ have cross product $(1,1,1)$, so the given order of vertices corresponds to the upward normal. For the graph, $\mathbf{N}=(-g_x,-g_y,1)=(1,1,1)$ ✓.
>
> **Step 3 — integrate.** $(\nabla\times\mathbf{F})\cdot\mathbf{N}=2z-1=2(1-x-y)-1=1-2x-2y$. Over $T$ (area $\tfrac12$, with $\iint_Tx\,dA=\iint_Ty\,dA=\tfrac16$):
> $$
> \oint_C\mathbf{F}\cdot d\mathbf{r}=\tfrac12-2\cdot\tfrac16-2\cdot\tfrac16=-\tfrac16 .
> $$
>
> **Check** by computing the three edge integrals directly (each a polynomial in $t$): they sum to $-\tfrac16$ ✓. One surface integral replaced three line integrals.

> [!example] Worked Example — the vector potential of a solenoid
> Inside a long solenoid of radius $a$ the magnetic flux density is uniform, $\mathbf{B}=B_0\hat{\mathbf{z}}$, and outside it is (nearly) zero. A vector potential with $\nabla\times\mathbf{A}=\mathbf{B}$ is sought in the form $\mathbf{A}=A_\phi(r)\,\hat{\boldsymbol{\phi}}$. Find it inside and outside.
>
> **Step 1 — Stokes on a circle of radius $r$.** $\oint\mathbf{A}\cdot d\mathbf{r}=2\pi rA_\phi(r)$, and this equals $\iint(\nabla\times\mathbf{A})\cdot d\mathbf{S}=\iint\mathbf{B}\cdot d\mathbf{S}=\Phi(r)$, the flux enclosed.
>
> **Step 2 — inside ($r<a$).** $\Phi=B_0\pi r^2$, so $A_\phi=B_0r/2$.
>
> **Step 3 — outside ($r>a$).** $\Phi=B_0\pi a^2$ (all the flux), so $A_\phi=\dfrac{B_0a^2}{2r}$.
>
> **Check with the cylindrical curl.** Inside: $\frac1r\partial_r(r\cdot B_0r/2)=B_0$ ✓. Outside: $\frac1r\partial_r(B_0a^2/2)=0$ ✓.
>
> The vector potential is non-zero outside even though $\mathbf{B}$ is zero there — it "remembers" the enclosed flux through its circulation. Faraday's law (Section 7.12) uses this: a changing enclosed flux drives an EMF in a loop that lies entirely in the field-free region.

```desmos-graph
left=0; right=4; bottom=0; top=0.7
---
y=\frac{x}{2}|blue|0<x<1
y=\frac{1}{2x}|blue|1<x<4
x=1|grey|dashed
```

*Vector potential $A_\phi(r)$ of an ideal solenoid ($B_0=1$, $a=1$): it grows linearly inside and decays as $1/r$ outside, where $\mathbf{B}=\mathbf{0}$ but the circulation of $\mathbf{A}$ still equals the enclosed flux.*

> [!warning] Common Pitfall
> Mismatched orientation. If the surface normal points down but you traverse $C$ counter-clockwise as seen from above, Stokes' theorem gives the wrong sign. Fix the curve direction first, then use the right-hand rule to choose $\hat{\mathbf{n}}$ — or vice versa — never independently.

> [!tip] Engineering Note
> Stokes' theorem is why a current probe (Rogowski coil) works: a coil wound around a closed path measures $\oint\mathbf{B}\cdot d\mathbf{l}=\mu_0I_{\text{enc}}$ regardless of where the conductor sits inside the loop. Likewise the voltage induced in a PCB loop depends only on the flux through the loop area, so minimising loop area (return path directly under the signal trace) minimises inductive pickup.

## 7.11 Divergence theorem

### Statement

> [!theorem] Divergence theorem (Gauss's theorem)
> Let $V$ be a bounded solid region whose boundary $S$ is a closed piecewise-smooth surface with **outward** unit normal. If $\mathbf{F}$ has continuous first partial derivatives on an open region containing $V$, then
> $$
> \oint_S\mathbf{F}\cdot d\mathbf{S}=\iiint_V\nabla\cdot\mathbf{F}\,dV .
> $$

In words: **the net outward flux through a closed surface equals the total source strength inside.** It is the three-dimensional analogue of the flux form of Green's theorem, and the mathematical content of Gauss's law.

### Proof sketches

*Sketch 1 — tiling.* Divide $V$ into tiny boxes. By the definition of divergence, the outward flux of each box is $\nabla\cdot\mathbf{F}\,\Delta V_k$. Summing over all boxes, the flux through every internal face is counted twice with opposite normals and cancels; only the faces on the outer boundary $S$ survive. ∎

*Sketch 2 — one component at a time.* Write $\mathbf{F}=F_x\hat{\mathbf{x}}+F_y\hat{\mathbf{y}}+F_z\hat{\mathbf{z}}$; it suffices to prove $\oint F_z\hat{\mathbf{z}}\cdot d\mathbf{S}=\iiint\partial_zF_z\,dV$ and its two companions. Let $V=\{(x,y)\in D,\ g_1(x,y)\le z\le g_2(x,y)\}$. Then

$$
\iiint_V\frac{\partial F_z}{\partial z}\,dV=\iint_D\Big[F_z(x,y,g_2)-F_z(x,y,g_1)\Big]dx\,dy .
$$

On the top surface the upward normal vector is $(-g_{2x},-g_{2y},1)\,dx\,dy$, so $F_z\hat{\mathbf{z}}\cdot d\mathbf{S}=F_z(x,y,g_2)\,dx\,dy$; on the bottom the outward normal points down, giving $-F_z(x,y,g_1)\,dx\,dy$; on any vertical side walls $\hat{\mathbf{z}}\cdot\hat{\mathbf{n}}=0$. The flux equals the volume integral. ∎

### Consequences

**Solenoidal fields.** If $\nabla\cdot\mathbf{F}=0$ in $V$, the flux into $V$ equals the flux out: what enters through one part of the surface leaves through another. For $\mathbf{B}$ this is "no magnetic monopoles"; it is why magnetic flux through a transformer core is the same at every cross-section.

**Conservation laws.** Let $\rho$ be the density of something (charge, say) and $\mathbf{J}$ its flux density. Conservation says the rate of decrease of the amount inside $V$ equals the outflow:

$$
-\frac{d}{dt}\iiint_V\rho\,dV=\oint_S\mathbf{J}\cdot d\mathbf{S}=\iiint_V\nabla\cdot\mathbf{J}\,dV .
$$

Since $V$ is arbitrary, the integrands must agree, giving the **continuity equation**

$$
\nabla\cdot\mathbf{J}+\frac{\partial\rho}{\partial t}=0 .
$$

The step "integral is zero for every region ⇒ integrand is zero" (valid for continuous integrands: if the integrand were positive at a point it would be positive on a small ball around it) is the device that converts every integral law into a differential one. We use it four more times in Section 7.12.

**Green's identities.** Apply the theorem to $\mathbf{F}=f\nabla g$, using identity 4 of Section 7.5, $\nabla\cdot(f\nabla g)=f\nabla^2g+\nabla f\cdot\nabla g$:

$$
\iiint_V\left(f\nabla^2g+\nabla f\cdot\nabla g\right)dV=\oint_Sf\,\nabla g\cdot d\mathbf{S}\qquad\text{(first identity)}.
$$

Swapping $f$ and $g$ and subtracting:

$$
\iiint_V\left(f\nabla^2g-g\nabla^2f\right)dV=\oint_S\left(f\nabla g-g\nabla f\right)\cdot d\mathbf{S}\qquad\text{(second identity)}.
$$

> [!theorem] Uniqueness of solutions of Poisson's equation
> If $\nabla^2V=-\rho/\varepsilon$ in $V$ and $V$ is prescribed on the boundary $S$, the solution is unique.

*Proof.* Let $V_1,V_2$ be two solutions and $w=V_1-V_2$. Then $\nabla^2w=0$ inside and $w=0$ on $S$. The first identity with $f=g=w$ gives $\iiint\lVert\nabla w\rVert^2dV=\oint w\nabla w\cdot d\mathbf{S}=0$, so $\nabla w=\mathbf{0}$, $w$ is constant, and it is zero because it vanishes on $S$. ∎ This theorem is what licenses the method of images and every "guess the solution and check it" approach in electrostatics ([Ch. 13](../part-4-differential-equations/13-partial-differential-equations.md)).

### Point sources and the delta function

The inverse-square field $\hat{\mathbf{R}}/R^2$ has zero divergence for $R\ne0$, yet its flux through any sphere about the origin is $4\pi$. The divergence theorem cannot be applied naively to a region containing the origin because the field is singular there. Two ways to handle it:

- **Excise the singularity** (as with the hole in Green's theorem): for any closed surface $S$ enclosing the origin, apply the theorem to the region between $S$ and a small sphere $S_\epsilon$. The divergence is zero there, so the flux through $S$ equals the flux through $S_\epsilon$, which is $4\pi$. If $S$ does not enclose the origin, the flux is $0$.
- **Use the Dirac delta** ([Ch. 14](../part-5-transforms/14-laplace-transforms.md), [Ch. 16](../part-5-transforms/16-fourier-transforms.md)): define $\delta(\mathbf{R})=\delta(x)\delta(y)\delta(z)$, with $\iiint_V\delta(\mathbf{R})\,dV=1$ if $V$ contains the origin and $0$ otherwise. Then the flux results above are summarised as
$$
\nabla\cdot\frac{\hat{\mathbf{R}}}{R^2}=4\pi\,\delta(\mathbf{R}),\qquad\nabla^2\frac1R=-4\pi\,\delta(\mathbf{R}).
$$

The second formula says that $G(\mathbf{R},\mathbf{R}')=\dfrac{1}{4\pi\lVert\mathbf{R}-\mathbf{R}'\rVert}$ is the **Green's function** of $-\nabla^2$: the potential of a unit point source. Superposing point sources gives the solution of Poisson's equation in free space, $V(\mathbf{R})=\dfrac1{4\pi\varepsilon_0}\iiint\dfrac{\rho(\mathbf{R}')}{\lVert\mathbf{R}-\mathbf{R}'\rVert}dV'$ — the same superposition-of-impulse-responses idea as convolution in circuit theory.

> [!example] Worked Example — verifying the theorem on a cube
> Let $\mathbf{F}=x\hat{\mathbf{x}}+y\hat{\mathbf{y}}+z\hat{\mathbf{z}}$ and $V=[0,1]^3$.
>
> **Volume side.** $\nabla\cdot\mathbf{F}=3$, so $\iiint3\,dV=3$.
>
> **Surface side.** On the face $x=1$, $\hat{\mathbf{n}}=\hat{\mathbf{x}}$ and $\mathbf{F}\cdot\hat{\mathbf{n}}=x=1$: flux $1$. On $x=0$, $\mathbf{F}\cdot(-\hat{\mathbf{x}})=-x=0$: flux $0$. The $y$ and $z$ faces behave the same way. Total $1+1+1=3$ ✓.

> [!example] Worked Example — a flux that is easier as a volume integral
> Find the outward flux of $\mathbf{F}=x^3\hat{\mathbf{x}}+y^3\hat{\mathbf{y}}+z^3\hat{\mathbf{z}}$ through the sphere $R=a$.
>
> **Direct route (hard).** $\mathbf{F}\cdot\hat{\mathbf{R}}=(x^4+y^4+z^4)/a$ on the sphere — an awkward angular integral.
>
> **Divergence theorem.** $\nabla\cdot\mathbf{F}=3(x^2+y^2+z^2)=3R^2$, so
> $$
> \oint_S\mathbf{F}\cdot d\mathbf{S}=\int_0^{2\pi}\int_0^\pi\int_0^a3R^2\cdot R^2\sin\theta\,dR\,d\theta\,d\phi=3\cdot\frac{a^5}{5}\cdot2\cdot2\pi=\frac{12\pi a^5}{5}.
> $$
> **Cross-check of the direct route.** By symmetry $\oint x^4\,dS=\oint y^4\,dS=\oint z^4\,dS$, and $\oint z^4\,dS=a^4\cdot a^2\int_0^{2\pi}\int_0^\pi\cos^4\theta\sin\theta\,d\theta\,d\phi=a^6\cdot\frac25\cdot2\pi=\frac{4\pi a^6}{5}$. So the flux is $\frac{3}{a}\cdot\frac{4\pi a^6}{5}=\frac{12\pi a^5}{5}$ ✓.

> [!example] Worked Example — closing an open surface
> Recover the flux $\pi/2$ of $\mathbf{F}=z\hat{\mathbf{z}}$ through the paraboloid cap $z=1-x^2-y^2$, $z\ge0$ (Section 7.7), using the divergence theorem.
>
> **Step 1 — close the surface** with the disk $z=0$, $r\le1$. The closed surface encloses the solid $V$ under the cap.
>
> **Step 2 — volume integral.** $\nabla\cdot\mathbf{F}=1$, so $\iiint_V1\,dV$ is the volume: $\int_0^{2\pi}\int_0^1(1-r^2)\,r\,dr\,d\phi=\pi/2$.
>
> **Step 3 — subtract the disk.** On the disk $z=0$, so $\mathbf{F}=\mathbf{0}$ and its flux is $0$.
>
> **Result.** Flux through the cap $=\pi/2-0=\pi/2$ ✓. "Close the surface, use the volume integral, subtract the easy lid" is a standard trick.

> [!example] Worked Example — Gauss's law for an arbitrary surface
> A point charge $q$ at the origin has $\mathbf{E}=\dfrac{q}{4\pi\varepsilon_0}\dfrac{\hat{\mathbf{R}}}{R^2}$. Find the flux of $\mathbf{E}$ out of (a) the cube $\lvert x\rvert,\lvert y\rvert,\lvert z\rvert\le1$; (b) the cube $1\le x\le3$, $\lvert y\rvert,\lvert z\rvert\le1$.
>
> **(a)** The cube encloses the origin. By the excision argument the flux equals that through a small sphere: $\dfrac{q}{4\pi\varepsilon_0}\cdot4\pi=\dfrac{q}{\varepsilon_0}$. (Consistent with the per-face result $q/(6\varepsilon_0)$ of Section 7.7.)
>
> **(b)** The origin lies outside; $\nabla\cdot\mathbf{E}=0$ throughout, so the flux is $0$ — the field lines that enter through the face $x=1$ all leave through the other five faces.

> [!warning] Common Pitfall
> Using an inward normal, or forgetting a piece of the boundary. The divergence theorem needs the **whole** closed boundary with **outward** normals. For a region with an internal cavity (a hollow sphere), the boundary includes the inner surface with its normal pointing *into the cavity* — away from the material.

> [!tip] Engineering Note
> The finite-volume method used in thermal simulators and in many TCAD semiconductor-device solvers is the divergence theorem applied cell by cell: each mesh cell's stored quantity changes by the sum of fluxes through its faces. Because each face flux is shared by two cells with opposite signs, the scheme conserves charge (or heat) exactly — the discrete version of the interior-cancellation argument.

## 7.12 Applications to electromagnetic field theory

Everything in this chapter was built for this section. Maxwell's equations are four statements about the fields $\mathbf{E},\mathbf{D},\mathbf{H},\mathbf{B}$; in **integral form** they are statements about line integrals around loops and fluxes through surfaces — the form in which they were discovered and are measured. Stokes' theorem and the divergence theorem turn them into **differential form** — local partial differential equations that hold at every point and can be solved. The mathematics is the subject here; the physics is the motivation.

### The field quantities

| Symbol | Name | SI unit |
|---|---|---|
| $\mathbf{E}$ | electric field intensity | V/m |
| $\mathbf{D}$ | electric flux density | $\mathrm{C/m^2}$ |
| $\mathbf{H}$ | magnetic field intensity | A/m |
| $\mathbf{B}$ | magnetic flux density | T $=\mathrm{Wb/m^2}$ |
| $\mathbf{J}$ | (free) current density | $\mathrm{A/m^2}$ |
| $\rho$ | (free) volume charge density | $\mathrm{C/m^3}$ |

In a linear, isotropic medium the **constitutive relations** are $\mathbf{D}=\varepsilon\mathbf{E}$, $\mathbf{B}=\mu\mathbf{H}$, $\mathbf{J}=\sigma\mathbf{E}$ (Ohm's law), with $\varepsilon=\varepsilon_r\varepsilon_0$, $\mu=\mu_r\mu_0$, $\varepsilon_0=8.854\times10^{-12}$ F/m and $\mu_0\approx4\pi\times10^{-7}$ H/m.

### Maxwell's equations in integral form

> [!theorem] Maxwell's equations — integral form
> For any closed surface $S$ enclosing a volume $V$, and any open surface $S$ bounded by a closed curve $C$ (orientations related by the right-hand rule):
> $$
> \begin{aligned}
> &\text{Gauss's law:} && \oint_S\mathbf{D}\cdot d\mathbf{S}=\iiint_V\rho\,dV=Q_{\text{enc}},\\
> &\text{no magnetic charge:} && \oint_S\mathbf{B}\cdot d\mathbf{S}=0,\\
> &\text{Faraday's law:} && \oint_C\mathbf{E}\cdot d\mathbf{l}=-\frac{d}{dt}\iint_S\mathbf{B}\cdot d\mathbf{S},\\
> &\text{Ampère–Maxwell law:} && \oint_C\mathbf{H}\cdot d\mathbf{l}=\iint_S\left(\mathbf{J}+\frac{\partial\mathbf{D}}{\partial t}\right)\cdot d\mathbf{S}.
> \end{aligned}
> $$

Read with the vocabulary of this chapter:

- **Gauss:** the outward flux of $\mathbf{D}$ through a closed surface is the enclosed charge.
- **No monopoles:** the net flux of $\mathbf{B}$ out of any closed surface is zero; flux lines of $\mathbf{B}$ close on themselves.
- **Faraday:** the circulation of $\mathbf{E}$ (the EMF) around a loop is minus the rate of change of magnetic flux through it.
- **Ampère–Maxwell:** the circulation of $\mathbf{H}$ around a loop is the total current — conduction plus displacement — through it.

### From integral form to differential form

The recipe is the same each time: (1) turn the boundary integral into an integral over the inside using Stokes or Gauss; (2) move everything to one side under a single integral; (3) since the region is arbitrary, the integrand vanishes.

**Gauss's law.** By the divergence theorem, $\oint_S\mathbf{D}\cdot d\mathbf{S}=\iiint_V\nabla\cdot\mathbf{D}\,dV$. Hence $\iiint_V(\nabla\cdot\mathbf{D}-\rho)\,dV=0$ for every $V$, so

$$
\nabla\cdot\mathbf{D}=\rho .
$$

**No magnetic charge.** The same argument gives $\nabla\cdot\mathbf{B}=0$.

**Faraday's law.** For a stationary surface the time derivative passes inside the integral. By Stokes' theorem $\oint_C\mathbf{E}\cdot d\mathbf{l}=\iint_S(\nabla\times\mathbf{E})\cdot d\mathbf{S}$, so $\iint_S\left(\nabla\times\mathbf{E}+\dfrac{\partial\mathbf{B}}{\partial t}\right)\cdot d\mathbf{S}=0$ for every surface $S$. Choosing tiny flat surfaces with each of the three normal directions shows that each component of the bracket vanishes:

$$
\nabla\times\mathbf{E}=-\frac{\partial\mathbf{B}}{\partial t}.
$$

**Ampère–Maxwell law.** Identically,

$$
\nabla\times\mathbf{H}=\mathbf{J}+\frac{\partial\mathbf{D}}{\partial t}.
$$

> [!theorem] Maxwell's equations — differential form
> $$
> \begin{aligned}
> \nabla\cdot\mathbf{D}&=\rho, & \nabla\times\mathbf{E}&=-\frac{\partial\mathbf{B}}{\partial t},\\
> \nabla\cdot\mathbf{B}&=0, & \nabla\times\mathbf{H}&=\mathbf{J}+\frac{\partial\mathbf{D}}{\partial t}.
> \end{aligned}
> $$
> For time-harmonic fields with phasor convention $\mathbf{E}(\mathbf{R},t)=\operatorname{Re}\{\mathbf{E}(\mathbf{R})e^{j\omega t}\}$, $\partial/\partial t\to j\omega$:
> $$
> \nabla\times\mathbf{E}=-j\omega\mathbf{B},\qquad\nabla\times\mathbf{H}=\mathbf{J}+j\omega\mathbf{D}.
> $$

| Law | Integral form | Theorem used | Differential form |
|---|---|---|---|
| Gauss | $\oint\mathbf{D}\cdot d\mathbf{S}=Q_{\text{enc}}$ | divergence | $\nabla\cdot\mathbf{D}=\rho$ |
| No monopoles | $\oint\mathbf{B}\cdot d\mathbf{S}=0$ | divergence | $\nabla\cdot\mathbf{B}=0$ |
| Faraday | $\oint\mathbf{E}\cdot d\mathbf{l}=-\dfrac{d\Phi}{dt}$ | Stokes | $\nabla\times\mathbf{E}=-\partial_t\mathbf{B}$ |
| Ampère–Maxwell | $\oint\mathbf{H}\cdot d\mathbf{l}=I_{\text{enc}}+\dfrac{d}{dt}\iint\mathbf{D}\cdot d\mathbf{S}$ | Stokes | $\nabla\times\mathbf{H}=\mathbf{J}+\partial_t\mathbf{D}$ |

The integral forms are best for problems with high symmetry (where the field can be pulled out of the integral) and for circuits and measurements; the differential forms are best for deriving wave equations, boundary-value problems and numerical methods. The two divergence equations are, by the Helmholtz theorem, exactly what is needed alongside the two curl equations to fix the fields.

### Why the displacement current is required

Take the divergence of Ampère's law without Maxwell's term, $\nabla\times\mathbf{H}=\mathbf{J}$. Since the divergence of a curl is zero (identity 2 of Section 7.5), this would force $\nabla\cdot\mathbf{J}=0$ — but the continuity equation (Section 7.11) says $\nabla\cdot\mathbf{J}=-\partial\rho/\partial t$, which is non-zero whenever charge accumulates (as on a charging capacitor plate). With Maxwell's term,

$$
0=\nabla\cdot(\nabla\times\mathbf{H})=\nabla\cdot\mathbf{J}+\frac{\partial}{\partial t}\nabla\cdot\mathbf{D}=\nabla\cdot\mathbf{J}+\frac{\partial\rho}{\partial t},
$$

which is precisely charge conservation. A vector identity, not an experiment, demanded the extra term — and that term is what makes electromagnetic waves possible.

### Potentials, Poisson's and Laplace's equations

Since $\nabla\cdot\mathbf{B}=0$, there is a vector potential with $\mathbf{B}=\nabla\times\mathbf{A}$ (Section 7.5). Substituting into Faraday's law gives $\nabla\times(\mathbf{E}+\partial_t\mathbf{A})=\mathbf{0}$, so the bracket is a gradient:

$$
\mathbf{E}=-\nabla V-\frac{\partial\mathbf{A}}{\partial t}.
$$

In electrostatics $\mathbf{E}=-\nabla V$; inserting into $\nabla\cdot(\varepsilon\mathbf{E})=\rho$ with constant $\varepsilon$ gives **Poisson's equation**, and in charge-free regions **Laplace's equation**:

$$
\nabla^2V=-\frac{\rho}{\varepsilon},\qquad\nabla^2V=0 .
$$

These are solved by separation of variables in [Ch. 13](../part-4-differential-equations/13-partial-differential-equations.md) and numerically in [Ch. 22](../part-7-numerical-computational/22-numerical-methods.md); uniqueness was proved in Section 7.11.

### The wave equation

In a source-free ($\rho=0$, $\mathbf{J}=\mathbf{0}$), linear, homogeneous medium, take the curl of Faraday's law and use Ampère–Maxwell:

$$
\nabla\times(\nabla\times\mathbf{E})=-\frac{\partial}{\partial t}\nabla\times\mathbf{B}=-\mu\frac{\partial}{\partial t}\left(\varepsilon\frac{\partial\mathbf{E}}{\partial t}\right)=-\mu\varepsilon\frac{\partial^2\mathbf{E}}{\partial t^2}.
$$

By identity 7, $\nabla\times(\nabla\times\mathbf{E})=\nabla(\nabla\cdot\mathbf{E})-\nabla^2\mathbf{E}=-\nabla^2\mathbf{E}$, because $\nabla\cdot\mathbf{E}=0$ here. Therefore

$$
\nabla^2\mathbf{E}=\mu\varepsilon\frac{\partial^2\mathbf{E}}{\partial t^2},\qquad u=\frac{1}{\sqrt{\mu\varepsilon}} .
$$

In vacuum $u=1/\sqrt{\mu_0\varepsilon_0}=2.998\times10^8$ m/s — the speed of light, predicted from two constants measured with capacitors and coils. In phasor form ($\partial_t\to j\omega$) it becomes the **Helmholtz equation** $\nabla^2\mathbf{E}+k^2\mathbf{E}=\mathbf{0}$ with wavenumber $k=\omega\sqrt{\mu\varepsilon}$. The same equation holds for $\mathbf{H}$.

### Boundary conditions from the integral forms

At an interface between two media the fields can jump, and the differential equations (which need derivatives) do not apply on the interface itself. The integral forms do.

- **Normal components — pillbox.** Apply Gauss's law to a thin pillbox straddling the interface, faces of area $\Delta S$ on either side, height $\to0$. The side-wall flux vanishes, leaving $(\mathbf{D}_1-\mathbf{D}_2)\cdot\hat{\mathbf{n}}\,\Delta S=\rho_s\Delta S$, where $\hat{\mathbf{n}}$ points from medium 2 into medium 1 and $\rho_s$ is surface charge. Likewise for $\mathbf{B}$ with zero right side.
- **Tangential components — thin loop.** Apply Faraday's law to a thin rectangular loop straddling the interface, sides $\Delta l$ parallel to it, height $\to0$. The enclosed flux vanishes with the area, leaving $(\mathbf{E}_1-\mathbf{E}_2)\cdot\hat{\mathbf{t}}\,\Delta l=0$. For $\mathbf{H}$ the right side is a surface current $J_s\Delta l$.

$$
\hat{\mathbf{n}}\cdot(\mathbf{D}_1-\mathbf{D}_2)=\rho_s,\quad\hat{\mathbf{n}}\cdot(\mathbf{B}_1-\mathbf{B}_2)=0,\quad\hat{\mathbf{n}}\times(\mathbf{E}_1-\mathbf{E}_2)=\mathbf{0},\quad\hat{\mathbf{n}}\times(\mathbf{H}_1-\mathbf{H}_2)=\mathbf{J}_s .
$$

These are exactly the "shrink the region onto the boundary" arguments by which we defined divergence and curl, now used on a surface where the fields are discontinuous.

### Energy flow: Poynting's theorem

Identity 6 of Section 7.5 and the two curl equations give

$$
\nabla\cdot(\mathbf{E}\times\mathbf{H})=\mathbf{H}\cdot(\nabla\times\mathbf{E})-\mathbf{E}\cdot(\nabla\times\mathbf{H})=-\mathbf{H}\cdot\frac{\partial\mathbf{B}}{\partial t}-\mathbf{E}\cdot\mathbf{J}-\mathbf{E}\cdot\frac{\partial\mathbf{D}}{\partial t}.
$$

For a linear medium $\mathbf{H}\cdot\partial_t\mathbf{B}=\partial_t(\tfrac12\mu H^2)$ and $\mathbf{E}\cdot\partial_t\mathbf{D}=\partial_t(\tfrac12\varepsilon E^2)$. Integrating over a volume and applying the divergence theorem:

$$
-\oint_S(\mathbf{E}\times\mathbf{H})\cdot d\mathbf{S}=\frac{d}{dt}\iiint_V\left(\tfrac12\varepsilon E^2+\tfrac12\mu H^2\right)dV+\iiint_V\mathbf{E}\cdot\mathbf{J}\,dV .
$$

Power flowing **in** through the surface (flux of the Poynting vector $\mathbf{S}=\mathbf{E}\times\mathbf{H}$, in $\mathrm{W/m^2}$) equals the rate of increase of stored field energy plus the power dissipated. It is the field version of "power in = $\frac{d}{dt}$(energy stored) + $I^2R$".

### Worked examples

> [!example] Worked Example — Gauss's law: the uniformly charged sphere
> A sphere of radius $a$ in free space holds charge $Q$ uniformly ($\rho_0=3Q/(4\pi a^3)$). Find $\mathbf{E}$ everywhere and check the differential form.
>
> **Step 1 — symmetry.** By spherical symmetry $\mathbf{E}=E_R(R)\hat{\mathbf{R}}$. On a Gaussian sphere of radius $R$, $\mathbf{D}\cdot d\mathbf{S}=\varepsilon_0E_R\,dS$ with $E_R$ constant, so $\oint\mathbf{D}\cdot d\mathbf{S}=4\pi R^2\varepsilon_0E_R$.
>
> **Step 2 — enclosed charge.** For $R\ge a$, $Q_{\text{enc}}=Q$. For $R<a$, $Q_{\text{enc}}=Q\,R^3/a^3$.
>
> **Step 3 — solve.**
> $$
> E_R=\begin{cases}\dfrac{QR}{4\pi\varepsilon_0a^3}, & R<a,\\[6pt]\dfrac{Q}{4\pi\varepsilon_0R^2}, & R\ge a.\end{cases}
> $$
> Both branches give $Q/(4\pi\varepsilon_0a^2)$ at $R=a$: the field is continuous (no surface charge).
>
> **Step 4 — differential check.** Inside, $\nabla\cdot\mathbf{E}=\frac1{R^2}\partial_R\big(R^2\cdot\frac{QR}{4\pi\varepsilon_0a^3}\big)=\frac{3Q}{4\pi\varepsilon_0a^3}=\rho_0/\varepsilon_0$ ✓. Outside, $\frac1{R^2}\partial_R(\text{const})=0$ ✓.
>
> **Numbers.** $Q=1$ nC, $a=1$ cm: peak field $E(a)=8.988\times10^9\times10^{-9}/10^{-4}=8.99\times10^4$ V/m $\approx90$ kV/m, about 3 % of the breakdown strength of air.

```desmos-graph
left=0; right=4; bottom=0; top=1.2
---
y=x|blue|0<x<1
y=\frac{1}{x^2}|blue|1<x<4
x=1|grey|dashed
```

*Field of a uniformly charged sphere, normalised to its surface value, against $R/a$: linear inside (where enclosed charge grows as $R^3$), inverse-square outside.*

> [!example] Worked Example — Ampère's law: the field in and around a thick wire
> A long straight wire of radius $a$ carries current $I$ uniformly ($J=I/(\pi a^2)$, along $\hat{\mathbf{z}}$). Find $\mathbf{H}$ and $\mathbf{B}$ ($\mu=\mu_0$), and verify $\nabla\times\mathbf{H}=\mathbf{J}$.
>
> **Step 1 — symmetry.** $\mathbf{H}=H_\phi(r)\hat{\boldsymbol{\phi}}$. Around a circle of radius $r$, $\oint\mathbf{H}\cdot d\mathbf{l}=2\pi rH_\phi$.
>
> **Step 2 — enclosed current.** $I_{\text{enc}}=I\,r^2/a^2$ for $r<a$, and $I$ for $r\ge a$ (static, so no displacement current).
>
> **Step 3 — solve.** $H_\phi=\dfrac{Ir}{2\pi a^2}$ inside, $\dfrac{I}{2\pi r}$ outside.
>
> **Step 4 — differential check.** Inside: $(\nabla\times\mathbf{H})_z=\frac1r\partial_r\big(r\cdot\frac{Ir}{2\pi a^2}\big)=\frac{I}{\pi a^2}=J$ ✓. Outside: $\frac1r\partial_r(I/2\pi)=0$ ✓ — the curl-free wire field of Sections 7.5–7.6, now seen as the field of a current that occupies only $r<a$.
>
> **Numbers.** $I=10$ A, $a=1$ mm: $B(a)=\mu_0I/(2\pi a)=2.0$ mT; $B(0.5\ \mathrm{mm})=1.0$ mT; $B(1\ \mathrm{cm})=0.20$ mT.

```desmos-graph
left=0; right=5; bottom=0; top=1.2
---
y=x|red|0<x<1
y=\frac{1}{x}|red|1<x<5
x=1|grey|dashed
```

*$\lVert\mathbf{B}\rVert$ of a thick wire normalised to its surface value, against $r/a$: linear inside, $1/r$ outside. Compare with the sphere: the outside decay is slower because the source is a line, not a point.*

> [!example] Worked Example — capacitance per metre of a coaxial cable
> A coaxial cable has inner radius $a$, outer radius $b$, and dielectric $\varepsilon$. Find the capacitance per metre and evaluate it for $a=0.45$ mm, $b=1.47$ mm, $\varepsilon_r=2.25$ (polyethylene).
>
> **Step 1 — Gauss.** Put charge $\rho_\ell$ per metre on the inner conductor. A coaxial Gaussian cylinder of radius $r$ and length $L$ ($a<r<b$) has flux only through its curved side: $\varepsilon E_r\cdot2\pi rL=\rho_\ell L$, so $E_r=\dfrac{\rho_\ell}{2\pi\varepsilon r}$.
>
> **Step 2 — line integral for the voltage.**
> $$
> V=-\int_b^a E_r\,dr=\frac{\rho_\ell}{2\pi\varepsilon}\int_a^b\frac{dr}{r}=\frac{\rho_\ell}{2\pi\varepsilon}\ln\frac ba .
> $$
>
> **Step 3 — capacitance.** $C'=\rho_\ell/V=\dfrac{2\pi\varepsilon}{\ln(b/a)}$.
>
> **Numbers.** $\ln(1.47/0.45)=1.1838$, $2\pi\times2.25\times8.854\times10^{-12}=1.2517\times10^{-10}$, so $C'=105.7$ pF/m — the right order for common 50 Ω cable (about 100 pF/m).
>
> Gauss's law (a surface integral) gave the field; a line integral gave the voltage. That pair of steps solves every symmetric capacitor problem.

> [!example] Worked Example — Faraday's law: induced EMF and induced field
> A uniform field $\mathbf{B}=B_0\sin(\omega t)\,\hat{\mathbf{z}}$ fills a region. (a) Find the EMF around a circular loop of radius $a$ in the $xy$-plane centred on the axis. (b) Find the induced electric field and check $\nabla\times\mathbf{E}=-\partial_t\mathbf{B}$. (c) Evaluate the peak EMF for $a=5$ cm, $B_0=10$ mT, $f=50$ Hz.
>
> **(a)** Flux $\Phi=\pi a^2B_0\sin\omega t$, so $\text{EMF}=\oint\mathbf{E}\cdot d\mathbf{l}=-d\Phi/dt=-\pi a^2B_0\omega\cos\omega t$.
>
> **(b)** By symmetry $\mathbf{E}=E_\phi(r)\hat{\boldsymbol{\phi}}$, so the EMF around radius $r$ is $2\pi rE_\phi=-\pi r^2\dot B$, giving $E_\phi=-\tfrac r2\dot B=-\tfrac r2B_0\omega\cos\omega t$. Check: $(\nabla\times\mathbf{E})_z=\frac1r\partial_r(rE_\phi)=\frac1r\partial_r\left(-\tfrac{r^2}{2}\dot B\right)=-\dot B$ ✓.
>
> **(c)** Peak EMF $=\pi(0.05)^2(0.01)(2\pi\cdot50)=24.7$ mV.
>
> This induced $\mathbf{E}$ has circular field lines and non-zero curl: it is **not** a gradient, and no scalar potential can describe it — Kirchhoff's voltage law fails for a loop around changing flux.

> [!example] Worked Example — displacement current in a capacitor
> A parallel-plate air capacitor has plate area $A=1\ \mathrm{cm^2}$ and gap $d=0.1$ mm, driven by $v(t)=V_0\sin\omega t$ with $V_0=10$ V, $f=1$ MHz. Show that the displacement current between the plates equals the conduction current in the leads.
>
> **Step 1.** Between the plates $E=v/d$, so the displacement current density is $\partial D/\partial t=\varepsilon_0\dot v/d$.
>
> **Step 2.** Its flux through a surface between the plates is $I_d=A\varepsilon_0\dot v/d=C\dot v$ with $C=\varepsilon_0A/d$ — exactly the lead current $i=C\,dv/dt$ ✓.
>
> **Step 3 — numbers.** $C=8.854\times10^{-12}\times10^{-4}/10^{-4}=8.854$ pF; peak $I=C\omega V_0=8.854\times10^{-12}\times2\pi\times10^6\times10=0.556$ mA.
>
> So a loop around the lead "sees" the same $\oint\mathbf{H}\cdot d\mathbf{l}$ whether the surface it bounds is pierced by the wire (conduction current) or bulges between the plates (displacement current) — Stokes' surface independence restored.

> [!example] Worked Example — Poynting's theorem for a resistive wire
> A straight wire of radius $a$, length $L$ and conductivity $\sigma$ carries a DC current $I$. Show that the Poynting flux into its surface equals $I^2R$.
>
> **Step 1 — fields at the surface.** $\mathbf{E}=\mathbf{J}/\sigma=\dfrac{I}{\sigma\pi a^2}\hat{\mathbf{z}}$; tangential $\mathbf{E}$ is continuous, so this is also the field just outside. From the Ampère example, $\mathbf{H}=\dfrac{I}{2\pi a}\hat{\boldsymbol{\phi}}$ at $r=a$.
>
> **Step 2 — Poynting vector.** $\mathbf{S}=\mathbf{E}\times\mathbf{H}=EH\,(\hat{\mathbf{z}}\times\hat{\boldsymbol{\phi}})=-EH\,\hat{\mathbf{r}}$: it points **into** the wire.
>
> **Step 3 — inward flux.** $P=EH\cdot2\pi aL=\dfrac{I}{\sigma\pi a^2}\cdot\dfrac{I}{2\pi a}\cdot2\pi aL=I^2\dfrac{L}{\sigma\pi a^2}=I^2R$ ✓.
>
> The energy dissipated in the wire flows in from the surrounding field, not along the inside of the conductor — a striking conclusion reached purely by a cross product and a surface integral.

> [!warning] Common Pitfall
> Using Gauss's or Ampère's law to "find" a field without enough symmetry. The integral forms are always true, but they determine the field only when symmetry lets you pull $E$ or $H$ out of the integral as a constant. For a finite wire, a square loop or a charged cube, $\oint\mathbf{H}\cdot d\mathbf{l}=I$ is still correct but does not tell you $\mathbf{H}$ at any point; use Biot–Savart, superposition or a numerical solver instead.

> [!warning] Common Pitfall
> Mixing time-domain and phasor forms. In the phasor equations $\mathbf{E}$ is a complex vector and $\partial_t$ has become $j\omega$; quantities like energy and power, which are quadratic, need $\tfrac12\operatorname{Re}\{\cdot\}$ and conjugates ($\langle\mathbf{S}\rangle=\tfrac12\operatorname{Re}\{\mathbf{E}\times\mathbf{H}^*\}$), not the time-domain formulas applied to phasors.

> [!tip] Engineering Note
> Circuit theory is the low-frequency limit of these equations. KCL is $\nabla\cdot\mathbf{J}=0$ (continuity with no charge accumulation at nodes, or with capacitor displacement currents counted as branch currents); KVL is $\oint\mathbf{E}\cdot d\mathbf{l}=0$ (Faraday with flux changes lumped into inductors); $v=L\,di/dt$ and $i=C\,dv/dt$ are Faraday's and the displacement term of Ampère–Maxwell confined to components. When a circuit's dimensions approach a tenth of a wavelength these lumped approximations fail and the full field equations — via transmission-line theory and the wave equation — take over.

## Summary

| Concept | Key formula |
|---|---|
| Dot product | $\mathbf{A}\cdot\mathbf{B}=AB\cos\theta=A_xB_x+A_yB_y+A_zB_z$ |
| Cross product | $\lVert\mathbf{A}\times\mathbf{B}\rVert=AB\sin\theta$; determinant with $\hat{\mathbf{x}},\hat{\mathbf{y}},\hat{\mathbf{z}}$ |
| Triple products | $\mathbf{A}\cdot(\mathbf{B}\times\mathbf{C})=$ volume; $\mathbf{A}\times(\mathbf{B}\times\mathbf{C})=\mathbf{B}(\mathbf{A}\cdot\mathbf{C})-\mathbf{C}(\mathbf{A}\cdot\mathbf{B})$ |
| Arc length, curvature | $s=\int\lVert\mathbf{r}'\rVert dt$, $\kappa=\lVert\mathbf{r}'\times\mathbf{r}''\rVert/\lVert\mathbf{r}'\rVert^3$ |
| Polar velocity | $\mathbf{v}=\dot r\hat{\mathbf{r}}+r\dot\phi\hat{\boldsymbol{\phi}}+\dot z\hat{\mathbf{z}}$ |
| Field lines | $dx/F_x=dy/F_y=dz/F_z$ |
| Gradient | $\nabla f=\sum\frac{1}{h_i}\partial_{u_i}f\,\hat{\mathbf{u}}_i$; $h=(1,r,1)$ cyl., $(1,R,R\sin\theta)$ sph. |
| Divergence (def.) | $\nabla\cdot\mathbf{F}=\lim\frac1{\Delta V}\oint\mathbf{F}\cdot d\mathbf{S}$ |
| Divergence (sph.) | $\frac1{R^2}\partial_R(R^2F_R)+\frac1{R\sin\theta}\partial_\theta(\sin\theta F_\theta)+\frac1{R\sin\theta}\partial_\phi F_\phi$ |
| Curl (def.) | $(\nabla\times\mathbf{F})\cdot\hat{\mathbf{n}}=\lim\frac1{\Delta S}\oint\mathbf{F}\cdot d\mathbf{r}$ |
| Curl $z$-comp. (cyl.) | $\frac1r[\partial_r(rF_\phi)-\partial_\phi F_r]$ |
| Laplacian (sph.) | $\frac1{R^2}\partial_R(R^2\partial_Rf)+\frac1{R^2\sin\theta}\partial_\theta(\sin\theta\,\partial_\theta f)+\frac1{R^2\sin^2\theta}\partial_\phi^2f$ |
| Zero identities | $\nabla\times\nabla f=\mathbf{0}$, $\nabla\cdot(\nabla\times\mathbf{A})=0$ |
| Curl curl | $\nabla\times\nabla\times\mathbf{A}=\nabla(\nabla\cdot\mathbf{A})-\nabla^2\mathbf{A}$ |
| Line integral | $\int_C\mathbf{F}\cdot d\mathbf{r}=\int_a^b\mathbf{F}(\mathbf{r}(t))\cdot\mathbf{r}'(t)\,dt$ |
| Fundamental theorem | $\int_C\nabla f\cdot d\mathbf{r}=f(\mathbf{b})-f(\mathbf{a})$ |
| Surface element | $d\mathbf{S}=\mathbf{r}_u\times\mathbf{r}_v\,du\,dv$; sphere $a^2\sin\theta\,d\theta\,d\phi\,\hat{\mathbf{R}}$ |
| Volume element | $dV=r\,dr\,d\phi\,dz=R^2\sin\theta\,dR\,d\theta\,d\phi$ |
| Green | $\oint(P\,dx+Q\,dy)=\iint(Q_x-P_y)\,dA$; $A=\frac12\oint(x\,dy-y\,dx)$ |
| Stokes | $\oint_C\mathbf{F}\cdot d\mathbf{r}=\iint_S(\nabla\times\mathbf{F})\cdot d\mathbf{S}$ |
| Divergence theorem | $\oint_S\mathbf{F}\cdot d\mathbf{S}=\iiint_V\nabla\cdot\mathbf{F}\,dV$ |
| Continuity | $\nabla\cdot\mathbf{J}+\partial_t\rho=0$ |
| Point source | $\nabla\cdot(\hat{\mathbf{R}}/R^2)=4\pi\delta(\mathbf{R})$, $\nabla^2(1/R)=-4\pi\delta(\mathbf{R})$ |
| Maxwell (differential) | $\nabla\cdot\mathbf{D}=\rho$, $\nabla\cdot\mathbf{B}=0$, $\nabla\times\mathbf{E}=-\partial_t\mathbf{B}$, $\nabla\times\mathbf{H}=\mathbf{J}+\partial_t\mathbf{D}$ |
| Potentials | $\mathbf{B}=\nabla\times\mathbf{A}$, $\mathbf{E}=-\nabla V-\partial_t\mathbf{A}$; $\nabla^2V=-\rho/\varepsilon$ |
| Wave / Helmholtz | $\nabla^2\mathbf{E}=\mu\varepsilon\,\partial_t^2\mathbf{E}$; $\nabla^2\mathbf{E}+k^2\mathbf{E}=\mathbf{0}$, $k=\omega\sqrt{\mu\varepsilon}$ |
| Poynting | $\mathbf{S}=\mathbf{E}\times\mathbf{H}$; $\langle\mathbf{S}\rangle=\frac12\operatorname{Re}\{\mathbf{E}\times\mathbf{H}^*\}$ |

## Further reading

- E. Kreyszig, *Advanced Engineering Mathematics*, 10th ed., Wiley — Ch. 9–10: vector differential and integral calculus, with the integral theorems.
- H. M. Schey, *Div, Grad, Curl, and All That: An Informal Text on Vector Calculus*, 4th ed., Norton — builds the operators from flux and circulation, exactly as in Section 7.5.
- D. K. Cheng, *Field and Wave Electromagnetics*, 2nd ed., Addison-Wesley — Ch. 2 (vector analysis in the coordinate conventions used here) and Ch. 7 (Maxwell's equations).
- D. J. Griffiths, *Introduction to Electrodynamics*, 4th ed., Cambridge University Press — Ch. 1 is a compact, careful treatment of vector calculus, including the Dirac delta and the Helmholtz theorem.
