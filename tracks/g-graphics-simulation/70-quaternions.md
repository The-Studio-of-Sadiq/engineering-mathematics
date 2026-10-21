# Chapter 70 — Quaternions and the Geometry of Rotation

[← Previous](69-mathematics-of-3d-worlds.md) · [Track G](README.md) · [Table of Contents](../../00-front/table-of-contents.md) · [Next →](71-ray-tracing.md)

> **Level 2 · Track G — Mathematics for Computer Graphics, Games, and Simulation**
> **Builds on:** [Ch. 2](../../part-1-foundations/02-complex-numbers.md), [Ch. 9](../../part-3-linear-algebra/09-vector-spaces-and-linear-transformations.md), [Ch. 10](../../part-3-linear-algebra/10-eigenvalues-and-eigenvectors.md), [Ch. 69](69-mathematics-of-3d-worlds.md) · **Leads to:** [Ch. 71](71-ray-tracing.md), [Ch. 72](72-physics-engines.md) · **Related:** [Ch. 56](../e-abstract-math/56-abstract-algebra.md), [Ch. 58](../e-abstract-math/58-topology.md), [Ch. 59](../e-abstract-math/59-differential-geometry.md)
> **Reading time:** ~25 min

## 70.1 The hook

July 1969. Apollo 11 is on its way to the Moon, and the spacecraft's sense of direction lives in a gyroscope-stabilised platform held inside three nested rings, or **gimbals**. Each gimbal lets the spacecraft turn about one axis while the platform stays fixed relative to the stars. Three rings, three angles: yaw, pitch, roll. It sounds complete.

It is not. If the spacecraft turns so that the middle ring is at 90°, the outer and inner rings line up, and their two axes become the same axis. The platform has lost one of its three directions of freedom. A turn in the missing direction now drags the platform with it, and the guidance computer's picture of "which way is up" is gone, to be rebuilt by hand from star sightings. This is **gimbal lock**. The crew were trained to steer clear of it, and the command module pilot, Michael Collins, is often quoted as asking mission control, half joking, for a fourth gimbal for Christmas.

Fifty years later every 3-D game, drone, phone and VR headset faces the same problem in software. Store an orientation as three angles, and there is some orientation where the angles jam. Store it as a $3\times3$ rotation matrix, as in [Ch. 69](69-mathematics-of-3d-worlds.md), and you carry nine numbers for three degrees of freedom; they drift in floating point, and there is no honest way to say what orientation lies *halfway* between two others. A camera that must swing smoothly from "facing north" to "facing straight up" has no good recipe.

The cure was found 126 years before Apollo, on 16 October 1843, by a man walking along a canal in Dublin. William Rowan Hamilton had spent years trying to multiply triples of numbers the way complex numbers multiply pairs. In a flash he saw that triples would never work, but *quadruples* would, and he carved the rule into the stone of Broom Bridge:

$$
\mathbf{i}^2=\mathbf{j}^2=\mathbf{k}^2=\mathbf{i}\mathbf{j}\mathbf{k}=-1.
$$

Why should rotating a 3-D world need a 4-D number? And why does the answer involve *half* the rotation angle?

## 70.2 The big idea

Start with something you already know. In the plane, multiplying by the unit complex number $e^{j\theta}$ rotates every point by $\theta$ ([Ch. 2 §2.7](../../part-1-foundations/02-complex-numbers.md#27-complex-plane-and-geometric-interpretation)). Rotations of the plane *are* unit complex numbers; composing rotations is multiplying them; and the halfway point between two orientations is just the point halfway along the arc of the unit circle. One number, no jamming, perfect interpolation.

Quaternions do the same for space, with two twists.

**Rotations have axes.** Whatever tangle of turns you apply to a rigid body, the net result is a single turn by some angle $\theta$ about some fixed axis $\hat{\mathbf{n}}$ (Euler's rotation theorem, proved below). So an orientation is really "an axis and an angle": four numbers with one constraint, which is exactly what a unit quaternion is. The unit quaternions form a 3-D sphere sitting in 4-D space, and **every point on that sphere is an orientation**. There are no poles, no special angles, no place where the coordinates jam.

**The angle is halved, and every rotation appears twice.** A unit quaternion $q$ does not rotate a vector by multiplying it once; it *sandwiches* it, $\mathbf{v}\mapsto q\mathbf{v}q^{*}$. Each factor contributes half the turn, so $q$ must carry $\theta/2$. The price is that $q$ and $-q$ give the same rotation: the sphere of quaternions wraps twice around the space of rotations.

Interpolation then becomes geometry on that sphere. To swing from orientation $q_0$ to $q_1$, walk along the great circle between them at constant speed. That walk is **SLERP**, and it is the reason a game camera can glide where Euler angles would lurch.

The analogy: Euler angles are like describing a location on Earth by latitude and longitude. They work almost everywhere, but at the North Pole longitude is meaningless and a small step can change it by 180°. Quaternions are like describing the location by a point in 3-D space on the globe's surface: no special places at all.

## 70.3 The mathematics

### Rotations and their axes

A **rotation** of space is a $3\times3$ matrix $\mathbf{R}$ that preserves lengths and handedness: $\mathbf{R}^\mathsf{T}\mathbf{R}=\mathbf{I}$ and $\det\mathbf{R}=1$ (orthogonal matrices are [Ch. 9 §9.7](../../part-3-linear-algebra/09-vector-spaces-and-linear-transformations.md#97-orthogonality)). The set of all of them is called $SO(3)$, the *special orthogonal group*. It is a group: products and inverses of rotations are rotations, and $\mathbf{R}^{-1}=\mathbf{R}^\mathsf{T}$. It is not commutative: tip a book forward and then turn it left, versus left then forward, and the book ends up in two different poses.

> [!theorem] Result 1: Euler's rotation theorem
> Every rotation $\mathbf{R}\in SO(3)$ fixes some unit vector $\hat{\mathbf{n}}$: $\mathbf{R}\hat{\mathbf{n}}=\hat{\mathbf{n}}$. So every rotation is a turn about an axis.

*Proof.* We need $1$ to be an eigenvalue of $\mathbf{R}$ ([Ch. 10 §10.2](../../part-3-linear-algebra/10-eigenvalues-and-eigenvectors.md#102-eigenvalues-and-eigenvectors)), that is, $\det(\mathbf{R}-\mathbf{I})=0$. Using $\mathbf{R}\mathbf{R}^\mathsf{T}=\mathbf{I}$, $\det\mathbf{R}=1$ and $\det\mathbf{M}^\mathsf{T}=\det\mathbf{M}$:

$$
\det(\mathbf{R}-\mathbf{I})=\det\big(\mathbf{R}(\mathbf{I}-\mathbf{R}^\mathsf{T})\big)=\det\mathbf{R}\,\det\big((\mathbf{I}-\mathbf{R})^\mathsf{T}\big)=\det(\mathbf{I}-\mathbf{R})=(-1)^3\det(\mathbf{R}-\mathbf{I}).
$$

A number equal to its own negative is zero. The factor $(-1)^3$ is where three dimensions matter: in the plane, $(-1)^2=1$ and the argument (rightly) fails, since a planar rotation fixes no direction. $\square$

The axis $\hat{\mathbf{n}}$ is the eigenvector for eigenvalue 1. The other two eigenvalues are $e^{\pm j\theta}$, and $\theta$ is the angle of turn; since the trace is the sum of the eigenvalues, $\operatorname{tr}\mathbf{R}=1+2\cos\theta$. Conversely, the rotation by $\theta$ about $\hat{\mathbf{n}}$ is given by **Rodrigues' formula**,

$$
\mathbf{R}\mathbf{v}=\mathbf{v}\cos\theta+(\hat{\mathbf{n}}\times\mathbf{v})\sin\theta+\hat{\mathbf{n}}(\hat{\mathbf{n}}\cdot\mathbf{v})(1-\cos\theta),
$$

which you can read straight off the figure: split $\mathbf{v}$ into its part along the axis, which does not move, and its part perpendicular to the axis, which turns by $\theta$ in the plane spanned by $\mathbf{v}_\perp$ and $\hat{\mathbf{n}}\times\mathbf{v}$.

![A vertical axis n-hat through the origin. A blue vector v from the origin has its tip on a dashed circle around the axis; an orange arc along the circle carries it through angle theta to a red vector R v. Green arrows from the circle's centre on the axis show the perpendicular parts of v and R v; the centre itself is the parallel part.](../../images/ch70-axis-rotation.svg)

*Rotation about an axis. The component $\mathbf{v}_\parallel$ along $\hat{\mathbf{n}}$ stays put; the component $\mathbf{v}_\perp$ swings by $\theta$ around the circle. Every rotation of space looks like this, for some axis (Result 1).*

### Euler angles and gimbal lock

The most popular way to write a rotation is as three turns about fixed coordinate axes, for example

$$
\mathbf{R}(\alpha,\beta,\gamma)=\mathbf{R}_z(\alpha)\,\mathbf{R}_y(\beta)\,\mathbf{R}_x(\gamma),
$$

with yaw $\alpha$, pitch $\beta$ and roll $\gamma$. Three numbers that a human can read: an excellent user interface. As a data structure, it has a hole.

> [!theorem] Result 2: gimbal lock
> At pitch $\beta=90^\circ$, the rotation depends only on $\alpha-\gamma$:
> $$
> \mathbf{R}_z(\alpha)\,\mathbf{R}_y(90^\circ)\,\mathbf{R}_x(\gamma)=\mathbf{R}_z(\alpha-\gamma)\,\mathbf{R}_y(90^\circ).
> $$

*Proof.* For any rotation $\mathbf{Q}$, the conjugate $\mathbf{Q}\mathbf{R}_{\hat{\mathbf{a}}}(\gamma)\mathbf{Q}^{-1}$ is the rotation by $\gamma$ about the turned axis $\mathbf{Q}\hat{\mathbf{a}}$ (turn the axis into place, rotate, turn back). Take $\mathbf{Q}=\mathbf{R}_y(90^\circ)$, which sends $\hat{\mathbf{x}}$ to $-\hat{\mathbf{z}}$. Then $\mathbf{R}_y(90^\circ)\mathbf{R}_x(\gamma)\mathbf{R}_y(90^\circ)^{-1}$ is a turn by $\gamma$ about $-\hat{\mathbf{z}}$, which is $\mathbf{R}_z(-\gamma)$. Hence $\mathbf{R}_y(90^\circ)\mathbf{R}_x(\gamma)=\mathbf{R}_z(-\gamma)\mathbf{R}_y(90^\circ)$, and multiplying on the left by $\mathbf{R}_z(\alpha)$ gives the claim. $\square$

So at that pitch, changing yaw and changing roll do *the same thing*. The map from three angles to rotations, which is smooth and three-to-three almost everywhere, collapses to two dimensions there: its Jacobian ([Ch. 6](../../part-2-calculus/06-multivariable-calculus.md)) drops rank. In a mechanical gimbal, that is the moment two rings become coplanar. In software, it shows up as a camera that suddenly spins wildly when it tries to follow a target passing overhead: a small change of orientation demands a huge change of yaw and roll.

![Two panels of three nested rings with their rotation axes. Left, at pitch 20 degrees: yaw (blue, vertical axis), pitch (green) and roll (red) axes point in three different directions. Right, at pitch 90 degrees: the red roll axis lies along the blue yaw axis and the red ring is parallel to the blue ring.](../../images/ch70-gimbal-lock.svg)

*Gimbal lock. Each ring turns about its own axis. At $90^\circ$ of pitch the roll axis lands on the yaw axis, so two rings produce the same motion and no ring can produce the third.*

The hole is not a flaw of this particular choice. Any description of all rotations by three smooth angles has such a place somewhere; the reason is topological (see §70.7). The escape is to use more numbers, but cleverly.

### Quaternion algebra

> [!definition] Quaternions
> A **quaternion** is $q=w+x\mathbf{i}+y\mathbf{j}+z\mathbf{k}$ with real $w,x,y,z$, multiplied using Hamilton's rule $\mathbf{i}^2=\mathbf{j}^2=\mathbf{k}^2=\mathbf{i}\mathbf{j}\mathbf{k}=-1$, which implies $\mathbf{i}\mathbf{j}=\mathbf{k}=-\mathbf{j}\mathbf{i}$, $\mathbf{j}\mathbf{k}=\mathbf{i}=-\mathbf{k}\mathbf{j}$, $\mathbf{k}\mathbf{i}=\mathbf{j}=-\mathbf{i}\mathbf{k}$. Write $q=(w,\mathbf{v})$ with scalar part $w$ and vector part $\mathbf{v}=(x,y,z)$. The **conjugate** is $q^*=(w,-\mathbf{v})$ and the **norm** is $\lVert q\rVert=\sqrt{w^2+x^2+y^2+z^2}$.

(The book's imaginary unit is $j$; the quaternion units are written bold, $\mathbf{i},\mathbf{j},\mathbf{k}$, and each of them squares to $-1$ like $j$ does. The complex numbers sit inside the quaternions as $w+x\mathbf{i}$.)

Multiplying out with the table and collecting terms gives the whole algebra in one line:

$$
(w_1,\mathbf{v}_1)(w_2,\mathbf{v}_2)=\big(w_1w_2-\mathbf{v}_1\cdot\mathbf{v}_2,\ \ w_1\mathbf{v}_2+w_2\mathbf{v}_1+\mathbf{v}_1\times\mathbf{v}_2\big).
$$

Both 3-D products appear, the dot and the cross. Hamilton's algebra is in fact where they came from: physicists later split quaternion multiplication into these two pieces and called the result vector analysis. The cross product is why quaternions do not commute, matching the fact that rotations do not.

Three facts follow quickly. $qq^*=(\lVert q\rVert^2,\mathbf{0})$, so every non-zero quaternion has an inverse $q^{-1}=q^*/\lVert q\rVert^2$. Conjugation reverses products, $(pq)^*=q^*p^*$. And the norm is multiplicative, $\lVert pq\rVert=\lVert p\rVert\,\lVert q\rVert$, since $\lVert pq\rVert^2=pq(pq)^*=p\,qq^*\,p^*=\lVert q\rVert^2\,pp^*$. So **products of unit quaternions are unit quaternions**: the unit sphere $S^3$ in 4-D space is closed under multiplication, just as the unit circle is for complex numbers.

### The sandwich

Identify a 3-D vector $\mathbf{v}$ with the "pure" quaternion $(0,\mathbf{v})$.

> [!theorem] Result 3: unit quaternions rotate by sandwiching
> Let $q=\big(\cos\tfrac{\theta}{2},\ \sin\tfrac{\theta}{2}\,\hat{\mathbf{n}}\big)$ with $\hat{\mathbf{n}}$ a unit vector. Then $q\,\mathbf{v}\,q^*$ is a pure quaternion, and it equals the rotation of $\mathbf{v}$ by angle $\theta$ about $\hat{\mathbf{n}}$.

*Proof.* Write $c=\cos\tfrac{\theta}{2}$, $s=\sin\tfrac{\theta}{2}$, and split $\mathbf{v}=\mathbf{v}_\parallel+\mathbf{v}_\perp$ along and across $\hat{\mathbf{n}}$. The map $\mathbf{v}\mapsto q\mathbf{v}q^*$ is linear, so treat the two parts separately.

*Parallel part.* $\mathbf{v}_\parallel=a\hat{\mathbf{n}}$ is a multiple of $q$'s own vector part, so it commutes with $q$ (the cross product term vanishes). Hence $q\mathbf{v}_\parallel q^*=\mathbf{v}_\parallel qq^*=\mathbf{v}_\parallel$: the axis does not move.

*Perpendicular part.* With $\hat{\mathbf{n}}\cdot\mathbf{v}_\perp=0$, the product rule gives

$$
q\,\mathbf{v}_\perp=(c,s\hat{\mathbf{n}})(0,\mathbf{v}_\perp)=\big(-s\,\hat{\mathbf{n}}\cdot\mathbf{v}_\perp,\ c\,\mathbf{v}_\perp+s\,\hat{\mathbf{n}}\times\mathbf{v}_\perp\big)=\big(0,\ c\,\mathbf{v}_\perp+s\,\hat{\mathbf{n}}\times\mathbf{v}_\perp\big),
$$

and the same computation shows $\mathbf{v}_\perp q^*=(0,\ c\,\mathbf{v}_\perp+s\,\hat{\mathbf{n}}\times\mathbf{v}_\perp)=q\mathbf{v}_\perp$. So $q\mathbf{v}_\perp q^*=q\,(q\mathbf{v}_\perp)=q^2\mathbf{v}_\perp$. Now $q^2=(c^2-s^2,\ 2cs\,\hat{\mathbf{n}})=(\cos\theta,\ \sin\theta\,\hat{\mathbf{n}})$ by the double-angle formulas, so applying the first computation again with $\theta$ in place of $\theta/2$:

$$
q\,\mathbf{v}_\perp\,q^*=\big(0,\ \mathbf{v}_\perp\cos\theta+(\hat{\mathbf{n}}\times\mathbf{v}_\perp)\sin\theta\big).
$$

That is $\mathbf{v}_\perp$ turned by $\theta$ in its plane. Adding the unmoved $\mathbf{v}_\parallel$ gives exactly Rodrigues' formula. $\square$

Now the half angle is no mystery: one factor $q$ on each side, each turning the perpendicular part by $\theta/2$. A check: rotate $\mathbf{v}=\hat{\mathbf{x}}$ by $90^\circ$ about $\hat{\mathbf{z}}$. Here $q=(\tfrac{1}{\sqrt2},0,0,\tfrac{1}{\sqrt2})$, and two applications of the product rule give $q\hat{\mathbf{x}}q^*=(0,0,1,0)$, the vector $\hat{\mathbf{y}}$. Correct.

Because the sandwich respects products, $(pq)\mathbf{v}(pq)^*=p\,(q\mathbf{v}q^*)\,p^*$, so **composing rotations is multiplying quaternions**, with the rightmost acting first as for matrices. A turn of $90^\circ$ about $\hat{\mathbf{x}}$ followed by $90^\circ$ about $\hat{\mathbf{y}}$ is

$$
q_y\,q_x=\left(\tfrac{1}{\sqrt2},0,\tfrac{1}{\sqrt2},0\right)\left(\tfrac{1}{\sqrt2},\tfrac{1}{\sqrt2},0,0\right)=\left(\tfrac12,\tfrac12,\tfrac12,-\tfrac12\right),
$$

whose scalar part $\cos\tfrac{\theta}{2}=\tfrac12$ says $\theta=120^\circ$, about the axis $(1,1,-1)/\sqrt3$. Two quarter-turns about perpendicular axes make a third of a turn about a diagonal: try it with a book.

Finally, since $(-q)\mathbf{v}(-q)^*=q\mathbf{v}q^*$, the antipodal points $q$ and $-q$ on $S^3$ are the same rotation. Unit quaternions **double cover** $SO(3)$. Rotating by $360^\circ$ about any axis gives $q=(\cos180^\circ,\mathbf{0})=-1$: the identity rotation, but not the identity quaternion. You must turn by $720^\circ$ to get back to $q=1$.

### Interpolation: SLERP

Given orientations $q_0,q_1$ (unit quaternions), what is "the orientation a fraction $t$ of the way from $q_0$ to $q_1$"? The natural answer is the point a fraction $t$ of the way along the shortest great-circle arc between them on $S^3$. Let $\cos\Omega=q_0\cdot q_1$ (the 4-D dot product). The arc lies in the plane spanned by $q_0$ and $q_1$; in that plane, take the unit vector $q_\perp=(q_1-q_0\cos\Omega)/\sin\Omega$ orthogonal to $q_0$, so the arc is $q_0\cos(t\Omega)+q_\perp\sin(t\Omega)$. Substituting and using $\sin\Omega\cos t\Omega-\cos\Omega\sin t\Omega=\sin\big((1-t)\Omega\big)$ gives

$$
\operatorname{slerp}(q_0,q_1;t)=\frac{\sin\big((1-t)\Omega\big)}{\sin\Omega}\,q_0+\frac{\sin(t\Omega)}{\sin\Omega}\,q_1.
$$

It moves at constant angular speed, and the physical rotation turns at the constant rate $2\Omega$ per unit $t$ about a fixed axis. One detail matters: $q_1$ and $-q_1$ are the same orientation but different points, and the arc to one is longer than to the other. If $q_0\cdot q_1<0$, replace $q_1$ by $-q_1$ first, or the camera goes "the long way round".

The cheap alternative, **normalised lerp**, draws the straight chord $(1-t)q_0+tq_1$ and pushes each point back onto the sphere. It follows the same path but at the wrong speed: slow at the ends, fast in the middle.

![Left: a unit circle with endpoints q0 and q1 separated by angle Omega of 150 degrees. Nine evenly spaced grey points on the chord project to red points on the arc that crowd near the ends and spread out in the middle; blue squares drawn just outside the circle are evenly spaced. Right: angle travelled against t, a straight blue line for slerp and an S-shaped red curve for normalised lerp.](../../images/ch70-slerp.svg)

*SLERP against normalised lerp, on a circle (the 2-D slice of $S^3$ containing $q_0$ and $q_1$), with the angle exaggerated to $150^\circ$. Evenly spaced points on the chord project to unevenly spaced points on the arc; at this angle the projected point moves at 2.85 times the slerp speed in the middle and 0.19 times at the ends.*

## 70.4 The algorithm

The core toolkit fits on a page.

```text
qmul((w1,v1),(w2,v2)) = (w1*w2 - v1.v2,  w1*v2 + w2*v1 + v1 x v2)   # 16 mul, 12 add
from_axis_angle(n, th) = (cos(th/2), sin(th/2) * n)                  # n must be unit
rotate(q, v)          = vector part of qmul(qmul(q, (0,v)), conj(q))
to_matrix(q)          = 3x3 matrix with entries 1-2(y^2+z^2), 2(xy-wz), ...   # for bulk use
slerp(q0, q1, t):
    d = dot(q0, q1)
    if d < 0: q1, d = -q1, -d            # take the short way round
    if d > 0.9995: return normalize((1-t)*q0 + t*q1)   # nearly equal: avoid 0/0
    W = acos(d)
    return (sin((1-t)*W)*q0 + sin(t*W)*q1) / sin(W)
```

**Cost.** Composing two rotations costs 16 multiplications as quaternions against 27 as matrices, and a quaternion is 4 numbers against 9. But rotating *one vector* by the sandwich costs more than a matrix–vector product (9 multiplies). So engines store and compose orientations as quaternions, then convert once to a matrix per object per frame before transforming thousands of vertices, exactly as in the pipeline of [Ch. 69 §69.4](69-mathematics-of-3d-worlds.md#694-the-algorithm).

**Rigid bodies.** A full pose is a rotation and a position, $\mathbf{x}\mapsto q\mathbf{x}q^*+\mathbf{t}$. Poses compose like the block matrices of Ch. 69: $(q_2,\mathbf{t}_2)\circ(q_1,\mathbf{t}_1)=(q_2q_1,\ q_2\mathbf{t}_1q_2^*+\mathbf{t}_2)$, where again the second pose rotates the first one's translation. Skeletal animation stores each bone as such a pair relative to its parent and multiplies down the chain from hip to fingertip. (Blending poses by blending the pairs separately makes twisting joints collapse like a candy wrapper, which is why animators use *dual quaternions*, a variant that interpolates rotation and translation together.)

**Failure cases.**

- *Drift.* Repeated multiplication in floating point slowly moves $q$ off the unit sphere, and a non-unit $q$ also scales whatever it rotates. The fix is one line, $q\leftarrow q/\lVert q\rVert$. A drifting matrix needs a full Gram–Schmidt re-orthonormalisation instead, and you must decide which column to trust.
- *Near-equal endpoints.* In SLERP, $\sin\Omega\to0$ as $q_0\to q_1$; switch to normalised lerp, which is then accurate.
- *Sign flips.* Data from sensors or animation files may jump between $q$ and $-q$ from frame to frame. Both are right, but averaging or filtering them naively gives nonsense near zero. Always align signs (by the dot product test) before blending.
- *Conventions.* Some libraries store $(x,y,z,w)$, some $(w,x,y,z)$; some multiply in the opposite order. Most quaternion bugs are convention bugs.

## 70.5 The experiment

First, see how uneven normalised lerp is. For two quaternions at angle $\Omega$ apart, the angle reached by normalised lerp at parameter $t$ (here $x$) is $\arctan\!\big(t\sin\Omega/(1-t+t\cos\Omega)\big)$; slerp reaches $t\Omega$. The graph uses $\Omega=1.5$ rad, close to the largest angle you meet after the sign test (a rotation of about $172^\circ$).

```desmos-graph
left=0; right=1; bottom=-0.1; top=1.6
---
y=1.5x|blue
y=\arctan\left(\frac{x\sin(1.5)}{1-x+x\cos(1.5)}\right)|red
y=10\left(\arctan\left(\frac{x\sin(1.5)}{1-x+x\cos(1.5)}\right)-1.5x\right)+0.8|grey|dashed
```

*Angle travelled by slerp (blue) and normalised lerp (red) for $\Omega=1.5$; the dashed curve is their difference magnified ten times (centred at 0.8). Replace every `1.5` by `0.5`: the curves nearly coincide, which is why games often get away with lerp for small steps. Keep $\Omega$ below about 1.57 so the denominator stays positive.*

At $\Omega=1.5$ the red curve starts at 0.665 times the slerp speed and reaches 1.24 times it in the middle; the largest angle error is about $3.5^\circ$, visible as a stutter in a slow camera move.

Second, watch numbers drift. The script below spins a body at constant angular velocity for $10^5$ small steps, once by multiplying rotation matrices and once by multiplying quaternions, and measures how far each has wandered from being a true rotation.

```python
import numpy as np
import matplotlib.pyplot as plt

DTYPE = np.float32               # KNOB: float32 (a game) or float64 (a simulator)
STEPS = 100_000                  # KNOB: number of frames
OMEGA = np.array([0.3, 1.0, 0.2])  # KNOB: angular velocity (rad/s)
DT = 1e-3                        # KNOB: time step (s)
RENORM_EVERY = 0                 # KNOB: renormalise the quaternion every N steps (0 = never)

def qmul(a, b):                  # Hamilton product, (w, x, y, z)
    w1, v1, w2, v2 = a[0], a[1:], b[0], b[1:]
    return np.r_[w1 * w2 - v1 @ v2, w1 * v2 + w2 * v1 + np.cross(v1, v2)].astype(DTYPE)

def qmat(q):                     # unit quaternion -> rotation matrix
    w, x, y, z = q
    return np.array([[1 - 2*(y*y + z*z), 2*(x*y - w*z), 2*(x*z + w*y)],
                     [2*(x*y + w*z), 1 - 2*(x*x + z*z), 2*(y*z - w*x)],
                     [2*(x*z - w*y), 2*(y*z + w*x), 1 - 2*(x*x + y*y)]], dtype=DTYPE)

th = np.linalg.norm(OMEGA) * DT  # one small step: angle th about OMEGA's direction
dq = np.r_[np.cos(th / 2), np.sin(th / 2) * OMEGA / np.linalg.norm(OMEGA)].astype(DTYPE)
dR = qmat(dq)
q, R = np.array([1, 0, 0, 0], DTYPE), np.eye(3, dtype=DTYPE)
log_t, err_R, err_q = [], [], []
for k in range(1, STEPS + 1):
    R = (R @ dR).astype(DTYPE)   # matrix: 27 multiplies per step
    q = qmul(q, dq)              # quaternion: 16 multiplies per step
    if RENORM_EVERY and k % RENORM_EVERY == 0:
        q = (q / np.linalg.norm(q)).astype(DTYPE)
    if k % 500 == 0:
        log_t.append(k * DT)
        err_R.append(np.abs(R.T @ R - np.eye(3)).max())     # how far from a rotation
        err_q.append(abs(float(np.linalg.norm(q)) - 1.0))  # how far from unit length
print(f"after {STEPS} steps: matrix error {err_R[-1]:.2e}, quaternion error {err_q[-1]:.2e}")
plt.semilogy(log_t, err_R, label="matrix: max |R^T R - I|")
plt.semilogy(log_t, err_q, label="quaternion: | |q| - 1 |")
plt.xlabel("simulated time (s)"); plt.legend(); plt.show()
```

With the defaults (single precision, no renormalisation) both drift by about $10^{-3}$ after $10^5$ steps: the matrix error is $8.6\times10^{-4}$ and the quaternion error $1.8\times10^{-3}$. Quaternions are *not* magically immune to rounding. The difference is what the repair costs. Now:

- Set `RENORM_EVERY = 100`. One square root and four divisions every hundred frames, and the quaternion error is gone, while the matrix keeps shearing.
- Set `DTYPE = np.float64`. Both errors drop to about $10^{-11}$. Simulators can afford this; GPUs and phone sensors often cannot.
- Raise `DT` to `0.05` and compare the final orientations of the two methods. They agree with each other, because each step is an exact rotation. But is each step the *right* rotation when $\boldsymbol\omega$ itself changes over time? That question belongs to physics engines.

## 70.6 The challenge

1. **Intuition.** Hold a cup of water flat on your palm and rotate your hand a full $360^\circ$ about the vertical, keeping the cup level (your arm will end up twisted). Keep going another $360^\circ$ in the same direction. What happens to your arm? Relate this to the fact that the rotation by $360^\circ$ is the quaternion $-1$, not $+1$.
2. **Derivation.** Show that $\operatorname{tr}\mathbf{R}=1+2\cos\theta$ for the rotation matrix built from $q=(\cos\tfrac{\theta}{2},\sin\tfrac{\theta}{2}\hat{\mathbf{n}})$, and use it to derive a formula recovering $q$ from $\mathbf{R}$ when $w=\cos\tfrac{\theta}{2}$ is not small. Which entries of $\mathbf{R}$ should you use instead when $\theta$ is close to $180^\circ$, and why does the naive formula then lose accuracy?
3. **Implementation.** Write a function that converts Euler angles $(\alpha,\beta,\gamma)$ to a quaternion and back. Feed it a camera that tracks an aircraft flying straight over your head, and plot the recovered yaw and roll against time. Watch them blow up near pitch $90^\circ$, then animate the same flight with SLERP between keyframe quaternions and confirm it stays smooth.
4. **Open investigation.** You have $N$ noisy measurements of the same orientation, as unit quaternions (some with flipped signs). What is their "average"? Try: (a) align signs and normalise the plain sum; (b) take the eigenvector for the largest eigenvalue of $\sum_i q_iq_i^\mathsf{T}$ (a $4\times4$ matrix); (c) iterate "average in the tangent space, then step along the sphere". Compare them as the noise grows from $1^\circ$ to $90^\circ$. When do they disagree, and which one would you trust in a drone?

## 70.7 Loose threads

Apollo's problem has an answer. Three angles always have a place where they jam (Result 2 shows where for one convention). Rotation matrices never jam but carry nine numbers. Unit quaternions carry four, compose with 16 multiplications, repair their drift with one division, and turn interpolation into walking along a great circle. Rotations have axes (Result 1), and the half-angle is the signature of the sandwich (Result 3).

What we did not tell you:

> [!mystery] Why the sphere wraps twice
> Every rotation appears twice on the quaternion sphere, as $q$ and $-q$. That is not a defect of quaternions but a fact about rotations themselves: a full $360^\circ$ turn traces a loop in $SO(3)$ that cannot be shrunk to nothing, yet going round it *twice* can be (the cup-and-arm trick of challenge 1 is a physical demonstration). This is also why no three smooth angles can describe every rotation without a jam. The tools that make "loops that cannot be shrunk" precise are in [Ch. 58](../e-abstract-math/58-topology.md).

> [!mystery] Why four, and not three?
> Hamilton spent a decade trying to multiply triples and failed. It was not lack of imagination: a number system in which you can multiply, divide and keep the norm multiplicative exists only in dimensions 1, 2, 4 and 8 (the reals, complexes, quaternions and the strange non-associative *octonions*). Why exactly those four dimensions? The algebraic side of the story starts in [Ch. 56](../e-abstract-math/56-abstract-algebra.md).

> [!mystery] Rotations that change in time
> A spinning body has an angular velocity $\boldsymbol\omega(t)$, and its orientation obeys $\dot q=\tfrac12(0,\boldsymbol\omega)\,q$. Integrate that with a textbook method and $q$ leaves the unit sphere; renormalise and you may still get the wrong spin, and a tumbling brick in a game can gain energy from nowhere. Which integrators respect the sphere and the physics is a question for [Ch. 72](72-physics-engines.md).

- **Lie groups.** $SO(3)$ and $S^3$ are curved spaces that are also groups. The "tangent space average" of challenge 4 and the exponential map $\boldsymbol\omega\mapsto e^{(0,\boldsymbol\omega)/2}$ are the first steps of Lie theory, the geometry of continuous symmetry; see [Ch. 59](../e-abstract-math/59-differential-geometry.md).
- **Smooth paths through many keys.** SLERP joins two orientations, but a path through five keyframes made of four SLERPs has kinks at the joins. Spherical splines (such as SQUAD) smooth them; making "curvature" meaningful on a sphere is again differential geometry.
- **Light.** We have placed objects and cameras in the world and turned them smoothly. Now the world needs illumination: [Ch. 71](71-ray-tracing.md) follows rays of light, and rotations reappear every time a ray reflects.

## Further reading

- J. B. Kuipers, *Quaternions and Rotation Sequences*, Princeton University Press, 1999. A patient, engineering-flavoured treatment, including aerospace sequences and gimbal lock.
- A. J. Hanson, *Visualizing Quaternions*, Morgan Kaufmann, 2006. The geometry of $S^3$ and the double cover, with many pictures.
- K. Shoemake, "Animating rotation with quaternion curves", *SIGGRAPH* 1985. The paper that brought SLERP to computer graphics.
- J. Stillwell, *Naive Lie Theory*, Springer, 2008. From quaternions and $SO(3)$ to Lie groups, using only undergraduate linear algebra.
