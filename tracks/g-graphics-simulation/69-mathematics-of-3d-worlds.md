# Chapter 69 — The Mathematics of 3D Worlds

[← Previous](README.md) · [Track G](README.md) · [Table of Contents](../../00-front/table-of-contents.md) · [Next →](70-quaternions.md)

> **Level 2 · Track G — Mathematics for Computer Graphics, Games, and Simulation**
> **Builds on:** [Ch. 8](../../part-3-linear-algebra/08-matrices-and-linear-systems.md), [Ch. 9](../../part-3-linear-algebra/09-vector-spaces-and-linear-transformations.md), [Ch. 7](../../part-2-calculus/07-vector-calculus.md) · **Leads to:** [Ch. 70](70-quaternions.md), [Ch. 71](71-ray-tracing.md), [Ch. 75](75-computational-geometry.md) · **Related:** [Ch. 22](../../part-7-numerical-computational/22-numerical-methods.md), [Ch. 58](../e-abstract-math/58-topology.md), [Ch. 59](../e-abstract-math/59-differential-geometry.md)
> **Reading time:** ~25 min

## 69.1 The hook

A first-year game programmer is building a toy solar system. She wants the Earth to spin on its axis and also sit 10 units from the Sun. She writes two lines: "rotate the Earth by the day angle" and "move the Earth 10 units out". It runs, and the Earth does something nobody asked for: instead of spinning in place, it swings around the Sun once a day like a stone on a string. She swaps the two lines, and now it spins in place. Same two operations, same numbers, opposite worlds.

Meanwhile her graphics card is doing something that should be impossible. Every frame, sixty times a second, it takes several million corners of triangles, each given as three numbers in the coordinates of its own little model, and lands each one on a pixel of your screen. Along the way every point must be

- **rotated** (the model faces some direction),
- **translated** (the model sits somewhere in the world, the camera sits somewhere else),
- **scaled** (the model was built in centimetres, the world is in metres),
- and **projected**, which means dividing by its distance from the eye, so that far things look small.

Here is the puzzle. A matrix is a *linear* map (see [Ch. 9 §9.5](../../part-3-linear-algebra/09-vector-spaces-and-linear-transformations.md#95-linear-transformations)): it always sends the zero vector to zero. A translation moves the origin, so no $3\times3$ matrix can translate anything. And projection divides by depth, which is not even close to linear. Yet the hardware does all of it with one $4\times4$ matrix per object and a single division per vertex.

How can a matrix with an extra row and column carry a whole world onto a flat screen? And why did swapping two lines move the Earth?

## 69.2 The big idea

Three ideas do almost all of the work.

**Every coordinate is relative to a frame.** "The cup is at $(2, 0, 1)$" means nothing until you say *measured from where, along which directions*. A frame is an origin plus three axes. A game has many: each model's own frame (where the artist built it), the world frame, the camera frame, and finally the screen. Rendering is nothing but a chain of changes of frame. The solar-system bug is a question of which frame each operation happens in: rotating *after* moving the Earth out swings it about the Sun's origin; rotating *before* spins it about its own centre.

**Add one coordinate, and translation becomes linear.** Write a 3-D point as $(x, y, z, 1)$. In four dimensions, the 3-D world is a flat slice at height $w=1$, and a translation of that slice is a *shear* of 4-D space, which is linear. Directions get $w=0$, so the shear leaves them alone, exactly as it should: moving a ship does not change which way its nose points.

**A camera is a box with a hole.** Light from a point passes through a pinhole and hits a screen. By similar triangles the image height is the true height divided by the distance. With the fourth coordinate in hand, we can make a matrix *copy* the depth into $w$; the hardware's final step "divide by $w$" then performs the projection.

![Left: a 3-D scene drawn by perspective projection, showing a blue cube on a ground grid with world axes at the origin, a camera with its own three axes to the right, an orange image plane between them, and grey rays from each cube corner to the camera crossing the image plane where a small red copy of the cube appears. Right: the red cube as the camera sees it in screen coordinates.](../../images/ch69-camera-frames.svg)

*Three frames in one picture. The world frame sits at the origin; the camera frame ($x_c$ right, $y_c$ up, $z_c$ along the line of sight) sits at the pinhole. Each ray from a cube corner to the pinhole pierces the image plane, and the pierce points (red) form the picture on the right. The left panel was itself drawn by the same projection, from a second camera.*

The analogy to hold onto: the GPU is a factory with one machine, "multiply by a $4\times4$ matrix", and one finishing step, "divide by the last number". Every trick in this chapter is a way of disguising a non-linear wish as a matrix.

## 69.3 The mathematics

### Points, vectors and frames

There are two kinds of 3-D objects, and confusing them causes real bugs. A **point** is a location. A **vector** is a displacement: the arrow from one point to another. Their algebra is different: point minus point is a vector, point plus vector is a point, vector plus vector is a vector, but "point plus point" means nothing (what is Paris plus Tokyo?).

> [!definition] Frame
> A **frame** $(\mathbf{o};\ \mathbf{e}_1,\mathbf{e}_2,\mathbf{e}_3)$ is an origin point $\mathbf{o}$ and a basis $\mathbf{e}_1,\mathbf{e}_2,\mathbf{e}_3$ of vectors. The coordinates of a point $\mathbf{p}$ in this frame are the unique numbers $(x,y,z)$ with $\mathbf{p}=\mathbf{o}+x\mathbf{e}_1+y\mathbf{e}_2+z\mathbf{e}_3$. The coordinates of a vector $\mathbf{v}$ are the numbers with $\mathbf{v}=x\mathbf{e}_1+y\mathbf{e}_2+z\mathbf{e}_3$; the origin plays no part.

When the basis is orthonormal, coordinates are just dot products: $x=(\mathbf{p}-\mathbf{o})\cdot\mathbf{e}_1$ and so on. For a camera at $\mathbf{c}$ with orthonormal axes $\hat{\mathbf{r}}$ (right), $\hat{\mathbf{u}}$ (up), $\hat{\mathbf{f}}$ (forward), a world point has camera coordinates

$$
x_c=(\mathbf{p}-\mathbf{c})\cdot\hat{\mathbf{r}},\qquad y_c=(\mathbf{p}-\mathbf{c})\cdot\hat{\mathbf{u}},\qquad z_c=(\mathbf{p}-\mathbf{c})\cdot\hat{\mathbf{f}}.
$$

This is a change of basis followed by a shift of origin; the change-of-basis part is [Ch. 9 §9.4](../../part-3-linear-algebra/09-vector-spaces-and-linear-transformations.md#94-basis-and-dimension), and the fact that the inverse of an orthogonal matrix is its transpose is [Ch. 9 §9.7](../../part-3-linear-algebra/09-vector-spaces-and-linear-transformations.md#97-orthogonality). How do we build $\hat{\mathbf{r}},\hat{\mathbf{u}},\hat{\mathbf{f}}$ from "camera at $\mathbf{c}$, looking at $\mathbf{t}$, world up is $\hat{\mathbf{z}}$"? Normalise $\mathbf{t}-\mathbf{c}$ to get $\hat{\mathbf{f}}$; then $\hat{\mathbf{r}}=\hat{\mathbf{f}}\times\hat{\mathbf{z}}/\lVert\hat{\mathbf{f}}\times\hat{\mathbf{z}}\rVert$ and $\hat{\mathbf{u}}=\hat{\mathbf{r}}\times\hat{\mathbf{f}}$. Two cross products ([Ch. 7 §7.2](../../part-2-calculus/07-vector-calculus.md#72-vector-operations)) and you have a camera. (Notice the failure case already: if you look straight up, $\hat{\mathbf{f}}\times\hat{\mathbf{z}}=\mathbf{0}$ and the camera has no idea which way is right.)

### Affine maps

> [!definition] Affine map
> A map $T:\mathbb{R}^3\to\mathbb{R}^3$ is **affine** if $T(\mathbf{x})=\mathbf{A}\mathbf{x}+\mathbf{t}$ for some $3\times3$ matrix $\mathbf{A}$ and vector $\mathbf{t}$. Rotations, scalings, shears, reflections and translations are all affine. It is **rigid** if $\mathbf{A}$ is a rotation ($\mathbf{A}^\mathsf{T}\mathbf{A}=\mathbf{I}$, $\det\mathbf{A}=1$).

Affine maps are not linear (unless $\mathbf{t}=\mathbf{0}$), but they preserve something almost as good. Call $\sum_i\lambda_i\mathbf{p}_i$ an **affine combination** of points if the weights sum to one, $\sum_i\lambda_i=1$. Only these combinations make sense for points: $\tfrac12\mathbf{p}+\tfrac12\mathbf{q}$ is the midpoint, wherever the origin is, while $\mathbf{p}+\mathbf{q}$ depends on the arbitrary choice of origin.

> [!theorem] Result 1: affine maps preserve affine combinations
> If $T(\mathbf{x})=\mathbf{A}\mathbf{x}+\mathbf{t}$ and $\sum_i\lambda_i=1$, then $T\left(\sum_i\lambda_i\mathbf{p}_i\right)=\sum_i\lambda_i\,T(\mathbf{p}_i)$.

*Proof.* Using linearity of $\mathbf{A}$ and then $\sum_i\lambda_i=1$ to write $\mathbf{t}=\sum_i\lambda_i\mathbf{t}$:

$$
T\Big(\sum_i\lambda_i\mathbf{p}_i\Big)=\sum_i\lambda_i\mathbf{A}\mathbf{p}_i+\mathbf{t}=\sum_i\lambda_i\mathbf{A}\mathbf{p}_i+\sum_i\lambda_i\mathbf{t}=\sum_i\lambda_i\left(\mathbf{A}\mathbf{p}_i+\mathbf{t}\right)=\sum_i\lambda_i T(\mathbf{p}_i).
$$

The step that used $\sum\lambda_i=1$ is exactly where a non-affine combination would leave a stray $\mathbf{t}$. $\square$

The consequences are what make 3-D graphics cheap. The segment $\{(1-s)\mathbf{p}+s\mathbf{q}\}$ maps to the segment between $T(\mathbf{p})$ and $T(\mathbf{q})$, so **lines map to lines** and you only ever transform the endpoints. Midpoints map to midpoints, parallel lines stay parallel, and a point inside a triangle with barycentric coordinates $(\lambda_1,\lambda_2,\lambda_3)$ keeps them. A whole triangle with a million interior points is described by three corners. Remember this; projection will break it.

### Homogeneous coordinates

Now the trick from §69.2. Attach a fourth coordinate: a point $\mathbf{p}$ becomes $(\mathbf{p},1)$ and a vector $\mathbf{v}$ becomes $(\mathbf{v},0)$. Then

$$
\begin{bmatrix}\mathbf{A}&\mathbf{t}\\ \mathbf{0}^\mathsf{T}&1\end{bmatrix}\begin{bmatrix}\mathbf{p}\\1\end{bmatrix}=\begin{bmatrix}\mathbf{A}\mathbf{p}+\mathbf{t}\\1\end{bmatrix},\qquad
\begin{bmatrix}\mathbf{A}&\mathbf{t}\\ \mathbf{0}^\mathsf{T}&1\end{bmatrix}\begin{bmatrix}\mathbf{v}\\0\end{bmatrix}=\begin{bmatrix}\mathbf{A}\mathbf{v}\\0\end{bmatrix}.
$$

A $4\times4$ matrix in this block form performs an affine map on points and its linear part on vectors, with no special cases. The bookkeeping rules for points and vectors fall out of the last coordinate: point minus point has $w=1-1=0$ (a vector), point plus vector has $w=1$ (a point), and point plus point would have $w=2$, which is not on the slice $w=1$, the algebra's way of saying "meaningless".

Composition is now matrix multiplication ([Ch. 8 §8.2](../../part-3-linear-algebra/08-matrices-and-linear-systems.md#82-matrix-multiplication)). Applying $T_1(\mathbf{x})=\mathbf{A}_1\mathbf{x}+\mathbf{t}_1$ and then $T_2$ gives

$$
\begin{bmatrix}\mathbf{A}_2&\mathbf{t}_2\\ \mathbf{0}^\mathsf{T}&1\end{bmatrix}\begin{bmatrix}\mathbf{A}_1&\mathbf{t}_1\\ \mathbf{0}^\mathsf{T}&1\end{bmatrix}=\begin{bmatrix}\mathbf{A}_2\mathbf{A}_1&\mathbf{A}_2\mathbf{t}_1+\mathbf{t}_2\\ \mathbf{0}^\mathsf{T}&1\end{bmatrix}.
$$

Look at the top-right block: the second map *rotates the first map's translation*. That is the solar-system bug, in one formula. Let $\mathbf{R}$ be a rotation by $90^\circ$ about $z$ and $\mathbf{T}$ a translation by $(10,0,0)$. Apply both to the Earth's centre, which starts at the origin:

- **rotate, then translate** ($\mathbf{T}\mathbf{R}$): the centre stays at the origin under $\mathbf{R}$, then moves to $(10,0,0)$. The Earth spins in place, out where it belongs.
- **translate, then rotate** ($\mathbf{R}\mathbf{T}$): the centre moves to $(10,0,0)$, then the rotation about the Sun's origin carries it to $(0,10,0)$. A quarter-day later, the Earth is a quarter of the way around its orbit.

Our programmer's first version applied the day rotation last. The rule of thumb, "matrices act right to left, the rightmost acts first, in the object's own frame", saves more debugging hours in graphics than any other sentence.

Inverting a rigid transform never needs a general inverse. If $\mathbf{x}'=\mathbf{R}\mathbf{x}+\mathbf{t}$, then $\mathbf{x}=\mathbf{R}^\mathsf{T}(\mathbf{x}'-\mathbf{t})$, so

$$
\begin{bmatrix}\mathbf{R}&\mathbf{t}\\ \mathbf{0}^\mathsf{T}&1\end{bmatrix}^{-1}=\begin{bmatrix}\mathbf{R}^\mathsf{T}&-\mathbf{R}^\mathsf{T}\mathbf{t}\\ \mathbf{0}^\mathsf{T}&1\end{bmatrix}.
$$

The **view matrix** of a camera is exactly this: the camera's own placement in the world is a rigid transform $(\mathbf{R}_c,\mathbf{c})$ whose rotation has columns $\hat{\mathbf{r}},\hat{\mathbf{u}},\hat{\mathbf{f}}$; the view matrix is its inverse, with rows $\hat{\mathbf{r}}^\mathsf{T},\hat{\mathbf{u}}^\mathsf{T},\hat{\mathbf{f}}^\mathsf{T}$, which reproduces the dot-product formulas for $x_c,y_c,z_c$ above. Moving the camera left is the same as moving the world right.

### Normals do not transform like points

Here is a trap that catches professionals. A surface's **normal** $\mathbf{n}$ is a vector perpendicular to it, used for lighting. Since normals are vectors, surely they transform by $\mathbf{A}$? Try a non-uniform stretch $\mathbf{A}=\operatorname{diag}(2,1,1)$ on the plane $x+y=0$, which has tangent $\mathbf{s}=(1,-1,0)$ and normal $\mathbf{n}=(1,1,0)$. After stretching, the tangent is $\mathbf{A}\mathbf{s}=(2,-1,0)$, and the naive normal $\mathbf{A}\mathbf{n}=(2,1,0)$ has dot product $4-1=3$ with it. The lighting is now computed with a normal that is not normal.

> [!theorem] Result 2: normals transform by the inverse transpose
> If tangent vectors transform by an invertible $\mathbf{A}$, then $\mathbf{n}'=\mathbf{A}^{-\mathsf{T}}\mathbf{n}$ is perpendicular to every transformed tangent.

*Proof.* Perpendicularity is $\mathbf{n}^\mathsf{T}\mathbf{s}=0$ for every tangent $\mathbf{s}$. Insert $\mathbf{I}=\mathbf{A}^{-1}\mathbf{A}$:

$$
0=\mathbf{n}^\mathsf{T}\mathbf{s}=\mathbf{n}^\mathsf{T}\mathbf{A}^{-1}\mathbf{A}\mathbf{s}=\left(\mathbf{A}^{-\mathsf{T}}\mathbf{n}\right)^\mathsf{T}(\mathbf{A}\mathbf{s}).
$$

So $\mathbf{A}^{-\mathsf{T}}\mathbf{n}$ is perpendicular to every transformed tangent $\mathbf{A}\mathbf{s}$. $\square$

In the example, $\mathbf{A}^{-\mathsf{T}}\mathbf{n}=(0.5,1,0)$ and $(0.5,1,0)\cdot(2,-1,0)=1-1=0$. For a rotation, $\mathbf{R}^{-\mathsf{T}}=\mathbf{R}$, which is why the bug hides until someone scales a model unevenly. The deeper lesson: tangents and normals are different kinds of object. Tangents are vectors; normals are really *covectors*, linear functions that eat vectors, and covectors transform "the other way".

### The pinhole and the projection matrix

![Side view of a pinhole camera. A green arrow of height Y stands at depth Z on the right. A blue ray from its tip passes through the pinhole at the origin. In front of the pinhole, at distance f, an orange virtual image plane is crossed by the ray at height fY/Z; behind the pinhole, the real film plane receives an inverted image at height minus fY/Z. Shaded similar triangles are marked.](../../images/ch69-pinhole.svg)

*The pinhole camera. The small shaded triangle (pinhole, $f$, $y_s$) and the large one (pinhole, $Z$, $Y$) are similar, so $y_s=fY/Z$. A real camera's film sits behind the hole and sees the world upside down; graphics puts a virtual screen in front, which gives the same picture the right way up.*

In camera coordinates (looking along $+z_c$), the similar triangles give the **perspective projection**

$$
x_s=\frac{f\,x_c}{z_c},\qquad y_s=\frac{f\,y_c}{z_c}.
$$

A point at $(2,1,4)$ with $f=1$ lands at $(0.5,0.25)$; move it twice as far, to $(4,2,8)$, and it lands at the same place: every point on a ray through the pinhole has the same image. That is the geometric content of projection, and it is the same as the defining property of homogeneous coordinates: $(x,y,z,w)$ and $(\lambda x,\lambda y,\lambda z,\lambda w)$ for $\lambda\ne0$ name the same point $(x/w,y/w,z/w)$. So we let a matrix put $z_c$ into the last slot and let the divide do the rest:

$$
\begin{bmatrix}f&0&0&0\\0&f&0&0\\0&0&a&b\\0&0&1&0\end{bmatrix}\begin{bmatrix}x_c\\y_c\\z_c\\1\end{bmatrix}=\begin{bmatrix}f x_c\\f y_c\\a z_c+b\\z_c\end{bmatrix}\ \xrightarrow{\ \div w\ }\ \left(\frac{f x_c}{z_c},\ \frac{f y_c}{z_c},\ a+\frac{b}{z_c}\right).
$$

The third row keeps a **depth** value for deciding which surface is in front. Choosing $a$ and $b$ so that the near plane $z_c=n$ maps to depth 0 and the far plane $z_c=F$ to depth 1 gives $a=F/(F-n)$ and $b=-nF/(F-n)$. (Conventions differ: OpenGL looks down $-z$ and maps depth to $[-1,1]$; the mathematics is the same.)

The matrix made the projection look linear, but the divide is not. What survives?

> [!theorem] Result 3: perspective maps lines to lines, but not midpoints to midpoints
> For points in front of the camera ($z_c>0$), the projection of a segment is a segment. Ratios along the segment are not preserved.

*Proof.* Let the 3-D segment be $\mathbf{p}(s)=(1-s)\mathbf{p}_0+s\mathbf{p}_1$ with depths $z_0,z_1>0$. Its screen image is $\boldsymbol\sigma(s)=f\,(x(s),y(s))/z(s)$. Write $\boldsymbol\sigma_0,\boldsymbol\sigma_1$ for the endpoint images, so $f(x_i,y_i)=z_i\boldsymbol\sigma_i$. Then

$$
\boldsymbol\sigma(s)=\frac{(1-s)z_0\,\boldsymbol\sigma_0+s\,z_1\,\boldsymbol\sigma_1}{(1-s)z_0+s z_1}=(1-u)\,\boldsymbol\sigma_0+u\,\boldsymbol\sigma_1,\qquad u=\frac{s z_1}{(1-s)z_0+s z_1}.
$$

As $s$ runs from 0 to 1, $u$ runs monotonically from 0 to 1 (both weights stay positive), so the image is exactly the screen segment from $\boldsymbol\sigma_0$ to $\boldsymbol\sigma_1$. But $u\ne s$ unless $z_0=z_1$. $\square$

Take $\mathbf{p}_0=(-1,0,2)$ and $\mathbf{p}_1=(1,0,6)$ with $f=1$. The endpoints land at $x_s=-0.5$ and $x_s=1/6$. The 3-D midpoint $(0,0,4)$ lands at $x_s=0$, but the screen midpoint is $-1/6$. The near half of the segment takes up three quarters of its screen length. This is foreshortening, and it is why a texture painted on a floor must not be interpolated linearly across the screen: some mid-1990s game consoles did exactly that, and their floor tiles visibly wobble and swim as you walk.

## 69.4 The algorithm

Every rasterising renderer, from a phone game to a film previsualisation tool, runs this pipeline on every vertex:

```text
for each object:
    M  = model matrix          (object frame  -> world frame; scale, rotate, translate)
    V  = view matrix           (world frame   -> camera frame; inverse of camera placement)
    P  = projection matrix     (camera frame  -> clip space; copies depth into w)
    MVP = P @ V @ M            (compose ONCE per object, 64 multiplies)
    for each vertex v = (x, y, z, 1):
        c = MVP @ v            (16 multiplies, 12 adds)
        clip triangles against the near plane (where w would approach 0)
        ndc = c.xyz / c.w      (the perspective divide)
        pixel = viewport(ndc)  (scale and shift to the window)
    normals: transform by inverse-transpose of the upper-left 3x3 of V @ M
```

**Complexity.** Composing three matrices costs a few hundred flops once per object; each vertex then costs one matrix–vector product. Precomposing is the whole point: the alternative, applying model, view and projection separately, triples the per-vertex work. With $10^6$ vertices at 60 frames per second, that is the difference between $10^9$ and $3\times10^9$ flops per second for this step alone.

**Failure cases.**

- *Points behind the eye.* For $z_c<0$ the divide flips the sign, and a vertex behind you appears in front, upside down. A triangle with one corner behind the camera would be smeared across the whole screen. This is why triangles are **clipped** against the near plane *before* the divide, in homogeneous coordinates where nothing has blown up yet.
- *Degenerate look-at.* Looking straight along "world up" makes $\hat{\mathbf{f}}\times\hat{\mathbf{z}}=\mathbf{0}$. Flight and space games cannot use a fixed world-up vector at all; they need a better description of orientation (Ch. 70).
- *Matrix drift.* An object rotated a little every frame, $\mathbf{R}\leftarrow\mathbf{R}\,\Delta\mathbf{R}$, slowly stops being a rotation in floating point: its columns lose length and orthogonality, and the model visibly shears. The repair (re-orthonormalising, Gram–Schmidt from [Ch. 9 §9.7](../../part-3-linear-algebra/09-vector-spaces-and-linear-transformations.md#97-orthogonality)) is a patch; Ch. 70 offers a cure.
- *Depth precision.* The stored depth $a+b/z_c$ is a function of $1/z_c$, not $z_c$. With $n=0.1$ and $F=1000$, depth already reaches $0.50005$ at $z_c=0.2$: half the whole depth range is spent on the first 10 cm in front of your nose. With a 24-bit depth buffer, two surfaces near $z_c=1000$ must be about 0.6 units apart to be told apart (at $z_c=100$, about 0.006). Closer than that, they flicker through each other, a defect called **z-fighting**.

> [!mystery] The clue we're holding back
> Pushing the near plane out helps z-fighting, but cuts off things close to the eye. Modern engines do something stranger: they flip the depth mapping so the *far* plane is 0 and the near plane is 1, store it as a floating-point number, and the precision problem nearly vanishes. Why does reversing a number line fix it? The answer lies in how floating-point numbers are spaced, in [Ch. 22 §22.1](../../part-7-numerical-computational/22-numerical-methods.md#221-floating-point-arithmetic).

## 69.5 The experiment

First, see foreshortening as a graph. A point slides along the 3-D segment from $(-1,0,2)$ to $(1,0,6)$; the horizontal axis is the 3-D parameter $s$ (here called $x$) and the curve is its screen position $x_s=(-1+2s)/(2+4s)$. The dashed line is what naive linear interpolation on the screen would give.

```desmos-graph
left=-0.1; right=1.1; bottom=-0.6; top=0.3
---
y=\frac{-1+2x}{2+4x}|blue|0<x<1
y=-0.5+\frac{2}{3}x|red|dashed|0<x<1
(0.5,0)|label:3-D midpoint lands here
(0.5,-0.1667)|label:screen midpoint
```

*Screen position (blue) against position along the 3-D segment. The blue curve bends: equal steps in the world are unequal on the screen. Change the depth slope `4` in the denominator to `0` (the segment is now parallel to the screen) and the curve becomes the straight line; make it `20` and watch the far half crowd into a sliver.*

Second, the whole pipeline in under forty lines. Two wireframe cubes, one near and one far, are carried by a single composed matrix and a perspective divide.

```python
import numpy as np
import matplotlib.pyplot as plt

DIST = 6.0           # KNOB: camera distance from the near cube's front face
FOCAL = 0.5 * DIST   # KNOB: focal length; FOCAL proportional to DIST = "dolly zoom"
YAW = 30.0           # KNOB: rotate both cubes about the vertical axis (degrees)
ORTHO = False        # KNOB: True = orthographic (no division by depth)

V = np.array([[x, y, z, 1.0] for x in (-1, 1) for y in (-1, 1) for z in (-1, 1)])
E = [(a, b) for a in range(8) for b in range(a + 1, 8) if np.sum(V[a] != V[b]) == 1]

def T(tx, ty, tz):                       # homogeneous translation
    M = np.eye(4); M[:3, 3] = (tx, ty, tz); return M

def Ry(deg):                             # rotation about the vertical (y) axis
    c, s = np.cos(np.radians(deg)), np.sin(np.radians(deg))
    return np.array([[c, 0, s, 0], [0, 1, 0, 0], [-s, 0, c, 0], [0, 0, 0, 1]])

P = np.array([[FOCAL, 0, 0, 0],          # perspective: w' = z, then divide
              [0, FOCAL, 0, 0],
              [0, 0, 1, 0],
              [0, 0, 1, 0]])
if ORTHO:
    P = np.diag([FOCAL / DIST, FOCAL / DIST, 1.0, 1.0])

view = T(0, 0, DIST + 1)                 # camera at origin looking down +z
cubes = {"near cube": T(0, 0, 0) @ Ry(YAW), "far cube": T(2.5, 0, 8) @ Ry(YAW)}
fig, ax = plt.subplots(figsize=(6, 4))
for (name, model), c in zip(cubes.items(), ["tab:blue", "tab:red"]):
    clip = (P @ view @ model @ V.T).T    # one matrix product carries every vertex
    s = clip[:, :2] / clip[:, 3:4]       # the perspective divide
    for a, b in E:
        ax.plot(s[[a, b], 0], s[[a, b], 1], color=c)
    ax.plot([], [], color=c, label=name)
ax.set_aspect("equal"); ax.set_xlim(-1.5, 1.5); ax.set_ylim(-1, 1)
ax.legend(); ax.set_title(f"DIST={DIST}, FOCAL={FOCAL:.2f}, ortho={ORTHO}")
plt.show()
```

Things to try:

- **The dolly zoom.** Because `FOCAL` is tied to `DIST`, the near cube's front face always has half-height $0.5$ on screen. Set `DIST = 2.0`: the far cube's front face shrinks to half-height $0.1$. Set `DIST = 20.0`: it grows to about $0.36$. The subject stays put while the background breathes. Hitchcock used this in *Vertigo* (1958) to make a staircase lurch.
- **Swap the order.** Replace `T(2.5, 0, 8) @ Ry(YAW)` by `Ry(YAW) @ T(2.5, 0, 8)` and vary `YAW`. The far cube now swings around the near cube (the world origin) instead of turning in place: the solar-system bug, reproduced.
- **Go orthographic.** Set `ORTHO = True`. Parallel edges stay parallel and the two cubes are the same size: an engineering drawing, not a photograph.
- **Walk through the eye.** Set `DIST = -1.0` and `FOCAL = 1.0`, so the near cube straddles the camera. Vertices behind the eye get negative $w$ and the drawing explodes. That is the clipping problem of §69.4, in your hands.

## 69.6 The challenge

1. **Intuition.** Railway tracks are parallel in the world, yet in a photograph they meet at a *vanishing point*. Using Result 3 and the fact that all lines with direction $\mathbf{d}$ share the homogeneous "point" $(\mathbf{d},0)$, explain where the vanishing point is on the screen. Why do lines parallel to the image plane have no vanishing point? Where do all the vanishing points of horizontal lines lie?
2. **Derivation.** Show that $1/z_c$ is an *affine* function of the screen coordinates $(x_s,y_s)$ across the image of any planar triangle (start from the plane equation $\alpha x_c+\beta y_c+\gamma z_c=\delta$). Deduce that a texture coordinate $q$ is interpolated correctly on the screen by interpolating $q/z_c$ and $1/z_c$ linearly and dividing.
3. **Implementation.** Extend the experiment script into a tiny **rasteriser**: fill each face of the cube with a flat colour, resolve visibility with a per-pixel depth buffer on a $160\times120$ grid instead of drawing faces in order, and compute flat shading from face normals transformed by the inverse transpose. Add a non-uniform scale to the model matrix and confirm that the wrong (non-inverse-transpose) normals visibly break the shading.
4. **Open investigation.** Real lenses are not pinholes. Wide-angle lenses bend straight lines near the edges (barrel distortion), often modelled by $r_d=r(1+k_1r^2+k_2r^4)$ in screen radius $r$. Distort your rasteriser's output with this model and then try to *undo* it from a photo of a checkerboard, estimating $k_1,k_2$ by least squares. How many squares do you need? How does the error behave near the corners? Camera calibration is a field in its own right; see how far you get with only this chapter.

## 69.7 Loose threads

The programmer's Earth now spins in place, because she knows that matrices act right to left and that the rightmost one acts in the object's own frame. The GPU's impossible trick has a name: homogeneous coordinates turn translation into a 4-D shear and projection into "copy depth into $w$, then divide". Three theorems carry the load: affine maps preserve affine combinations, normals transform by the inverse transpose, and perspective keeps lines straight while distorting distances along them.

What we did not tell you:

- **Rotation is the weak link.** We stored orientation as nine numbers in a matrix, though a rotation has only three degrees of freedom. Nine numbers drift, cannot be averaged, and give no sensible way to turn smoothly from one orientation to another.

> [!mystery] Halfway between two orientations
> A camera faces north; a second later it must face straight up and slightly east. Average the two rotation matrices entrywise and you do not get a rotation at all: the result shrinks and shears the world. What *is* halfway between two orientations, and why does the honest answer need a four-dimensional number system invented on a Dublin bridge in 1843? That is [Ch. 70](70-quaternions.md).

- **Points at infinity.** The homogeneous points with $w=0$ behave like genuine points: parallel lines meet there, and a projective transformation can move them into view as vanishing points. Adding them turns 3-D space into *projective space*, where every two distinct lines in a plane meet in exactly one point.

> [!mystery] A world with one side
> Collapse each line through the origin of $\mathbb{R}^3$ to a single point and you get the *real projective plane*: the world of homogeneous coordinates $(x,y,w)$ up to scale. It is a closed surface that cannot be built in ordinary 3-D space without crossing itself, and it has only one side. Why, and what does that mean for a creature living on it? The answer is in [Ch. 58](../e-abstract-math/58-topology.md).

- **The inverse problem.** Projection throws away depth. Recovering it from two or more photographs (stereo, structure from motion) is an inverse problem with its own geometry, the *epipolar* constraint, and its own instabilities; see [Ch. 91](../j-frontier/91-inverse-problems.md).
- **The other way to draw.** The pipeline here pushes triangles toward the eye. The opposite approach pulls light backward from the eye into the scene, one ray per pixel, and gets shadows, mirrors and glass almost for free: [Ch. 71](71-ray-tracing.md). Clipping, culling and finding what a ray hits quickly are algorithms on shapes, in [Ch. 75](75-computational-geometry.md).

## Further reading

- J. F. Hughes, A. van Dam, M. McGuire, D. F. Sklar, J. D. Foley, S. K. Feiner, K. Akeley, *Computer Graphics: Principles and Practice*, 3rd ed., Addison-Wesley, 2013. Chapters on transformations and the camera, with the point/vector distinction treated carefully.
- S. Marschner, P. Shirley, *Fundamentals of Computer Graphics*, 5th ed., CRC Press, 2021. Compact derivations of the viewing pipeline and perspective-correct interpolation.
- E. Lengyel, *Mathematics for 3D Game Programming and Computer Graphics*, 3rd ed., Cengage, 2011. Homogeneous coordinates, normals and depth precision from a game engine's point of view.
- R. Hartley, A. Zisserman, *Multiple View Geometry in Computer Vision*, 2nd ed., Cambridge University Press, 2004. The pinhole camera, projective geometry and calibration, taken all the way.
