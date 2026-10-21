# Chapter 25 — Mathematical Modeling

[← Previous](24-computational-mathematics.md) · [Table of Contents](../00-front/table-of-contents.md) · [Next →](../00-front/table-of-contents.md)

> **Prerequisites:** Ch. 1 (Mathematical Language and Algebra), Ch. 4 (Differential Calculus), Ch. 6 (Multivariable Calculus), Ch. 8 (Matrices and Linear Systems), Ch. 11–12 (Ordinary and Linear Differential Equations), Ch. 14 (Laplace Transforms), Ch. 20 (Statistical Analysis), Ch. 22–23 (Numerical Methods, Optimization)
> **Learning outcomes:** After this chapter you can …
> - run the modeling cycle deliberately: state the purpose, choose the system boundary, write conservation laws and constitutive relations, and decide what to neglect;
> - separate state variables, parameters, inputs, disturbances and outputs, and put a lumped model into state-space form;
> - apply the Buckingham Pi theorem to reduce the number of variables in a problem and to obtain scaling laws without solving any equation;
> - nondimensionalize a model (RLC circuit, thermal network, heat equation, nonlinear oscillator) so that its behaviour is governed by a minimal set of dimensionless groups such as $\zeta$, $\omega_n$, Biot and Fourier numbers;
> - linearize a nonlinear model about an operating point, compute small-signal parameters, and estimate the amplitude range over which the linear model is trustworthy;
> - validate a model with residual analysis, goodness-of-fit measures and cross-validation, identify its parameters from step responses or by least squares, test identifiability, and choose between analytical and numerical solution.

Every chapter of this book has taken a mathematical problem as given: *solve this differential equation*, *invert this transform*, *minimize this function*. In engineering practice the problem is almost never given. Someone hands you a circuit board that oscillates when it should not, a sensor whose reading drifts with temperature, a transmission line whose losses are higher than the datasheet suggests, and the first and hardest step is to decide **which mathematics describes it**. That step is *mathematical modeling*.

Modeling is a craft as much as a science, but it is not mysterious. It has a recognisable cycle (Section 25.1), a vocabulary (25.2), a set of tools that extract information *before* any equation is solved (dimensional analysis and scaling, 25.3–25.4), a standard way to make hard problems tractable (linearization, 25.5), a discipline of stating and checking assumptions (25.6), and a way to confront the model with data — validating its structure (25.7) and estimating its numbers (25.8). The chapter ends with the choice every modeller faces between closed-form and computed solutions (25.9), and a short synthesis that ties the methods of the whole book together.

Because this is the final chapter, it leans on almost everything before it. Cross-references are given where a technique is used; the emphasis here is on *when* and *why* to use each one, not on re-deriving it.

## 25.1 Model formulation

### What a model is, and what it is for

> [!definition] Mathematical model
> A **mathematical model** of a system is a set of mathematical relations — algebraic equations, differential or difference equations, probability distributions, inequalities — together with a statement of the variables they relate, the assumptions under which they hold, and the purpose for which they are intended.

The last clause is the one most often forgotten. A model is a tool built for a job, and "Is the model correct?" is the wrong question; the right one is "Is the model *adequate for this purpose* over *this range of conditions*?" The same resistor is

- an ideal $v = Ri$ for a DC bias calculation,
- $v = R(T)\,i$ with $R(T) = R_0\,[1+\alpha(T-T_0)]$ for a precision current shunt,
- $R$ in series with a lead inductance $L_s$ and in parallel with a pad capacitance $C_p$ for a 2 GHz matching network,
- $R$ plus a Johnson noise source with power spectral density $4kTR$ for a low-noise amplifier front end (Ch. 21).

None of these is "the" model of a resistor. Each is the simplest model that answers a particular question. George Box's often-quoted remark that all models are wrong but some are useful expresses exactly this point.

> [!tip] Engineering Note — the purpose fixes the fidelity
> Before writing an equation, write one sentence: *"This model must predict ___ to within ___ for inputs in the range ___."* For example: "predict the junction temperature of Q3 to within 5 °C for dissipations from 0 to 25 W, on time scales longer than 1 s." That sentence decides the boundary, the time scale, the variables to keep and the accuracy needed from every parameter. Models built without it tend to be either too crude to answer the question or too elaborate to be understood, calibrated or trusted.

### The modeling cycle

Modeling is iterative. A widely used description has six stages:

1. **Problem statement and purpose.** What question must the model answer? What accuracy, what operating range, what time scale?
2. **Conceptualization.** Choose the system boundary, identify the inputs crossing it and the outputs of interest, and decide which physical effects matter. Sketch an idealized system (a circuit diagram, a lumped thermal network, a free-body diagram).
3. **Mathematical formulation.** Write the governing equations: conservation laws plus constitutive relations, with initial and boundary conditions.
4. **Analysis and solution.** Simplify (dimensional analysis, scaling, linearization), then solve — analytically or numerically.
5. **Validation.** Compare predictions with independent measurements. Are the residuals consistent with measurement noise? Does the model work on data not used to build it?
6. **Use and refinement.** Use the model for design or prediction; when it fails, return to stage 2 or 3 with what the failure taught.

$$
\text{purpose}\;\to\;\text{concept}\;\to\;\text{equations}\;\to\;\text{solution}\;\to\;\text{validation}\;\to\;\text{use}\;\circlearrowleft
$$

The loop back is essential. The first model is rarely the last; but a *documented* first model, with explicit assumptions, is the only thing that can be systematically improved.

### Building blocks: conservation laws and constitutive relations

Almost every physical model in electrical engineering is assembled from two kinds of equation.

> [!definition] Balance (conservation) law and constitutive relation
> A **balance law** states that the rate of change of an extensive quantity inside a boundary equals its net inflow plus its internal production:
> $$
> \frac{d}{dt}(\text{amount inside}) = (\text{rate in}) - (\text{rate out}) + (\text{rate generated}).
> $$
> A **constitutive relation** characterises a particular material or component: it links the *flow* of the conserved quantity to the *effort* (potential difference) that drives it, or links the stored amount to the effort.

Balance laws are universal and exact (within classical physics); constitutive relations are empirical, approximate and valid over limited ranges. That asymmetry tells you where model error usually lives: in the constitutive relations and in what was left out, rather than in the conservation laws themselves.

| Domain | Conserved quantity | Effort (potential) | Flow | Storage relation | Dissipation relation |
|---|---|---|---|---|---|
| Electrical | charge $q$ | voltage $v$ | current $i$ | $q = Cv$ | $v = Ri$ |
| Magnetic (flux linkage) | flux linkage $\lambda$ | — | — | $\lambda = Li$ | — |
| Thermal | heat $Q$ | temperature $T$ | heat flow $P$ | $Q = C_{\text{th}}T$ | $\Delta T = R_{\text{th}}P$ |
| Mechanical (translation) | momentum $p$ | velocity $v$ | force $F$ | $p = mv$ | $F = bv$ |
| Mechanical (rotation) | angular momentum | angular velocity $\omega$ | torque $\tau$ | $J\omega$ | $\tau = b\omega$ |
| Fluid | volume $V$ | pressure $p$ | volume flow $Q$ | $V = C_f\,p$ | $p = R_f\,Q$ |

Kirchhoff's current law is charge conservation at a node; Kirchhoff's voltage law is energy conservation around a loop (in the quasi-static, lumped regime of Section 25.6). The *analogy* across domains is not decorative: it lets an electrical engineer write the thermal model of a heat sink or the mechanical model of a motor shaft as an equivalent circuit, and reuse every circuit-analysis technique of Ch. 12 and Ch. 14.

### Levels of description

Models can be classified along several independent axes. Choosing a position on each axis is part of formulation.

| Axis | Choices | Mathematics that results |
|---|---|---|
| Spatial | lumped / distributed | ODEs (Ch. 11–12) / PDEs (Ch. 13) |
| Time | static / dynamic; continuous / discrete | algebraic / differential / difference equations (Ch. 17) |
| Linearity | linear / nonlinear | superposition, transforms (Ch. 14–17) / iteration, simulation (Ch. 22) |
| Time dependence of parameters | time-invariant / time-varying | constant-coefficient / variable-coefficient equations |
| Certainty | deterministic / stochastic | ODEs / random processes (Ch. 21) |
| Origin | white-box (physics) / grey-box / black-box (data) | first-principles equations / fitted structure / regression, ARX, neural nets |

A **white-box** model is derived entirely from physical laws with parameters known from geometry and material data. A **black-box** model is a flexible input–output relation fitted to data with no physical interpretation of its coefficients. Most useful engineering models are **grey-box**: the structure comes from physics (a first-order thermal lag, a second-order mechanical resonance), and the parameter values come from measurement (Section 25.8).

> [!example] Worked Example 25.1 — An RC charging model from first principles
> **Problem.** A sensor node draws its supply from a capacitor $C = 100\ \mathrm{\mu F}$ charged through a resistor $R = 1\ \mathrm{k\Omega}$ from a source $V_s = 5\ \mathrm{V}$ that is switched on at $t=0$. The load draws a constant current $I_L = 1\ \mathrm{mA}$. Formulate a model for the capacitor voltage $v(t)$ and find its steady state and time constant.
>
> **Step 1 — purpose.** Predict $v(t)$ to a few percent over the first second after switch-on, to check the node can start when $v > 3\ \mathrm{V}$.
>
> **Step 2 — concept.** Boundary: the capacitor node. Inputs: $V_s$ (switched), $I_L$. Neglect capacitor leakage and equivalent series resistance (justified later if the leakage current is $\ll I_L$).
>
> **Step 3 — equations.** Charge balance at the node (KCL):
> $$
> C\frac{dv}{dt} = \frac{V_s - v}{R} - I_L, \qquad v(0) = 0.
> $$
> Rearranged:
> $$
> RC\,\frac{dv}{dt} + v = V_s - RI_L .
> $$
>
> **Step 4 — solution.** This is the standard first-order linear equation of Ch. 11 with time constant $\tau = RC = 10^3\times 10^{-4} = 0.1\ \mathrm{s}$ and steady state $v_\infty = V_s - RI_L = 5 - 1 = 4\ \mathrm{V}$:
> $$
> v(t) = v_\infty\left(1 - e^{-t/\tau}\right) = 4\left(1-e^{-t/0.1}\right)\ \mathrm{V}.
> $$
>
> **Step 5 — use.** The node reaches $3\ \mathrm{V}$ when $1-e^{-t/\tau} = 0.75$, i.e. $t = \tau\ln 4 = 0.139\ \mathrm{s}$.
>
> **Dimensional check.** $[RC] = \mathrm{(V/A)(A\,s/V)} = \mathrm{s}$ and $[RI_L] = \mathrm{V}$, so every term has the right units.

The example is elementary on purpose: notice that the mathematics (Step 4) occupied one line, while purpose, boundary and assumptions occupied the rest. That proportion is typical of real modeling.

> [!example] Worked Example 25.2 — A lumped thermal model of a power transistor
> **Problem.** A MOSFET dissipating $P = 20\ \mathrm{W}$ is mounted on a heat sink in ambient air at $T_a = 40\ ^\circ\mathrm{C}$. Datasheets give the junction-to-case resistance $R_{jc} = 1.5\ \mathrm{K/W}$, the interface (thermal pad) resistance $R_{cs} = 0.5\ \mathrm{K/W}$ and the sink-to-ambient resistance $R_{sa} = 2.0\ \mathrm{K/W}$. The heat sink (aluminium, $m = 0.25\ \mathrm{kg}$, $c_p = 900\ \mathrm{J/(kg\,K)}$) has much larger heat capacity than the transistor package. Formulate a model for the junction temperature after power is applied.
>
> **Concept.** Heat flows junction → case → sink → ambient through three series resistances. Heat is *stored* mainly in the sink, with thermal capacitance $C_s = mc_p = 0.25\times900 = 225\ \mathrm{J/K}$. Neglect the package's own capacitance (it charges in seconds, while the sink takes minutes).
>
> **Equations.** Energy balance on the sink, with sink temperature $T_s$:
> $$
> C_s\frac{dT_s}{dt} = P - \frac{T_s - T_a}{R_{sa}},
> $$
> and, because the package stores negligible heat, the full $P$ flows through $R_{jc}+R_{cs}$ at every instant:
> $$
> T_j = T_s + (R_{jc}+R_{cs})P .
> $$
>
> **Solution.** With $T_s(0) = T_a$, the sink equation is first order with $\tau = R_{sa}C_s = 2.0\times225 = 450\ \mathrm{s}$ and
> $$
> T_s(t) = T_a + PR_{sa}\left(1-e^{-t/\tau}\right).
> $$
> Hence
> $$
> T_j(t) = T_a + (R_{jc}+R_{cs})P + R_{sa}P\left(1-e^{-t/\tau}\right) = 40 + 40 + 40\left(1-e^{-t/450}\right)\ ^\circ\mathrm{C}.
> $$
> At switch-on the junction jumps (on the model's time scale) to $80\ ^\circ\mathrm{C}$, and it approaches $120\ ^\circ\mathrm{C}$ with a 7.5-minute time constant.
>
> **What the model omits.** The instantaneous jump is an artefact of neglecting the package capacitance; on a millisecond time scale the junction rises smoothly. That is acceptable for the stated purpose (steady state and minutes-scale warm-up) but would be wrong for a 10 ms overload pulse — there the package capacitance is the *only* thing that matters. Same device, different purpose, different model.

> [!example] Worked Example 25.3 — An electromechanical model: the DC motor
> **Problem.** Formulate a model for a permanent-magnet DC motor driven by voltage $v_a(t)$ and loaded by torque $\tau_L(t)$, suitable for designing a speed controller.
>
> **Concept.** Two coupled domains. Electrical: armature resistance $R_a$, inductance $L_a$, and back-EMF proportional to speed. Mechanical: rotor and load inertia $J$, viscous friction $b$. Coupling: the torque constant $K_t$ (torque per ampere) and back-EMF constant $K_e$ (volts per rad/s). Energy conservation across the air gap (electrical power $e\,i$ converted equals mechanical power $\tau\omega$) forces $K_t = K_e \equiv K$ in SI units.
>
> **Equations.** KVL round the armature and Newton's law for rotation:
> $$
> \begin{aligned}
> L_a\frac{di}{dt} &= v_a - R_a i - K\omega,\\
> J\frac{d\omega}{dt} &= K i - b\omega - \tau_L .
> \end{aligned}
> $$
> If shaft position $\theta$ is also wanted, add $\dot\theta = \omega$.
>
> **Check.** Units: $[K\omega] = (\mathrm{V\,s/rad})(\mathrm{rad/s}) = \mathrm{V}$; $[Ki] = (\mathrm{N\,m/A})(\mathrm{A}) = \mathrm{N\,m}$; and $\mathrm{V\,s} = \mathrm{N\,m/A}$ since $1\ \mathrm{V\,A} = 1\ \mathrm{W} = 1\ \mathrm{N\,m/s}$. Steady state with $\tau_L=0$: $\omega_\infty = Kv_a/(R_ab + K^2)$, which is the familiar no-load speed when $b$ is small, $\omega_\infty \approx v_a/K$.
>
> This model is linear only because friction was taken viscous. Real gearboxes add Coulomb friction $\tau_c\operatorname{sgn}(\omega)$, which is nonlinear and discontinuous at $\omega = 0$ — the reason speed controllers behave poorly near standstill. Whether to include it depends, once again, on the purpose.

> [!warning] Common Pitfall — modeling the device instead of the question
> A frequent failure is to start from "let us model everything about the motor": winding temperature, magnetic saturation, brush arcing, cogging torque, shaft flexibility. The result has dozens of parameters, most unknown, and no-one can tell which effects drive the answer. Start with the simplest model that could possibly answer the question; add an effect only when validation (Section 25.7) shows the simple model fails *in a way that effect explains*.

## 25.2 State variables and parameters

### Five kinds of quantity

Once a model is written down, every symbol in it plays one of a small number of roles. Getting the roles straight is what turns a page of equations into something that can be simulated, identified and controlled.

> [!definition] Roles of quantities in a model
> - **State variables** $\mathbf{x}(t)\in\mathbb{R}^n$: the minimal set of quantities whose present values, together with future inputs, determine all future behaviour. They carry the system's *memory*.
> - **Inputs** $\mathbf{u}(t)$: externally imposed signals that the modeller (or a controller) chooses — a drive voltage, a commanded current, a heater power.
> - **Disturbances** $\mathbf{w}(t)$: external signals that act on the system but are *not* chosen — ambient temperature, load torque, supply ripple, noise.
> - **Parameters** $\mathbf{p}$: quantities that characterise the system and are constant (or vary much more slowly than the states) over the time scale of interest — $R$, $L$, $C$, $J$, $K$, $R_{\text{th}}$.
> - **Outputs** $\mathbf{y}(t)$: quantities that are measured or are of interest; functions of the states and inputs.
>
> The general lumped, continuous-time model is
> $$
> \dot{\mathbf{x}} = \mathbf{f}(\mathbf{x},\mathbf{u},\mathbf{w};\mathbf{p},t), \qquad \mathbf{y} = \mathbf{h}(\mathbf{x},\mathbf{u};\mathbf{p},t), \qquad \mathbf{x}(t_0) = \mathbf{x}_0 .
> $$

Constants that are *universal* — $\pi$, the Boltzmann constant $k$, the elementary charge $q$ — are sometimes called **constants** to distinguish them from parameters, which belong to a particular device.

The distinction between a parameter and a state is a matter of *time scale*, not of physics. The temperature of a resistor is a parameter in a circuit simulation lasting a microsecond, because it cannot change appreciably in that time; it is a state in a thermal simulation lasting an hour. The same goes for battery state of charge, capacitor ageing, or a slowly drifting oscillator frequency. A useful rule:

> [!tip] Engineering Note — parameter or state?
> If a quantity changes by less than the required model accuracy over the time window of interest, treat it as a parameter. If it changes by more, it must be a state with its own equation, or the model must be re-identified as it drifts (adaptive or gain-scheduled models).

The distinction between an input and a disturbance is a matter of *control*. Ambient temperature is a disturbance for an amplifier, but an input for the environmental chamber used to test it.

### Why states: the state-space form

The state-space form $\dot{\mathbf{x}} = \mathbf{f}(\mathbf{x},\mathbf{u})$ is the standard because

- it is first order, so every numerical ODE solver of Ch. 22 accepts it directly;
- its linear version $\dot{\mathbf{x}} = \mathbf{A}\mathbf{x} + \mathbf{B}\mathbf{u}$ is analysed by the eigenvalue methods of Ch. 10 (modes, stability) and the matrix exponential;
- the initial condition $\mathbf{x}(t_0)$ is exactly the information needed — no more, no less.

For linear time-invariant (LTI) systems,

$$
\dot{\mathbf{x}} = \mathbf{A}\mathbf{x} + \mathbf{B}\mathbf{u}, \qquad \mathbf{y} = \mathbf{C}\mathbf{x} + \mathbf{D}\mathbf{u},
$$

with solution

$$
\mathbf{x}(t) = e^{\mathbf{A}t}\mathbf{x}(0) + \int_0^t e^{\mathbf{A}(t-\tau)}\mathbf{B}\mathbf{u}(\tau)\,d\tau
$$

and transfer function (Ch. 14), for zero initial state,

$$
\mathbf{G}(s) = \mathbf{C}(s\mathbf{I}-\mathbf{A})^{-1}\mathbf{B} + \mathbf{D}.
$$

The poles of $\mathbf{G}(s)$ are among the eigenvalues of $\mathbf{A}$; they coincide with all of them when the realization is minimal (controllable and observable).

### Choosing state variables

> [!theorem] Natural state variables for lumped physical systems
> In a lumped network of ideal elements, a natural choice of states is one *energy variable per independent energy-storage element*: the capacitor voltages (or charges) and inductor currents (or flux linkages) in a circuit; the masses' velocities and springs' deflections in a mechanism; the node temperatures of thermal capacitances. The stored energy is then a quadratic form in the states, e.g. $E = \tfrac12\sum C_k v_k^2 + \tfrac12\sum L_m i_m^2$.

*Why it works (sketch).* The constitutive relations of storage elements are the only ones involving derivatives ($i = C\,dv/dt$, $v = L\,di/dt$). Resistive elements and Kirchhoff's laws are algebraic. Solving the algebraic equations for the capacitor currents and inductor voltages in terms of the capacitor voltages, inductor currents and sources gives exactly $\dot{\mathbf{x}} = \mathbf{f}(\mathbf{x},\mathbf{u})$.

The word *independent* matters. Two capacitors directly in parallel share one voltage and contribute one state; a loop consisting only of capacitors and voltage sources, or a cut-set consisting only of inductors and current sources, removes a state. Such configurations are called **degenerate**; the order of the network is

$$
n = n_C + n_L - n_{\text{C-loops}} - n_{\text{L-cutsets}},
$$

where $n_{\text{C-loops}}$ counts independent loops of capacitors (and voltage sources only) and $n_{\text{L-cutsets}}$ counts independent cut-sets of inductors (and current sources only).

States are not unique. Any invertible linear transformation $\mathbf{z} = \mathbf{T}\mathbf{x}$ gives an equally valid state vector, with

$$
\dot{\mathbf{z}} = \mathbf{T}\mathbf{A}\mathbf{T}^{-1}\mathbf{z} + \mathbf{T}\mathbf{B}\mathbf{u}, \qquad \mathbf{y} = \mathbf{C}\mathbf{T}^{-1}\mathbf{z} + \mathbf{D}\mathbf{u}.
$$

The eigenvalues ($\mathbf{T}\mathbf{A}\mathbf{T}^{-1}$ is similar to $\mathbf{A}$, Ch. 10) and the transfer function are unchanged: $\mathbf{C}\mathbf{T}^{-1}(s\mathbf{I}-\mathbf{T}\mathbf{A}\mathbf{T}^{-1})^{-1}\mathbf{T}\mathbf{B} = \mathbf{C}(s\mathbf{I}-\mathbf{A})^{-1}\mathbf{B}$. What changes is interpretability and numerical conditioning — physical states with similar magnitudes are usually the best-conditioned choice.

> [!example] Worked Example 25.4 — State-space model of the series RLC circuit
> **Problem.** A series RLC circuit with $R = 10\ \Omega$, $L = 1\ \mathrm{mH}$, $C = 1\ \mathrm{\mu F}$ is driven by $v_s(t)$. The output is the capacitor voltage. Write the state-space model, find the eigenvalues and interpret them.
>
> **States.** One capacitor and one inductor, no degeneracy: $n = 2$, $\mathbf{x} = [v_C,\ i]^\mathsf{T}$.
>
> **Equations.** Capacitor: $C\dot v_C = i$. KVL: $L\dot i = v_s - Ri - v_C$. Hence
> $$
> \frac{d}{dt}\begin{bmatrix} v_C\\ i\end{bmatrix} = \begin{bmatrix} 0 & 1/C\\ -1/L & -R/L\end{bmatrix}\begin{bmatrix} v_C\\ i\end{bmatrix} + \begin{bmatrix} 0\\ 1/L\end{bmatrix}v_s, \qquad y = \begin{bmatrix}1 & 0\end{bmatrix}\mathbf{x}.
> $$
>
> **Eigenvalues.** $\det(\lambda\mathbf{I}-\mathbf{A}) = \lambda^2 + (R/L)\lambda + 1/(LC) = 0$. Numerically $R/L = 10^4\ \mathrm{s^{-1}}$ and $1/(LC) = 10^9\ \mathrm{s^{-2}}$, so
> $$
> \lambda = -5000 \pm j\sqrt{10^9 - 2.5\times10^7} = -5000 \pm j\,31\,225\ \mathrm{s^{-1}}.
> $$
> The real part gives an envelope decay time $1/5000 = 0.2\ \mathrm{ms}$; the imaginary part gives a ringing frequency $31\,225/(2\pi) = 4.97\ \mathrm{kHz}$.
>
> **Transfer function.** $G(s) = \mathbf{C}(s\mathbf{I}-\mathbf{A})^{-1}\mathbf{B} = \dfrac{1/(LC)}{s^2 + (R/L)s + 1/(LC)}$, the standard second-order low-pass form whose dimensionless parameters are extracted in Section 25.4.

> [!example] Worked Example 25.5 — Counting states in a degenerate circuit
> **Problem.** A voltage source $v_s$ feeds a resistor $R_1$ into node A. Node A has capacitor $C_1$ to ground and capacitor $C_2$ to node B; node B has $C_3$ to ground and $R_2$ to ground. How many states? Now suppose $R_1$ is removed so the source connects directly to node A. How many states?
>
> **Case 1.** Three capacitors. Is there a loop of capacitors only? The loop ground–$C_1$–A–$C_2$–B–$C_3$–ground contains only capacitors, so $n_{\text{C-loops}} = 1$ and $n = 3 - 1 = 2$. Indeed $v_{C_1} = v_{C_2} + v_{C_3}$ (KVL round the loop), so only two capacitor voltages are independent; choose $\mathbf{x} = [v_{C_1}, v_{C_3}]^\mathsf{T}$.
>
> **Case 2.** With the source wired straight to node A, $v_{C_1} = v_s$ is now forced: the loop source–$C_1$ is a second loop of capacitors and voltage sources only. So $n = 3 - 2 = 1$, with state $v_{C_3}$ (and $v_{C_2} = v_s - v_{C_3}$).
>
> **Consequence.** In case 2 a step in $v_s$ produces an *impulse* of current into $C_1$ and $C_2$ — the ideal model predicts infinite current. A circuit simulator will either fail or insert a tiny series resistance. The degeneracy is a symptom of an idealization (zero source resistance) that is invalid on short time scales, and is a reminder that assumptions must be checked (Section 25.6).

> [!example] Worked Example 25.6 — Changing state coordinates for the DC motor
> **Problem.** For the DC motor of Example 25.3 with $\tau_L = 0$, the states $\mathbf{x} = [i,\ \omega]^\mathsf{T}$ give
> $$
> \mathbf{A} = \begin{bmatrix} -R_a/L_a & -K/L_a\\ K/J & -b/J\end{bmatrix}, \quad \mathbf{B} = \begin{bmatrix}1/L_a\\ 0\end{bmatrix}.
> $$
> Use instead the momentum-like states $\mathbf{z} = [\lambda,\ h]^\mathsf{T} = [L_a i,\ J\omega]^\mathsf{T}$ (flux linkage and angular momentum). Find the new matrices and confirm the eigenvalues are unchanged. Use $R_a = 2\ \Omega$, $L_a = 5\ \mathrm{mH}$, $K = 0.05\ \mathrm{V\,s/rad}$, $J = 2\times10^{-5}\ \mathrm{kg\,m^2}$, $b = 10^{-5}\ \mathrm{N\,m\,s/rad}$.
>
> **Transformation.** $\mathbf{T} = \operatorname{diag}(L_a, J)$, so
> $$
> \mathbf{T}\mathbf{A}\mathbf{T}^{-1} = \begin{bmatrix} -R_a/L_a & -K/J\\ K/L_a & -b/J\end{bmatrix}, \qquad \mathbf{T}\mathbf{B} = \begin{bmatrix}1\\0\end{bmatrix}.
> $$
> In the new coordinates the input enters with unit gain ($\dot\lambda = v_a - \dots$), which is the energy-based (port-Hamiltonian) form preferred in some control designs.
>
> **Eigenvalues.** Trace $= -R_a/L_a - b/J = -400 - 0.5 = -400.5\ \mathrm{s^{-1}}$; determinant $= (R_a b + K^2)/(L_a J) = (2\times10^{-5} + 2.5\times10^{-3})/(10^{-7}) = 25\,200\ \mathrm{s^{-2}}$. Both matrices share trace and determinant, so both give
> $$
> \lambda^2 + 400.5\lambda + 25\,200 = 0 \;\Rightarrow\; \lambda = \frac{-400.5 \pm \sqrt{400.5^2 - 4(25\,200)}}{2} = \frac{-400.5 \pm 244.1}{2},
> $$
> that is $\lambda_1 \approx -78.2\ \mathrm{s^{-1}}$ and $\lambda_2 \approx -322.3\ \mathrm{s^{-1}}$.
> The slow mode (time constant $12.8\ \mathrm{ms}$) is essentially mechanical, related to $JR_a/K^2 = 16\ \mathrm{ms}$; the fast mode ($3.1\ \mathrm{ms}$) is essentially electrical, related to $L_a/R_a = 2.5\ \mathrm{ms}$. Because the domains are coupled through $K$, neither eigenvalue equals the naive single-domain estimate — a small but genuine lesson about decoupled "back-of-envelope" models.

> [!warning] Common Pitfall — treating a derivative of an input as a state
> In a model like $\ddot y + a\dot y + by = c\dot u + du$, it is tempting to take $\mathbf{x} = [y, \dot y]^\mathsf{T}$. Then $\ddot y$ contains $\dot u$, and the "state equation" needs the derivative of the input — which does not exist for a step input. The fix is a proper realization (e.g. controllable canonical form, Ch. 14) in which the states are internal variables and $u$ enters without differentiation; the output is then $y = \mathbf{C}\mathbf{x}$ with $\mathbf{C} = [d,\ c]$.

### Algebraic variables and DAEs

Many models contain variables that obey *algebraic* rather than differential equations: node voltages with no capacitor, currents through resistors, a junction temperature determined by an instantaneous heat flow (Example 25.2). The full model is then a **differential-algebraic equation** (DAE)

$$
\dot{\mathbf{x}} = \mathbf{f}(\mathbf{x},\mathbf{z},\mathbf{u}), \qquad \mathbf{0} = \mathbf{g}(\mathbf{x},\mathbf{z},\mathbf{u}),
$$

with differential states $\mathbf{x}$ and algebraic variables $\mathbf{z}$. If $\partial\mathbf{g}/\partial\mathbf{z}$ is invertible (an *index-1* DAE), the implicit function theorem (Ch. 6) lets us eliminate $\mathbf{z}$ locally and recover an ODE. Circuit simulators of the SPICE family solve the DAE form directly using modified nodal analysis, because eliminating every algebraic variable by hand would destroy the sparsity of the network equations.

### Discrete-time states

Sampled-data systems and digital filters have states that evolve by a difference equation,

$$
\mathbf{x}[n+1] = \mathbf{A}_d\mathbf{x}[n] + \mathbf{B}_d\mathbf{u}[n], \qquad \mathbf{y}[n] = \mathbf{C}\mathbf{x}[n] + \mathbf{D}\mathbf{u}[n].
$$

For a continuous plant whose input is held constant between samples (zero-order hold) at period $T$, the exact discretization is $\mathbf{A}_d = e^{\mathbf{A}T}$ and $\mathbf{B}_d = \int_0^T e^{\mathbf{A}\sigma}\,d\sigma\,\mathbf{B}$; the eigenvalues map as $z = e^{sT}$ (Ch. 17). The states of a digital filter are its delay-line contents; the parameters are its coefficients.

## 25.3 Dimensional analysis

### Motivation: information for free

Before solving a single equation — sometimes before even knowing the equations — one can deduce the *form* of a physical relation from the requirement that it be dimensionally consistent. Dimensional analysis

- detects errors (an equation that adds volts to amperes is wrong, full stop);
- reduces the number of variables an experiment or simulation must sweep, often dramatically;
- yields **scaling laws** ("double the length, quadruple the radiated power") that guide design before detailed analysis;
- organises results into dimensionless groups that are universal across sizes and materials.

Chapter 1 introduced dimensional consistency as an error check. This section turns it into a method.

### Dimensions, units and homogeneity

> [!definition] Dimension, base dimensions, dimensional formula
> A **dimension** is a physical kind of quantity (length, mass, time…); a **unit** is a particular standard amount of it (metre, kilogram, second). The SI uses seven **base dimensions**; for electrical engineering four or five suffice: mass $\mathrm{M}$, length $\mathrm{L}$, time $\mathrm{T}$, current $\mathrm{I}$, and temperature $\Theta$. Every derived quantity $Q$ has a **dimensional formula**
> $$
> [Q] = \mathrm{M}^{a_1}\mathrm{L}^{a_2}\mathrm{T}^{a_3}\mathrm{I}^{a_4}\Theta^{a_5},
> $$
> and the exponent vector $(a_1,\dots,a_5)$ is its **dimension vector**. A quantity is **dimensionless** if all exponents are zero.

Some dimension vectors used constantly in EE:

| Quantity | SI unit | $\mathrm{M}$ | $\mathrm{L}$ | $\mathrm{T}$ | $\mathrm{I}$ | $\Theta$ |
|---|---|---|---|---|---|---|
| Voltage $V$ | V | 1 | 2 | −3 | −1 | 0 |
| Resistance $R$ | Ω | 1 | 2 | −3 | −2 | 0 |
| Capacitance $C$ | F | −1 | −2 | 4 | 2 | 0 |
| Inductance $L$ | H | 1 | 2 | −2 | −2 | 0 |
| Permittivity $\varepsilon$ | F/m | −1 | −3 | 4 | 2 | 0 |
| Permeability $\mu$ | H/m | 1 | 1 | −2 | −2 | 0 |
| Conductivity $\sigma$ | S/m | −1 | −3 | 3 | 2 | 0 |
| Power $P$ | W | 1 | 2 | −3 | 0 | 0 |
| Thermal resistance $R_{\text{th}}$ | K/W | −1 | −2 | 3 | 0 | 1 |
| Thermal conductivity $k$ | W/(m K) | 1 | 1 | −3 | 0 | −1 |
| Heat-transfer coefficient $h$ | W/(m² K) | 1 | 0 | −3 | 0 | −1 |

> [!theorem] Principle of dimensional homogeneity
> A physically meaningful equation must have the same dimensions in every additive term. Consequently, the arguments of transcendental functions ($\exp$, $\ln$, $\sin$, …) must be dimensionless, and any relation among physical quantities can be rearranged so that it relates dimensionless combinations only.

The first sentence is the error check of Ch. 1. The second is the seed of the Pi theorem: a law of nature cannot depend on whether we measure in metres or feet, so it must be expressible in quantities that do not change when units change — dimensionless groups.

### The Buckingham Pi theorem

> [!theorem] Buckingham Pi theorem
> Suppose a physical relation involves $n$ quantities $Q_1,\dots,Q_n$,
> $$
> f(Q_1, Q_2, \dots, Q_n) = 0,
> $$
> and is dimensionally homogeneous. Let $\mathbf{D}$ be the $k\times n$ **dimension matrix** whose column $i$ is the dimension vector of $Q_i$, and let $r = \operatorname{rank}\mathbf{D}$. Then the relation is equivalent to one among $n-r$ independent dimensionless products
> $$
> F(\Pi_1, \Pi_2, \dots, \Pi_{n-r}) = 0, \qquad \Pi_m = Q_1^{a_{1m}}Q_2^{a_{2m}}\cdots Q_n^{a_{nm}},
> $$
> where the exponent vectors $\mathbf{a}_m = (a_{1m},\dots,a_{nm})^\mathsf{T}$ form a basis of the null space of $\mathbf{D}$.

*Proof sketch.* A product $\Pi = \prod Q_i^{a_i}$ has dimension vector $\mathbf{D}\mathbf{a}$, so it is dimensionless exactly when $\mathbf{D}\mathbf{a} = \mathbf{0}$. By the rank–nullity theorem (Ch. 9) the null space has dimension $n - r$, so there are exactly $n-r$ independent dimensionless products; any other dimensionless product is a product of powers of these. Now pick $r$ of the quantities whose dimension vectors are linearly independent (the **repeating variables**). Changing the size of each base unit rescales every $Q_i$ by a factor, and because the repeating variables span the dimension space, units can be chosen so that each repeating variable has numerical value 1. In those units the relation $f = 0$ depends only on the remaining $n - r$ quantities, whose numerical values are now precisely the $\Pi$ groups. Since the physical relation cannot depend on the choice of units, it must be a relation among the $\Pi$'s alone. $\square$

The theorem's power lies in the count: a relation among $n$ quantities becomes one among $n - r$ groups. A function of one variable needs one curve to describe it; a function of zero variables is just a constant.

> [!tip] Procedure — applying the Pi theorem
> 1. **List** every quantity that could influence the result, including the result itself. Omitting a relevant quantity is the main way the method fails.
> 2. **Write** the dimension matrix $\mathbf{D}$ and find its rank $r$.
> 3. **Choose** $r$ repeating variables with independent dimensions, preferably those that are "given" (geometry, material constants), not the unknown.
> 4. **Form** each remaining quantity into a $\Pi$ by multiplying by powers of the repeating variables and solving the linear equations for the exponents that make it dimensionless. (Equivalently: compute a null-space basis of $\mathbf{D}$.)
> 5. **Write** the result as $\Pi_1 = \Phi(\Pi_2,\dots,\Pi_{n-r})$ with the unknown in $\Pi_1$, and determine $\Phi$ by theory, simulation or experiment.

> [!example] Worked Example 25.7 — The period of a pendulum
> **Problem.** A simple pendulum has bob mass $m$, length $\ell$, release angle $\theta_0$, in gravity $g$. Find the form of the period $T_p$.
>
> **List.** $T_p, m, \ell, g, \theta_0$: $n = 5$.
>
> **Dimension matrix** (rows $\mathrm{M}, \mathrm{L}, \mathrm{T}$):
> $$
> \mathbf{D} = \begin{array}{c|ccccc} & T_p & m & \ell & g & \theta_0\\ \hline \mathrm{M} & 0 & 1 & 0 & 0 & 0\\ \mathrm{L} & 0 & 0 & 1 & 1 & 0\\ \mathrm{T} & 1 & 0 & 0 & -2 & 0\end{array}
> $$
> The columns of $m$, $\ell$, $g$ are independent, so $r = 3$ and there are $5-3 = 2$ groups.
>
> **Groups.** $\theta_0$ (an angle, radians, is already dimensionless): $\Pi_2 = \theta_0$. For $T_p$, try $\Pi_1 = T_p\,m^a\ell^b g^c$: M: $a = 0$; L: $b + c = 0$; T: $1 - 2c = 0$. So $c = \tfrac12$, $b = -\tfrac12$:
> $$
> \Pi_1 = T_p\sqrt{\frac{g}{\ell}}, \qquad T_p = \sqrt{\frac{\ell}{g}}\;\Phi(\theta_0).
> $$
>
> **Conclusions without solving anything.** (i) The period is independent of the mass — no dimensionless group can contain $m$, because $m$ is the only quantity carrying $\mathrm{M}$. (ii) It scales as $\sqrt{\ell}$: a pendulum four times longer swings half as fast. (iii) Only one curve $\Phi(\theta_0)$ need be measured. Solving the ODE (Section 25.5) gives $\Phi(0) = 2\pi$ and $\Phi(\theta_0) = 4K\!\left(\sin(\theta_0/2)\right)$ in terms of the complete elliptic integral $K$; dimensional analysis cannot find the constant $2\pi$, but it found everything else.

> [!example] Worked Example 25.8 — Dimensionless groups of the RLC circuit
> **Problem.** The free response of a series RLC circuit is characterised by a frequency $\omega$ (any characteristic frequency). Which dimensionless groups govern it?
>
> **List.** $\omega, L, C, R$: $n = 4$. **Dimension matrix** (rows $\mathrm{M,L,T,I}$):
> $$
> \mathbf{D} = \begin{array}{c|cccc} & \omega & L & C & R\\ \hline \mathrm{M} & 0 & 1 & -1 & 1\\ \mathrm{L} & 0 & 2 & -2 & 2\\ \mathrm{T} & -1 & -2 & 4 & -3\\ \mathrm{I} & 0 & -2 & 2 & -2\end{array}
> $$
> Rows M, L and I are proportional (each component appears only through the combination $\mathrm{M L^2 I^{-2}}$), so $r = 2$ and there are $4 - 2 = 2$ groups even though four base dimensions appear. This is why the rank, not the number of base dimensions, must be used.
>
> **Null space.** Solving $\mathbf{D}\mathbf{a} = \mathbf{0}$ gives basis vectors $(2,1,1,0)$ and $(0,-1,1,2)$ for $(\omega,L,C,R)$ — the second equivalent, after taking a square root, to $(0,-\tfrac12,\tfrac12,1)$. Hence
> $$
> \Pi_1 = \omega^2 LC \;\to\; \omega\sqrt{LC}, \qquad \Pi_2 = R\sqrt{\frac{C}{L}} .
> $$
> Any characteristic frequency is therefore $\omega = \Phi\!\left(R\sqrt{C/L}\right)/\sqrt{LC}$. Section 25.4 shows these are exactly $\omega_n = 1/\sqrt{LC}$ and $\zeta = \tfrac12 R\sqrt{C/L}$; the damped frequency $\omega_d = \omega_n\sqrt{1-\zeta^2}$ is one instance of the general form.

> [!example] Worked Example 25.9 — Skin depth from dimensions alone
> **Problem.** At angular frequency $\omega$, current in a conductor of conductivity $\sigma$ and permeability $\mu$ concentrates in a surface layer of thickness $\delta$. Assuming $\delta$ depends only on $\omega$, $\mu$, $\sigma$, find its form and estimate the scaling with frequency.
>
> **List and matrix.** $n = 4$: $\delta, \omega, \mu, \sigma$. Rows $\mathrm{M, L, T, I}$:
> $$
> \mathbf{D} = \begin{array}{c|cccc} & \delta & \omega & \mu & \sigma\\ \hline \mathrm{M} & 0 & 0 & 1 & -1\\ \mathrm{L} & 1 & 0 & 1 & -3\\ \mathrm{T} & 0 & -1 & -2 & 3\\ \mathrm{I} & 0 & 0 & -2 & 2\end{array}
> $$
> Rank $r = 3$, so a single group exists. Try $\delta^2\omega^a\mu^b\sigma^c$. M: $b - c = 0$; I: $-2b + 2c = 0$; L: $2 + b - 3c = 0 \Rightarrow c = 1$, $b = 1$; T: $-a - 2b + 3c = 0 \Rightarrow a = 1$. So
> $$
> \Pi = \delta^2\omega\mu\sigma = \text{const} \quad\Longrightarrow\quad \delta = \frac{c_0}{\sqrt{\omega\mu\sigma}} .
> $$
> With one group and no other, $\Phi$ is a pure constant. Solving the diffusion equation for the field in the conductor (Ch. 13) gives $c_0 = \sqrt2$, i.e. $\delta = \sqrt{2/(\omega\mu\sigma)}$.
>
> **Numbers.** For copper, $\sigma = 5.8\times10^{7}\ \mathrm{S/m}$, $\mu = \mu_0$: $\delta = 8.53\ \mathrm{mm}$ at 60 Hz, $66\ \mathrm{\mu m}$ at 1 MHz, $2.09\ \mathrm{\mu m}$ at 1 GHz — every factor of 100 in frequency shrinks $\delta$ by 10, exactly as $\delta\propto\omega^{-1/2}$ demands.

> [!example] Worked Example 25.10 — Radiated power of a short dipole (two groups)
> **Problem.** A short wire antenna of length $\ell$ carries a uniform sinusoidal current of amplitude $I_0$ at angular frequency $\omega$ in free space ($\mu_0$, $\varepsilon_0$). Find the form of the time-averaged radiated power $P$.
>
> **List.** $P, I_0, \ell, \omega, \mu_0, \varepsilon_0$: $n = 6$. **Matrix** (rows $\mathrm{M, L, T, I}$):
> $$
> \mathbf{D} = \begin{array}{c|cccccc} & P & I_0 & \ell & \omega & \mu_0 & \varepsilon_0\\ \hline \mathrm{M} & 1 & 0 & 0 & 0 & 1 & -1\\ \mathrm{L} & 2 & 0 & 1 & 0 & 1 & -3\\ \mathrm{T} & -3 & 0 & 0 & -1 & -2 & 4\\ \mathrm{I} & 0 & 1 & 0 & 0 & -2 & 2\end{array}
> $$
> Rank $r = 4$, so $6 - 4 = 2$ groups.
>
> **Null space.** A basis is $\mathbf{a}_1 = (-1, 2, 1, 1, 1, 0)$ and $\mathbf{a}_2 = (1, -2, 1, 1, 0, 1)$, i.e.
> $$
> \Pi_a = \frac{I_0^2\,\ell\,\omega\,\mu_0}{P}, \qquad \Pi_b = \frac{P\,\ell\,\omega\,\varepsilon_0}{I_0^2}.
> $$
> These are correct but not illuminating. Any products of powers are equally valid groups, so recombine: $\Pi_a\Pi_b = \omega^2\ell^2\mu_0\varepsilon_0 = (k\ell)^2$ with $k = \omega\sqrt{\mu_0\varepsilon_0} = 2\pi/\lambda$, and $\sqrt{\Pi_b/\Pi_a} = (P/I_0^2)\sqrt{\varepsilon_0/\mu_0} = P/(\eta_0 I_0^2)$, with $\eta_0 = \sqrt{\mu_0/\varepsilon_0}\approx 377\ \Omega$. Therefore
> $$
> P = \eta_0 I_0^2\;\Phi(k\ell).
> $$
>
> **Physics fills in $\Phi$.** For a short dipole ($k\ell\ll1$) the field solution gives $\Phi = (k\ell)^2/(12\pi)$, i.e. $P = \eta_0 I_0^2(k\ell)^2/(12\pi) = 40\pi^2(\ell/\lambda)^2I_0^2$ watts. For $\ell = 1\ \mathrm{m}$, $f = 10\ \mathrm{MHz}$ ($\lambda = 30\ \mathrm{m}$), $I_0 = 1\ \mathrm{A}$: $P = 40\pi^2/900 = 0.439\ \mathrm{W}$. Dimensional analysis alone already told us the answer is $\eta_0I_0^2$ times a function of electrical length, and that a geometrically similar antenna scaled together with the wavelength radiates the same power for the same current.
>
> **Choice of groups.** The lesson of the recombination step: the null space is unique, its basis is not. Choose groups that have recognised physical meanings (electrical length $k\ell$, power normalised to $\eta_0I_0^2$).

> [!example] Worked Example 25.11 — Forced-air cooling of a heat sink: Nusselt, Reynolds and Prandtl
> **Problem.** The convective heat-transfer coefficient $h$ of a heat-sink fin in an air stream depends on air speed $V$, fin length $\ell$, and the air's thermal conductivity $k$, density $\rho$, viscosity $\mu$ and specific heat $c_p$. Reduce the problem.
>
> **List.** $n = 7$; base dimensions $\mathrm{M, L, T}, \Theta$, and the dimension matrix has rank $r = 4$. So there are $7 - 4 = 3$ groups. A null-space computation gives the basis (after recombination):
> $$
> \mathrm{Nu} = \frac{h\ell}{k}, \qquad \mathrm{Re} = \frac{\rho V\ell}{\mu}, \qquad \mathrm{Pr} = \frac{\mu c_p}{k}.
> $$
> Check one: $[h\ell/k] = (\mathrm{M T^{-3}\Theta^{-1}})(\mathrm{L})/(\mathrm{M L T^{-3}\Theta^{-1}}) = 1$.
>
> **Result.** $\mathrm{Nu} = \Phi(\mathrm{Re}, \mathrm{Pr})$. A seven-variable problem has become a two-variable one; for air $\mathrm{Pr}\approx0.71$ is essentially fixed, so a *single curve* $\mathrm{Nu}(\mathrm{Re})$ — measured once in a wind tunnel on any convenient size of fin — serves every fin and every fan speed. Laminar flat-plate theory gives $\mathrm{Nu} = 0.664\,\mathrm{Re}^{1/2}\mathrm{Pr}^{1/3}$, from which $h\propto\sqrt{V/\ell}$: doubling the fan speed raises $h$ by only $\sqrt2$.

> [!warning] Common Pitfall — the missing variable
> Dimensional analysis is only as good as the list of variables. Leave out $g$ in the pendulum and no group containing $T_p$ can be formed — the method then "proves" that the period is meaningless. Worse, leave out a variable whose dimensions happen to be expressible through the others (say, the wire radius $a$ in the dipole problem, which has the same dimension as $\ell$) and the method silently returns a wrong, incomplete answer: the true result would contain the additional group $a/\ell$. Physical insight decides what goes in the list; the algebra only processes it.

> [!warning] Common Pitfall — dimensionless is not the same as negligible
> Angles, ratios such as $a/\ell$, and numbers like $\mathrm{Pr}$ are dimensionless but may matter enormously. Dimensional analysis never tells you which groups can be dropped; that comes from scaling (Section 25.4) or experiment.

> [!tip] Engineering Note — scaling laws in electronics
> Dennard's MOSFET scaling rules are dimensional reasoning applied to a device: scale all lengths and the supply voltage by $1/\kappa$ and doping by $\kappa$, so that every electric field — and therefore every dimensionless group in the device equations — is unchanged. Then current scales as $1/\kappa$, capacitance as $1/\kappa$, gate delay $CV/I$ as $1/\kappa$, and power density stays constant. When threshold-voltage and leakage effects (which introduce the thermal voltage $kT/q$, a quantity that does *not* scale) became important, the rules broke down — a dimensional argument explaining why the "free lunch" of scaling ended.

## 25.4 Scaling and nondimensionalization

### From dimensionless groups to dimensionless equations

Dimensional analysis (Section 25.3) works on a *list of variables* and tells us which combinations can appear. **Nondimensionalization** works on the *equations themselves*: every dependent and independent variable is divided by a characteristic scale, and the equation is rewritten in the new variables. The parameters then collapse into a handful of dimensionless groups, each of which can be read as a ratio of competing effects. When the scales are chosen well, the size of each group tells us which terms dominate and which can be neglected.

> [!definition] Characteristic scale, scaled variable
> For a variable $q$ with characteristic magnitude $q_c$ (a constant built from the model's parameters, inputs or geometry), the **scaled variable** is $q^* = q/q_c$. For time, $t^* = t/t_c$, and derivatives transform by the chain rule:
> $$
> \frac{d q}{dt} = \frac{q_c}{t_c}\frac{dq^*}{dt^*}, \qquad \frac{d^2q}{dt^2} = \frac{q_c}{t_c^2}\frac{d^2q^*}{dt^{*2}} .
> $$
> A scaling is **good** when, for the behaviour of interest, the scaled variables and their derivatives are of order one, so that the size of each dimensionless coefficient directly measures the importance of its term.

The procedure is mechanical:

1. Identify dependent variables, independent variables and parameters.
2. Introduce scaled variables with *unknown* scales $q_c$, $t_c$, $x_c$.
3. Substitute and divide through by the coefficient of one chosen term (usually the one that must always be present).
4. Choose the scales to make as many coefficients as possible equal to 1.
5. The coefficients that remain are the governing dimensionless groups.

The number of groups that survive equals the number predicted by the Pi theorem for the quantities in the equation — the two methods are consistent, but scaling also tells you *where* in the equation each group sits.

### The series RLC circuit: $\omega_n$ and $\zeta$

The capacitor voltage of a series RLC circuit driven by $v_s(t)$ obeys

$$
LC\frac{d^2v}{dt^2} + RC\frac{dv}{dt} + v = v_s(t).
$$

Let $v = V_0\,y$ (with $V_0$ the size of the input, e.g. the step height) and $t = t_c\,\tau$:

$$
\frac{LC}{t_c^2}\,y'' + \frac{RC}{t_c}\,y' + y = u(\tau), \qquad u = v_s/V_0,
$$

with primes denoting $d/d\tau$. Choosing $t_c$ so that the first coefficient is 1 gives $t_c = \sqrt{LC}$, and the middle coefficient becomes $RC/\sqrt{LC} = R\sqrt{C/L}$. Writing it as $2\zeta$:

> [!theorem] Standard form of a second-order system
> With $\omega_n = 1/\sqrt{LC}$, $\tau = \omega_n t$ and
> $$
> \zeta = \frac{R}{2}\sqrt{\frac{C}{L}} = \frac{1}{2Q},
> $$
> the series RLC circuit becomes
> $$
> y'' + 2\zeta y' + y = u(\tau).
> $$
> Every series RLC circuit — every mass–spring–damper, every second-order servo — with the same $\zeta$ has the *same* scaled response; $\omega_n$ only stretches the time axis and $V_0$ only stretches the voltage axis.

Three parameters $(R, L, C)$ have become one ($\zeta$) plus two scales. Equivalently, the physical step response is

$$
v(t) = V_0\,S_\zeta(\omega_n t),
$$

where $S_\zeta$ is the universal unit-step response of the standard form. For $0<\zeta<1$,

$$
S_\zeta(\tau) = 1 - e^{-\zeta\tau}\left[\cos\!\left(\sqrt{1-\zeta^2}\,\tau\right) + \frac{\zeta}{\sqrt{1-\zeta^2}}\sin\!\left(\sqrt{1-\zeta^2}\,\tau\right)\right],
$$

for $\zeta = 1$, $S_1(\tau) = 1 - (1+\tau)e^{-\tau}$, and for $\zeta>1$ it is the sum of two decaying real exponentials with rates $\zeta\pm\sqrt{\zeta^2-1}$.

```desmos-graph
left=0; right=15; bottom=-0.1; top=1.8
---
y=1-e^{-0.2x}(\cos(0.9798x)+0.2041\sin(0.9798x))|blue|0<x
y=1-e^{-0.5x}(\cos(0.8660x)+0.5774\sin(0.8660x))|green|0<x
y=1-(1+x)e^{-x}|red|0<x
y=1-(1.0774e^{-0.2679x}-0.0774e^{-3.7321x})|purple|0<x
y=1|grey|dashed
```

*Universal step responses $S_\zeta(\tau)$ of $y''+2\zeta y'+y=u$ for $\zeta = 0.2$ (blue), $0.5$ (green), $1$ (red, critically damped) and $2$ (purple, overdamped); the horizontal axis is $\tau=\omega_n t$.*

The dimensionless figures of merit follow directly: percent overshoot $\mathrm{PO} = 100\,e^{-\pi\zeta/\sqrt{1-\zeta^2}}$, peak time $\tau_p = \pi/\sqrt{1-\zeta^2}$, and 2 % settling time $\tau_s \approx 4/\zeta$. Each is a function of $\zeta$ only; physical times are recovered by dividing by $\omega_n$.

> [!example] Worked Example 25.12 — Reading the RLC circuit of Example 25.4 in scaled form
> **Problem.** For $R = 10\ \Omega$, $L = 1\ \mathrm{mH}$, $C = 1\ \mathrm{\mu F}$, find $\omega_n$, $\zeta$, $Q$, overshoot, peak time and settling time of the capacitor voltage after a 1 V step. What resistance gives critical damping?
>
> **Scales.** $\omega_n = 1/\sqrt{10^{-3}\times10^{-6}} = 1/\sqrt{10^{-9}} = 31\,623\ \mathrm{rad/s}$ ($f_n = 5.03\ \mathrm{kHz}$). $\zeta = \tfrac{10}{2}\sqrt{10^{-6}/10^{-3}} = 5\times0.031\,62 = 0.158$; $Q = 1/(2\zeta) = 3.16$.
>
> **Figures of merit.** $\sqrt{1-\zeta^2} = 0.9874$.
> $$
> \mathrm{PO} = 100\exp\!\left(-\frac{\pi\times0.158}{0.9874}\right) = 100\,e^{-0.503} = 60.5\,\% .
> $$
> Peak time $t_p = \pi/(\omega_n\sqrt{1-\zeta^2}) = \pi/31\,225 = 100.6\ \mathrm{\mu s}$; settling time $t_s\approx4/(\zeta\omega_n) = 4/5000 = 0.8\ \mathrm{ms}$.
>
> **Critical damping.** $\zeta = 1$ requires $R = 2\sqrt{L/C} = 2\sqrt{1000} = 63.2\ \Omega$. The quantity $\sqrt{L/C} = 31.6\ \Omega$ is the **characteristic impedance** of the resonator — the natural resistance scale of the circuit — and $\zeta$ is simply $R$ measured in units of $2\sqrt{L/C}$.
>
> **Check against Example 25.4.** There the eigenvalues were $-5000\pm j\,31\,225$: real part $-\zeta\omega_n = -0.158\times31\,623 = -5000$ ✓, imaginary part $\omega_n\sqrt{1-\zeta^2} = 31\,225$ ✓.

### Thermal models: lumped scaling

The sink equation of Example 25.2, $C_s\dot T_s = P - (T_s-T_a)/R_{sa}$, scales with temperature-rise scale $\Delta T_c = PR_{sa}$ and time scale $t_c = R_{sa}C_s$. Setting $\theta = (T_s-T_a)/(P_0R_{sa})$ and $\tau = t/(R_{sa}C_s)$,

$$
\frac{d\theta}{d\tau} + \theta = p(\tau), \qquad p = P/P_0 .
$$

There is *no* dimensionless parameter left: every first-order lag is the same system, and its only "parameters" are the two scales. This is why one plot of $1 - e^{-\tau}$ serves RC filters, thermometers and heat sinks alike.

### Distributed thermal models: Biot and Fourier numbers

The lumped model assumes the body is at a uniform temperature. Whether that is reasonable is itself answered by scaling. Consider a plate of half-thickness $\ell$, conductivity $k$, diffusivity $\alpha = k/(\rho c_p)$, initially at $T_0$, suddenly exposed on both faces to a fluid at $T_\infty$ with heat-transfer coefficient $h$. By symmetry, on $0\le x\le\ell$:

$$
\frac{\partial T}{\partial t} = \alpha\frac{\partial^2T}{\partial x^2}, \qquad \frac{\partial T}{\partial x}(0,t) = 0, \qquad -k\frac{\partial T}{\partial x}(\ell,t) = h\,\bigl(T(\ell,t)-T_\infty\bigr), \qquad T(x,0) = T_0 .
$$

Scale with $\xi = x/\ell$, $\theta = (T-T_\infty)/(T_0-T_\infty)$ and an as-yet-free time scale $t_c$:

$$
\frac{\ell^2}{\alpha t_c}\frac{\partial\theta}{\partial\tau} = \frac{\partial^2\theta}{\partial\xi^2}, \qquad \frac{\partial\theta}{\partial\xi}(1,\tau) = -\frac{h\ell}{k}\,\theta(1,\tau).
$$

Choosing $t_c = \ell^2/\alpha$ (the diffusion time) makes the scaled time the **Fourier number** and leaves a single group in the boundary condition, the **Biot number**:

$$
\mathrm{Fo} = \frac{\alpha t}{\ell^2}, \qquad \mathrm{Bi} = \frac{h\ell}{k} = \frac{\ell/k}{1/h} = \frac{\text{internal conduction resistance}}{\text{surface convection resistance}} .
$$

> [!theorem] The lumped-capacitance criterion
> If $\mathrm{Bi}\ll1$ (in practice $\mathrm{Bi}<0.1$), the temperature inside the body is uniform to within about 5 % of the overall temperature difference, and the distributed model reduces to the lumped first-order model
> $$
> \frac{d\theta}{dt} = -\frac{hA}{\rho c_p V}\,\theta, \qquad \tau_{\text{lumped}} = \frac{\rho c_p V}{hA} = \frac{\ell^2}{\alpha\,\mathrm{Bi}} \quad(\text{plate}: V/A = \ell).
> $$

*Sketch.* With $\mathrm{Bi}\to0$ the boundary condition says $\partial\theta/\partial\xi\approx0$ at both faces, so to leading order $\theta$ is independent of $\xi$. Integrating the PDE over $0\le\xi\le1$ and using the boundary conditions gives $d\bar\theta/d\mathrm{Fo} = -\mathrm{Bi}\,\theta(1)\approx-\mathrm{Bi}\,\bar\theta$, i.e. decay on the slow time $\mathrm{Fo}\sim1/\mathrm{Bi}$, which is the lumped model. The internal profile relaxes on the fast time $\mathrm{Fo}\sim1$; when $\mathrm{Bi}\ll1$ these two time scales are well separated. $\square$

> [!example] Worked Example 25.13 — Is the heat sink really lumped?
> **Problem.** The heat sink of Example 25.2 has a 5 mm thick aluminium base ($k = 200\ \mathrm{W/(m\,K)}$, $\alpha = 8.4\times10^{-5}\ \mathrm{m^2/s}$), cooled by natural convection, $h\approx 10\ \mathrm{W/(m^2K)}$ (effective, per unit fin area, $h_{\text{eff}}\approx 100\ \mathrm{W/(m^2K)}$ referred to the base area). A second design uses a 1.6 mm FR-4 board as the "sink" ($k = 0.3\ \mathrm{W/(m\,K)}$). Compute Biot numbers and diffusion times.
>
> **Aluminium.** Heated from one side, the relevant length is the full thickness $\ell = 5\ \mathrm{mm}$:
> $$
> \mathrm{Bi} = \frac{100\times0.005}{200} = 2.5\times10^{-3}\ll0.1,
> $$
> and the diffusion time $\ell^2/\alpha = 2.5\times10^{-5}/8.4\times10^{-5} = 0.30\ \mathrm{s}$, versus a lumped time constant of 450 s. Lumped model: excellent.
>
> **FR-4.** $\mathrm{Bi} = 100\times0.0016/0.3 = 0.53$. The lumped model is *not* valid: the board surface far from the component runs much cooler than the copper under it, and the heat spreads laterally over a distance comparable to the board thickness times $1/\sqrt{\mathrm{Bi}}$. A distributed (PDE or finite-element) model, or at least a multi-node RC ladder, is required.
>
> **Moral.** The *same equation* (first-order lag) is a superb model in one case and a poor one in the other; the dimensionless number tells which, before any simulation is run.

### Time-scale separation and small parameters

A good scaling often exposes a **small parameter** $\varepsilon\ll1$ multiplying one term. Neglecting that term gives a **reduced model**; the conditions under which this is legitimate are the subject of perturbation theory.

- If $\varepsilon$ multiplies a term that is *not* the highest derivative, setting $\varepsilon = 0$ is a **regular** perturbation: the solution changes smoothly, by $O(\varepsilon)$.
- If $\varepsilon$ multiplies the *highest* derivative, setting $\varepsilon=0$ lowers the order of the equation and some initial conditions can no longer be satisfied: a **singular** perturbation. The neglected dynamics act in a fast **boundary layer** of duration $O(\varepsilon)$ in scaled time, after which the reduced model is accurate.

The RLC circuit with very small inductance is the canonical example. Scaling time on the *slow* scale $t_c = RC$ instead of $\sqrt{LC}$:

$$
\varepsilon\,y'' + y' + y = u, \qquad \varepsilon = \frac{L/R}{RC} = \frac{L}{R^2C} = \frac{1}{4\zeta^2}.
$$

When $\zeta\gg1$, $\varepsilon\ll1$ and the circuit behaves as an RC circuit with time constant $RC$, except for a brief initial interval of order $L/R$ during which the inductor current builds up. The two time constants are approximately $RC$ and $L/R$, and their ratio is $\varepsilon$ — a direct measure of the **stiffness** that will matter for numerical solution (Section 25.9).

> [!example] Worked Example 25.14 — A nonlinear oscillator reduced to one parameter
> **Problem.** A negative-resistance LC oscillator (for example, a tunnel-diode or cross-coupled-pair oscillator) is modelled by an inductor $L$ and capacitor $C$ in parallel with a nonlinear element drawing current $i_N(v) = -av + bv^3$ ($a>0$ is the small-signal negative conductance, $b>0$ models gain compression). Nondimensionalize and predict amplitude and frequency for $L = 10\ \mathrm{\mu H}$, $C = 1\ \mathrm{nF}$, $a = 1\ \mathrm{mS}$, $b = 1\ \mathrm{mA/V^3}$.
>
> **Equations.** KCL at the top node with inductor current $i_L$: $C\dot v + i_L + i_N(v) = 0$, $L\,\dot i_L = v$. Differentiate the first and substitute the second:
> $$
> C\ddot v - a\left(1 - \frac{3b}{a}v^2\right)\dot v + \frac{v}{L} = 0 .
> $$
>
> **Scaling.** Let $t = \sqrt{LC}\,\tau$ and $v = V_c\,x$. Multiplying by $L$:
> $$
> x'' - a\sqrt{\frac{L}{C}}\left(1 - \frac{3bV_c^2}{a}x^2\right)x' + x = 0 .
> $$
> Choose $V_c = \sqrt{a/(3b)}$ to make the bracket $1-x^2$. The result is the **Van der Pol equation**
> $$
> x'' - \mu\,(1-x^2)\,x' + x = 0, \qquad \mu = a\sqrt{\frac{L}{C}} .
> $$
> Five physical parameters $(L,C,a,b)$ plus time and voltage units have collapsed into *one* number $\mu$, the ratio of the negative conductance to the tank's characteristic admittance $\sqrt{C/L}$.
>
> **Numbers.** $\sqrt{L/C} = \sqrt{10^{-5}/10^{-9}} = 100\ \Omega$, so $\mu = 10^{-3}\times100 = 0.1$. $V_c = \sqrt{10^{-3}/(3\times10^{-3})} = 0.577\ \mathrm{V}$. For small $\mu$ the Van der Pol limit cycle is nearly sinusoidal with amplitude $x\approx2$ and scaled angular frequency $\approx1$; numerical integration gives amplitude $2.000$ and period $1.0006\times2\pi$ for $\mu = 0.1$. Hence $v(t)\approx 1.155\cos(\omega_0 t)\ \mathrm{V}$ with $f_0 = 1/(2\pi\sqrt{LC}) = 1.59\ \mathrm{MHz}$.
>
> **The parameter's meaning.** For $\mu = 1$ the scaled amplitude is $2.01$ and the period $1.06\times2\pi$; for $\mu = 5$ the oscillation becomes a *relaxation* oscillation with period $1.85\times2\pi$ and sharp transitions — a square-ish rather than sinusoidal waveform. A designer wanting a clean sine wave therefore wants $\mu\ll1$, i.e. just enough negative conductance to sustain oscillation; this design rule is invisible in the five-parameter form and obvious in the one-parameter form.

```desmos-graph
left=-2; right=2; bottom=-3; top=3
---
y=-x+x^{3}|blue
y=-x|red|dashed
x=1.155|grey|dotted
x=-1.155|grey|dotted
```

*Nonlinear element current $i_N(v)=-av+bv^3$ in mA against $v$ in volts ($a=1$ mS, $b=1$ mA/V³): the small-signal slope $-a$ (red dashed) supplies energy near $v=0$, the cubic term absorbs it at large swing, and the dotted lines mark the predicted oscillation amplitude $2V_c=1.155$ V.*

> [!warning] Common Pitfall — scaling with the wrong scale
> A scaling is a *hypothesis* about which balance dominates. Scaling the RLC circuit with $t_c=\sqrt{LC}$ when $\zeta = 50$ produces a coefficient $2\zeta = 100$ in front of $y'$ and hides the fact that the response has two widely separated time scales; scaling with $t_c = RC$ instead reveals $\varepsilon = 1/(4\zeta^2) = 10^{-4}$. If a scaled solution produces derivatives far from order one, the chosen scales were wrong for that regime — rescale.

> [!tip] Engineering Note — normalized design charts
> Filter tables (Butterworth, Chebyshev), Smith charts, universal resonance curves and heat-sink design charts are all products of nondimensionalization. A Butterworth table lists poles for $\omega_c = 1\ \mathrm{rad/s}$ and $R = 1\ \Omega$; *frequency scaling* and *impedance scaling* then produce any physical design. Recognising that a design problem is a scaled version of a tabulated one is often the fastest route to an answer.

## 25.5 Linearization

### Why linearize

Almost every real device is nonlinear: diodes and transistors are exponential or square-law, magnetic cores saturate, pendulums follow $\sin\theta$, power is quadratic in voltage. Yet the overwhelming majority of engineering analysis — transfer functions, Bode plots, phasors, eigenvalue stability, Kalman filters, linear control — is linear. The bridge is **linearization**: replacing a nonlinear model by its first-order Taylor approximation about an **operating point**, and treating the deviations from that point as small signals. Small-signal analysis of amplifiers, the "AC analysis" of every circuit simulator, and the stability analysis of every power converter and control loop rest on it.

The price is a **regime of validity**: the linear model is accurate only for deviations small enough that the neglected second-order terms are negligible. An honest linear model states how small "small" is.

### Operating points and the Jacobian

> [!definition] Operating point and linearized model
> For the model $\dot{\mathbf{x}} = \mathbf{f}(\mathbf{x},\mathbf{u})$, $\mathbf{y} = \mathbf{h}(\mathbf{x},\mathbf{u})$, an **operating point** (equilibrium) is a pair $(\bar{\mathbf{x}},\bar{\mathbf{u}})$ with $\mathbf{f}(\bar{\mathbf{x}},\bar{\mathbf{u}}) = \mathbf{0}$; the corresponding output is $\bar{\mathbf{y}} = \mathbf{h}(\bar{\mathbf{x}},\bar{\mathbf{u}})$. Writing $\mathbf{x} = \bar{\mathbf{x}}+\delta\mathbf{x}$, $\mathbf{u} = \bar{\mathbf{u}}+\delta\mathbf{u}$, $\mathbf{y} = \bar{\mathbf{y}}+\delta\mathbf{y}$, the **linearized (small-signal) model** is
> $$
> \delta\dot{\mathbf{x}} = \mathbf{A}\,\delta\mathbf{x} + \mathbf{B}\,\delta\mathbf{u}, \qquad \delta\mathbf{y} = \mathbf{C}\,\delta\mathbf{x} + \mathbf{D}\,\delta\mathbf{u},
> $$
> with Jacobian matrices evaluated at the operating point:
> $$
> \mathbf{A} = \left.\frac{\partial\mathbf{f}}{\partial\mathbf{x}}\right|_{\bar{\mathbf{x}},\bar{\mathbf{u}}},\quad \mathbf{B} = \left.\frac{\partial\mathbf{f}}{\partial\mathbf{u}}\right|_{\bar{\mathbf{x}},\bar{\mathbf{u}}},\quad \mathbf{C} = \left.\frac{\partial\mathbf{h}}{\partial\mathbf{x}}\right|_{\bar{\mathbf{x}},\bar{\mathbf{u}}},\quad \mathbf{D} = \left.\frac{\partial\mathbf{h}}{\partial\mathbf{u}}\right|_{\bar{\mathbf{x}},\bar{\mathbf{u}}}.
> $$

*Derivation.* Taylor's theorem in several variables (Ch. 6) gives

$$
\mathbf{f}(\bar{\mathbf{x}}+\delta\mathbf{x},\bar{\mathbf{u}}+\delta\mathbf{u}) = \underbrace{\mathbf{f}(\bar{\mathbf{x}},\bar{\mathbf{u}})}_{=\mathbf{0}} + \mathbf{A}\,\delta\mathbf{x} + \mathbf{B}\,\delta\mathbf{u} + \mathbf{r}(\delta\mathbf{x},\delta\mathbf{u}),
$$

and since $\dot{\bar{\mathbf{x}}} = \mathbf{0}$, $\delta\dot{\mathbf{x}} = \dot{\mathbf{x}}$. Dropping the remainder $\mathbf{r}$ gives the linear model. If $\mathbf{f}$ is twice continuously differentiable, the remainder is bounded by

$$
\lVert\mathbf{r}\rVert \le \tfrac12 M\left(\lVert\delta\mathbf{x}\rVert + \lVert\delta\mathbf{u}\rVert\right)^2,
$$

where $M$ bounds the second derivatives of $\mathbf{f}$ near the operating point. The linear terms are of size $\lVert\mathbf{A}\rVert\lVert\delta\mathbf{x}\rVert$, so the *relative* error of linearization grows in proportion to the signal amplitude.

> [!theorem] Linearization and local stability (Lyapunov's indirect method)
> Let $\bar{\mathbf{x}}$ be an equilibrium of $\dot{\mathbf{x}} = \mathbf{f}(\mathbf{x})$ with $\mathbf{f}$ continuously differentiable and Jacobian $\mathbf{A}$.
> 1. If every eigenvalue of $\mathbf{A}$ has negative real part, $\bar{\mathbf{x}}$ is locally asymptotically stable for the nonlinear system.
> 2. If some eigenvalue has positive real part, $\bar{\mathbf{x}}$ is unstable for the nonlinear system.
> 3. If no eigenvalue has positive real part but some lie on the imaginary axis, linearization is **inconclusive**; the nonlinear terms decide.
>
> Moreover (Hartman–Grobman), when no eigenvalue lies on the imaginary axis (a *hyperbolic* equilibrium), the nonlinear trajectories near $\bar{\mathbf{x}}$ are a continuous deformation of the linear ones.

The theorem is the formal licence for small-signal stability analysis — and case 3 is its fine print.

### Small-signal model of a diode

> [!example] Worked Example 25.15 — The diode small-signal resistance and its range of validity
> **Problem.** A silicon diode obeys the Shockley equation $i = I_S\left(e^{v/(nV_T)}-1\right)$ with $n = 1$ and $V_T = kT/q = 25.85\ \mathrm{mV}$ at 300 K. It is biased at $I_D = 1\ \mathrm{mA}$. Find the small-signal resistance and determine how large a voltage excursion the linear model tolerates for 2 %, 10 % and 20 % accuracy in the current change.
>
> **Linearization.** At the bias point $\bar v = V_D$, with $I_D\gg I_S$,
> $$
> g_d = \left.\frac{di}{dv}\right|_{V_D} = \frac{I_S}{nV_T}e^{V_D/(nV_T)} \approx \frac{I_D}{nV_T}, \qquad r_d = \frac{1}{g_d} = \frac{nV_T}{I_D} = \frac{25.85\ \mathrm{mV}}{1\ \mathrm{mA}} = 25.85\ \Omega .
> $$
> Remarkably, $r_d$ depends only on the bias current, not on $I_S$ — an unknown, process-dependent parameter has dropped out.
>
> **Exact versus linear.** The exact current change for a voltage step $\delta v$ is $\delta i = I_D\left(e^{\delta v/V_T}-1\right)$; the linear prediction is $\delta i_{\text{lin}} = \delta v/r_d = I_D\,\delta v/V_T$. With $x = \delta v/V_T$, the relative error is
> $$
> \frac{\delta i_{\text{lin}} - \delta i}{\delta i} = \frac{x}{e^x-1} - 1 \approx -\frac{x}{2} \quad(\lvert x\rvert\ll1).
> $$
>
> | $\delta v$ | exact $\delta i$ (µA) | linear $\delta i$ (µA) | error |
> |---|---|---|---|
> | $+1\ \mathrm{mV}$ | 39.44 | 38.68 | −1.9 % |
> | $+5\ \mathrm{mV}$ | 213.4 | 193.4 | −9.4 % |
> | $+10\ \mathrm{mV}$ | 472.3 | 386.8 | −18.1 % |
> | $-10\ \mathrm{mV}$ | −320.8 | −386.8 | +20.6 % |
> | $+20\ \mathrm{mV}$ | 1167.6 | 773.6 | −33.7 % |
>
> **Regime of validity.** Since the error is about $\lvert\delta v\rvert/(2V_T)$, an accuracy $e$ requires $\lvert\delta v\rvert\lesssim 2eV_T$: about $1\ \mathrm{mV}$ for 2 %, $5\ \mathrm{mV}$ for 10 %, $10\ \mathrm{mV}$ for 20 %. The often-quoted rule "small signal means $\lvert v_d\rvert\ll V_T$" is this calculation in words. Note the asymmetry between $\pm10\ \mathrm{mV}$: the linear model *underestimates* forward excursions and *overestimates* reverse ones, which is the origin of even-order distortion.

The second-order Taylor term gives that distortion quantitatively. For $\delta v = \hat V\cos\omega t$,

$$
\delta i \approx I_D\left[\frac{\hat V}{V_T}\cos\omega t + \frac12\left(\frac{\hat V}{V_T}\right)^2\cos^2\omega t\right] = I_D\left[\frac{\hat V}{V_T}\cos\omega t + \frac{\hat V^2}{4V_T^2} + \frac{\hat V^2}{4V_T^2}\cos2\omega t\right],
$$

so the second-harmonic distortion is $\mathrm{HD}_2 \approx \hat V/(4V_T)$: $0.97\,\%$ at $\hat V = 1\ \mathrm{mV}$ and $4.8\,\%$ at $5\ \mathrm{mV}$ (the exact Bessel-function result, $I_2(\hat V/V_T)/I_1(\hat V/V_T)$, agrees to three digits). The constant term $\hat V^2/(4V_T^2)$ is a shift of the operating point — **rectification** — and is exactly what the linear model misses entirely.

```desmos-graph
left=-60; right=60; bottom=-1; top=4
---
y=e^{x/25.85}-1|blue
y=x/25.85|red|dashed
y=0|grey
```

*Diode current change $\delta i/I_D=e^{\delta v/V_T}-1$ (blue) and its tangent-line linearization $\delta v/V_T$ (red dashed), with $\delta v$ in mV; they agree only within a few millivolts of the operating point.*

### The nonlinear pendulum (and every oscillator with a $\sin$)

The pendulum $\ddot\theta + \omega_0^2\sin\theta = 0$, $\omega_0^2 = g/\ell$, is the prototype for a whole class of EE models: the swing equation of a synchronous generator ($M\ddot\delta + D\dot\delta = P_m - P_{\max}\sin\delta$), the phase-locked loop with a sinusoidal phase detector, and the Josephson junction all have the same $\sin$ nonlinearity.

With $\mathbf{x} = [\theta,\ \dot\theta]^\mathsf{T}$, $\mathbf{f} = [\dot\theta,\ -\omega_0^2\sin\theta]^\mathsf{T}$, equilibria satisfy $\sin\bar\theta = 0$: $\bar\theta = 0$ (hanging) and $\bar\theta = \pi$ (inverted). The Jacobian is

$$
\mathbf{A}(\bar\theta) = \begin{bmatrix} 0 & 1\\ -\omega_0^2\cos\bar\theta & 0\end{bmatrix}.
$$

> [!example] Worked Example 25.16 — Linearizing the pendulum about both equilibria
> **Hanging, $\bar\theta = 0$.** $\mathbf{A} = \begin{bmatrix}0&1\\-\omega_0^2&0\end{bmatrix}$, eigenvalues $\pm j\omega_0$. Linearization predicts a centre (undamped oscillation at $\omega_0$) — but by the theorem this is case 3, *inconclusive*. Here energy conservation, $E = \tfrac12\dot\theta^2 + \omega_0^2(1-\cos\theta)$, settles it: level curves of $E$ near the bottom are closed, so the nonlinear system also has a centre.
>
> **Inverted, $\bar\theta = \pi$.** $\cos\pi = -1$, so $\mathbf{A} = \begin{bmatrix}0&1\\\omega_0^2&0\end{bmatrix}$, eigenvalues $\pm\omega_0$: a saddle, hyperbolic, so the nonlinear equilibrium is genuinely unstable, with small deviations growing as $e^{\omega_0 t}$. For $\ell = 1\ \mathrm{m}$, $\omega_0 = 3.13\ \mathrm{s^{-1}}$: an error doubles every $\ln2/\omega_0 = 0.22\ \mathrm{s}$, which sets the bandwidth an inverted-pendulum balancing controller must have.
>
> **How good is the hanging linearization?** The linear model predicts a period $T_0 = 2\pi/\omega_0$ independent of amplitude. The exact period is $T = (2/\pi)\,K\!\left(\sin(\theta_0/2)\right)T_0$, with $K$ the complete elliptic integral of the first kind; its series is $T/T_0 = 1 + \theta_0^2/16 + \dots$
>
> | Amplitude $\theta_0$ | exact $T/T_0$ | $1+\theta_0^2/16$ | linear-model error |
> |---|---|---|---|
> | 10° | 1.0019 | 1.0019 | 0.19 % |
> | 20° | 1.0077 | 1.0076 | 0.77 % |
> | 30° | 1.0174 | 1.0171 | 1.7 % |
> | 45° | 1.0400 | 1.0386 | 4.0 % |
> | 90° | 1.1803 | 1.1542 | 18 % |
> | 150° | 1.7622 | 1.4284 | 76 % |
>
> The "small-angle approximation" is good to 1 % up to about 23° and fails completely as $\theta_0\to180°$, where the period diverges.

```desmos-graph
left=-3.5; right=3.5; bottom=-1.5; top=1.5
---
y=\sin(x)|blue
y=x|red|dashed
y=-(x-3.14159)|green|dashed
```

*The pendulum nonlinearity $\sin\theta$ (blue) and its linearizations about $\bar\theta=0$ (red, slope $+1$) and $\bar\theta=\pi$ (green, slope $-1$): the sign of the slope decides between oscillation and exponential divergence.*

### A multivariable example: magnetic levitation

> [!example] Worked Example 25.17 — Linearizing an electromagnet suspension
> **Problem.** A steel ball of mass $m$ hangs below an electromagnet. With gap $x$ (measured downward from the magnet face) and coil current $i$, the magnetic force is approximately $k\,i^2/x^2$ upward, so
> $$
> m\ddot x = mg - k\frac{i^2}{x^2}.
> $$
> For $m = 0.05\ \mathrm{kg}$, $k = 10^{-4}\ \mathrm{N\,m^2/A^2}$, find the current that holds the ball at $x_0 = 10\ \mathrm{mm}$, linearize, and classify the equilibrium.
>
> **Operating point.** $mg = k i_0^2/x_0^2 \Rightarrow i_0 = x_0\sqrt{mg/k} = 0.01\sqrt{0.4905/10^{-4}} = 0.700\ \mathrm{A}$.
>
> **Jacobian.** With $\mathbf{x} = [x,\ \dot x]^\mathsf{T}$ and input $i$, $f_2 = g - (k/m)\,i^2/x^2$:
> $$
> \frac{\partial f_2}{\partial x} = \frac{2k i_0^2}{m x_0^3} = \frac{2g}{x_0}, \qquad \frac{\partial f_2}{\partial i} = -\frac{2k i_0}{m x_0^2} = -\frac{2g}{i_0},
> $$
> where the equilibrium condition was used to simplify. So
> $$
> \delta\ddot x = \frac{2g}{x_0}\,\delta x - \frac{2g}{i_0}\,\delta i = 1962\,\delta x - 28.0\,\delta i .
> $$
>
> **Classification.** Open loop ($\delta i = 0$), $\lambda = \pm\sqrt{2g/x_0} = \pm44.3\ \mathrm{s^{-1}}$: an unstable saddle with time constant $1/44.3 = 23\ \mathrm{ms}$. The linearized model, with transfer function $\delta X/\delta I = -28.0/(s^2 - 1962)$, is the starting point for designing the stabilizing feedback loop.
>
> **Validity.** The second-order term in $\delta x$ is $\tfrac12(\partial^2f_2/\partial x^2)\delta x^2 = -(3g/x_0^2)\,\delta x^2$; compared with the linear term $(2g/x_0)\delta x$, the ratio is $1.5\,\lvert\delta x\rvert/x_0$. For 10 % accuracy the ball must stay within about $\pm0.7\ \mathrm{mm}$ of the set point — the linear controller design is valid only if the controller keeps it there.

### Beyond equilibria

Two extensions are worth knowing by name:

- **Linearization about a trajectory.** If $\bar{\mathbf{x}}(t)$ is a nominal solution (e.g. a switching converter's periodic steady state), the deviations obey a *time-varying* linear system $\delta\dot{\mathbf{x}} = \mathbf{A}(t)\delta\mathbf{x} + \mathbf{B}(t)\delta\mathbf{u}$, with Jacobians evaluated along $\bar{\mathbf{x}}(t)$. For periodic trajectories this leads to Floquet theory and to state-space averaging in power electronics.
- **Quasi-linearization (describing functions).** For a sinusoidal signal of amplitude $A$ through a static nonlinearity, replace the nonlinearity by the gain of its fundamental Fourier component (Ch. 15), $N(A)$. This amplitude-dependent "linearization" predicts limit cycles, as in the Van der Pol oscillator, where it gives amplitude $x = 2$ exactly.

> [!warning] Common Pitfall — linearizing what is not differentiable
> Saturation, dead zones, ideal switches, Coulomb friction $\tau_c\operatorname{sgn}(\omega)$ and absolute-value detectors have no derivative at the critical point. A Jacobian taken "on one side" describes only that side, and small-signal analysis at the corner is meaningless: an amplifier biased exactly at the edge of saturation has *no* well-defined small-signal gain. Either bias away from the corner or use describing-function, piecewise-linear or simulation methods.

> [!tip] Engineering Note — what SPICE "AC analysis" actually does
> A circuit simulator's AC analysis first solves the nonlinear DC operating-point equations by Newton's method (Ch. 22), then replaces every nonlinear device by its Jacobian — $g_m$, $r_\pi$, $r_o$, junction capacitances evaluated at the bias point — and solves the resulting *linear* phasor equations at each frequency. An AC analysis therefore reports a gain of, say, 40 dB for a 10 V input signal without complaint; it is the user's responsibility to remember that the result is valid only for signals small enough to keep the operating point fixed.

## 25.6 Model assumptions

### Assumptions are part of the model

Every model in this chapter so far rested on assumptions: the capacitor had no leakage (Example 25.1), the package stored no heat (25.2), friction was viscous (25.3), the source had zero resistance (25.5), the diode signal was small (25.15). An assumption is not a weakness to be hidden; it is the *specification* of the model's domain of validity. A model without stated assumptions cannot be checked, cannot be trusted outside the conditions it was tested in, and cannot be improved systematically.

> [!definition] Assumption and regime of validity
> A **modeling assumption** is a statement that some effect, quantity or interaction is negligible, constant, ideal or of a particular form. Its **regime of validity** is the set of operating conditions — usually expressible as an inequality on a dimensionless group — under which the error it introduces is below the accuracy required by the model's purpose.

The key move is the second sentence: turning a qualitative assumption ("the wire is short") into a quantitative inequality ("$\ell<\lambda/10$"). Sections 25.3–25.5 provide the tools; this section collects the standard assumptions of EE modeling and the inequalities that police them.

### A catalogue of standard assumptions

| Assumption | What it buys | Regime of validity (typical) | When it fails |
|---|---|---|---|
| Lumped (quasi-static) circuit | ODEs instead of PDEs; Kirchhoff's laws | size $\ell\ll\lambda$, say $\ell<\lambda/10$; or edge rise time $t_r > 6\,t_{pd}$ | long cables, fast digital edges, RF |
| Lumped thermal body | one temperature per body | $\mathrm{Bi} = h\ell/k < 0.1$ | thin, poorly conducting boards (Ex. 25.13) |
| Linearity | superposition, transforms, eigenvalues | deviations small vs. nonlinearity scale, e.g. $\lvert v_d\rvert\ll V_T$ | large signals, saturation, clipping |
| Time invariance | constant coefficients, transfer functions | parameters drift less than accuracy over the window | warm-up, ageing, battery discharge |
| Ideal op-amp | virtual short, infinite input impedance | loop gain $\lvert A\beta\rvert\gg1$, i.e. $f\ll f_{\text{GBW}}\beta$ | high frequency, high closed-loop gain |
| Ideal capacitor | $Z = 1/(j\omega C)$ | $f$ well below self-resonance $1/(2\pi\sqrt{L_{\text{ESL}}C})$ | decoupling at hundreds of MHz |
| Constant resistance | $R$ fixed | $\alpha\,\Delta T\ll$ accuracy | windings, filaments, current shunts |
| Deterministic | no noise terms | signal ≫ noise (SNR large) | low-level sensing, receivers |
| Neglect of a time constant | lower-order model | that time constant ≪ time scale of interest (Sec. 25.4) | stiff, closely spaced modes |
| Continuum / Ohm's law | $\mathbf{J} = \sigma\mathbf{E}$ | dimensions ≫ mean free path; fields moderate | nanoscale devices, breakdown |

### Making an assumption quantitative

The pattern is always the same: identify the neglected term, estimate its size relative to the retained terms using the scales of Section 25.4, and require the ratio to be below the tolerated error.

> [!example] Worked Example 25.18 — Is a PCB trace a lumped wire or a transmission line?
> **Problem.** A 10 cm trace on FR-4 (effective relative permittivity $\varepsilon_{\text{eff}}\approx 4$) carries (a) a 100 MHz sinusoid, (b) a digital signal with 0.5 ns rise time. Is the lumped assumption valid?
>
> **Propagation.** $v_p = c/\sqrt{\varepsilon_{\text{eff}}} = 1.5\times10^8\ \mathrm{m/s}$, so the one-way delay is $t_{pd} = 0.1/1.5\times10^8 = 0.67\ \mathrm{ns}$.
>
> **(a) Sinusoid.** $\lambda = v_p/f = 1.5\ \mathrm{m}$, $\ell/\lambda = 0.067 < 0.1$. The phase shift along the trace is $360°\times0.067 = 24°$ — borderline. For a rough amplitude estimate the lumped model is acceptable; for a phase-critical path it is not.
>
> **(b) Digital edge.** The relevant comparison is rise time versus delay. With $t_r = 0.5\ \mathrm{ns} < t_{pd} = 0.67\ \mathrm{ns}$, the edge has fully risen at the driver before it reaches the load: reflections will be seen as distinct steps. Equivalently, the signal's "knee" frequency $0.35/t_r = 700\ \mathrm{MHz}$ has wavelength 0.21 m, and $\ell/\lambda\approx0.5$. The trace must be modelled as a transmission line (Ch. 13, telegrapher's equations) and terminated.
>
> **Moral.** "Is the circuit lumped?" has no answer until the *signal* is specified. The assumption belongs to the combination of system and input.

> [!example] Worked Example 25.19 — When is an op-amp ideal?
> **Problem.** A non-inverting amplifier with ideal gain $1/\beta = 100$ uses an op-amp with DC gain $A_0 = 10^5$ and a single pole giving gain–bandwidth product $f_{\text{GBW}} = 1\ \mathrm{MHz}$ (pole at 10 Hz). Over what frequency range is the ideal-op-amp model (gain exactly 100) accurate to 1 %? To 3 dB?
>
> **Model.** The closed-loop gain is
> $$
> G = \frac{A}{1 + A\beta} = \frac{1}{\beta}\cdot\frac{1}{1 + 1/(A\beta)}, \qquad A(f) = \frac{A_0}{1 + jf/10\ \mathrm{Hz}} .
> $$
> The ideal model sets $1/(A\beta) = 0$; the neglected term is the reciprocal of the loop gain.
>
> **DC.** $A_0\beta = 1000$, error $\approx 0.1\,\%$: excellent.
>
> **Frequency dependence.** For $f\gg10\ \mathrm{Hz}$, $A\approx -j f_{\text{GBW}}/f$, so $1/(A\beta) \approx j f/(f_{\text{GBW}}\beta) = jf/10\ \mathrm{kHz}$, and
> $$
> \lvert G\rvert \approx \frac{100}{\sqrt{1 + (f/10\ \mathrm{kHz})^2}} .
> $$
> A 1 % magnitude error requires $(f/10\ \mathrm{kHz})^2 < 0.02$, i.e. $f < 1.4\ \mathrm{kHz}$. The 3 dB point is $f = f_{\text{GBW}}\beta = 10\ \mathrm{kHz}$. A direct computation gives $\lvert G\rvert = 99.41$ at 1 kHz, $89.4$ at 5 kHz and $70.7$ (−3.0 dB) at 10 kHz, confirming the estimate.
>
> **Moral.** The "ideal op-amp" assumption is not a property of the op-amp; it is a condition, $\lvert A\beta\rvert\gg1$, on the op-amp *in this circuit at this frequency*. An op-amp that is ideal for an audio buffer is far from ideal at gain 100 in the same band.

> [!example] Worked Example 25.20 — Hidden parasitics: the capacitor that becomes an inductor
> **Problem.** A 100 nF ceramic decoupling capacitor has equivalent series inductance $L_{\text{ESL}} = 1\ \mathrm{nH}$ and equivalent series resistance $R_{\text{ESR}} = 10\ \mathrm{m\Omega}$. Over what frequency range is it "a capacitor"?
>
> **Model.** $Z(\omega) = R_{\text{ESR}} + j\left(\omega L_{\text{ESL}} - \dfrac{1}{\omega C}\right)$, a series RLC. Its self-resonant frequency is
> $$
> f_{\text{SRF}} = \frac{1}{2\pi\sqrt{L_{\text{ESL}}C}} = \frac{1}{2\pi\sqrt{10^{-9}\times10^{-7}}} = 15.9\ \mathrm{MHz}.
> $$
> At 1 MHz, $\lvert Z\rvert = 1.585\ \Omega$ versus the ideal $1/(\omega C) = 1.592\ \Omega$ — the ideal model is within 0.4 %. At 15.9 MHz the reactances cancel and $\lvert Z\rvert = R_{\text{ESR}} = 10\ \mathrm{m\Omega}$. At 100 MHz, $\lvert Z\rvert = 0.61\ \Omega$ and the part is *inductive*: an ideal-capacitor model would predict $0.016\ \Omega$, wrong by a factor of 38.
>
> **Regime.** The ideal model holds for $f\ll f_{\text{SRF}}$; the inequality comes from the ratio of the neglected to the retained reactance, $\omega^2 L_{\text{ESL}}C = (f/f_{\text{SRF}})^2$.

### Checking assumptions a posteriori

Many assumptions can only be checked after the model has been solved, because the size of the neglected term depends on the solution. The discipline is:

1. Solve the simplified model.
2. Use the solution to *evaluate* every neglected term.
3. Confirm each is small compared with the retained terms, by the margin the purpose requires.
4. If not, restore the term and re-solve.

> [!example] Worked Example 25.21 — A posteriori check of a constant-resistance assumption
> **Problem.** In the DC motor model of Example 25.3, $R_a = 2\ \Omega$ was taken constant. The copper winding has temperature coefficient $\alpha = 0.003\,93\ \mathrm{K^{-1}}$. A thermal model of the motor shows that at rated load the winding settles $40\ \mathrm{K}$ above the temperature at which $R_a$ was measured. The speed model must be accurate to 5 %. Is the assumption acceptable?
>
> **Evaluate the neglected effect.** $R_a(T) = R_a(1 + \alpha\Delta T) = 2(1 + 0.157) = 2.31\ \Omega$ — a 15.7 % change.
>
> **Effect on the output.** At load torque $\tau_L$ the steady speed (with $b$ small) is $\omega \approx v_a/K - R_a\tau_L/K^2$. The *speed droop* term is proportional to $R_a$, so it is underestimated by 15.7 %. If the droop is 10 % of no-load speed, the speed error is $1.6\,\%$: acceptable. If the motor runs at stall-like loads where droop is 50 %, the error is 7.9 %: not acceptable, and $R_a$ must become temperature-dependent — making temperature a *state* coupled to an electrical model (Section 25.2).
>
> **Moral.** Whether an assumption is valid depends on how sensitive the *output of interest* is to it, which is a question about the model's sensitivity, not only about the size of the neglected effect.

> [!tip] Engineering Note — write the assumptions down
> Professional models (a thermal design report, a link budget, a SPICE netlist submitted for review) should begin with a numbered list of assumptions, each with its validity criterion and, ideally, the computed margin ("lumped: $\ell/\lambda = 0.03$ at the highest frequency of interest"). When the design changes — a faster edge, a hotter enclosure — the list tells you immediately which parts of the model must be revisited.

> [!warning] Common Pitfall — validity by habit
> Assumptions inherited from a textbook or a previous project ("op-amps are ideal", "wires are short", "noise is white", "the load is resistive") are the most dangerous, because they are never questioned. The circuit that worked at 1 MHz does not work at 100 MHz for precisely this reason. Each assumption should be re-justified for every new operating range.

### Conservative and non-conservative assumptions

Some assumptions err in a known direction. Neglecting a heat path *over*-predicts temperature (conservative for reliability); neglecting a parasitic capacitance may *under*-predict oscillation risk (non-conservative). When an assumption cannot be justified precisely, choosing the conservative direction lets the model still give a safe answer — a bound rather than an estimate. Worst-case and corner analysis (setting each uncertain parameter at the extreme that hurts most) is the systematic version of this idea; Monte Carlo analysis (Ch. 24) is the statistical version.

## 25.7 Model validation

### Verification, calibration, validation

Three activities are often confused:

> [!definition] Verification, calibration and validation
> - **Verification** asks *"Are we solving the equations right?"* — is the code or derivation a correct solution of the stated mathematical model? (Checked against analytical solutions, convergence tests, conservation checks.)
> - **Calibration** (parameter identification, Section 25.8) asks *"What parameter values make the model match these data?"*
> - **Validation** asks *"Are we solving the right equations?"* — does the calibrated model reproduce the behaviour of the real system, to the accuracy required, on data that were **not** used for calibration?

A model can be perfectly verified and perfectly calibrated and still be invalid: a first-order model fitted by least squares to a second-order system will be solved exactly and fitted optimally, and still be wrong. Validation is the step that catches structural error, and it requires data the model has not seen.

### Residuals

> [!definition] Residual
> Given measurements $y_k$ ($k = 1,\dots,N$) and model predictions $\hat y_k$ at the same conditions, the **residuals** are
> $$
> e_k = y_k - \hat y_k .
> $$
> If the model is structurally correct and the measurements contain independent zero-mean noise of standard deviation $\sigma$, the residuals should look like that noise: zero mean, standard deviation $\approx\sigma$, no trend against time, input or predicted value, and no correlation between neighbours.

Residual analysis is the single most informative validation tool, because *structure in the residuals points to the missing physics*. A bowed residual plot means a missing curvature (a quadratic term, a saturation); residuals growing with $\hat y$ mean multiplicative rather than additive error; residuals that oscillate in time mean an unmodelled mode; a step in the residuals means an unmodelled event.

### Goodness-of-fit measures

For a model with $p$ fitted parameters and $N$ data points, with $\bar y = \frac1N\sum y_k$:

| Measure | Definition | Interpretation |
|---|---|---|
| Sum of squared errors | $\mathrm{SSE} = \sum e_k^2$ | raw misfit; depends on units and $N$ |
| Root-mean-square error | $\mathrm{RMSE} = \sqrt{\mathrm{SSE}/N}$ | typical error, in units of $y$ |
| Residual standard error | $s = \sqrt{\mathrm{SSE}/(N-p)}$ | unbiased estimate of $\sigma$ (Ch. 20) |
| Coefficient of determination | $R^2 = 1 - \mathrm{SSE}/\mathrm{SST}$, $\ \mathrm{SST} = \sum(y_k-\bar y)^2$ | fraction of variance "explained" |
| Adjusted $R^2$ | $1 - \dfrac{\mathrm{SSE}/(N-p)}{\mathrm{SST}/(N-1)}$ | penalises extra parameters |
| Chi-square | $\chi^2 = \sum (e_k/\sigma_k)^2$ | with known noise $\sigma_k$ |
| Reduced chi-square | $\chi^2_\nu = \chi^2/\nu$, $\ \nu = N - p$ | $\approx1$ for an adequate model |
| Normalized fit | $100\left(1 - \dfrac{\lVert\mathbf{y}-\hat{\mathbf{y}}\rVert}{\lVert\mathbf{y}-\bar y\mathbf{1}\rVert}\right)\%$ | used in system identification |

> [!theorem] Reduced chi-square test
> If the model is correct, the noise is Gaussian with known standard deviations $\sigma_k$, and $p$ linear parameters are fitted, then $\chi^2$ follows a chi-square distribution with $\nu = N-p$ degrees of freedom (Ch. 20). Hence $E[\chi^2_\nu] = 1$ and $\operatorname{Var}(\chi^2_\nu) = 2/\nu$. A value of $\chi^2_\nu$ far above $1 + 2\sqrt{2/\nu}$ indicates model inadequacy (or underestimated noise); a value far below 1 indicates overestimated noise or overfitting.

> [!warning] Common Pitfall — a high $R^2$ does not validate a model
> $R^2$ measures how much of the *variation in the data* the model reproduces. When the data span a wide range (a sensor calibrated from 0 to 400 °C), almost any monotonic model achieves $R^2 > 0.999$, while its errors can still be fifty times the measurement noise. $R^2$ cannot detect structural error; residual plots and $\chi^2_\nu$ can.

> [!example] Worked Example 25.22 — Validating a linear RTD model
> **Problem.** A platinum resistance thermometer is calibrated at nine temperatures; the resistance meter has standard uncertainty $\sigma = 0.02\ \Omega$.
>
> | $T$ (°C) | 0 | 50 | 100 | 150 | 200 | 250 | 300 | 350 | 400 |
> |---|---|---|---|---|---|---|---|---|---|
> | $R$ (Ω) | 100.04 | 119.41 | 138.51 | 157.29 | 175.85 | 194.09 | 212.05 | 229.70 | 247.09 |
>
> Fit (a) a linear model $R = c_0 + c_1T$ and (b) a quadratic model $R = c_0 + c_1T + c_2T^2$. Validate each.
>
> **(a) Linear.** Least squares (Section 25.8) gives $c_0 = 101.362\ \Omega$, $c_1 = 0.367\,65\ \Omega/^\circ\mathrm{C}$. Then
> $$
> \mathrm{SSE} = 6.242\ \Omega^2, \quad R^2 = 0.999\,69, \quad \mathrm{RMSE} = 0.833\ \Omega, \quad \chi^2_\nu = \frac{6.242/0.02^2}{7} = 2229 .
> $$
> $R^2$ looks superb. But the residuals are
> $$
> -1.32,\ -0.34,\ +0.38,\ +0.78,\ +0.96,\ +0.82,\ +0.39,\ -0.34,\ -1.33\ \ \Omega,
> $$
> a symmetric bow, up to sixty-five times the noise level, and $\chi^2_\nu\approx2229$ is absurdly larger than 1. The linear model is **invalid** at this accuracy: it is missing curvature. As a temperature sensor it would err by up to $1.33/0.368\approx3.6\ ^\circ\mathrm{C}$.
>
> **(b) Quadratic.** $c_0 = 100.034\ \Omega$, $c_1 = 0.390\,43\ \Omega/^\circ\mathrm{C}$, $c_2 = -5.694\times10^{-5}\ \Omega/^\circ\mathrm{C}^2$, with
> $$
> \mathrm{SSE} = 1.15\times10^{-3}\ \Omega^2, \quad \mathrm{RMSE} = 0.011\ \Omega, \quad \chi^2_\nu = \frac{2.87}{6} = 0.48 .
> $$
> Residuals ($+0.006, -0.003, +0.003, -0.026, +0.009, +0.009, +0.013, -0.008, -0.004$ Ω) show no pattern and are at the noise level. With $\nu = 6$, the 95 % range of $\chi^2_\nu$ is roughly $0.21$ to $2.41$, so $0.48$ is consistent with a correct model (the data were also rounded to 0.01 Ω, which slightly reduces the apparent noise).
>
> **Physics check.** The fitted $c_1/c_0 = 3.903\times10^{-3}$ and $c_2/c_0 = -5.69\times10^{-7}$ per °C² are close to the Callendar–Van Dusen standard coefficients $A = 3.9083\times10^{-3}$, $B = -5.775\times10^{-7}$ — independent confirmation that the quadratic structure is physical, not just a better curve.

### Residual tests for dynamic models

For time-series data, residuals must also be **uncorrelated in time**. The sample autocorrelation of the residuals,

$$
\hat r(k) = \frac{\sum_{n=1}^{N-k}(e_n-\bar e)(e_{n+k}-\bar e)}{\sum_{n=1}^{N}(e_n-\bar e)^2},
$$

is approximately $\mathcal{N}(0, 1/N)$ for each $k\ge1$ if the residuals are white (Ch. 21). Values outside $\pm1.96/\sqrt N$ at more than about 5 % of lags indicate unmodelled dynamics. A complementary test is the cross-correlation between residuals and past inputs, which should also vanish if the model has captured the input–output dynamics.

> [!example] Worked Example 25.23 — Detecting a missing time constant
> **Problem.** A sensor's step response is sampled at 251 points over 5 s, with noise of standard deviation 0.01. Unknown to the modeller, the true response is second order with time constants 1.0 s and 0.2 s. Fit a first-order model $K(1-e^{-t/\tau})$ and a second-order model, and compare residual autocorrelations.
>
> **First-order fit.** Nonlinear least squares gives $K = 1.037$, $\tau = 1.334\ \mathrm{s}$, RMSE $= 0.030$ — only three times the noise, and a plot of model against data "looks good". But the residual autocorrelation is
> $$
> \hat r(1) = 0.88,\quad \hat r(2) = 0.88,\quad \hat r(5) = 0.84,\quad \hat r(10) = 0.75,
> $$
> against a whiteness bound $1.96/\sqrt{251} = 0.124$. The residuals are a slow, smooth curve — the signature of a missing mode. Note also that the fitted $\tau$ (1.33 s) is not either true time constant: an incorrect structure biases the parameters.
>
> **Second-order fit.** $K = 0.999$, $\tau_1 = 0.991\ \mathrm{s}$, $\tau_2 = 0.205\ \mathrm{s}$, RMSE $= 0.0096$ (matching the noise), $\hat r(1) = -0.06$, $\hat r(2) = -0.04$, $\hat r(5) = -0.08$, $\hat r(10) = 0.02$: all within the bound. The model is accepted, and its parameters are now close to the true ones.

### Overfitting and the need for independent data

Adding parameters always reduces SSE on the data used for fitting: a polynomial of degree $N-1$ passes through all $N$ points exactly. A model that fits the noise rather than the system is **overfitted**: it reproduces the calibration data beautifully and predicts new data badly. The expected prediction error of a model decomposes (for squared error) as

$$
E\left[(y_{\text{new}} - \hat y)^2\right] = \underbrace{\sigma^2}_{\text{irreducible noise}} + \underbrace{\bigl(\text{bias}\bigr)^2}_{\text{structural error}} + \underbrace{\text{variance}}_{\text{parameter uncertainty}},
$$

the **bias–variance trade-off**. Too simple a model has large bias; too complex a model has large variance because its many parameters are poorly determined by finite data. The goal is the minimum of the sum, and it can only be found by testing on data that were not used for fitting.

> [!definition] Cross-validation
> - **Hold-out validation:** split the data into a *training* set (for fitting) and a *validation* set (for evaluating prediction error).
> - **$k$-fold cross-validation:** split the data into $k$ parts; fit $k$ times, each time leaving one part out and computing the prediction error on it; average the $k$ errors.
> - **Leave-one-out cross-validation (LOOCV):** $k = N$. For linear least squares it needs no refitting: with hat matrix $\mathbf{H} = \mathbf{X}(\mathbf{X}^\mathsf{T}\mathbf{X})^{-1}\mathbf{X}^\mathsf{T}$ and leverages $h_{kk}$, the leave-one-out residual is $e_k/(1-h_{kk})$, and the **PRESS** statistic is
> $$
> \mathrm{PRESS} = \sum_{k=1}^{N}\left(\frac{e_k}{1-h_{kk}}\right)^2 .
> $$

*Why the shortcut holds (sketch).* Removing observation $k$ changes the normal-equation matrix by a rank-one term $-\mathbf{x}_k\mathbf{x}_k^\mathsf{T}$; the Sherman–Morrison formula for the inverse of a rank-one update gives the new prediction at $\mathbf{x}_k$ in closed form, and simplifying yields $y_k - \hat y_{k,(-k)} = e_k/(1-h_{kk})$. Points with high leverage ($h_{kk}$ near 1, typically at the ends of the range) have their leave-one-out errors amplified the most, which is exactly where overfitted polynomials misbehave.

> [!example] Worked Example 25.24 — Choosing the polynomial order by cross-validation
> **Problem.** For the RTD data of Example 25.22, compare degrees 1, 2 and 3 by LOOCV.
>
> **Computation** (PRESS via leverages):
>
> | Degree | $p$ | Training RMSE (Ω) | LOOCV RMSE $\sqrt{\mathrm{PRESS}/N}$ (Ω) | AICc |
> |---|---|---|---|---|
> | 1 | 2 | 0.833 | 1.192 | — (structurally inadequate) |
> | 2 | 3 | 0.0113 | **0.0161** | −69.9 |
> | 3 | 4 | 0.0097 | 0.0175 | −65.5 |
>
> The training RMSE keeps falling as parameters are added, as it always must. The leave-one-out error is smallest for the quadratic: the cubic term reduces the training error by fitting noise, and that "improvement" does not survive prediction of withheld points. The small-sample Akaike criterion, $\mathrm{AICc} = N\ln(\mathrm{SSE}/N) + 2p + 2p(p+1)/(N-p-1)$, agrees. Both select the physically motivated quadratic model.

> [!warning] Common Pitfall — cross-validating time series by random shuffling
> For dynamic data, neighbouring samples are correlated. Randomly assigning samples to folds lets the model "see" the neighbours of every validation point, and the cross-validation error is optimistically biased. Use **blocked** folds (contiguous segments of time), or better, validate on an entirely separate experiment with a different input signal.

> [!tip] Engineering Note — validate against the purpose
> A model is validated *for a purpose*. A thermal model intended to predict peak junction temperature should be validated on the peak, under the worst-case load profile, not on the average temperature under a gentle one. Choose validation experiments that stress the conditions where the model will be used, and especially those where the assumptions of Section 25.6 are closest to failing.

## 25.8 Parameter identification

### The identification problem

A grey-box model has a structure from physics and parameters from data. **Parameter identification** (estimation, calibration, fitting) is the inverse problem: given input–output measurements, find the parameter vector $\mathbf{p}$ that makes the model reproduce them.

> [!definition] Output-error parameter estimation
> Given data $\{(u_k, y_k)\}_{k=1}^N$ and a model that predicts $\hat y_k(\mathbf{p})$, the **least-squares estimate** is
> $$
> \hat{\mathbf{p}} = \arg\min_{\mathbf{p}} J(\mathbf{p}), \qquad J(\mathbf{p}) = \sum_{k=1}^N w_k\bigl(y_k - \hat y_k(\mathbf{p})\bigr)^2,
> $$
> with weights $w_k = 1/\sigma_k^2$ when measurement uncertainties differ. If the errors are independent Gaussian, this is also the maximum-likelihood estimate (Ch. 20).

Three families of method are used in practice, in increasing generality: graphical/characteristic-point methods (read parameters off a step response), linear least squares (when the model is linear in the parameters) and nonlinear least squares (everything else). All of them must be followed by an identifiability check: are the data actually capable of determining the parameters?

### Identification from a step response

The step response is the most common identification experiment in EE: switch on a heater, a motor voltage, a reference input, and record the output. Standard model structures have parameters that can be read from characteristic features.

> [!theorem] First-order step response
> For $\tau\dot y + y = Ku$ with a step of size $\Delta u$ at $t=0$ from equilibrium, $y(t) - y(0) = K\Delta u\left(1-e^{-t/\tau}\right)$. Hence
> - **gain:** $K = \Delta y_\infty/\Delta u$;
> - **time constant:** $\tau$ is the time to reach $1 - e^{-1} = 63.2\,\%$ of the final change; equivalently, the initial slope is $K\Delta u/\tau$;
> - **log-linear form:** $\ln\!\left(1 - \dfrac{y(t)-y(0)}{\Delta y_\infty}\right) = -\dfrac{t}{\tau}$, a straight line through the origin with slope $-1/\tau$.

The log-linear form is preferable to reading a single point, because it uses all the data and its straightness is itself a validation test: curvature at early times reveals a second, faster time constant or a dead time. For a first-order-plus-dead-time (FOPDT) model $K e^{-\theta s}/(\tau s + 1)$, a classic two-point method uses the times $t_{28}$ and $t_{63}$ at which the response reaches 28.3 % and 63.2 %: $\tau = 1.5\,(t_{63} - t_{28})$, $\theta = t_{63} - \tau$ (these follow from $1 - e^{-1/3} = 0.283$ and $1 - e^{-1} = 0.632$).

> [!example] Worked Example 25.25 — Identifying a first-order system from step data
> **Problem.** The supply node of Example 25.1 is tested: after switch-on, its voltage settles at $v_\infty = 4.00\ \mathrm{V}$, and samples are
>
> | $t$ (s) | 0.02 | 0.05 | 0.10 | 0.15 | 0.20 | 0.30 |
> |---|---|---|---|---|---|---|
> | $v$ (V) | 0.743 | 1.571 | 2.524 | 3.081 | 3.459 | 3.798 |
>
> Estimate $\tau$.
>
> **Single point.** $0.632\times4.00 = 2.528\ \mathrm{V}$ is reached at about $t = 0.10\ \mathrm{s}$, so $\tau\approx0.10\ \mathrm{s}$.
>
> **Log-linear fit.** Compute $z = \ln(1 - v/4.00)$: $-0.2055, -0.4988, -0.9970, -1.4708, -2.0006, -2.9858$. The model $z = -t/\tau$ is a line through the origin, so the least-squares slope is
> $$
> \hat m = \frac{\sum t_kz_k}{\sum t_k^2} = -9.947\ \mathrm{s^{-1}}, \qquad \hat\tau = -1/\hat m = 0.1005\ \mathrm{s}.
> $$
> The points lie on the line without systematic curvature (e.g. $z/t = -10.3, -9.98, -9.97, -9.81, -10.0, -9.95$), supporting the first-order structure. Combined with $R = 1\ \mathrm{k\Omega}$ known, $C = \tau/R = 100.5\ \mathrm{\mu F}$.
>
> **Caution.** The log transform amplifies noise near the end of the transient, where $1 - v/v_\infty$ is small; points beyond about $3\tau$ should be down-weighted or omitted, and an error in $v_\infty$ bends the whole line. Nonlinear least squares on the untransformed data avoids both problems.

> [!theorem] Underdamped second-order step response
> For $y'' + 2\zeta\omega_n y' + \omega_n^2 y = \omega_n^2Ku$ with $0<\zeta<1$, the fractional overshoot $M_p$ and the peak time $t_p$ satisfy
> $$
> M_p = e^{-\pi\zeta/\sqrt{1-\zeta^2}}, \qquad t_p = \frac{\pi}{\omega_d} = \frac{\pi}{\omega_n\sqrt{1-\zeta^2}} .
> $$
> Inverting,
> $$
> \zeta = \frac{-\ln M_p}{\sqrt{\pi^2 + \ln^2M_p}}, \qquad \omega_n = \frac{\pi}{t_p\sqrt{1-\zeta^2}} .
> $$
> For a free (or impulse) response with successive peaks $x_k$, $x_{k+1}$ one period apart, the **logarithmic decrement** $\delta = \ln(x_k/x_{k+1}) = 2\pi\zeta/\sqrt{1-\zeta^2}$ gives $\zeta = \delta/\sqrt{4\pi^2+\delta^2}$; using peaks $m$ periods apart, $\delta = \frac1m\ln(x_k/x_{k+m})$, is more accurate.

> [!example] Worked Example 25.26 — Identifying an RLC network from overshoot and peak time
> **Problem.** A filter section modelled as a series RLC circuit with known $C = 1\ \mathrm{\mu F}$ shows, in its capacitor-voltage step response, 25 % overshoot at $t_p = 1.2\ \mathrm{ms}$. Identify $\zeta$, $\omega_n$, $L$ and $R$.
>
> **Damping.** $\ln 0.25 = -1.3863$:
> $$
> \zeta = \frac{1.3863}{\sqrt{9.8696 + 1.9218}} = \frac{1.3863}{3.4339} = 0.404 .
> $$
> **Natural frequency.** $\omega_d = \pi/t_p = 2618\ \mathrm{rad/s}$; $\sqrt{1-\zeta^2} = 0.9149$, so $\omega_n = 2862\ \mathrm{rad/s}$ ($f_n = 455\ \mathrm{Hz}$).
>
> **Element values.** From $\omega_n = 1/\sqrt{LC}$: $L = 1/(\omega_n^2C) = 1/(8.19\times10^6\times10^{-6}) = 0.122\ \mathrm{H}$. From $\zeta = \tfrac R2\sqrt{C/L}$: $R = 2\zeta\sqrt{L/C} = 2(0.404)(349.5) = 282\ \Omega$.
>
> **Identifiability note.** Without knowing $C$, the step response determines only $\omega_n$ and $\zeta$ — two numbers — and not the three element values. We return to this below.

### Linear least squares

If the model output is linear in the parameters, $\hat y_k = \boldsymbol{\varphi}_k^\mathsf{T}\mathbf{p}$ (the **regressors** $\boldsymbol{\varphi}_k$ may be nonlinear functions of the inputs), stacking gives $\hat{\mathbf{y}} = \mathbf{X}\mathbf{p}$ with $\mathbf{X}$ the $N\times p$ regression matrix.

> [!theorem] Least-squares estimate and its covariance
> The minimizer of $\lVert\mathbf{y}-\mathbf{X}\mathbf{p}\rVert^2$ satisfies the **normal equations** $\mathbf{X}^\mathsf{T}\mathbf{X}\hat{\mathbf{p}} = \mathbf{X}^\mathsf{T}\mathbf{y}$; if $\mathbf{X}$ has full column rank,
> $$
> \hat{\mathbf{p}} = (\mathbf{X}^\mathsf{T}\mathbf{X})^{-1}\mathbf{X}^\mathsf{T}\mathbf{y}.
> $$
> If $\mathbf{y} = \mathbf{X}\mathbf{p}_0 + \boldsymbol{\epsilon}$ with $E[\boldsymbol{\epsilon}] = \mathbf{0}$ and $\operatorname{Cov}(\boldsymbol{\epsilon}) = \sigma^2\mathbf{I}$, then $\hat{\mathbf{p}}$ is unbiased and
> $$
> \operatorname{Cov}(\hat{\mathbf{p}}) = \sigma^2(\mathbf{X}^\mathsf{T}\mathbf{X})^{-1},
> $$
> with $\sigma^2$ estimated by $s^2 = \mathrm{SSE}/(N-p)$. (Gauss–Markov: no other linear unbiased estimator has smaller variance.)

*Derivation.* $J = (\mathbf{y}-\mathbf{X}\mathbf{p})^\mathsf{T}(\mathbf{y}-\mathbf{X}\mathbf{p})$ has gradient $-2\mathbf{X}^\mathsf{T}(\mathbf{y}-\mathbf{X}\mathbf{p})$ (Ch. 6); setting it to zero gives the normal equations, and the Hessian $2\mathbf{X}^\mathsf{T}\mathbf{X}$ is positive definite under full rank, so the stationary point is the global minimum (Ch. 23). Substituting $\mathbf{y}$: $\hat{\mathbf{p}} - \mathbf{p}_0 = (\mathbf{X}^\mathsf{T}\mathbf{X})^{-1}\mathbf{X}^\mathsf{T}\boldsymbol{\epsilon}$, whose mean is zero and whose covariance is $(\mathbf{X}^\mathsf{T}\mathbf{X})^{-1}\mathbf{X}^\mathsf{T}(\sigma^2\mathbf{I})\mathbf{X}(\mathbf{X}^\mathsf{T}\mathbf{X})^{-1} = \sigma^2(\mathbf{X}^\mathsf{T}\mathbf{X})^{-1}$. $\square$

Numerically, the normal equations square the condition number of $\mathbf{X}$; QR factorization or the SVD (Ch. 22, Ch. 24) should be used for anything but small, well-conditioned problems.

Many nonlinear models become linear in transformed parameters. The diode equation, the Arrhenius law, power laws $y = ax^b$ and exponential decays all linearize under a logarithm. The transformation also transforms the noise: additive noise of constant variance on $y$ becomes noise of variance $\sigma^2/y^2$ on $\ln y$, so an unweighted log-domain fit over-weights small values. Conversely, *multiplicative* (percentage) noise becomes additive with constant variance in the log domain — then the log fit is exactly right.

> [!example] Worked Example 25.27 — Extracting diode parameters by linear least squares
> **Problem.** A small-signal diode is measured at 300 K ($V_T = 25.85\ \mathrm{mV}$); the ammeter error is about 2 % of reading.
>
> | $V$ (V) | 0.40 | 0.45 | 0.50 | 0.55 | 0.60 | 0.65 |
> |---|---|---|---|---|---|---|
> | $I$ (A) | $1.744\times10^{-5}$ | $5.186\times10^{-5}$ | $1.653\times10^{-4}$ | $4.737\times10^{-4}$ | $1.441\times10^{-3}$ | $4.479\times10^{-3}$ |
>
> Estimate $I_S$ and $n$ with standard errors.
>
> **Linearize the model.** For $I\gg I_S$, $\ln I = \ln I_S + V/(nV_T)$: linear in $p_1 = \ln I_S$ and $p_2 = 1/(nV_T)$, with regressors $\boldsymbol\varphi_k = [1,\ V_k]^\mathsf{T}$. Because the noise is a percentage of reading, it is additive with constant variance on $\ln I$, so an unweighted fit is appropriate.
>
> **Normal equations.** $\sum1 = 6$, $\sum V = 3.15$, $\sum V^2 = 1.6975$, giving $\mathbf{X}^\mathsf{T}\mathbf{X} = \begin{bmatrix}6 & 3.15\\3.15 & 1.6975\end{bmatrix}$. Solving with $\mathbf{X}^\mathsf{T}\mathbf{y}$ from $y_k = \ln I_k$:
> $$
> \hat p_1 = -19.820, \qquad \hat p_2 = 22.153\ \mathrm{V^{-1}} .
> $$
> **Parameters.** $\hat n = 1/(\hat p_2V_T) = 1/(22.153\times0.025\,85) = 1.746$; $\hat I_S = e^{-19.820} = 2.47\times10^{-9}\ \mathrm{A}$.
>
> **Uncertainty.** Residual standard error $s = 0.0237$ (i.e. 2.4 % — consistent with the stated 2 % meter error). $\operatorname{Cov}(\hat{\mathbf{p}}) = s^2(\mathbf{X}^\mathsf{T}\mathbf{X})^{-1}$ gives standard errors $0.060$ for $\hat p_1$ and $0.113\ \mathrm{V^{-1}}$ for $\hat p_2$. Propagating (Ch. 20): $\sigma_n = \sigma_{p_2}/(p_2^2V_T) = 0.009$, and since $\sigma_{\ln I_S} = 0.060$, $I_S$ is known to about 6 %. Result: $n = 1.746\pm0.009$, $I_S = (2.47\pm0.15)\ \mathrm{nA}$.
>
> **Correlation.** The correlation coefficient between $\hat p_1$ and $\hat p_2$ is $-0.987$. The intercept $\ln I_S$ is an extrapolation to $V = 0$, far outside the data (0.40–0.65 V); any tilt of the line moves the intercept. This strong correlation is a mild form of the identifiability problem discussed next. Re-parameterizing around the data centre, $\ln I = \ln I(V_0) + (V-V_0)/(nV_T)$ with $V_0 = 0.525\ \mathrm{V}$, makes the two estimates uncorrelated.
>
> **Validation.** At higher currents the series resistance $R_s$ makes the curve bend below the straight line; residuals trending negative at the top of the range would be the signal to add $R_s$ to the model ($V = nV_T\ln(I/I_S) + IR_s$, no longer linear in the parameters).

### Nonlinear least squares

When the model is not linear in its parameters — a step response in $\tau$, a resonance in $\omega_n$ and $\zeta$, a diode with series resistance — $J(\mathbf{p})$ is minimized iteratively by the **Gauss–Newton** or **Levenberg–Marquardt** methods of Ch. 23. Linearizing the model about the current estimate $\mathbf{p}^{(m)}$,

$$
\hat{\mathbf{y}}(\mathbf{p}) \approx \hat{\mathbf{y}}(\mathbf{p}^{(m)}) + \mathbf{S}\,\Delta\mathbf{p}, \qquad S_{ki} = \frac{\partial\hat y_k}{\partial p_i}\bigg|_{\mathbf{p}^{(m)}},
$$

turns each step into a linear least-squares problem for $\Delta\mathbf{p}$:

$$
(\mathbf{S}^\mathsf{T}\mathbf{S} + \lambda\mathbf{I})\,\Delta\mathbf{p} = \mathbf{S}^\mathsf{T}\bigl(\mathbf{y} - \hat{\mathbf{y}}(\mathbf{p}^{(m)})\bigr),
$$

with $\lambda = 0$ for Gauss–Newton and $\lambda>0$ adapted for robustness in Levenberg–Marquardt. The **sensitivity matrix** $\mathbf{S}$ plays the role of $\mathbf{X}$; at convergence, $\operatorname{Cov}(\hat{\mathbf{p}})\approx\sigma^2(\mathbf{S}^\mathsf{T}\mathbf{S})^{-1}$. Good starting values — often from the graphical methods above — are essential, because $J$ may have several local minima.

### Identifiability

> [!definition] Structural and practical identifiability
> A parameter is **structurally identifiable** if, with perfect noise-free data from the chosen experiment, it can be determined uniquely (globally) or up to a finite set of values (locally). It is **practically identifiable** if it can be determined with acceptable precision from the finite, noisy data actually available.

Structural identifiability is a property of the *model and experiment*, independent of the data quality. The test for a linear model: compute the transfer function from the measured input to the measured output; its coefficients are all the data can ever reveal. Each parameter is identifiable only if it can be recovered from those coefficients.

> [!example] Worked Example 25.28 — Which RC parameters can a step test identify?
> **Problem.** In the circuit of Example 25.1, a voltage step drives $R$ and $C$, and only the capacitor voltage is measured. Can $R$ and $C$ both be identified? What if the source current is measured too?
>
> **Voltage only.** The transfer function is $V_C(s)/V_s(s) = 1/(RCs + 1)$. Its single coefficient is $RC$. Every pair $(R, C)$ with the same product gives *identical* data: $R$ and $C$ are individually **not structurally identifiable**; only $\tau = RC$ is. No amount of data or cleverness changes this.
>
> **Current as well.** $I(s)/V_s(s) = Cs/(RCs+1)$. The initial current after the step is $V_s/R$, giving $R$; combined with $\tau$, $C = \tau/R$. Both are now identifiable. The remedy for structural non-identifiability is a *different experiment* — another sensor, another input — not a better optimizer.
>
> The same reasoning explains Example 25.26: the capacitor-voltage step response of a series RLC circuit has transfer function $1/(LCs^2 + RCs + 1)$, i.e. two coefficients for three unknowns, so one element value must be known independently.

Practical identifiability is diagnosed from the **Fisher information matrix** $\mathbf{F} = \mathbf{S}^\mathsf{T}\mathbf{S}/\sigma^2$, whose inverse approximates the covariance of the estimates (the Cramér–Rao bound). Warning signs are a large condition number of $\mathbf{F}$, correlation coefficients between parameter estimates near $\pm1$, and standard errors comparable to the parameters themselves.

> [!example] Worked Example 25.29 — Practical identifiability depends on the experiment duration
> **Problem.** A first-order step response $y = K(1 - e^{-t/\tau})$ with true $K = 1$, $\tau = 1\ \mathrm{s}$ is sampled at 50 equally spaced points over a window $[0, W]$, with noise $\sigma = 0.01$. Compute the standard errors of $\hat K$ and $\hat\tau$ and their correlation for $W = 0.2\ \mathrm{s}$, $1\ \mathrm{s}$ and $5\ \mathrm{s}$.
>
> **Sensitivities.** $\partial y/\partial K = 1 - e^{-t/\tau}$, $\partial y/\partial\tau = -(Kt/\tau^2)e^{-t/\tau}$. Form $\mathbf{S}$ (50 × 2), $\mathbf{F} = \mathbf{S}^\mathsf{T}\mathbf{S}/\sigma^2$, and invert:
>
> | Window $W$ | $\sigma_{\hat K}$ | $\sigma_{\hat\tau}$ (s) | corr$(\hat K,\hat\tau)$ | cond$(\mathbf{F})$ |
> |---|---|---|---|---|
> | $0.2\tau$ | 0.63 | 0.68 | 0.9998 | 9 600 |
> | $1\tau$ | 0.028 | 0.041 | 0.993 | 320 |
> | $5\tau$ | 0.0026 | 0.0097 | 0.75 | 36 |
>
> **Interpretation.** For $t\ll\tau$, $y\approx Kt/\tau$: the early response depends only on the *ratio* $K/\tau$, so the two parameters are almost perfectly correlated and individually undetermined (63 % and 68 % uncertainty despite 1 % noise!). Recording until the response has settled separates them, since the final value fixes $K$ alone. The model is structurally identifiable in all three cases; practically, only the long experiment identifies it.
>
> **Design rule.** Plan identification experiments so that each parameter has a distinct, strong effect on the measured output somewhere in the record: run step tests to at least $4$–$5\tau$, excite all frequencies of interest (Ch. 21: a persistently exciting input, e.g. a pseudo-random binary sequence), and measure the variables that break structural ambiguities.

### Discrete-time (ARX) identification

For sampled data, a powerful shortcut is to fit a difference equation directly. The first-order ARX model

$$
y[n] = a\,y[n-1] + b\,u[n-1] + \epsilon[n]
$$

is *linear* in $(a, b)$, so linear least squares applies with regressors $\boldsymbol\varphi_n = [y[n-1],\ u[n-1]]^\mathsf{T}$. The continuous parameters follow from the ZOH relations $a = e^{-T/\tau}$, $b = K(1-a)$: $\hat\tau = -T/\ln\hat a$, $\hat K = \hat b/(1-\hat a)$. Because past noisy outputs appear as regressors, the estimate is biased when the measurement noise is not white in the equation-error sense; instrumental-variable or output-error methods remove the bias. Higher-order ARX, ARMAX and subspace methods generalize the idea and are the workhorses of black- and grey-box system identification.

> [!warning] Common Pitfall — fitting parameters that the model cannot carry
> If the residuals of Section 25.7 show structure, fitted parameters are biased in ways that cannot be fixed by more data (Example 25.23: the first-order fit gave $\tau = 1.33\ \mathrm{s}$, matching neither true time constant). Parameters are only meaningful once the structure has been validated; identification and validation must iterate together.

## 25.9 Analytical versus numerical solutions

### Two ways to answer the same question

Once a model is formulated, scaled, simplified and calibrated, it must be solved. There are two broad routes. An **analytical** (closed-form) solution expresses the answer as a formula in the parameters — $v(t) = v_\infty(1-e^{-t/RC})$, $\delta = \sqrt{2/(\omega\mu\sigma)}$, $r_d = nV_T/I_D$. A **numerical** solution computes the answer for specific parameter values by an algorithm — an ODE integrator, a finite-element solver, Newton iteration, a Monte Carlo simulation (Ch. 22, Ch. 24). Neither is universally better. The choice is a modeling decision, made with the purpose in mind.

| Criterion | Analytical | Numerical |
|---|---|---|
| Scope | linear, low-order, simple geometry, special nonlinearities | essentially unrestricted: nonlinear, high-order, complex geometry |
| Parameter dependence | explicit — sensitivities, scaling laws and design rules read off directly | one parameter set per run; trends need sweeps |
| Insight | high: shows *why* and which groups matter | low unless combined with scaling and visualization |
| Accuracy | exact for the model (but the model may be oversimplified to permit it) | approximate: discretization, round-off, convergence tolerances |
| Cost | cheap to evaluate once found; finding it may be impossible | cheap to set up; can be expensive to run (stiff, 3-D, many runs) |
| Failure modes | algebra errors; hidden assumptions made "to get a formula" | instability, non-convergence, under-resolved meshes, silently wrong answers |
| Verification | check limits, dimensions, special cases | compare with analytical cases, refine steps, check conservation |

The most effective practice combines both: use analysis to understand structure, find scales, derive limiting cases and check codes; use computation to handle the full model.

### When a closed form exists — and when it does not

Closed forms are available for linear constant-coefficient ODEs and PDEs on simple domains (Chs. 11–16), for some nonlinear first-order equations (separable, Bernoulli, Riccati with a known solution), and for special nonlinear problems whose solutions are **special functions** (elliptic integrals for the pendulum, Bessel functions for cylindrical problems, the Lambert $W$ function for diode circuits). Outside these islands — nonlinear circuits of order two or more, coupled multi-physics, irregular geometry — numerical methods are the only route.

Between the two lie **semi-analytical** methods: perturbation expansions in a small parameter (Section 25.4), averaging, describing functions (Section 25.5), and series or transform solutions evaluated numerically. They often give the best of both worlds — formulas that are accurate in a stated regime.

> [!example] Worked Example 25.30 — A diode circuit: special function versus root-finding
> **Problem.** A source $V_s = 5\ \mathrm{V}$ drives the diode of Example 25.27 ($I_S = 2.47\ \mathrm{nA}$, $n = 1.746$) through $R = 1\ \mathrm{k\Omega}$. Find the operating point (a) numerically, (b) in closed form.
>
> **Equation.** KCL: $\dfrac{V_s - v}{R} = I_S\left(e^{v/(nV_T)}-1\right)$, a transcendental equation with no solution in elementary functions.
>
> **(a) Numerical.** Bracketing on $[0, 1]\ \mathrm{V}$ and applying Brent's method (or Newton's method, Ch. 22) gives $v = 0.6492\ \mathrm{V}$, $I = 4.351\ \mathrm{mA}$ in a few iterations.
>
> **(b) Closed form.** Rearranging to the form $we^w = z$ gives, with $a = V_s + RI_S$,
> $$
> v = a - nV_T\,W\!\left(\frac{RI_S}{nV_T}\,e^{a/(nV_T)}\right),
> $$
> where $W$ is the Lambert $W$ function (the inverse of $we^w$). Evaluating gives $v = 0.6492\ \mathrm{V}$, identical.
>
> **Which is better?** For one operating point, either. The closed form, however, shows the structure at a glance — $v$ depends on $R I_S/(nV_T)$ and $V_s/(nV_T)$ only, two dimensionless groups — and can be differentiated to obtain sensitivities in closed form. In practice $W$ is itself evaluated by a Newton iteration inside a library, so the "analytical" answer is numerical one level down. The real distinction is not formula versus algorithm, but whether the parameter dependence is *explicit*.

### Numerical solutions must be verified

A numerical solution is an approximation whose error must be controlled. Analytical solutions of simplified cases are the main verification tool, and analytical insight (eigenvalues, time scales) predicts the step sizes a numerical method needs.

> [!theorem] Stability of forward Euler for a linear system
> Applied to $\dot{\mathbf{x}} = \mathbf{A}\mathbf{x}$ with step $h$, forward Euler gives $\mathbf{x}_{k+1} = (\mathbf{I} + h\mathbf{A})\mathbf{x}_k$, which is stable if and only if $\lvert1 + h\lambda_i\rvert < 1$ for every eigenvalue $\lambda_i$ of $\mathbf{A}$ (for diagonalizable $\mathbf{A}$; Ch. 22). For $\lambda = -\sigma \pm j\omega_d$, this requires
> $$
> h < \frac{2\sigma}{\sigma^2+\omega_d^2} = \frac{2\zeta}{\omega_n} .
> $$

> [!example] Worked Example 25.31 — Verifying an integrator on the RLC circuit
> **Problem.** Simulate the step response of the RLC circuit of Example 25.4 ($\omega_n = 31\,623\ \mathrm{rad/s}$, $\zeta = 0.158$) to $t = 1\ \mathrm{ms}$ by forward Euler and by fourth-order Runge–Kutta, and compare with the exact solution $v_C(1\ \mathrm{ms}) = 0.993\,59\ \mathrm{V}$ for a 1 V step.
>
> **Predicted stability limit.** $h_{\max} = 2\zeta/\omega_n = 2(0.158)/31\,623 = 10.0\ \mathrm{\mu s}$ — even though the ringing period is $201\ \mathrm{\mu s}$.
>
> **Results.**
>
> | Method | $h$ (µs) | $v_C(1\ \mathrm{ms})$ (V) | error (V) |
> |---|---|---|---|
> | Forward Euler | 1 | 0.9892 | $4.4\times10^{-3}$ |
> | Forward Euler | 2 | 0.9821 | $1.1\times10^{-2}$ |
> | Forward Euler | 5 | 0.9208 | $7.3\times10^{-2}$ |
> | Forward Euler | 10 | 0.0039 | $\approx1$ (no decay: at the stability limit) |
> | Forward Euler | 12 | −1.698 | diverging |
> | RK4 | 10 | 0.993 59 | $3.4\times10^{-6}$ |
>
> **Interpretation.** The Euler error roughly doubles when $h$ doubles (first-order convergence) until the stability limit predicted analytically is reached, where the answer becomes meaningless. RK4 at $h = 10\ \mathrm{\mu s}$ is a thousand times more accurate than Euler at $1\ \mathrm{\mu s}$. Without the exact solution to compare against, the $h = 5\ \mathrm{\mu s}$ Euler answer (7 % low) might well have been accepted. The habit to cultivate: *before trusting a simulation of the full model, run it on a special case with a known answer*, and *refine the step until the answer stops changing*.

### Stiffness: when analysis tells you which solver to use

> [!definition] Stiff system
> A system is **stiff** when its eigenvalues (time constants) span many orders of magnitude and the behaviour of interest is governed by the slow ones. Explicit methods are then forced, by stability rather than accuracy, to use steps comparable to the *fastest* time constant over a horizon set by the *slowest*.

> [!example] Worked Example 25.32 — A stiff RLC circuit
> **Problem.** A series circuit has $R = 1\ \mathrm{k\Omega}$, a parasitic inductance $L = 1\ \mathrm{\mu H}$, and $C = 1\ \mathrm{\mu F}$. Characterise its time scales and choose a solver for a 5 ms simulation.
>
> **Scaling (Section 25.4).** $\zeta = \tfrac{R}{2}\sqrt{C/L} = 500$, so $\varepsilon = 1/(4\zeta^2) = 10^{-6}$. The eigenvalues are the roots of $\lambda^2 + 10^9\lambda + 10^{12} = 0$: $\lambda_1 \approx -10^3\ \mathrm{s^{-1}}$ ($RC = 1\ \mathrm{ms}$) and $\lambda_2\approx-10^9\ \mathrm{s^{-1}}$ ($L/R = 1\ \mathrm{ns}$). The stiffness ratio is $10^6$.
>
> **Explicit solver.** Forward Euler needs $h < 2/10^9 = 2\ \mathrm{ns}$ for stability, i.e. at least $2.5\times10^6$ steps for 5 ms — all spent resolving a nanosecond transient that ended long ago.
>
> **Implicit solver.** Backward Euler, the trapezoidal rule or BDF methods (Ch. 22) are stable for any step on this problem; accuracy on the slow mode requires only $h\sim10$–$50\ \mathrm{\mu s}$, about 100–500 steps. Circuit simulators use implicit methods for precisely this reason.
>
> **Or remove the stiffness analytically.** Since $\varepsilon\ll1$ and the purpose concerns the millisecond behaviour, the reduced model $RC\dot v + v = v_s$ (Section 25.4) is accurate to about $\varepsilon$ after the first few nanoseconds, and has a closed-form solution. Analysis has not only chosen the solver; it has made the solver unnecessary.

> [!warning] Common Pitfall — believing the simulation because it is detailed
> A detailed numerical model inspires confidence it may not deserve. Mesh too coarse, step too long, tolerance too loose, a unit typed as mV instead of V, a boundary condition on the wrong face — all produce smooth, plausible, wrong results. The analytical checks of this chapter (dimensions, limiting cases, scaling laws, conservation, special-case solutions) are the defence. A simulation that has never been compared with anything is a hypothesis, not a result.

> [!tip] Engineering Note — a practical decision rule
> 1. Always derive what can be derived cheaply: equilibria, scales, dimensionless groups, linearized eigenvalues, limiting cases.
> 2. If the question can be answered from those, stop — the formula is the best deliverable.
> 3. Otherwise, solve numerically, using the analytical results to choose the method and step, to verify the code, and to interpret the output in terms of the dimensionless groups.
> 4. Present numerical results in scaled form wherever possible, so that one set of curves covers a family of designs.

## Synthesis: the book as a modeling toolkit

This chapter has been, in effect, a map of how the rest of the book is used. Every method in the earlier chapters enters the modeling cycle at a particular stage:

| Stage of modeling | Methods from this book |
|---|---|
| Formulating relations | algebra, functions and dimensional consistency (Ch. 1); complex numbers and phasors (Ch. 2–3); derivatives as rates (Ch. 4); integrals as accumulation (Ch. 5); vector calculus for fields and balance laws (Ch. 7) |
| Organizing the equations | matrices and linear systems (Ch. 8); vector spaces, null spaces and rank — the heart of the Pi theorem (Ch. 9); ODEs and PDEs as the language of dynamics and distributed systems (Ch. 11–13) |
| Simplifying | Taylor expansion and Jacobians for linearization (Ch. 4, 6); eigenvalues for modes, stability and time-scale separation (Ch. 10); scaling and dimensional analysis (Ch. 1, 25) |
| Solving | Laplace, Fourier and Z-transforms for linear time-invariant models (Ch. 14–17); complex analysis for inversion and stability (Ch. 18); numerical methods, optimization and computational mathematics when closed forms run out (Ch. 22–24) |
| Confronting data | probability, statistics and random processes for noise, residuals, estimation, identification and uncertainty (Ch. 19–21); least squares and optimization for calibration (Ch. 20, 23) |

Three ideas recur throughout. **Linearity** — superposition, transforms, eigen-decomposition — is the great simplifier, and linearization is how nonlinear reality is brought within its reach. **Dimensionless thinking** — scales, groups, regimes of validity — tells us which effects matter and when a simplification is legitimate. **Comparison with data**, through residuals, likelihoods and cross-validation, is what turns a mathematical exercise into an engineering model. The mathematics of this book is not a collection of separate techniques but a single toolkit, and modeling is the craft of choosing the right tool for each stage of the question.

## Summary

| Topic | Key formula or rule |
|---|---|
| Modeling cycle | purpose → concept → equations → solution → validation → use, iterated |
| Balance law | $\dfrac{d}{dt}(\text{stored}) = \text{in} - \text{out} + \text{generated}$; constitutive relations close the system |
| State-space model | $\dot{\mathbf{x}} = \mathbf{f}(\mathbf{x},\mathbf{u};\mathbf{p})$, $\ \mathbf{y} = \mathbf{h}(\mathbf{x},\mathbf{u};\mathbf{p})$ |
| LTI transfer function | $\mathbf{G}(s) = \mathbf{C}(s\mathbf{I}-\mathbf{A})^{-1}\mathbf{B} + \mathbf{D}$; invariant under $\mathbf{z} = \mathbf{T}\mathbf{x}$ |
| Circuit order | $n = n_C + n_L - n_{\text{C-loops}} - n_{\text{L-cutsets}}$ |
| Buckingham Pi | $n$ quantities, dimension matrix rank $r$ $\Rightarrow$ $n - r$ groups; exponents span $\operatorname{null}(\mathbf{D})$ |
| Skin depth | $\delta = \sqrt{2/(\omega\mu\sigma)}$ (form from dimensions alone) |
| Short dipole | $P = \eta_0 I_0^2\,\Phi(k\ell)$, $\ \Phi = (k\ell)^2/(12\pi)$ for $k\ell\ll1$ |
| Convection | $\mathrm{Nu} = h\ell/k = \Phi(\mathrm{Re},\mathrm{Pr})$, $\ \mathrm{Re} = \rho V\ell/\mu$, $\ \mathrm{Pr} = \mu c_p/k$ |
| Second-order standard form | $y'' + 2\zeta y' + y = u$, $\ \tau = \omega_n t$, $\ \omega_n = 1/\sqrt{LC}$, $\ \zeta = \tfrac R2\sqrt{C/L} = 1/(2Q)$ |
| Step-response metrics | $M_p = e^{-\pi\zeta/\sqrt{1-\zeta^2}}$, $\ t_p = \pi/(\omega_n\sqrt{1-\zeta^2})$, $\ t_s\approx4/(\zeta\omega_n)$ |
| Biot and Fourier numbers | $\mathrm{Bi} = h\ell/k$ (lumped if $<0.1$), $\ \mathrm{Fo} = \alpha t/\ell^2$ |
| Singular perturbation (RLC) | $\varepsilon y'' + y' + y = u$, $\ \varepsilon = L/(R^2C) = 1/(4\zeta^2)$ |
| Van der Pol scaling | $x'' - \mu(1-x^2)x' + x = 0$, $\ \mu = a\sqrt{L/C}$, $\ v = x\sqrt{a/(3b)}$ |
| Linearization | $\mathbf{A} = \partial\mathbf{f}/\partial\mathbf{x}$, $\ \mathbf{B} = \partial\mathbf{f}/\partial\mathbf{u}$ at $(\bar{\mathbf{x}},\bar{\mathbf{u}})$; stability from eigenvalues unless on the $j\omega$-axis |
| Diode small signal | $r_d = nV_T/I_D$; relative error $\approx\lvert\delta v\rvert/(2nV_T)$; $\ \mathrm{HD}_2\approx\hat V/(4nV_T)$ |
| Pendulum | $T/T_0 = \tfrac{2}{\pi}K(\sin\tfrac{\theta_0}{2}) \approx 1 + \theta_0^2/16$ |
| Lumped-circuit validity | $\ell < \lambda/10$, or $t_r \gtrsim 6\,t_{pd}$ |
| Ideal op-amp validity | $\lvert A\beta\rvert\gg1$; closed-loop bandwidth $\approx f_{\text{GBW}}\beta$ |
| Capacitor self-resonance | $f_{\text{SRF}} = 1/(2\pi\sqrt{L_{\text{ESL}}C})$ |
| Goodness of fit | $R^2 = 1 - \mathrm{SSE}/\mathrm{SST}$; $\ s = \sqrt{\mathrm{SSE}/(N-p)}$; $\ \chi^2_\nu = \sum(e_k/\sigma_k)^2/(N-p)\approx1$ |
| Residual whiteness | $\lvert\hat r(k)\rvert < 1.96/\sqrt N$ for most lags |
| LOOCV (linear LS) | $\mathrm{PRESS} = \sum\bigl(e_k/(1-h_{kk})\bigr)^2$, $\ \mathbf{H} = \mathbf{X}(\mathbf{X}^\mathsf{T}\mathbf{X})^{-1}\mathbf{X}^\mathsf{T}$ |
| First-order identification | $K = \Delta y_\infty/\Delta u$; $\tau$ at 63.2 %; $\ln(1 - \Delta y/\Delta y_\infty) = -t/\tau$ |
| Second-order identification | $\zeta = -\ln M_p/\sqrt{\pi^2+\ln^2M_p}$; $\ \omega_n = \pi/(t_p\sqrt{1-\zeta^2})$; $\ \zeta = \delta/\sqrt{4\pi^2+\delta^2}$ |
| Linear least squares | $\hat{\mathbf{p}} = (\mathbf{X}^\mathsf{T}\mathbf{X})^{-1}\mathbf{X}^\mathsf{T}\mathbf{y}$, $\ \operatorname{Cov}(\hat{\mathbf{p}}) = \sigma^2(\mathbf{X}^\mathsf{T}\mathbf{X})^{-1}$ |
| Nonlinear LS step | $(\mathbf{S}^\mathsf{T}\mathbf{S} + \lambda\mathbf{I})\Delta\mathbf{p} = \mathbf{S}^\mathsf{T}(\mathbf{y}-\hat{\mathbf{y}})$ |
| Identifiability | structural: recoverable from transfer-function coefficients; practical: $\mathbf{F} = \mathbf{S}^\mathsf{T}\mathbf{S}/\sigma^2$ well conditioned |
| ARX first order | $y[n] = ay[n-1] + bu[n-1]$; $\ \tau = -T/\ln a$, $\ K = b/(1-a)$ |
| Forward Euler stability | $\lvert1 + h\lambda\rvert < 1$; for $\lambda = -\zeta\omega_n\pm j\omega_d$: $h < 2\zeta/\omega_n$ |
| Stiffness ratio | $\lvert\lambda_{\max}\rvert/\lvert\lambda_{\min}\rvert\gg1$ $\Rightarrow$ implicit solver or reduced model |

## Further reading

- E. A. Bender, *An Introduction to Mathematical Modeling*, Dover, 2000 — short, readable treatment of the modeling cycle, dimensional analysis and model simplification.
- C. C. Lin and L. A. Segel, *Mathematics Applied to Deterministic Problems in the Natural Sciences*, SIAM Classics, 1988 — the standard reference on scaling, nondimensionalization and perturbation methods.
- G. I. Barenblatt, *Scaling*, Cambridge University Press, 2003 — dimensional analysis and similarity in depth, including the Pi theorem and its limitations.
- L. Ljung, *System Identification: Theory for the User*, 2nd ed., Prentice Hall, 1999 — parameter estimation, model validation, identifiability and ARX/output-error methods for dynamic systems.
- K. J. Åström and R. M. Murray, *Feedback Systems: An Introduction for Scientists and Engineers*, 2nd ed., Princeton University Press, 2021 — state-space modeling, linearization and the use of models in control (freely available online).
