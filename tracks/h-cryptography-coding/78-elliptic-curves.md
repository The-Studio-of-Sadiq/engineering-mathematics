# Chapter 78 — Elliptic Curves and Modern Cryptographic Mathematics

[← Previous](77-cryptography.md) · [Track H](README.md) · [Table of Contents](../../00-front/table-of-contents.md) · [Next →](79-error-correcting-codes.md)

> **Level 2 · Track H — Cryptography, Coding, and Digital Security**
> **Builds on:** [Ch. 57](../e-abstract-math/57-number-theory.md), [Ch. 56](../e-abstract-math/56-abstract-algebra.md), [Ch. 77](77-cryptography.md), [Ch. 1](../../part-1-foundations/01-mathematical-language-and-algebra.md) · **Leads to:** [Ch. 79](79-error-correcting-codes.md), [Ch. 82](82-zero-knowledge-proofs.md), [Ch. 94](../j-frontier/94-algebraic-geometry.md) · **Related:** [Ch. 32](../b-probability-information/32-counting-random-worlds.md), [Ch. 26](../a-algorithms/26-the-mathematics-of-algorithms.md), [Ch. 81](81-pseudorandomness.md)
> **Reading time:** ~25 min

## 78.1 The hook

The PlayStation 3 was, for four years, the console nobody could crack. Every program it ran had to carry a digital signature from Sony, made with Sony's private key using the Elliptic Curve Digital Signature Algorithm (ECDSA). The mathematics behind it is the same that protects most of today's secure web connections, and nobody has ever broken it. Without the private key, no unofficial software would run.

In December 2010, a group of hackers calling themselves fail0verflow stood on stage at a computer congress in Berlin and showed that they had the key. They had not solved a hard problem. They had noticed that two different signed programs shared one number in their signatures. ECDSA requires a fresh secret random number for every signature, and Sony had used the same one every time. From two signatures and a few lines of algebra, high-school algebra modulo a prime, the private key fell out. Anyone could now sign software as Sony.

Three years later the same mistake emptied Bitcoin wallets: a flaw in a random-number generator on Android phones made some wallets repeat their signing numbers, and thieves who watched the public blockchain for repeats took the coins.

The cryptosystem was not wrong. The mathematics was flawless, the implementation was not. To understand both the strength and the fragility, we need to know what these signatures are made of. The answer is surprising: a *curve*, a rule for "adding" points on it that is built from straight lines, and arithmetic in a finite field where the curve stops looking like a curve at all.

## 78.2 The big idea

Public-key cryptography needs a group where going forward is easy and going back is hard. In [Ch. 57](../e-abstract-math/57-number-theory.md) the group was $\mathbb{Z}_p^\times$: computing $g^a$ is fast, recovering $a$ from $g^a$ (the **discrete logarithm**) is believed slow. Elliptic curves supply a different group with the same shape and, as far as anyone knows, no shortcuts.

The construction is geometric. Draw the cubic curve $y^2=x^3+ax+b$. Pick two points on it, draw the straight line through them, and the line hits the curve in exactly **one more point**. Reflect that point in the $x$-axis: that is the "sum" of the two. To add a point to itself, use the tangent line. Astonishingly, this rule obeys the laws of a group, so you can add a point $G$ to itself again and again, $G,\,2G,\,3G,\dots$, and the multiples behave like the powers $g,g^2,g^3,\dots$ of a number.

Now do the same arithmetic modulo a prime $p$. The curve dissolves into a scattered cloud of points on a $p\times p$ grid, but the formulas still work, and the multiples $kG$ hop around the cloud with no visible pattern. Computing $kG$ from $k$ is fast, by doubling and adding. Recovering $k$ from $kG$ is the **elliptic-curve discrete logarithm problem**, and the best known general attack is no better than a sophisticated form of guessing.

An analogy: imagine a billiard table with strange, perfectly deterministic walls. You shoot a ball from a known spot, and it bounces $k$ times and stops. Anyone can simulate $k$ bounces quickly if they know $k$. But seeing only where the ball stopped, finding $k$ means re-simulating bounce after bounce, because the table remembers nothing about how many bounces got it there.

## 78.3 The mathematics

### The fields underneath

We need a field: a number system where we can add, subtract, multiply and divide by anything nonzero. Over the real numbers we can draw pictures. For cryptography we use the finite field $\mathbb{F}_p=\mathbb{Z}_p$ of integers mod a prime $p$, where division means multiplying by a modular inverse from the extended Euclidean algorithm. [Ch. 57](../e-abstract-math/57-number-theory.md#573-the-mathematics) builds these fields, and the prime-power fields $\mathbb{F}_{2^m}$, and we use them without re-deriving anything. One caution: the formulas below divide by 2 and by 3, so we assume the field's characteristic is not 2 or 3. (Binary curves exist but use a different equation.)

> [!definition] Elliptic curve
> Let $K$ be a field with $1+1\neq0$ and $1+1+1\neq0$, and let $a,b\in K$ satisfy $4a^3+27b^2\neq0$. The **elliptic curve** $E$ is the set
> $$
> E(K)=\{(x,y)\in K^2:\ y^2=x^3+ax+b\}\ \cup\ \{\mathcal{O}\},
> $$
> where $\mathcal{O}$ is an extra **point at infinity**.

The condition $4a^3+27b^2\neq0$ says the cubic $x^3+ax+b$ has no repeated root, which keeps the curve smooth: no cusps, no self-crossings, and a well-defined tangent everywhere. The point $\mathcal{O}$ is best pictured as sitting infinitely far up (and down) every vertical line; every vertical line meets the curve there.

### Result 1: the chord-and-tangent law

![Two panels showing the curve y squared equals x cubed plus 17 over the real numbers. Left: the orange line through P = (-2, 3) and Q = (2, 5) meets the curve again at R = (1/4, 33/8); a dashed vertical line drops to its mirror image P + Q = (1/4, -33/8). Right: the green tangent at P = (-1, 4) meets the curve again at R, and reflecting gives 2P, approximately (2.14, -5.18).](../../images/ch78-chord-tangent.svg)

*Addition and doubling on $y^2=x^3+17$. All four named points have rational coordinates: rational points add to rational points, because the formulas use only field operations.*

> [!theorem] The addition formulas
> Let $P=(x_1,y_1)$ and $Q=(x_2,y_2)$ lie on $E$ with $Q\neq-P$, where $-P=(x_1,-y_1)$. Put
> $$
> \lambda=\begin{cases}\dfrac{y_2-y_1}{x_2-x_1}, & P\neq Q,\\[2ex] \dfrac{3x_1^2+a}{2y_1}, & P=Q.\end{cases}
> $$
> Then the line of slope $\lambda$ through $P$ meets $E$ in a third point $R=(x_3,-y_3)$, and $P+Q:=(x_3,y_3)$ is given by
> $$
> x_3=\lambda^2-x_1-x_2,\qquad y_3=\lambda(x_1-x_3)-y_1 .
> $$

*Proof.* The line through $P$ with slope $\lambda$ is $y=\lambda(x-x_1)+y_1$. In the case $P\neq Q$ it passes through $Q$ by the choice of $\lambda$. Substitute it into the curve:

$$
\bigl(\lambda(x-x_1)+y_1\bigr)^2=x^3+ax+b
\quad\Longleftrightarrow\quad
x^3-\lambda^2x^2+(\cdots)x+(\cdots)=0 .
$$

This is a monic cubic in $x$, and its roots are the $x$-coordinates where line and curve meet. Two roots are known: $x_1$ and $x_2$. A cubic over a field with two roots in the field has a third, because dividing out $(x-x_1)(x-x_2)$ leaves a linear factor. By Vieta's formula the three roots sum to minus the coefficient of $x^2$:

$$
x_1+x_2+x_3=\lambda^2\quad\Longrightarrow\quad x_3=\lambda^2-x_1-x_2 .
$$

The third intersection point lies on the line, so its $y$-coordinate is $\lambda(x_3-x_1)+y_1$; reflecting changes its sign, giving $y_3=\lambda(x_1-x_3)-y_1$. For $P=Q$, implicit differentiation of $y^2=x^3+ax+b$ gives $2y\,y'=3x^2+a$, so the tangent slope is $(3x_1^2+a)/(2y_1)$; the tangent touches the curve at $P$ with multiplicity two, so $x_1$ counts as a double root and the same Vieta argument applies with $x_2=x_1$. $\blacksquare$

Check it on the figure. For $P=(-2,3)$, $Q=(2,5)$ on $y^2=x^3+17$: $\lambda=\tfrac{2}{4}=\tfrac12$, $x_3=\tfrac14+2-2=\tfrac14$, $y_3=\tfrac12(-2-\tfrac14)-3=-\tfrac{33}{8}$. And indeed $(\tfrac14)^3+17=\tfrac{1089}{64}=(\tfrac{33}{8})^2$. Doubling $P=(-1,4)$: $\lambda=\tfrac38$, and $2P=(\tfrac{137}{64},-\tfrac{2651}{512})$. The denominators grow quickly, which is one reason cryptographers prefer finite fields.

Nothing in the proof used real numbers. It used only $+,-,\times,\div$ and the fact that a cubic with two roots in the field has a third. So **the same formulas define addition on $E(\mathbb{F}_p)$**.

### The group

> [!theorem] $E(K)$ is an abelian group
> With the chord-and-tangent rule, $P+\mathcal{O}=P$, and $P+(-P)=\mathcal{O}$, the set $E(K)$ is an abelian group with identity $\mathcal{O}$.

Most axioms are immediate. The identity: the vertical line through $P$ meets the curve at $P$, $-P$ and $\mathcal{O}$, so "the third point" of $P$ and $\mathcal{O}$ is $-P$, reflected back to $P$. Inverses: the line through $P$ and $-P$ is vertical, its third point is $\mathcal{O}$. Commutativity: the line through $P$ and $Q$ is the line through $Q$ and $P$.

**Associativity is the miracle.** Why should $(P+Q)+R=P+(Q+R)$ for a rule made of lines? One can grind it out from the formulas (a computer algebra system does it in seconds, a human in many pages). A cleaner argument uses a classical fact about cubics: *if two cubic curves meet in exactly nine points, then any cubic passing through eight of them also passes through the ninth* (the Cayley–Bacharach theorem). Multiply the equations of three lines, the line through $P$ and $Q$, the line through $P+Q$ and $R$, and the vertical line through $Q+R$, to get a cubic $C_1$. Likewise the line through $Q$ and $R$, the line through $P$ and $Q+R$, and the vertical line through $P+Q$ give a cubic $C_2$. Each meets $E$ in nine points, and eight are shared: $P,Q,R,\mathcal{O},\pm(P+Q),\pm(Q+R)$. The ninth intersection of $C_1$ with $E$ is $-\bigl((P+Q)+R\bigr)$, and Cayley–Bacharach forces $C_2$ through it, which means it equals $-\bigl(P+(Q+R)\bigr)$. (The special cases where some of the points coincide need a limiting argument.)

> [!mystery] Why lines make a group
> The nine-point argument shows associativity is *true*, but not *why* a cubic curve, of all things, should carry a group law. Conics do not; curves of degree four do not in this simple way. The deep reason is that an elliptic curve is secretly a doughnut: over $\mathbb{C}$ it is the plane modulo a lattice, $\mathbb{C}/\Lambda$, where addition is just addition of complex numbers, and the coordinates $x,y$ are doubly periodic functions of the kind [Ch. 18](../../part-5-transforms/18-complex-analysis.md) hints at. Over any field, the group is a shadow of a more general object, the group of "divisors" on a curve of genus one. That story is [Ch. 94](../j-frontier/94-algebraic-geometry.md).

### Curves over finite fields

![Left: the 18 affine points of y squared = x cubed + 2x + 2 over F_17 on a 17 by 17 grid, connected in the order G, 2G, 3G and so on for G = (5, 1), each labelled with its multiple k; the path jumps erratically across the grid, and the points are symmetric about the dotted line y = 8.5. Right: the same equation over F_263, a cloud of 281 scattered points with the same mirror symmetry.](../../images/ch78-finite-field-points.svg)

*Over a finite field the curve becomes a cloud. The group law still works: the orange path visits all 18 affine points as $G,2G,\dots,18G$ and returns to $\mathcal{O}$ at $19G$. Each $x$ gives either no points or a mirror pair $(x,\pm y)$, the only trace left of the real picture.*

How many points does $E(\mathbb{F}_p)$ have? Each $x\in\mathbb{F}_p$ gives two points if $x^3+ax+b$ is a nonzero square mod $p$, one if it is zero, and none otherwise. Roughly half of the nonzero residues are squares, so we expect about $p$ affine points plus $\mathcal{O}$.

> [!theorem] Hasse's theorem (1933)
> The number of points satisfies $\lvert\,\#E(\mathbb{F}_p)-(p+1)\,\rvert\le2\sqrt p$.

On $y^2=x^3+2x+2$ the count is $19$ over $\mathbb{F}_{17}$, within $1$ of $p+1=18$, and $282$ over $\mathbb{F}_{263}$, within $18$ of $264$; the bound allows $8.2$ and $32.4$. We will not prove Hasse's theorem here; it is the first case of a vast generalisation (the Weil conjectures) that took most of the twentieth century to prove.

Since $19$ is prime, Lagrange's theorem ([Ch. 56](../e-abstract-math/56-abstract-algebra.md#563-the-mathematics)) leaves every point except $\mathcal{O}$ only one possible order, 19: each one generates the whole group, and $E(\mathbb{F}_{17})\cong\mathbb{Z}_{19}$. Cryptographic curves are chosen the same way, with a group of prime order $n$ (or a small multiple of it) near $2^{256}$.

### The discrete logarithm, and how fast it can be attacked

> [!definition] Elliptic-curve discrete logarithm problem (ECDLP)
> Given a point $G$ of order $n$ and a point $Q$ in the group it generates, find $d\in\{0,\dots,n-1\}$ with $Q=dG$.

The forward direction is fast. **Double-and-add** reads the bits of $d$ from the top: double the running total at each bit and add $G$ when the bit is 1. That computes $dG$ in at most $2\log_2 d$ group operations, about 510 for a 256-bit $d$. Brute force in reverse tries $G,2G,3G,\dots$ and needs up to $n$ steps. Can we do better than brute force without any special structure? Yes, but only to the square root.

> [!theorem] Baby-step giant-step (Shanks)
> The ECDLP in a group of order $n$ can be solved with $O(\sqrt n)$ group operations and $O(\sqrt n)$ memory.

*Proof.* Let $m=\lceil\sqrt n\,\rceil$. Every $d$ with $0\le d<n\le m^2$ can be written uniquely as $d=im+j$ with $0\le i,j<m$ (division with remainder). Then

$$
Q=dG=imG+jG\quad\Longleftrightarrow\quad Q-i(mG)=jG .
$$

**Baby steps:** compute $jG$ for $j=0,\dots,m-1$ and store them in a hash table ($m$ additions). **Giant steps:** compute $Q,\ Q-mG,\ Q-2mG,\dots$ and look each one up ($\le m$ additions). Some $i<m$ must produce a hit, because the true $d$ has such a decomposition; the hit gives $d=im+j$. Total work is at most $2m=O(\sqrt n)$. $\blacksquare$

For a group of order near $2^{256}$, $\sqrt n\approx2^{128}\approx3.4\times10^{38}$ operations: the same scale as brute-forcing a 128-bit AES key in [Ch. 77](77-cryptography.md#773-the-mathematics). The memory problem is removed by **Pollard's rho** method, which takes a pseudorandom walk through the group and waits for it to collide with itself. By the birthday bound of [Ch. 32](../b-probability-information/32-counting-random-worlds.md#the-birthday-bound-made-general), a collision arrives after about $\sqrt{\pi n/2}$ steps, using almost no memory.

Two further facts shape how curves are chosen.

- **Pohlig–Hellman.** If the group order $n$ factors into small primes, the ECDLP splits by the Chinese remainder theorem into small problems, one per prime power, so its cost is governed by the *largest prime factor* of $n$, not $n$ itself. Hence the insistence on prime (or nearly prime) group order.
- **Generic lower bound.** Shoup proved in 1997 that any algorithm that treats the group as a black box, only adding and comparing points, needs $\Omega(\sqrt n)$ operations. Baby-step giant-step and rho are optimal among such algorithms.

**The state of play.** For $\mathbb{Z}_p^\times$ there are faster, non-generic attacks: **index calculus** factors many powers $g^k$ over a base of small primes and solves linear equations for their logarithms, giving sub-exponential running time. That is why finite-field Diffie–Hellman and RSA need 3072-bit numbers for 128-bit security. Points on an elliptic curve have no known notion of "small" that a factor base could exploit. So, as far as anyone knows, the best attack on a well-chosen curve is generic, and 256-bit curves give roughly 128-bit security. That is the whole reason elliptic-curve keys are twelve times shorter.

"Well-chosen" carries weight. Some curves *do* have shortcuts: on **supersingular** curves a pairing maps the ECDLP into a finite field where index calculus applies (the MOV attack, 1993), and on **anomalous** curves with exactly $p$ points the ECDLP can be solved in polynomial time (Smart; Satoh–Araki; Semaev, 1998–99). Standard curves are checked against every known weakness, and whether some cleverer index calculus exists for general curves is an open question.

### Pairings: a weakness that became a tool

A **pairing** is a map $e:\ \langle P\rangle\times\langle Q\rangle\to\mathbb{F}_{p^k}^\times$ from pairs of points of prime order $n$ into the multiplicative group of an extension field, which is **bilinear**, $e(aP,bQ)=e(P,Q)^{ab}$, and **non-degenerate**, $e(P,Q)\neq1$. The Weil and Tate pairings, built from rational functions on the curve, have these properties for every elliptic curve, and Miller's algorithm (1986) evaluates them with $O(\log n)$ operations, in a sense the double-and-add of pairings. The catch is the field they land in. The values are $n$-th roots of unity, which live in $\mathbb{F}_{p^k}$ only once $n$ divides $p^k-1$; the least such $k$ is the **embedding degree**. For a random curve $k$ is about as large as $n$ itself, so the target field has some $2^{2^{256}}$ elements and the pairing cannot even be written down. Only special curves have small $k$: supersingular curves have $k\le6$, and families such as BN and BLS12 are constructed on purpose to have $k=12$.

Small $k$ cuts both ways. It is exactly the MOV attack: $e(P,dP)=e(P,P)^d$ moves the discrete logarithm into $\mathbb{F}_{p^k}^\times$, where index calculus works, so pairing-friendly curves must be made larger to stay safe. But it also gives a capability no ordinary group has: anyone can test a *multiplicative* relation between hidden exponents. Given $aG$, $bG$ and $cG$, the check $e(aG,bG)=e(G,cG)$ holds exactly when $ab\equiv c\pmod n$, yet reveals none of the three. **BLS signatures** use this directly: the signature on $m$ is $\sigma=d\,H(m)$ for a hash-to-curve point $H(m)$, and it verifies when $e(\sigma,G)=e(H(m),dG)$. Signatures from many signers on many messages can be added into one point, which is why blockchains with thousands of validators use them.

### Using the group: key exchange and signatures

**ECDH** is Diffie–Hellman with points: public $G$; Alice sends $aG$, Bob sends $bG$, and both compute $abG$. **ECDSA** signs a message hash $h$ (an integer mod $n$) with private key $d$ and public key $Q=dG$:

1. pick a secret, fresh, uniformly random **nonce** $k\in\{1,\dots,n-1\}$;
2. compute $r=x(kG)\bmod n$ (the $x$-coordinate of $kG$, reduced mod $n$);
3. compute $s=k^{-1}(h+rd)\bmod n$. The signature is $(r,s)$.

To verify, compute $w=s^{-1}$, $u_1=hw$, $u_2=rw$ and check that $x(u_1G+u_2Q)\bmod n=r$. This works because $u_1G+u_2Q=(h+rd)wG=kG$.

### Result 3: one repeated nonce reveals the key

> [!theorem] Nonce reuse breaks ECDSA
> If two signatures $(r,s_1)$ on $h_1$ and $(r,s_2)$ on $h_2\neq h_1$ were made with the same nonce $k$, then anyone can compute
> $$
> k=\frac{h_1-h_2}{s_1-s_2}\bmod n,\qquad d=\frac{s_1k-h_1}{r}\bmod n .
> $$

*Proof.* The same $k$ gives the same $kG$, hence the same $r$; that is how the repeat is spotted. From $s_ik\equiv h_i+rd\pmod n$ for $i=1,2$, subtract: $(s_1-s_2)k\equiv h_1-h_2$. Since $h_1\neq h_2$ and $k$ is invertible, $s_1\neq s_2$, so divide to get $k$. Then solve $s_1k=h_1+rd$ for $d$ ($r\neq0$ is required by the standard). Everything happens in the field $\mathbb{Z}_n$, since $n$ is prime. $\blacksquare$

That is the PlayStation 3 key and the Android wallets. The theorem needs no hard problem to be solved; it is two linear equations in two unknowns. Worse, it is not only exact repetition that kills: if the nonces are merely *biased*, say their top few bits are often zero, lattice methods recover the key from a few hundred signatures. This is why modern implementations derive $k$ deterministically from the key and the message (RFC 6979) or use schemes like Ed25519 that do so by design.

> [!mystery] From one multiplication to a whole computation
> A pairing checks a single product of hidden exponents, $ab=c$. Yet pairing-based SNARKs convince a verifier, with three curve points, that a computation of a billion steps was carried out correctly, while revealing nothing about its inputs. How does one hidden multiplication grow into a check of an entire program, and what does it mean for a verifier to "learn nothing"? See [Ch. 82](82-zero-knowledge-proofs.md).

## 78.4 The algorithm

Here is the whole pipeline, from point addition to the PlayStation-style key recovery, on the 19-point curve of the figure. Every value is small enough to check by hand against the labelled multiples in the left panel.

```python
# Toy ECDSA on y^2 = x^3 + 2x + 2 over F_17 (group order 19). NEVER use for real data.
p, a, b = 17, 2, 2
G, n = (5, 1), 19

def add(P, Q):
    if P is None: return Q
    if Q is None: return P
    (x1, y1), (x2, y2) = P, Q
    if x1 == x2 and (y1 + y2) % p == 0: return None            # P + (-P) = O
    if P == Q: lam = (3 * x1 * x1 + a) * pow(2 * y1, -1, p) % p
    else:      lam = (y2 - y1) * pow(x2 - x1, -1, p) % p
    x3 = (lam * lam - x1 - x2) % p
    return (x3, (lam * (x1 - x3) - y1) % p)

def mul(k, P):                                  # double-and-add, NOT constant time
    R = None
    for bit in bin(k)[2:]:
        R = add(R, R)
        if bit == "1": R = add(R, P)
    return R

def sign(d, h, k):                              # k MUST be fresh and secret for every signature
    r = mul(k, G)[0] % n
    s = pow(k, -1, n) * (h + r * d) % n
    return r, s

def verify(Q, h, r, s):
    w = pow(s, -1, n)
    X = add(mul(h * w % n, G), mul(r * w % n, Q))
    return X is not None and X[0] % n == r

d = 11; Q = mul(d, G)                            # private key, public key
r1, s1 = sign(d, 5, k=10)                        # two messages, SAME nonce
r2, s2 = sign(d, 12, k=10)
print(Q, (r1, s1), (r2, s2), verify(Q, 5, r1, s1), verify(Q, 12, r2, s2))

# Eve sees only public data: r1 == r2 betrays the repeated nonce
k = (5 - 12) * pow(s1 - s2, -1, n) % n
d_found = (s1 * k - 5) * pow(r1, -1, n) % n
print("recovered k =", k, " recovered d =", d_found)
```

The output is `(13, 10) (7, 12) (7, 7) True True`, then `recovered k = 10  recovered d = 11`. Both signatures verify, both share $r=7$, and Result 3 hands Eve the private key.

**Complexity.** One addition costs one modular inversion and a few multiplications. Real implementations avoid the inversion by working in **projective coordinates** $(X:Y:Z)$, with $x=X/Z$, $y=Y/Z$, which turn each addition into about a dozen multiplications and postpone the single division to the end. Scalar multiplication is $O(\log n)$ additions, so $O(\lambda^3)$ bit operations for a $\lambda$-bit curve, on par with RSA at a far smaller $\lambda$.

**Failure cases, mathematical.** If $4a^3+27b^2=0$, the curve is singular, the group law degenerates, and its "discrete logarithm" collapses into one in the additive or multiplicative group of the field, which is easy or much easier. If $n$ is not prime, Pohlig–Hellman splits the problem into pieces. If a received point is not checked to lie on the curve, an attacker can send a point on a *different* curve with a weak group (the formulas never use $b$!) and harvest the private key bit by bit. This **invalid-curve attack** has broken real TLS implementations.

**Failure cases, engineering.** `mul` branches on the bits of the secret scalar, so its running time and power draw leak them. Production code uses constant-time ladders (the Montgomery ladder does one addition and one doubling for every bit, whatever its value), complete addition formulas without special cases, and libraries such as libsodium or a platform's TLS stack. Nonces must come from a cryptographically secure generator, or be derived deterministically as in RFC 6979 or Ed25519.

> [!warning] Never roll your own crypto
> The 30 lines above are correct mathematics and a catastrophically insecure implementation: a 19-element group, a branching scalar multiplication, no point validation, a caller-supplied nonce and no hashing of the message. Production elliptic-curve code is written by specialists, verified, fuzzed, audited and attacked for years. Use a vetted library's high-level API (for example X25519 for key exchange and Ed25519 for signatures), standard curves, and current guidance; never invent curves, nonce schemes or "simplified" protocols.

> [!mystery] The generator with a hidden key
> Where do the "random" numbers of cryptography come from? In 2006 a standards body published Dual_EC_DRBG, a random-number generator built from two elliptic-curve points $P$ and $Q$: each output is roughly the $x$-coordinate of a multiple of $Q$, and the internal state advances via $P$. Researchers soon pointed out that whoever knew a secret $e$ with $P=eQ$ could predict every future output from one sample. The points' origin was never explained, and documents leaked in 2013 suggested that the generator had been deliberately backdoored. How can a generator look perfectly random to everyone but one person, and what definition of "random" would have ruled this out? See [Ch. 81](81-pseudorandomness.md).

## 78.5 The experiment

First, drag the geometry. The block draws $y^2=x^3+17$ with the chord through $P=(-2,3)$ and $Q=(2,5)$ and its reflection.

```desmos-graph
left=-4; right=6; bottom=-10; top=10
---
y^2=x^3+17|blue
y=0.5x+4|orange
x=0.25|grey|dashed
(-2,3)|label:P
(2,5)|label:Q
(0.25,4.125)|label:R
(0.25,-4.125)|red|label:P+Q
```

*The chord through $P$ and $Q$ on $y^2=x^3+17$ meets the curve at $R=(\tfrac14,\tfrac{33}{8})$; its mirror image is $P+Q$.*

- Replace the orange line by the tangent at $(-1,4)$, `y=0.375x+4.375`, and find where it meets the curve again (near $x=2.14$).
- Change the curve to `y^2=x^3-3x+2`. Here $4a^3+27b^2=-108+108=0$: the curve crosses itself at $(1,0)$, where no tangent is defined. Then try `y^2=x^3-3x+3` and `y^2=x^3-3x+1`; the second splits into two pieces, an oval and an infinite branch, and the group law still works across both.

Then go finite. The script counts the points of a curve mod $p$, checks Hasse's bound, and attacks a discrete logarithm with baby-step giant-step.

```python
import numpy as np
import matplotlib.pyplot as plt

P_MOD, A, B = 1009, 2, 3      # KNOBS: a prime p > 3 and the curve y^2 = x^3 + A x + B
SECRET = 777                  # KNOB: Alice's private scalar d

sq = {}                                            # y^2 mod p -> list of y
for y in range(P_MOD): sq.setdefault(y * y % P_MOD, []).append(y)
pts = [(x, y) for x in range(P_MOD) for y in sq.get((x**3 + A*x + B) % P_MOD, [])]
N = len(pts) + 1                                   # + point at infinity (None)
print(f"N = {N}, p + 1 = {P_MOD + 1}, Hasse bound 2*sqrt(p) = {2*np.sqrt(P_MOD):.1f}")

def add(P, Q):
    if P is None: return Q
    if Q is None: return P
    (x1, y1), (x2, y2) = P, Q
    if x1 == x2 and (y1 + y2) % P_MOD == 0: return None
    lam = ((3*x1*x1 + A) * pow(2*y1, -1, P_MOD) if P == Q else (y2 - y1) * pow(x2 - x1, -1, P_MOD)) % P_MOD
    x3 = (lam * lam - x1 - x2) % P_MOD
    return (x3, (lam * (x1 - x3) - y1) % P_MOD)

def mul(k, P):                                     # double-and-add
    R = None
    for bit in bin(k)[2:]:
        R = add(R, R)
        if bit == "1": R = add(R, P)
    return R

G = pts[1]
n = next(k for k in range(1, N + 1) if N % k == 0 and mul(k, G) is None)   # order of G
Q = mul(SECRET % n, G)                             # Alice's public key
m = int(np.ceil(np.sqrt(n)))                       # baby-step giant-step
baby = {}; R = None
for j in range(m): baby.setdefault(R, j); R = add(R, G)
step, R = mul(n - m, G), Q                         # step = -mG
for i in range(m):
    if R in baby: print(f"order of G = {n}; recovered d = {i*m + baby[R]} after {m + i} group ops"); break
    R = add(R, step)
plt.scatter(*np.array(pts).T, s=2); plt.title(f"{N} points mod {P_MOD}"); plt.show()
```

With the defaults the curve has $N=1068$ points, within 58 of $p+1=1010$ and inside the Hasse bound of 63.5. The chosen $G$ has order only $267$, so the script "recovers" $d=243$: the secret $777$ only matters modulo the order, and $777-2\cdot267=243$. Things to try:

- **Factor $N$.** $1068=2^2\cdot3\cdot89$. A real attacker would not run baby-step giant-step on the full group but on the pieces of size 4, 3 and 89, as Pohlig–Hellman suggests. Can you find $A,B$ for $p=1009$ that make $N$ prime?
- **Scale $p$.** Try $p=10007$ and $p=100003$, keeping a prime-order or near-prime-order curve. The number of group operations grows like $\sqrt N$: it is about 320 for $p=100003$. Extrapolate to $N\approx2^{256}$.
- **Break the curve.** Set $A=0$, $B=0$, which makes $4A^3+27B^2=0$. The count comes out as exactly $p+1$, and the script then crashes while looking for the order of $G$. Which point is to blame, and what group do the remaining points form?

## 78.6 The challenge

No answers are given.

1. **Intuition.** On the real curve, a line through two points of an elliptic curve meets it in a third. Why must a line through two points of a *conic* (a circle, say) fail to give a useful "third point", and what replaces the group law there? Then explain, using only the pictures, why the points of order 2 on $y^2=x^3+ax+b$ are exactly those with $y=0$, and why there are at most three of them.
2. **Derivation.** Show that a point $P=(x,y)$ has order 3 exactly when it is an inflection point of the curve. Derive the polynomial whose roots are the $x$-coordinates of the points of order 3, and check that it has degree 4. Then prove that the map $x\mapsto x^3+ax+b$ cannot be a square for *every* $x\in\mathbb{F}_p$, so that $\#E(\mathbb{F}_p)<2p+1$.
3. **Implementation.** Implement Pollard's rho for the ECDLP (partition the group into three sets by $x\bmod3$ and walk by $R\mapsto R+G$, $2R$, or $R+Q$; detect a cycle with Floyd's tortoise and hare). Measure the number of steps for prime group orders $n$ between $10^3$ and $10^7$, and compare the average with $\sqrt{\pi n/2}$.
4. **Open investigation.** For $p=1009$, compute $N-(p+1)$ for every pair $(a,b)$ with $4a^3+27b^2\neq0$, and normalise by $2\sqrt p$. Plot the histogram. Does it look like a semicircle? Repeat for other primes and for a fixed curve over many primes. (The Sato–Tate conjecture, proved for most curves only in 2011, predicts a distribution here; find out which one by experiment.)

## 78.7 Loose threads

The hook asked how one repeated number broke an unbroken system. The **group** is real and strong: the chord-and-tangent law makes the points of a cubic into an abelian group, the same formulas work over $\mathbb{F}_p$, and the best known attack on the discrete logarithm in a well-chosen curve is generic, needing $\sqrt n$ steps, which is why 256-bit curves rival 3072-bit RSA. The **signature** is fragile: ECDSA hides the private key behind a fresh secret nonce, and repeating that nonce turns the hard problem into two linear equations. Sony's mathematics was perfect and its implementation reused $k$.

What we did not tell you:

- **Why the group law exists at all.** The first mystery points to [Ch. 94](../j-frontier/94-algebraic-geometry.md), where curves of genus one, complex tori and divisors live.
- **How to count points on a 256-bit curve.** Brute force is hopeless at that size, yet curve designers need the exact group order. Schoof's algorithm (1985) computes it in polynomial time by computing $N$ modulo many small primes, using the points of small order, then gluing the answers with the Chinese remainder theorem.
- **Whether randomness can be trusted.** The second mystery points to [Ch. 81](81-pseudorandomness.md).
- **What pairings make possible.** The third mystery points to [Ch. 82](82-zero-knowledge-proofs.md).
- **What quantum computers do here.** Shor's period-finding algorithm ([Ch. 77](77-cryptography.md#773-the-mathematics)) solves the ECDLP too, and with fewer qubits than factoring a comparable RSA modulus. Elliptic curves are not a refuge from quantum computers; lattices may be.
- **Fermat's last theorem.** The proof by Wiles in 1994 runs through elliptic curves over $\mathbb{Q}$ and the modularity of their point counts: the same numbers $N-(p+1)$ that you histogrammed in the challenge, arranged into a single object.

## Further reading

- J. H. Silverman and J. Tate, *Rational Points on Elliptic Curves*, 2nd ed., Springer, 2015. The geometry and arithmetic of the group law, for undergraduates.
- L. C. Washington, *Elliptic Curves: Number Theory and Cryptography*, 2nd ed., CRC Press, 2008. Point counting, discrete-log attacks and pairings, with full proofs.
- D. Hankerson, A. Menezes and S. Vanstone, *Guide to Elliptic Curve Cryptography*, Springer, 2004. The implementation side: coordinates, fast arithmetic, side channels.
- J. Hoffstein, J. Pipher and J. H. Silverman, *An Introduction to Mathematical Cryptography*, 2nd ed., Springer, 2014. Elliptic curves alongside lattices and the rest of public-key cryptography.
