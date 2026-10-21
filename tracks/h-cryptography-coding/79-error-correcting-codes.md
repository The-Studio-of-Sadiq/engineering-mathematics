# Chapter 79 — Error-Correcting Codes: Recovering the Lost

[← Previous](78-elliptic-curves.md) · [Track H](README.md) · [Table of Contents](../../00-front/table-of-contents.md) · [Next →](80-coding-meets-information-theory.md)

> **Level 2 · Track H — Cryptography, Coding, and Digital Security**
> **Builds on:** [Ch. 8](../../part-3-linear-algebra/08-matrices-and-linear-systems.md), [Ch. 9](../../part-3-linear-algebra/09-vector-spaces-and-linear-transformations.md), [Ch. 57](../e-abstract-math/57-number-theory.md), [Ch. 56](../e-abstract-math/56-abstract-algebra.md) · **Leads to:** [Ch. 80](80-coding-meets-information-theory.md), [Ch. 82](82-zero-knowledge-proofs.md) · **Related:** [Ch. 19](../../part-6-probability-statistics/19-probability.md), [Ch. 33](../b-probability-information/33-information-theory.md), [Ch. 52](../d-optimization/52-dynamic-programming.md), [Ch. 60](../e-abstract-math/60-real-analysis.md)
> **Reading time:** ~25 min

## 79.1 The hook

Look at a QR code on a restaurant menu or a concert poster. Many have a logo pasted right over the middle: a coffee cup, a band's face, a company emblem. Those pixels are simply gone, replaced by something that is not data at all. Point your phone at it anyway. It scans on the first try.

Scratch a music CD with a key and it will usually still play. The disc can lose about 4,000 consecutive bits, some 2.5 mm of track, and reproduce the sound *exactly*, not approximately. In 1989 Voyager 2 sent pictures of Neptune from 4.4 billion kilometres away, with a transmitter of about 20 watts, roughly a refrigerator light bulb. By the time the signal reached Earth it was far weaker than the noise in the receiving antennas, and yet the images arrived essentially without a wrong pixel.

None of these systems has a way to ask for the data again. The QR code cannot re-print itself; Neptune is not taking requests. The missing information is *recovered*, computed from what remains. That sounds impossible: how can you reconstruct what was destroyed?

The answer is that nothing useful was destroyed. Before sending, the data was dressed in carefully designed redundancy, so that any small enough damage leaves a pattern that points back, unambiguously, to the original. This chapter is about how to design that redundancy: what makes it efficient, how the receiver finds and fixes errors it cannot see, and where the limits lie.

## 79.2 The big idea

The naive way to protect a bit is to repeat it: send `000` for 0 and `111` for 1, and decode by majority vote. One flipped bit is corrected. But you have tripled the cost, and if you want to survive two flips you must send five copies. There must be a better way.

The key idea is **distance**. Think of every possible received word as a point in a space, and the words you actually send, the **codewords**, as a few special points. Noise moves a codeword a short distance. If the codewords are spread so far apart that the little "clouds" of damaged versions around them never overlap, then any received word lies in exactly one cloud, and the receiver simply snaps it back to the centre. Error correction is **nearest-neighbour decoding**, and code design is the art of packing as many non-overlapping clouds as possible into the space.

An analogy: the NATO spelling alphabet. Over a crackling radio, "B" and "D" and "E" sound alike, but "Bravo", "Delta" and "Echo" do not. The longer words carry redundancy, but it is *arranged* so that the words are far apart in sound. A burst of static that blurs one syllable still leaves a word closer to the right one than to any other.

Three questions organise the chapter.

1. **Geometry.** How far apart must codewords be to correct $t$ errors, and how many codewords can fit?
2. **Algebra.** Searching for the nearest of $2^{k}$ codewords is hopeless when $k=200$. Codes built from linear algebra, and then from polynomials, let the receiver compute the error instead of searching for it.
3. **Erasures versus errors.** A logo on a QR code is damage *at a known place*. A scratch the reader cannot detect is damage at an unknown place. Knowing where damage is turns out to be worth exactly a factor of two.

## 79.3 The mathematics

### Hamming distance

> [!definition] Codes and distance
> A **block code** of length $n$ over an alphabet $A$ is a set $C\subseteq A^n$ of codewords. The **Hamming distance** $d(x,y)$ between two words is the number of positions where they differ. The **minimum distance** of $C$ is $d=\min\{d(x,y):x,y\in C,\ x\neq y\}$. A code with $\lvert C\rvert=q^k$ codewords over a $q$-letter alphabet has **rate** $k/n$.

Hamming distance is a genuine metric ([Ch. 60](../e-abstract-math/60-real-analysis.md#metric-spaces-distance-abstracted) lists it among its examples): it is symmetric, zero only for equal words, and satisfies the **triangle inequality** $d(x,z)\le d(x,y)+d(y,z)$, because any position where $x$ and $z$ differ must be one where $x$ differs from $y$ or $y$ differs from $z$. The **Hamming ball** of radius $t$ around $c$ is the set of all words within distance $t$ of $c$: the "cloud" of the big idea.

![Two drawings of the 3-dimensional cube of bit strings. Left: the repetition code with codewords 000 and 111 at opposite corners, filled in; the three neighbours of 000 are ringed in blue and the three neighbours of 111 in red, so the two balls of radius 1 split the cube exactly. Right: the even-parity code with filled codewords 000, 011, 101 and 110, and the four odd-weight words in grey; every grey word is adjacent to three codewords.](../../images/ch79-hamming-cube.svg)

*Codes as points on the cube $\{0,1\}^3$, where neighbours differ in one bit. Left: distance 3; the radius-1 balls around the two codewords tile the cube, so every single error is corrected. Right: distance 2; a single error always lands on a non-codeword (detected), but that word is equally close to three codewords, so it cannot be corrected.*

### Result 1: distance decides everything

> [!theorem] Detection and correction
> A code with minimum distance $d$ can **detect** any pattern of up to $d-1$ errors, and **correct** any pattern of up to $t=\lfloor (d-1)/2\rfloor$ errors by nearest-codeword decoding.

*Proof.* *Detection.* If $c$ is sent and between 1 and $d-1$ positions change, the received word $r$ satisfies $0<d(c,r)\le d-1$. Since distinct codewords are at least $d$ apart, $r$ is not a codeword, so the receiver sees that something is wrong.

*Correction.* Suppose $c$ is sent and $r$ is received with $d(c,r)\le t$. Let $c'\neq c$ be any other codeword. By the triangle inequality,

$$
d\le d(c,c')\le d(c,r)+d(r,c')\quad\Longrightarrow\quad d(r,c')\ge d-t\ge 2t+1-t=t+1>t\ge d(c,r).
$$

So $c$ is strictly the closest codeword to $r$, and nearest-codeword decoding returns it. $\blacksquare$

Both bounds are sharp. With $d$ errors, $c$ can be turned into another codeword, and nobody can tell. With $t+1$ errors, $r$ can land at distance $t+1$ from $c$ and $d-t-1\le t+1$ from some $c'$, and the decoder may be fooled. The parity code on the right of the figure has $d=2$, so it detects one error and corrects none; the repetition code on the left has $d=3$ and corrects one.

### Result 2: the sphere-packing bound

Large distance costs codewords. How many can a binary code of length $n$ correcting $t$ errors have?

> [!theorem] Hamming bound
> A binary code of length $n$ that corrects $t$ errors has
> $$
> \lvert C\rvert\sum_{i=0}^{t}\binom{n}{i}\ \le\ 2^n .
> $$

*Proof.* A ball of radius $t$ around a word contains the word itself, the $\binom n1$ words differing in one position, the $\binom n2$ differing in two, and so on: $V=\sum_{i\le t}\binom ni$ words. The proof of Result 1 shows that if the code corrects $t$ errors then balls of radius $t$ around distinct codewords are disjoint (a word in two balls would be decoded two ways). So $\lvert C\rvert$ disjoint balls of $V$ words each fit inside the $2^n$ words of the whole space. $\blacksquare$

A code meeting the bound with equality is **perfect**: its balls tile the space with nothing left over. The repetition code of length 3 is perfect: $2\cdot(1+3)=8=2^3$. So is the code of the next section, with $16\cdot(1+7)=128=2^7$, and so is the remarkable binary **Golay code** of length 23 with $2^{12}$ codewords correcting 3 errors: $2^{12}\bigl(1+23+253+1771\bigr)=2^{12}\cdot2^{11}=2^{23}$. Apart from repetition codes and these families, there are no other perfect binary codes; that was a hard theorem of the 1970s.

### Linear codes: let the algebra search for you

Use the alphabet $\mathbb{F}_2=\{0,1\}$ with addition mod 2 (XOR), so that words of length $n$ form the vector space $\mathbb{F}_2^n$ of [Ch. 9](../../part-3-linear-algebra/09-vector-spaces-and-linear-transformations.md).

> [!definition] Linear code
> A **linear $[n,k,d]$ code** is a $k$-dimensional subspace $C\subseteq\mathbb{F}_2^n$ with minimum distance $d$. It has a $k\times n$ **generator matrix** $\mathbf{G}$ whose rows form a basis, so the codewords are $\mathbf{m}\mathbf{G}$ for messages $\mathbf{m}\in\mathbb{F}_2^k$, and an $(n-k)\times n$ **parity-check matrix** $\mathbf{H}$ whose null space is $C$: $\mathbf{c}\in C\iff\mathbf{H}\mathbf{c}^\mathsf{T}=\mathbf{0}$.

Linearity buys two things at once. First, the **minimum distance equals the minimum weight** (number of 1s) of a nonzero codeword: $d(\mathbf x,\mathbf y)=\operatorname{wt}(\mathbf x-\mathbf y)$, and $\mathbf x-\mathbf y$ is again a codeword. Second, the **syndrome** of a received word $\mathbf r=\mathbf c+\mathbf e$ depends only on the error:

$$
\mathbf{s}=\mathbf{H}\mathbf{r}^\mathsf{T}=\mathbf{H}\mathbf{c}^\mathsf{T}+\mathbf{H}\mathbf{e}^\mathsf{T}=\mathbf{H}\mathbf{e}^\mathsf{T}.
$$

The receiver never needs to know what was sent. If $\mathbf e$ has a single 1, in position $j$, then $\mathbf s$ is column $j$ of $\mathbf H$. So: if all columns of $\mathbf H$ are distinct and nonzero, any single error can be located by reading the syndrome, and $d\ge3$. (The null-space picture of [Ch. 9 §9.6](../../part-3-linear-algebra/09-vector-spaces-and-linear-transformations.md#96-null-space-and-column-space) is exactly what is happening.)

Richard Hamming turned this into the most economical possible design in 1950. Take as the columns of $\mathbf H$ *all* the nonzero binary vectors of length 3, in counting order:

$$
\mathbf{H}=\begin{pmatrix}1&0&1&0&1&0&1\\0&1&1&0&0&1&1\\0&0&0&1&1&1&1\end{pmatrix}.
$$

Column $j$ is the number $j$ in binary (least significant bit on top). This is the **Hamming $[7,4,3]$ code**: 7 bits carry 4 data bits, at rate $4/7$ instead of the repetition code's $1/3$, and the syndrome of a single error *is the binary address of the wrong bit*. Positions 1, 2 and 4 (one 1 in their column) hold parity bits $p_1,p_2,p_3$; positions 3, 5, 6, 7 hold data $d_1,\dots,d_4$. Each row of $\mathbf H$ is one parity rule, and the three rules are three overlapping circles.

![Two Venn diagrams of three overlapping circles A, B and C, with the seven bits of a Hamming(7,4) codeword placed in the seven regions: d4 in the centre, d1 in A and B only, d2 in A and C only, d3 in B and C only, and the parity bits p1, p2, p3 in the parts of A, B and C alone. Left: the codeword 0110011 for data 1011, where all three circles contain an even number of 1s. Right: the same word with d3 flipped to 0, shown in red; circles B and C are dashed red and labelled odd, circle A is even.](../../images/ch79-hamming74-venn.svg)

*The Hamming code as three parity circles. Each circle must contain an even number of 1s. A single flipped bit upsets exactly the circles that contain it, and each of the seven regions lies inside a different combination of circles, so the pattern of failures names the culprit. Here B and C fail and A passes: the error is $d_3$, at position $6=110_2$.*

### Cyclic codes and CRCs: codes as polynomials

Write a word $(c_0,\dots,c_{n-1})$ as a polynomial $c(x)=c_0+c_1x+\dots+c_{n-1}x^{n-1}$ over $\mathbb{F}_2$. A **cyclic shift** of the word is multiplication by $x$, followed by reducing $x^n$ to $1$, that is, working in the ring $\mathbb{F}_2[x]/(x^n-1)$ of [Ch. 56](../e-abstract-math/56-abstract-algebra.md#563-the-mathematics). A linear code is **cyclic** if every cyclic shift of a codeword is a codeword. Algebraically, cyclic codes are exactly the ideals of that ring, and every one has the form

$$
C=\{\,m(x)\,g(x)\ :\ \deg m<n-\deg g\,\},
$$

for a **generator polynomial** $g(x)$ that divides $x^n-1$. Over $\mathbb{F}_2$, $x^7+1=(x+1)(x^3+x+1)(x^3+x^2+1)$, and $g(x)=x^3+x+1$ generates a $[7,4,3]$ cyclic code equivalent to Hamming's. Encoding is now polynomial multiplication, which hardware does with a shift register and a few XOR gates.

The same division trick gives the **cyclic redundancy check** (CRC) on every Ethernet frame, USB packet and ZIP file. To protect a message $m(x)$, append the remainder $r(x)=x^{r}m(x)\bmod g(x)$, where $r=\deg g$; the transmitted word $x^rm(x)+r(x)$ is then divisible by $g(x)$. With $g=x^3+x+1$ (bits `1011`) and message `1101`, the remainder of `1101000` is `001`, and `1101001` is sent. The receiver divides by $g$ and accepts only remainder zero. If $g$ has a nonzero constant term, every **burst** of errors confined to $r$ consecutive positions is detected: such an error is $x^ib(x)$ with $\deg b<r$ and $b\neq0$, and $g$ can divide neither $x^i$ (it is not a power of $x$) nor $b$ (whose degree is too small). CRC-32 detects every burst up to 32 bits long. It *detects*; it does not correct. For that we need bigger fields.

### Reed–Solomon codes: data as a polynomial

Bits flip one at a time, but scratches, dust and logos destroy whole **symbols**: bytes, or blocks of pixels. Reed and Solomon (1960) chose a large alphabet, a finite field $\mathbb{F}_q$ such as $\mathbb{F}_{2^8}$ whose elements are bytes ([Ch. 57](../e-abstract-math/57-number-theory.md#573-the-mathematics)), and a beautifully simple idea.

> [!definition] Reed–Solomon code
> Fix distinct points $\alpha_1,\dots,\alpha_n\in\mathbb{F}_q$ (so $n\le q$) and $k\le n$. A message $(m_0,\dots,m_{k-1})\in\mathbb{F}_q^k$ becomes the polynomial $f(x)=m_0+m_1x+\dots+m_{k-1}x^{k-1}$, and its codeword is the list of values
> $$
> \bigl(f(\alpha_1),\,f(\alpha_2),\,\dots,\,f(\alpha_n)\bigr)\in\mathbb{F}_q^n .
> $$

Send $n$ samples of a curve that is determined by $k$ of them. The extra $n-k$ samples are the redundancy.

![A smooth cubic curve f(x) through seven evaluation points x = 1 to 7. Five points are solid blue dots on the curve. The point at x = 3 is hollow with a question mark: erased, at a known position. The point at x = 6 is a red cross well above the curve, with a dashed arrow down to the true value: corrupted, at an unknown position.](../../images/ch79-reed-solomon.svg)

*Reed–Solomon, drawn over the real numbers for intuition (real codes use a finite field). The message is a cubic ($k=4$) and seven values are sent. Any four correct values pin the cubic down, so up to three erasures can be filled in; an undetected wrong value, like the red cross, costs twice as much to fix.*

### Result 3: Reed–Solomon codes are as good as possible

How large can the minimum distance of *any* code be? There is a simple ceiling.

> [!theorem] Singleton bound, and Reed–Solomon meets it
> (a) Any code of length $n$ with $q^k$ codewords over a $q$-letter alphabet has $d\le n-k+1$.
> (b) The Reed–Solomon code has $d=n-k+1$ exactly.

*Proof.* (a) Delete the last $d-1$ positions of every codeword. Two distinct codewords differ in at least $d$ places, so they still differ after deleting $d-1$ of them: the shortened words, of length $n-d+1$, are all distinct. There are $q^k$ of them and only $q^{\,n-d+1}$ possible words of that length, so $k\le n-d+1$.

(b) The code is linear (adding messages adds polynomials and their values), so $d$ is the minimum weight of a nonzero codeword. A nonzero polynomial $f$ of degree at most $k-1$ has at most $k-1$ roots in a field, so at most $k-1$ of the values $f(\alpha_i)$ are zero, and at least $n-k+1$ are nonzero. Hence every nonzero codeword has weight $\ge n-k+1$, and with (a), $d=n-k+1$. $\blacksquare$

Codes meeting the Singleton bound are called **maximum distance separable** (MDS). Combined with Result 1 and a small extension of it to erasures, Result 3 says what a Reed–Solomon decoder can do: with $s$ erasures (positions known) and $e$ errors (positions unknown), it recovers the message whenever

$$
2e+s\le n-k .
$$

An erasure costs one redundant symbol, an error costs two: knowing *where* damage is halves its price. With only erasures, any $k$ surviving values determine $f$ by Lagrange interpolation, which is the uniqueness of the interpolating polynomial in [Ch. 22 §22.5](../../part-7-numerical-computational/22-numerical-methods.md#225-interpolation) carried over to a finite field.

That answers the hook. A QR code splits its data into blocks, each a Reed–Solomon codeword over $\mathbb{F}_{2^8}$. At the highest error-correction level, H, the smallest code has 26 bytes of which 17 are redundant, enough to fix 8 wrong bytes, about 30 % of the block. A designer pasting a logo over the centre is spending part of that budget on purpose. The CD uses two Reed–Solomon codes over bytes with the data shuffled between them, so that a long scratch becomes many short, correctable damages spread across many codewords. Deep-space links use a $(255,223)$ Reed–Solomon code, correcting 16 byte errors in every 255, wrapped around a convolutional code.

> [!mystery] How good can a code be?
> The repetition code drives errors to zero only by driving the rate to zero. Hamming's code has rate $4/7$ but corrects only one error in seven. Is there a trade-off you cannot beat? Shannon proved in 1948 that for a channel flipping bits with probability $p$, every rate below $1-H_b(p)$ (the binary entropy of [Ch. 33](../b-probability-information/33-information-theory.md)) is achievable with error probability *as small as you like*, and no rate above it is. His proof used random codes and was not constructive at all. What that capacity is, why random codes reach it, and how LDPC and polar codes finally got there in practice, is [Ch. 80](80-coding-meets-information-theory.md).

## 79.4 The algorithm

Two decoders, one per idea. The Hamming decoder computes a syndrome and flips one bit. The Reed–Solomon decoder works over the prime field $\mathbb{F}_{929}$ (the field used by the PDF417 barcode on many boarding passes and driving licences) and recovers the message by interpolation.

```python
import numpy as np
from itertools import combinations

# --- Hamming(7,4): column j of H is j in binary, positions p1 p2 d1 p3 d2 d3 d4 ---
H = np.array([[(j >> i) & 1 for j in range(1, 8)] for i in range(3)])
DATA_POS = [2, 4, 5, 6]                       # 0-based positions of d1..d4

def ham_encode(d):
    c = np.zeros(7, dtype=int); c[DATA_POS] = d
    for i, pos in enumerate([0, 1, 3]):        # parity bit i fixes row i of H
        c[pos] = H[i] @ c % 2
    return c

def ham_decode(r):
    s = H @ r % 2
    j = s[0] + 2 * s[1] + 4 * s[2]             # syndrome = position of the error
    r = r.copy()
    if j: r[j - 1] ^= 1
    return r[DATA_POS]

c = ham_encode([1, 0, 1, 1]); r = c.copy(); r[5] ^= 1
print(c, r, ham_decode(r))                     # [0 1 1 0 0 1 1] [0 1 1 0 0 0 1] [1 0 1 1]

# --- Reed-Solomon over the prime field GF(929): evaluate, then interpolate ---
P = 929
def rs_encode(msg, n):                         # msg = k coefficients, codeword = f(1..n)
    return [sum(m * pow(x, i, P) for i, m in enumerate(msg)) % P for x in range(1, n + 1)]

def interpolate(pts, x0):                      # Lagrange: value at x0 of the poly through pts
    total = 0
    for xi, yi in pts:
        num = den = 1
        for xj, _ in pts:
            if xj != xi: num, den = num * (x0 - xj) % P, den * (xi - xj) % P
        total += yi * num * pow(den, -1, P)
    return total % P

def rs_decode_bruteforce(word, k):             # try k-subsets; exponential, for illustration only
    n = len(word)
    known = [(x, y) for x, y in zip(range(1, n + 1), word) if y is not None]
    t = (len(known) - k) // 2                  # errors we can afford: 2e + erasures <= n - k
    for sub in combinations(known, k):
        cand = [interpolate(sub, x) for x in range(1, n + 1)]
        if sum(c != y for c, y in zip(cand, word) if y is not None) <= t:
            return cand
    return None

cw = rs_encode([3, 1, 4, 1], n=7)              # k = 4, n = 7, d = 4
print(cw)                                      # [9, 29, 69, 135, 233, 369, 549]
bad = cw.copy(); bad[2] = None; bad[5] = (bad[5] + 100) % P   # one erasure, one error
print(rs_decode_bruteforce(bad, 4) == cw)      # True: 2*1 + 1 <= 7 - 4
```

The Hamming part reproduces the Venn figure: data `1011` becomes `0110011`, the flipped sixth bit gives syndrome $6$, and the data comes back. The Reed–Solomon part erases one symbol and corrupts another, a damage pattern with $2e+s=3=n-k$, right at the limit, and still recovers the codeword exactly.

**Complexity.** Hamming decoding is a matrix–vector product: $O(n\log n)$ bit operations for a length-$n$ Hamming code. The Reed–Solomon decoder above is deliberately naive. It tries up to $\binom{n}{k}$ subsets, which is fine for $n=7$ and absurd for the $(255,223)$ deep-space code, where $\binom{255}{223}\approx5\times10^{40}$. Real decoders find the errors algebraically. The **Berlekamp–Welch** idea fits in three lines. Let $E(x)$ be the unknown **error-locator** polynomial, monic of degree $e$, vanishing exactly at the error positions. Then at *every* evaluation point,

$$
y_i\,E(\alpha_i)=f(\alpha_i)\,E(\alpha_i)=Q(\alpha_i),\qquad Q:=fE,\ \deg Q<k+e,
$$

because where $y_i$ is wrong, $E(\alpha_i)=0$ kills both sides. These are $n$ *linear* equations in the $k+2e\le n$ unknown coefficients of $E$ and $Q$. Solve them by Gaussian elimination ([Ch. 8 §8.7](../../part-3-linear-algebra/08-matrices-and-linear-systems.md#87-gaussian-elimination)) in $O(n^3)$, then divide: $f=Q/E$. The nonlinear problem "where are the errors?" has become linear algebra. Production decoders (Berlekamp–Massey with Chien search and Forney's formula) do the same job in $O(n^2)$ or better.

**Failure cases.** Every decoder has a breaking point, and its behaviour past that point matters.

- *Too many errors.* With $t+1$ errors, a nearest-codeword decoder may confidently return the *wrong* codeword. For the Hamming code, two errors always produce a nonzero syndrome pointing at a third, innocent bit; the "correction" makes three errors. A system that must know when it failed adds an overall parity bit (the extended $[8,4,4]$ code corrects one error *and* detects two).
- *Burst errors.* A scratch destroys many neighbouring bits. Bit-level codes are helpless against it; symbol-level codes plus **interleaving** (spreading each codeword across the medium) turn one long burst into many small, correctable damages.
- *The wrong noise model.* Codes designed for independent flips can do badly on channels where noise comes in bursts or fades. The design starts from the channel, not from the code.
- *Detection is not security.* A CRC catches accidents, but anyone can recompute it after tampering. Protection against an adversary needs the keyed authentication tags of [Ch. 77](77-cryptography.md#773-the-mathematics), not a checksum.

## 79.5 The experiment

First, the formulas. A 4-bit message sent unprotected over a channel that flips each bit with probability $x$ survives only if all four bits survive. Sent as a Hamming codeword, it survives if at most one of the seven bits flips (using the binomial distribution of [Ch. 19 §19.6](../../part-6-probability-statistics/19-probability.md#196-probability-distributions)).

```desmos-graph
left=0; right=0.25; bottom=0; top=0.7
---
y=1-(1-x)^4|red
y=1-(1-x)^7-7x(1-x)^6|blue
y=1-((1-x)^3+3x(1-x)^2)^4|green|dashed
```

*Message error probability against bit-flip probability $x$: uncoded (red), Hamming(7,4) (blue), triple repetition of each bit (green, dashed). At $x=0.01$ the Hamming code cuts the failure rate from 3.9 % to 0.2 %.*

Change the exponents to explore. Replace the blue curve by $1-(1-x)^{15}-15x(1-x)^{14}$, the $[15,11,3]$ Hamming code, which carries 11 bits instead of 4; compare it with $1-(1-x)^{11}$. The longer code is cheaper (rate $11/15$ instead of $4/7$), but its curve rises much faster: one correction per 15 bits is spread thinner than one per 7. Where on the $x$-axis does that difference start to matter?

Then simulate. The script sends many random 4-bit messages through a binary symmetric channel, decodes them three ways, and plots the measured error rates.

```python
import numpy as np
import matplotlib.pyplot as plt

FLIP_PROBS = np.linspace(0.001, 0.2, 25)   # KNOB: channel bit-flip probabilities to test
BLOCKS = 20000                             # KNOB: 4-bit messages sent per probability
SEED = 79

rng = np.random.default_rng(SEED)
H = np.array([[(j >> i) & 1 for j in range(1, 8)] for i in range(3)])   # Hamming(7,4)
DATA = [2, 4, 5, 6]
G = np.zeros((4, 7), dtype=int)            # generator: encode each unit data vector
for row, pos in enumerate(DATA):
    c = np.zeros(7, dtype=int); c[pos] = 1
    c[[0, 1, 3]] = H @ c % 2               # parity bits are the syndrome of the data part
    G[row] = c

def run(p):
    d = rng.integers(0, 2, (BLOCKS, 4))
    raw = (rng.random(d.shape) < p).any(axis=1).mean()                     # uncoded: any flip
    rep = d.repeat(3, axis=1) ^ (rng.random((BLOCKS, 12)) < p)             # repetition x3
    rep_dec = rep.reshape(BLOCKS, 4, 3).sum(axis=2) >= 2
    rep_err = (rep_dec != d).any(axis=1).mean()
    r = (d @ G % 2) ^ (rng.random((BLOCKS, 7)) < p)                        # Hamming(7,4)
    j = (H @ r.T % 2).T @ np.array([1, 2, 4])                             # syndrome -> position
    r[np.arange(BLOCKS)[j > 0], j[j > 0] - 1] ^= 1
    ham_err = (r[:, DATA] != d).any(axis=1).mean()
    return raw, rep_err, ham_err

res = np.array([run(p) for p in FLIP_PROBS])
for lab, col in zip(["uncoded (rate 1)", "repetition x3 (rate 1/3)", "Hamming(7,4) (rate 4/7)"], res.T):
    plt.semilogy(FLIP_PROBS, np.maximum(col, 1e-5), "o-", label=lab)
plt.xlabel("bit-flip probability p"); plt.ylabel("message error rate")
plt.legend(); plt.grid(True, which="both"); plt.show()
print("at p = %.3f:" % FLIP_PROBS[5], res[5])
```

At $p\approx0.042$ the run prints roughly `0.155 0.022 0.031`: uncoded messages fail 15.5 % of the time, Hamming-coded ones 3.1 %, and triple repetition 2.2 %. Repetition wins on reliability, but it sends 12 bits per message where Hamming sends 7. Things to try:

- **Equalise the cost.** Compare codes at the same *rate*, not the same message. Add a $[15,11]$ Hamming code (column $j$ of $\mathbf H$ is $j$ in 4-bit binary) and send 11-bit messages. Which wins per transmitted bit?
- **Add a burst.** Instead of independent flips, flip three *consecutive* bits of each block with probability `p`. Watch the Hamming code fall apart, then interleave: send bit 1 of ten codewords, then bit 2 of each, and so on.
- **Find the wall.** Push `FLIP_PROBS` up to 0.5. All three curves converge on the same failure rate as $p\to\tfrac12$. Why can no code, however clever, help at $p=\tfrac12$, and what happens if you push past it to $p=0.9$?

## 79.6 The challenge

No answers are given.

1. **Intuition.** A QR code at level H can survive about 30 % damage, yet a logo covering only 20 % of the symbol can still make it unreadable if placed badly. Using erasures, errors and the structure of the code (data split into several Reed–Solomon blocks, plus fixed finder patterns in three corners), explain why *where* the damage falls matters as much as how much there is.
2. **Derivation.** Prove the **Gilbert–Varshamov bound**: there exists a binary code of length $n$ and minimum distance $d$ with at least $2^n/\sum_{i=0}^{d-1}\binom ni$ codewords. (Hint: build the code greedily; when you can add no more words, the balls of radius $d-1$ around the codewords must cover the whole space.) Compare with the Hamming bound for $n=15$, $d=3$.
3. **Implementation.** Implement Berlekamp–Welch over $\mathbb{F}_{929}$: set up the $n\times(k+2e)$ linear system, solve it by Gaussian elimination mod 929, and divide polynomials to get $f$. Test it on the $(7,4)$ example above, then on $n=100$, $k=50$ with 25 random errors, and compare its running time with the brute-force decoder.
4. **Open investigation.** For $n$ up to about 12, search exhaustively (or with a clever backtracking search) for the largest binary code with minimum distance 3, and compare with the Hamming bound and the Gilbert–Varshamov bound. How fast does the exhaustive search blow up? Can a random-restart greedy search match the known optimum for $n=10$ and $n=11$ (72 and 144 codewords)?

## 79.7 Loose threads

The hook asked how a QR code with a logo, a scratched CD or a faint signal from Neptune can recover data that is literally gone. The answer is geometry and algebra. Codewords are placed far apart in Hamming distance, so that damaged words still point home (Result 1), at a cost bounded by sphere packing (Result 2). Linear codes compute the error from a syndrome instead of searching for it. Reed–Solomon codes turn a block of data into a polynomial and send extra samples of it; because a low-degree polynomial has few roots, any $k$ intact samples determine it, and the code meets the Singleton bound exactly (Result 3). Erasures, like a logo at a known place, cost half as much as errors.

What we did not tell you:

- **The limit.** The first mystery points to [Ch. 80](80-coding-meets-information-theory.md): Shannon's capacity, and the codes that approach it.
- **Convolutional codes and Viterbi decoding.** Instead of blocks, a convolutional encoder passes an endless bit stream through a shift register. Its decoder finds the most likely path through a trellis of states, a shortest-path dynamic program ([Ch. 52](../d-optimization/52-dynamic-programming.md)).
- **Decoding beyond half the distance.** Result 1 seems to forbid correcting more than $(d-1)/2$ errors. But if the decoder may output a short *list* of candidates, Guruswami and Sudan (1999) showed that Reed–Solomon codes can be decoded far beyond that radius, using curves through the received points.

> [!mystery] Are random codes the best possible?
> The Gilbert–Varshamov bound of the challenge says that, as $n\to\infty$ with $d/n=\delta$ fixed, there exist binary codes of rate at least $1-H_b(\delta)$, and a random linear code achieves it. Seventy years later, nobody has found an explicit family of binary codes that beats it, and nobody has proved that none exists. Upper bounds from linear programming leave a wide gap. Over larger alphabets, codes built from algebraic curves *do* beat the corresponding bound, a stunning result from 1982. For binary codes, the true best rate is an open problem.

> [!mystery] Checking a million symbols by reading twenty
> A Reed–Solomon codeword has a rigid structure: its values come from a low-degree polynomial. Suppose a prover sends you a list of a million field values and claims it is (close to) such a codeword. Could you check that claim by reading only a few dozen randomly chosen positions? Surprisingly, yes, and this "proximity test" is the engine of modern succinct proofs, which let a blockchain verify a long computation in milliseconds. How a few random queries can certify a huge object, and how this connects to proving things without revealing them, is [Ch. 82](82-zero-knowledge-proofs.md).

## Further reading

- R. W. Hamming, *Coding and Information Theory*, 2nd ed., Prentice Hall, 1986. The inventor's own gentle introduction.
- R. Roth, *Introduction to Coding Theory*, Cambridge University Press, 2006. Linear, cyclic and Reed–Solomon codes with complete algebraic decoders.
- V. Guruswami, A. Rudra and M. Sudan, *Essential Coding Theory*, draft textbook, freely available online. Bounds, list decoding and the modern algorithmic view.
- F. J. MacWilliams and N. J. A. Sloane, *The Theory of Error-Correcting Codes*, North-Holland, 1977. The classic encyclopedic reference.
