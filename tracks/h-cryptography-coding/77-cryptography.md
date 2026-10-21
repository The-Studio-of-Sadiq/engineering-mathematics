# Chapter 77 — Cryptography: Mathematics Against Adversaries

[← Previous](README.md) · [Track H](README.md) · [Table of Contents](../../00-front/table-of-contents.md) · [Next →](78-elliptic-curves.md)

> **Level 2 · Track H — Cryptography, Coding, and Digital Security**
> **Builds on:** [Ch. 19](../../part-6-probability-statistics/19-probability.md), [Ch. 57](../e-abstract-math/57-number-theory.md), [Ch. 30](../a-algorithms/30-limits-of-computation.md), [Ch. 33](../b-probability-information/33-information-theory.md) · **Leads to:** [Ch. 78](78-elliptic-curves.md), [Ch. 81](81-pseudorandomness.md), [Ch. 82](82-zero-knowledge-proofs.md) · **Related:** [Ch. 46](../c-game-theory/46-adversarial-decision-making.md), [Ch. 32](../b-probability-information/32-counting-random-worlds.md), [Ch. 56](../e-abstract-math/56-abstract-algebra.md)
> **Reading time:** ~25 min

## 77.1 The hook

In the winter of 1942, with German armies outside Moscow, the Soviet office that printed one-time pads ran short of time. A one-time pad is a booklet of random numbers. Sender and receiver hold identical copies; each page is added to one message and then destroyed. Used correctly it is the only cipher in history that has been *proved* unbreakable. Under wartime pressure, some pages were printed twice and issued to different users.

The American codebreakers of the Venona project, working from 1943 on, noticed. Over the next four decades they read parts of nearly three thousand Soviet messages, among them traffic that exposed spies inside the Manhattan Project. The unbreakable cipher broke. Nobody beat its mathematics; someone *used it twice*.

The lesson keeps returning. Ciphers in Wi-Fi routers, payment cards and game consoles are rarely defeated by mathematics; they fall to a reused number, a predictable seed or a timing leak. And yet mathematics is what makes modern security possible at all. Your browser has never met your bank and shares no secret with it, yet over wires anyone may tap, the two agree within milliseconds on a key the eavesdropper cannot compute.

So there are three questions, and this chapter is about all of them:

- What does it even *mean* for a cipher to be secure, stated precisely enough to prove?
- How can a short key protect a long message, when the one provably perfect cipher needs a key as long as the message?
- How can two strangers agree on a secret in public?

## 77.2 The big idea

Cryptography is **mathematics with an opponent**. In most of this book the world is indifferent: noise is random, a circuit does not try to fool you. Here the noise is a clever person, Eve, who reads everything you send, knows exactly which algorithm you use, and will exploit any pattern you leave. A security claim is therefore a statement about *every possible strategy* Eve could follow, which is why it has to be framed as a game.

The game has three rules that took cryptographers a century to agree on.

1. **Eve knows the system.** Only the key is secret. This is **Kerckhoffs' principle** (1883), the cryptographer's version of publishing your mixed strategy in [Ch. 46](../c-game-theory/46-adversarial-decision-making.md). A cipher whose security rests on a secret design is broken the day the design leaks.
2. **Say exactly what Eve can do.** Only listen? Get messages of her choice encrypted? Submit forged ciphertexts and watch your reaction? Each is a different attack model.
3. **Say exactly what "winning" means.** If Eve learns *one bit* about the message, say "buy" versus "sell", she has won.

With the game fixed, there are three levels of security, and they form the spine of this chapter.

- **Perfect secrecy.** The ciphertext carries no information at all about the message, even for an Eve with infinite computing power. Achievable, but the key must be as long as everything you will ever send.
- **Computational security.** Eve could in principle learn something, but only by a computation so large it will not finish before the Sun burns out. This is how short keys protect long messages.
- **Public-key security.** Some functions are easy to compute but hard to invert unless you hold a hidden "trapdoor". Bob publishes the lock and keeps the key; anyone can snap the padlock shut, only Bob can open it.

![Two panels. Top: symmetric encryption, where Alice encrypts with a key k that was delivered in advance along a dashed purple path, sends the ciphertext over a public channel that Eve taps, and Bob decrypts with the same k. Bottom: public-key encryption, where Bob posts his public key pk in a public directory, Alice encrypts with pk, Eve sees both the ciphertext and pk, and Bob decrypts with a private key sk that never leaves him.](../../images/ch77-alice-bob-eve.svg)

*The two architectures. In symmetric cryptography the hard problem is getting k to Bob in the first place. Public-key cryptography removes that problem by splitting the key in two, at the price of much slower arithmetic; real systems use public keys only to agree on a symmetric key.*

The analogy to keep in mind: a symmetric cipher is a **safe with two identical keys**, one for Alice and one for Bob. A public-key cipher is a **mailbox with a slot**: anyone can drop a letter in, but only the owner has the key to the back door.

> [!warning] Never roll your own crypto
> This chapter separates two things that are easy to confuse. The *mathematics* of a cryptosystem can be correct, proved and beautiful while an *implementation* of it is trivially broken: by a reused nonce, a weak random-number generator, a non-constant-time comparison, a padding error message, or a key left in memory. Every code sample in this chapter uses toy parameters and textbook constructions chosen for clarity, and every one of them is insecure. For real systems use a vetted, maintained library (libsodium, your platform's TLS stack, a FIPS-validated module) through its highest-level interface, and follow current guidance from bodies such as NIST. Designing ciphers, protocols or even "small tweaks" is a job for specialists, and their work is trusted only after years of public attack.

## 77.3 The mathematics

### Encryption schemes

> [!definition] Symmetric encryption scheme
> A scheme is three algorithms. $\operatorname{Gen}$ outputs a key $k$ from a key space $\mathcal{K}$. $\operatorname{Enc}_k(m)$ turns a message $m\in\mathcal{M}$ into a ciphertext $c$; it may use fresh randomness. $\operatorname{Dec}_k(c)$ returns a message. **Correctness** demands $\operatorname{Dec}_k(\operatorname{Enc}_k(m))=m$ for every $k$ and $m$.

Treat the message as a random variable $M$ (Eve has some prior belief about what Alice might say) and the key as an independent random variable $K$, chosen by $\operatorname{Gen}$. Then the ciphertext $C=\operatorname{Enc}_K(M)$ is random too, and we can ask what seeing it does to Eve's beliefs. That is a conditional-probability question, and [Ch. 19](../../part-6-probability-statistics/19-probability.md#193-bayes-theorem) has the tools.

> [!definition] Perfect secrecy (Shannon, 1949)
> A scheme is **perfectly secret** if for every prior distribution of $M$, every message $m$ and every ciphertext $c$ with $P(C=c)>0$,
> $$
> P(M=m\mid C=c)=P(M=m).
> $$
> Seeing the ciphertext does not change Eve's beliefs at all. Equivalently, $M$ and $C$ are independent, so their mutual information ([Ch. 33](../b-probability-information/33-information-theory.md)) is zero.

### Result 1: the one-time pad is perfectly secret

Let messages, keys and ciphertexts all be $n$-bit strings, let $\operatorname{Gen}$ choose $K$ uniformly from all $2^n$ strings, and encrypt by bitwise XOR (addition mod 2, written $\oplus$):

$$
\operatorname{Enc}_k(m)=m\oplus k,\qquad \operatorname{Dec}_k(c)=c\oplus k .
$$

Correctness is one line: $(m\oplus k)\oplus k=m\oplus(k\oplus k)=m\oplus 0=m$.

> [!theorem] The one-time pad is perfectly secret
> For any distribution of $M$ independent of the uniform key $K$, $P(M=m\mid C=c)=P(M=m)$.

*Proof.* Fix $m$ and $c$. Given $M=m$, the ciphertext equals $c$ exactly when $K=m\oplus c$, a single key. Because $K$ is uniform and independent of $M$,

$$
P(C=c\mid M=m)=P(K=m\oplus c)=2^{-n},
$$

the same number for every $m$. By the law of total probability, $P(C=c)=\sum_{m'}P(C=c\mid M=m')P(M=m')=2^{-n}$. Bayes' theorem then gives

$$
P(M=m\mid C=c)=\frac{P(C=c\mid M=m)\,P(M=m)}{P(C=c)}=\frac{2^{-n}\,P(M=m)}{2^{-n}}=P(M=m). \qquad\blacksquare
$$

Read the proof again and see where it used its assumptions. It needed $K$ **uniform** (otherwise $P(K=m\oplus c)$ depends on $m$), **independent of $M$**, and **used once**. Break the last condition and the pad cancels: if $c_1=m_1\oplus k$ and $c_2=m_2\oplus k$, then

$$
c_1\oplus c_2=m_1\oplus m_2 ,
$$

which no longer involves the key at all. The XOR of two English texts is full of structure, and guessing a common word in one reveals the other at the same position. That is all Venona needed, and it is the experiment in §77.5.

### Result 2: perfect secrecy needs a key as long as the message

Could a cleverer scheme achieve perfect secrecy with a shorter key? Shannon proved that it cannot.

> [!theorem] Shannon's key-length bound
> If a scheme is perfectly secret, then $\lvert\mathcal{K}\rvert\ge\lvert\mathcal{M}\rvert$.

*Proof.* Take $M$ uniform on $\mathcal{M}$ and pick any ciphertext $c$ that can occur. Let $D(c)=\{\operatorname{Dec}_k(c):k\in\mathcal{K}\}$ be the set of messages that $c$ could decrypt to under some key. It has at most $\lvert\mathcal{K}\rvert$ elements. Suppose $\lvert\mathcal{K}\rvert<\lvert\mathcal{M}\rvert$. Then some message $m^\ast$ lies outside $D(c)$, and by correctness $c$ can never be an encryption of $m^\ast$, so $P(M=m^\ast\mid C=c)=0$. But $P(M=m^\ast)=1/\lvert\mathcal{M}\rvert>0$. This contradicts perfect secrecy. $\blacksquare$

The counting argument covers every cipher anyone will ever invent: a 128-bit key cannot perfectly protect a 129-bit message. That is why the one-time pad lived on only in diplomatic bags and hotlines, and why the rest of cryptography had to change the rules of the game.

### Computational security: limiting Eve's budget

The way out is to stop defending against an Eve with unlimited computing power. Instead we ask that every Eve who runs for at most $t$ steps succeeds with probability at most $\varepsilon$ above blind guessing. The scheme is then called **$(t,\varepsilon)$-secure**. Asymptotically, with a security parameter $\lambda$ (roughly the key length), we demand that every adversary running in time polynomial in $\lambda$ has an advantage that is **negligible**: smaller than $1/\lambda^c$ for every constant $c$ once $\lambda$ is large enough, as $2^{-\lambda}$ is.

How large is large enough? A 128-bit key has $2^{128}\approx3.4\times10^{38}$ values. A machine testing $10^{18}$ keys per second, far beyond any that exists, would need about $1.1\times10^{13}$ years to try them all, some 800 times the age of the universe. In 1998 a purpose-built machine searched the 56-bit key space of the old DES standard in 56 hours. Key length is not a detail; each extra bit doubles Eve's work.

Two caveats make computational security a different kind of mathematics from Results 1 and 2.

- **It rests on assumptions.** We cannot prove that any practical cipher needs $2^{128}$ steps to break, only that "if Eve breaks this scheme, she also solves problem X", where X has resisted the world for decades. This is a **reduction**, the tool [Ch. 30](../a-algorithms/30-limits-of-computation.md#reductions-undecidability-spreads) uses to spread undecidability, now spreading hardness.
- **It needs average-case hardness.** NP-hardness ([Ch. 30](../a-algorithms/30-limits-of-computation.md#efficient-computation-p-and-np)) says that *some* instances are hard. A cipher needs *random* instances, those made by a random key, to be hard almost always. The minimal object that delivers this is a **one-way function**: a function $f$ that is easy to compute, but such that for a random $x$ every efficient algorithm, given $f(x)$, finds some preimage only with negligible probability. If one-way functions exist then $\mathrm{P}\neq\mathrm{NP}$, because inverting $f$ is an NP search problem; the converse is not known. This is the precise sense in which modern cryptography is a bet against $\mathrm{P}=\mathrm{NP}$, and a stronger bet than that.

The engine of symmetric cryptography is the **pseudorandom generator** (PRG): a deterministic function $G$ that stretches a short random seed $k$ into a long string $G(k)$ that no efficient test can tell apart from truly random bits. Replace the one-time pad by $G(k)$ and you have a **stream cipher**: $c=m\oplus G(k)$. Its security proof is a reduction. If some Eve could learn anything about $m$ from $c$, then running her on $m\oplus s$ for a challenge string $s$ would tell you whether $s$ came from $G$ or from fair coins, because against truly random $s$ we already proved she learns nothing. So a break of the cipher *is* a break of the generator.

> [!mystery] Do one-way functions exist?
> Every computational guarantee in this chapter is conditional. Remarkably, a single assumption would suffice for almost all of symmetric cryptography: Håstad, Impagliazzo, Levin and Luby proved in 1999 that a pseudorandom generator can be built from *any* one-way function. But nobody can prove that a single one-way function exists, and doing so would settle $\mathrm{P}\neq\mathrm{NP}$ along the way. What "indistinguishable from random" means exactly, and how the construction squeezes randomness out of mere one-wayness, is [Ch. 81](81-pseudorandomness.md).

### Block ciphers, and why encryption must be randomised

Most data is encrypted with a **block cipher**, a keyed permutation $E_k$ of fixed-size blocks. AES, standardised in 2001, permutes 128-bit blocks under a 128-, 192- or 256-bit key. The security goal is that $E_k$, for a random hidden $k$, is indistinguishable from a permutation chosen uniformly at random from all $(2^{128})!$ of them, a **pseudorandom permutation**.

A block cipher alone is not an encryption scheme for long messages; you need a **mode of operation**. The obvious one, encrypting each block separately (electronic codebook, ECB), is a famous mistake. Equal plaintext blocks give equal ciphertext blocks, so every repetition in the data survives encryption.

![Three panels of a 64 by 48 pixel image with a disc and a square with a hole. Left: the plaintext. Middle: ECB encryption with a toy keyed permutation, where the background is replaced by vertical stripes but the outlines of the disc and the square are still plainly visible. Right: counter-mode encryption with the same toy cipher, which looks like uniform noise.](../../images/ch77-ecb-vs-ctr.svg)

*A toy 16-bit "block cipher" (a random permutation of $2^{16}$ values) applied to an image four pixels at a time. In ECB mode every block of uniform colour maps to the same ciphertext block, so the shapes survive. In counter (CTR) mode each block is XORed with $E_k(\text{counter})$, a fresh pseudorandom pad, and the structure disappears.*

The image shows a general fact, and the right definition of security exposes it. In the game of **indistinguishability under chosen-plaintext attack** (IND-CPA), Eve may ask for encryptions of any messages she likes. Then she submits two messages $m_0,m_1$ of equal length, receives $c^\ast=\operatorname{Enc}_k(m_b)$ for a secret fair coin $b$, may ask for more encryptions, and finally guesses $b$. The scheme is IND-CPA secure if no efficient Eve guesses correctly with probability non-negligibly above $\tfrac12$.

> [!theorem] Deterministic encryption is never IND-CPA secure
> If $\operatorname{Enc}_k$ is a deterministic function of $m$, there is an Eve who wins the IND-CPA game with probability 1.

*Proof.* Eve picks any two distinct messages $m_0\neq m_1$ and asks the encryption oracle for $c_0=\operatorname{Enc}_k(m_0)$. She submits $(m_0,m_1)$ and receives $c^\ast$. Since encryption is deterministic, $c^\ast=c_0$ if $b=0$. If $b=1$, then $c^\ast=\operatorname{Enc}_k(m_1)\neq c_0$, because by correctness two messages with the same ciphertext under one key would decrypt to the same thing. So Eve outputs 0 if $c^\ast=c_0$ and 1 otherwise, and is always right. $\blacksquare$

This three-line proof condemns ECB, and it equally condemns the textbook RSA of [Ch. 57](../e-abstract-math/57-number-theory.md#573-the-mathematics), where $c=m^e\bmod n$ depends on nothing but $m$. Secure encryption must use fresh randomness or a never-repeating **nonce** every time. That is why counter mode feeds a counter into the cipher, why real RSA encryption pads the message with random bytes (OAEP), and why, as Venona showed, *reusing* the fresh value destroys everything.

Confidentiality is also not enough. CTR mode is **malleable**: flipping ciphertext bit $i$ flips plaintext bit $i$, so Eve can turn "pay 100" into "pay 900" blind. Real systems use **authenticated encryption** (AES-GCM, ChaCha20-Poly1305), which adds a tag no forger can compute and rejects tampered ciphertexts.

### Public-key cryptography

The key-distribution problem in the figure is solved by functions with a trapdoor. You have met the two classical ones in [Ch. 57](../e-abstract-math/57-number-theory.md#573-the-mathematics), which proves why they work; here we look at them as security engineers.

- **Diffie–Hellman** (1976). In a cyclic group with generator $g$, Alice sends $g^a$, Bob sends $g^b$, both compute $g^{ab}$. Getting $g^{ab}$ from $g,g^a,g^b$ is the **computational Diffie–Hellman** problem; the stronger **decisional** assumption says $g^{ab}$ even *looks* random to Eve. Unauthenticated, the exchange falls to a "man in the middle" who runs one exchange with each party; signatures fix that.
- **RSA** (1977). Bob publishes $n=pq$ and $e$, keeps $d=e^{-1}\bmod\varphi(n)$. Inverting $m\mapsto m^e\bmod n$ without $d$ is the **RSA problem**; it is no harder than factoring $n$, and is believed to be no easier. RSA also gives **signatures**: Bob signs a hash $h$ of a document as $s=h^d\bmod n$, and anyone checks $s^e\equiv h$.

![A flow diagram of toy RSA. Key generation uses primes 101 and 113, giving n = 11413 and phi(n) = 11200. The public key is (11413, 3) and the private key is d = 7467, since 3 times 7467 = 22401 is 1 mod 11200. Alice encrypts m = 4242 to c = 303, which crosses the open network, and Bob computes 303 to the power 7467 mod 11413 = 4242.](../../images/ch77-toy-rsa.svg)

*Toy RSA from key generation to decryption. Everything Eve needs is public except $\varphi(n)$, and computing $\varphi(n)$ is as hard as factoring $n$. With 4-digit $n$ that takes microseconds; with a 617-digit $n$ the best known classical algorithm, the number field sieve, would take longer than any realistic budget.*

Public-key arithmetic is thousands of times slower than AES, so practice is always **hybrid**: a public-key step agrees on a fresh symmetric key, and authenticated symmetric encryption carries the data. That is what happens in the handshake when your browser opens a bank's page.

### Result 3: what a quantum computer would change

Every public-key system above rests on factoring or on discrete logarithms. Both reduce to a single task: finding the **period** of a function. Here is the reduction for factoring, which is completely classical.

> [!theorem] Factoring reduces to order finding
> Let $n$ be odd and composite, and let $a$ be coprime to $n$ with multiplicative order $r$ (the least $r\ge1$ with $a^r\equiv1\pmod n$). If $r$ is even and $a^{r/2}\not\equiv-1\pmod n$, then $\gcd(a^{r/2}-1,\,n)$ is a non-trivial factor of $n$.

*Proof.* Let $x=a^{r/2}\bmod n$. Then $x^2\equiv1$, so $n$ divides $x^2-1=(x-1)(x+1)$. Now $n$ does not divide $x-1$, because $x\equiv1$ would mean $a^{r/2}\equiv1$ with $r/2<r$, contradicting the minimality of $r$. And $n$ does not divide $x+1$ by hypothesis. So $n$ divides the product without dividing either factor, which is only possible if $n$ shares a proper factor with each. Hence $1<\gcd(x-1,n)<n$. $\blacksquare$

For $n=15$ and $a=7$: the powers of 7 mod 15 are $7,4,13,1$, so $r=4$; then $x=7^2\bmod 15=4$, and $\gcd(3,15)=3$, $\gcd(5,15)=5$. One can show that when $n$ has at least two distinct odd prime factors, a random $a$ satisfies the hypotheses with probability at least $\tfrac12$; for the toy modulus $n=11413$ the fraction is about 91 %. So factoring is no harder than order finding.

Classically, order finding is as hard as factoring. But $k\mapsto a^k\bmod n$ is a *periodic sequence*, and finding a period is what Fourier transforms do ([Ch. 16](../../part-5-transforms/16-fourier-transforms.md)). Shor's 1994 algorithm prepares a quantum superposition over all $k$, computes $a^k\bmod n$ in superposition, and applies a quantum Fourier transform, whose output peaks near multiples of $1/r$; continued fractions then recover $r$. A large, error-corrected quantum computer would thus break RSA, Diffie–Hellman and the elliptic-curve systems of [Ch. 78](78-elliptic-curves.md) in polynomial time. None exists yet, but traffic recorded today can be decrypted later, so the migration has begun: in 2024 NIST standardised **post-quantum** algorithms built on finding short vectors in high-dimensional lattices ("learning with errors"), for which no efficient quantum attack is known. ML-KEM (FIPS 203) replaces Diffie–Hellman key agreement and ML-DSA (FIPS 204) replaces RSA and elliptic-curve signatures; SLH-DSA (FIPS 205) is a conservative backup whose signatures rest only on the security of a hash function. The price is size: post-quantum keys and signatures run to kilobytes instead of tens of bytes. Symmetric ciphers survive with longer keys, since Grover's quantum search only square-roots brute force.

## 77.4 The algorithm

The code below is a teaching model. It implements textbook RSA and the classical half of Shor's algorithm with numbers small enough to check by hand. Python's built-in three-argument `pow` does modular exponentiation by square-and-multiply ([Ch. 57 §57.4](../e-abstract-math/57-number-theory.md#574-the-algorithm)), and `pow(e, -1, phi)` computes a modular inverse by the extended Euclidean algorithm.

```python
from math import gcd
import random

# --- toy RSA: illustration only, NEVER for real data ---
def keygen(p, q, e=3):
    n, phi = p * q, (p - 1) * (q - 1)
    assert gcd(e, phi) == 1
    return (n, e), pow(e, -1, phi)          # public key, private exponent

def enc(pk, m):  n, e = pk; return pow(m, e, n)
def dec(pk, d, c): n, _ = pk; return pow(c, d, n)

pk, d = keygen(101, 113)
c = enc(pk, 4242)
print(pk, d, c, dec(pk, d, c))              # (11413, 3) 7467 303 4242

# --- factoring from order finding (the classical half of Shor's algorithm) ---
def order(a, n):                            # brute force: the step a quantum computer replaces
    r, x = 1, a % n
    while x != 1:
        x, r = x * a % n, r + 1
    return r

def factor_via_order(n, tries=20):
    for _ in range(tries):
        a = random.randrange(2, n - 1)
        g = gcd(a, n)
        if g > 1:
            return g                        # lucky guess
        r = order(a, n)
        if r % 2 == 0 and pow(a, r // 2, n) != n - 1:
            return gcd(pow(a, r // 2, n) - 1, n)
    return None

print(factor_via_order(15), factor_via_order(11413))   # e.g. 3 101
```

**Complexity.** Encryption and decryption cost $O(\log d)$ multiplications of $\lambda$-bit numbers, roughly $O(\lambda^3)$ bit operations: fast for $\lambda=2048$. The function `order` may loop $r$ times, and $r$ can be close to $\varphi(n)/2$, around $2^{2047}$ for a 2048-bit modulus. That one loop is the gap between "encrypt in a millisecond" and "break never". Shor's algorithm replaces exactly that loop, and nothing else changes.

**Failure cases, mathematical.** `keygen` fails when $\gcd(e,\varphi(n))>1$; with $p=7$, $q=13$ and $e=3$, $\varphi(n)=72$ is divisible by 3 and $d$ does not exist. Decryption fails if $m\ge n$, since only $m\bmod n$ survives. `factor_via_order` fails for prime powers $n=p^k$ (try $n=9$), which need a separate test.

**Failure cases, engineering.** These are the ones that break real systems, and the code above has all of them.

- *Determinism.* `enc` is deterministic, so by the theorem of §77.3 it is not IND-CPA secure. With $e=3$ and a message $m<n^{1/3}$, $m^3$ never wraps around $n$, and an ordinary integer cube root recovers $m$.
- *Malleability.* $\operatorname{Enc}(m_1)\operatorname{Enc}(m_2)=\operatorname{Enc}(m_1m_2)\bmod n$, so Eve can transform ciphertexts meaningfully without the key.
- *Side channels.* `pow` takes time, and draws current, that depends on the bits of $d$; such measurements have recovered keys from smart cards and remote servers. Production code is constant-time and blinded.
- *Bad randomness.* Devices with poor entropy at boot have produced moduli that share a prime; one $\gcd(n_1,n_2)$ then factors both.
- *Wrong error messages.* A server that reports "bad padding" differently from other errors lets Eve decrypt any ciphertext by adaptive queries (Bleichenbacher, 1998).

None of these is a mistake in the number theory. Each is a gap between a mathematical object and a running program, which is why the warning in §77.2 is not a formality.

> [!mystery] Why are elliptic-curve keys so short?
> To get security comparable to a 128-bit symmetric key, current guidance asks for an RSA modulus of about 3072 bits, but an elliptic-curve key of only 256 bits. Both rest on "easy forward, hard backward" arithmetic in a finite group. Why does one group need twelve times fewer bits than the other, and what exactly is the group made of? The answer is a curve, a chord, and the absence of a shortcut, in [Ch. 78](78-elliptic-curves.md).

## 77.5 The experiment

The script below plays Eve against a **two-time pad**. Two messages are encrypted with the same random pad, exactly the Venona mistake. Eve XORs the ciphertexts, which cancels the pad, and then **drags a crib**: she guesses that a common word appears in one message, slides it along $c_1\oplus c_2$, and at each position XORs it in. Wherever the guess is right, a fragment of the *other* message appears in plain text.

```python
import numpy as np
import matplotlib.pyplot as plt

M1 = b"meet the courier at the north gate at dawn"   # KNOB: first secret message
M2 = b"the shipment leaves the harbour on tuesday"   # KNOB: second secret message
CRIB = b" the "                                      # KNOB: a word Eve guesses is in one message
REUSE = True                                         # KNOB: False gives each message a fresh pad

rng = np.random.default_rng(2024)
L = min(len(M1), len(M2))
m1 = np.frombuffer(M1[:L], dtype=np.uint8)
m2 = np.frombuffer(M2[:L], dtype=np.uint8)
k1 = rng.integers(0, 256, L, dtype=np.uint8)
k2 = k1 if REUSE else rng.integers(0, 256, L, dtype=np.uint8)
c1, c2 = m1 ^ k1, m2 ^ k2                # what Eve intercepts
x = c1 ^ c2                              # the pad cancels if it was reused

crib = np.frombuffer(CRIB, dtype=np.uint8)
score = []
for pos in range(L - len(crib) + 1):     # crib dragging: slide the guess along
    frag = x[pos:pos + len(crib)] ^ crib
    ok = np.mean([(32 <= b < 127) and (chr(b).isalpha() or b == 32) for b in frag])
    score.append(ok)
    if ok == 1.0:
        print(f"pos {pos:2d}: other message reads {bytes(frag)!r}")

plt.bar(range(len(score)), score)
plt.xlabel("crib position"); plt.ylabel("fraction of letters/spaces")
plt.title("crib dragging on c1 XOR c2"); plt.show()
```

With the defaults, six positions produce readable-looking fragments. One of them, at position 4, reads `'shipm'`: a genuine piece of the second message, exposed because the first message has " the " there. Another, at position 19, reads `' the '` because *both* messages happen to contain " the " at the same place. The rest are look-alikes such as `'mridr'` that a human rejects at a glance. Things to try:

- **Extend a hit.** Guess that `shipm` continues as `shipment`, set `CRIB = b"shipment"` and look at position 4 again: you now read eight characters of the first message. Real codebreakers bootstrapped whole messages this way, one guessed word at a time.
- **Set `REUSE = False`.** Every position now gives noise, and no fragment passes the filter. Result 1 is at work: with independent pads, $c_1\oplus c_2$ is itself a uniformly random string, and no guess can be confirmed.
- **Change the language model.** The filter accepts any lowercase letters, which is crude. Score by English letter frequencies instead (English carries only about 1 to 1.5 bits per letter out of a possible $\log_2 26\approx4.7$, [Ch. 33](../b-probability-information/33-information-theory.md)). How short a crib still finds the true position?

## 77.6 The challenge

No answers are given.

1. **Intuition.** A friend proposes "double encryption": encrypt with a one-time pad $k_1$ and then again with an independent pad $k_2$, for twice the security. Explain in words why this gains nothing. Then explain why encrypting twice with a *block cipher* under two keys (as in "double DES") also gains far less than doubling the key length suggests. (Think of meeting in the middle.)
2. **Derivation.** Prove the converse half of Shannon's theory: if $\lvert\mathcal{K}\rvert=\lvert\mathcal{M}\rvert=\lvert\mathcal{C}\rvert$ and the scheme is perfectly secret, then every key is used with probability $1/\lvert\mathcal{K}\rvert$, and for every pair $(m,c)$ there is exactly one key with $\operatorname{Enc}_k(m)=c$. Then show that the definition of perfect secrecy is equivalent to: for all $m_0,m_1,c$, $P(C=c\mid M=m_0)=P(C=c\mid M=m_1)$.
3. **Implementation.** Extend the toy RSA code with the Chinese-remainder speed-up for decryption (compute $c^d$ mod $p$ and mod $q$ separately and combine). Measure the speed-up for 512-bit primes generated with a Miller–Rabin test. Then implement a toy version of the common-factor attack: generate 1000 moduli where, through a deliberately weak random generator, a few share a prime, and find them all with pairwise gcds. Keep everything in a sandbox: this is to understand the failure, not to build anything.
4. **Open investigation.** In the experiment, replace the letter filter by a statistical score and automate crib dragging over a dictionary of the 200 most common English words. For two messages of length $L$, how does the fraction of plaintext you can recover automatically grow with $L$? Is there a threshold length below which the two-time pad is "practically safe", and how does it depend on the language (try a language with a different letter distribution, or compressed text)?

## 77.7 Loose threads

The hook asked three questions. **Security** means winning a precisely defined game against every allowed adversary; perfect secrecy and IND-CPA are two such games, and the second showed in three lines that deterministic ciphers fail. **Short keys** protect long messages only against a bounded Eve, through pseudorandom generators and block ciphers resting on unproved assumptions; Shannon's counting shows nothing else is possible. **Strangers agree on a secret** through trapdoor arithmetic whose hardness rests on factoring and discrete logarithms, which a quantum computer would turn into period finding.

The Venona lesson runs under all of it: mathematics proves a scheme secure *under stated conditions*, and practice fails by quietly violating one.

What we did not tell you:

- **Where good randomness comes from**, and what "looks random" must mean so that a stream cipher inherits the one-time pad's security. The first mystery points to [Ch. 81](81-pseudorandomness.md).
- **What the elliptic-curve group is.** The second mystery points to [Ch. 78](78-elliptic-curves.md), and the finite fields underneath are introduced in [Ch. 57](../e-abstract-math/57-number-theory.md).
- **Hash functions.** Signatures, passwords and blockchains all depend on functions where finding two inputs with the same output is infeasible. The birthday bound of [Ch. 32](../b-probability-information/32-counting-random-worlds.md#the-birthday-bound-made-general) is why a 256-bit hash gives only 128 bits of collision resistance.
- **Why lattice problems resist quantum computers.** Nobody knows. Ajtai (1996) and Regev (2005) showed that random lattice instances are as hard as the worst ones, a guarantee factoring never had, but a quantum algorithm for short lattice vectors has not been ruled out.

> [!mystery] Proving you know a secret without revealing it
> A signature convinces Bob that Alice holds a private key, but the signature itself is a lasting piece of evidence that Bob can show to anyone. Is it possible to convince Bob that you know the factorisation of $n$, or the password, or the solution to a puzzle, so that he learns *nothing* except that the claim is true, not even enough to convince a third party? It is, and the definition of "learns nothing" is one of the strangest ideas in computer science: anything Bob could see, he could have simulated alone. See [Ch. 82](82-zero-knowledge-proofs.md).

## Further reading

- J. Katz and Y. Lindell, *Introduction to Modern Cryptography*, 3rd ed., CRC Press, 2020. The standard text for security definitions, reductions and proofs.
- D. Boneh and V. Shoup, *A Graduate Course in Applied Cryptography*, freely available online. Rigorous, current, and full of real attacks.
- C. E. Shannon, "Communication theory of secrecy systems", *Bell System Technical Journal* 28 (1949), 656–715. Where perfect secrecy and the key-length bound began.
- N. Ferguson, B. Schneier and T. Kohno, *Cryptography Engineering*, Wiley, 2010. The engineering side: what goes wrong between a correct algorithm and a secure system.
