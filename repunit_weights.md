# Prime-length repunit weights

Write $R_n(X)=1+X+\cdots+X^{n-1}$. Let $n\ge3$ be prime and $Q=p^a$ with $p$ prime and $a\ge1$. The following results concern **every choice of three positive trinomial rows** at $N=R_n(Q)$.

* The conjecture holds for all $Q$ whenever $n$ is one of the prime lengths in the finite table below.
* For every prime length it holds whenever $Q>1+n(n+1)3^{3(n-1)/4}$. In particular the range of lengths covered at once grows at least proportionally to $\log Q$.
* If $n-1=dr$ with $r$ prime and $d\le6$, a possible exception requires $Q=O_d(n^{1+3d/4})$. For prime lengths $n=2r+1$ with $r$ prime this gives $Q=O(n^{5/2})$. These statements do not require infinitely many primes of any specified shape.

## What a single digit subset must satisfy

By the [endpoint reduction](all_weights_endpoint_reduction.md), a counterexample contains $(1,1,N-2)$, whose coefficient is $N(N-1)$. Distribute the full powers of $N$ and $N-1$ between its two other rows according to which row avoids the underlying prime. A row avoids a prime when its coefficient is coprime to it. Choose the first row to avoid $p$, which divides $N-1$. Write $UV=N$ and $ST=N-1$ for the resulting assignments, with $U,S$ assigned to the first row. Its parts are multiples of $U$. The other row has parts $Vy_1,Vy_2,Vy_3$ with positive $y_i$ summing to $U$, so $U\ge3$. At each full power in $T$, at least two $y_i$ are zero modulo that power. Thus $27T^2\le U^3$. At each full power in $S$ the residues of the first row are $0,0,1$. These facts follow from carrying at the given prime power. [Assignment proof](endpoint_factors.md)

The base-$p$ digits of $N$ consist of ones at positions $0,a,\ldots,(n-1)a$. By Kummer's criterion the first row partitions the $n$ powers of $Q$ into three nonempty subsets. Take a subset of smallest cardinality $k\le n/3$, with polynomial $F(X)=\sum_{i\in B}X^i$. It cannot be a singleton since $U$ is coprime to $Q$. With $m=n-1$, the necessary divisibilities are
$$
 U\mid R_n(Q),\quad U\mid F(Q),\qquad R_m(Q)\mid T F(Q)(F(Q)-1). \tag{1}
$$
The last one follows from $S\mid F(Q)(F(Q)-1)$ and $N-1=Q R_m(Q)$. The same divisibilities hold for each of the three parts.

Every prime $s\mid U$ has $\operatorname{ord}_s(Q)=n$. Indeed $Q^n\equiv1\pmod s$. If $Q\equiv1\pmod s$, the two divisibilities force $s\mid n$ and $s\mid k$, impossible. In particular $s\equiv1\pmod n$, and no part can have exactly two terms either. Such a part would give $Q^h\equiv-1\pmod s$ for $0<h<n$, contradicting the odd order $n$.

## The two resultant bounds

For any integer polynomial $G$ with $U\mid G(Q)$ the integer resultant $\operatorname{Res}(R_n,G)$ is divisible by $U$, by the integral Bézout identity for resultants. For a nonempty proper subset polynomial $G$ with $v$ terms, this resultant is nonzero because $R_n=\Phi_n$ is irreducible. At the $n$th roots of unity Parseval's identity gives
$$
 \sum_{\zeta^n=1,\,\zeta\ne1}|G(\zeta)|^2=nv-v^2.
$$
The arithmetic-geometric mean inequality gives
$$
 U\le|\operatorname{Res}(R_n,G)|\le
 \left(\frac{v(n-v)}{n-1}\right)^{(n-1)/2}\le v^{(n-1)/2}. \tag{2}
$$

At the nontrivial $m$th roots call a root bad for $F$ when $F(\zeta)\notin\{0,1\}$. There is at least one. After folding the coefficients of $F$ modulo $m$, the sum of their squares is at least $k$. Parseval then gives at least $mk-k^2>m-1$ for the sum of the squares of the absolute values over nontrivial roots, since $(k-1)(m-k-1)>0$. If they were all good the sum would be at most $m-1$.

Take any nonempty union of bad Galois orbits, of total degree $e$, and let $H$ be the product of their cyclotomic polynomials. As $H\mid R_m$, (1) gives
$$
 (Q-1)^e\le |H(Q)|\le T\,|\operatorname{Res}(H,F(F-1))|\le T c^e,
 \qquad
 c=\min\left(k(k+1),\frac{m(k+3)}e\right). \tag{3}
$$
For the middle inequality, $H(Q)/\gcd(H(Q),F(Q)(F(Q)-1))$ divides $T$ and the gcd divides the nonzero resultant. The triangle inequality gives the first choice of $c$. For the second, only the coefficients at positions 0 and $m$ can collide upon folding. The folded sums of squares for $F$ and $F-1$ are at most $k+2$ and $k+3$. Parseval, followed by arithmetic-geometric mean applied separately to $|F|^2$ and $|F-1|^2$ on the chosen roots, gives the bound.

## A uniform exponential bound on the base

Use all bad roots of a smallest part and let $D$ be the lcm of their orders. Thus $D\mid m$. All orders $s\mid m$ with $s\nmid D$ have the same good value, either 0 or 1. To see this, follow a chain from $s$ to $m$ by multiplying by one prime $l$ at a time. Each order still fails to divide $D$. The polynomials $\Phi_s$ and $\Phi_{sl}$ have a common root in characteristic $l$. Indeed if $s=l^a u$ with $l\nmid u$, reducing $\Phi_{l^a u}$ modulo $l$ gives $\Phi_u^{\varphi(l^a)}$ for $a\ge1$, directly from the cyclotomic factorization of $X^{ul^a}-1$. If the good values differed, evaluating $F$ at the common root would equate 0 and 1.

Choose $G=F$ for common value 0, and $G=R_n-F$ for common value 1, since $R_n(\zeta)=1$ at all nontrivial $m$th roots. For $D=m$ take $G=F$. We have
$$
 G=L A,\qquad L=\frac{X^m-1}{X^D-1}=1+X^D+\cdots+X^{m-D}.
$$
The quotient $A$ is integral of degree at most $D$. Its coefficients in degrees below $D$ are those of $G$, and its coefficient in degree $D$ is that of $G$ in degree $m$. Hence $A$ has only zero and one coefficients. Also $L(Q)\mid R_m(Q)$ and $N=1+Q R_m(Q)$, so $U\mid A(Q)$. The nonempty support of $A$ is proper in $\{0,\ldots,n-1\}$, as $D<m$ or $A=G$. Applying (2) gives $U\le(D+1)^{(n-1)/2}$.

For each bad order $s$, we have $s+1\le3^{\varphi(s)}$. In fact $\varphi(s)^2/s\ge1/2$ by its prime factorization, and $2b^2+1\le3^b$ for integers $b\ge1$ by induction. Multiplying over these orders yields $D+1\le3^e$. Now (3), $T\le U^{3/2}$ and $c\le n(n+1)$ give
$$ Q\le1+n(n+1)3^{3(n-1)/4}.$$

## A polynomial bound when the predecessor of the length has a large prime factor

Suppose $m=dr$ with $r$ prime and $d\le6$. Some part has a bad primitive root of an order $s\mid m$ not dividing $d$. Otherwise apply the chain argument just given to each part. Each has a constant value at all these orders, and the three values sum to 1. The two parts with value 0 equal $L A_i$ with $L=1+X^d+\cdots+X^{(r-1)d}$ and $A_i$ a zero-one polynomial of degree at most $d$. Since $U\mid A_i(Q)$, neither $A_i$ can be a monomial or binomial, by the order-$n$ condition. Each of the two corresponding parts therefore contains at least $3r$ digits. The remaining part has a number of digits congruent to 1 modulo $r$. It cannot be a singleton, so it has at least $r+1$ digits. This would give $n\ge7r+1$, contrary to $n=dr+1$.

Here $r\mid s$, so that orbit has $e=\varphi(s)\ge r-1$. Apply (3) to the part with this orbit, using its size at most $n$, and (2) to the smallest part. The explicit bound for an exception is
$$
 Q\le1+\frac{(n-1)(n+3)}{r-1}
             (n/3)^{3(n-1)/(4(r-1))}.
$$
For fixed $d$ this is $O_d(n^{1+3d/4})$.

## Finite certificates

For each smallest subset shape the verifier computes its algebraic norm as the absolute determinant of multiplication by $F$ in $\mathbb Z[X]/R_n$. It factors the integer and checks the factorization. At a candidate prime $s\mid U$ only roots of order $n$ modulo $s^j$ are possible. Starting with a primitive $n$th root $z$ modulo $s$, the powers of $z^{s^{j-1}}$ give all such roots modulo $s^j$, by lifting or the binomial theorem. The exponent $j$ need not exceed its multiplicity in the norm. The largest possible $s^j$ for each $s$ gives a product upper bound on $U$. For $s^j>1000$ the two unused parts must also split the remaining powers into two zero sums of at least $k$ terms modulo $s^j$. Testing this is optional for the proof and reduces the enumeration.

Affine maps on the positions modulo $n$ preserve the possible prime powers and transform their root classes explicitly. The verifier checks that the orbit sizes total $\binom nk$. For every actual subset it obtains a lower bound on $e$ in (3) by testing primitive roots in a finite field of characteristic $1\bmod m$. A value outside $\{0,1\}$ proves that the entire complex orbit is bad. When such a modular test is inconclusive there is a second finite bound $Q\le1+(k+4)T$. To see it, reduce $F^2-F$ modulo $X^m-1$ to cyclic coefficients $b_j$. They are between $-2$ and $k+2$ by Cauchy-Schwarz for the convolution. They cannot all be equal since there is a bad root. From (1), $Q^m-1$ divides $T\sum_j (b_{j-1}-b_j)Q^j$ with indices taken cyclically. If $Q-1>(k+4)T$ this integer has absolute value less than $Q^m-1$ and must vanish. Successive reduction modulo $Q$ forces all differences to vanish, a contradiction.

For every remaining integer base and subset, put $f=F(Q)$, $u=\gcd(N,f)$ and $t=(N-1)/\gcd(N-1,f(f-1))$. A counterexample needs $u\ge3$ and $27t^2\le u^3$. Where a subset passes, the program tests the full prime-power assignments by the Chinese remainder theorem and Kummer's criterion. In this last step the row labels may be chosen so that $S\ge V$. Indeed $S<V$ and $T<U$ together would give $N-1=ST\le(V-1)(U-1)<N-1$. Each of the three positive normalized parts, being less than $V$, then needs an exceptional residue for at least one power in $S$. [CRT row algorithm](search_endpoint.py)

| Prime length | Subsets of sizes $2$ through $\lfloor n/3\rfloor$ | Subsets with a possible weight prime | Largest base cutoff | Base and subset pairs tested |
|---:|---:|---:|---:|---:|
| 3 | 0 | 0 | 0 | 0 |
| 5 | 0 | 0 | 0 | 0 |
| 7 | 21 | 0 | 0 | 0 |
| 11 | 220 | 110 | 26 | 71 |
| 13 | 1,079 | 624 | 84 | 130 |
| 17 | 9,384 | 7,616 | 268 | 469 |
| 19 | 43,776 | 30,438 | 131 | 967 |
| 23 | 390,632 | 348,887 | 178 | 107,563 |
| 29 | 16,489,516 | 14,124,740 | 6,107 | 4,746,792 |
| 31 | 75,973,157 | 52,713,795 | 3,771 | 2,502,482 |

Only $n=29$, $Q=2$ required the final CRT step. There the weight is $536870911$ and the row receiving $S\ge V$ already fails in every assignment.

The independent finite-base verifier enumerates **every prime power** up to the indicated cutoff. It forms the subset sums of two halves and looks for pairs with total divisible by a prime of the weight. At lengths through 23 it factors the weight directly. At lengths 29 and 31 it computes a catalog of all possible weight primes independently. It enumerates binary necklaces under rotation, groups them by exponent multiplication, and checks coverage against the binomial counts. For each representative it multiplies its evaluations at all nontrivial $n$th roots in one finite field whose characteristic exceeds (2). This recovers the exact positive integer norm, since the complex roots occur in conjugate pairs. All prime factors of that norm congruent to 1 modulo $n$ are included. Where a full predecessor power $w$ satisfies $2w+1>U_{\max}$, it must belong to $S$, since otherwise two of the complementary row's normalized parts would be at least $w$. The verifier can use $f\equiv0,1\pmod w$ to accelerate the subset search. Its CRT fallback reconstructs each part directly and uses factorial floors.

This second calculation covers 840 prime-power bases at length 29 and 564 at length 31. It finds the same single base requiring a full row check. Both programs use Python with python-flint. The primary verifier's output also gives the largest remaining weight prime and the number of distinct CRT fallbacks. Use the base cutoffs in the table for the independent checks.

~~~sh
python3 check_repunit_weights.py
python3 check_repunit_weights_cross.py 11 26
python3 check_repunit_weights_cross.py 13 84
python3 check_repunit_weights_cross.py 17 268
python3 check_repunit_weights_cross.py 19 131
python3 check_repunit_weights_cross.py 23 178
python3 check_repunit_weights_cross.py 29 6107 --catalog
python3 check_repunit_weights_cross.py 31 3771 --catalog
~~~
