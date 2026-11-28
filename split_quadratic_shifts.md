# Split quadratic blocks at unrestricted outer exponents

Fix a prime $p$, positive exponents $a,e$, $X=p^a$, $Q=p^e$, and the weight $N=1+QG(X)$. For **120 of the 121** polynomials $G=X^2+bX+c$ with $0\le b,c\le10$, every three positive trinomial rows of weight $N$ have a common prime. For the last polynomial the result holds except possibly when
$$
G=(X+2)(X+4),\qquad p=3,\qquad e+1=an,\qquad n\text{ odd}.
$$
Even there any failure must include both the endpoint and the row
$$
(3Q(X+1),\ Q(X+1)(X+2),\ 3Q+1).
$$
The factor certificate below excludes every outer exponent for $1\le a\le300$. More generally, for **every fixed monic quadratic** with nonnegative integer coefficients and every $K>0$, all $e\le K a^2$ work once $p^a$ is sufficiently large, uniformly in $p$.

There is an unrestricted-exponent theorem for larger coefficients. If $d=p^m$ with $m\ge0$, $k\ge1$, $k\ne d$, $G=(X+d)(X+k)$, and $X>(3+d+k+dk)^{60}$, every possible failure has the form (3) below. In particular **$G=(X+1)(X+k)$ works for every outer exponent whenever $2\le k<2p+2$ and $X>(4+2k)^{60}$.**

## The two rows outside the endpoint

Write $H=3+d+k+dk$. Assign each prime of $N(N-1)$ to a nonendpoint row avoiding it, with $p$ assigned to $A$. Let $U$ be the product of the full prime powers assigned to $A$ in $N$, and $QB$ their counterpart in $N-1$. Thus $G(X)=BT$. The companion row is $(N/U)(y_1,y_2,y_3)$ with positive $y_i$ summing to $U$. Each full prime power of $T$ divides two $y_i$. Their product and AM–GM give
$$
U\ge3,\qquad 27G(X)^2\le U^3 B^2.                                      \tag{1}
$$
When $X>\max(d+k,dk)$, the parts of $A$ are $Qx(X),Qy(X),1+Qz(X)$ for three nonzero polynomials of nonnegative integral coefficients summing to $G$. Kummer's criterion requires the splitting within each coefficient block to be carry-free in base $p$. We have $U\mid x(X),y(X)$ and $\gcd(U,pG(X))=1$. At any full prime power of $QB$, the parts of $A$ have residues $0,0,1$ in some order. In particular each full power of $B$ divides at least one unshifted value $x(X),y(X)$, including at $p$ by cancellation of $Q$. [Endpoint reduction](all_weights_endpoint_reduction.md) · [Normalization and shifted congruences](shifted_monic_blocks.md)

Write $B=B_d B_k$, where $B_i$ is the part of $X+i$ supported on primes of $B$. This accounts for multiplicities even when the evaluated factors share primes. If $x(-i)y(-i)\ne0$, then $B_i$ divides $\operatorname{lcm}(|x(-i)|,|y(-i)|)\le H^6$.

## When one factor has a prime-power constant term

Assume $d=p^m$, $X>H^{60}$. At least one of $-d,-k$ is a zero of neither unshifted polynomial. Otherwise distinct unshifted polynomials contain the two roots, since a quadratic containing both would exhaust $G$. If a low-degree part is $\alpha(X+k)$, carry freedom at the constant block $dk$ forces $d\mid\alpha$, and positivity gives $\alpha\le d$. Then $U\mid d$ is impossible. Otherwise the monic quadratic containing $X+k$ is $(X+k)(X+j)$ with $0\le j<d$. The constant block forces $d\mid j$, hence $j=0$. The other unshifted part is $\alpha(X+d)$, and it must take the whole linear coefficient of the remainder $d(X+k)$ since $d$ consists of one digit unit. Again $U\mid d$.

Consequently $B\le2XH^6$. Equation (1) excludes $U\le H^3$. The lower-degree unshifted part is therefore $x=\alpha X+\beta$ with both coefficients positive. Put $r=\beta/\alpha$. The other unshifted polynomial satisfies $y(-r)=0$, since otherwise $U$ divides the nonzero integer $\alpha^{\deg y}y(-r)$ of magnitude at most $H^3$. The root parameter $r$ differs from $d,k$ by coprimality of $U$ with $G(X)$. If neither factor of $G$ divides $y$, then $B\le H^{12}$ and $U\le HX$ contradict (1). Thus $y$ is monic quadratic, with roots $-r$ and one of $-d,-k$. In particular $r$ is integral. Taking the second root to be $-k$ would give $0<r<d$ and $d\mid r$ from the constant block. Hence
$$
x=\alpha(X+r),\quad y=(X+r)(X+d),\quad z=(D-\alpha)X+dD-\alpha r,\qquad D=k-r>0.
                                                                    \tag{2}
$$
Here $\alpha\le D$, $\alpha r\le dD$, and $d\mid\alpha r$. We have $U\mid X+r$ since $U\mid y(X)$ and $\gcd(U,X+d)=1$. Also $U\mid Q(d-r)D+1$. If $r=d+\alpha$, the reduced sum of the two unshifted parts is $(X+r)/\gcd(\alpha,X+d)$, a multiple of $U$. For positive $v,w$, $(v+w)/\gcd(v,w)$ divides $\binom{v+w}v$ by the two adjacent binomial identities. This contradicts avoidance of every prime of $U$ by $A$.

Put $C=\alpha(r-d)$ and $L=D\operatorname{lcm}(\alpha,|k-d|)$. At $-k$ the evaluations give $B_k\mid L$. At $-d$ the product of assigned powers whose exceptional residue is outside slot $x$ divides $|C|$. The powers switching to slot $x$ divide $CQ-1$. With $A_0=a-m>0$ and $g=\gcd(p^{A_0}+1,Cp^e-1)$ we get $B_d\le d|C|g$. Using $U\le X+r$ in (1) yields $27X<2L^2 d^2 C^2g^2$.

We use a gcd estimate for any signed nonzero integer $C$. If $|C|$ is not a power of $p$, choose $j=1$ or 2 and $\ell\ge0$ with $|je-\ell A_0|\le A_0/3$ by comparing the three circular gaps between $0,e,2e$ modulo $A_0$. Put $v=je-\ell A_0$. Then $g$ divides $(-1)^\ell C^j p^{\max(v,0)}-p^{\max(-v,0)}$. This is nonzero, so $g\le2|C|^2p^{A_0/3}$. For $C=\sigma p^\nu$, $\sigma=\pm1$, the full gcd occurs precisely when $e+\nu=n A_0$ and $(-1)^n=\sigma$. In every other case $g\le p^{A_0/3}+1$. Indeed if $g>2$, Bézout at $b=\gcd(A_0,e+\nu)$ gives $p^b\equiv\pm1\pmod g$. It must be $-1$ with $A_0/b$ odd. If $b=A_0$ the signs force the full gcd, and otherwise $A_0/b\ge3$. These arguments also prove the asserted criterion for the full gcd.

The nonfull bounds in the preceding paragraph, together with $L\le H^3$, $d\le H$ and $|C|\le H^2$, give $X\le H^{60}$. Therefore
$$
 r=d+\sigma p^s>0,\qquad \alpha=p^t,\qquad r\ne d+\alpha,\qquad
 e+s+t=n(a-m),\qquad (-1)^n=\sigma,\quad n\ge1,\quad s,t\ge0.       \tag{3}
$$
Neither $D$ nor $r$ can be a power of $p$ (including 1). If $D$ is such a power, the remainder $D(X+d)$ has a single digit unit in each coefficient. The positive line $x$ takes both, leaving $z=0$. If $r$ is such a power and $D$ is not, apply the nonpower gcd estimate to $X/r+1$ and $-Q(d-r)D-1$. It bounds $U\le2|(d-r)D|^2 X^{1/3}$. But $27(X+k)^2\le L^2U^3$ by (1) and $B_d\le X+d$, giving $X<H^{18}$. This proves the claimed restrictions. For $d=1$, positivity forces $r=1+p^s$, and $s\ne t$. Whether $s\ge1$ or $t\ge1$, the inequality $k\ge r(1+\alpha)$ gives $k\ge2p+2$.

The digit conditions sharpen the $d=1$ range to $k<27$ at $p=2$ and $k<20$ at $p=3$. At 2, $s=0$ gives a forbidden power $r=2$. For $s>0$, $t=0$ repeats the same nonzero constant block twice. If $s=1$, the linear coefficient of $y$ is 4, so $t=2$ is impossible as well as $t=s$. This gives $k\ge3(1+8)$. If $s\ge2$, its linear coefficient $2+2^s$ forbids $t=1$, so distinct $s,t\ge2$ give $k\ge45$. At 3, $s\ge1,t=0$ creates a units carry in the linear block. Distinct positive $s,t$ give $k\ge40$. For $s=0$ we have $r=2$, $t\ge1$. If $t=1$, the two constant contributions $2,6$ require the bottom two ternary digits of $k$ to be 22. Then $k+1$ has zero tens digit, incompatible with the two linear contributions $3,3$. The remaining $t\ge2$ gives $k\ge20$. At the endpoints $k=27,20$ respectively, $(s,t)=(1,3),(0,2)$ pass the coefficient digit tests. For $p\ge5$, $k=2p+2$, $(s,t)=(1,0)$ passes them. These examples establish sharpness for this digit test alone.

For each fixed shape in (3), substitution at $X\equiv-r\pmod U$ gives
$$
U\mid D r^n-\alpha d^n,\qquad U\mid X+r,\qquad U^3\ge27(X+k)^2/L^2. \tag{4}
$$
For two possible indices $n_1<n_2$ at the same $a$, the common divisor of their two $U$'s is $\gg_{d,k} X^{1/3}$ since their lcm divides $X+r$. By (4) that divisor divides the nonzero integer $r^{n_2-n_1}-d^{n_2-n_1}$, using coprimality with $p,r,D$. Thus $n_2-n_1\gg_{p,d,k}a$ at large $a$. Summing over the finitely many shapes, every interval of length $E$ contains $O_{p,d,k}(1+E/a^2)$ possible outer exponents.

## Every fixed quadratic throughout any quadratic range of exponents

**Theorem.** For any fixed nonnegative integers $b,c$ and any real $K>0$, the conjecture holds at every $1+p^eG(p^a)$ with $G=X^2+bX+c$, $e\le Ka^2$, and $p^a$ sufficiently large in terms of $b,c,K$. The threshold is uniform in the prime and need not be effective.

The [explicit quadratic bound](shifted_monic_blocks.md) settles all primes $p\ge H^{\max(212,24K)}$ for $H=b+c+3$. It suffices to consider one fixed smaller prime at a time. The same bound settles every sequence with bounded $e/a$ and $a\to\infty$. An infinite sequence of failures at this $p$ would have a subsequence with $a\to\infty$, $e/a\to\infty$ (in fact $e\gg_{p,H}a^2$). The irreducible and square cases already hold at all large inputs. If $G=X(X+k)$, $k>0$, the [successor bound](successor_powers.md) requires $27(X+k)^2<k^5$ in a failure. Thus assume $G=(X+h)(X+k)$ with positive distinct $h,k$. We may fix the coefficient split along a subsequence.

**Gcd lemma.** Along such a sequence, let $j$ be a fixed positive integer, $c_0$ a fixed nonzero integer. If $g=\gcd(p^a+j,c_0p^e-1)>p^{\epsilon a}$ for some fixed $\epsilon>0$ infinitely often, both $j$ and $|c_0|$ are powers of $p$, including 1.

Write $e=\ell a+v$ with $0\le v<a$. The rational $S$-units
$$u=-p^a/j,\qquad v'=c_0 p^v(-j)^\ell$$
are both congruent to 1 modulo $g$ since $\gcd(g,j)=1$. Take $S$ to contain infinity and the prime divisors of $pjc_0$. Both logarithmic heights are $O_{p,j,c_0,K}(a)$ because $\ell\le Ka$. The sum over places of $\log^-\max(|u-1|_w,|v'-1|_w)$ is at most $-\log g$. Proposition 2 of Corvaja and Zannier therefore gives, apart from finitely many pairs, finitely many possible relations $u^i=(v')^t$ with $(i,t)\ne(0,0)$. Along a subsequence fix the relation. Since $|u|\to\infty$, $t\ne0$. Any prime $q\ne p$ dividing $j$ would give $-i v_q(j)=t(v_q(c_0)+\ell v_q(j))$, impossible as $\ell\to\infty$. Now $j$ is a power of $p$, and any $q\ne p$ dividing $c_0$ would give $0=t v_q(c_0)$. This proves the lemma, including when $v'$ is a root of unity. [Corvaja and Zannier, Proposition 2](https://arxiv.org/pdf/math/0311030)

If $U\le H^3$ along a subsequence, (1) gives $B\gg_H X^2$, hence $B_h,B_k\gg_H X$. For each $j=h,k$, if both unshifted evaluations at $-j$ are nonzero then $B_j$ is bounded. If just one is nonzero, call it $c_0$. The exceptional slot must switch to that part except for a group dividing $|c_0|$, so $B_j\le |c_0|\gcd(X+j,c_0 Q-1)$. The lemma forces $j$ to be a power of $p$. If this happens at both roots, $G(X)$ is a sum of four powers of $p$ counting multiplicity, so has at most four digit units. Among three nonempty polynomial parts there are at least two singleton parts, one unshifted. This would force $U\mid p^t$ for some $t$. Otherwise all parts vanish at one root $-j'$, and the other parameter $j$ is a power of $p$ (both roots cannot divide a nonzero linear part). The two lower-degree parts must be positive constant multiples $v_1(X+j'),v_2(X+j')$. Their multipliers sum to at most $j$ by dividing the polynomial identity by $X+j'$ and evaluating at zero. Carry freedom in the constant block $jj'$ implies $j\mid v_1,v_2$, another impossibility.

Finally, for $U>H^3$ the common-root argument preceding (2) works for arbitrary positive $h,k$. It gives (2) with $d$ one of the two root parameters and $k$ there denoting the other. Indeed a common root $-r$ distinct from both factor roots is forced by the linear part. If neither factor divides the other unshifted polynomial, $B\le H^{12}$ contradicts (1) at large $X$. Hence that part is monic quadratic with roots $-r,-d$, making $r$ a positive integer. Here $d$ is an arbitrary parameter. Still $U\mid X+r,\ U\mid -(d-r)DQ-1$, and $U\gg_H X^{2/3}$ by (1) and $B\le L(X+d)$. Also $B_d\le |C|\gcd(X+d,CQ-1)$ gives $\gcd(X+d,CQ-1)\gg_H X^{1/2}$, where $C=\alpha(r-d)$. Apply the gcd lemma at $j=r$ and $j=d$. We obtain that $r,d$ and $|(d-r)D|$ are powers of $p$, hence so is $D$. As in (2), $x$ would have to take both digit units of $D(X+d)$, leaving $z=0$. This excludes the last case and proves the theorem.

## Exhaustion of the coefficient box

The [irreducible and square cases](shifted_monic_blocks.md) cover 100 of the quadratics. The ten other cases with $c=0$ reduce to $1+p^{e+a}(p^a+b)$ and follow from the [successor theorem](successor_exponents.md). For the eleven distinct positive factor pairs $h<k$, label the shifted part and enumerate all nonzero coefficient splits of $(X+h)(X+k)$, identifying the two unshifted slots under exchange. There are 19,330 shapes. When $X\le\max(h+k,hk)$ the programs enumerate the actual digit allocations directly.

Here $h\le2$. If both factor roots appear in the unshifted parts, they appear in different parts. Two lines $\alpha(X+h),\beta(X+k)$ have $U\mid\alpha,\beta$, so their constant sum is at least $U(h+k)$ and at most $hk$. For a line $\alpha(X+j)$ and a quadratic $(X+j')(X+l)$, where $\{j,j'\}=\{h,k\}$, subtraction gives $\alpha\le j-l$ and $\alpha j\le(j-l)j'$. Hence $U\le\alpha\le\min(j,j')$. Both alternatives contradict $U\ge3$.

When a root $-j$ belongs to neither part, write $m_j=\operatorname{lcm}(|x(-j)|,|y(-j)|)$. If $U$ is bounded by an integer $u$ from a constant part or a nonzero resultant, then $X+j=bt$ where $b\mid m_j$ and $27t^2\le u^3$. Factors $X,X+h,X+k$ can be removed from divisibility tests for $U$ by coprimality. When both $B_h,B_k$ are bounded we enumerate $B\mid m_hm_k$ and the arithmetic progressions $G(X)\equiv0\pmod B$. Inequality (1), using the lower-degree part to bound $U$, gives an exact polynomial cutoff in $X$.

The remaining shapes have (2), now with $d$ either root parameter and $k$ in (2) standing for the other. As before $U\mid X+r$, $B_k\mid L$, $B_d\le |\alpha(r-d)|\gcd(X+d,\alpha(r-d)Q-1)$. After the reduced-pair-sum test, each of the 24 remaining shapes has $r=1$ or $d=1$. In the former case take $c_0=-(d-r)D$, $g=\gcd(X+1,c_0Q-1)$ and $U\le g$, giving $27X^2\le L^2g^3$. Otherwise take $c_0=\alpha(r-d)$, $g=\gcd(X+1,c_0Q-1)$ and $27X<2L^2 c_0^2g^2$. The enumeration checks $1<|c_0|\le4$ and $L\le256$. The gcd estimate gives $g\le2|c_0|^2 X^{1/3}$ outside a full resonance at $|c_0|=p^\nu$. Hence $X<2^{81}$. Also $g\le |c_0|^a+1$ by raising $c_0Q\equiv1$ to power $a$. These two bounds make it practical to enumerate the possible $p,a$ and the multiplicative orbit of $Q$ at $X+1$.

For each remaining numerical choice the first implementation enumerates full prime-power assignments, uses the normalized-denominator restrictions, and reconstructs all possible companion rows by CRT. At a forced prime, a carry or a divisor of a reduced pair sum eliminates the row. The second uses a determinant for the two unshifted parts, solves linear congruences separately at each full power, and constructs the companion by searching allowed coordinates. Both use the previously proved [normalized-denominator exclusions](shifted_directions.md). The finite searches leave only the full resonance displayed at the start, with $r=1,d=2,k=4,\alpha=3,c_0=-3$.

## The remaining ternary shape

In this section put $X=3^a$, $s=X+1$, $F=X+2$ and $M=X+4$. At $-4$ the unshifted evaluations are $-9,6$. As $M$ is coprime to 6, all its primes go to the companion. At $-2$ the linear part evaluates to $-3$, a unit at every prime of $F$. Thus all powers assigned from $F$ switch to the linear part, requiring $3Q\equiv-1\pmod B$. Since $3Q=X^n$ and $n$ is odd,
$$
 U\mid s,\qquad B\mid\gcd(F,2^n-1),\qquad B>\sqrt{27F}.                \tag{5}
$$
The last inequality follows from (1) with $U\le s$. If $n_1<n_2<n_3$ survive at fixed $a$, each pairwise gcd of their $B_i$ divides $2^{n_j-n_i}-1$. All $B_i$ divide $F$, so prime by prime their product is at most $F\prod_{i<j}\gcd(B_i,B_j)$. By (5), $F^{3/2}<F 2^{2(n_3-n_1)}$. Thus any interval of outer exponents of length $\tfrac a4\log_2(3^a+2)$ contains **at most two** possible exceptions. For $e\le Ka^2$, applying Proposition 2 cited above to $(-3^a/2,2^n)$ gives finiteness, since the common numerator divisor exceeds $\sqrt F$, the heights are $O_K(a)$, and the two numbers are multiplicatively independent by their 3-adic valuations.

There are two exponent-independent carry tests. Throughout, $q^v\parallel h$ means $q^v$ is the full power of the prime $q$ in $h$.

1. Let $w=q^v\parallel s$, $t=s/w$. The first two entries of $A$ are divisible by $w$. After division by $w$ their residues modulo $w$ are $-t,-t/3$, since $X^n\equiv-1$ and $F\equiv1\pmod s$. If at any level $q^j\mid w$ their least residues sum to at least $q^j$, a carry in the two entries forbids $q\mid U$. In particular $2\nmid U$. Also each prime of $F$ at which 2 has even multiplicative order must go to $T$ by (5). These tests give $U\mid U_{\max}$ and $T_{\min}\mid T$. Whenever $U_{\max}^3<27T_{\min}^2$ the shape is impossible.
2. Now $s\mid N$ by $X^n\equiv-1\pmod s$ and $FM\equiv3\pmod s$. Write $S=N/s$. The companion has entries $S z_i$ with $\sum z_i=s$ and every $z_i$ a multiple of $s/U$. For each $w=q^v\parallel M$, the residues of $z_i$ are $s,0,0$ in some order. Indeed $Ss=N\equiv1\pmod M$. As $M>s$, each complete assignment specifies all three coordinates by CRT. If a reduced pair sum has a common prime with $M$, its binomial factor forces a carry. Further, $3S\equiv-1\pmod M$ since $s=M-3$. For the two zero slots at $w$, the next blocks of $Sz_i$ are $-z_i/(3w)\bmod w$. A sum of least residues at least $w$ rules out the assignment. All these tests are independent of $n$.

One quick consequence is that a full prime power $w\parallel s$ with $w^2\ge4s$ suffices to certify every outer exponent. Then $w\ge4t$ and the residues at level $w$ in test 1 are $w-t$ and at least $(w-t)/3\ge t$. Thus $U\le s/w$, contradicting (1) with $T\ge M$. Another consequence is that $M$ must have at least three distinct prime divisors, since every positive $z_i<M$ needs an exceptional slot.

## Reproduction

Run the following commands in the repository's Python environment with python-flint and SymPy installed.

```sh
python check_split_quadratic_blocks.py
python check_split_quadratic_blocks_cross.py
python check_ternary_quadratic_resonance.py
python check_ternary_quadratic_resonance_cross.py
```

The first two programs independently enumerate 19,330 shapes and return $(p,\alpha,r,d,k,c_0)=(3,3,1,2,4,-3)$ as the only residual shape. The first carries 4,455 numerical assignments to the existing companion tests. The second performs 172 direct companion tests and also verifies the signed gcd bounds by complete multiplicative orbits at five primes and 32 signed constants through inner exponent 15.

The [factor witnesses](ternary_quadratic_factors.txt) give one record per inner exponent. Each line lists primes forbidding factors of $U$, primes forcing factors of $T$, and optionally all primes of $M$. A dash denotes an empty list. Both checkers prove primality with FLINT, recompute full powers and verify every exclusion. At the 300 certified indices the size bound settles 236. The first checker eliminates two more using an individual prime power and examines 29,500 CRT rows over the other 62 indices. The second independently constructs coordinates using CRT idempotents and tests carries by integer division on the representative $n=3$. The congruences above justify that representative at both required levels for every odd $n$. Factors can be sought locally with the first checker's `--generate LIMIT --deep` option. Some factors used as hints were retrieved from [FactorDB](https://factordb.com/).
