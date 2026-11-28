# Ratios of cuts in different rows

For a row of three positive integers of sum $N$, a *cut* is one part or the sum of two parts. Its binomial coefficient divides the row's trinomial coefficient. Suppose two rows $A,B$, together with the endpoint row $(1,1,N-2)$, have coefficient gcd one. Choose any cut $a$ of $A$ and any cut $b$ of $B$, and write $a:b=x:y$ in lowest terms. Then

$$\min(x,y)>\max\!\left(4096,\ (3125N/108)^{1/5}\right). \tag{1}$$

The fifth-root bound is analytic. The constant 4096 comes from a finite enumeration over ratios, allowing all weights. Together with the [endpoint reduction](all_weights_endpoint_reduction.md), (1) applies to every possible counterexample. Neither nonendpoint row can have a part 1. Such a cut has coefficient $N$, forcing the other row to avoid every prime of $N$. The full prime powers of $N$ would then divide each of its parts.

## A quantitative bound for two binomial coefficients

Take any integers $2\le a,b<N$, put $C_a=\binom Na$, $C_b=\binom Nb$, and set
$$D=\gcd(N,C_a,C_b),\qquad J=\gcd(N-1,C_a,C_b).$$
Write $a=gx,b=gy$ with $\gcd(x,y)=1$. Then
$$3125N<108 J^2 c^5 \min(x,y)^5,\qquad
c=\gcd\!\left(\frac{N}{\gcd(N,a)},\frac{N}{\gcd(N,b)}\right)\mid D. \tag{2}$$
In particular the common divisor $DJ$ of $N(N-1),C_a,C_b$ exceeds $(3125N/108)^{1/5}/\min(x,y)$. The bound is uniform in the larger term of the ratio.

Write $a/N=A/(cu)$ and $b/N=B/(cv)$ in lowest terms, so $\gcd(u,v)=1$. Each reduced denominator divides the corresponding binomial coefficient by $aC_a=N\binom{N-1}{a-1}$, giving $c\mid D$. There is a positive integer $w$ with
$$N=wcuv,\quad a=wvA,\quad b=wuB.$$
Write $A=fr,B=fs$ with $\gcd(r,s)=1$. Reducedness gives $g=wf$, $x=vr$, $y=us$. Put $m=frs$. Then
$$mN=cgxy,\quad a/N=m/(cy),\quad b/N=m/(cx),\quad 1\le m<\min(cx,cy).$$
Also $\gcd(m,x)=r$ and $\gcd(m,y)=s$.

For $p^e\parallel N-1$ write $\ell=\min(e,v_p(C_a),v_p(C_b))$. At modulus $p^{e-\ell}>1$, a cut attaining this minimum is congruent to 0 or 1. Indeed, if the cut residue exceeds 1 there, the two binomial parts have a carry at that layer and all subsequent layers through $e$, giving at least $\ell+1$ factors of $p$. Since $N$ and the denominators of the cuts are units at $p$, the retained power $p^{e-\ell}$ divides $A(cu-A)$ or $B(cv-B)$. Multiplying over $p$ gives
$$\frac{N-1}{J}\mid\operatorname{lcm}\big(A(cu-A),B(cv-B)\big)
  \mid f(cx-m)(cy-m).$$
For the second divisibility take the product of the two entries divided by their common factor $f$, and use $cx-m=r(cv-B)$ and $cy-m=s(cu-A)$.

Define
$$t=cx-m,\qquad h=\frac{Jmt(cy-m)}{N-1},\qquad
k=\frac{gh}{m},\qquad L=Jmt-kx.$$
The divisibility just proved implies $rs\mid h$ and $h>0$, hence $k$ is a positive integer since $g=wf$. Elimination gives
$$hN=ckxy,\qquad cLy=Jm^2t-h. \tag{3}$$
Now $cy=mN/a$ and $a>1$, so $cy-m<m(N-1)$. Hence $h<Jm^2t$ and $L$ is a positive integer. Multiplying $kx<Jmt$ by $cLy<Jm^2t$ in (3) yields
$$hLN<J^2m^3t^2\le\frac{108}{3125}J^2(cx)^5. \tag{4}$$
The maximum of $m^3t^2$ at fixed $m+t$ occurs at $m=3(m+t)/5$. Interchanging the cuts proves (2). Two other bounds used below follow from (3) and $\max_{0<u<1}u^2(1-u)=4/27$:
$$y<\frac{4Jc^2x^3}{27L},\qquad g<\frac{4Jc^3x^2}{27h}.$$

[Lean verification of the elimination, the sharp polynomial inequality and the strict fifth-power bound from the integer parameters](CrossCutRatios.lean)

## Enumeration without a weight cutoff

For $D=J=1$ order the cuts so $x\le y$. Then $c=1$, and all possibilities occur among the following parameters:
$$1\le m<x,\quad 1\le k<\frac{m(x-m)}x,\quad
L=m(x-m)-kx,\quad \gcd(x,m)\mid h,\quad h\mid km,$$
$$y=\frac{m^2(x-m)-h}{L},\quad g=km/h,\quad N=kxy/h,\quad a=gx,\quad b=gy.$$
Retain integral values with $\gcd(x,y)=1$, $y\ge x$ and $2\le a\le b<N$. The case $x=1$ is impossible since $1\le m<x$. Every prime of $N(N-1)$ divides one of $k,x,y,m,x-m,y-m$, by (3) and the definition of $h$. A base $p$ divides $C_a$ precisely when some digit of $a$ exceeds the corresponding digit of $N$.

| Range with $2\le x\le4096$ | Parameter systems | Passing base two | Passing all endpoint primes |
| --- | ---: | ---: | ---: |
| All weights | 7,432,387 | 99,489 | 162 |
| $N>10^9$ | 1,780,607 | 5,274 | 0 |

The largest candidate weight was 36,800,424,470,974,099. The largest weight among the 162 surviving pairs was 55,784,447, with cuts 65,722 and 14,718,090. The smallest surviving reduced term was 10, at $N=46$ with cuts 10 and 23.

For every parity survivor a second program factors $N$ and $N-1$ directly and checks the two factorial valuations at each prime. To check extensions of each surviving cut pair, assign each endpoint prime to a row that must avoid it. A divisibility already detected in one cut forces the prime to the other row. Across all 162 pairs there are 164 assignments. To extend a cut, split either side of size $s$ into $i,s-i$ with $1\le i\le s/2$. Its row avoids all assigned primes precisely when the cut avoids them and $\binom si$ avoids them. Enumerate the indices digitwise in one assigned base, then use factorial valuations at the others. Every assignment fails on at least one row. This proves the computational part of (1) without using an earlier bounded-weight result.

A separate enumeration for $x\le32$ dispenses with $h,k,L$. Write $b/N=j/v$ in lowest terms, so $v\mid x$ since $b/N=m/x$. For each $g<4x^2/27$ put $a=gx$. A predecessor power assigned to the first cut divides $a(a-1)$ and one assigned to the second divides $j(v-j)$. Since $v\mid N$ and $v\mid a$, it follows that
$$N-1\mid (a/v)(a-1)j(v-j).$$
Enumerating divisors on the right gives exactly the same three surviving cut pairs. An exact-binomial check also verifies (2) and (3) for all 161,700 unordered cut pairs with repetition at $3\le N\le100$.

## Reproduction

The C++ program accepts a maximum reduced term through 4096 and a strict lower weight cutoff between zero and one billion. All intermediate products fit signed 64-bit integers in this range. For example, the numerator $kx\,m^2(x-m)$ has $kx\le x^2/4$ and $m^2(x-m)\le4x^3/27$, giving $x^5/27<2^{63}$. The denominator $L(\mathrm{cutoff}+1)$ is at most $x^2(10^9+1)/4$. The reconstructed $y$ is at most $m^2(x-m)$.

~~~sh
mkdir -p .work
c++ -std=c++17 -O3 check_cross_cut_ratios.cpp -o .work/check_cross_cut_ratios
.work/check_cross_cut_ratios 4096 0 --audit > .work/cross_cut_full_audit.log
.work/check_cross_cut_ratios 4096 1000000000
python3 check_cross_cut_audit.py --audit-log .work/cross_cut_full_audit.log
lean CrossCutRatios.lean
~~~

The Python audit uses python-flint. The Lean file uses Lean 4.36.0 and its standard library. Further results for ratios near the fifth-root bound appear in [the case of small $hL$](cross_cut_critical.md).
