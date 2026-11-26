# Adjacent cuts stay away from the center and every fixed rational proportion

Let $1\le a<a+1\le N/2$, put $M=N-1$, $d=N-2a-1$, and suppose
$$\gcd\left(NM,\binom Na,\binom N{a+1}\right)=1. \tag{1}$$
Then
$$9d^6>4N^4(N-1). \tag{2}$$
Moreover, for each fixed $C>0$, only finitely many such pairs have $d\le C N^{5/6}$. For any fixed rational $0<r/t<1/2$ in lowest terms, all sufficiently large pairs satisfying (1) obey $|a-rN/t|\ge c_{r,t}N^{2/3}$ for some $c_{r,t}>0$.

Through $N=10^9$, (1) has **no solutions at any position**. Thus throughout this range any two adjacent interior binomial coefficients share a prime divisor of $N(N-1)$. Reflection covers the upper half of the row. The two equal middle coefficients at odd weight $N$ are each divisible by $N$, by $N\mid a\binom Na$ and $\gcd(N,(N-1)/2)=1$.

In a [trinomial counterexample](all_weights_endpoint_reduction.md), any two adjacent cuts of the two nonendpoint rows must satisfy (1), since the cut binomials divide the row coefficients. The estimates therefore apply at every weight to such a pair.

## Four integer parameters

Each full prime power of $N$ divides $a$ or $a+1$ by the low digits of an avoiding binomial. For each full power $w$ of $M$, the cut in its avoiding binomial is $0$ or $1$ modulo $w$, so $a$ is $-1,0$ or $1$. Write
$$g=\gcd(a,M),\quad i=M/g,\quad a=gc,\quad M-a=gf.$$
Here $c,f$ are positive coprime integers with sum $i$. If a prime of $M$ divides $a$, its full power divides $a$. Every remaining full power divides $U=a+1$ or $V=N-a$. Thus $i\mid UV$. Each full power of $N$ also divides $U$ or $V$, so $Ni\mid UV$. We obtain positive integers $\ell,j$ satisfying
$$UV=Ni\ell=N+g^2cf,\qquad cf=Nj,\qquad \ell i=1+g^2j. \tag{3}$$
In particular $\gcd(g,\ell)=1$ and $i\ge2\sqrt{Nj}$. The difference $\delta=g-4\ell$ satisfies
$$d^2=(V-U)^2=3N+1+Ni\delta. \tag{4}$$
A negative $\delta$ would make the right-hand side negative since $i\ge4$ and $N\ge4$.

## Why equality is impossible

If $\delta=0$, (4) gives $3N=d^2-1$, and $4\mid M$ since $g=4\ell$. In particular $d$ is even. An odd binomial requires each side of its cut to have ones only in positions where $N$ has ones. At least one of the two binomials must be odd by (1).

The following table uses $3N=d^2-1$ modulo 16 and $2a=N-d-1$. It lists the residues of $d$ for which at least one of the two cuts passes the lowest three binary places.

| $d\bmod16$ | 0 | 2 | 4 | 6 | 8 | 10 | 12 | 14 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| $N\bmod8$ | 5 | 1 | 5 | 1 | 5 | 1 | 5 | 1 |
| $a\bmod8$ | 2 | 7 | 0 | 5 | 6 | 3 | 4 | 1 |
| Either cut passes | no | yes | yes | no | no | no | yes | yes |

Choose $r\in\{2,4,-4,-2\}$ in the surviving class of $d$. The weights $N_0=(r^2-1)/3$ are 1 or 5. With $a_0=(N_0-r-1)/2$, each of the two cuts $a_0,a_0+1$ has at least one side $c_0$ in $[0,5]$, allowing the other side to be negative at this formal representative. If $d\ne r$, let $k=v_2(d-r)\ge4$. From $3(N-N_0)=(d-r)(d+r)$ we get $N\equiv N_0\bmod 2^{k+1}$. Every side $c$ of either actual cut satisfies
$$c\equiv c_0+2^{k-1}\pmod{2^k}$$
for its corresponding representative, by the formula for $2a$. Choose the side with $0\le c_0\le5$. It has a one in place $k-1$, where $N$ has a zero. Both binomials are therefore even. If $d=r$, direct substitution gives no positive pair in the stated range. This proves $\delta\ge1$ for every solution of (1).

## The exponent and the finite set at its boundary

Reducing $\ell i=1+g^2j$ modulo $\ell$ gives the positive integer $e=(1+\delta^2j)/\ell$. Consequently
$$g=4\ell+\delta\le 4+4\delta^2j+\delta\le9\delta^2j,\qquad
j=\frac{cf}{N}\le\frac{M^2}{4Ng^2}.$$
It follows that $g^3\le 9\delta^2 M^2/(4N)$. Cubing $d^2>NM\delta/g$ now gives (2), including the stronger factor $\delta$ on its right.

Suppose $d\le C N^{5/6}$. The stronger bound makes $\delta$ bounded. Equation (4) gives $g\gg_C N^{1/3}$, so $j=O_C(N^{1/3})$ and $e=(1+\delta^2j)/\ell=O_C(1)$. Fix the positive integers $\delta,e$. Let $h=f-c$ and $z=4\delta j+e$. Equations (3) give
$$h^2=i(i-4gj)-4j,\qquad i=4gj+4\delta j+e,$$
and substitution yields
$$\delta e h^2=\delta z^3+(z-e)(4z-e). \tag{5}$$
The cubic on the right has discriminant $e^2(144+140v-27v^2)$ with $v=\delta e$. The factor in parentheses is positive for integer $1\le v\le6$ and negative for $v\ge7$. Thus (5) is nonsingular and Siegel's theorem gives only finitely many integral pairs $(z,h)$. Each $z$ determines $j,\ell,g,i,N$. This proves the finiteness claim, without giving a largest weight. [Siegel's integral-point theorem](https://doi.org/10.1007/978-88-7642-520-2)

## Other rational proportions

Fix $0<r/t<1/2$, put $A=r(t-r)$, $B=t^2$, $\theta=A/B$, and $\Delta=Ag-B\ell$. If $\Delta\ne0$, (3) gives $\ell\mid A^2+\Delta^2j$. When $|\ell/g-\theta|\le\theta/2$,
$$\frac{\theta g}{2}\le\ell\le (A^2+1)\Delta^2 j,\qquad
\left|\frac{\ell}{g}-\theta\right|^3\ge\frac{2\theta N|\Delta|}{B^3(A^2+1)M^2}.$$
On writing $y=a/N$, the exact difference $\ell/g-y(1-y)$ equals $(1-y^2)/M$. Also $|y(1-y)-\theta|\le|y-r/t|$. These estimates give the asserted separation for $\Delta\ne0$, including the case $|\ell/g-\theta|>\theta/2$ for large $N$.

When $\Delta=0$ we have $UV=\theta NM$. Put $q=t-2r>0$ and $E=t^2+2A$. Expanding $d^2=(N+1)^2-4\theta NM$ gives
$$(q^2N+E-qt d)(q^2N+E+qt d)=4A(2t^2+A).$$
Both factors are positive integers, so only finitely many $N$ occur at fixed $r,t$.

## Finite verification

Run `python3 check_adjacent_cuts.py`. From (3) and $\delta\ge1$, the program enumerates $5\le g<\sqrt N/2$, $1\le\ell\le(g-1)/4$ coprime to $g$, and positive $j$ in the class $g^2j\equiv-1\bmod\ell$. The cutoff $N\le X$ gives $g^3j\le(X-1)\ell-g$. For each choice $i=(1+g^2j)/\ell$ is fixed. Requiring $i(i-4gj)-4j$ to be a square recovers the only possible $c\le f$, with their parity and coprimality checked directly.

At $X=10^9$ there are 35,767,295 parameter triples and 192 reconstructed index pairs, of which 183 pass all the low-power checks. Each of the 192 fails the binary carry test. An independent scan of the actual prime powers and cut indices through weight 4,096 checks the enumeration in that range. The program also computes 65,280 gcds directly from exact adjacent binomial coefficients through weight 512.

[AdjacentCuts.lean](AdjacentCuts.lean) checks the residue inequality at the lifting step and the cubic identity (5). Run `lean AdjacentCuts.lean` with Lean 4.36.0. The unrestricted adjacent-binomial assertion beyond $10^9$, and the full trinomial conjecture, remain open here.
