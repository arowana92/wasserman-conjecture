# Reciprocal residues at the second prime-power layer

Suppose $M=N-1=w_1w_2w_3$ is the product of full powers of three distinct primes. In a counterexample the [three-factor CRT row](endpoint_factors.md) has parts
$$A_i=k_iM/w_i,\qquad 1\le k_i<w_i,\qquad k_i(M/w_i)\equiv1\pmod{w_i}.$$
It avoids all three predecessor primes. It also avoids the primes in a unitary divisor $U\ge3$ of $N$ dividing every $k_i$.

**The three multipliers $k_i$ cannot all be equal. If $U$ is a prime power, at most one of the three multipliers can equal $U$.** With $w_1<w_2<w_3$, the latter gives the numerical restriction
$$U\left(\frac1{w_1}+\frac2{w_2}+\frac2{w_3}\right)\le1+\frac1M.$$
This follows by minimizing $\sum k_i/w_i=N/M$ subject to the multiplicity restriction. The two smallest allowed multiples are $U$ and $2U$.

Write $[t]_w$ for the least nonnegative residue and $\{x\}$ for fractional part. Define
$$\beta_i=\frac{[M/w_i]_{w_i}}{w_i},\qquad \alpha_{ij}=\frac{[A_j/w_i]_{w_i}}{w_i}\quad (i\ne j).$$
At the second base-$w_i$ digit, avoiding its prime base requires $\sum_{j\ne i}\alpha_{ij}\le\beta_i$. The other part also has a nonnegative digit. If $k_i=k_j$, then $w_i\alpha_{ij}$ is the inverse of $w_j$ modulo $w_i$. Indeed $A_j/w_i=k_jM/(w_iw_j)$, and multiplication by $w_j$ gives $k_i M/w_i\equiv1$. The two positive inverses satisfy
$$\alpha_{ij}+\alpha_{ji}=1+\frac1{w_iw_j}.$$
After multiplication by $w_iw_j$ the left numerator is 1 modulo both powers, strictly greater than 1, and less than $2w_iw_j$. If all three multipliers coincide, summing over three unordered pairs gives more than 3, whereas $\sum_i\beta_i<3$.

## A restriction for unequal multipliers

Let $u\ge2$ be any prime power dividing $N$ and all three $k_i$, whose prime base the row avoids. Set $a_i=k_i/u$. Then
$$\sum_i\{a_i^2\beta_i\}<1. \tag{1}$$
For a pair of different indices $i,j$, suppose $a_i^2\beta_i<1$ and $a_j^2\beta_j<1$. A necessary condition is
$$a_i a_j\in w_i\mathbb Z_{\ge0}+w_j\mathbb Z_{\ge0}. \tag{2}$$
In particular $a_i a_j<\min(w_i,w_j)$ excludes this case. Taking $a_i=a_j=1$ proves the second assertion on equal multipliers.

To prove (1), let $d_i=w_i\beta_i$ and write $u a_i d_i=1+h_i w_i$. Then $h_i\equiv-w_i^{-1}\pmod u$. Also $A_i/u=a_iM/w_i\equiv-a_i w_i^{-1}\pmod u$. It follows that
$$\{a_i^2\beta_i\}=\frac{[A_i/u]_u}{u}+\frac{a_i}{u w_i}.$$
The right side is less than 1 since $a_i<w_i$. The low digits of the normalized row $A/u$ sum to at most $u-1$, by carry freedom in the prime base of $u$. Using $\sum_i a_i/w_i=(1+1/M)/u$ gives (1).

For (2), set $x=w_i\{a_i^2\alpha_{ij}\}$ and $y=w_j\{a_j^2\alpha_{ji}\}$. These are nonnegative integers. The congruence defining $k_i$ gives $xw_j\equiv a_i a_j\pmod{w_i}$, and the corresponding congruence for $k_j$ gives $yw_i\equiv a_i a_j\pmod{w_j}$. Thus $xw_j+yw_i\equiv a_i a_j\pmod{w_i w_j}$. Under the bounds in the hypothesis, $x/w_i\le\{a_i^2\beta_i\}$ and $y/w_j\le\{a_j^2\beta_j\}$ by the second-digit tests. Equation (1) puts $xw_j+yw_i$ below $w_iw_j$, as is $a_i a_j$. They must be equal.

For an example where the scale prime is essential, take $N=18{,}801{,}054$ and $M=113\cdot379\cdot439$. The row $(17{,}969{,}148,446{,}463,385{,}443)$ has multipliers 9 at 379 and 439. Its coefficient avoids both of those primes, while its 3-adic valuation is 10. Here $9\mid N$.

## Searching assignments with a prescribed scale

The search fixes odd primes $u<p<q$ and lets the third factor $R>q$ vary. It tests the case where the first row avoids $u,p,q$ and $u$ is a full power in $N=pqR+1$. The two multipliers for $p,q$ have the form $ua,ub$ with $a,b>0$. Setting $D=pq-u(aq+bp)$ gives the row
$$(uaqR,\ ubpR,\ DR+1).$$
We need $D>0$ and $upq\mid DR+1$ for the third part. Conversely, these conditions give precisely the required three low congruences at $p,q,R$, and divisibility of the row and weight by $u$. We have $\gcd(D,upq)=1$, since $ua<p$, $ub<q$. Thus $R$ runs through a single progression $R_0\pmod{upq}$.

Three necessary low-digit tests are constant on this progression. At $u$ require $[aqR]_u+[bpR]_u<u$. At $p$ require $[ubR]_p\le[qR]_p$, and at $q$ require $[uaR]_q\le[pR]_q$. The first tests the normalized row at its units digit and the others compare a part's second digit to the total's second digit. They can eliminate a progression at every $R$.

For $103\le u\le800$ and $u<p<q\le3659$, exactly 21,334,937 of 223,692,214 progressions survive these tests. In the interval $q<R<12upq$, the full set contains 2,684,305,628 parameter choices. After the three fixed tests, exact unitarity and a complete base-$u$ addition leave 24,823. Of those, 101 also avoid $p$ and **none avoids $q$**. This excludes the prescribed assignment throughout the indicated range even when $R$ is allowed to be composite. The largest possible weight is below $12\cdot800\cdot3659^4+1<2^{63}$. Other assignments of the prime factors at a given weight are outside this search.

Two explicit three-prime predecessors give rows where exactly one carry prevents the proposed assignment. In each row the coefficient has valuations $(0,0,0,1)$ at $(u,p,q,R)$, respectively. Both $R$ are prime, and the three parts have the required $1,0,0$ residues at each predecessor prime.

| $N$ | $(u,p,q,R)$ | CRT row |
|---:|---|---|
| 5,480,060,281,140 | $(79,1777,7507,410801)$ | $(1461760592718,2710465408801,1307834279621)$ |
| 42,564,932,000,094,920 | $(83,23911,24767,71875487)$ | $(15957176632045956,21682043933992712,4925711434056252)$ |

In each case the remaining row assigned $V=N/u$ has normalized sum $u<2R+1$. It cannot avoid $R$, which would require two positive multiples of $R$ in that sum. The endpoint and the two rows therefore share $R$ for this attempted assignment.

The independent verifier scans $R$ directly in a smaller box, computes the two inverse multipliers, and checks Legendre's factorial floors. Separately it checks 16,827 reciprocal pairs, including 228 where both multipliers equal the tested prime power. Twelve pairs already force a predecessor carry from the modular criterion with (1). A separate positive control at $N=1{,}050{,}469{,}803{,}524$ passes the three tested bases $101,1063,1597$. Its remaining factor is $618793=7\cdot109\cdot811$, and direct factorial floors confirm all three zero valuations.

```sh
mkdir -p .work
c++ -O3 -std=c++17 check_three_factor_progressions.cpp -o .work/check_three_factor_progressions
.work/check_three_factor_progressions
python3 check_crt_reciprocity.py .work/check_three_factor_progressions
```
