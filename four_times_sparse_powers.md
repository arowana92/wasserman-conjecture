# Composite prime-power cofactors at a = 4

Let $b=q^m$ with $q$ prime and $m\ge2$, and let $n=4b+1$. Write $L_j$ for the exponent of the unit group modulo $\operatorname{lcm}(1,\ldots,j)$. For all but finitely many such $b$, Wasserman's three-row conjecture holds at **every** weight
$$
 N=1+4b\,t^e,\qquad t\ge1,\quad \omega(t)\le2,\quad e\ge1,\quad L_{n^3}\mid e .
$$
This includes **every prime $q\le 1{,}000{,}000$ and every $m\ge2$**. The exponents of the two additional prime powers may also be chosen independently as positive multiples of $L_{n^3}$. The finiteness proof gives an absolute effective upper bound for $m$. The cutoff on $q$ in the computation is explicit, while the finiteness argument gives no effective final cutoff on $q$.

## An exact binary description

The [static cofactor theorem](static_power_transfer.md) already handles $q=2$. For odd $b$ its only possible exceptions have $b\equiv2\bmod 3$ and a carry-free binary sum $(b+1)+3b=4b+1$. These conditions together are equivalent to
$$
 b=2^h(2^{6j+3}+1)/3-1,\qquad h\ge2,\quad j\ge0. \tag{1}
$$
In particular this binary test leaves exactly $\sum_{j\ge0}\max(0,L-6j-3)$ odd integers $b<2^L$, before testing primality or any other prime.

For the equivalence, write $x\preceq y$ if every one-bit of $x$ occurs in $y$. The cut passes exactly when $b+1\preceq4b+1$. Bit one of the total is zero, so $4\mid b+1$. Write $b+1=2^h u$, $h\ge2$, $u$ odd. After dividing the two words first by $4$ (rounding the total down) and then by $2^{h-2}$, the test is $u\preceq4u-1$. The two low bits of $4u-1$ are ones, and its higher word is $u-1$. Thus bit two of $u$ is zero and each bit $u_i$ for $i\ge3$ is at most $u_{i-2}$. The odd positions occupied are precisely $1,3,\ldots,2k-1$ for some $k\ge0$. Consequently $u=(2\cdot4^k+1)/3$. Since $4^k\equiv1+3k\bmod9$, the condition $3\mid u$ is $k\equiv1\bmod3$. Setting $k=3j+1$ proves (1) in both directions. Here $u$ has $6j+2$ bits and $2^h u-1$ has $h+6j+2$ bits, which gives the count.

## Uniformly bounding the exponent

The following argument works even for an arbitrary integer base $q>1$. If $q^m$ satisfies (1), then $q\equiv11\bmod12$, $m$ is odd, and $h=v_2(q+1)$ because $(q^m+1)/(q+1)$ is odd. For $m\ge3$ put $A=h+6j+3$. Then
$$
 3q^m=2^A+2^h-3,\qquad 4\le2^h\le q+1.
$$
Divide $A=mk+r$, $0\le r<m$. Folding the leading dyadic exponent into the rational base gives
$$
 0<|\Lambda|=\left|(2^k/q)^m 2^r 3^{-1}-1\right|
 =\frac{2^h-3}{3q^m}<q^{1-m}.
$$
Here $2^k<2q$, so the logarithmic height of $2^k/q$ is at most $2\log q$. Apply the real rational case of Matveev's theorem with the three positive rationals $2^k/q,2,3$ and exponent bound $m$. Omit $2$ when $r=0$. Using natural logarithms, the coefficient $1.4\cdot30^6\cdot3^{4.5}\cdot2\log2\log3$ is below $10^{13}$, as is the corresponding two-logarithm coefficient. Thus
$$
 (m-1)\log q <10^{13}\log q\,(1+\log m).
$$
The ratio $(m-1)/(1+\log m)$ is increasing for $m>1$ and exceeds $10^{13}$ at $10^{16}$. This proves $m<10^{16}$ for all bases. [Matveev's paper](https://www.mathnet.ru/eng/im314) · [Explicit multiplicative formulation, Theorem 2.1](https://doi.org/10.1016/j.jnt.2020.05.005)

Fix one of these finitely many $m$ and one residue $r$. The positive algebraic number $\gamma=(2^r/3)^{1/m}$ is irrational, since the 3-adic valuation of its $m$th power is $-1$. The equation gives
$$
 0<(q/2^k)^m-\gamma^m=\frac{2^h-3}{3\cdot2^{km}}.
$$
For fixed $m,r$, the ratio $q/2^k$ is bounded above and below by positive constants since $2q^m<2^A<3q^m$. Factoring the difference of powers gives $|q/2^k-\gamma|\ll_m 2^{-(m-1)k}$. These are reduced dyadic fractions. Ridout's theorem gives just finitely many reduced dyadic rational approximations to an algebraic irrational with error less than $(2^k)^{-3/2}$. Since $m\ge3$, this applies for all sufficiently large $k$. Each remaining bounded $k$ allows only finitely many $q$. Taking the finite union proves the assertion. [Ridout, *Rational approximations to algebraic numbers*](https://doi.org/10.1112/S0025579300001182)

## Checking every exponent through a specified prime base

For $3\le m<128$ odd the verifier tests the binary cut on the exact integer $q^m$. For $m\ge128$ we have $A\ge h+128$ from $2^A>2q^m$, $2^h\le q+1$ and $q\ge11$. Consequently
$$
 (-q)^m\equiv1-2^h/3\pmod{2^{h+128}} .
$$
Here $-q\equiv1\bmod 2^h$ with $v_2(-q-1)=h\ge2$. Its order at this modulus is $2^{128}$: each squaring increases $v_2(x-1)$ by exactly one when $x\equiv1\bmod4$. Thus there is exactly one possible exponent modulo $2^{128}$, recovered one bit at a time. Independent modular exponentiation verifies its value and the exact order for every base.

Run **python3 check_four_times_sparse_powers.py**. The binary description agrees with direct enumeration for all odd $b<2^{20}$, including all 33 survivors. The prime sieve finds 19,654 bases $q\equiv11\bmod12$ through one million. None passes at a small odd exponent $m\ge3$. The smallest lifted exponent residue over these bases is $9{,}319{,}914{,}061{,}748{,}364{,}423{,}010{,}636{,}459{,}303{,}371$, already greater than the analytic bound. [Verifier](check_four_times_sparse_powers.py)
