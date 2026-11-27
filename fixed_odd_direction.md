# Prescribing the common primes while one primitive triple stays fixed

Fix a positive primitive triple $x=(x_1,x_2,x_3)$ whose binary ones are disjoint across its three entries. Let $u=\sum x_i$, and let $A_x(X)$ count $1\le v\le X$ for which the three entries of $vx$ still have disjoint binary ones. Fix any cutoff $B$ containing all prime divisors of $u$.

**A positive proportion of these binary-admissible multipliers $v$, uniformly for all sufficiently large $X$, have the following property.** There is a primitive positive triple $y$ of sum $v$, with $\gcd(u,v)=1$, such that at weight $N=uv$ the rows $vx,uy$ pass the lowest-layer tests at every prime of $N(N-1)$, and their common gcd $G$ with the endpoint satisfies
$$
                        \prod_{p\le B} p^{v_p(G)}=u.
$$
The implied proportion can depend on $x,B$. The same fixed triple works however large the cutoff becomes. The common gcd of the six actual entries is one.

For $x=(1,2,4)$, at every fixed $B\ge7$ the number of distinct companion sums $v\le X$ with these properties and an odd first coefficient is $\Theta_B(X^\delta)$, where $\delta=\log_2\lambda=0.55146\ldots$ and $\lambda^3=\lambda^2+1$. Throughout this family 7 is a common prime of exact valuation one and every other common prime exceeds $B$. The [more general prescription](prescribed_common_prime.md) allows an arbitrary prescribed small-prime part by choosing a direction in terms of the cutoff.

## Constructing the companion

Primitivity and the bit condition imply that exactly one entry of $x$ is odd, so $u$ is odd. Choose pairwise coprime odd integers $d,e>u$, each coprime to every prime through $B$, with
$$
                2d\equiv 2e\equiv u/p^h\pmod{p^h}
                    \qquad (p^h\parallel u).
$$
For existence, take the indicated unit classes at the powers of $u$ and the class one at the remaining primes through $B$ and at 2. Choose $d$ in this CRT class. Choose $e$ in the same class and congruent to one modulo $d$.

Set $D=de$ and determine the unit $r_0\bmod D$ by $u e r_0\equiv1\bmod d$ and $u d r_0\equiv1\bmod e$. Impose the following constraints on an odd, binary-admissible $v$:
$$
 uv\equiv1+2Dr_0\pmod{D^2},\qquad v\equiv1\pmod u,\qquad
 uv\equiv-1\pmod\ell\quad(3\le\ell\le B\text{ prime},\ \ell\nmid u).
$$
They specify one residue class modulo an odd integer $L$. Now $r=(uv-1)/(2D)$ is an integer congruent to $r_0\bmod D$. Take $v$ large enough that $r>0$. The second actual row is
$$
                  C=(uer,\ udr,\ (2de-u(d+e))r+1).
$$
It is positive and sums to $uv$. Its coordinate gcd is exactly $u$: all three are multiples of $u$, the first two have gcd $ur$, and the third is one modulo $r$. Thus $C=uy$ for a primitive $y$ of sum $v$. The full weight powers divide all entries of one row by $\gcd(u,v)=1$. The odd predecessor powers lie wholly in $d,e$ or $r$. In each of these three groups one actual entry is one and the others zero modulo the group. Binary carry freedom in $vx$ handles 2.

For $p^h\parallel u$, put $c=u/p^h$. Strip $p^h$ from every entry of $C$, which preserves its coefficient valuation at $p$ by Legendre's digit-sum formula. The first two remaining entries are $cer,cdr$. Since $2de r\equiv-1\bmod p^h$, they are congruent to $-c/(2d)$ and $-c/(2e)$ respectively, both $-1\bmod p^h$. Hence at each of the $h$ layers $w=p^j$ their residues alone sum to $2w-2\ge w$. Legendre's formula gives $v_p(K_C)\ge h$. Also $u\mid K_{vx}$ by primitivity and the factorial decrement identity. The endpoint has valuation exactly $h$ because $v\equiv1\bmod u$. Every other odd prime through $B$ misses the endpoint, and $K_{vx}$ is odd.

## Counting the scales in the prescribed class

Choose $b$ with $2^b>u$. For every admissible $0\le w\le Y$, start with $v=1+2^b w$. Take an integer $k$ with $2^b(uY+1)\le2^k<2^{b+1}(uY+1)$. Above the free block allow optional ones at $k,k+b,\ldots,k+(L-2)b$. The products with the three $x_i$ form separate binary blocks: a copy of $x$, a copy of $wx$, and separated copies of $x$ at the high positions. The full row adds without carries. Each high bit is a unit modulo $L$. Adding the option to use such a bit increases the covered residue set by at least one until it is the whole group: invariance under translation by a unit would make any nonempty set full. Thus the high bits can realize the required class. The low block recovers $w$, giving an injection into admissible odd $v\le C_{x,B}(Y+1)$. Only finitely many small $v$ need exclusion for positivity.

For completeness, $A_x(X)\ll_{c,x}1+A_x(X/c)$ for each fixed $c>1$. Read the bits of a multiplier from low to high, keeping the three multiplication carries. The carry for $x_i$ lies in $[0,x_i)$, and a transition is allowed only when the three output bits are disjoint. A state accepts a finished word if its residual carries have disjoint binary ones. From any state visited along an accepted word an accepting state is accessible in at most $s=x_1x_2x_3$ further input bits, by removing loops. Take length $h=\lceil\log_2(X+1)\rceil$, and discard the last $s+\lceil\log_2 c\rceil+2$ high bits of every padded accepted word, for $X$ sufficiently large. Each retained prefix occurs only a bounded number of times. Completing it by at most $s$ bits gives distinct accepted integers at most $X/c$ for distinct prefixes. This proves the estimate. The admissible count grows without bound, already by separated copies of the bit one. Applying the estimate to the constant in the injection proves the asserted positive proportion.

For $(1,2,4)$ the ones of $v$ must be at least three places apart. The count $F_h$ of padded $h$-bit words satisfies $F_h=F_{h-1}+F_{h-3}$ for $h\ge3$ and $F_0=1,F_1=2,F_2=3$. Bounding the initial values above and below by constant multiples of $\lambda^h$ gives the same bounds for every $h$ by induction. This proves the sharp exponent in the first section.

## The exponent for the lowest layers alone

For **every** fixed primitive positive triple $x$, allowing arbitrary digit overlap in $x$, the number of primitive companion directions $y$ of sum through $X$ that pass the lowest layers at some common weight is $X^{1+o(1)}$. Already at common scale $g=1$ a positive proportion of all integer sums $v\le X$ admit a partner. In particular the linear exponent for the number of directions is sharp on every fixed triple.

For the lower bound let $u=\sum x_i$, $d=2u+1$, $e=2ud+1$, $D=de$. Impose $v\equiv1\bmod u$ and $uv\equiv1+Dr_0\bmod D^2$ with the same inverses defining $r_0$ as above. These select an arithmetic progression of positive density. For sufficiently large $v$ set $r=(uv-1)/D$ and use the second row
$$
                     (uer,\ udr,\ (de-u(d+e))r+1).
$$
It has coordinate gcd $u$. All weight powers are assigned because $\gcd(u,v)=1$ and all predecessor powers belong to one of $d,e,r$, with actual residues $1,0,0$ in the corresponding order. Distinct sums give distinct primitive companions.

For the upper bound write $H_x=\gcd_{i<j}(x_i x_j)$, and likewise for $y$. The [single-scale criterion](primitive_row_scales.md) gives $N-1\le H_x H_y\ll_x X^{3/2}$. The assigned predecessor product $T$ of $x$ divides the fixed integer $H_x$; put $S=(N-1)/T$. The other row has normalized sum $v\le X$ and assigned predecessor product $S$. This is exactly the data counted in the [two-congruence weight bound](two_predecessor_congruences.md#recovering-the-row-from-a-product) with its parameter $K=vT\le H_x X$. On $Q/2\le N<Q$ it gives $O_x(X Q^{o(1)})$ possible weights, using the trivial count if $H_x X\ge Q$. At each such weight, choose $T\mid H_x$ and distribute the full powers of $S$ among the three exceptional positions in $3^{\omega(S)}=X^{o(1)}$ ways. Each actual part of the companion is fixed modulo $S$, with $O_x(1)$ choices below $N$ since $N/S\le H_x+1$. Its primitive direction is then determined. Summing over the dyadic weight intervals gives the claimed upper bound.

Run the [checker](check_fixed_odd_direction.py) with standard-library Python. It tests a mixed prime-power modulus at $u=21$, a repeated prime power at $u=27$, the coverage of 3,024 binary block repairs, and the first-layer construction for arbitrary fixed directions. For $x=(1,2,4)$ one explicit example has $v=2,267,760,042,017$, $N=15,874,320,294,119$ and full common endpoint gcd $16,319,685,886,490,987$, whose part on primes through 7 is exactly 7.

~~~sh
python3 check_fixed_odd_direction.py
~~~
