# Relations of height six in the two nonendpoint rows

Let $A,B$ be positive trinomial rows of weight $N$. Assume each satisfies a nonzero homogeneous integer relation with all coefficients of magnitude at most six. Then:

* If at least one of the two chosen relations uses all three parts, the coefficients share a prime dividing $N(N-1)$ at **every** weight.
* The same conclusion holds without that condition at every $N\le10^{18}$.
* Over all weights, at most $O(X^{0.359})$ integers $N\le X$ can fail the assertion under these height-six hypotheses. For the particular pair of within-row ratios $1:1$ and $5:6$ the count is $O_\epsilon(X^\epsilon)$ for every $\epsilon>0$.

Consequently, through $10^{18}$ the full three-row conclusion holds whenever each of the three rows has such a relation. The [endpoint reduction](all_weights_endpoint_reduction.md) deals with triples without the endpoint. The unrestricted finite verification is [through $10^9$](full_range_1b.md).

## A relation involving all three parts

Assume an endpoint failure. Write $N=UV$, $M=N-1=ST$ for the [prime-power assignment](endpoint_pair_bound.md); $U,V\ge3$. At every full predecessor power, the avoiding row has residues $0,0,1$. For a row on a homogeneous relation with all three coefficients nonzero this full power divides one of those coefficients. If both relations use all three parts, $M\mid60$, settling those weights.

Otherwise one row has two parts $ag,bg$ in a reduced ratio $a:b$ of height at most six. Write the companion as $V(x,y,z)$ with sum $x+y+z=U$. Its relation is $rx+sy=tz$ after permuting and changing sign, with $1\le r\le s\le6$, $1\le t\le6$, $\gcd(r,s,t)=1$. Put
$$
 D=\gcd(M,\operatorname{lcm}(a,b,r,s,t)).
$$
The [ratio assignment](affine_endpoint_gap.md) shows $M/D\mid g$ and $U\mid g$. Since $U$ is coprime to $M$ and the third part is at least $U$, this gives $(a+b)U<D$. Enumerate $D$ over divisors of the indicated lcm, $U\ge3$ coprime to $D$ below this bound, and the solutions of the displayed relation and sum. Against **every** pair of parts of $(x,y,z)$ in reduced ratio $i:j$ we must have $\operatorname{lcm}(a,b,i,j)>(a+b)(i+j)$. The only remaining tuple, up to ordering of the companion, is
$$
           (a,b,D,U,x,y,z)=(3,5,60,7,1,2,4).
$$
The first row is even in coefficient because $3g$ and $5g$ have the same lowest binary bit. In the companion $V$ is odd. If $V,2V,4V$ add without carries, both lowest bits of $N$ are one. This contradicts $60\mid N-1$. Two enumeration orders check the finite reduction in the programs below.

## Two ratios

Write the pairs as $ag,bg$ and $ch,dh$, and $L=\operatorname{lcm}(a,b,c,d)$. In an endpoint failure,
$$
                    N(N-1)\mid\operatorname{lcm}(g,h,L).               \tag{1}
$$
Every full power in $N$ divides $g$ or $h$. A full predecessor power divides the scale of its assigned row if the exceptional position is outside the pair; otherwise the scale is a unit and the power divides the term at the vanishing position. This proves (1).

Enumerate pairs of reduced ratios through six. The [sharp constant test](part_ratio_predecessor.md) leaves the following table. Just the improvement $C(3,4)=6$ over the product 12 is enough to give this same table.

| Second ratio $c:d$ | First ratios $a:b$ | Bases required in the second row |
|---|---|---|
| $2:5$ | $1:3$ | 2, 5 |
| $3:5$ | $1:2$, $1:4$ | 3, 5 |
| $4:5$ | $1:1$, $1:3$, $2:3$ | 2, 5 |
| $5:6$ | $1:1$, $1:4$ | 3, 5; also 2 against $1:1$ |

For every entry $L<2(a+b)(c+d)$. Positivity and (1) therefore force $Lgh=N(N-1)$ and $\operatorname{lcm}(g,h,L)=Lgh$. Thus $g,h,L$ are pairwise coprime. Since $N\mid gh$, each full power $w=p^e$ of $L$ is also the full power in $N-1$. A row avoiding $p$ must have one ratio term divisible by $w$, and its unit term times the scale must be $1\bmod w$. In all eight pairs exactly one row can take each such prime.

For the row with scale $h$ it is necessary that $ch,dh$ add without a carry at both displayed bases, with the specified unit congruences. For each candidate $h$ one can enumerate all possible weights. Write $h=vt$ with $v=\gcd(N,h)$ a unitary divisor of $h$. Then $N=1+Ltz$ for the unique candidate $0<z<v/(a+b)$ satisfying $Ltz\equiv-1\bmod v$ and $g=Nz/v$. This follows by splitting the coprime factors $N,N-1$ in $Lgh$ and using $(a+b)g<N$.

| Second ratio | Two-base scale candidates through $N=10^{18}$ | Row pairs after primary congruences | Row pairs in independent CRT enumeration |
|---|---:|---:|---:|
| $2:5$ | 194,432 | 1,910 | 16,996 |
| $3:5$ | 94,202 | 7,700 | 44,779 |
| $4:5$ | 379 | 78 | 496 |
| $5:6$ | 49,099 | 5,791 | 26,452 |

The [first verifier](check_small_ratio_pairs.py) intersects carry-graph prefixes, splits the unitary divisors, and checks the prime assignments. The [independent scale sieve](check_two_base_ratio_scales.cpp) uses factorial floors. At any violating power $q$ the violation persists until $\lfloor ch/q\rfloor$ or $\lfloor dh/q\rfloor$ increases, giving a provable skip. Its [weight verifier](check_small_ratio_pairs_cross.py) instead enumerates all roots of $N(N-1)=0\bmod h$ and imposes $N=1\bmod L$. Since $0<N\le Lh$, there is at most one weight per root. Every reconstructed pair in the larger enumeration shares an endpoint prime in $\{2,3,5\}$. This last certificate concerns the row pairs left by the scale and product conditions.

## The count beyond the finite range

The carry graph at base $p$ for $ch,dh$ records the multiplication carries for a common input digit of $h$ and retains only transitions with output digit sum below $p$. For block length $B=4096$ let $W$ be the total number of paths of length $B$ from all reachable states, with arbitrary endpoints. The first verifier checks $W^{1000}<p^{B u}$ by exact integer arithmetic for the following numerators $u$.

| Ratio | Base 2 | Base 3 | Base 5 |
|---|---:|---:|---:|
| $2:5$ | 631 | | 728 |
| $3:5$ | | 644 | 710 |
| $4:5$ | 574 | | 699 |
| $5:6$ | 608 | 644 | 693 |

At each base use the compact set of numbers with base-$p^B$ digits among the block labels, reversing low-to-high reading to place the digits in their usual numerical order. The zero block is allowed. This set is invariant under multiplication by $p^B$ modulo one and has dimension bounded strictly by $u/1000$. For $h\le X$ normalize using the least power $p^{Bn}>X$. The point $h/p^{Bn}$ belongs to this compact set by padding the input word and its fractional tail with zeros. Therefore $h/X$ is in a diagonal rescaling of this set with scale bounded by the fixed number $p^B$.

Corso and Shmerkin's uniform intersection estimate applies because the prime-power bases used are pairwise multiplicatively independent. At resolution $1/X$ the number of compatible $h$ is $O_\epsilon(X^{\max(0,\sum s_p-(j-1))+\epsilon})$ for $j$ required bases of dimensions $s_p$. Each $h$ yields at most $2^{\omega(h)}=O_\epsilon(X^\epsilon)$ weights. The table gives strict exponent below $0.359$ in every two-base case, and dimension sum below two in the three-base case. [Corso and Shmerkin, Corollary 1.17](https://arxiv.org/pdf/2409.04608)

There is an elementary count for **any two fixed ratios**, regardless of their terms. Retain $L$ from (1); then $Lgh=kN(N-1)$ with $1\le k<L/((a+b)(c+d))$. For given $k,g$, the modulus $g/\gcd(g,k)$ divides $N(N-1)$ and $N\le1+Lg/(k(c+d))$. The CRT thus gives $O_{a,b,c,d,\epsilon}(X^\epsilon)$ possibilities per $g$ through weight $X$. At least one row must be odd. The [binary multiple bound](direction_carry_counts.md) gives $O_{a,b}(X^{\log_2\varphi})$ choices for $g$ when that row is odd, with $\varphi=(1+\sqrt5)/2$; interchange the rows in the other case. Hence the exceptional count is $O_{a,b,c,d,\epsilon}(X^{\log_2\varphi+\epsilon})$. The same exponent works for two arbitrary fixed homogeneous relations: with two nonzero-coordinate normals $N$ is bounded, and in the mixed case the bound on $U$ reduces the companion to finitely many ratios.

Run from the repository root with Python, python-flint and a C++17 compiler.

```sh
python3 check_small_ratio_pairs.py
c++ -O3 -std=c++17 check_two_base_ratio_scales.cpp -o /tmp/check_two_base_ratio_scales
python3 check_small_ratio_pairs_cross.py /tmp/check_two_base_ratio_scales
```
