# A spacing theorem for six base-p digit units

Let $p$ be any prime dividing $N-1$, and let $s_p(N)$ be the sum of the base-$p$ digits of $N$. Count a digit of size $b$ as $b$ units of the same exponent.

**Theorem.** If $s_p(N)\le5$, the full three-row conjecture holds. If $s_p(N)=6$, a counterexample requires a pairing of the six units with parts
$$1+p^u,\qquad p^i+p^{i+d},\qquad p^j+p^{j+e},$$
where $u,i,j,d,e$ are positive. Both $d$ and $e$ divide $u$ with **odd quotients**. Let $H$ be the difference between the greatest and least of the five positive exponents, counting multiplicities. Then
$$H<d+e,\qquad H<3\min(d,e)+\gcd(d,e).$$
Thus one can exclude a weight by testing at most 15 pairings of its digits, without factoring either $N$ or $N-1$.

In particular, the full conjecture holds at every
$$1+p^k+p^{3k}+p^{6k}+p^{7k}+p^{8k},\qquad k\ge1.$$
After dividing the exponents by $k$, the only triples $(u,d,e)$, up to swapping $d,e$, with all three having the same 2-adic valuation are $(1,3,1),(1,5,1),(3,5,1),(3,7,1),(7,5,5),(7,7,3)$. Each violates one of the two divisibilities.

**Counting consequence.** For a fixed prime $p$, at most $O(L^3)$ weights below $p^{L+1}$ with $p\mid N-1$ and $s_p(N)=6$ can be exceptions. The implied constant is absolute. The total number of such weights is of order $L^5$ as $L$ tends to infinity. In base two this improves the earlier $O(L^4)$ exceptional-set estimate for odd six-bit weights and gives an elementary test on the exponents. [Earlier sparse-weight estimate](odd_sparse_endpoint.md)

## What the two remaining rows require

The [global endpoint reduction](all_weights_endpoint_reduction.md) places $E=(1,1,N-2)$ in every counterexample. Since its coefficient is $N(N-1)$, another row $A$ must avoid $p$. Assign $p$ to $A$ and distribute the other prime bases of $N(N-1)$ between $A$ and the last row $C$, assigning each to a row whose coefficient it does not divide. Write $U,V$ for the products of the full powers assigned from $N$ to $A,C$. Then $UV=N$, $U,V\ge3$, and $U$ divides every part of $A$. Also $p\nmid U$. [Factor assignment](endpoint_reduction.md)

By Kummer, the parts of $A$ partition the base-$p$ units of $N$. A part consisting of one unit is coprime to $U$, which is impossible. This proves the assertion for at most five units. At six units each part uses two. The units digit of $N$ is 1. A pair of repeated exponents in either remaining part would force $U\mid 2p^i$ and hence $U\le2$. This proves the displayed form of the pairing.

Put $k=\min(u,i,j)$ and $m=(N-1)/p^k$. Valid base-$p$ digits give $p\nmid m$ and $m\ge p^H$. At a full power of $m$ assigned to $A$, two parts of $A$ are 0 modulo that power. If they were the last two, the first would be 1 modulo the power, which would then divide $p^u$. Thus one of the two zero positions must be the first. Write $S_d,S_e$ for the products of these powers according as the other zero is in the second or third part. Set $R=m/(S_dS_e)$. Every full prime power of $R$ divides two parts of $C$ and is coprime to $V$. The three positive parts of $C/V$ sum to $U$, so
$$R^2\mid\prod_{z\in C}(z/V),\qquad 27m^2\le U^3(S_dS_e)^2. \tag{1}$$

Define $g_d=\gcd(u,d)$, $g_e=\gcd(u,e)$, $g=\gcd(u,d,e)$ and
$$G=\gcd(p^u+1,p^d+1,p^e+1).$$
We have $U\mid G$, so $G\ge3$. The gcd of any of these two or three binomials is $p^h+1$ when their gcd exceeds 2, where $h$ is the gcd of their exponents. Indeed Bézout's identity gives $p^h\equiv\pm1$ modulo their numerical gcd. The positive sign would make that gcd divide 2, and the negative sign requires all quotient exponents to be odd. The reverse divisibility follows by factoring an odd power plus one. In particular $G=p^g+1$, and $d/g_d,e/g_e,u/g_d,u/g_e$ are odd.

Since $G$ divides $N$, it is coprime to $m$. The two zero positions defining each product therefore give
$$S_d\mid (p^{g_d}+1)/G,\qquad S_e\mid(p^{g_e}+1)/G. \tag{2}$$
Together (1) and (2) imply
$$p^H\le m\le \frac{(p^{g_d}+1)(p^{g_e}+1)}{3^{3/2}\sqrt G}<p^{g_d+g_e}. \tag{3}$$
Here $p^a+1\le(3/2)p^a$ for $a\ge1$ and $G\ge3$.

The case $p^{g_d}=2$ is impossible. Then $G=3$, $S_d=1$, and (1) gives $m\le S_e\le(2^e+1)/3<2^e\le m$. The same holds with $d,e$ interchanged. Thus $p^{g_d},p^{g_e}\ge3$. Keeping $S_e$ in (1) and using (2) gives the stronger one-sided estimate
$$S_e>3p^{H-g_d-g/2}. \tag{4}$$
Indeed $m\le S_e (p^{g_d}+1)\sqrt G/3^{3/2}$, and the factor multiplying $S_e p^{g_d+g/2}$ is at most $(4/3)\sqrt{3/2}/3^{3/2}=4/(9\sqrt2)<1/3$. The estimate for $S_d$ is symmetric.

## Rotating three exponents

We have $S_e\mid p^e+1$ and $S_e\mid p^i+p^{i+d}+p^u$ by subtracting the third part from $N-1$. On the circle of circumference $e$, place the residues of $i,i+d,u$. If an empty arc between consecutive points has length $b$, rotate the three residues into $[0,e-b]$. Algebraically this multiplies the trinomial by a power of $p$ and uses $p^e\equiv-1\pmod {S_e}$. Thus $S_e$ divides a sum of three signed powers with exponents in $[0,e-b]$. If $S_e>3p^{e-b}$ this signed sum vanishes as an integer. That is impossible for odd $p$ by parity. For $p=2$, one power must equal the other two summed, so their exponents are $a,a,a+1$. In particular $d\equiv0,1,$ or $-1\pmod e$.

Suppose $d$ fails to divide $u$. Then $d/g_d$ is odd and at least 3, so $g_d\le d/3$. By (3), $g_e>H-g_d\ge2H/3$. Since $e\le H$, we have $e=g_e\mid u$ and $d<3e/2$. We can rotate into an interval of length at most $H-d/2$:

- If $d\le e$, the arc of length $d$ from $i$ to $i+d$ has at most one further point in its interior. An empty arc has length at least $d/2$, leaving length at most $e-d/2$.
- If $e<d<3e/2$, the reverse arc between these residues has length $2e-d$. Its longer empty piece has length at least $(2e-d)/2$, leaving length at most $d/2$.

Since $g\le g_d$, (4) is strictly greater than $3p^{H-d/2}$. The rotated sum must therefore vanish. Only $p=2$ is possible. Here $d\ge3$ and $e>2H/3\ge2$, so $d<2e-1$. The congruence leaves $d=1,e-1,e,$ or $e+1$. The values 1 and $e$ divide $u$. Consecutive exponents cannot both have odd quotient on division by their gcd, so $d=e\pm1$ contradicts $G\ge3$. This proves $d\mid u$, and symmetry proves $e\mid u$.

Now $H<d+e$ by (3). Suppose $e\le d$ and $H\ge3e+g$. Then $d>H-e\ge2e+g$. For $S_d$, rotate $j,j+e,u$ modulo $d$. The arc of length $d-e$ between the first two is split in at most two pieces, leaving a span of at most $(d+e)/2\le(H+e)/2$. The symmetric (4) gives $S_d>3p^{H-e-g/2}\ge3p^{(H+e)/2}$. Vanishing would require $p=2$ and $e\equiv0,\pm1\pmod d$. Since $1\le e<d-1$, this forces $e=1$, already excluded by $p^{g_e}\ge3$. This proves the second span bound.

## Two groups of three digit units

Suppose the six unit positions, with their valid base-$p$ multiplicities, are
$$0,a,b,\quad t,t+c,t+d,\qquad 1\le a\le b\le h,\quad 0\le c\le d\le h,\quad t\ge3h.$$
**The full conjecture holds unless $(c,d)=(a,b)$.** In that remaining case $N=(p^t+1)(1+p^a+p^b)$, and the only possible row avoiding $p$ in a counterexample is $(p^t+1)(1,p^a,p^b)$ up to permutation.

Indeed $H=t+d-a\ge2h$. If the units digit were paired within the lower group, its partner exponent would be at most $h$. Both other gaps divide this exponent, contradicting their sum exceeding $H$. Thus the units digit is paired across the groups. The two remaining pairs must either both cross or both stay within groups. Internal gaps are again too small. Every cross gap is at least $t-h>(t+h)/3$, while the partner exponent $u$ is at most $t+h$. Its odd quotient by each gap must equal one. Thus the upper positions are translates of the lower positions by $u$, which is necessarily $t$. For this pairing, (2) also gives $S_d=S_e=1$ in the notation of the proof. The complementary row must avoid every prime of $(N-1)/p^{v_p(N-1)}$.

The matched groups explain a limit of the size test. Take $20\le t\le L/2$ and $1\le a<b\le t/4$ with upper positions $t,t+a,t+b$. Here $G=p^t+1$ and $m\le5p^{5t/4}$, so $27m^2\le675p^{5t/2}\le p^{3t}\le G^3$. These weights pass the numerical bound (1) with the upper bounds $G,1,1$ substituted for $U,S_d,S_e$. The ordered exponents recover $t,a,b$ uniquely, yielding $\gg L^3$ distinct weights. Further savings over the cubic count require a condition beyond this size test.

## Keeping the multiplier fixed

Suppose $m$ is a fixed positive integer coprime to $p$ with exactly five base-$p$ units. The five positive exponents of $N=1+p^km$ are the exponents of those units shifted by $k$, so their span $H$ is independent of $k$. If four of them fit in an interval of length at most $H/4$, the full conjecture holds for every $k\ge1$. At least one of the two pairs with gaps $d,e$ must lie entirely in that interval, whereas the span bound gives $H<4\min(d,e)$.

More generally, at each pairing let $c$ be the exponent in $m$ matched to the units digit of $N$, so $u=k+c$. The theorem requires $d,e$ to have the same 2-adic valuation. With $D=\operatorname{lcm}(d,e)$ it then requires $k+c\equiv D\pmod{2D}$. Since $2D\ge d+e>H$, in any interval of $K$ consecutive exponents at most $15(1+K/H)$ values of $k$ can admit a counterexample when $H>0$. All these statements are uniform in $m$ and $p$.

## Counting the weights

Fix $u\le L$. There are at most the divisor pairs $d,e\mid u$ to choose from. The starting positions $i,j$ lie within $H<4\min(d,e)$ of $u$, giving $O(\min(d,e)^2)=O(de)$ choices. With $\sigma(u)$ the sum of divisors of $u$, the weight count is at most a constant times
$$\sum_{u\le L}\sigma(u)^2\le L^3\sum_{a,b\ge1}\frac1{ab\operatorname{lcm}(a,b)}=O(L^3).$$
Indeed $\sigma(u)=u\sum_{a\mid u}1/a$; count the multiples of $\operatorname{lcm}(a,b)$ and use $\operatorname{lcm}(a,b)\ge\sqrt{ab}$ for convergence. For the total class choose five distinct positive positions out of $L$ for a lower bound, or allow all multisets of five for an upper bound.

## Finite certificates

The full conjecture holds at every weight below the indicated bound with $p\mid N-1$ and $s_p(N)=6$.

| Base $p$ | Weight bound | Pairings after the size test | After assigned-base carry tests | Complementary CRT searches | Complementary rows checked |
|---:|---:|---:|---:|---:|---:|
| 2 | $2^{160}$ | 178,067 | 16,562 | 12,013 | 77 |
| 3 | $3^{64}$ | 11,364 | 1,058 | 790 | 5 |
| 5 | $5^{40}$ | 2,757 | 261 | 167 | 0 |
| 7 | $7^{40}$ | 2,979 | 494 | 337 | 0 |

The binary range contains $\binom{159}{5}=794{,}747{,}031$ distinct weights. The search uses (1) and (2) directly, independently of the rotation argument. For fixed $u,d,e$ they bound $m$, hence $H$, so the program enumerates every feasible pair of starts. For a particular row it assigns to $A$ every full power of $N$ and of $m$ whose base it actually avoids, giving the weakest requirements on $C$. The base $p$ itself is already covered by $A$. The complementary scales and residues now determine a finite [CRT search](endpoint_reduction.md) with all positive lifts of the normalized parts. A second enumeration lists all admissible single parts by subset sums of CRT residues and matches three parts by their ordinary integer sum. It gives exactly the same complementary row sets. Every resulting row fails a required prime-base test.

Run `python3 check_six_digit_gaps.py` with python-flint installed. The checker certifies every factorization by multiplication and primality tests. A separate enumeration of digit multisets and pairings through exponent 16 reproduces all 506 positive pairings surviving the size test across the four bases. Exact binomial products also agree with the prime-base tests in 29,788 comparisons, including 2,635 prime avoidances.
