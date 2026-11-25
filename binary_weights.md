# Weights of the form $m2^k+1$

**Theorem 1.** For every integer $k\ge0$, any three positive trinomial coefficients of weight $5\cdot2^k+1$ have a common prime divisor.

**Theorem 2.** Fix any positive odd integer $m$. For all sufficiently large $k$, any counterexample of weight $N=m2^k+1$ must include $(1,1,N-2)$, up to permutation. For each of the following choices of $m$, the full three-row conjecture holds for all sufficiently large $k$.

- The binary expansion of $m$ has at most four ones.
- $m$ is a power of an odd prime.

The threshold in Theorem 2 may depend on $m$. Theorem 1 has no exceptional exponents. Infinitely many of its weights have no universal prime. For example, $k\equiv2\pmod6$ gives $N\equiv3\pmod9$ and $N>3$, so $N$ is not a prime power. In an odd prime base, an odd $N$ cannot have digit sum two. In base two this $N$ has three ones. In each base there are therefore enough digit units to form a carry-free triple of positive parts.

## Removing one row from the weighted-product bound

Suppose three rows have multinomial gcd one. Say that a row *avoids* $p$ if its coefficient is coprime to $p$. By Kummer the sum of its parts then has no carries in base $p$. In particular, the residues of its parts modulo $w=p^e$ sum to the residue of $N$ modulo $w$. Every prime must be avoided by some row.

Take $k\ge1$. If $m=1$ then $N$ has only two binary digit units, so 2 divides every coefficient. Otherwise take a row $A$ avoiding 2. Its parts have the form
$$A_i=2^k x_i+\delta_i,\qquad \sum_i x_i=m,\qquad \sum_i\delta_i=1,\qquad \delta_i\in\{0,1\}.$$
The $x_i$ have disjoint binary digits. At least two are positive because all three parts of $A$ are positive. Put $d_0=\gcd(x_1,x_2,x_3)$, $z_i=x_i/d_0$, and $d=m/d_0$. Then $d$ is odd, $d\ge3$, $\sum z_i=d$, $\gcd(z_1,z_2,z_3)=1$, and
$$d A_i=z_i(N-1)+d\delta_i.$$

Fix $H$, and use offsets $J=\{0\le h<H:d\nmid h-1\}$. Suppose $w=p^e>H$ divides $N-h$ for $h\in J$ and $A$ avoids $p$. The residues $u_i=A_i\bmod w$ sum to $h$. Consequently
$$w\mid z_i(h-1)+d(\delta_i-u_i)\quad\text{for every }i.$$
For some $i$, $d\nmid z_i(h-1)$ by the gcd condition on the $z_i$. That integer on the right is nonzero and its absolute value is at most $2m(H+1)$. Thus every layer $w>2m(H+1)$ at these offsets must be covered by one of the other two rows $B,C$.

Write $W(t)=t(t+1)/2$ and $B_t(x)=\prod_{j=0}^{t-1}(x-j)^{t-j}$. For each of the rows $B,C$ choose nonnegative cutoffs $t_i\le x_i'$ on its actual parts $x_i'$, summing to $H$ in that row. Write $P$ for the product of the six resulting $B_{t_i}(x_i')$. Put $L=\sum_{h\in J}(H-h)$. For $N\ge2H$,
$$N^L\ll_{m,H}P. \tag{1}$$
Indeed, let $K$ be the lcm of all integers through $2m(H+1)$. Remove $\gcd(N-h,K)$ from each $N-h$. Each layer left over occurs at a unique $h$. In a covering row with residues $u_i$, the factors $x_i'-u_i$ with $u_i<t_i$ supply at least $\sum_i\max(0,t_i-u_i)\ge H-h$ copies at that layer. Hence the product of the remaining factors raised to $H-h$ divides $P$. Since $N-h\ge N/2$, (1) follows. This is the layer argument of the [weighted divisibility lemma](weighted_obstruction.md).

## Proof of Theorem 2

Call a row primitive if its three parts have gcd one. For $p\mid N$, a row avoiding $p$ must have every part divisible by $p^{v_p(N)}$. If two rows are primitive, the third must have every part divisible by $N$, which is impossible. An imprimitive row at our odd weight has every part at least 3.

Sort each row increasingly. In a counterexample the sum of the smallest two parts is either 2 or at least 11 by Bergman’s Proposition 10. A sum of 2 gives precisely the row $(1,1,N-2)$, denoted by $E$. Let $r$ be the largest prime at most $N$ and $g=N-r$. For large $N$, $r>N/2$, so a row avoiding $r$ has two smaller parts summing to at most $g$. Guth–Maynard’s Corollary 1.3 gives $g=O(N^{57/100})$. With $m$ fixed we therefore have $g<2^k$ eventually. The row $A$ has at least two parts $\ge2^k$ and cannot avoid $r$. Designate a row among $B,C$ that does. [Bergman, Proposition 10](https://arxiv.org/abs/0806.0607) · [Guth–Maynard, Corollary 1.3](https://ora.ox.ac.uk/objects/uuid%3Aad11b8bf-ad2b-4ebf-a627-647f023c378f/files/s73666708v)

Use $H=12$. For any $d\ge3$ the omitted offsets $1,1+d,\dots$ cost at most $11+8+5+2=26$ out of $W(12)=78$, giving $L\ge52$. The following costs bound the product for one row. They use $B_t(x)\le x^{W(t)}$, treating $B_a(a)$ as an absolute constant when $a\le4$. For the designated row the second part is at most $g$.

| Smallest part $a$ | Cutoffs in ordinary row | Cost | Cutoffs in designated row | Cost |
|---|---|---|---|---|
| $a\ge5$ | $(4,4,4)$ | $N^{30}$ | $(5,5,2)$ | $g^{30}N^3$ |
| $a=4$ | $(4,4,4)$ | $N^{30}$ | $(4,5,3)$ | $O(g^{15}N^6)$ |
| $a=3$ | $(3,4,5)$ | $O(N^{25})$ | $(3,6,3)$ | $O(g^{21}N^6)$ |
| $a=1,2$, row different from $E$ | $(a,5,7-a)$ | $O(N^{36})$ | $(a,8,4-a)$ | $O(g^{36}N^6)$ |

The lower bound of 11 on the two-part sum ensures every displayed cutoff fits. If $A$ is primitive, $B,C$ both have minimum at least 3. The largest total exponent of $N$ in the first three lines, using $g=O(N^{0.57})$, is $30+3+30(0.57)=50.1<52$. This contradicts (1) for large $N$.

If $A$ is imprimitive, every $x_i$ is positive, since $x_i=0$ would give a part of 1. Its three positive $x_i$ are distinct by disjointness of their binary digits. Thus the $z_i$ are distinct positive integers, and their odd sum $d$ is at least 7. This improves $L$ to at least $78-11-4=63$. Suppose neither $B$ nor $C$ is $E$. At most one of them contains a 1 or 2, since such a row is primitive at odd $N$. If that row is designated the exponent is at most $30+6+36(0.57)=56.52$. If it is ordinary the exponent is at most $36+3+30(0.57)=56.1$. The first three lines suffice if all parts are at least 3. All cases contradict (1). This proves the required presence of $E$.

If $m$ has at most four binary ones, $N$ has at most five. Three positive parts formed without a carry therefore include a part consisting of just one binary unit. This power of two forces $A$ to be primitive because $N$ is odd. The primitive case above already excluded a counterexample. If $m$ is an odd prime power then $N-1$ has just two distinct prime factors. The [endpoint reduction](endpoint_reduction.md) excludes $E$ at every such weight, proving the other assertion.

## A uniform tail bound for $m=5$

For $k\ge1$ the only row avoiding 2 is $A=(1,2^k,4\cdot2^k)$ up to permutation. It is primitive, so in a counterexample $B,C$ have all parts at least 3. Here $d=5$. Use $H=15$ and $J=\{0,\ldots,14\}\setminus\{1,6,11\}$, giving $L=120-(14+9+4)=93$.

For each $h\in J$ set
$$K_h=\operatorname{lcm}\big(1,\ldots,15,\ |h-1-5u|\ (0\le u\le h)\big).$$
Each integer inside the absolute values is nonzero. Indeed $h\not\equiv1\pmod5$. If $w>15$ divides $N-h$ and $A$ avoids its underlying prime, apply the congruence above to the part $2^k$ with $z_i=1,\delta_i=0$. Then $w$ divides $K_h$. The layer argument gives the exact divisibility
$$F(N)=\prod_{h\in J}\left(\frac{N-h}{\gcd(N-h,K_h)}\right)^{15-h}\ \mid\ P. \tag{2}$$

For $N\ge30$ each of $B,C$ costs at most $N^{45}$. If all parts are at least 5, take $(5,5,5)$. For minimum 3 take $(3,6,6)$, costing at most $108N^{42}$, and for minimum 4 take $(4,5,6)$, costing at most $27648N^{36}$. The two-part sum bound ensures these cutoffs fit, and both constants are absorbed for $N\ge30$. Thus a counterexample requires $F(N)\le N^{90}$. In particular, using $N-h\ge N/2$, it requires
$$N^3\le D,\qquad D=2^{93}\prod_{h\in J}K_h^{15-h}.$$
The checker computes $D$ with integers and finds that its bit length is 2,216, giving cutoff 739. For $k\ge\lceil\operatorname{bitlength}(D)/3\rceil$, $N^3>D$. Below this cutoff, testing (2) leaves only 81 positive exponents, the largest 124. This argument for $m=5$ uses no prime-gap estimate.

## Checking the remaining exponents of Theorem 1

For each remaining exponent, choose a certified prime $4\cdot2^k<r\le N$ and put $g=N-r$. The binary row $A$ cannot avoid $r$, so another row must be $(i,j,N-i-j)$ with $1\le i\le j$ and $i+j\le g$. Write $u=\gcd(N,i,j)$ and $v=N/u$. The last row must have all parts divisible by $v$. Indeed, if $p^e\parallel N$ and $p$ divides $v$, the middle row has some part not divisible by $p^e$. It cannot avoid $p$, and neither can $A$ since it contains 1. The last row must have every part divisible by $p^e$. Taking all such primes gives the stated multiple. All candidates are thus covered by looping over sorted positive triples $(a,b,c)$ summing to $u$ and testing $(va,vb,vc)$.

The loop examines 8,203 pairs of rows. In every case one of 3 and 5 divides all three coefficients. Both a carry test and Legendre’s factorial-valuation formula check each pair. For $k=0$ the weight is 6, which has only two binary digit units. This completes Theorem 1.

Run the standard-library verifier from this directory.

```sh
python3 check_binary_ray.py
```

The [prime data](binary_ray_primes.txt) are certificates of primality. For each entry $p$ the verifier factors $p-1$ completely using previously certified primes, then checks a witness of multiplicative order exactly $p-1$ modulo $p$. This forces $\varphi(p)\ge p-1$ and hence $p$ prime. Trial division handles entries below $2^{16}$. The choice of the nearby primes and the checks of the bounds and row pairs can all be reproduced offline.
