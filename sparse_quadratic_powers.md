# The quadratic weights $P(3P+2)$

Let $P=p^a$ for a prime $p$ and $a\ge1$, and put $N=P(3P+2)$. Each of the following hypotheses suffices for **every three positive rows** of weight $N$ to have a common prime divisor.

* Some full prime power $w$ in $P+1$ satisfies $3w^2\ge4(P+1)$. This gives infinitely many weights even with $P$ prime. Baker and Harman's Theorem 2 with shift $-1$ gives infinitely many primes $P$ for which the largest prime factor of $P+1$ exceeds $P^{0.677}$. [Baker and Harman, *Shifted primes without large prime factors*, p. 332](https://matwbn.icm.edu.pl/ksiazki/aa/aa83/aa8342.pdf)
* After removing all powers of 2 and 3 from $P+1$, the remaining integer is a prime power times a factor at most 13, or is 1.
* The odd part of $3P-1$ is at most **1,000,000**.
* $3P-1=DQ$, where $Q$ is a prime power and $D\le615$. Up to $D\le2400$ the conclusion still holds if either $P$ or $Q$ is a square.

The [leading-digit and exponent theorem](sparse_quadratic_exponents.md) settles all twelve remaining shapes through $P=10^{1000}$ and at unrestricted exponents when either prime base is at most one million.

The last two conclusions are proved by finite carry certificates uniform in the size of the prime powers. The first two have elementary arithmetic proofs. For $D\le2400$ the table below restricts possible failures even when both powers are nonsquares.

## The forced row

Write $H=3P-1$ and $V=3P+2$, so $N-1=(P+1)H$. The [endpoint reduction](all_weights_endpoint_reduction.md) assigns each prime of $N(N-1)$ to one of two rows avoiding it. A row avoids a prime if its three parts add without a carry in that base. The [five-unit weight theorem](sparse_weight_factors.md) gives the following deductions.

For $p=3$, the base-$p$ weight has three digit units and full $p$-part $P$. Its necessary inequality $27(N-1)^2\le P^3$ is impossible. For $p=2$ and $P\ge4$ there are three binary units and the full 2-part is $2P$, which also violates the inequality. At $P=2$ the weight is 16.

Take $p\ge5$. Its five units are $P,P,P^2,P^2,P^2$. The row $A$ avoiding $p$ has a singleton, its assigned weight factor is exactly $P$, and the companion is $C=V(y_1,y_2,y_3)$ with $\sum y_i=P$. If $S$ is the product of full predecessor powers assigned to $A$, the other row must satisfy
$$
                     27(N-1)^2\le P^3S^2.                         \tag{1}
$$
In particular $S>2$. Two singleton parts would force $S=1$, as both are units at every predecessor prime. Thus the other two parts must each contain two digit units. The singleton must be the residue-one position at every full power in $S$, so each of the other parts is divisible by $S$. A pair of equal units would give $S\mid2$ after cancelling a power of $p$. Both pairs are therefore mixed, giving the unique shape
$$
              A=(P^2,P(P+1),P(P+1)).
$$
The repeated part forces an even coefficient. Every prime of $H$ must be avoided by $C$, since $S\mid P+1$, $\gcd(H,P+1)\mid4$, and $A$ cannot take 2. We may use the earlier verified range $P\le200{,}000{,}000$ whenever needed. [Five-unit proof and finite verification](sparse_weight_factors.md)

## A restriction on every prime in $S$

If $r^e=w\parallel P+1$ for odd $r$, write $P=bw-1$. After dropping the first $e$ base-$r$ places, the entries of $A$ are
$$
             b^2w-2b,\quad b^2w-b,\quad b^2w-b.
$$
Only the first part had a nonzero discarded block, equal to 1. At each $m=r^j$ with $1\le j\le e$, put $t=b\bmod m$. The three displayed entries have residues $[-2t]_m,[-t]_m,[-t]_m$, where brackets denote least nonnegative residues. Their sum is below $m$ exactly when $4t>3m$. Indeed the repeated residues are $m-t$ and the first is $m-2t$ or $2m-2t$ according as $2t\le m$ or $2t>m$. In particular $r=3$ cannot be assigned to $A$, and every $w$ assigned to $A$ satisfies
$$
             b>3w/4,\qquad 3w^2<4(P+1).
$$
The residue equivalence is also machine checked. [Lean proof](SparseQuadraticLayers.lean)

Equation (1) gives $S^2>243(P+1)$, since $(P+1)(3P-1)^2-9P^3=3P^2-5P+1>0$. If a full power $w$ is omitted from $S$, then $S\le(P+1)/w$, so $w^2<(P+1)/243$. This proves the large-power criterion. If the part coprime to 6 is a prime power times a factor at most 13, taking the full power of that prime leaves an even smaller cofactor. Whether it is assigned or omitted, $S\le13\sqrt{4(P+1)/3}<\sqrt{243(P+1)}$. The same conclusion holds if the part coprime to 6 is 1.

## Cutting at a full power of $H$

Take $Q$ to be a full prime power in $H$, with $H=DQ$. At its base, two of the $y$'s vanish modulo $Q$, since $V$ is a unit. Write them $Qi,Qj$ with $1\le i\le j$. Then $3(i+j)<D$ by positivity and $3\nmid D$. Every full power in $D$ divides at least one of $i,j$, so $D\mid ij$. For each unshifted part we also have $D\mid i(3Qi-1)$, since the actual part $VQi$ is 0 or 1 modulo each full power and $V\equiv3\bmod D$.

**Odd $Q$.** Here $D$ is even. For $P>200{,}000{,}000$ and $D\le2400$, $Q>2^{12}$. The direction of $C$ differs coordinatewise by at most $1/(DQ+1)$ from
$$
                     s/D=(3i,3j,D-3(i+j))/D .
$$
An odd coefficient requires different binary orders for $i,j$ and disjoint binary digits throughout $C$. To test the first twelve columns, start with $z=0$. At each column allow any part-digit vector $u$ of weight at most one, with weight one at the first column. Update $z\gets2z+Du-s\sum u$ and keep $s_\ell-D\le z_\ell\le s_\ell$.

To see the bounds, divide the parts by the least power of two strictly above $N$. At depth $h$ their scaled prefixes have integer parts $v_\ell$ and tails $e_\ell$ summing to $\tau<1$, by disjointness. The coordinate error in $D(v_\ell+e_\ell)-s_\ell(\sum v+\tau)$ is at most $D2^h/(DQ+1)<1$. The tail contribution $s_\ell\tau-De_\ell$ is between $s_\ell-D$ and $s_\ell$. Rounding the prefix discrepancy, which is an integer, gives the stated bounds. [General direction argument](rational_proportion_divisors.md)

For a base $r\mid D$ and $m=r^j$, the residue $x=Q\bmod m$ is a unit. Write $v=Dx+3$. Carry freedom requires
$$
            [vxi]_m+[vxj]_m\le [(Dx+1)v/3]_m ,
$$
where 3 is inverted modulo $m$. The verifier lifts residues one digit at a time and rejects a shape if one base runs out of extensions. Every exclusion is checked independently by constructing positive integer rows and evaluating factorial floors at the layers. Truncating a search retains its shape.

Among 62,935 pairs with $D\mid ij$, 24,302 pass the binary order and prefix tests, and 24,290 of these have a finite local exclusion. The twelve remaining possibilities are listed below. In the final column every unit square for $Q$ violates the displayed inequality at that modulus.

| $D$ | $(i,j)$ | Square exclusion modulus |
|---:|---|---:|
| 616 | (77,128) | 8 |
| 872 | (32,109) | 8 |
| 1432 | (179,288) | 179 |
| 1528 | (191,288) | 8 |
| 1640 | (32,205) | 8 |
| 1672 | (32,209) | 8 |
| 1808 | (113,320) | 8 |
| 2072 | (64,259) | 7 |
| 2168 | (271,448) | 8 |
| 2224 | (128,417) or (192,417) | 8 |
| 2344 | (96,293) | 8 |

Every remaining $D$ is divisible by 8, forcing $P=(DQ+1)/3\equiv3\bmod8$. This excludes square $P$. The statements also allow a representation in which the displayed power is not full. Enlarging an odd $Q$ to the full power would multiply it by at least 5, since $3\nmid H$, and would leave a cofactor at most $2400/5<616$.

**Even $Q$.** Now $Q=2^e$ and $D$ is the odd part of $H$. If $Q>2D$, set
$$
 h=(D^2-1)/3,\qquad l=(Q+4D)/3.
$$
These are integers, $l<Q$, and $N=1+lQ+hQ^2$. The first two parts of $C$ are $3iQ+DiQ^2$ and $3jQ+DjQ^2$, with $3(i+j)<D$. Thus $Di,Dj$ must have disjoint bit positions inside $h$, and $3i,3j$ must have disjoint bit positions inside $l$. As $3l\equiv4D\bmod Q$, the latter inclusion only reads the low bits of $4D/3$ computed by modular inversion. These tests apply to every admissible $Q>2D$ simultaneously.

For all odd $D\le10^6$ coprime to 3 the bit conditions retain 15,433 single entries across the cofactors. The three pairs $(D,i,j)$ are $(20063,1728,3456)$, $(246973,56,69632)$ and $(450871,224,7168)$. Each violates $D\mid ij$.

For $Q\le2D$ enumerate each power of two with integral $P>200{,}000{,}000$. The 1,471,770 resulting scales include composite inputs $P$ as well. For every full power dividing $D$ choose independently whether $3Qi$ is zero or one. The Chinese remainder theorem fixes $i\bmod D$. An admissible $0<3i<D$ would make $Vi$ a binary submask of $(N-1)/Q=D(DQ+4)/3$. There are zero such single entries. Taking the full 2-part only decreases the odd cofactor bound if a given power of two was not full.

## Reproduction

The Python verifier uses modular bit lifting for the separated blocks. The independent C++ verifier selects an actual power of two representing each separated case and tests the integer products directly. For the overlapping cases C++ constructs $i$ through complementary divisors of $D$ and uses binary digit sums. The Python controls include 10,572 middle-layer comparisons and 1,424 exact trinomial parity comparisons. The Lean file checks the middle-layer equivalence and the block identity.

~~~sh
python3 check_sparse_quadratic_powers.py
mkdir -p .work
c++ -std=c++17 -O3 check_sparse_quadratic_blocks.cpp -o .work/check_sparse_quadratic_blocks
.work/check_sparse_quadratic_blocks
lean SparseQuadraticLayers.lean
~~~
