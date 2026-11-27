# Two freely placed prime powers and an arbitrary shifted coefficient

Let $p$ be a prime, $H,L,Q$ powers of $p$, $H\ge L\ge p$, $Q\ge1$, and $m$ a positive integer. Set
$$N=H+mQ+L.$$
Assume $p\mid N(N-1)$. For $Q\ge p$ this is automatic. For $Q=1$ it means $m\equiv0,1\bmod p$.

**Cutoff theorem.** If $27(N-1)^2>(Lm)^5$, Wasserman's three-row conjecture holds at $N$. The exponents of $H$ and $Q$ may occur in either order. The base-$p$ digits of the three summands may overlap.

**Quantitative form.** For $Q\ge p$, take any positive trinomial row $A$ of sum $N$ whose coefficient $K_A$ is prime to $p$, and any other positive row $B$ of the same sum. Then
$$\gcd\big(N(N-1),K_A,K_B\big)^3(Lm)^5\ \ge\ 27(N-1)^2.$$
For example, at $p=2$, $L=2$, $m=65535$, $H=2^{1000}$ and $Q=2^{400}$ this forces the endpoint gcd above $2^{638}$ for every such pair. There are odd coefficients here since the binary digits of the three summands are disjoint and can be assigned to three nonempty parts.

**All-exponent computation.** The conjecture holds under the stated hypotheses whenever $Lm\le131072$. This covers 453,582 distinct triples $(p,L,m)$. In particular $2^a+m2^b+2$ works for every $a\ge1$, every $b\ge0$ and every $1\le m\le65536$.

## The assigned prime powers

By the [endpoint reduction](all_weights_endpoint_reduction.md), it suffices to exclude a pair $A,B$ with $\gcd(N(N-1),K_A,K_B)=1$. In such a pair one row avoids $p$, and we designate that row $A$. For the quantitative result start with the designated $p$-avoiding row. Put $D=\gcd(N,K_A,K_B)$, $J=\gcd(N-1,K_A,K_B)$. Use the [partial prime-power assignment](endpoint_pair_bound.md), assigning the power at $p$ to $A$, and call its two assigned products $U\mid N$ and $S\mid N-1$. The other row gives
$$27(N-1)^2\le D^3J^2 U^3S^2. \tag{1}$$
All parts of $A$ are multiples of $U$. Each prime power exactly dividing $S$ has one residue 1 and two residues 0 in $A$. When $DJ=1$, the parts of $B/(N/U)$ are positive integers with sum $U$, so $U\ge3$.

Kummer's carry criterion lets us distribute the single units $H,L$ and the digits of $m$ among the three parts of $A$. When original digit units overlap, merge $p$ equal units into a larger one until the digits of $N$ remain. Give the original units underlying each final unit to its recipient in $A$. At least one part is pure, of the form $z=kQ$ with $1\le k\le m$. Its digits fit inside $N$ and the digits of $k$ fit inside $m$.

When $Q\le L$ use $U\mid z$ and $S\mid (z/U)(z-1)$, which follows from the residue condition and $\gcd(S,U)=1$. If $z>1$, then
$$U^3S^2\le U z^2(z-1)^2<z^5\le (Lm)^5. \tag{2}$$
For the qualitative theorem $z\ge U\ge3$. The quantitative theorem has $z\ge Q\ge p$.

## Eliminating the free displacement

Suppose $Q>L$ and either $H>L$ or $p$ is odd. The full power of $p$ in $N$ is $L$ (when $H=L$ the low coefficient is 2). Thus $U=Lu$ with $p\nmid u$ and $u\mid k$. Every factor of $S$ is prime to $p$.

First suppose $H,L$ go to different parts. The other part containing $L$ is $lQ+L$ with $0\le l\le m-k$. The product of the powers of $S$ for which $kQ$ is zero divides $k/u$. At every remaining power $w$,
$$kQ\equiv1,\qquad lQ+L\equiv0\pmod w.$$
Multiplying the latter by $k$ and using the former shows that $w\mid l+Lk$. These powers are pairwise coprime. If the two singletons go to the same part, write the two pure parts as $kQ,lQ$. The powers at which $kQ$ is one now divide $l/u$. Either placement gives
$$uS\le k\big(m+(L-1)k\big),\qquad
 U^3S^2\le L^3u k^2\big(m+(L-1)k\big)^2\le (Lm)^5. \tag{3}$$
The central elimination and the final numerical implication in the split case are also [checked in Lean](TwoSingletonPowers.lean).

The remaining case has $p=2$, $H=L<Q$. Merge the identical powers to the single unit $2L$ before distributing. Two slots are pure, $kQ,lQ$ with $k,l>0$ and $k+l\le m$. If $Q>2L$, then $U=2Lu$ with $u$ odd dividing $k,l$. Each power of the odd number $S$ divides one of $k/u,l/u$. Consequently
$$ U^3S^2 \le 8L^3 (kl)^2/u\le L^3m^4/2\le(Lm)^5. \tag{4}$$
If $Q=2L$, write $g=\gcd(k,l)$. Then $U\le Qg$ and $S\mid\operatorname{lcm}(k,l)=kl/g$. This gives $U^3S^2\le8L^3g(kl)^2\le L^3m^5/4$. This proves both theorems using (1) and $D^3J^2\le(DJ)^3$.

## Finite verification

For each $(p,L,m)$ in the range $Lm\le131072$, the strict cutoff leaves finitely many pairs $(H,Q)$. The [contiguous verification](full_range_1b.md) handles $N\le10^9$. For every remaining parameter choice the C++ program generates the possible digit subsets $k$ in decreasing order. It skips a part unless $U_0=\gcd(N,kQ)\ge3$, $N-kQ\ge6$, and
$$U_0^3\,\gcd\big(N-1,kQ(kQ-1)\big)^2\ge27(N-1)^2.$$
This is necessary by (1) at $DJ=1$. It also applies (2), (3), or (4) where applicable. In (3) the largest possible $u$ for a particular part is $\gcd(N/L,k)$. The search stops descending through $k$ only when a bound monotone in $k$ fails even with $u=k$. All powers and products used in the comparisons are computed exactly.

| Base and ordering | Parameter choices above one billion and within cutoff | Choices with an eligible part | Distinct surviving weights | Greatest weight within cutoff tested above one billion |
|---|---:|---:|---:|---:|
| Binary, $Q\le H$ | 31,853,223 | 180 | 113 | 1,168,231,235,584 |
| Binary, $Q>H$ | 15,870,318 | 147 | 117 | 1,168,231,235,584 |
| Odd prime, either order | 14,895,878 | 15,104 | 15,039 | 1,195,442,245,376 |

There are 131,071 binary families and 322,511 odd-base families. Across the 62,619,419 parameter choices the union of the survivor lists contains 15,255 distinct weights. Of the complementary assignments, 959,723 pass the factor and size filters. None supplies a heavy row avoiding all its assigned primes. The heavy row is the one whose assigned predecessor product is at least the other row's assigned weight product, as in the [endpoint enumeration](endpoint_reduction.md).

The Python check factors both endpoints of each surviving weight, verifies the prime factors and their products, and enumerates complementary prime-power assignments using the separate [CRT and Kummer verifier](check_binary_gap_finite.py). It also tests the quantitative and assigned-product inequalities using exact trinomial coefficients at weights through 48. These direct tests cover 53,570 choices, including 4,721 with partially retained powers and 452 with coincident binary singletons above which $Q$ is placed.

A direct Python enumeration over individual integers $k$ independently reproduced all three candidate lists restricted to $L\ge1024$. It tested 276,476 parameter choices and retained exactly 33, 47 and 14,626 entries in the three rows of the table, respectively. At $Lm\le8500$ the same method tested 30,708 choices and retained none.

~~~sh
mkdir -p .work
c++ -O3 -std=c++17 check_two_singletons.cpp -o .work/check-two-singletons
for mode in binary-ordered binary-above odd
do
  .work/check-two-singletons 131072 "$mode" > ".work/two-singleton-$mode.rows"
done
python3 check_two_singletons.py --cross-cap 131072 --cross-min-low 1024 .work/two-singleton-{binary-ordered,binary-above,odd}.rows
lean TwoSingletonPowers.lean
~~~

The final finite check uses Python with python-flint. Lean was run with version 4.36.0.
