# Three rows at the weights $(2P+1)(5P+1)$

Let $P=p^k$ be a prime power with $k\ge1$, and put $H=10P+7$. Wasserman's three-row conjecture holds at $N=1+PH$ for every $P\le100{,}000{,}000$. The bound is on the input $P$ and allows weights up to $100{,}000{,}000{,}700{,}000{,}001$.

The following criterion has no upper bound on $P$. A failure requires a factorization $H=d_1d_2d_3$ into pairwise coprime integers $d_i>1$ with
$$
                        1/d_1+1/d_2+1/d_3<1/4 .
$$
In particular the conjecture holds if $H$ has at most two distinct prime factors, or has exactly three and one occurs to full power 3. Here “full power” means the entire power of that prime in $H$.

## Classifying the first row

By the [all-weight endpoint reduction](all_weights_endpoint_reduction.md), a failure contains $(1,1,N-2)$. Assign each prime of $N(N-1)$ to one of the other two rows avoiding it. Let $A$ take $p$, and let $C$ denote the other row. Their assigned products of full powers of $N$ are $U,V$ respectively. We have $UV=N$, $\gcd(U,V)=1$, both scales are at least 3, and each divides every part of its row. One of the two coefficients must be odd since $2\mid N(N-1)$. [Assignment proof](endpoint_reduction.md)

The [one-billion verification](full_range_1b.md) covers $P\le3981$. For larger $P$ the blocks $10P^2,7P,1$ are disjoint in base $p$. By Kummer, the two unshifted parts of $A$ have the form $P(aP+b),P(cP+e)$, with $a+c\le10$ and $b+e\le7$. Order the pairs lexicographically. Both lines are nonzero, and so is the remaining line with coefficients $(10-a-c,7-b-e)$, since otherwise a part would be zero or one. As $U$ is coprime to $P$, it divides both unshifted lines. Consequently $U$ divides
$$
                 u_* = \gcd(a^2-7ab+10b^2,c^2-7ce+10e^2,ae-bc).
$$
Indeed $a^2N$ reduces to the first expression modulo $U$ using $aP\equiv-b$; the second and third expressions follow in the same way.

Let $S$ be the product of the full powers in $H$ whose bases are assigned to $A$. Each has two zero positions modulo its power in $N-1$. For a base other than $p$, at least one unshifted line vanishes there, and $10P\equiv-7$ gives $10b-7a$ or $10e-7c$ divisible by its power. At $p$ the shifted part is 1 modulo $P$ and both unshifted parts are multiples of the full $p$-part of $PH$. Cancel $P$ to get the same conclusion. Thus $S$ divides $S_* = \operatorname{lcm}(|10b-7a|,|10e-7c|)$. Neither entry can vanish since $(10,7)$ is primitive and all three lines are nonzero. Two parts of $C/V$ are divisible by each full power in $H/S$, so $(H/S)^2$ divides their three-part product, which is at most $(U/3)^3$.

Of the 1,069 line pairs, seven have $u_*=0$. For each of the rest with $u_*\ge3$, $27(10P+7)^2\le u_*^3S_*^2$ implies $P\le3981$. The [short integer enumeration](check_seven_quadratic_shapes.py) checks these bounds. In the seven remaining cases both lines are multiples of the same linear factor, giving
$$
 T=1+tP,\quad W=1+uP,\quad N=TW,\qquad
 A=W(Pi,Pj,1+(t-i-j)P),
$$
where $\{t,u\}=\{2,5\}$ and $1\le i\le j$, $i+j\le t$. Let $B=\gcd(U,T)$. As $U\mid W\gcd(i,j)$ and every prime assigned to $U$ occurs there to its full exponent in $TW$, $B\mid\gcd(i,j)$. The only possible $B>1$ is 2 at $(i,j)=(2,2)$. The repeated parts make the coefficient even, contradicting the assignment of 2 to $A$. Hence $B=1$ and $T\mid V$.

When $t=5,j=3=t-u$, the first and third parts sum to $W^2$. Here $i=1$ or 2, so $\gcd(i,j)=1$ and $U\mid W$. If $\ell\mid U$ and $\ell^c\parallel W$, the sum of these two parts is zero modulo $\ell^{2c}$. They add without a carry in base $\ell$ since $A$ avoids $\ell$. Both parts must vanish modulo $\ell^{2c}$, so $\ell^c\mid i$. As $U\mid W$ this gives $U\mid i$, impossible for $U\ge3$.

Now $i,j$ are powers of two, including 1. An odd prime of $H$ assigned to $A$ would have to divide $i$ or $j$: the two zero positions meet an unshifted part and $P,W$ are units there except at $p$, where we first cancel $P$. Thus every prime of $H$ is assigned to $C$. In particular $p\mid H$ is impossible in a failure.

## Reconstructing the companion

Write $C=T(f_1,f_2,f_3)$ with positive integers of sum $W<H$. Assign each full power of $H$ to the position at which the actual part is 1 modulo that power. The three group products $d_i$ are pairwise coprime. Each exceeds 1 since a part in an empty group would be a multiple of $H$ even after dividing out $T$. Since $TW\equiv1\pmod H$, the parts are uniquely determined by the least nonnegative residues
$$
            f_i=(H/d_i)\,[W(H/d_i)^{-1}]_{d_i}.
$$
For $t=2$, $A$ has $i=j=1$ so $C$ must be odd. Also $P$ cannot be even in this case, since $A$ is assigned $p$. For odd $P$, $T$ is odd and $W$ even. Carry freedom forces each $f_i$ to be even since $N$ is even. Thus $W\ge2H\sum_i1/d_i$ and $2W<H$. For $t=5$, $5W<H$ directly. Both cases give the stated reciprocal bound. When $t=5$ and $P$ is odd, $T$ is even and divides $V$, so $C$ must again be odd. An even prime power is a power of two.

The first search sieves $H$ into prime powers and enumerates all partitions into three nonempty unlabeled groups. It covers every positive integer $P$ through the bound for $t=2$, and every odd integer and every power of two through the bound for $t=5$. Where an odd coefficient is required it tests the three pairwise intersections of the binary digits of the actual parts.

| $t$ | Inputs with at least three prime factors in $H$ | Partitions | Positive rows of the required sum | Rows surviving parity |
|---|---:|---:|---:|---:|
| 2 | 52,964,235 | 225,926,579 | 25,059,885 | 0 |
| 5 | 26,482,458 | 112,962,943 | 1,429,771 | 0 |

The second search chooses the smallest group $d$ first. It satisfies $d^3\le H$. Its part is $T(H/d)n$, where $\gcd(n,d)=1$ and $0<tn<d$ since $tW=H-u$. The congruence $tn H/d\equiv-u\pmod d$ fixes $P$ modulo $d^2$ since $d$ is coprime to 10. After testing the binary cut between that part and its complement, the program tries every ordinary divisor pair $ab=H/d$ with $d\le a\le b$ and $\gcd(a,b)=1$. It recovers one more part modulo $H$, obtains the last by subtraction and tests the remaining binary intersection. The searches use distinct enumeration orders and implementations of modular inversion. A third small check tests $H\mid(Tf)(Tf-1)$ directly for every positive $f<W-1$ through $P=2500$. It finds 58 rows at $t=2$ and one at $t=5$.

| $t$ | Smallest-group cuts tested | Surviving first cut | Surviving complete rows |
|---|---:|---:|---:|
| 2 | 93,858,743 | 73 | 0 |
| 5 | 15,977,839 | 14 | 0 |

Run from the repository root with Python 3.9 or later and a C++17 compiler.

```sh
python3 check_seven_quadratic_shapes.py
mkdir -p .work
c++ -std=c++17 -O3 check_seven_quadratic_row.cpp -o .work/check_seven_quadratic_row
c++ -std=c++17 -O3 check_seven_quadratic_cut.cpp -o .work/check_seven_quadratic_cut
.work/check_seven_quadratic_row 100000000
.work/check_seven_quadratic_cut 100000000
```
