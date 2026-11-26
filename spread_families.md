# A lower bound on spread and weights above three times a prime power

**Every nonconstant nonendpoint row of a counterexample has spread greater than 1,000,000.** Spread means largest part minus smallest part. A row with three equal parts remains possible under this theorem. An independent analytic bound gives $w^5>27(N-1)^2$ for every positive such spread $w$. [Gap bounds and proof](row_spreads.md)

Write $M=N-1$. For every prime $p\ge5$, power $Q=p^k$ with $k\ge1$, and $0<c<Q$ with $p\mid M$ and $N=3Q+c$, the **full conjecture holds** under either condition below.

- $c\le1{,}000{,}000$ or $c^5\le27 M^2$.
- All digits of $c$ in base $p$ are 0 or 1. If $P$ is its largest digit unit and $D=c-P$, then $D^4c^5\le27M^4$.

The second condition proves, for example, all $2^m$ weights
$$N=3p^{16m}+p^{12m}+1+\sum_{i\in I}p^i,\qquad I\subseteq\{1,\ldots,m\},\quad m\ge1.$$
The gap from $3p^{16m}$ has order $N^{3/4}$. The number of lower digit units is unrestricted. Below we also give an optimal four-prime list for one member of this family.

## Proof for the prime-power weights

A counterexample contains the [endpoint row](all_weights_endpoint_reduction.md), whose coefficient is divisible by $p$. Some other row $A$ must therefore avoid $p$. Its parts distribute the three upper digit units $Q$ and the digit units of $c$ with no carries. In the assignment $UV=N$, $U\ge3$ divides every part of $A$ and $\gcd(U,p)=1$ because $p\mid M$. This excludes a part exactly equal to $Q$.

Under the first condition every part of $A$ exceeds $c$ by the [individual-part bounds](million_part_bound.md). Thus each part receives one $Q$. The spread is at most $c$, so the spread bounds force all three parts to agree. Their three residues cannot sum to 1 at the units digit of $N$ in base $p$. This proves the first case. The same case also holds when $p\mid N$ and $3\nmid N$, since an equal row then has nonintegral parts.

For the second condition, the three lower portions in $A$ have disjoint base-$p$ digit supports. Two equal positive lower portions, or two in the ratio $1:2$, would share their lowest occupied digit position since $p$ is odd. Also $cD\le27^{1/4}M/c^{1/4}< 9M/2$. There are three ways to distribute the upper units up to permutation.

- With all three in one part, the other two have product at most $Dc$ since at least one misses $P$. But (4) of the [gap bounds](row_spreads.md#another-obstruction-from-the-other-two-parts) gives product at least $9M$ for every two parts.
- With two in one part and one in another, write the row as $(a,Q+b,2Q+b')$. If $P$ is in $a$ then $|2b-b'|\le2D$. Otherwise $a\le D$. Thus $a|2b-b'|\le2cD<9M$. By the same obstruction this forces $2(Q+b)=2Q+b'$. Disjoint supports then give $b=b'=0$, producing the excluded part $Q$.
- With one upper unit in each part, all three lower portions must be positive to avoid a pure $Q$. They are distinct, and the two missing $P$ have sum at most $D$. In increasing order the first gap is at most $D$, the other gap is less than the spread, and the spread is at most $c$. This contradicts $w(df)^4\ge27M^4$ under the stated inequality.

For the displayed example, $D\le(5/4)p^m$ and $c\le2p^{12m}$, while $M\ge3p^{16m}$. These inequalities imply the second condition at every stated exponent.

## Four primes for a 34-digit weight

Take $p=11$, $m=2$, $I=\{2\}$. At
$$N=6334133033455498532444456857341726$$
an optimal certificate consists of the four primes
$$11,\qquad q=3604276009509249793981575479,\qquad
s=72204423293878581162091272241,\qquad r=N-109.$$
Every triple of rows has a common prime from this list. Indeed the 30 unordered positive rows avoiding 11 are obtained by allocating the three copies of $11^{32}$ and the units $11^{24},11^2,1$. Exact factorial valuations show every one carries at $q,s,r$. The two factorizations are
$$N=2\cdot 3^2\cdot89\cdot1097\cdot q,\qquad
N-1=5^2 11^2\cdot29\cdot s.$$
A row avoiding $q$ has all parts multiples of $q$. Avoiding $s$ requires two parts multiples of $s$. The two conditions together require weight at least $q(2s+1)>N$. A row avoiding $r>N/2$ has two parts of sum at most 109, incompatible with avoidance of either $q$ or $s$. Thus the four avoided-prime row sets are pairwise disjoint.

Every universal prime would divide the endpoint coefficient. For each prime in the two factorizations the digit sum of $N$ in that base is at least three, so the digit units can be distributed to three positive parts without a carry. Hence there is no universal prime and three rows can avoid any prescribed list of at most three primes between them. The verifier checks the factorizations, digit sums, the 30 rows and recursive integer primality certificates for the listed factors and $r$.

## Exhaustion at bounded spread

The full result is [already verified through one billion](full_range_1b.md). Above that threshold put $W=10^6$, and suppose $A$ has positive width at most $W$. A binary carry-free row has width greater than $N/4$, so the companion must be assigned 2. Thus $U$ and $S$ are odd. The normalized companion parts of sum $U$ have three distinct 2-adic valuations because multiplying by their common scale $V$ preserves the differences of valuations. Hence $U\ge1+2+4=7$.

Write $A=(a,a+d,a+d+f)$ with $w=d+f\le W$. Allow one gap to vanish. If $g=\gcd(d,f)$, then $U\mid g$. For the predecessor bins $S_i$, we have $S_i\mid B_i$ with
$$B_1=\gcd(f,d+1),\quad B_2=\gcd(w,d-1),\quad B_3=\gcd(d,f-1). \tag{1}$$
This follows from the three possible residue vectors $(1,0,0),(0,1,0),(0,0,1)$. The $B_i$ are pairwise coprime since a prime cannot give two different such gap patterns. They are positive here since $g\ge7$. Put $S=\prod S_i$. The remaining possibilities satisfy
$$T\le\sqrt{U^3/27},\quad ST\equiv-1\bmod U,\quad N=ST+1,\quad
a=(N-2d-f)/3.$$
Require integral positive parts divisible by $U$, coprime weight scales and coprime predecessor products. Then check all three actual residues at each $S_i$ and carry freedom at every prime of $U S$.

Two enumerations generate these assignments.

1. Fix $g$ and use $d=gr$, $f=g(h-r)$ with $0<r<h$ and $\gcd(r,h)=1$ in the distinct case. The repeated cases have $d=0$ or $f=0$. For speed, $B_1B_2,B_2B_3\le2w+1$ and $B_1B_3\le w+1$ since these products divide $d+2f+1$, $2d+f-1$, and $f-d-1$, respectively, in absolute value. None vanishes for $g\ge7$. Also $B_1B_2B_3\le h^3/4$ for distinct parts by coprimality with $g$ and $r(h-r)h\le h^3/4$. For $T_{\max}=\lfloor\sqrt{\operatorname{odd}(g)^3/27}\rfloor$, at least one odd part of a $B_i$ is at least $(10^9/T_{\max})^{1/3}$, where $\operatorname{odd}(g)$ removes all factors of 2 from $g$. Generate the largest such bin from the divisors of $w+1,h,w-1$, respectively. It fixes $r$ modulo that bin. Recompute (1) to resolve ties and proper divisors. Then enumerate $U\mid g$ and the three bin divisors.
2. Independently choose $w$ and its odd divisor $U\ge7$. Take odd $S_1\mid w+1$, $S_2\mid w/U$ and $S_3\mid w-1$. With $x=d/U$ the requirements $x\equiv0\bmod S_3$, $Ux\equiv-1\bmod S_1$ and $Ux\equiv1\bmod S_2$ give one class modulo $S$. Enumerate every member in $[0,w/U]$.

Both use integer bounds, including square and cube roots corrected by exact products. The table counts assignments after the scale and integrality tests and before checking all the assigned row digits.

| Row $A$ | Assignments | Passing its carry tests | Largest candidate weight |
|---|---:|---:|---:|
| Two equal parts | 426,608,408 | 1,711,987 | 95,999,713,000,191 |
| Three distinct parts | 19,113 | 0 | 246,213,009,120 |

The sorted files of passing assignments from the two searches have the same SHA-256, `01203109ded8b19f978e2ef9022e579f53426f4933d23675d522c54a8e6bcc76`. A separate program verifies the valuations on every emitted row. Of these assignments 1,136,662 pass the test $2q+1\le U$ at every full prime power $q$ of $T$. The companion search distributes these powers to exceptional positions, lifts the resulting residues to a positive composition of $U$, and tests all assigned primes.

A second companion check works one part at a time. Each normalized part $y$ is either 0 or $U$ at every full power of $T$, and $Vy$ can occupy only bits of $N$. Only 63 assignments admitted even one such putative smallest part with $0<y\le U/3$. None gave a binary-free row with the required predecessor residues after testing the middle part and reconstructing the third. Thus every passing assignment is excluded.

## Reproduction

Build with a C++17 compiler supporting 128-bit integers. The Python audit and four-prime check use python-flint. The emitted files can be regenerated locally.

~~~sh
mkdir -p .work
c++ -std=c++17 -O3 check_row_spreads.cpp -o .work/row_spreads
.work/row_spreads 1000000 1000000000 --by-span --emit > .work/spread-span.txt
.work/row_spreads 1000000 1000000000 --emit > .work/spread-primitive.txt
LC_ALL=C sort .work/spread-span.txt | shasum -a 256
LC_ALL=C sort .work/spread-primitive.txt | shasum -a 256
python3 check_row_spread_companions.py .work/spread-span.txt
.work/row_spreads 1000 4000 --by-span --emit > .work/spread-control.txt
python3 check_row_spread_companions.py .work/spread-control.txt --control
python3 check_spread_bounds.py
~~~

The control enumerates weight factors and actual normalized rows directly at $35{,}000\le N\le55{,}000$ and matches all six passing first rows of spread at most 1,000. For the displayed prime-power family the last script tests every carry-free allocation at $p=5,7$, $1\le m\le5$ and every subset $I$. The 211,044 positive labelled allocations at 124 weights are all excluded by the row gcd, the bound (4) in the [proof](row_spreads.md), or the distinct-gap inequality.
