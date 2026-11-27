# Common primes across arithmetic progressions of binary exponents

Let $a_j=b+jh$ with $b,h\ge1$. In **every 28 terms** of this progression there is an exponent $a_j$ such that Wasserman's conjecture holds at $1+d2^{a_j}$ **simultaneously for every odd $d\le65535$**. The following upper bounds count indices at which even one of the 32,768 cofactors could fail. The starting exponent and the step have no upper bound.

| Number of terms | Possibly exceptional indices |
|---:|---:|
| 28 | 27 |
| 50 | 31 |
| 100 | 42 |
| 1,000 | 81 |
| 1,000,000 | 879 |
| 1,000,000,000 | 10,234 |

For a **single** odd cofactor in the same range, some term in every **four-term** exponent progression satisfies the full conjecture. **Three terms suffice** at every cofactor except possibly $33661,56279,60829,61423$. In particular at the first unresolved cofactor, 12031, its exceptional exponents cannot contain a three-term arithmetic progression **or a three-term geometric progression**, including a nonintegral positive ratio giving distinct integer exponents.

## At the first differing digit

Use the previously enumerated necessary binary rows $A=(xq,yq,zq+1)$, with $q=2^a$, and the forced carry base $p$ associated with each. A counterexample must have $a\equiv a_0\pmod H$. In every case $2^H\equiv1\pmod p$, so the input's units digit in base $p$ is fixed. At 56279 there are two rows, with $a\equiv717$ or 801 modulo 1743. In a progression of at least three terms consisting entirely of these two residues the residue is constant, since $2(801-717)\not\equiv0\pmod{1743}$. Thus it suffices to consider one row at a time. [Bounded-cofactor theorem and row table](binary_cofactor_rays.md)

Suppose $m$ successive terms use a single row. Then $H\mid h$. Put $t=v_p(2^h-1)\ge1$ and $2^h=1+p^t u$ with $p\nmid u$. All $q_i=2^{b+ih}$ have the same low $t$ digits and reach the same carry state there. A state is called live if some finite sequence of later input digits gives a carry-free completed row. For $t\ge k$, reduce the geometric sequence modulo $p^{t+k}$. After removing its common low digits the next $k$ digits form an arithmetic progression modulo $p^k$ with unit difference $q_0u$. This follows by expanding $(1+p^tu)^i$ modulo $p^{t+k}$.

The verifier finds a $k$ for which **no** reachable live state permits $m$ such words. For $t<k$ it enumerates the common live prefixes of length $t$ and checks the actual geometric progression. At each new layer the first two words determine every subsequent one by their unit ratio. The row table below records the number of prefixes tested over $1\le t<k$ and the maximum additional number of digits required there.

| $d$ | initial $a_0$ | $p$ | $m$ | $k$ | shorter common prefixes | extra digits |
|---:|---:|---:|---:|---:|---:|---:|
| 12031 | 907 | 3 | 3 | 1 | 0 | 0 |
| 15871 | 6002 | 5 | 3 | 12 | 243 | 9 |
| 33661 | 410 | 5 | 4 | 1 | 0 | 0 |
| 56279 | 801 | 7 | 4 | 13 | 1054 | 5 |
| 56279 | 717 | 7 | 4 | 8 | 128 | 5 |
| 60829 | 35049 | 7 | 4 | 5 | 7 | 1 |
| 61423 | 7030 | 11 | 4 | 8 | 314 | 4 |
| 62463 | 5401 | 3 | 3 | 5 | 7 | 1 |
| 62809 | 15529 | 3 | 3 | 4 | 3 | 1 |
| 65095 | 312 | 5 | 3 | 12 | 78 | 4 |
| 65243 | 14788 | 3 | 3 | 1 | 0 | 0 |

At $d=12031$ every one of the 765 reachable live states after the units digit permits at most two next digits. More generally, three exceptional exponents $a,b,c$ cannot have equal values of $v_3(a-b),v_3(b-c),v_3(a-c)$. Since the exponents have the same parity, lifting the exponent gives $v_3(2^a-2^b)=1+v_3(a-b)$, and likewise for the other pairs. Equal valuations would require three different digits immediately after a common prefix. The exponent refinement below gives $a\equiv1\pmod3$. If $a,ar,ar^2$ were three distinct exceptional integer exponents with positive rational $r$, then $r\equiv1\pmod3$ as a 3-adic unit. The three differences have equal valuation since $3\nmid r+1$. This proves the geometric assertion.

## Counting indices for all cofactors together

First strengthen each necessary exponent progression. At each forced prime use complete power orbits through the largest prime power at most 65536, with at least one layer beyond any assigned weight power. In the tabulated carry base require the prefix to reach a live state. In the other bases test the factorial floors. If the retained indices in a full orbit of length $n$ are $i_1,\ldots,i_s$, they all lie in the class of $i_1$ modulo $\gcd(n,i_2-i_1,\ldots,i_s-i_1)$. Use this to refine the exponent class before the next prime. At 12031 the result is $a\equiv769251571\pmod{2043279108}$. In particular $19\mid N$, since $18$ divides the period and $12031\cdot2^{769251571}+1\equiv0\pmod{19}$.

For one row, restrict $b+jh$ to its refined class modulo $H$. Its indices form a progression of step $H/\gcd(h,H)$. On this progression write the exponents $e_0+iD$ with $H\mid D$. The scalars $2^{e_0+iD}$ have a common prefix of $t=v_p(2^D-1)$ digits. At the next digit the $p$ index classes modulo $p$ occupy different digits, because the valuation of the difference of the two scalars is exactly $t+v_p(i-j)$. The same statement applies inside each child, replacing $D$ by $pD$. This handles steps divisible by an arbitrarily high power of $p$.

From a live state $s$ let $F_s(n)$ bound the number retained among $n$ consecutive indices under arbitrary digit permutations at these splits. Put $F_s(0)=0$, $F_s(1)=1$. For $n=mp+r$ each child gets $m$ or $m+1$ indices, and exactly $r$ of the $p$ digits get the larger count. Sum $F_v(m)$ over live digit edges $s\to v$, then add the $r$ largest gains $F_v(m+1)-F_v(m)$, including a zero for each missing edge. The gains are zero or one. Indeed the resulting maximum over all tree permutations counts points in the first $n$ indices, and extending that interval by one changes any one permutation's count by at most one. Equivalently, the same property follows inductively from the recurrence, including at the boundary $r=p-1$. The maximum over the reachable live states bounds the initial packet.

For a subset $I$ of the eleven rows that actually occur, write $r_i,H_i$ for their refined progressions. Compatibility requires
$$\gcd(h,H_i,H_j)\mid r_i-r_j\qquad(i,j\in I).$$
Starting with $G=\operatorname{lcm}_{i\in I}H_i$, for each pair decrease every prime exponent violating this condition to the largest allowed exponent. Then $\gcd(h,H_i)\mid\gcd(G,H_i)$ simultaneously for all $i$. In $K$ input terms the packet for row $i$ therefore has length at most $\lceil K\gcd(G,H_i)/H_i\rceil$. Sum the corresponding $F$ bounds and maximize over the 2,047 nonempty subsets. This gives the first table. It counts a possibly bad index even when further prime tests would exclude it.

Run `python3 check_binary_exponent_progressions.py` for both tables. Run `python3 check_binary_exponent_progressions_cross.py` for an independent verification of all eleven progression exclusions. The second program starts with the partition of $[0,1)$ by the rational breakpoints of the three slope floors, allowing a one-unit offset on the last coordinate. It builds edges using factorial-floor equality and searches backward from the zero state. It tests words by lifting a starting residue and a difference, or a starting residue and a geometric ratio. Both constructions obtain the same reachable state counts.
