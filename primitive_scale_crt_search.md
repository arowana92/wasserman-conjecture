# A search for rows that avoid a large product of predecessor primes

The conjecture holds for **all choices of three rows** at each of 7,364 weights identified below. Of these, 6,690 exceed one billion. The largest is $10{,}940{,}257{,}758{,}976$. These weights are selected by a carry condition and do not form a consecutive range.

Write a positive row as $Vy$ with $y$ primitive, $u=\sum y_i$ and $N=uV$. Take powers $d,e,f$ of three distinct primes and put $t=def$. We search for rows satisfying
$$
  2\le V\le2000,\qquad \max(d,e,f)\le2000,\qquad
  t\mid N-1,\qquad t>u,\qquad \gcd(K_{Vy},t)=1.
$$
The powers can be proper divisors of their full powers in $N-1$. There is no other restriction on $u$ or $N$. Such rows concentrate large missing factors in one coefficient and are relevant to the [endpoint assignment](endpoint_reduction.md).

## Enumerating the row

At each power the actual row has residues $1,0,0$ in some order. Every primitive coordinate must be the exceptional coordinate at one power, since one vanishing at all three would be divisible by $t$ and at least $t>u$. Thus there is only one candidate up to ordering for any $(d,e,f,V)$ with $\gcd(V,t)=1$. The coordinate exceptional at $d$ is
$$
             ef\,((Vef)^{-1}\bmod d),
$$
and cyclically at $e,f$, using least positive residues. Check the sum, primitivity, and carries at all three prime bases.

| Stage | Count |
|---|---:|
| Unordered triples of powers with an admissible unit scale | 11,118,969,242 |
| Primitive candidate rows with sum less than $t$ | 1,548,724,646 |
| Rows passing all three carry tests, counted with their powers | 7,370 |
| Such records using a nonprime power | 16 |
| Distinct weights | 7,364 |
| Records with $\gcd(N(N-1),K)<N$ | 11 |

For an independent comparison through power 257 and scale 8, start at $u=V^{-1}\bmod t$. Multiply $u$ by each of the three ordinary CRT idempotents modulo $t$, check that the sum of positive residues is $u$ and evaluate the three factorial valuations. All four records in this smaller box match. The Python verifier also recomputes the required valuations of every record in the full output.

## Checking all rows at these weights

Factor $N$ and $N-1$ into full prime powers, with primality checks and exact recomposition. By the [all-weight endpoint reduction](all_weights_endpoint_reduction.md), it suffices to rule out a pair of rows whose coefficients have no common prime dividing $N(N-1)$. Enumerate all allocations $UV=N$, $ST=N-1$ of the full powers, assigning $U,S$ to one avoiding row. Name the sides so that $S\ge V$. Such a side exists since $(U-1)(V-1)<N-1$. Its normalized parts sum to $V$, two are divisible by each assigned predecessor power, and each receives an exceptional residue. Thus $S$ has at least three prime bases, $2w+1\le V$ for each of its powers, and $27S^2\le V^3$. The CRT rows for both sides and the carry tests are described in the [assignment proof](endpoint_reduction.md). This checks all row choices even when the originally identified row is absent.

There are 692,169 assignments after the elementary tests. Exactly 11 admit a first row and none admit a row on the complementary side. Additionally, for each of the 11 recorded rows with endpoint gcd below $N$, the verifier independently enumerates positive partitions for its required companion, using factorial floors and no CRT. Every prime divisor of $N$ present in the recorded coefficient forces its **full weight power** into each part of a companion. After division by this product the remaining sum is at most $V$, since $u$ divides the recorded coefficient. This makes direct verification short even at the largest weights.

## A small endpoint gcd

At $N=4{,}217{,}316{,}279$ the row $(180377267,1044164635,2992774377)$ has endpoint gcd $522{,}499{,}362$. The coefficient valuations at the complete endpoint factorization are as follows.

| Integer | Prime factors, in order | Coefficient valuations |
|---|---|---|
| $N$ | $3,11,113,1130951$ | $16,10,0,1$ |
| $N-1$ | $2,7,163,1321,1399$ | $22,9,0,0,0$ |

Both endpoint factors are squarefree. The product of the avoided endpoint primes is $34{,}039{,}767{,}101$, more than eight times the weight. A companion avoiding the five present primes must have each part divisible by $3\cdot11\cdot1130951=37321383$. Its three normalized parts sum to 113, and none of their partitions gives an avoiding row. Thus this example disproves the single-row bound $\gcd(N(N-1),K)\ge N$ at this weight while satisfying Wasserman's conjecture for every selection.

## Reproduction

Use C++17, Python and python-flint. The search emits rows in deterministic order. The SHA-256 of the complete row file is `f489a4d31ca1a72e91147c391c94e1ed578775b08680f32091df76f213c66aef`.

~~~sh
mkdir -p .work
c++ -O3 -std=c++17 check_three_factors.cpp -o .work/three_factor_scales
.work/three_factor_scales 2000 --primitive-scales 2000 > .work/power_scale2000.rows
python3 check_heavy_predecessor_crt.py .work/power_scale2000.rows --cross
~~~

[C++ enumeration](check_three_factors.cpp) · [Independent valuations, companion and full-weight checks](check_heavy_predecessor_crt.py)
