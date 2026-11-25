# The remaining obstruction at three-bit weights

**Theorem.** Write $N=2^r+2^s+2^t$, $0\le r<s<t$, and set $P=2^r$, $m=N/P$. For all sufficiently large $N$, any counterexample at $N$ has even weight and consists of the three rows
$$E=(1,1,N-2),\qquad A=(2^r,2^s,2^t),\qquad C=(my_1,my_2,my_3),\qquad \sum y_i=P.$$
The only remaining requirement on $C$ is that its coefficient avoid every odd prime of $N(N-1)$. Conversely any such $C$ at an even three-bit weight gives an actual counterexample with $E,A$. In particular **the full conjecture holds** at sufficiently large three-bit weights whenever $27(N-1)^2>P^3$. The threshold here is uniform over the exponents.

The [odd three-bit theorem](three_binary_digits.md) already treats odd $N$. For even $N$, the unique row avoiding 2 is $A$, since the three binary digit units must go to three different parts. At any odd prime dividing $N(N-1)$ its three residues modulo the prime are nonzero and cannot sum to 0 or 1. Thus $A$ avoids none of these odd primes. If $E$ occurs, its coefficient is $N(N-1)$, giving the stated equivalence. The primes of $m$ then force all parts of $C$ to be multiples of $m$, and avoiding $N-1$ implies $(N-1)^2\mid y_1y_2y_3\le(P/3)^3$. [Endpoint size calculation](even_endpoint.md)

## Excluding a counterexample without the endpoint

Assume $r\ge1$ and $E$ is absent. Call the two rows other than $A$ the ordinary rows. Write each row in increasing order. Their second parts tend uniformly to infinity by the [height 6 bound](power_intervals.md). No two rows can both have a part at most 2, since any full prime power of $N$ avoided by either row divides one of those two parts. Their product is at most 2, which would leave the last row with three positive normalized parts of sum at most 2.

We use the following form of the [weighted lemma](weighted_obstruction.md). At a fixed height $H$ retain offsets $0\le h<H$ in a set $J$ and let $L=\sum_{h\in J}(H-h)$. Suppose an integer $K$ includes all layers through $H$ and is divisible by every $w=p^e>H$ with $w\mid N-h$, $h\in J$, for which $A$ avoids $p$. The weighted product of the other two rows, with cutoffs totalling $H$ on each row, is at least $(N/(2K))^L$ when $N\ge2H$. Indeed, it is divisible by the product of $((N-h)/\gcd(N-h,K))^{H-h}$ over $J$.

Fix $H=300$. If the two ordinary rows both have minima at most $H$, the full powers of primes of $N$ avoided by them have product at most $H^2$. The complementary product divides every part of $A$, hence divides $P$. Thus $m\le H^2$. We handle these finitely many multipliers below. Otherwise one ordinary row can take cutoffs $(100,100,100)$ and the other $(0,150,150)$ for all sufficiently large $N$. The combined upper bound on their weighted product is $N^{37{,}800}$.

## Retaining offsets under two power relations

Let $u,v,z$ be candidate residues of $2^r,2^s,2^t$ respectively, all nonnegative with $u+v+z=h<H$. For each of $e=r,s$ write $x=u,v$ respectively and form a constraint.

- If $e\le\sqrt t$, use $x=2^e$.
- Otherwise the pigeonhole argument for the odd three-bit case gives coprime $1\le a\le b\le\lceil\sqrt t\rceil$ with $D=be-at$ satisfying $|D|\le\sqrt t$. Use $2^{\max(-D,0)}x^b=2^{\max(D,0)}z^a$. [Short-relation proof](three_binary_digits.md)

Any carry-free layer $w>H$ dividing $N-h$ satisfies both constraints modulo $w$ at its actual residues. Omit an offset if both constraints hold as exact integer equalities for some candidate triple. At every retained offset, take the product of one nonzero difference from the two constraints for each triple, using absolute values. Multiply all these products by $\operatorname{lcm}(1,\ldots,H)$ to obtain $K$. There are $O_H(1)$ differences and each has logarithmic size $O_H(\sqrt t)$, including the first type. Hence $K=N^{o(1)}$ uniformly.

Here is the count of omitted offsets. If the first type of constraint occurs with $2^e\ge H$, none are omitted. Both exponents cannot satisfy $2^e<H$ for sufficiently large $N$, by the second-part bound for $A$. In the remaining cases at least one power relation is present.

If a relation has $b\ge2$, it permits at most $1+\lfloor\sqrt H\rfloor=18$ possibilities for $z$. Indeed its positive integer solutions have $z=2^{j_0}d^b$ for a fixed $j_0$ and positive integers $d$, by taking prime valuations and using coprimality of $a,b$. The common zero gives the one other possibility. For fixed $z$ each constraint determines its $x$ uniquely on the nonnegative reals. There are therefore at most 18 omitted offsets.

Otherwise each power relation has $a=b=1$. With two such relations the offset is a multiple of $1+2^{s-r}+2^{t-r}=m$. For $m\ge17$ this omits at most 18 offsets, and smaller $m$ is among the fixed-multiplier cases below. With just one power relation we have $q=2^r<H$ fixed in a finite list and the offsets lie in $q+(1+2^{t-s})v$. If $t-s\ge4$ there are at most 18 of them. The remaining fixed shifts are treated below.

In all the cases with at most 18 omissions, $L\ge45{,}150-18\cdot300=39{,}750>37{,}800$. The lower bound $N^{L-o(1)}$ contradicts the upper bound from the two ordinary rows.

## The bounded parameters

For these cases take a prime $\ell\le N$ with $g=N-\ell=O(N^{0.57})$. Guth–Maynard’s Corollary 1.3 gives such a prime for every sufficiently large $N$. In each case covered below $A$ has two parts bounded below by a fixed positive multiple of $N$, so $A$ cannot avoid $\ell$. One ordinary row therefore has its smaller two parts summing to at most $g$. [Guth–Maynard, Corollary 1.3](https://ora.ox.ac.uk/objects/uuid%3Aad11b8bf-ad2b-4ebf-a627-647f023c378f/files/s73666708v)

For fixed $m\le300^2$ we have $m\ge7$. Use height 12 and omit offsets divisible by $m$. The identity $m2^r=N$ implies $w\mid h-mu$ at a carry-free layer of $A$. This is a nonzero integer bounded in absolute value by a constant depending on $m$ at every retained offset. Thus a constant suffices for $K$. The retained degree is at least $78-12-5=61$. The [height 12 cutoff table](binary_weights.md) bounds the combined cost of the two ordinary rows by $O(N^{56.52})$: it uses only the smaller-part sum of at least 11 for a non-endpoint row, at most one row with a part at most 2, and the bound on $g$. This is a contradiction.

It remains to cover $N=q+(1+2^d)2^s$ with $q=2^r<300$, $r\ge1$, and $d=t-s\in\{1,2,3\}$. When $q\ge4$, use height 3. For each $h=0,1,2$ and each odd $w=p^e>3$ dividing $N-h$, the three residues of $A$ are nonzero and cannot sum to $h$. The 2-parts of $N,N-1,N-2$ are $q,1,2$ respectively. Consequently all three offsets are retained with constant $K$, giving lower degree 6. Cutoffs $(1,1,1)$ on each ordinary row cost at most $N^4g^2=O(N^{5.14})$.

When $q=2$, both ordinary rows have minimum at least 3. Use height 12 and omit offsets $h\equiv2\pmod{1+2^d}$. The relation $(1+2^d)2^s=N-2$ again gives bounded layers on the retained offsets. The largest loss is $10+7+4+1=22$ at $d=1$, leaving degree at least 56. The first three lines of the height 12 cutoff table bound the combined cost by $O(N^{50.1})$. This final contradiction proves the theorem.

## Every three-bit weight below $2^{128}$

**Finite theorem.** The full conjecture holds at all $\binom{128}{3}=341{,}376$ weights with three binary ones below $2^{128}$.

The endpoint test uses the exact CRT reduction above. After the odd-weight exclusion and the bound $27(N-1)^2\le P^3$, just 31,862 weights require factoring. In this checker the three groups $d_i$ are formed directly. The normalized part for group $d$ is $(M/d)\big((mM/d)^{-1}\bmod d\big)$, with the inverse taken in $[1,d-1]$ and $M=N-1$. Retain the values smaller than $P$ summing to $P$. Only 56 unordered CRT rows survive across the whole range, and each has positive factorial valuation at a prime of $N-1$.

For the independent case without $E$, write $A=(a,b,c)$. If $a+b\le10$, its sum is at least 3 and Proposition 10 applies. Otherwise, for $h=0,\ldots,5$, compute
$$D_h=\operatorname{lcm}\big(\gcd(N-h,60),\ \gcd(N-h,a-i,b-j)\text{ for }i,j\ge0,\ i+j\le h\big).$$
This removes every small layer and every larger layer which $A$ could cover, by the same two-residue condition as in the [one-million proof](full_range_1m.md). Therefore $F=\prod_h((N-h)/D_h)^{6-h}$ is at most the total height 6 cost of the ordinary rows. Write $B=\lfloor N^9/3^9\rfloor$ and $C=\lfloor64N^9/3^9\rfloor$. At most one row contains 1, so the ordinary product is bounded by $BC$.

If $F\le BC$, look for a proved prime $c<\ell\le N$. The binary row cannot avoid it. An ordinary row avoiding it has its two smaller parts summing to at most $g=N-\ell$. If $g\le10$ there is no such non-endpoint row. Otherwise the row costs at most $g^{15}$ by cutoffs $(1,5,0)$ for minimum 1 or $(2,4,0)$ for minimum at least 2. The remaining ordinary row costs at most $C$, so $F>Cg^{15}$ excludes the selection.

| Deciding step without the endpoint | Weights |
|---|---:|
| Binary row has smaller-part sum at most 10 | 624 |
| $F>BC$ | 246,944 |
| Proved prime gives $g\le10$ | 15,389 |
| $F>Cg^{15}$ | 78,419 |

Run `python3 check_three_bit_range.py 128` with python-flint installed. All comparisons use integer arithmetic. The verifier confirms each prime factor returned by the factorization routine and uses a proving primality test for the nearby primes. [FLINT primality tests](https://flintlib.org/doc/fmpz/fmpz_is_prime.html)
