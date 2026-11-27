# Two freely placed binary powers at offsets 83 and 2579

**For $c=83$ or $2579$ and all integers $a,h,n$ with $a\ge1$ and $n>h\ge1$, every three positive trinomial rows at**
$$N=1+2^a(2^n+2^h+c)$$
**have a common prime divisor.** All three exponents vary independently subject to these inequalities. The offset 2579 has five binary ones, so the family includes weights with eight distinct binary ones.

The occupied positions of 83 are $0,1,4,6$, and those of 2579 are $0,1,4,9,11$. In either list all positive pairwise distances are different. The [two-block theorem for one shared distance](distinct_block_ratios.md#one-shared-distance-can-also-work) therefore covers every $2^n>(2c)^{64}$. This leaves $n\le472$ and $n\le789$, respectively. The calculation below covers every outer exponent at each remaining pair.

## Finite reduction

Set $P=2^h$, $u=n-h$, $B=1+2^u$, $D=BP+c$. Discard the case $27D^2\ge(P+c)^5$ by the [single-high-power inequality](successor_powers.md). For $P>c$ apply the two-block test outright if $27D^2>(Bc)^7$ and no distance in $c$ equals $u$. If there is a match at $v,v+u$, delete those two binary units from $c$, leaving $R$. The same lower bound together with $729D^4>P^3(R+1)^6R^4$ excludes this case by the one-shared-distance theorem. Each test uses integers only.

In a remaining counterexample the row avoiding 2 must have parts $2^a x,2^a y,2^a z+1$, with $x,y,z>0$ a binary carry-free split of $D$. Take $x\le y$. An assigned weight scale $U$ divides $\gcd(x,y)$ and is coprime to $2D$. If $S$ is the predecessor factor from $D$ assigned to the same row, then $S\mid\gcd(D,xy)$ and $27D^2\le U^3 S^2$. With $G$ the part of $\gcd(x,y)$ coprime to $2D$ and $S_0=\gcd(D,xy)$, discard $G^3 S_0^2<27D^2$. [Assignment proof](distinct_block_ratios.md#three-determinants)

The first verifier distributes the two large powers and the lower units among the three parts. When summing the original units in $D$ produces a carry, every true split is still among these allocations: merge the original units to form the final ones, and give each merged unit and all its constituents to its chosen part. At the last stage require the actual bits of $x,y,z$ to be disjoint.

Since $U\mid 2^aD+1$ and $U\mid G$, if the powers of two modulo $G$ return to 1 within $4n+64$ steps, exhaust this complete cycle. At a residue $v$ replace $G$ by $\gcd(G,vD+1)$ and repeat the product test. If the cycle has not closed, keep the full $G$. For each remaining scale cap, the companion CRT test in the [binary block computation](distinct_block_ratios.md#every-binary-gap-in-a-finite-set-of-block-pairs) excludes the case. This includes all positive outer exponents regardless of the search cutoff on the cycle.

| $c$ | Remaining $(n,h)$ | Allocations tested | Passing product and carry checks | Closed cycle | Excluded by the cycle | Rows reaching companion checks |
|---:|---:|---:|---:|---:|---:|---:|
| 83 | 37,349 | 10,084,279 | 114 | 114 | 114 | 0 |
| 2579 | 105,770 | 95,510,427 | 511 | 505 | 501 | 10 |

## Independent verification

The second verifier generates each pair from the union of its bits inside the actual binary expansion of $D$. It fixes the lowest bit of the union in one part. It uses a different elimination of the outer exponent. If the current cap is $G$ and $m=\gcd(G,2^j-1)$, then
$$U\mid \frac{G}{m}\gcd\big(m,(-D)^j-1\big).$$
Indeed $\gcd(U,m)$ divides $(-D)^j-1$ because $(-D)2^a\equiv1\bmod U$ and $2^j\equiv1\bmod m$. Also $U/\gcd(U,m)$ divides $G/m$. Updating the cap at $1\le j\le2n+64$ gives a necessary condition at every outer exponent. The remaining cases use the possible single CRT coordinates of the normalized companion. The remaining subset tests number 3 at offset 83 and 15 at offset 2579. All fail. Trial division suffices for every factorization used by either program.

Run the two verifiers from the repository root.

~~~sh
python3 check_two_free_binary_powers.py
python3 check_two_free_binary_powers_cross.py
~~~
