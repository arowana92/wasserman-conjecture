# Exponents that work simultaneously for every prime and every cofactor through 1000

Let $1\le d\le1000$, let $q$ be any prime, and let $e\ge1$. **Wasserman's full trinomial conjecture holds at $1+dq^e$ if $e$ is divisible by at least one of**
$$
 2,\quad 21,\quad 51,\quad 595,\quad 2523,\quad 17017,\quad 24871,\quad 83657 .
$$
These exponent classes have natural density $11490087461/21558946695>0.5329$. The same classes work for all the primes and cofactors simultaneously.

For **$3\mid e$**, a possible failure requires $d=986$, $e$ odd, $7\nmid e$, $17\nmid e$, and $29^2\nmid e$. If also $29\mid e$, only one distinguished row remains, namely $(203Q,714Q,69Q+1)$ for $Q=q^e$.

## The two rows left at a cube

The [primitive denominator theorem](shifted_directions.md) and the [exclusions of 364 and 880](denominator_residue_hulls.md) leave two cube-eligible rows, both at 986. These restrictions also hold for an arbitrary integer cube $W$ used as the full primitive scale $\gcd(x,y,z-1)$ of a nonendpoint row. At a positive power $q^e$ with $d\le1000$, the earlier prime-power reduction shows that a possible failure must have $d$ as its primitive denominator and $W=q^e$. All even exponents were excluded there.

| $(a,b,c)$ in $(aW,bW,cW+1)$ | $W\bmod6902$ | $W\bmod17$ | $W\bmod29$ |
|---|---:|---:|---:|
| $(119,203,664)$ | 4215 | 16 | 10 |
| $(203,714,69)$ | 3501 | 16 | 21 |

Both rows must avoid 17 and 29. The only nonzero seventh powers modulo 29 are $1,12,17,28$, excluding a simultaneous seventh power in both cases. For a prime $p$ and $h\ge1$, if $W=x^{p^h}$ and $W\equiv s\not\equiv0\pmod p$, then
$$
 W\equiv s^{p^h}\pmod{p^{h+1}} .
$$
Indeed Fermat's theorem gives $x\equiv s\pmod p$, and raising congruent units modulo $p^k$ to the $p$th power preserves the congruence modulo $p^{k+1}$, by the binomial theorem.

The table gives the resulting overflowing layers. Each residue sum exceeds its layer, which suffices for a carry by Legendre's formula.

| Slopes | Required power | Layer | $W$ at the layer | Three part residues | Sum |
|---|---:|---:|---:|---|---:|
| $(119,203,664)$ | 17 | 289 | 288 | $170,86,204$ | 460 |
| $(203,714,69)$ | 17 | 289 | 288 | $86,153,221$ | 460 |
| $(119,203,664)$ | 29 | 841 | 416 | $726,348,377$ | 1451 |
| $(203,714,69)$ | $29^2$ | 24389 | 9388 | $3422,20446,13659$ | 37527 |

This proves the conclusions for exponents divisible by $21,51$ or $2523=3\cdot29^2$. Lean 4.36.0 independently checks the seventh-power exclusions and all three prime-power lifts for arbitrary natural roots. [Kernel-checked modular implications](CofactorPowerResidues.lean)

## Combining different obstructions

For an odd prime $s$ and a remaining shape with $W\equiv R\pmod M$, first test whether $R$ is an $s$th power modulo each full prime power of $M$. If any fails, this shape excludes all $s$th powers. Also test the preceding lift when $s$ is a base the distinguished row must avoid and $M$ fixes a nonzero residue modulo $s$. Try every possible lift through the third base-$s$ digit and reject only if all lifts produce a carry. Python computes each carry both by factorial layers and by individual digit sums.

For example, $595=5\cdot7\cdot17$ works at every scale parity. Among the 19 remaining combinations of denominator, shape and parity, one may use the following primes. Each listed prime excludes every corresponding row at that denominator.

| Power obstruction | Denominators |
|---:|---|
| 5 | 380, 589, 710, 738, 874, 946, 950, 962 |
| 7 | 703, 780, 850, 986 |
| 17 | 629, 646 |

For completeness, when the scale is odd the inclusion-minimal sets of odd primes through 79 meeting every shape's exclusion list in this test are
$$
 \{3,7\},\ \{3,17\},\ \{5,7,17\},\ \{7,11,13,17\},\
 \{7,11,17,19\},\ \{7,17,19,37\}.
$$
Their products give six of the displayed exponent classes. For a prime-power scale, the even-scale case $q=2$ was already [settled for all exponents and all cofactors through 1000](prime_base_density.md#from-trinomial-rows-to-residue-classes-of-the-base). For unrestricted integer scales of both parities, the first four sets still work. The other inclusion-minimal sets from this test are $\{7,11,17,19,41\}$, $\{7,13,17,19,37\}$, and $\{7,17,19,37,41\}$. These minimality statements refer to the specified residue and digit tests.

To reconstruct the lists, start with the empty set and process each shape in turn. Adjoin an element of its exclusion list in every possible way and discard strict supersets. Induction gives exactly the minimal sets meeting all processed lists. The density in the theorem follows by inclusion-exclusion over the eight positive integers displayed there, using reciprocal least common multiples for intersections.

Run the enumeration with Python and python-flint, and the modular verification with Lean 4.36.0.

~~~~sh
python3 check_power_obstruction_cover.py
lean CofactorPowerResidues.lean
~~~~
