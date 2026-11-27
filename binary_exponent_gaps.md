# Gaps between possible exceptional binary exponents

Write $N=1+d2^a$, $a\ge1$. Consider the following sets of cofactors.

- $S$ contains all 32,768 positive odd integers through 65,535.
- $T$ contains $5Q$ for every odd prime power $Q<2^{23}$.

An exponent is exceptional for a set if Wasserman's full conjecture fails at even one of its cofactors. Throughout, the starting exponents are unrestricted.

| Cofactors | Terms sufficient to include an exponent good at every cofactor, in any positive-step arithmetic progression | Terms sufficient at step one |
|:---|---:|---:|
| $S$ | 5 | 3 |
| $T$ | 3 | 3 |
| $S\cup T$ | 12 | 4 |

For $S$, four terms already suffice unless the step is $125,155,265$ or $295$ modulo 420. At every single exponent, **at most one cofactor in $S$ can fail**. The following six disjoint residue classes determine which cofactor it could be. Together they occupy 70 of the 420 exponent residues. They follow by reducing the necessary exponent progressions.

| Possible cofactor | Necessary exponent class |
|---:|:---|
| 12031 | $7\pmod{12}$ |
| 56279 | $3\pmod{21}$ |
| 60829 | $189\pmod{210}$ |
| 61423 | $30\pmod{140}$ |
| 65095 | $60\pmod{140}$ |
| 65243 | $28\pmod{60}$ |

**Separation.** In any interval of diameter $10^{15}$ there are at most **four** distinct exceptional exponents for $S\cup T$. Within that interval at most **one** exponent can fail at any specified cofactor. At $d=12031$ the one-exponent assertion holds in every interval of diameter $10^{48}$, including for prime-power exponents across all of $S$, since no other cofactor in $S$ can fail at a prime-power exponent. For the whole of $T$ there is at most **one** exceptional exponent in every interval of diameter $10^{128}$, even if the cofactor may differ between the two putative exceptions.

The proof uses the higher-layer exclusions which leave seven rows at six cofactors in $S$, and two rows for $T$. The only possible carry bases needed below are grouped in the table. The last column bounds the growth exponent of the number of possible exponent residues in that base. [Higher-layer exclusions](binary_perfect_power_exponents.md) · [Small-cofactor carry graphs](binary_cofactor_rays.md) · [Factor-five carry graphs](factor_five_exponents.md)

| Carry base $p$ | Possible cofactors | Growth exponent $\delta_p$ |
|---:|:---|---:|
| 3 | $12031,65243,5\cdot1467779,5\cdot2823523$ | $.14$ |
| 5 | $65095$ | $.22$ |
| 7 | $56279$ (two rows), $60829$ | $.29$ |
| 11 | $61423$ | $.36$ |

## Ruling out a whole progression

Label a putatively exceptional exponent by one necessary row, with exponent class $a_i$ modulo $H_i$. Suppose the exponents are $b+jh$. For a prescribed sequence of labels, the CRT in $b$ requires, and is implied by, the pairwise congruences
$$
(j-i)h\equiv a_j-a_i\pmod{\gcd(H_i,H_j)}.
$$
Solving these single-variable congruences gives either an empty set or one class of steps. A fixed step can be substituted directly.

There are additional obstructions for repeated labels. At each small-cofactor row an earlier carry-graph calculation excludes a three- or four-term arithmetic progression even at unrestricted exponents. At the two surviving small ternary rows, each reachable state admitting finite completion has at most two possible next digits after the fixed units digit. Three exponents with equal pairwise 3-adic difference valuations would force three distinct next digits after a common prefix, because same-parity exponents satisfy $v_3(2^a-2^{a'})=1+v_3(a-a')$. The same exclusion holds for each factor-five row. Its productive graph starts after the first two ternary digits and has maximum outdegree two. If three scalars first split at the second digit, one gives $9\mid N$, which the factor-five reduction excludes. Equality of the three valuations on positions $i,j,k$ is preserved after multiplying the differences by $h$. [Short progression tests and individual row bounds](binary_exponent_progressions.md) · [Factor-five initial digits and state degrees](factor_five_exponents.md)

Appending a label, testing these obstructions and merging the step congruences gives the following counts of surviving ordered label strings, starting at length one and ending at the first empty length.

| Set and step | Counts |
|:---|:---|
| $S$, arbitrary | $7,49,62,8,0$ |
| $S$, step one | $7,10,0$ |
| $T$, arbitrary | $2,4,0$ |
| $S\cup T$, arbitrary | $9,81,228,242,132,134,144,160,152,132,84,0$ |
| $S\cup T$, step one | $9,19,2,0$ |

A second verifier checks all pairs of starting residue and step directly. For $S$ it uses the six classes modulo 420, relying on the four-term exclusion for the cofactor 56279 across its two rows. For $S\cup T$ it suffices to project these six classes to modulus 84 and add the two factor-five labels there. These smaller rings independently prove the progression lengths in the first table. The full CRT calculation restricts the exceptional four-term steps to the four classes given above.

For $T$ there is also a short direct argument. Its two labels have exponents respectively $3$ and $2$ modulo 4. Three terms in those residues must have step divisible by 4, giving a single label, where the ternary obstruction applies.

## Ruling out two exponents at any height

Take any two rows using the same prime $p$, with periods $H_1,H_2$. Both force their scalar $q=2^a$ modulo $p$. Start with $P=\gcd(H_1,H_2)$, which satisfies $2^P\equiv1\pmod p$ in all pairs used. Let $w$ be the largest power of $p$ dividing $2^P-1$. An exponent residue modulo $P$ determines its scalar modulo $w$. At the next stage lift a residue $e$ to $e+iP$, $0\le i<p$, keeping its congruence modulo $\gcd(pP,H_j)$ for the corresponding row $j$. These are all possibilities modulo $pP$. The new scalar is $2^e(2^P)^i$. Its next digit is an affine permutation of $i$ since $(2^P-1)/w$ is a unit modulo $p$. The next width is $pw$, by binomial expansion at the odd prime $p$.

The first verifier retains a lift only when the new digit follows a carry-free transition to a state admitting finite completion. For pairs involving $T$, the other verifier checks the individual factorial floors from the part residues at each new ternary layer, with the required $9\nmid N$. The small ternary rows in these mixed checks also satisfy $9\nmid N$. Both verifiers maintain sorted lists of possible absolute exponent residues.

Once $P>B$, two distinct exponents separated by at most $B$ must reduce to distinct residues at cyclic distance at most $B$ modulo $P$. Remove residues with no such partner in the other list, or in their own list for a repeated row. It suffices to inspect the strictly preceding and following residues, including wraparound. Removing a residue also removes all its later lifts. Thus exhaustion at a finite stage certifies **every** pair of positive exponents of gap at most $B$, including pairs above every prescribed initial bound.

| Checks | Gap bound $B$ | Maximum refinement depth | Largest retained list size, summing both lists for distinct rows |
|:---|---:|---:|---:|
| Small rows at the same prime, apart from the extended 12031 check | $10^{15}$ | 41 | 1,493,201 |
| 12031 against itself | $10^{48}$ | 107 | 237,108 |
| One small and one factor-five ternary row | $10^{15}$ | 37 | 92 |
| The three unordered factor-five row pairs, allowing repetition | $10^{128}$ | 290 | 638,163 |

A set of five distinct exceptional exponents in an interval of diameter $10^{15}$ would have two of the same carry base. This proves the four-exponent bound.

## Most steps give a uniform bound at every starting exponent

Let $K\ge2$ and $L\ge2$ be integers. Among any $K$ consecutive positive steps, at most
$$O(LK^{0.72})$$
can admit **any** length-$L$ exponent progression containing more than four exceptional indices for $S\cup T$. Thus one may take $L=\lfloor K^{0.27}\rfloor$ and discard only $O(K^{0.99})$ steps for sufficiently large $K$. For $T$ alone the bound for steps admitting two or more exceptional indices is $O(LK^{52/625})$, allowing $L=\lfloor K^{0.9}\rfloor$ with a power saving.

Indeed the number of exponent residues modulo $p^n$ arising from a row in the table is $O(p^{n\delta_p})$. To transfer the cited productive-path bounds from scalars to exponents, put $t=v_p(2^{\operatorname{ord}_p(2)}-1)$. The scalar modulo $p^{n+t}$ determines the exponent modulo $p^n$ by the order-lifting formula. With only the two rows of $T$, take $\delta_3=26/625$. The corresponding set of differences for two exponents in the same group occupies $O(p^{2n\delta_p})$ residues, uniformly in $n$. In particular this difference set has zero density even when neither exponent is bounded.

Write $\beta=2\delta_p<1$ and take $p^{n-1}<K\le p^n$. For positive $m<L$, let $v=v_p(m)$. On any $K$ consecutive steps, $h\mapsto mh\bmod p^{n+v}$ is injective. At most $O(K^\beta p^{\beta v})$ steps give such an exceptional difference. Summing in $m$ costs $O(L)$ times the bound $K^\beta$, since at most $L/p^j$ indices are multiples of $p^j$ and $\sum_{j\ge0}p^{(\beta-1)j}$ converges. Sum over the four groups. Every progression with more than four exceptional indices contains a same-group pair with difference $mh$. For $T$ there is just one group.

Run the following from the repository root with Python 3.9 or later.

```sh
python3 check_binary_exponent_patterns.py
python3 check_binary_exponent_patterns_cross.py
python3 check_binary_exponent_gaps.py
python3 check_factor_five_gaps.py
```
