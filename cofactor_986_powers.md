# Separating the exceptional exponents at cofactor 986

Let $q$ be any prime and $E$ any positive integer. **Among the one million weights**
$$
 1+986q^{Ej},\qquad 1\le j\le10^6,
$$
**at most one can fail Wasserman's three-row conjecture.** Equivalently, two failing exponents $a<b$ at the same prime would satisfy $b>10^6\gcd(a,b)$.

Consequently, among **all one billion choices** $1+dq^{3Ej}$ with $1\le d\le1000$ and $1\le j\le10^6$, at most one can fail, and its cofactor must be 986. The [earlier cube theorem](denominator_residue_hulls.md) supplies every other cofactor.

There is an additive restriction as well. For each $q$, the exceptional exponents at 986 **contain no four distinct terms in arithmetic progression**. If they contain three terms in arithmetic progression, the common difference is divisible by $68\cdot7^{18}$. The starting exponent and the common difference have no size bound.

## A single row and exponent class for each prime

The [primitive row classification](shifted_directions.md) and the [prime-base reduction](prime_base_density.md#from-trinomial-rows-to-residue-classes-of-the-base) give the following necessary conditions at $1+986q^n$. The base 2 is settled. For any other prime, $n$ is odd. Put $W=q^n$. A failing weight requires one of two rows $(aW,bW,cW+1)$, with slopes $(119,203,664)$ or $(203,714,69)$. Both must add without carries in bases $7,17,29$. In their respective order $W$ is $4215$ or $3501$ modulo $6902$. These are all the odd-scale rows at this denominator, including at exponents outside the cubic case.

The table uses low digits that can extend to a complete carry-free row.

| Slopes | $W\bmod49$ | $W\bmod289$ | $W\bmod29$ |
|---|---|---|---|
| $(119,203,664)$ | 22 or 29 | 84 | 10 |
| $(203,714,69)$ | 22 | 84 | 21 |

To verify extendibility, start the multiplication carries at $(0,0,1)$. For input digit $d$ compute $v_i=s_i+a_i d$, discard it if $\sum_i(v_i\bmod p)\ge p$, and otherwise pass to $(\lfloor v_i/p\rfloor)_i$. The states are bounded. Delete states from which no accepting state can be reached, then enumerate the indicated prefixes. A state accepts if its entries add without carries after the input ends. The independent checker instead tries each possible prefix, checks its factorial floors, and searches the remaining graph for the all-zero state. A completed word can reach that state by padding with zeros. At the second row this also gives the single residue $5608\bmod 7^5$.

Both residues at 29 are generators, with $21=10^{19}\bmod29$. Thus $q$ generates the units modulo 29 and $n$ is coprime to 28. At 7 write $W=1+7t$ modulo 49, with $t=3$ or 4. If $q^6=1+7u$ modulo 49, comparison of $W^6$ with $(q^6)^n$ gives $un\equiv6t\bmod7$. In particular $u\ne0\bmod7$. For two exceptional exponents $n,n'$ at this $q$, the ratio $n'/n\bmod7$ is therefore $1$ or $-1$. On the other hand modulo 28 the ratio is $1$, $19$ or $19^{-1}=3$, according to the two residues at 29. Only $1$ passes the mod-7 condition. Thus the row stays the same and $n'\equiv n\bmod28$.

At 17, $q^n\equiv-1\bmod17$ and $n$ is odd, so $q\equiv-1\bmod17$. Write $q=-1+17v$. The common residue $84=-1+5\cdot17$ gives $vn\equiv5\bmod17$. Hence all possible exponents are units in a **single residue class modulo 476**. With the second row they even occupy one class modulo $68\cdot7^4$, since $v_7(q^6-1)=1$ and $q^{6n}$ is fixed modulo $7^5$. The binomial expansion for principal units and the finite implication from the two ratios are also checked in [Lean](CofactorPowerRatios.lean).

## The million multipliers

If some $n=Ej$ fails, then $E$ and $j$ are units modulo 476. At all failing indices the same row occurs and the indices occupy one class modulo 476. The unit $x=q^E$ has order 1 or 3 modulo 7 since $x^j\equiv1\bmod7$ with $j$ odd. In the latter case every failing $j$ is divisible by 3. Replace $x$ by $x^3$ and $j$ by $j/3$ throughout the comparison, preserving the common class property and the bound. It suffices to allow **any** $x\equiv1\bmod7$.

For each row and index class, grow the base-seven digits of $x$. At a current prefix $x_0\bmod 7^k$, keep the indices $m$ whose power $x_0^m$ has readable low digits up to this level, together with its multiplication-carry state. Initially $x_0=1$, $k=1$. For each digit $d$,
$$
 (x_0+d7^k)^m\equiv x_0^m + m d 7^k\pmod{7^{k+1}},
$$
because $x_0\equiv1\bmod7$. Thus each allowed next output digit determines a unique next input digit, using $7\nmid m$. Discard a child as soon as fewer than two distinct indices remain there. Exhaustion proves there are no two simultaneous indices even at an unbounded input. The two searches visited 22,740,603 and 561,222 nodes. Their largest examined depths were 27 and 18.

## Independent verification by rational reconstruction

The second verifier directly enumerates $x=W\bmod7^{14}$ with $W\equiv1\bmod7$ that pass the first 14 factorial-floor inequalities. It retains 4,136 and 168 words for the two rows, including words that might never finish. Write $x=8^z\bmod7^{14}$; the base 8 has order $7^{13}$ by the binomial expansion. Each $z$ is a unit at 7. For two words $x,y$ with logs $z,z'$ and reduced positive exponents $a,b$, the required identity $W_1^b=W_2^a$ gives
$$
 a/b\equiv z/z'\pmod{7^{13}},\qquad a/b\equiv1\pmod{68}.
$$
Reduction by the gcd preserves the second congruence since both original exponents are units modulo 476. Equal words would force $a=b$ within the bound. For distinct words combine the congruences to a residue $r\bmod M$, $M=68\cdot7^{13}>2(10^6)^2$.

There is at most one reduced positive ratio with numerator and denominator at most $B=10^6$ in this class, by the cross determinant. Run the Euclidean algorithm on $(M,r)$ until the first remainder at most $B$, tracking its signed coefficient of $r$. To see why this finds every solution, $(rb-a)/(Mb)$ differs from $r/M$ by at most $B/(Mb)<1/(2b^2)$ and has reduced denominator $b$. Here $b$ is coprime to $M$ and $\gcd(a,b)=1$. Thus it is a continued-fraction convergent. Its remainder has magnitude $a$ and its signed coefficient is $b$ when the remainder is positive. After the first small remainder, the next coefficient magnitude must exceed $B$. Otherwise both adjacent remainders and both coefficient magnitudes would be at most $B$ by monotonicity, and the determinant identity for neighboring stages would give $M\le2B^2$. Hence it suffices to test the first small remainder.

For a reconstructed pair lift the congruence $x^b=y^a$ one digit at a time. Modulo the next power, the first added digit $i$ determines the second $j$ by
$$
 aj\equiv (x^b-y^a)/7^k +bi\pmod7.
$$
Test both new words by factorial floors. This eliminates all 409,986 reconstructed prefix pairs in the first row, and all 407 in the second row, within another 16 and 6 digits respectively. This verification uses neither the carry states nor the common-base search tree of the first computation.

## Three or four equally spaced exponents

If $W_i=q^{n+iL}$ all fail, they use the same row and coincide modulo 7. Put $h=v_7(q^L-1)$; then $h=1+v_7(L)$, since $v_7(q^6-1)=1$ and $7\nmid6$. The valuation follows by the binomial expansion of powers of a unit congruent to 1 modulo 7. The $W_i$ have the same prefix $r\bmod7^h$, and the first two differ at digit $h$.

For four terms and $h\ge2$, write $W_i=r+7^h y_i$ and $q^L=1+7^h v$. Then $y_i\equiv y_0+i vW_0\bmod49$. The step is a unit at 7. From every productive carry state reachable after the initial digit 1, the checker lists the readable two-digit words. None of these lists contains a four-term progression modulo 49 with step a unit at 7. At $h=1$, test the actual geometric recurrence through two more digits as follows.

For three terms and each $1\le h\le18$, generate every readable common prefix $r\bmod7^h$ starting in digit 1. Extend the first two values $x,y$ by all readable digits, requiring them to differ at their first extension. Modulo every new layer the ratio is $y/x$ and the third term is $y^2/x$. Require its new digit to be readable as well. For four terms also test $y^3/x^2$. Every search exhausts. At most 42 digits past the common prefix were needed for any three-term case, and the four-term case at $h=1$ exhausts within two digits. Thus a failing triple requires $v_7(L)\ge18$ as well as $476\mid L$. Four terms are impossible at every step.

Run from the repository root. Both Python programs use standard-library arithmetic.

~~~~sh
python3 check_cofactor_986_powers.py
python3 check_cofactor_986_powers_cross.py
lean CofactorPowerRatios.lean
~~~~
