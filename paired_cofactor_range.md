# An all-cofactor bound from the remaining binomial cut

Under the [two sparse-layer hypotheses](three_unit_layers.md), the trinomial conjecture holds for **every weight through $10^{1000}$ at all pairs of primes $p,q<257$**. It also holds for **every cofactor $D$ without a size limit** at these pairs whenever both moving prime powers are at most $10^{334}$.

## An elementary alternative at arbitrary bases

Use $P=p^a$, $Q=q^b$. In an exception the common remainder must be $s\mid D$. Put $E=D/s$, $W=EPQ$ and $w=s-1$. One row contains $wW$, and its two remaining parts sum to $W+1$. All three are divisible by a unitary factor $U\ge3$ of the weight with $U\mid w$.

**If $P+Q+1\le s^4E$, then $27E^2\le U^3$.** Every prime factor of $E$ must be assigned to the companion row at its full power in the predecessor. Thus whenever $E>((s-1)/3)^{3/2}$, a failure requires
$$
                              s^4E<P+Q+1.                                \tag{1}
$$

To prove this, interchange the powers so that $P>Q$ and set $\Delta=P-Q$. The [two-base binomial cut](three_unit_layers.md#three-simultaneous-binomial-cuts) forces one of the remaining parts to be $B=(W-Q)/\Delta$. Its complement is $C=W+1-B$, and
$$
                           (\Delta-1)B-C=-Q-1.                            \tag{2}
$$
For clarity, at the chosen layers one has $W/P=1+t_p L_p$ and $W/Q=1+t_q L_q$, with $0<t_p<Q/s^2$ and $0<t_q<P/s^2$. The assumed inequality gives $t_pt_q(P+Q+1)<W$, the hypothesis of the binomial result. Both parts are at least $U$, so the cut is interior.

Take a prime $r\mid E$ and its full power $H\parallel N-1$ (that is, $H=r^{v_r(N-1)}$). If this prime goes to the grouped row, its three parts are $0,0,1$ in some order modulo $H$. The part $wW$ is divisible by $r$, so the one must occur at $B$ or $C$. In the first case (2) gives $H\mid Q+\Delta=P$. In the second it gives $H\mid Q$. Both are impossible since $E$ is coprime to $pq$. Thus the companion takes every prime of $E$, including primes shared by $s$ and $E$. Its three normalized positive entries sum to $U$. For each $H$ just considered, two normalized entries are divisible by $H$ because the normalizing scale divides $N$. Their triple product is therefore a multiple of $E^2$ and is at most $(U/3)^3$, proving the claim. The Lean kernel checks the implication from (2) and the two possible residue assignments to divisibility of $P$ or $Q$. [Lean lemma](ResidualPowerClasses.lean)

## Finite enumeration and its extension

The classification at bases below 257 left only 25 triples $(p,q,w)$ and gives $7\le s=w+1\le2211$. For each one we enumerate exactly the positive triples $(E,a,b)$ for which $N=1+sEp^a q^b\le10^{1000}$ and the paired remainders equal $s$. This is equivalent to
$$
                   Eq^b\equiv1\pmod{L},\qquad Ep^a\equiv1\pmod{M},
$$
where $L$ and $M$ are the smallest powers of $p$ and $q$, respectively, at least $s^2E$. These congruences also force $E$ to be prime to $pq$.

Fix $b$, set $H=(N_{\max}-1)/(spq^b)$ and let $A$ be the largest power of $p$ at most $s^2$. In any interval $[A^{j-1},A^j)$ with $E\le H$, the first modulus must be at least $A^j$. Hence the only possible $E$ is the least residue of $q^{-b}$ modulo $A^j$. Test it at its actual modulus $L$. In binary the equivalent implementation lists candidate prefixes of the modular inverse by locating a leading one followed by enough zero bits. For the second modulus a table of inverse powers of $p$ gives every possible $a$ up to the weight bound. Once a complete period would exceed that bound, lookup at the smaller modulus suffices to get a single candidate exponent, tested again at $M$.

The program repeats the enumeration with the roles of $p,q$ reversed and requires the full sets of exponent and cofactor triples to agree. An additional small check loops directly over $E,a,b$ and extracts the valuations of $Eq^b-1$ and $Ep^a-1$. Of the 25 cases only the five below yield parameter triples. Together they give 2,245. Exactly 269 admit the required unitary divisor, 264 at $(2,3)$ and five at $(2,5)$. Every one fails either the unitary divisor condition $3\le U\mid\gcd(w,W+1)$, with $U\mid N$ and $\gcd(U,N/U)=1$, or the required carry freedom of $wW,W+1$ at $p$ or $q$.

| $p$ | $q$ | $s$ | Parameter triples |
|---:|---:|---:|---:|
| 2 | 3 | 11 | 636 |
| 2 | 5 | 7 | 1,487 |
| 2 | 11 | 13 | 23 |
| 2 | 17 | 19 | 7 |
| 3 | 11 | 13 | 92 |

For the all-cofactor conclusion, put $H=10^{334}$ and $P,Q\le H$. In case (1), $N<1+H^2(2H+1)/7^3<10^{1000}$. Otherwise $E\le ((s-1)/3)^{3/2}<s^2$ and $P,Q<s^4 E$, giving $N<1+s^{15}<10^{51}$. Both ranges are covered. The same reasoning shows that (1) is necessary for **every** possible exception at these base pairs, without any bound on the exponents.

## Bases 2 and 3 through $10^{10000}$

In this case $s=11$, $U=5$. Apply the [40 necessary cofactor classes modulo 1320](residual_power_classes.md#the-cofactor-classes-at-bases-2-and-3) before looking up the ternary exponent. Exactly 71 tuples remain through $10^{10000}$. A second implementation reconstructs the cofactor from the ternary digits of $2^{-a}$, then looks up $b$ in the binary orbit of 3. It compares the entire set of triples with the binary-prefix search. [Independent reverse search](check_paired_ternary_prefixes.py) At 11 the grouped row requires $W\equiv-1\pmod{11}$. At 5 it requires $W\equiv4\pmod{25}$. Indeed, if $5^h\parallel W+1$, then $10W\equiv-10\pmod{5^{h+1}}$. For $10W,W+1$ to add without a carry the positive residue of $W+1$ at that modulus must be less than ten. Hence $h=1$ and the residue is five.

Only two tuples pass both congruences. The table gives exact valuations of $\binom{11W+1}{W+1}$ by Legendre's factorial formula. The grouped trinomial coefficient is a multiple of this binomial coefficient and must avoid all four primes. Thus the conjecture holds throughout the paired family through $10^{10000}$. Applying the all-cofactor bound with $s=11$ also covers **all cofactors whenever both $2^a$ and $3^b$ are at most $10^{3334}$**.

| $E$ | $a$ | $b$ | $v_2$ | $v_3$ | $v_5$ | $v_{11}$ |
|---:|---:|---:|---:|---:|---:|---:|
| 43 | 1976 | 5473 | 4423 | 619 | 2252 | 1545 |
| 43 | 1976 | 15713 | 12540 | 619 | 5851 | 3892 |

## Counting across all cofactors

For each fixed pair $p,q$, the number of possible exceptional weights through $X$ under the two-layer hypotheses is $O_{p,q}((\log(2X))^2)$, allowing arbitrary cofactors and both exponents to vary. The two-power sum argument gives only finitely many $s$ for the fixed bases. Fix one and write $H=\log(2X)$, so $a,b=O(H)$. At a given $E$, the first congruence restricts $b$ to one class modulo the order of $q$ at $L$. This order is $\gg_{p,q}L$ because lifting at odd $p$ gives $p^{v_p(q^j-1)}\ll_{p,q} j$ whenever $p\mid q^j-1$, and the same bound holds at 2 by the odd and even cases of binary lifting. Likewise the period for $a$ is $\gg_{p,q} M$. Since $L,M\ge s^2E$, the sum over $E\le H$ is $O(\sum_{E\le H}(1+H/E)^2)=O(H^2)$. For $E>H$, each fixed $b$ specifies at most one candidate $E$ per interval $[A^{j-1},A^j)$ and therefore $O(H)$ in total. Each such $E$ admits $O(1)$ possible values of $a$. Summing over $b$ finishes the count.

Run `python3 check_paired_cofactor_range.py 1000` and `python3 check_paired_cofactor_range.py 10000 --binary-ternary` to reproduce the two weight ranges. [Verifier](check_paired_cofactor_range.py)

The optional independent search uses `python-flint`. Run `python3 check_paired_ternary_prefixes.py` in an environment with that package installed.
