# Individual weights with two bounded predecessor powers and an unrestricted third

Let $w=2^j$, $j\ge1$, let $P$ be an odd prime power and put $D=wP$. For **every prime power $Q$**, the coefficients of any three positive trinomial rows of weight $1+DQ$ have a common prime divisor, under any of the following conditions:

* $D\le15{,}423$.
* $w\in\{8,16\}$ and $D\le40{,}000$.
* $w\le4$, with arbitrary $P$.

The conclusion concerns each individual weight, and permits any base or exponent for $Q$. The last case follows from the [small predecessor-power theorem](three_factors_small_power.md). If $Q$ shares a prime with $wP$, the [two-factor theorem](all_weights_endpoint_reduction.md#consequences) applies. Assume henceforth the three bases are distinct.

## An exact arithmetic progression for the third power

The [three-factor assignment](endpoint_factors.md#when-n-1-has-three-distinct-prime-factors) gives a single row $A$ avoiding all three predecessor bases and all primes in a unitary factor $U\ge3$ of $N$. Since $N$ is odd, choose an odd prime $u\mid U$. Order the exceptional positions by $w,P,Q$. For some positive $a<w$ and $b<P$, the parts are
$$
 A=(aPQ,\;bwQ,\;cQ+1),\qquad c=(w-a)P-bw>0.
$$
In fact $aPQ\equiv1\pmod w$ and $bwQ\equiv1\pmod P$. The last part is a positive multiple of $U$, so it cannot equal one. All parts are multiples of $u$, giving $u\mid a,b$. The distinguished congruences and $u\mid wPQ+1$ fix a unique residue $R$ of $Q$ modulo $M=wuP$:
$$
 aPR\equiv1\pmod w,\qquad bwR\equiv1\pmod P,\qquad wPR\equiv-1\pmod u.
$$
If $a=u=w-1$, then $U\mid a$ and $U\ge3$ imply $U=u$, excluded by the [one-below-a-power carry](three_factors_small_power.md).

The same row must avoid the primes $2,u,p$ where $P$ is a power of $p$. Two short tests discard many assignments. The digits of $aR\bmod P$ must be contained in the digits of $wR\bmod P$ in base $p$: the distinguished part $aPQ$ is shifted by $P$ and the total above the shift is $wQ$. After dividing the entire row by $u$, the units digits of its first two parts are $(a/u)PR$ and $(b/u)wR$ modulo $u$. Their least residues have sum below $u$.

## Carry languages and the prime-power obstruction

For the slopes $(aP,bw,c)$ and offsets $(0,0,1)$, read digits of $Q$ in a fixed base $s$. The standard multiplication carries give a finite graph of the inputs whose three entries are carry-free. After the prescribed low residue modulo $w$, compute the [smallest arithmetic lattice](denominator_residue_hulls.md#computing-the-strongest-single-congruence) containing the binary inputs. Intersect it with $R+M\mathbb Z$.

If this intersection is $L+T\mathbb Z$, the positive integer $\gcd(L,T)$ divides $Q$. Two distinct divisors which are primes exclude any prime-power $Q$. If a prime $q$ divides this gcd, a possible $Q$ must be a power of $q$. Thus the three slopes themselves have to add without a carry in base $q$: multiplying them by $Q$ shifts all their digits, and the added unit lies below the shift. The verifier also handles singleton and empty languages.

For example, at $D=32\cdot383$ a surviving preliminary choice was $u=7,a=21,b=28$. Its slopes are $(8043,896,3317)$. With the fixed low residue 3 modulo 32 the accepted binary inputs are all $163\pmod{42{,}729{,}472}$. The modulus is divisible by 163. The actual positive word $Q=163$ passes that binary test, but the slope residues modulo 163 are 56, 81 and 57, whose sum is 194. This excludes the whole assignment, regardless of the exponent of $Q$.

For the remaining possibilities the verifier tests carry freedom along the full progression. It reads either the digits of $k$ in $Q=R+Mk$ or the high digits of $Q$ after its forced low power of the base. In the second method, traverse edges backwards from all stopping states with high-part residue zero. An edge labelled $d$ transforms a suffix residue $h$ into $sh+d$ modulo the remaining coprime modulus. Failure to reach the target residue at the initial state excludes the assignment. A search stopped at an implementation resource limit leaves that particular test inconclusive. Every shape in the claimed ranges has a conclusive exclusion.

| First conclusive test | Assignments excluded |
|---|---:|
| Fixed digit at $p$ or $u$ | 597,562 |
| No accepted binary input after the prescribed low digits | 229,583 |
| Binary lattice incompatible with the progression | 116 |
| Forced divisor of a prime power | 22 |
| Direct progression graph | 22 |
| Suffix residues | 7 |
| **Total** | **827,312** |

A second enumeration starts with $D$, determines its full power of 2 and tests the odd part for being a single prime power by trial division. It enumerates integer pairs $a,b$ by the positive-slope bound, then factors $\gcd(a,b)$ to recover every candidate $u$. It checks the complete assignment list and all congruences against the primary generator. Exact factorial-floor positive controls verify accepted words, including $Q=163$ in the example. Neither enumeration truncates $Q$ in the proof.

~~~~sh
python3 check_even_two_power_cofactors.py
python3 check_even_two_power_cross.py
~~~~
