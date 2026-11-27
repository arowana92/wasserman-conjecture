# The first cofactor where every finite digit depth can miss the assigned carries

Consider the [bounded-digit endpoint test](cofactor_index_principle.md) at weights $N=1+dq^e$, with $q>d+1$ prime. A possible counterexample must contain the endpoint $(1,1,N-2)$ and two rows of the specified shapes, with one row avoiding $q$ and each full weight prime power assigned to an avoiding row.

**For every $d\le246$ the local tests eliminate every candidate, uniformly for all positive exponents. The first cofactor allowing a row shape that survives every finite depth is exactly $247$.** Through $d=247$ the only such shape, up to the ordering of the two unshifted entries and the companion entries, is
$$
 P=q^e,\quad A=(114P,117P,16P+1),\quad C=(V,V,V),\quad 3V=1+247P.
$$
It survives arbitrarily long simultaneous local tests exactly when **none of $2,3,13,361$ divides $e$**. In particular, the finite divisibility list of the general principle is $\{1\}$ for $D\le246$ and $\{2,3,13,361\}$ for $D=247$. This exact list concerns the stated local tests. The [primitive denominator bound](shifted_directions.md) proves the full conjecture at every $1+247q^e$ for every prime $q$ and positive exponent.

For each surviving exponent, any prime cutoff $B$ and any depth $k\ge1$, a reduced progression gives prime bases $q>\max(B,248)$ for which no prime through $B$ has a carry in all three displayed rows (including the endpoint) through the first $k$ digits. Nevertheless 2 is always common. Its first carry in $C$ occurs at **exactly digit $k+1$**, numbering the units digit as 1.

## Exhausting the smaller cofactors

Write the two candidate rows as $(aP,bP,cP+1)$ and $(fV,gV,hV)$, where $a+b+c=d$, $U=f+g+h$, $UV=1+dP$, $U\mid a,b$ and $(U,d)=1$. All slopes are positive. At a full prime power $w\parallel d$ the residue of the weight is 1. An avoiding row has three residues $0,0,1$ modulo $w$. Since $P$ is a unit there, the first row can supply this pattern only if the two unshifted slopes are both 0 modulo $w$, or one is 0 and the other a unit. Otherwise at least two companion slopes must be multiples of $w$, because $V$ is a unit there.

Two programs enumerate the positive slopes with those constraints, one starting at $U$ and one at $a,b$. They agree on every tuple for every $3\le d\le247$, giving 607 tuples in total. At $e=1$ the scalar at each prime may be any unit. For each of 606 tuples the finite-state search finds an excluding prime, and a separate calculation confirms it by lifting **all** surviving unit residues one digit at a time. The latter checks the sums of the three actual residues at each layer, preserving the choice of avoiding row along a prefix. All 606 exclusions need at most 24 digits at the chosen prime. The initial $0,0,1$ constraints also require at most 7 digits here. At $d=1,2$ there are no positive slope triples. Higher exponents only restrict the available unit scalars. These checks prove the lower bound and uniqueness.

After the initial residue conditions the eight primes $2,3,5,7,13,17,23,47$ suffice for all 606 exclusions, still within 24 digits. This set is necessary for the remaining local test at $e=1$. For each of its eight primes the verifier finds a tuple for which that is the only prime with an exclusion. To check this last fact it tests the primes of $dU(d+1)$ for a carry-free infinite continuation. Every other prime admits scalar 1 and misses the endpoint. This is a minimal set of primes for excluding the enumerated shapes under their assignments.

## The four exponent obstructions

At 3 the first row is assigned the weight prime, so $P\equiv2\pmod3$. This rules out all even exponents. At 13 and 19, which divide $247$, the three equal companion residues cannot be $0,0,1$. The first row must avoid both primes. The required residues are $P\equiv4\pmod{13}$ and $P\equiv13\pmod{19}$. The first is not a cube modulo 13, excluding multiples of 3.

If $13\mid e$, the only possible lift of 4 modulo $13^2$ is 147. The first-row residues there are $27,130,156$, whose sum exceeds 169. If $19^2\mid e$, the only possible lift of 13 modulo $19^3$ is 5903. The residues $760,4751,5282$ sum to more than 6859. The uniqueness in both cases follows by raising the mod-$r$ residue to $r$ or $r^2$ respectively. Such powering gains one digit of precision at each step. Both residue sums force a carry.

## Local scalars proving the converse

Suppose $2,3,13,361$ all miss $e$. At 2 take the unit scalar $P=-1/247$ locally. This makes $V=0$ in the binary test. To prescribe the exact first carry instead, take $1+247P\equiv2^k\pmod{2^{k+1}}$. Every odd residue is locally an $e$th power since $e$ is odd. The actual $V$ then has binary valuation $k$. The endpoint and $A$ both carry in the units column, and the three identical entries of $C$ first carry at column $k+1$.

The following scalars handle the other forced primes.

| Prime | Local scalar $P$ | Check |
|---:|---:|:---|
| 3 | $-16/121$ | $v_3(1+247P)=1$. The first row adds without a carry at every digit. |
| 13 | $43$ | $A(43)$ adds without a carry in base 13. |
| 19 | $127$ | $A(127)$ adds without a carry in base 19 and $127^{18}\equiv1\pmod{361}$. |

At the ternary scalar, the first four triples of digits of $A$ are $(0,0,0),(1,0,0),(0,2,0),(0,0,1)$. The repeating block is $(0,1,0),(0,0,1),(0,2,0),(2,0,0),(0,1,1)$. This can be checked by subtracting the units digit from each rational entry and dividing by 3. The verifier checks the return to the identical three rational remainders. Since $e$ is coprime to 6 and to $12\cdot13$, all units are locally attainable at 3 and 13. At 19 one has $(e,18)=1$ and $v_{19}(e)\le1$, so the congruence modulo 361 makes 127 attainable there as well.

At 31 use $q\equiv-1\pmod{31}$, since $e$ is odd and $31\nmid1-247$. At other primes through $B$ outside $2,3,13,19,31$ use $q\equiv1$. Such primes divide neither 247 nor 248, so they miss the endpoint. Choose root residues for the four forced primes to the required precisions and combine by the Chinese remainder theorem. Dirichlet's theorem gives arbitrarily large primes in the resulting reduced progression. For each of them the rows are positive and integral, and $A$ avoids $q$ by shifting its slope digits. This proves the sharpness assertion.

Run both programs with Python 3.9 or later. The first also accepts a positive exponent and a cofactor limit as command line options.

```sh
python3 check_cofactor_index_depth.py
python3 check_cofactor_index_depth_cross.py
```
