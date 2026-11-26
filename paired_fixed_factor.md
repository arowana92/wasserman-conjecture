# Four of five weights when one prime-power factor is at most 500,000

**Theorem.** Let $P,Q$ be full powers of two distinct primes, with $P<Q$ and $P\le500{,}000$. For every such $Q$, at least four of
$$5PQ+1,\quad7PQ+1,\quad8PQ+1,\quad9PQ+1,\quad16PQ+1$$
satisfy Wasserman's three-row conjecture. There is no upper bound on $Q$ or on its exponent. The [five-weight theorem](five_related_weights.md) already covers all five-weight inputs through $10^{1000}$ in the sense of this four-out-of-five conclusion. When the two underlying primes coincide, all five predecessors have at most two distinct prime factors.

## Equality of parts at different weights

Suppose two distinguished rows have a pair of complementary parts summing to the same $C=dt+1$ at both weights, where $t=PQ$. Match a part $x$ in one pair with $x'$ in the other by their residues at $P,Q$, namely $(1,0)$ or $(0,1)$. Thus $t\mid x'-x$. Suppose each pair adds without carries in each prime base $r$ in a set $S$. Put $g=\prod_{r\in S}r^{v_r(C)}$. At a layer $r^{v_r(C)}$ dividing $C$, the residues of both parts must be zero, since they add to the residue of $C$ without overflow. Hence $g\mid x,x'$, and $\gcd(g,t)=1$. **If $g\ge d$ and all parts of the two pairs exceed one, the two pairs coincide.** Indeed $|x'-x|<dt\le gt$, while $gt$ divides the difference. Equality of the second parts follows by subtraction.

After this identification the common pair must avoid every prime base avoided by either of the two rows. Let $G$ be the full part of $C$ supported on this union of primes. Then **$G(P+Q)\le C$**: label the shared parts to be multiples of $P$ and $Q$, respectively. Both are multiples of $G$, which is coprime to $PQ$. Lean checks the layer divisibility, the equality under the range condition, and the scale bound. [Lean lemmas](SharedParts.lean)

For the 9/16 pair, $C=5t+1$ and $16\mid C$ since $t\equiv1267\pmod{1584}$. Both rows avoid 2. The criterion above forces their complementary pairs to agree, and the common pair avoids 2, 3 and 11. Thus the full $\{2,3,11\}$-part of $5t+1$ times $P+Q$ is at most $5t+1$. In particular each prime power is greater than $1584/5$.

## Fixing one factor removes the numerical cutoff on the other

Only simultaneous failures at multipliers 7/16 and 9/16 need to be excluded. The [shared-part restrictions](five_related_weights.md#multipliers-7-and-16) turn the rows into affine expressions in $Q$. In each row, place the part exceptional modulo $P$ in its second position, where it is a multiple of $Q$ and congruent to 1 modulo $P$.

* For 7/16 the two rows are $(6PQ,cQ,(P-c)Q+1)$ and $(3PQ,eQ,(13P-e)Q+1)$. Here $0<c<P$, $42\mid c$, $48\mid e$, and $e=c$ or $c+12P$. Indeed the two matched shifts are 0 and $12t$, and the two row pairs of nondistinguished parts are multiples of 42 and 48, respectively. The residue class $t\equiv587\pmod{1008}$ ensures $\gcd(42,Q)=1$ and $\gcd(48,Q)=1$. Thus $c\equiv0$ modulo 336 in the first case and $c\equiv84P\pmod{336}$ in the second. The second congruence follows by writing $c=42h$ and solving $42h+12P\equiv0\pmod{48}$.
* For 9/16 the rows are $(4PQ,cQ,(5P-c)Q+1)$ and $(11PQ,cQ,(5P-c)Q+1)$, with $0<c<5P$ and $1584\mid c$. Here $t\equiv1267\pmod{1584}$ ensures $\gcd(Q,1584)=1$.

In both cases $\gcd(c,P)=1$ and $cQ\equiv1\pmod P$. All six parts at the two weights must add without binary carries within each row. For each fixed $P,c,e$ this is a finite-state test on the digits of $Q$, with no limit on word length. For the six slopes $s_i$ and the offsets $b_i=0,0,1,0,0,1$, first try each low byte $q$ allowed by the two binary cuts and the necessary $PQ\bmod16$. Retain the byte if the residues of $s_iq+b_i$ modulo 256 are bitwise disjoint within each row. Their quotients by 256 form its state. Subsequent bits $d$ produce the next output bits $(s_i d+u_i)\bmod2$ from carries $u_i$. Accept a state if the three remaining carries in each row are pairwise disjoint. All states reached from the permitted bytes are explored until acceptance or exhaustion. Zeros may be appended at the high end of a finite input.

| Pair | Shapes tested | Shapes accepting any binary word | All accepted values of $Q$ |
|---|---:|---:|---|
| 7/16 | 59,125,413 | 47 | 1 for 40 shapes, 49 for 7 shapes |
| 9/16 | 31,344,480 | 1 | 1 |

None of these values is an admissible partner. In the seven shapes accepting $Q=49$ one has $P=456523$, and the 7/16 residue class also excludes $7\mid Q$.

The C++ generator enumerates all prime powers through the stated bound by a sieve. Python independently enumerates every shape through $P=5000$ by trial division. It also reconstructs the complete accepted language at each of the 48 surviving shapes through $500000$, starting from four low bits. To obtain canonical words, a nonempty high suffix must end with a one and stop immediately. The verifier discards states unable to reach such a final one and checks that the remaining paths from the initial state are acyclic. Enumerating these paths gives the exact integer values in the table. This distinguishes arbitrarily long zero padding from genuinely larger partners.

~~~~sh
mkdir -p .work
c++ -std=c++17 -O3 check_paired_fixed_factor.cpp -o .work/paired_fixed_factor
python3 check_paired_fixed_factor.py .work/paired_fixed_factor
lean SharedParts.lean
~~~~
