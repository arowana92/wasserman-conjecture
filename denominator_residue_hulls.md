# Two exclusions from the congruences satisfied by an entire carry language

**In every possible counterexample, for each ordering $(x,y,z)$ of either nonendpoint row,**
$$
 \frac{N-1}{\gcd(x,y,z-1)}\ge380 .
$$
**The quotient cannot equal 880 either.** There is no condition on the common scale $W=\gcd(x,y,z-1)$.

Consequently every weight $1+dq^e$ satisfies the full conjecture for $1\le d<380$, and also for $d=880$, for every prime $q$ and every $e\ge1$. For weights $1+dq^{3e}$ the full conclusion now holds throughout $d\le1000$ except possibly $d=986$. These consequences follow from the [primitive denominator theorem and its prime-power reduction](shifted_directions.md#simultaneous-tests-and-the-prime-power-consequence): a failing weight with $d\le1000$ requires primitive denominator $d$ itself, and the only candidates for a cube had been 880 and 986.

The [classification through 1000](shifted_directions.md) already excludes every denominator below 380 except 364. For 364 and 880 it specifies just the following possibilities. The first row is $(aW,bW,cW+1)$, and the companion is $V(f,g,h)$ with $UV=N$.

| Denominator | $(a,b,c)$ | $U$ | $(f,g,h)$ | Required scale | Base the first row must avoid |
|---:|---|---:|---|---|---:|
| 364 | $(76,247,41)$ | 19 | $(5,7,7)$ | $W\equiv747\pmod{988}$ | 2 |
| 880 | $(272,561,47)$ | 17 | $(2,5,10)$ | $W\equiv2801\pmod{2992}$ | 2 |

At 364 take just $W\equiv3\pmod4$ and the binary carry condition. All such $W$ are congruent to 135 modulo 69,120. In particular $3\mid W$ and $3\mid N-1$. Write $W=3^j K$ with $3\nmid K$. At the first column above the shift the row has digits $76K,247K,41K$ modulo 3. They sum to at least 3 for both unit classes of $K$. The companion also has a carry: $3\nmid V$ and the digits of $5V,7V,7V$ sum to at least 3. The endpoint coefficient $N(N-1)$ and both row coefficients thus have the common prime 3.

At 880 use $W\equiv1\pmod{16}$. Binary carry freedom forces $W\equiv241\pmod{15360}$, in particular $W\equiv1\pmod5$. Since $17V=880W+1$, one has $V\equiv18\pmod{25}$. The residues of $2V,5V,10V$ modulo 25 are $11,15,5$ and sum to 31. This forces a carry at 5. The first row also carries at 5 since its three units digits are $2,1,3$. Again 5 divides $N-1$, so 5 is common to all three coefficients.

## Computing the strongest single congruence

For a finite carry graph in base $p$, let $L_s$ be the set of nonnegative suffix integers readable from state $s$ up to acceptance. If it is nonempty, its smallest containing arithmetic lattice is $h_s+d_s\mathbb Z$, where $h_s$ can be any actual member and $d_s$ is the gcd of all differences of members. The convention $d_s=0$ describes a singleton.

Start each stopping state at $(0,0)$ and all others empty. If an edge with digit $e$ goes from $s$ to a state with summary $(h,d)$, propagate $(e+ph,pd)$ to $s$. To merge into a nonempty summary $(h_0,d_0)$ keep $h_0$ and replace $d_0$ by $\gcd(d_0,pd,e+ph-h_0)$. This computes the hull of the union. The representatives remain actual words. The procedure terminates: once a state changes from empty to a singleton or a positive modulus, every subsequent change after the first positive modulus is to a proper divisor. At stabilization every finite accepted word belongs to its state lattice by induction on its length. Conversely every propagated lattice lies in the true hull, by induction on updates. This proves equality, including for graphs with cycles.

For a fixed low input residue $l$ modulo a power $B$ of $p$, check those low columns and start with the ordinary multiplication carries there. A suffix summary $(h,d)$ yields $l+Bh+Bd\mathbb Z$. The binary computation gives both stated congruences. An independent forward traversal tracks the actual input residue and place instead and confirms that the only residues modulo 135 and 5 are 0 and 1, respectively. The inputs 135 and 241 themselves pass the respective low-bit and binary-row tests.

~~~~sh
python3 check_two_denominator_hulls.py
~~~~
