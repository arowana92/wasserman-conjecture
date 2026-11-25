# Denominators of rows at arbitrary scales

In any counterexample to Wasserman's three-row conjecture, take either row other than the endpoint $(1,1,N-2)$. For **every** ordering $(x,y,z)$ of this row, set
$$W=\gcd(x,y,z-1),\qquad m=(N-1)/W.$$

**Theorem.** Always $m\ge364$. If $m\le1000$, then $W>10^{1000}$, $W$ is neither a square nor a 105th power, and the only possible denominators are
$$364,380,589,629,646,703,710,738,780,850,874,880,946,950,962,986.$$
A cube $W$ requires $m=880$ or $986$. For every $\epsilon>0$, the number of possible counterexample weights $N\le X$ having any such $m\le1000$ is $O_\epsilon(X^\epsilon)$. The scale $W$ can have arbitrarily many prime factors and can share factors with $m$.

**Consequence for weights.** For every prime $p$, every $e\ge1$, and every $1\le d\le1000$, the full conjecture holds at both $d p^{2e}+1$ and $d p^{105e}+1$. It holds at $d p^{3e}+1$ under the same hypotheses except possibly when $d=880$ or $986$. These conclusions include $p\mid d$.

## The arithmetic that bounds the enumeration

The [all-weight endpoint theorem](all_weights_endpoint_reduction.md) supplies the endpoint row. Assign every full prime power of $N$ and $N-1$ to one of the two remaining rows that adds without carries in that base. Write $UV=N$, with $U$ assigned to the selected row $A$ and $V$ to the other row $C$. Their parts are divisible by $U,V$, respectively, and $U,V\ge3$. The integer $W$ divides $N-1$, so it is coprime to $U,V$. We obtain
$$A=(aW,bW,cW+1),\quad a+b+c=m,\quad \gcd(a,b,c)=1,\qquad C=V(f,g,h),\quad f+g+h=U.$$
All six multipliers are positive since the shifted part is divisible by $U\ge3$. We may order $a\le b$ and $f\le g\le h$. Also $U\mid a,b$ and $W\equiv-m^{-1}\pmod U$.

Let $w=r^j$ be a full power in $m$. If $C$ takes $r$, then two of $f,g,h$ are divisible by the entire $r$-part of $mW$, because the actual residues at that power are $0,0,1$. If $A$ takes $r$ and $r\mid W$, only its shifted part can have residue 1. Cancelling the $r$-part of $W$ from each unshifted part gives $w\mid a,b$. Then $r\mid c$ as well, violating primitivity. Thus if $A$ takes $r$, the scale $W$ is a unit at $r$. Both $a,b$ cannot vanish modulo $w$ for the same reason. Exactly one vanishes modulo $w$, the other is a unit modulo $r$, and $(a+b)W\equiv1\pmod w$.

In particular **every full prime power $w$ in the denominator satisfies $3w<m$**, even when $m>1000$. If $A$ takes it, $Uw$ divides $a$ or $b$. If $C$ takes it, $2w+1\le U$ and $m\ge2U+1$. More generally, partition the full powers of $m$ into products $F_1,F_2,D$, placing those assigned to $A$ at the unshifted slot they divide and those assigned to $C$ in $D$. Then
$$U(F_1+F_2)<m,\qquad 27D^2\le U^3.$$
Indeed $D^2\mid fgh$. The entire product of the full powers in $mW$ whose bases divide both $m,W$ also has its square dividing $fgh$.

## Unbounded local tests

For each $m\le1000$, enumerate the positive multiples $a,b$ of $U$ and the positive partitions $(f,g,h)$ of $U$. The routine shared with the [prime-power cofactor calculation](prime_power_cofactors.md) initially allows both unshifted slopes to vanish modulo $w$. This is a superset of the required shapes. Impose primitivity, both possibilities for each $w$ just proved, and the joint congruences on $W$ wherever $C$ cannot take its base. Write the resulting progression, including parity, as $W=R+Mn$. Also $V=(mR+1)/U+(mM/U)n$.

Test every base dividing $2mU$. At a prime of $U$ the row $A$ must avoid that prime and $V$ must be a unit. At a prime of $m$ at least one eligible row must avoid it. Every weight has even endpoint coefficient. If $m,W,U$ are odd then $V$ is even and $C$ takes 2. When $m,W$ are even, $C$ must have two slopes divisible by twice the 2-part of $m$. Two slopes of $C$ with equal 2-adic order always cause a binary carry. The same holds for the slopes of $A$ when $W$ is even, since the added unit then occupies the unused bottom bit.

For two unshifted parts $vT,wT$, the coefficient is divisible by $(v+w)/\gcd(v,w)$. This follows by factoring off their binomial coefficient and using
$$\frac{v}{v+w}\binom{(v+w)T}{vT}=\binom{(v+w)T-1}{vT-1}.$$
Apply this to $a,b$ and all pairs of $f,g,h$. The [affine carry automaton](carry_intersections.md) decides whether each required row can avoid its base for some nonnegative input $n$. The first input digit tests unitarity at the primes of $U$. Failure on the larger initial progression modulo $U$ also rules out the full progression. Each prime can use its own witness at this stage, so every rejection is necessary for an actual counterexample.

Of 41,029 preliminary shapes, 21 pairs of a shape and a parity remain. The checker prints all their rows, progressions and forced bases. The even-scale denominators are $589,629,703,738,986$. Testing all residues modulo the 21 progression moduli gives no square or 105th power. Only the two stated denominators admit a cube. The separate checker verifies this conclusion by taking roots modulo each prime-power factor of each modulus.

## Simultaneous tests and the prime-power consequence

Each surviving shape forces $A$ to avoid at least two specified primes. The two smallest upper bounds for its missing-digit dimensions sum to at most $86/100$. These are integer certificates from length-64 block counts. The [two-base intersection argument](carry_intersections.md#counting-at-unbounded-weights) applies to $W$ itself and gives the stated $O_\epsilon(X^\epsilon)$ bound after summing over the finite list.

Two independent algorithms cover $2\le W\le10^{1000}$ for all 21 progressions. One intersects high-digit intervals of the carry graphs, visiting 745,367 nodes. The other uses 326,658 intervals on which Legendre's factorial floors force a carry. Both return the empty set. At $W=1$ the weight is already in the [verified initial interval](full_range_1b.md).

Finally suppose $N=d p^k+1$, $1\le d\le1000$. Choose $A$ avoiding $p$ and choose its shifted slot by the residues modulo $p^k$. Then $W=t p^k$ with $mt=d$. As $m\ge364$, $t\le2$. If $t=2$, then $m\le500$ and $W$ is even. Both available denominators in this range, 364 and 380, require odd $W$. Hence $t=1$, $W=p^k$, and $m=d$. The tests for perfect powers give the consequence. In particular if $p\mid d$, the full $p$-part of $N-1$ would force $t>1$, so all exponents are excluded in that case.

## Reproduction

Run the following from the repository root with Python and python-flint. The second command independently enumerates 506 primitive shapes at denominators through 100 and at 364, 380, 738 and 986. It compares the local arithmetic with exact multinomials in 6,948 prime tests, including 2,047 cases where the denominator and scale share a prime. The two interval routines also recover 24 positive controls.

```sh
python3 check_shifted_directions.py
python3 check_shifted_directions_cross.py
```
