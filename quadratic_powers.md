# Quadratic polynomials at prime powers

Let $h,d$ be positive integers, $P=p^a$ with $p$ prime and $a\ge1$, and $\epsilon\in\{0,1\}$. Put $N=P(hP+d)+\epsilon$, $K=\max(2,h,d)$ and $\Delta_\epsilon=d^2+4h(1-2\epsilon)$.

**Theorem.** Either of the following hypotheses suffices for Wasserman's three-row conjecture at every such weight with $P>8K^{15}$.

1. $\gcd(h,d)\le2$ and $\Delta_\epsilon$ is not the square of an integer. Negative discriminants qualify.
2. At least one of $h,d$ is 1 or a positive power of $p$.

For **$1\le h,d\le15$**, the conclusion under the same hypotheses holds for **every $P$**, including the initial powers. For example, $13P^2+15P$ and $13P^2+15P+1$ both satisfy the conjecture for every prime power $P$, since their discriminants are 277 and 173. The proof below gives sharper effective bounds for each pair of coefficients.

## The three linear polynomials

Assume $P>d$. By the [endpoint reduction](all_weights_endpoint_reduction.md), assign all primes of $N(N-1)$ to the two other rows of a counterexample, each prime to a row avoiding it. Write $A$ for the row taking $p$. Kummer's criterion writes its parts as $P L_i(P)+\delta_i$, where
$$L_i(X)=a_iX+b_i,\qquad \sum a_i=h,\quad\sum b_i=d,\quad\sum\delta_i=\epsilon.$$
Indeed $hP+d$ has disjoint high and low blocks in base $p$ and both blocks split without carries. If $\epsilon=1$ label the shifted slot 3. Each $L_i$ is nonzero. Otherwise a part is zero or one. A part equal to one would mean the scale assigned to $A$ from $N$ divides 1, leaving the complementary row with three positive normalized parts of sum one.

For each line set $D_i=h b_i-d a_i$. Write $M_{ij}=a_i b_j-a_j b_i$ and
$$E_i=\begin{cases}|h b_i^2-d a_i b_i+(2\epsilon-1)a_i^2|,&a_i>0,\\ b_i,&a_i=0.\end{cases}$$
If an integer divides both $L_i(P)$ and $hP^2+dP+(2\epsilon-1)$, it divides $E_i$. Multiply the quadratic by $a_i^2$ and substitute $a_iP\equiv-b_i$ modulo the common divisor. Likewise any common divisor of $L_i(P),L_j(P)$ divides $M_{ij}$, and a common divisor of $L_i(P),hP+d$ divides $D_i$.

## The unshifted weight

Let $N=P(hP+d)$. Write $Pv$ for the product of the full powers of $N$ assigned to $A$, and $S$ for those assigned to it from $N-1$. Then $v$ divides $hP+d$ and all three $L_i(P)$, including when $p$ divides $v$. Group the full powers making up $S$ by their exceptional slot, where the residue is 1. The group at $k$ divides the two other $L_i(P)$ since it is coprime to $P$. These three groups are pairwise coprime. Thus
$$v\mid v_* = \gcd(D_1,D_2,D_3,M_{12},M_{13},M_{23}),\qquad
S\mid S_* = \operatorname{lcm}_{i<j}\gcd(E_i,E_j,M_{ij}).$$
The complementary row normalized by its share of $N$ has sum $Pv$. Two of its parts are divisible by each full power in $(N-1)/S$, so
$$27(N-1)^2\le (Pv)^3 S^2,\qquad 27h^2 P<v^3 S^2. \tag{1}$$

## The shifted weight

Let $N=P(hP+d)+1$, write $Q=hP+d$, and let $U$ be the scale assigned to $A$ from $N$. It divides both unshifted $L_1(P),L_2(P)$ and satisfies $U\ge3$. Therefore
$$U\mid U_* = \gcd(E_1,E_2,M_{12}).$$
Write $Q=BT$, where the full powers assigned to $A$ from the predecessor have product $PB$. Split $B$ into three coprime groups $B_k$ by exceptional slot. At a full power $w$ in $B_k$ one has $P L_i(P)+\delta_i\equiv [i=k]\pmod w$. Using $hP\equiv-d\pmod w$ gives
$$B_k\mid h^2(\delta_i-[i=k])-dD_i\qquad(1\le i\le3).$$
Also $B_k\mid D_i$ when $i$ is unshifted and $i\ne k$, and for all $i$ when $k=3$. In these cases $w\mid L_i(P)$ by cancelling $P$ if $w$ is coprime to $p$. At the base $p$, use the entire power $Pw$ in $PQ$ and cancel $P$ from the corresponding congruence at that larger modulus.

Let $b_k$ be the gcd of the listed integers for slot $k$, and $B_*=\operatorname{lcm}(b_1,b_2,b_3)$. Then $B\mid B_*$. The normalized complementary row now has sum $U$ and two multiples of each full power in $T$. We obtain
$$27(hP+d)^2\le U^3B^2. \tag{2}$$

## Why the bounds are finite

Under hypothesis 1 all $E_i$ are positive. A vanishing $E_i$ for a nonconstant line would make $-b_i/a_i$ a rational root of the quadratic, requiring its integral discriminant to be a square. Put $g=\gcd(h,d)$ and $(h',d')=(h/g,d/g)$. If every $D_i$ vanished, the three pairs $(a_i,b_i)$ would be positive integer multiples of $(h',d')$, whose multipliers sum to $g\le2$. Thus $v_*>0$, and $b_3>0$ by the same argument.

If $k\ne3$ and $b_k=0$, the neutral slot different from $k,3$ has $D_i=0$ and uses at least $(h',d')$. This is impossible for $g=1$. For $g=2$ it uses exactly $(h',d')$, leaving total $(h',d')$ in the other slots. At the shifted slot the vanishing equation gives
$$d'(h'b_i-d'a_i)=h'^2.$$
Coprimality forces $d'=1$. The remaining coefficients satisfy $|h'b_i-a_i|\le h'$, so $h'=1$ too. Equality then forces the other two lines to be $1$ and $X$, one of them unshifted. Both evaluate to units modulo $U$, contradicting $U\ge3$. Hence every shape that could give a counterexample has $B_*>0$.

Under hypothesis 2 the single digit unit in $h$ or $d$ goes wholly to one slot. No individual $D_i$ can vanish. A nonzero line proportional to $(h,d)$ using that unit would consume both totals. Thus every $b_k>0$ by the unshifted, nonexceptional slot. At most one line has both coefficients positive, so at most one $E_i$ can vanish. This still gives $U_*>0$ and $S_*>0$.

For an explicit uniform threshold, $|D_i|,|M_{ij}|\le K^2$ and $E_i\le3K^3$. If all $E_i$ are positive, $S^2\mid E_1E_2E_3$ because every full power in $S$ divides two of them. If one $E_i$ vanishes under hypothesis 2, each power divides at least one of the remaining two, both at most $K^2$. Thus $S^2\le27K^9$ in either case, and $v\le K^2$. Equation (1) gives $P<K^{15}$. In (2), $U\le3K^3$ and $B\le(2K^3)^3$, since each nonzero group bound is at most $K^3+K^2$. This gives $P\le8K^{27/2}$. Together these prove the claimed threshold.

## Checking all initial powers for coefficients through 15

For the finite box, $P\le d$ is already inside the [verified interval through one billion](full_range_1b.md). Enumerate the three coefficient pairs of the lines, ordering all slots when $\epsilon=0$ and the two unshifted slots when $\epsilon=1$. For pairs covered only by hypothesis 2, keep shapes splitting $h,d$ without a carry in an eligible base. The bounds on $v,S,U,B$ depend only on the shape.

For $\epsilon=0$, enumerate $v\mid v_*$ and $S\mid S_*$ coprime. Equation (1) ends the progression at $\lfloor(v^3S^2-1)/(27h^2)\rfloor$. Enumerate $P$ modulo $v$ by vanishing of all three $L_i(P)$ and modulo $S$ by $P(hP+d)\equiv1$ together with $S\mid\operatorname{lcm}_{i<j}\gcd(L_i(P),L_j(P))$. Combine these congruences, retain prime powers and unitary assignments, and test the first row in every assigned base. For $\epsilon=1$, (2) gives a short interval of $P$ to check directly, using the divisors of $U_*$ and the full powers in $Q$.

| Computation | $\epsilon=0$ | $\epsilon=1$ |
|:---|---:|---:|
| Coefficient pairs covered at every prime base | 161 | 159 |
| Additional pairs covered at selected bases | 50 | 52 |
| Shapes checked | 64,180 | 181,975 |
| Shapes whose upper bound reached beyond the initial interval | 915 | 10 |
| First rows passing the assignment and carry tests above the initial interval | 0 | 0 |

The unshifted congruences produced 36,266 progression points before primality and full-power assignment tests. This proves the finite-box result for all exponents.

Run the following with Python 3 and python-flint.

```sh
python3 check_quadratic_powers.py
python3 check_quadratic_powers_cross.py
```

The independent small-weight control enumerates normalized parts directly and computes each coefficient as two exact binomials. Taking every prime avoided by that coefficient gives the maximal first-row assignments. All 3,202 unshifted and 9,341 shifted carry-free rows tested satisfy the divisibilities for those assignments, including nontrivial controls for each of the four bounds.
