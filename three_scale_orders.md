# Three independent carry bases and two unrestricted prime powers

Let $t=p^a q^b$ for primes $p,q$ and integers $a,b\ge1$. A binary order of $t$ is the integer $\lfloor\log_2 t\rfloor$. In **every** interval of $H$ consecutive orders, the following bounds count orders containing even one possible failure of the full three-row conjecture.

| Weight | Condition on $t$ | Bound on exceptional orders |
|---|---|---:|
| $5t+1$ | All inputs | $O(H^{0.927})$ |
| $7t+1$ | All inputs | $O(H^{0.980})$ |
| $9t+1$ | $t\not\equiv115\pmod{144}$ | $O(H^{0.940})$ |
| $9t+1$ | $t\equiv115\pmod{144}$, and $(9t+1)/4$ has a prime factor $\le2000$ | $O(H^{0.980})$ |
| $9t+1$ | Same residue, prime factor $\le100000$ | $O(H^{0.9997})$ |
| $9t+1$ | Same residue, prime factor $\le110000$ | $O(H^{0.9999})$ |

The constants are uniform in the starting order. In particular every sufficiently long block contains a run of $\gg H^{0.020}$ consecutive binary orders in which the conjecture holds for both cofactors five and seven, and also cofactor nine except possibly in the stated residue with no prime factor through 2000 in $(9t+1)/4$. The same run can include all [prime-power cofactors through 1000](two_scale_orders.md) by adding their exceptional orders before taking gaps.

## Counting a rotation in the phase plane

Let $b_1,b_2,b_3$ be pairwise multiplicatively independent integer bases and $B_i=b_i^L$ for a fixed $L$. For input $t$, suppose a positive integer $T_i=c_i t+O(1)$ has all its base-$B_i$ blocks in a fixed digit alphabet, including padded initial and zero fractional blocks, where $c_i>0$ and the errors are uniformly bounded. Let $K_i$ be the compact set of infinite expansions in that alphabet. Write $s_i=\log_{B_i}|\mathcal D_i|$, $S=\sum_i s_i<2$ and $\tau=\max(S-1,s_1,s_2,s_3)<1$. For each $\epsilon>0$ the number of occupied base-$b_1$ orders in an arbitrary interval of length $H$ is
$$
                       O_\epsilon(H^{(2+\tau)/3+\epsilon}).
$$
Bounds on the $s_i$ can replace their exact values in this formula.

Normalize $T_i$ to $x_i=T_i/B_i^{n_i}$ with $n_i=1+\lfloor\log_{B_i}T_i\rfloor$. The product $K=\prod K_i$ then contains $x$ in $D=\prod[B_i^{-1},1]$. Put
$$
 \Phi(x)=\big(\log_{B_2}(x_2/x_1),\ \log_{B_3}(x_3/x_1)\big),\qquad
 \theta=(\log_{b_2}b_1,\log_{b_3}b_1).
$$
Modulo $\mathbb Z^2$, $n_1\theta$ is within $O(1/t)$ of $\Phi(x)$ translated by the two constants $-\log_{B_i}(c_i/c_1)$. This follows by substituting $x_i/x_1=(T_i/T_1)B_1^{n_1}/B_i^{n_i}$. It suffices to count a window of $O(H)$ consecutive $n_1$; each occupied $n_1$ gives only $O(1)$ ordinary orders of $t$.

We need two local covering facts. All constants may depend on the fixed bases and alphabets.

1. The inverse image in $D\cap K$ of a disk of radius $w\le1$ under $\Phi$ meets at most $O_\xi(w^{-\xi})$ grid cubes of side $w$, for every $\xi>0$. If nonempty, this inverse image lies in an $O(w)$-tube about a line $(u,\beta u,\gamma u)$ whose direction has all coordinates bounded away from zero after normalization. The relevant ratios $\beta,\gamma$ are in compact positive ranges since $D$ is bounded away from the coordinate planes. In the normal two-plane every unit vector thus has at most one coordinate of magnitude less than a fixed small constant. For the natural product measure on $K$, the proof of Theorem 1.15 of Corso and Shmerkin gives a uniform bound $O_\xi(w^{S-\xi})$ on the mass projected onto any disk of radius $O(w)$ in these normal planes, since $S<2$. Every ball of radius $w$ centered on $K$ has mass $\gg w^S$. Packing disjoint balls in an enlarged tube proves the count. [Corso and Shmerkin, proof of Theorem 1.15, especially the projection estimates after (6.3)](https://arxiv.org/pdf/2409.04608)
2. Fix any cube of side $w\le1$ and a line in the phase plane. The points of $D\cap K$ in that cube for which $\Phi(x)$ is within $h\le w$ of the line meet $O_\xi((w/h)^{\tau+\xi})$ cubes of side $h$, uniformly in the initial cube and the line. Here is the local argument. Write the line condition as $|g(x)-c|\le h$ with $g=v\cdot\Phi$ for a unit normal $v$. Its gradient has the form $(-(a+b)/x_1,a/x_2,b/x_3)$, with $a=v_1/\log B_2$ and $b=v_2/\log B_3$. Uniformly on $D$ at least two components have magnitude bounded away from zero, and all second derivatives are bounded on a neighborhood of $D$. In a cube of side $R$ the product is covered by boundedly many digit rectangles with sides comparable to $R$; these are translated diagonal images of $K$. After normalization to a unit cube their scale factors lie in compact positive intervals. For any sequence with $R\to0$, $h/R\to0$, each subsequential Hausdorff limit of the normalized pieces lies in a finite union of these product copies intersected by a single plane. Indeed Taylor's remainder after normalization is $O(R)$ and the normalized tolerance is $O(h/R)$. Pass to a convergent subsequence of normals, offsets and product copies, using a point of each nonempty piece to bound the plane offsets. The limiting plane has at least two nonzero normal coordinates, also after the diagonal rescaling. If three are nonzero its intersection has upper box dimension at most $\max(S-1,0)$. With only two nonzero the bound is $s_k+\max(s_i+s_j-1,0)\le\tau$ by the two-coordinate theorem, with coordinate $k$ free. [Corso and Shmerkin, Theorem 1.15](https://arxiv.org/pdf/2409.04608)

   For completeness this limit statement gives a bound down to the thickness $h$. There exist $R_0,A>0$ and an integer $m_0$ so that every nonempty dyadic piece of side $R\le R_0$ with $R\ge Ah$ has at most $2^{m(\tau+\xi)}$ occupied children at **some** relative dyadic depth $1\le m\le m_0$. Otherwise choose a sequence failing at all depths through a growing bound, with $R\to0$ and $h/R\to0$. The limit just described has upper box dimension at most $\tau$. At some sufficiently large fixed depth it is covered by $O(2^{m(\tau+\xi/2)})$ cells, and adjoining a bounded number of neighbors covers the convergent pieces, a contradiction. Start with boundedly many dyadic cubes comparable to $w$ (if $w>R_0$ subdivide first at the fixed scale $R_0$). At each stage use the available depth until reaching scale below $Ah$. The sum of the $(\tau+\xi)$-powers of side lengths never increases. Terminal sides are comparable to $h$, which proves the bound, also when $w$ is already comparable to $h$. This argument allows finite thickness and the initial cube can be arbitrarily small. [Earlier exact-line argument](three_base_orders.md#a-uniform-bound-for-straight-lines-in-the-phase-plane)

Now take $Q=\lfloor H^{1/3}\rfloor$. By placing $Q^2+1$ orbit points in $Q^2$ squares, there is a step $1\le m\le Q^2$ whose torus displacement has sup norm $0<\rho\le 1/Q$. For some constant $K_0$ we have $\rho\gg m^{-K_0}\gg H^{-K_0}$ by the two-logarithm estimate for the first coordinate. Discarding $t\ll 1/\rho$ costs $O(1+\log H)$ initial ordinary orders. All remaining phase errors are $O(\rho)$. [Two-logarithm bound and its hypotheses](logarithmic_intervals.md#two-leading-prefixes-and-a-rotation)

Split the $H$ indices into progressions modulo $m$, then into segments of at most $\lfloor1/\rho\rfloor$ points. A segment of $M$ points lifts to a straight segment with spacing $\rho$ and length at most $w=M\rho\le1$. The allowed product points map within $O(\rho)$ of its line and within a disk of radius $O(w)$. Boundedly many integer translates suffice on the torus. The first cover leaves $O_\xi(w^{-\xi})$ cubes at scale $w$, and the second cover leaves $O_\xi(M^{\tau+\xi})$ small cubes per first cube. Each small cube projects into a region of diameter $O(\rho)$ by Lipschitz continuity, hence can account for $O(1)$ indices in this segment even with the approximation error. Fixed factors in the radii change only the constants. Since $w\ge\rho\gg H^{-K_0}$, the loss factors $w^{-\xi}M^\xi$ are bounded by $O(H^{(K_0+1)\xi})$.

There are $J=O(m+H\rho)=O(H^{2/3})$ segments in all. Their lengths add to $H$, so concavity gives $\sum M^\tau\le H^\tau J^{1-\tau}$. Take $\xi$ sufficiently small for the desired $\epsilon$. This proves the counting theorem at every starting index, including when the two components of $\theta$ satisfy a joint rational relation.

## Which digit restrictions apply

The [five and seven reductions](five_seven_pairs.md) and the [factor-nine theorem](factor_nine_pairs.md) give the following necessary cuts. In the first five rows $kt,(d-k)t+1$ add without carries in each of the three indicated bases. The numbers are certified upper bounds for block dimensions, in units of $1/10000$. For each graph choose $L=4096$. All internal words of length $L$ are bounded by the total number $C$ of paths, with $C^{10000}\le b^{Lj}$ for the indicated bound $j$. High zeros can be padded on an accepting integer, and an all-zero input block occurs in the graph.

| $d$ | $k$ | Bases | Dimension numerators |
|---:|---:|---|---|
| 5 | 3 | 2,3,5 | 5979,6788,5033 |
| 7 | 3 | 2,3,7 | 5522,6472,6109 |
| 7 | 5 | 2,5,7 | 6311,7281,5807 |
| 9 | 5 | 2,3,5 | 5737,5453,6982 |
| 9 | 7 | 2,3,7 | 5647,4840,7573 |
| 9 | 8 | 2,3 | 6951,68 |
| 9 | 4 | 2,3 | 5737,5453 |

The remaining multiplier at cofactor seven is $6$. Its monotone seven-digit argument and binary cut give the [two-scale estimate](two_scale_orders.md#even-weights-at-cofactor-seven) whenever this multiplier is required, at either parity. At cofactor nine and even $t$ there are only $k=5,7$. At odd $t$, $k=8$ is covered by the two-base theorem since the sum is $0.7019$, and $k=5$ is impossible. For odd $t$ and $k=7$, use $V=(9t+1)/7$ in binary: $V,2V,4V$ must add without carries and their dimension numerator is 5519. The other two coordinates still use $t$ at bases 3 and 7. This is the affine case of the preceding theorem.

At odd $t$ the remaining $k=4$ requires $t\equiv115\bmod144$. Each prime $r\mid V=(9t+1)/4$ is at least five, and every base-$r$ digit of $V$ is at most $\lfloor(r-1)/4\rfloor$ because $V,V,2V$ add freely. Thus the dimension in the third base is at most $\log_r(1+\lfloor(r-1)/4\rfloor)$. It is at most $0.8200,0.8800,0.8806$ through $r=2000,100000,110000$, respectively. For each such exponent $s$, $(r+3)/r^s$ increases for $r>3s/(1-s)$, so it suffices to check the upper endpoint and the small primes directly. Sum over the finitely many eligible $r$. All exponents in the statement leave positive room for the $\epsilon$ in the counting theorem.

The checker compares enumeration of paths with an independently formed matrix power for every two-part cut, and verifies all exponent inequalities by integer arithmetic.

~~~~sh
python3 check_three_scale_orders.py
~~~~
