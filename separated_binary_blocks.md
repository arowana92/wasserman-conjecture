# A uniform theorem for two separated binary blocks

**Theorem.** For all sufficiently large $N=d2^k+c$ with $d,c\ge1$ and $dc\le N^{1/200}$, every three-row counterexample must contain $E=(1,1,N-2)$. If also
$$\gcd(d,c)\le5,\qquad \gcd(d,c-1)\le5,$$
the full conjecture holds at $N$. The threshold is absolute. The number of binary ones in either block may grow. For instance all offsets $1\le c\le(3\cdot2^k)^{1/200}/3$ above $3\cdot2^k$ are covered for sufficiently large $k$.

Take a binary-free row $A$. Since $Q=2^k=(N-c)/d>c$ for large $N$, its parts have the form $x_i=a_iQ+b_i$. Here the $a_i$ partition the bits of $d$ and the $b_i$ partition those of $c$. At least two $a_i$ are positive: otherwise the second-smallest part is at most $c\le N^{1/200}$, contrary to the [second-part bound](power_intervals.md) for a non-endpoint row. The binary-free row cannot be $E$, whose coefficient is even. Write $a_i=g z_i$ with $g=\gcd(a_1,a_2,a_3)$ and $D=\sum z_i=d/g$. Its positive slopes are distinct by the disjointness of the bits. Thus $D\ge3$ and three positive slopes would give $D\ge6$.

## Removing layers uniformly

For any fixed height $H$, retain offsets $0\le h<H$ with $h\not\equiv c\pmod D$. For each such $h$ choose $i$ with $D\nmid z_i(h-c)$, possible because the $z_i$ have gcd one. At any layer $w=p^e>H$ dividing $N-h$ where $A$ avoids $p$, the residues of its parts sum to $h$. For its residue $j=x_i\bmod w$ we have
$$w\mid z_i(h-c)+D(b_i-j),\qquad 0\le j\le h.$$
This follows from $D x_i=z_i(N-c)+D b_i$. The right side is a nonzero integer of size $O_H(dc)$. Use the product of these $h+1$ nonzero integers in absolute value, multiplied by $\operatorname{lcm}(1,\ldots,H)$, as $K_h$. By the [weighted layer argument](weighted_obstruction.md),
$$N^L\ll_H P\prod_{h\in J}K_h^{H-h}\ll_H P (dc)^{M_H},\qquad
 M_H=\sum_{h=0}^{H-1}(h+1)(H-h)=\binom{H+2}{3}.$$
Here $J$ is any subset of the retained offsets, $L=\sum_{h\in J}(H-h)$, and $P$ is the actual weighted product for the two remaining rows. Indeed the layers in $(N-h)/\gcd(N-h,K_h)$ must all be covered by these rows. If removing two rows, multiply their corresponding $K_h$ at offsets retained for both.

Assume for now that the selection omits $E$. The [Baker–Harman–Pintz bound](https://doi.org/10.1112/plms/83.3.532) supplies a prime $r\le N$ with $G=N-r=O(N^{21/40})$. A row avoiding $r$ has two small parts with sum at most $G$. Since $A$ has at least two parts at least $Q\gg N^{199/200}$, one of the others is such a near row.

For $D\ge5$, at $H=12$ omitting the single class leaves $L\ge78-(12+7+2)=57$. The [ordinary and near cutoff bounds](prime_power_cofactors.md) give $P\ll N^{54.75}$. The additional loss has exponent at most $364/200$, so $54.75+364/200<57$ excludes this case.

For $D=3,4$ use $H=24$, $M_H=2600$. If both remaining rows have minimum at least 24, the ordinary cutoffs $(8,8,8)$ cost $N^{108}$. The near cutoffs $(9,10,5)$ cost $G^{100}N^{15}\ll N^{67.5}$. Here $L\ge300-108=192$ and $108+67.5+2600/200=188.5<192$.

Otherwise call a row with minimum below 24 $B$. Since $D\le4$, one of the $a_i$ vanishes, so $A$ has a part at most $c$. Assign full powers of $N$ to avoiding rows. The product assigned to $A$ and $B$ is at most $23c$, as it divides the product of these two small parts. Hence each part of the third row $C$ is at least $N/(23c)$, and $B$ is the near row. Write $C_i=z'_i N/D'$ with the direction divided by its gcd, so $\gcd_i z'_i=1$, $3\le D'=\sum z'_i\le23c$ and $D'\mid N$. Except at multiples of $D'$, a layer covered by $C$ divides one of the nonzero integers $z'_i h-D' j$ for a suitable fixed $i$ and $0\le j\le h$. Thus its removal costs $O_H(c^{M_H})$ as above.

Retain offsets outside $c\pmod D$ and $0\pmod{D'}$. The classes are compatible since both contain $N$. The retained degree is at least 132. If both moduli are at least 4, each loses at most 84 of the total 300. If one is 3 and the other is at least 6, their losses are at most 108 and 60. Equal moduli 3 omit only one class. For the remaining ordered pairs $(D,D')=(3,4),(3,5),(4,3)$ the minima over compatible classes are 144, 144 and 132, by summing $24-h$ for $0\le h<24$.

On $B$ take $(1,15,8)$, using that its second part tends uniformly to infinity. Its minimum is below 24, so this costs $O(G^{120}N^{36})=O(N^{99})$. The two removals add at most $2\cdot2600/200=26$ to the exponent. Since $125<132$, the case is impossible. This proves the uniform endpoint reduction.

## Excluding the endpoint under the two gcd conditions

At an endpoint counterexample, [each part](endpoint_factors.md) of the other rows exceeds $27^{1/5}(N-1)^{2/5}$. All three $a_i$ of the binary row are therefore positive, so $D\ge6$. As above $D\mid d$. Write $U,S$ for the full powers of $N,N-1$ assigned to $A$, and $T=(N-1)/S$.

Not all determinants $a_i b_j-a_j b_i$ can vanish: otherwise $b_i=z_i c/D$ for all $i$. Integrality and coprimality would imply $D\mid c$, contrary to $\gcd(d,c)\le5$. Since $U$ divides every part, it divides each determinant by eliminating $Q$. A nonzero one gives $U\le dc$.

For each pair of positions $i,j$, at least one of $d b_i-a_i(c-1)$ and $d b_j-a_j(c-1)$ is nonzero. Otherwise those two $b$'s equal $z_i(c-1)/D,z_j(c-1)/D$, and the remaining one equals $1+z_l(c-1)/D$ by summing. This would force $D\mid c-1$, contrary to the other gcd hypothesis. Select one such nonzero integer for each pair. A full power in $S$ divides two actual parts as well as $N-1$, hence divides the selected integer for that pair by $d x_i-a_i(N-1)=d b_i-a_i(c-1)$. Thus $S\le(dc)^3$.

Now $T^2\le (U/3)^3$ from the [endpoint factor assignment](endpoint_reduction.md). This gives $N-1=ST\le(dc)^{9/2}/\sqrt{27}$, impossible. In fact this last endpoint argument by itself applies whenever $dc\le N^{2/9-\epsilon}$ for fixed $\epsilon>0$: the individual part bound still forces all slopes to be positive in that larger range.
