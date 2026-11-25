# Weights with bounded binary digit sum

Write $s_2(N)$ for the number of ones in the binary expansion of $N$, and $E=(1,1,N-2)$.

**Theorem.** For every fixed $s$, all sufficiently large weights with $s_2(N)\le s$ satisfy the following assertion. Any selection of three rows whose coefficients have gcd one must contain $E$, up to permutation.

**Corollary.** Wasserman’s full three-row conjecture holds for all sufficiently large **odd** $N$ with $s_2(N)\le5$. The threshold is uniform in the bit positions. The proof does not give its value. An independent finite calculation below covers all such weights below $2^{64}$.

For the corollary, a counterexample must have a row $A$ avoiding 2. Its parts partition the binary ones into three nonempty sets, so with at most five ones some part is a power of two. Since $N$ is odd, $A$ is primitive (its parts have gcd one). $E$ is another primitive row, and its coefficient $N(N-1)$ is even. Two primitive rows in a counterexample are impossible. Indeed, assign each full power $p^{v_p(N)}$ to a row avoiding $p$. It divides every part of that row by Kummer, so the third row would have every part divisible by $N$.

## The estimates for two remaining rows

Suppose to the contrary that an infinite sequence of distinct counterexample weights omits $E$ and has $s_2(N)\le s$. Pass to a subsequence tending to infinity. Designate a binary carry-free row $A$ and denote the two others by $B,C$. Disjointness of the binary digits of the three parts of $A$ will be important. The second part of each row, sorted increasingly, tends to infinity by the [second-part bound](power_intervals.md). Every row's two smallest parts sum to at least 11 by Bergman’s Proposition 10. A row avoiding the largest prime $r\le N$ has two parts of sum at most $g=N-r$, with $g=O(N^{21/40})$. [Bergman](https://arxiv.org/abs/0806.0607) · [Baker–Harman–Pintz](https://doi.org/10.1112/plms/83.3.532)

Use the [weighted divisibility lemma](weighted_obstruction.md) at height $H$. Set $W(t)=t(t+1)/2$. Cutoffs $t_i$ on a row sum to $H$, fit its three parts $x_i$, and cost at most $\prod x_i^{W(t_i)}$. The actual product consists of the weighted falling factorials in the lemma. For offsets $J\subseteq\{0,\ldots,H-1\}$ write $L=\sum_{h\in J}(H-h)$. If all layers $w=p^e>H$ at these offsets where $A$ avoids $p$ divide some $K$, then the actual product $P$ for $B,C$ satisfies
$$P\ge (N/(2K'))^L,\qquad K'=\operatorname{lcm}(K,1,\ldots,H),\quad N\ge2H.$$
This follows by removing $\gcd(N-h,K')$ from each $N-h$. We will also use this argument with two rows removed when both of them have bounded layers. Constants and heights below can depend on the chosen subsequence.

Call a row **on a fixed line** if its actual triples on the subsequence are integer points of one fixed affine line. We first treat this case for $A$.

## When the binary row is on a fixed line

An affine line through infinitely many distinct integral points has a primitive integer direction. Write $A_i=a_i t+b_i$ with fixed integers, $\gcd(a_1,a_2,a_3)=1$ and $t$ integer tending to $+\infty$ after orienting the line. All $a_i\ge0$ and at least two are positive by the second-part bound. Two positive slopes cannot agree. Indeed, if the greater of two binary-disjoint parts has highest bit $2^j$, the other has at most $s$ bits in smaller positions, and the difference is at least $2^{j-s}$. Two unbounded parts therefore cannot differ by a constant. It follows that $D=\sum a_i\ge3$, and that three positive slopes would give $D\ge6$. Put $b_0=\sum b_i$.

At offsets $h\not\equiv b_0\pmod D$ this row can cover only bounded layers. At such a layer its residues $u_i$ sum to $h$. Choose integers $l_i$ with $\sum l_i a_i=1$. Then
$$w\mid h-b_0-D\sum_i l_i(u_i-b_i).$$
The integer on the right is nonzero and bounded for fixed $H$. The prime $r$ is avoided by $B$ or $C$, since $A$ has at least two parts comparable to $N$. Call a row avoiding $r$ the near row.

If $D\ge5$, use $H=12$. Omitting one residue class costs at most $12+7+2$, leaving $L\ge57$. The combined cost of the near and ordinary rows is $O(N^{54.75})$ by the [height-12 bounds](prime_power_cofactors.md). This contradicts the lower bound.

If $D=3,4$, there is a fixed part in $A$. Pass to a subsequence on which each of the minima of $B,C$ is fixed or tends to infinity. If both tend to infinity, give the ordinary row cutoffs as equal as possible among three parts, and give the near row cutoffs as equal as possible among its two small parts. The combined exponent in $N$ is at most
$$H^2(1/6+(21/40)/4)+O(H).$$
Retaining all but one class modulo $D$ gives $L\ge H^2/3-O(H)$. The strict inequality between the quadratic coefficients gives a contradiction on choosing $H$ large and fixed.

Otherwise suppose $B$ has a fixed minimum. The prime powers of $N$ assigned to $A,B$ have bounded product $k$, since this product divides the product of their fixed parts. Pass to a subsequence with $k$ fixed. The parts of $C$ are $N/k$ times a fixed positive triple of sum $k$, after another subsequence. Reducing this direction by its gcd gives a primitive direction of sum $D'\ge3$; thus $C$ is a line for which the exceptional class is $0\pmod{D'}$. Both $A,C$ have at least two parts comparable to $N$, so $B$ is near. The omitted classes $b_0\pmod D$ and $0\pmod{D'}$ are compatible since $N$ belongs to both. Their union has density at most $1/2$. For two moduli at least 4 this follows from the sum of their densities. For one modulus 3, the other cases to check are 3, 4 and 5, where inclusion–exclusion gives $1/3,1/2,7/15$; moduli at least 6 give the bound directly. Periodicity now gives $L\ge H^2/4-O_{D,D'}(H)$. On $B$ put cutoffs $0,2H/3,H/3$ in increasing order of parts, with integer rounding preserving the sum. Its exponent is at most
$$H^2\big((21/40)\,2/9+1/18\big)+O(H)<H^2/4-O_{D,D'}(H)$$
for sufficiently large $H$. This excludes the last line case.

We can consequently assume at least one of the minima of $B,C$ tends to infinity. Otherwise pass to two fixed minima. The prime powers of $N$ assigned to $B,C$ have bounded product, forcing $A$ to lie on one of finitely many lines through the origin by the same normalization. Along a subsequence one line is fixed. With at least one minimum tending to infinity we can give that row a three-part split and split the other row’s cutoffs between its two larger parts. For every fixed $H$ this gives
$$P\le N^{5H^2/12+O(H)} \tag{1}$$
eventually, with an absolute implied constant in the exponent.

## Restricting the exponents to an algebraic torus

Pass to a fixed number of bits and a fixed partition of their positions among the parts of $A$. Write the increasing exponent vector as $e$. Among rational affine subspaces of the exponent space containing infinitely many of these vectors, choose one of smallest dimension $d$ and retain the vectors in it. We have $d\ge1$. Each proper rational affine subspace inside it now contains only finitely many sequence terms. Its integer points have the form
$$e=e_0+M v,\qquad v\in\mathbb Z^d,$$
where the columns of the fixed integer matrix $M$ are a lattice basis. In particular $\|v\|_\infty=O(1+\max e_i)=O(\log N)$. Put $x_j=2^{v_j}$. The actual parts of $A$ are evaluations of three fixed rational Laurent polynomials $F_i(x)$, obtained by adding the monomials for their bits. Let $Y$ be the Zariski closure over $\mathbb C$ of the image of $F=(F_1,F_2,F_3)$ on the torus. This is an irreducible variety. Its dimension is positive, as the weights are unbounded.

A nonzero Laurent polynomial has only finitely many zeros on this sequence. At a zero, its two largest monomial exponents evaluated at $2^v$ must differ by a bounded amount, since otherwise the leading term outweighs the sum of the others. The bound depends only on its finitely many nonzero coefficients. The difference is a nonzero integer linear form in $v$, and each possible value gives a proper affine subspace.

For any fixed integral $u=(u_1,u_2,u_3)$ for which the three $F_i-u_i$ have no common nonunit factor over the complex Laurent polynomial ring, we have
$$\log\gcd_{\rm odd}(A_1-u_1,A_2-u_2,A_3-u_3)=o(\log N). \tag{2}$$
Here the gcd takes absolute values and ignores the factor 2. To justify (2), take one nonzero polynomial and a rational integer linear combination of the three avoiding its finitely many irreducible factors. Such a choice exists since each factor excludes a proper linear subspace of choices. If a polynomial is a Laurent unit, its value has bounded odd part and the assertion is immediate. Otherwise these two polynomials are coprime. Clear powers of 2 in their coefficients, and clear negative exponents while stripping coordinate monomials from each polynomial individually. This gives coprime ordinary polynomials and does not change odd valuations at $x$. The odd gcd of the three values divides the resulting pair's odd gcd, up to fixed constants.

Apply the outside-$S$ gcd theorem with $S=\{2,\infty\}$ to the pair. For any $\epsilon>0$ its bound for the log gcd is $\epsilon\max h(x_i)$ off a finite union of proper translated algebraic subgroups. Each such translate meeting our points imposes a nonzero character equation $x^m=c$, hence $m\cdot v=$ constant at these points. This is a proper rational affine subspace, containing only finitely many of our terms by minimality. The coordinate heights are $|v_j|\log2=O(\log N)$, so arbitrary $\epsilon$ gives (2). For $d=1$ the same conclusion follows directly from the polynomial Bézout identity. [Levin–Wang, Theorem 1.13(a), stating the Corvaja–Zannier and Levin gcd theorem](https://doi.org/10.1515/crelle-2019-0033)

## Counting the exceptional residues

If $\dim Y\ge2$, two of the coordinate polynomials are algebraically independent. In characteristic zero some two-by-two minor of their Jacobian is a nonzero Laurent polynomial $J_0$. Every irreducible factor common to all $F_i-u_i$ divides $J_0$: modulo such a factor the differentials of the shifted polynomials are multiples of the differential of the factor, so the determinant vanishes. There are only finitely many irreducible factors of $J_0$. Each fixes at most one triple $u$, because divisibility at two different values of one coordinate would imply divisibility of their nonzero constant difference. Thus only finitely many integral residue triples fail the coprimality hypothesis of (2), independently of $H$.

If $\dim Y=1$ and $Y$ is a line, all actual rows lie on that line, already excluded. For a nonlinear curve $Y$, all triples failing the hypothesis of (2) belong to $Y$: a nonunit common factor has a zero in the complex torus. Only $o(H)$ integer points of $Y$ lie in $[0,H]^3$. Indeed, choose a nonconstant coordinate. Some projection onto it and one other coordinate has nonlinear irreducible plane curve as its closure, since otherwise every coordinate would be affine linear in the chosen one. Bombieri–Pila give $O_{Y,\epsilon}(H^{1/2+\epsilon})$ integral points on that projection. Each lifts to at most $\deg Y$ points by intersection with a plane fixing the chosen nonconstant coordinate. [Bombieri–Pila, Theorem 5](https://people.maths.ox.ac.uk/pila/Ovals.pdf)

Fix $H$. Omit all offsets $0\le h<H$ that are a sum of three nonnegative residues in one of these exceptional triples. Omit also the one offset, if present, for which $N-h$ is divisible by a power of 2 exceeding $H$. There is at most one since such an offset must be congruent to $N$ modulo the smallest power of two greater than $H$. The retained degree satisfies $L\ge W(H)-o(H^2)$ as $H$ grows, uniformly on the subsequence. At any retained odd layer covered by $A$, the actual residues sum to $h$, and the layer divides the gcd in (2) for that residue triple. For fixed $H$ multiply those odd gcds over all nonexceptional triples of sum less than $H$. By (2) this provides $K=N^{o(1)}$ as $N$ grows, besides the bounded layers. Choose $H$ large enough that the retained degree exceeds the exponent in (1) by a fixed positive margin. The lower bound $P\ge N^{L-o(1)}$ now contradicts (1). This completes the proof.

## Finite verification for odd weights

The full conjecture holds at all 333,375 odd four-bit weights below $2^{128}$ and all 595,665 odd five-bit weights below $2^{64}$. For fewer bits use the [three-bit verification](even_three_binary_digits.md) and the universal divisor 2 when there are at most two ones.

The calculation first uses the existing [one-million theorem](full_range_1m.md). At the other weights the endpoint is impossible by primitivity. With $H=6$, the weighted lemma gives the lower bound
$$F=\prod_{h=0}^5\left((N-h)/\gcd(N-h,60)\right)^{6-h}$$
for the product over all three rows. The row avoiding the proved prime $r$ has its small-part sum bounded by $g$. A gap of at most 10 excludes the selection; otherwise use the cutoffs $(1,5,0)$ or $(2,4,0)$ according to whether its minimum is 1 or larger. These cost at most $g^{15}$. The two other rows cost at most $64 N^{18}/3^{18}$ together by the one-million proof. This excludes 331,354 and 583,644 weights respectively, after accounting for 969 and 3,876 previously covered weights.

For the remaining 1,052 and 8,145 weights the verifier enumerates the 6 or 25 partitions of the bits into three nonempty parts. Each possible binary row is excluded by the [individual-row height-six and prime-gap checks](even_three_binary_digits.md). The verifier uses python-flint for proved primality and integer arithmetic for all comparisons. A separate run using individual-row tests throughout gave the same conclusion for four bits below $2^{64}$ and five bits below $2^{40}$.

```sh
python3 check_sparse_odd.py
```
