# Almost every pair of binary blocks works at every gap

For positive integers $d,c$ let $N_a=d2^a+c$, $a\ge1$. Among $1\le d,c\le X$, at most **$O(X^{1.78})$ pairs** can fail Wasserman's three-row conjecture at even one exponent with $2^a>c$. If all positive exponents are included, the bound is **$O(X^{1.869})$**.

The quantifier over gaps is uniform. For example, fix an integer $L$ and any residue $b$. Among pairs with $d\equiv0$ and $c\equiv b\bmod L$, a proportion tending to one gives only good weights at every positive exponent, all in the class $b$. This proportion tends to one uniformly for $L\le X^{0.065}$. For separated blocks alone one can take $L\le X^{0.109}$. Indeed each box contains $\gg(X/L)^2$ pairs in these two classes. By the Chinese remainder theorem one can prescribe any finite disjoint lists of prime divisors of the weights and of their predecessors.

## Small and intermediate gaps

Enlarge the box to $1\le d,c<X=2^n$, at a cost of at most a factor two in its side. Write $Q=2^a$. For $Q<X$ with $c<Q$ each weight has a unique representation at that $Q$, and the weight is below $2XQ$. The [global exceptional-set estimate](two_predecessor_congruences.md) gives $O((XQ)^{0.8689})$ possibilities. Summing over the powers of two gives $O(X^{1.7378})$. With unrestricted $c$, a weight has at most $1+X/Q$ representations. Summation then gives $O(X^{1.8689})$.

Consider $X\le Q<X^4$. Set
$$
 \tau=79/114,\qquad \gamma=85/109,\qquad \kappa=1+\gamma=194/109.
$$
We show that $X^{\kappa+o(1)}$ pairs suffice in this range. In any counterexample take the [endpoint](all_weights_endpoint_reduction.md) and name the other rows $A,C$, with $A$ avoiding 2. Assign the full prime powers of $N$ and $N-1$ to rows avoiding their prime bases, writing $UV=N$, $ST=N-1$ with $U,S$ assigned to $A$. All parts of $A$ are divisible by $U\ge3$, and at each power in $S$ their residues are $0,0,1$ in some order. Thus any part $y$ of $A$ gives $U\mid y$ and $S\mid y(y-1)$. These specify at most $X^{o(1)}$ pairs $(U,S)$ for fixed $y$, by the divisor bound and $3\le y<2X^5$. Each fixes $N$ in one class modulo $K=US$.

For $K\le X^\kappa$, the [reconstruction bound on the complementary row](two_predecessor_congruences.md#recovering-the-row-from-a-product) gives $O(X^{\kappa+o(1)})$ possible weights. Apply it on dyadic intervals up to $2X^5$, counting the weights below $X^\kappa$ directly. Each weight fixes its pair at this $Q$.

There are $2n$ possible binary positions. The two blocks are disjoint. At most $2^{2nH(\tau)}$ pairs have more than $2n\tau$ ones, where $H(t)=-t\log_2t-(1-t)\log_2(1-t)$. For the remaining pairs one part $y$ receives at most $2n\tau/3$ ones, so at most $X^{2H(\tau/3)}$ distinct values of $y$ occur. Both estimates follow by bounding binomial tails with the binomial theorem. The exact inequalities $2H(\tau)<\kappa$ and $H(\tau/3)<\gamma$ are checked by integer powers below.

For fixed $U,S$, the pairs lie in an affine translate of the lattice
$$
 \Lambda_K=\{(i,j)\in\mathbb Z^2:Qi+j\equiv0\bmod K\}.
$$
It has determinant $K$. If $\lambda$ is the sup norm of a shortest nonzero vector, every translate has $O(1+X/\lambda+X^2/K)$ points in the box. For $\lambda\ge X$ there is at most one. Otherwise a shortest vector $v$ is primitive in the lattice. Parallel lattice lines in its direction are spaced $K/\|v\|_2$ apart, with spacing $\|v\|_2$ along each line. Since $\lambda\le\|v\|_2\le\sqrt2\lambda$, the count is $O((1+X\lambda/K)(1+X/\lambda))$, giving the bound.

For $K>X^\kappa$ and $\lambda>X^\gamma$, summing over the choices supplied by $y$ gives
$$
 X^{2H(\tau/3)+o(1)}(1+X^{1-\gamma}+X^{2-\kappa})\ll X^{\kappa+o(1)}.
$$
For $\lambda\le X^\gamma$ put $P\le\lambda<2P$ in dyadic intervals. There are $O(P^2)$ possible shortest vectors $(i,j)$. The integer $Qi+j$ is nonzero, since $\max(|i|,|j|)<X\le Q$ and the vector is nonzero. Its absolute value is polynomially bounded. Consequently its divisors yield only $X^{o(1)}$ possibilities for $K,U,S$ at a fixed vector. The lattice count now gives per interval
$$
 X^{o(1)}(P^2+XP+P^2 X^{2-\kappa})\ll X^{\kappa+o(1)}.
$$
There are $O(\log X)$ intermediate exponents, absorbed by the margin $\kappa<1.78$.

## Arbitrarily long gaps

When $Q\ge X^4$, the inequality $d^3c^5\le X^8\le Q^2$ implies $(dc)^5<27(N-1)^2$. By the [determinant theorem](projected_rows.md#two-blocks-of-digits-in-any-prime-base), the upper three parts of the binary row must be $gz$, where $z$ is a positive primitive integer triple, and its lower parts must be $mz+e$. Here $g\ge1$, $m\ge0$, and $e$ is zero or one coordinate unit vector. Both triples must add without binary carries. Writing $s=\sum z_i$, we have $s\mid d$ and $s\mid c$ or $s\mid c-1$. None of these necessary conditions depends on $a$.

The pairs with $c=1$ contribute $O(X)$. For $c>1$ and $s>Y=X^{2/9}$ the divisibilities bound the total by $O(\sum_{Y<s\le X}(X/s)^2)=O(X^{16/9})$. For $s\le Y$ there are $O(Y^3)$ directions and four choices of $e$. For any such choice there are $O(X^\alpha)$ possibilities for $g$ and independently $O(X^\alpha)$ for $m$, with $\alpha=\log_2\lambda_3$ and $\lambda_3^3=\lambda_3^2+1$, $\lambda_3>1$. Indeed both inputs are at most $X/s$. Use the [affine triple count](binary_affine_counts.md) at $L=\lceil\log_2(1+X/s)\rceil$ with delay spread at most $\log_2 s$. This contributes $O(X^{2\alpha+2/3})=O(X^{16/9})$ since $\alpha<5/9$. Together the three gap ranges prove the theorem.

## Finite families with every gap verified

A direct count marks each pair for which an unshifted or shifted direction exists. In the shifted case the two unchanged coordinates of $z$ must have gcd at least three. The column headings record which necessary patterns exist.

| Upper limit on both blocks | Neither | Unshifted only | Shifted only | Both |
|---:|---:|---:|---:|---:|
| 31 | 913 | 28 | 20 | 0 |
| 63 | 3,791 | 98 | 80 | 0 |
| 255 | 63,457 | 827 | 741 | 0 |
| 4,095 | 16,699,343 | 37,220 | 32,461 | 1 |
| 16,383 | 267,985,169 | 223,234 | 194,272 | 14 |

The first count column certifies all separated weights satisfying $(dc)^5<27(N-1)^2$. For $d,c\le255$ the [shifted-row exclusion](shifted_binary_prefixes.md) also applies to every shifted pattern. If $d$ is even, move its power of two into $Q$. This decreases $dc$ and leaves the primitive upper directions unchanged. Thus **64,198 of the 65,025 pairs through 255 satisfy the conjecture at every positive exponent**. Their exact description is the absence of an unshifted pattern. The remaining 827 pairs pass every gap outside the long-gap determinant regime as well.

For this last assertion, weights up to $10^9$ use the [contiguous verification](full_range_1b.md). At each larger weight before the determinant criterion applies, enumerate the [full prime-power assignments](endpoint_reduction.md) and test positive rows by CRT and digit addition. In 146,015 parameter triples, up to weight 136,902,082,815, none of the 1,313,156 size-feasible assignments has a row on its larger-modulus side. At and beyond the first separated gap satisfying the determinant inequality, it holds for every larger gap.

Run the following from the repository. The finite-gap script uses python-flint. The independent Python pattern enumeration compares upper and lower parts by cross multiplication and reproduces the table through 63.

~~~sh
python3 check_binary_gap_bounds.py
c++ -std=c++17 -O3 check_binary_gap_templates.cpp -o /tmp/wasserman-gap-templates
for n in 31 63 255 4095 16383; do /tmp/wasserman-gap-templates "$n"; done
python3 check_binary_gap_finite.py
~~~
