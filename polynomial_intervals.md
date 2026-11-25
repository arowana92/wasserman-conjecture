# Consecutive certified weights on polynomial families

Let $G$ be monic of degree $n\ge1$ with nonnegative integer coefficients, $P=p^a$ a prime power, and $c\ge2$. Write $F=XG+c$ and $R_c=(c^2+\sum_i G_i^2)^{1/2}$.

**Interval criterion.** If $F$ and $F-1$ are irreducible over $\mathbb Q$, $p\mid c(c-1)$, and $P>R_c^6$, the full trinomial conjecture holds at $N=PG(P)+c$. For $p=2$ the base condition holds at every offset.

Indeed a nonendpoint row $A$ avoiding $p$ splits the base-$P$ coefficient blocks of $F$ into nonnegative polynomials $F_1,F_2,F_3$ summing to $F$. Each is nonzero, proper, and coprime over $\mathbb Q$ to both $F$ and $F-1$. A multiple of either of the same degree would equal it, since $F$ is monic. In particular equality with $F-1$ would leave two nonzero nonnegative polynomials summing to 1, an impossibility. Write $d=n+1$. In the [retained-power inequality](polynomial_prime_powers.md#proof-by-specialization), $U$ divides a nonzero resultant of $F$ with one part. Also $S^2$ divides the product of the three nonzero resultants of $F-1$ with the parts. All these resultants have absolute value at most $R_c^{2d}$. Consequently
$$D^3J^2\ge 27(N-1)^2/R_c^{12d}.$$
Since $N-1\ge P^d$ this rules out a counterexample under the stated hypotheses.

## Almost every coefficient vector gives whole intervals

Fix $n\ge1$. Take each of the $n$ lower coefficients of $G$ in $[0,H]$ and let $1\le L\le H$, $H\ge2$. Apart from
$$O_n\!\left(H^{n-1}\big(\log(2H)+(L+1)\log(2L)\big)\right)$$
coefficient vectors, **all** weights $PG(P)+c$ for $0\le c\le L$ and $p\mid c(c-1)$ satisfy the conjecture whenever $P>(H\sqrt{n+2})^{15}$. In particular, if $L=o(H/\log(2H))$, a fraction tending to one of the vectors gives the entire interval of $L+1$ weights at **every** sufficiently large power of two. On the same vectors the first two weights are covered at every prime base.

To see this it suffices to exclude reducible $G$ and reducible $XG+t$ for $t\in\{-1,1,\ldots,L\}$. Use the [paired polynomial theorem](polynomial_prime_powers.md) for offsets 0 and 1. Every polynomial in the count has Mahler measure at most its coefficient Euclidean norm, hence $O_n(H)$. For a monic degree-$n$ polynomial the reducible count is $O_n(H^{n-1}\log(2H))$, including the empty count at degree one. For monic degree $n+1$ with specified nonzero constant $t$, the bound is $O_n(\tau(|t|)H^{n-1})$. These are the Mahler-measure counts in Propositions 10.3 and 10.4 of Grizzard and Gunther. Summing uses $\sum_{t\le L}\tau(t)=\sum_{k\le L}\lfloor L/k\rfloor\ll L\log(2L)$. [Grizzard and Gunther, *Slicing the stars*](https://sites.lafayette.edu/grizzarr/files/2017/10/units.pdf)

For completeness, the integer powers of $H$ come directly from factoring. For monic $A$ of degree $k$ and Mahler measure $M$, its coefficient at $X^{k-i}$ has absolute value at most $\binom{k}{i}M$ by the elementary symmetric functions in the roots. There are $O_k(M^{k-1})$ such polynomials with a specified nonzero constant. In a monic factorization of $XG+t$, the two constants can be chosen in at most $2\tau(|t|)$ ways, and the two Mahler measures multiply. For factor degrees $k,l$ with $k+l=n+1$, summing over dyadic sizes of one factor gives $O_n(H^{\max(k-1,l-1)})$ if $k\ne l$. For $k=l\ge2$ the extra logarithm is absorbed by $H^{n-1}$, and for $k=l=1$ both factors are fixed by their constants.

## Explicit intervals of 1,021 weights

**For every $a\ge1$ and $P=2^a$, every one of the 1,021 weights**
$$P^4+1024P^3+P^2+P+c,\qquad 0\le c\le1020,$$
**satisfies Wasserman's conjecture.**

More generally take $A=2^b\ge8$ in place of 1024 and $0\le c\le A-4$. All these weights are covered whenever $P\ge(2A)^6$. For $2\le c\le A-4$ the term $AX^3$ in $X^4+AX^3+X^2+X+c$ strictly dominates all the other terms together on the complex unit circle, as it does at offset $c-1$. Rouché's theorem gives three roots strictly inside and one outside the circle, counting multiplicity. A nontrivial factorization into monic integer polynomials would have one factor with all roots of modulus less than one. Its nonzero integral constant makes that impossible. The interval criterion applies, with $R_c^2=A^2+3+c^2<(2A)^2$.

The two offsets 0 and 1 work for all $a,b$. The unshifted weight has at most four binary ones and its full two-part is $P$. A positive binary carry-free row must have at least two parts that are pure powers of two. Hence its assigned scale from $N$ is exactly $P$, and it can avoid no prime of $N-1$ since two slots must vanish modulo each such assigned prime. The complementary inequality $27(N-1)^2\le P^3$ fails. For the shifted weight there are at most five binary ones, so a binary carry-free triple has a part that is a power of two. Its assigned scale from the odd weight divides this part, and would equal 1, contradicting the required scale of at least 3.

For $A=1024$ the bulk criterion applies already at $P>1446^6$, since $1024^2+3+1020^2<1446^2$. It remains to check $1\le a\le62$. After the [initial interval through one billion](full_range_1b.md), this leaves 57,064 weights. For any part $z$ of the required binary row, its set bits form a subset of those of $N$. It must satisfy $3\le z\le N-6$ and
$$27(N-1)^2\le \gcd(N,z)^3\gcd(N-1,z(z-1))^2.$$
This follows from the complementary size inequality, since its assigned $U$ divides $z$ and $S$ divides $z(z-1)$. Over all the remaining weights **no single part** meets these conditions. The independent check uses the fact that one part has at most a third of the occupied bits. It enumerates combinations of bit positions and computes the second gcd by its two coprime factors. It tests 8,286,018 sparse parts over exactly the same 57,064 weights.

```sh
python3 check_polynomial_intervals.py
python3 check_polynomial_powers_cross.py
```
