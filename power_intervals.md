# Intervals above powers of two

**Theorem.** Fix $0<\alpha<2/5$. For every sufficiently large $j$, the three-row conjecture holds at **every** integer
$$2^j\le N\le 2^j+2^{\alpha j}.$$
More generally, for sufficiently large $N$, the same conclusion holds whenever a prime power $Q=p^j\le N$ satisfies $N-Q\le N^\alpha$ and $p\mid N(N-1)$.

An explicit, smaller interval is valid without an asymptotic threshold. For $N\ge12$ it suffices that $p\mid N(N-1)$ and
$$0\le N-Q < 120^{-7/5}N^{1/5}.$$

## A size bound for the second part of every row

The following restriction does not require $N$ to be close to a prime power. Write each row in increasing order. In a counterexample of sufficiently large weight, **every row other than $E=(1,1,N-2)$ has second part greater than $N^\alpha$** for each fixed $0<\alpha<2/5$.

If $E$ is present this follows from the stronger bound $a^5>27(N-1)^2$ on each part $a$ of the other two rows. We prove the claim when $E$ is absent. The two smaller parts of each row then sum to at least 11 by Bergman’s Proposition 10. [Endpoint bound](endpoint_factors.md) · [Bergman](https://arxiv.org/abs/0806.0607)

We use the weighted inequality with height $H$. Its left side has degree $W(H)=H(H+1)/2$ in $N$. A cutoff $t$ on a part $x$ costs at most $x^{W(t)}$, the cutoffs of each row sum to $H$, and a constant depending on $H$ multiplies the right side. All heights in this proof are fixed independently of $N$. [Weighted divisibility lemma](weighted_obstruction.md)

First, in **any** selected row $(a,b,c)$ we have $b\ge120^{-7/5}N^{1/5}$ for $N\ge12$. Apply the weighted inequality with $H=6$. Each of the other two rows costs at most $N^9$. Indeed, take $(2,2,2)$ when its minimum is at least 2, or $(1,2,3)$ when the minimum is 1. On the selected row take $(1,5,0)$ if $a=1$, or $(2,4,0)$ otherwise. These cost at most $b^{15}$. Thus $N^{21}\le120^{21}N^{18}b^{15}$. In particular every second part eventually exceeds any fixed height.

Suppose now that a selected row $D=(a,b,c)$ has $b\le N^\alpha$. Fix a height $H\ge18$, to be chosen sufficiently large in terms of $\alpha$. The minima of the other two rows cannot both be at most $H$ for large $N$. If they were $x,y\le H$, let $U$ be the product of the full powers of primes dividing $N$ avoided by either of their rows. Each such power divides $x$ or $y$, so $U\mid\operatorname{lcm}(x,y)$. All parts of $D$ must be multiples of $N/U\ge N/H^2$, contrary to $b\le N^\alpha$. This observation also proves that two different rows cannot both have a part at most 2, since then $U\le2$ and the last row cannot have three positive parts.

For the exponent estimates below, an even division of $H$ cutoffs among three parts costs $N^{H^2/6+O(H)}$. Using just two parts costs $N^{H^2/4+O(H)}$. Integer rounding changes each cutoff by $O(1)$, accounting for the $O(H)$ terms. The second and third parts already exceed $H$ for large $N$.

- **If $a\ge H$**, give the first two parts of $D$ cutoffs differing by $O(1)$ from $H/(2+\alpha)$, and its third the remaining cutoff. Their combined cost has exponent $H^2\alpha/(2(2+\alpha))+O(H)$ in $N$. One other row has all parts $>H$ and can split its cutoffs over three parts. The last row can use its two larger parts. The sum of the quadratic coefficients is
  $$\frac{\alpha}{2(2+\alpha)}+\frac16+\frac14<\frac12$$
  exactly for $\alpha<2/5$. Choosing $H$ large gives an exponent strictly below $W(H)$, a contradiction for large $N$.

- **If $3\le a<H$**, apply the height 6 bound again. Give $D$ the cutoffs $(3,2,1)$, costing $O_H(b^3N)$. The other rows cost at most $N^{18}$ together. This gives $N^2\ll_H b^3$, inconsistent with $b\le N^\alpha$.

- **If $a\le2$ and both other minima exceed $H$**, use their three-part splits. On $D$ give its second part a cutoff within 1 of $H/(1+\alpha)$, and its third the remainder. The sum of quadratic coefficients is
  $$\frac{\alpha}{2(1+\alpha)}+\frac13<\frac12$$
  for $\alpha<1/2$. Again the exponent is too small.

- **If $a\le2$ and another row has minimum at most $H$**, that minimum is at least 3 and the last row has minimum $>H$. Use height 15. Put $(1,10,4)$ on $D$, $(3,6,6)$ on the other row of bounded minimum, and $(5,5,5)$ on the last. These fit because all second parts eventually exceed $H$. The costs are $O(b^{55}N^{10})$, $O_H(N^{42})$, and $N^{45}$. Thus $N^{120}\ll_H N^{97+55\alpha}$. The exponent on the right is at most 119 even at $\alpha=2/5$.

These cases prove the size bound.

## Applying a prime-power digit

Suppose a counterexample exists and select a row avoiding $p$. When $Q=p^j$ and $Q\le N<2Q$, a carry-free row must allocate the single digit at position $j$ to one part. The other two positive parts have sum at most $N-Q$. The row cannot be $E$ because $p\mid N(N-1)$, its coefficient. Its second part is therefore greater than $N^\alpha$ for sufficiently large $N$, contradicting $N-Q\le N^\alpha$. Taking $p=2$ automatically gives the required divisibility, proving the interval result.

For the explicit bound, the same digit argument applies since the stated distance is less than $N/2$ for $N\ge12$. Without $E$ the height 6 calculation already gives the contradiction. With $E$, every other part exceeds $27^{1/5}(N-1)^{2/5}$, which is greater than $120^{-7/5}N^{1/5}$ for $N\ge12$.
