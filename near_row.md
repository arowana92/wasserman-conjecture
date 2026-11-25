# One row must contain almost all the weight in a counterexample

For a positive trinomial row of weight $N$, let $s$ be the sum of its two smallest parts. Put $E=(1,1,N-2)$ and identify permutations of a row. The following restrictions hold for every fixed $0<\delta<5/28$ and all sufficiently large $N$.

1. In any counterexample there is a row $R$ with $s_R\le N^{1/2+\delta}$.
2. In any counterexample omitting $E$, the two other rows each satisfy $s>N^{6/7-\delta}$. In particular, for each prime power $p^j\in[N-N^{6/7-\delta},N]$, this same row $R$ must avoid $p$.
3. In a counterexample omitting $E$, exactly one row has a part at most 3. Its smallest part is 1 or 2. If this is $R$, then $s_R>1.33\sqrt N$. Otherwise its minimum is 1, and its two remaining parts have gcd $\gg_\delta N^{1/8-7\delta/4}$.

Here a row avoids $p$ if $p$ does not divide its coefficient. The interval restriction on $\delta$ keeps $1/2+\delta<6/7-\delta$. Choosing a smaller positive $\delta$ gives all the assertions about larger allowable exponents in the first bound or smaller ones in the second.

## Four nearby multiples

Let $p_i u_i=N-h_i$ for four distinct primes $p_i$ and positive integers $u_i$. Suppose $0\le h_i\le H<p_i$, the $u_i$ are pairwise coprime and $2Hu_i u_j<N$ for all $i\ne j$. Then any selection of three rows has a common divisor among the four primes, or it contains a row with $s\le H$.

Indeed, if each prime is avoided by some row, one row avoids two of them. At either such prime write each part as $x=k p_i+r$, $0\le r<p_i$. The three remainders sum to $h_i$ by Kummer, hence the three quotients sum to $u_i$. Thus
$$\left|\frac xN-\frac{k}{u_i}\right|=\frac{|r-kh_i/u_i|}{N}\le\frac{h_i}{N}.$$
The fractions for the two primes differ by less than $1/(u_i u_j)$ and must agree. Their denominators are coprime, so the common fraction is an integer. Across the three parts these integers are nonnegative and sum to one. One part is therefore at least $p_i u_i=N-h_i$. In fact $s\le\min(h_i,h_j)$.

At $N=305$, the products $2\cdot151,3\cdot101,5\cdot61,7\cdot43$ have offsets $3,2,0,4$. The largest prime below 305 has offset 12. Exact coefficients show that only $(1,1,303)$ and $(1,2,302)$ avoid any pair of the listed primes.

## Producing four coprime cofactors

We use the quantitative bound in §§3 and 6 of Jia and Liu, *On the largest prime factor of integers*, Acta Arithmetica 95 (2000), 17–48. Fix an arbitrarily small admissible $\epsilon>0$ with $K=(8/\epsilon)(1/26+\epsilon/2)$ an integer. Write $Y=x^{1/2+\epsilon}$ and $P=x^{\epsilon/8}$. Their quantity $\Phi$ counts ordered tuples of primes with
$$x<pq_1\cdots q_K\le x+Y,\qquad P<q_j\le2P,\qquad p\asymp_\epsilon x^{25/26-\epsilon/2}.$$
Equation (3) defines the count, equation (5) gives the scale $X=Y(\sum_{P<q\le2P}1/q)^K$, and §6 bounds $\Phi$ below by a fixed positive multiple of $X/\log x$. The prime number theorem on $[P,2P]$ therefore gives $\Phi\gg_\epsilon Y/(\log x)^{K+1}$. [Jia and Liu, full text](https://matwbn.icm.edu.pl/ksiazki/aa/aa95/aa9512.pdf)

Each represented integer has a unique largest prime and at most $K!$ orderings. Once at most three cofactors have been selected, at most $3K$ primes have been used in them. The number of integers in the interval divisible by any of these primes is $O_K(Y/P+1)$, smaller than the count of distinct represented integers for sufficiently large $x$. Choose the next cofactor avoiding them all. This gives four pairwise coprime cofactors $u_i\le (2P)^K\ll_\epsilon x^{1/26+\epsilon/2}$.

Take $x=N-2N^{1/2+\epsilon}$ and $H=2N^{1/2+\epsilon}$. Then $p_i>H$ and $2H u_i u_j<N$ for sufficiently small fixed $\epsilon$. The selected integers are distinct and lie within an interval shorter than each $p_i$, so the large primes are distinct. Taking $\epsilon<\delta$ proves assertion 1.

## Which rows can have small parts?

A 1 and a 3 cannot occur in distinct rows of a counterexample at any weight. To see this, call these rows $A,B$, and the third $C$. Assign each full prime power of $N$ to a row avoiding its prime. All parts of that row are multiples of the assigned power. The product assigned to $A$ is 1 and that assigned to $B$ divides 3. The product assigned to $B$ must be at least 3, since the three positive parts of $C$ are multiples of $N$ divided by that product. Thus $3\parallel N$ and $C=(N/3,N/3,N/3)$.

Row $C$ cannot avoid a prime of $N-1$, since three equal residues cannot sum to 1. Each full power of $N-1$ assigned to $B$ divides 3 or 2 by its part 3, and $3\nmid N-1$. Their product is at most 2. The remaining product, assigned to $A$, divides both companions of its part 1. They sum to $N-1$, so the product assigned to $B$ is exactly 2 and these companions are both $(N-1)/2$. In particular $N$ is odd. Now a prime $p\mid N-2$ exists and is at least 5, since $N$ is an odd multiple of 3 admitting a part 3. To avoid $p$, the three residues modulo $p$ must sum to 2. In $B$ one residue is already 3. In $C$ the sum is a multiple of 3. In $A$ it is odd. This is a common prime, a contradiction.

For sufficiently large $N$ a counterexample contains a part 1 or 2 by the [weighted factorial argument](weighted_obstruction.md). Two rows cannot both have minimum at most 2 since the product of their assigned powers of $N$ would divide 2, leaving fewer than three units for the remaining row. A 2 and a 3 in distinct rows are excluded at large $N$ by the [product bound for two parts in distinct rows](parts_in_different_rows.md). Together with the previous paragraph this proves the claimed uniqueness of the row with a small part.

In a counterexample omitting $E$, every part greater than 3 grows at least as a positive constant times $N^{1/7}$ by [the part-size bound](critical_small_parts.md). In the exceptional row at most one part is bounded, by the lower bound there on $s$. By the [bound for a specified part in another row](endpoint_pair_bound.md), when $s_R\le N^{1/2+\epsilon}$ every part $m\ge2$ outside $R$ satisfies $m^{10}>6^6(N-1)^4/s_R^7$. Apply assertion 1 with $\epsilon<1/14$ to exclude a 2 outside $R$ (choosing $R$ with minimal $s$). For a 1 outside $R$, equation (6) of that bound gives the gcd in assertion 3.

## Separating the other two sums

Let $R$ minimize $s$ and suppose a second row $A$ has $s_A\le N^\alpha$, with $1/2<\alpha<6/7$. By the first assertion, $s_R\le N^\gamma$ for any fixed $\gamma>1/2$ at sufficiently large weights. Apply the [weighted divisibility lemma](weighted_obstruction.md) with a large fixed height $H$. Its left side has degree $W(H)=H(H+1)/2$ in $N$. A part bounded by $N^z$ with cutoff $t$ costs at most $N^{zW(t)}$.

For a row with $s\le N^z$ and all parts growing, use cutoffs in proportions $(1,1,z)/(2+z)$. Its degree divided by $H^2$ is at most $f_2(z)+O(1/H)$, where $f_2(z)=z/(2(2+z))$. For a row with one bounded part give it cutoff zero and the others proportions $(1,z)/(1+z)$. The corresponding bound is $f_1(z)=z/(2(1+z))$. Integer cutoffs summing to $H$ within one of the prescribed values suffice. In the third row put $z=1$ and use $f_2(1)=1/6$ or $f_1(1)=1/4$ as applicable. All cutoffs fit for large $N$.

According to the location of the one bounded part, the total leading degree divided by $H^2$ is at most one of
$$f_2(\gamma)+f_2(\alpha)+1/4,\quad
 f_2(\gamma)+f_1(\alpha)+1/6,\quad
 f_1(\gamma)+f_2(\alpha)+1/6.$$
At $\gamma=1/2$ these are strictly less than $1/2$ for $\alpha<6/7,\ 7/8,\ 1$, respectively. Fix $\gamma$ close enough to $1/2$ to preserve strictness, then $H$ large enough to absorb rounding, then $N$ large enough to dominate the constant in the lemma. This contradicts the left side. Taking $\alpha=6/7-\delta$ proves the second sum bound. If $p^j\in[N-N^\alpha,N]$, then $N<2p^j$ for large $N$. In a row avoiding $p$ the quotients of the parts by $p^j$ sum to 1, so its largest part is at least $p^j$. Only $R$ can do this.

## If the small part is in the minimizing row

Suppose $s_R\le c\sqrt N$ for fixed $c$ and the minimum of $R$ is 1 or 2. Use $t=\lceil(\log N)^2\rceil$, $H=6t$ and the [sharper weighted inequality and loss estimate](critical_small_parts.md). Give $R$ cutoffs $(1,4t,2t-1)$ and the other rows $(2t,2t,2t)$. The growing parts are $\gg N^{1/7}$, so these fit. The product on the right, before the loss factor $D(H)$, is at most
$$2s_R^{W(4t)}N^{W(2t-1)}(N/3)^{6W(2t)}.$$
Substitute $s_R\le c\sqrt N$. The degree in $N$ exceeds $W(H)$ by $O(t)$, contributing $o(1)$ after taking logarithms and dividing by $H^2$. The bound $\limsup\log D(H)/H^2\le\pi/(6\sqrt3)$ then implies
$$0\le\frac{\pi}{6\sqrt3}+\frac29\log c-\frac13\log3.$$
This excludes every $c<\exp((9/2)((\log3)/3-\pi/(6\sqrt3)))$, including $c=1.33$.

## Reproduction

Run `python3 check_near_row.py` to verify the finite four-prime example by direct binomial products, the local obstruction for a 1 and a 3 through weight 100, and the rational degree comparisons. The asymptotic results use the published count and the preceding proofs.
