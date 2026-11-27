# Exceptional weights when an odd prime divides the weight or its predecessor

For $p=3$, the number of possible failures of Wasserman's full three-row conjecture at $N\le X$ with $N\equiv0,1\pmod3$ is **$O(X^{0.8178})$**. The count includes every choice of three positive rows. More generally use the following exponents whenever the indicated odd prime $p$ divides $N(N-1)$.

| $p$ | Upper exponent |
|---|---:|
| 3 | .8178 |
| 5 | .8293 |
| 7 | .8369 |
| 11 | .8462 |
| 13 | .8493 |
| 17,19,23 | .8581 |
| 29 through 53, prime | .8665 |
| 59 through 73, prime | .8685 |
| 79 | .86885 |

The union of these conditions is periodic of natural density $1-\prod_{3\le p\le79,\ p\ {\rm prime}}(1-2/p)>0.959$. Consequently $O(X^{0.86885})$ covers their union, while the [unrestricted bound](two_predecessor_congruences.md) is $O(X^{0.8689})$. At $p=3$ the subcase where one nonendpoint row avoids both 2 and 3 has the stronger bound $O(X^{0.7611})$.

## A weighted count in two bases

Fix an odd prime $p$ and $0<u,v<1$. Write $D_b(w)=\log_b(1+w+\cdots+w^{b-1})$ and
$$
 A=\log_2(1/u),\quad B=(p-1)\log_p(1/v),\quad
 \Delta=\max(0,D_2(u)+D_p(v)-1).
$$
Let $M_b=b^{L_b}$ with both $M_2,M_p$ comparable to $X$. Uniformly in the integer shift $a$,
$$
 \sum_{\substack{0\le z<M_2\\0\le z+a<M_p}} u^{s_2(z)}v^{s_p(z+a)}
             \ \ll_\epsilon X^{\Delta+\epsilon}.                 \tag{1}
$$
Here $s_b$ is the base-$b$ digit sum.

To prove (1), take independent Bernoulli digits in each base with probabilities $w^d/(1+\cdots+w^{b-1})$, and let $\mu$ be the product of the two resulting measures on $[0,1]$. A pair of specified length-$L_b$ words gives a product rectangle of mass exactly the corresponding summand divided by $M_2^{D_2(u)}M_p^{D_p(v)}$. Rectangles for different summands have disjoint interiors; their boundaries have measure zero. All lie in one strip of width $O(X^{-1})$ about the line $y=(M_2/M_p)x+a/M_p$. The slope is bounded above and away from zero.

For any $q>1$, the $L^q$ dimension of the measure in base $b$ is $-\log_b(\sum_d\Pr(d)^q)/(q-1)$ by direct cylinder counts, and it is at least $D_b(w)=-\log_b\max_d\Pr(d)$. For projection to the normal of the strip the dimension is the smaller of 1 and the sum of the two $L^q$ dimensions. This is the convolution theorem for homogeneous self-similar measures at incommensurable logarithmic scales. The two grids have exponential separation at every level; distinct translations at level $n$ are separated by at least $b^{-n}$. Reflection, scalar dilation and translation preserve the hypotheses. The contraction logarithms have irrational ratio since 2 and $p$ are distinct primes. [Shmerkin, Theorem 7.2](https://arxiv.org/pdf/1609.07802)

The finite-scale lower dimension estimate is uniform for these normals. Specifically, the product model has $\mu$ at its identity parameter, and Proposition 3.2 of Corso and Shmerkin gives the lower bound for the finite-scale $L^q$ spectrum locally uniformly in the projection and uniformly in the model parameter. Use the dimension formula above at each noncoordinate direction, or directly their Theorem 5.4, and take a finite cover of the compact set of normals. Thus each interval of length $h$ in the projected line has mass at most $O_\epsilon(h^{\min(1,D_2(u)+D_p(v))-\epsilon})$: bound its mass by a constant times the $q$-norm of the dyadic masses at this scale, and take $q$ so large that the factor $(q-1)/q$ costs less than $\epsilon/2$. Apply at $h\asymp X^{-1}$ and multiply back the two rectangle normalizations to get (1). [Corso and Shmerkin, Lemma 5.1, Proposition 3.2 and Theorem 5.4](https://arxiv.org/pdf/2409.04608)

In particular, for any random integer $Z$ supported on the summation range, entropy against the probability proportional to the summands gives
$$
 \frac{\mathcal H(Z)}{\log X}
 \le \Delta+\frac{\mathbb E[s_2(Z)]\log(1/u)+\mathbb E[s_p(Z+a)]\log(1/v)}
 {\log X}+o(1).                                                   \tag{2}
$$
The little-oh in this upper bound means a limsup bound with arbitrarily small fixed loss. It is uniform in the shift and the law of $Z$.

## Assigning the two primes

On $[X,2X)$ fix a canonical failing pair $C,F$ with the endpoint, using the [global endpoint reduction](all_weights_endpoint_reduction.md). Each full power dividing $N$ or $N-1$ is assigned to a row whose coefficient avoids its base. Let $K_i=U_iS_i$ be the two products of these assignments, split by weight and predecessor. Then $K_1K_2=N(N-1)$. Every part $z$ of row $i$ satisfies $U_i\mid z$, $S_i\mid z(z-1)$; all six parts exceed one. By the [two-congruence reconstruction](two_predecessor_congruences.md#recovering-the-row-from-a-product), weights with either $K_i\le X^k$ number $O(X^{k+o(1)})$ for $0<k<1$, using the complementary trinomial for the reconstruction. For the remaining weights both products exceed $X^k$, and the larger is at least $\sqrt{N(N-1)}$.

Work in a stratum on which the indicated rows and orientations are chosen deterministically (the actual parts can vary). Take the uniform law on the weights. Write $H=\mathcal H(N)/\log X$. Knowing $(U_i,S_i)$ leaves $O(1+X/K_i)$ possibilities for $N$ when a lower bound for $K_i$ is given. Since every part determines this label up to $X^{o(1)}$ possibilities and the row determines $N$, if $K_i\gg X^b$ then
$$ H\le\frac{1}{3\log X}\sum_{z\in{\rm row}\ i}\mathcal H(z)
          +\tfrac23(1-b)_++o(1).                                 \tag{3}
$$
The sum here counts all three positions. [Entropy proof](common_factor_count.md#information-that-every-individual-part-supplies)

Take $L_b=\lceil\log_b(2X+1)\rceil$. If a row avoids $b$, its three digit sums add to $s_b(N)$ by Kummer's criterion. The complement $M_b-1-N$ then makes four digit sums totaling $(b-1)L_b$.

**Both primes on one row.** Apply (2) to each of its three parts with shift zero. Apply it once more to the complement $M_2-1-N$ with shift $M_p-M_2$. Its entropy equals that of $N$. Adding the four bounds and using the digit sums gives
$$ H+\sum_{z\in{\rm row}}\mathcal H(z)/\log X\le 4\Delta+A+B+o(1). $$
Combine with (3) at $b=k$. This leaves $H\le\Delta+(A+B)/4+(1-k)/2+o(1)$. Balancing with the small-assignment estimate gives the exponent
$$                    \frac{1+2\Delta+(A+B)/2}{3}                \tag{4} $$
up to an arbitrarily small loss.

**Primes on different rows.** If the binary row has the larger assigned product, the [single-base four-word bound](four_digit_sums.md) gives exponent $h(1/4)+o(1)<.8113$. Otherwise the base-$p$ row has the larger product. Write $r=\mathbb E[s_2(N)]/L_2$ and $s=\mathbb E[s_p(N)]/((p-1)L_p)$. With $0<x,y<1$, (3) gives the first two bounds below; (2) applied to the two complements gives the third. Each has an additive $o(1)$:
$$
\begin{aligned}
 H&\le D_2(x)+(r/3)\log_2(1/x)+\tfrac23(1-k),\\
 H&\le D_p(y)+((p-1)s/3)\log_p(1/y),\\
 H&\le \Delta+(1-r)A+(1-s)B .
\end{aligned}
$$
In the one-base estimates, entropy of a digit is bounded by the log partition function plus its expected value times $\log(1/x)$ or $\log(1/y)$, then summed over digit positions. Take weights $a,b,1-a-b>0$ on these three inequalities. Let $J$ be the maximum of the resulting affine expression in $(r,s)$ over the four corners of $[0,1]^2$ after removing the term $2a(1-k)/3$. The bound is $J+2a(1-k)/3+o(1)$. Balancing with the small-assignment count costs exponent
$$                        \frac{J+2a/3}{1+2a/3} .                \tag{5}
$$
If both rows could take a given prime, either permissible assignment suffices. When $p\mid N(N-1)$ both 2 and $p$ must be taken, so the cases cover every possible failure.

For example at $p=3$ the verifier uses $(u,v)=(.4030,.5020)$ in (4) and obtains an upper bound below $.760957$. For (5) it uses $(a,b)=(.1158,.6570)$, $(u,v)=(.7793,.4443)$ and $(x,y)=(.2305,.4309)$, obtaining less than $.817550$. It checks (4) and all four corners in (5) with outward rational log bounds for every prime in the table. The strict margins absorb all little-oh terms, the finitely many orientations and the dyadic summation. The density follows by the Chinese remainder theorem: the excluded residues modulo each prime are exactly the $p-2$ classes different from zero and one.

~~~~sh
python3 check_prime_residue_counts.py
~~~~
