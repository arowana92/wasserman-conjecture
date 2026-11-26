# A small cut and the proportions of the opposite row

Let $N\ge3$, $M=N-1$, and $2\le a<N$. Let $B=(b_1,b_2,b_3)$ be a positive row of weight $N$. Write $K=\binom Na$, $L=N!/\prod_i b_i!$, $D=\gcd(N,K,L)$ and $J=\gcd(M,K,L)$. This formulation allows the other row to be specified by just one cut. If that binomial and $L$ share an endpoint prime, so do all trinomial refinements of the cut and $L$.

Put $\beta_i=b_i/N$, $\lambda=a(a-1)/M$, and take $\beta$ as a column vector. Define
$$F=\lambda\operatorname{diag}(\beta)+(a-\lambda)\beta\beta^{\mathsf T}.$$

**Theorem.** The matrix $D^2JF$ has integer entries and is positive definite. In particular, when $D=J=1$, both $F$ and its adjugate are positive definite integral matrices. They give the following two obstructions.

1. For every cut $b$ of $B$, every integer $u$, and every positive integer $t$,
   $$1\le a\left(\frac{tb}{N}-u\right)^2+
          \frac{a(a-1)}M t^2\frac bN\left(1-\frac bN\right). \tag{1}$$
   For instance, if $a(a-1)t^2\le2M$, every cut in the opposite row satisfies $|tb-uN|\ge N/\sqrt{2a}$. The other two entries need not be close to any prescribed proportions.
2. For every nonzero integer vector $v$ of length three, with $e=\sum v_i$,
   $$1\le a\lambda\sum_i \beta_j\beta_k(v_i-e\beta_i)^2
           +\lambda^2 e^2\beta_1\beta_2\beta_3, \tag{2}$$
   where $j,k$ are the complementary indices in each summand. Thus if a part in one row is at most $h$ and all three opposite positions satisfy $|v_iN-eb_i|\le Y$, a common endpoint prime is guaranteed whenever
   $$\frac{h^2(h-1)Y^2}{3MN^2}+\frac{h^2(h-1)^2e^2}{27M^2}<1. \tag{3}$$

In particular, for fixed $e$, a cut of size $a\asymp N^\alpha$ with $0<\alpha<1/2$ rules out an opposite row within $o(N^{3(1-\alpha)/2})$ of a fixed positive rational direction in the endpoint-coprime case. For a single opposite cut, (1) gives the scale $N^{1-\alpha/2}$. The bounds have content beyond the [individual-part threshold](extremal_parts.md) at $\alpha=2/5$.

## Integrality and the adjugate

Use the [retained-power argument](endpoint_pair_bound.md) at the binomial cut and the trinomial row. With $n=N/D$, $m=M/J$, the weight powers imply $n\mid ab_i$, and also $n\mid (N-a)b_i$. At every predecessor power assigned to the cut, one of $a,N-a$ vanishes. At a predecessor power assigned to $B$, at least one of $b_i,b_j$ vanishes for $i\ne j$. Consequently
$$n^2m\mid a(N-a)b_i b_j \qquad(i\ne j)$$
since $n,m$ are coprime. The off-diagonal entry $F_{ij}$ is exactly this numerator divided by $N^2M$. The row sums of $F$ are $a\beta_i$, whose denominators divide $D$, proving the diagonal assertion as well. The same proof works with two positive parts in $B$ for (1).

For real $c$ the quadratic form is
$$c^{\mathsf T}Fc=a(\beta^{\mathsf T}c)^2+\lambda\left(\sum_i\beta_i c_i^2-(\beta^{\mathsf T}c)^2\right).$$
It is strictly positive for $c\ne0$ since $\lambda>0$, all $\beta_i>0$, and only constant vectors have zero variance. For (1), put the label $t-u$ on the indices forming $b$ and $-u$ on the others. The stated bound follows from $\beta(1-\beta)\le1/4$ for $0\le\beta\le1$.

A diagonal matrix plus a rank-one matrix gives
$$\det F=a\lambda^2\beta_1\beta_2\beta_3,\qquad
F^{-1}=\lambda^{-1}\left(\operatorname{diag}(1/\beta_i)-(1-\lambda/a)\mathbf1\mathbf1^{\mathsf T}\right).$$
Multiplying these expressions and using $\sum_i(v_i-e\beta_i)^2/\beta_i=\sum_i v_i^2/\beta_i-e^2$ gives exactly the right side of (2) as $v^{\mathsf T}\operatorname{adj}(F)v$. This is a positive integer when $D=J=1$. Use $\sum\beta_j\beta_k\le1/3$ and $\prod\beta_i\le1/27$ for (3). For arbitrary $D,J$, the corresponding lower bounds in (1) and (2) are $1/(D^2J)$ and $1/(D^2J)^2$ respectively.

## A whole-weight test at unequal prime-power remainders

Let $P=p^i,Q=q^j$ have distinct prime bases dividing $NM$ and positive exponents. Suppose
$$N=dP+h=eQ+k,\qquad 0\le h<P,\quad 0\le k<Q,\quad d,e\ge1,\quad s_p(d)\le2.$$
Assume $h\ge2$, use $Y=ek$ in (3), and require the strict inequality there. Then the full conjecture holds provided either of the following tests holds.

- There are powers of $p,q$ dividing $N$ whose product exceeds $h$.
- $\gcd(d,e)=1$, $de(h+k)<N$, and $h^7\le6^6 M^4$.

For $h<2$ the conjecture already follows from the individual-part bound. In a counterexample the [endpoint reduction](all_weights_endpoint_reduction.md) leaves two nonendpoint rows covering the bases of $NM$. At $P$, a carry-free row distributes at most two upper digit units. It has a positive part $a\le h$, necessarily at least two. Under the first test this row cannot also avoid $q$ since all parts would be divisible by the specified product. Under the second, if it avoided both bases, write $r_l,v_l$ for its upper blocks at $P,Q$. Then $|d A_l-r_lN|\le dh$, $|e A_l-v_lN|\le ek$, giving $|(er_l-dv_l)N|\le de(h+k)<N$. Coprimality concentrates all upper blocks in one position; the two other parts sum to at most $h$, contrary to the [pair-sum bound](endpoint_pair_bound.md).

The opposite row now avoids $q$. Its integer upper blocks $v_i$ sum to $e$, and $|v_iN-eb_i|\le ek$ since its lower blocks sum to $k$. Condition (3) excludes this pair of rows. For fixed $d,e$, the inequality holds throughout $h=o(\sqrt N)$, $k=o((N/h)^{3/2})$ for $h\ge2$.

## Infinite families of consecutive weights

Fix a positive integer $e$ coprime to three and any fixed modulus $L\ge1$. For every sufficiently large $X$, there are $\gg_{e,L} X^{11/120}$ distinct starting weights $n\le X$, with $L\mid n$, for which the conjecture holds at **both $n$ and $n+1$**. Each pair has the same two upper powers $P=2^i$ and $Q=q^3$. At both weights the remainder above $3P$ has order $n^{5/12}$ and the remainder above $eQ$ has order $n^{17/20}$.

Choose $B=3\cdot 2^i$ in $[X/4,X/2]$. Use primes $q$ with
$$B-2B^{17/20}\le e q^3\le B-B^{17/20}.$$
Their interval has length $\asymp_e B^{11/60}$ at size $\asymp_e B^{1/3}$. The prime-gap theorem of Baker, Harman and Pintz gives a prime in each subinterval of length $C_e B^{7/40}$ for a sufficiently large fixed $C_e$. Taking separated subintervals gives $\gg_e B^{1/120}$ primes. [Baker, Harman and Pintz, Theorem 1](https://doi.org/10.1112/plms/83.3.532)

For each $q$, take all multiples $n$ of $qL$ with $B^{5/12}\le n-B\le2B^{5/12}$. There are $\gg_{e,L}B^{1/12}$ per prime for large $B$. A weight has at most three of these large prime divisors since the product of four exceeds $n$. For each pair $n,n+1$ the prime $q$ divides the weight or its predecessor, and the binary upper quotient is three with two digit units. The second same-row test holds at both weights. In (3) the first term is $O_{e}(B^{-1/20})$ and the second is $O_e(B^{-1/3})$. This proves the count and the bounds on both remainders.

## Integer examples and verification

Set $P=2^{420}$, $B=3P$, $e=10$, $H=2^{175}$, and
$$q=933053352955195767836684502058363471888193.$$
For every $n$ with $B+H\le n<B+2H$, $210q\mid n$, and $143\mid n+1$, both $n,n+1$ satisfy the conjecture. There are **1,709,179 disjoint consecutive pairs**. The verifier checks the extrema in (3) with integer arithmetic and certifies the primality of $q$ recursively.

For a smaller example put $P=2^{192}$ and
$$n+1=18831305206160042291507368269623000702000273884576598888581.$$
The same test works at both weights with $q=12348858302118435931$. At $N=n+1$ the largest prime power dividing $N-1$ is $q$, and $2q<h$. Thus the [single predecessor power test](two_power_cut_windows.md) does not cover the binary remainder at this weight.

The verifier checks 34,206 binomial–trinomial pairs through weight 36. It recovers the common endpoint valuations from exact coefficients, compares the quadratic forms with direct cofactor arithmetic, and checks the adjugate identity as a polynomial identity over the integers. Run it with Python, SymPy and python-flint.

~~~sh
python3 check_small_cut_forms.py
~~~
