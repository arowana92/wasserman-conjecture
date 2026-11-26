# Two close prime powers and a third near their product

Let $P<Q=P+u$ and $R=PQ-k$ be powers of three distinct primes, with $0<k<u^2$. Set $F=PQR$. For **every** positive integer $d$ with $200du^3\le P$, any two positive trinomial coefficients share one of the three base primes at either weight $dF$ or $dF+1$. The [endpoint reduction](all_weights_endpoint_reduction.md) consequently gives the three-row conjecture at both weights. The third power here has order $P^2$.

For example,
$$P=100000000003,\quad Q=100000000019,\quad R=10000000002200000000023$$
are prime, with $u=16$ and $k=34$. The same three primes work at both weights for all $1\le d\le122070$.

## The two smaller bases

It suffices to forbid every pair of bases for positive binomial rows at $n=dF$, by the [transfer lemma](clustered_predecessors.md). At $P$, the quotient $n/P$ is the value at $X=P$ of
$$d(X+u)\big(X(X+u)-k\big)=d\big(X^3+2uX^2+(u^2-k)X-ku\big).$$
The constant is negative and $u^2-k\ge1$. All coefficients have magnitude less than $P/4-1$, so a single borrow normalizes the constant. At $Q$, the quotient $n/Q$ is the value at $Y=Q$ of
$$d(Y-u)\big(Y(Y-u)-k\big)=d\big(Y^3-2uY^2+(u^2-k)Y+ku\big).$$
Here only the quadratic coefficient is nonpositive. In the [two-power digit test](signed_power_pairs.md), the bound from the nonconstant digits at $P$ is
$$L\le3\big(1+u+d(u^2(1+u)+2u(1+u)^2+(1+u)^3)\big)\le60du^3,$$
using $1+u\le2u$ and $u\ge2$. Thus $2L+2<Q$, and all the signed coefficients at $Q$ satisfy the test. This excludes the pair $P,Q$.

## The two pairs with the larger power

Write $W=P$ or $Q$ and suppose a positive binomial row avoids the bases of $W$ and $R$. Both parts are divisible by $WR$. After dividing by $R$ they sum to $dPQ=dR+dk$. Here $dk<W<R$ and $d<R$. Carry freedom at $R$ expresses their quotients as $e_iR+j_i$ with nonnegative $e_i,j_i$, $\sum e_i=d$ and $\sum j_i=dk$. Divisibility by $W$ requires $j_i\equiv ke_i\bmod W$ since $R=PQ-k$. As $0\le j_i,ke_i\le dk<W$, these are ordinary equalities. The original parts must therefore be $e_iF$, with $e_i>0$.

For $W=P$ the quotients $e_iF/P$ have a strictly negative constant coefficient $-e_i ku$ in the first polynomial above with $e_i$ as multiplier. For $W=Q$ the quotients $e_iF/Q$ have a strictly negative quadratic coefficient $-2e_i u$ in the second polynomial. All coefficients have absolute value less than $W/4-1$ since $e_i\le d$. In either case both parts have a digit greater than $W/2$ at the same place. This forces a carry. All three pairs are now excluded.

The example and small comparisons of this argument with factorial-floor valuations are reproduced by `python3 check_powers_near_a_product.py` with python-flint installed.
