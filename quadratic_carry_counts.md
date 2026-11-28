# A sharper count at every positive quadratic

Fix positive integers $M,h$ and let $P$ range over all prime powers of positive exponent. The number of $P\le X$ for which Wasserman's trinomial conjecture can fail at $MP^2+hP+1$ is
$$
 O_{M,h}\!\left(X^{(2+\beta)/3}\log(2X)\right),\qquad
 \beta=\log_2\varphi,\quad \varphi=(1+\sqrt5)/2 .
$$
Here $(2+\beta)/3=0.89808\ldots$. This improves the [previous exponent](quadratic_prime_power_density.md) of $(5+\beta)/6=0.94904\ldots$. There is no restriction on the factorization of $MP+h$.

For an interval of $K$ consecutive inputs in $[X,2X]$ put $R=K/X^{2/3}$. The count is $O_{M,h}(X^{2/3}(1+R)^\beta\log(2+R))$. In particular it is $o(K)$ uniformly whenever $R\to\infty$.

If $h^2-4M$ is not a nonnegative square, there are at most $O_{M,h}(K^\alpha\log(2K))$ exceptions in **every arithmetic progression of $K$ positive inputs with odd common difference**, where $\alpha=\log_2\lambda=0.551463\ldots$ and $\lambda^3=\lambda^2+1$, $\lambda>1$. This count is independent of where the progression starts.

## Uniform dependence on the quadratic coefficients

Take two or three integral polynomials of degree at most two and include their sum as a further polynomial. All must be nonconstant. Write their derivatives as
$$
            F_i'(x)=2^{e_i}(u_i x+v_i),\qquad \gcd(u_i,v_i)\ \text{odd}.
$$
For $n\ge0$, the number of residues modulo $2^n$ extending to any integer where the parts are nonnegative and add without a binary carry is at most
$$
 C(n+1)\lambda_k^{\,n+\max e_i-\min e_i}.                 \tag{1}
$$
Here $C$ is absolute, $k$ is the number of parts, $\lambda_2=\varphi$ and $\lambda_3=\lambda$. The coefficients and offsets may have either sign.

At depth $l\ge2$, stop a residue ball $x+2^l\mathbb Z$ when every $s_i=v_2(u_i x+v_i)$ is less than $l-1$, using $v_2(0)=+\infty$. At most $2(k+1)$ balls per depth can still require splitting. Indeed, a primitive affine expression with even slope is odd everywhere. One with odd slope vanishes on precisely one class modulo $2^{l-1}$, hence two classes modulo $2^l$. Start with the four balls at depth two, treating $n<2$ directly. There are at most $4(k+1)$ newly stopped balls per depth.

For two distinct inputs $x_1,x_2$ in a stopped ball, the normalized derivative valuations stay fixed. Writing $a_i$ for the quadratic coefficient, $2a_i=2^{e_i}u_i$, so
$$
 v_2\big(F_i(x_2)-F_i(x_1)\big)=v_2(x_2-x_1)+e_i+s_i .
$$
In fact $v_2(a_i(x_2-x_1))\ge e_i+l-1>e_i+s_i$, including $a_i=0$. This proves the digit-toggle property on the ball. The delay spread of the parts and their total is at most $\max e_i-\min e_i+l-2$. The [weighted prefix argument](affine_predecessor_triples.md#a-binary-count-for-polynomial-rows) now bounds the number of passing prefixes in the ball by $\lambda_k^{n-l+\max e_i-\min e_i+l-2}$. For two parts the undelayed count is two and there is one delayed bit. The free-column sum is at most $1+\varphi^{-1}=\varphi$. For three parts the bound is $1+\lambda^{-2}=\lambda$, or at most one choice if all four bits are undelayed. Adding the bounded number of unstopped balls at depth $n$ proves (1). An odd-step progression of length $K$ has distinct residues modulo $2^{\lceil\log_2 K\rceil}$, so (1) also bounds its count.

## Applying the estimate

Use the [quadratic row classification](quadratic_prime_power_density.md#the-first-row). The row assigned the base of $P$ is one of finitely many polynomial triples with coefficients depending only on $M,h$. Its parts have positive degree. The two unshifted parts are positive multiples of $P$, and a constant shifted part would be 1, incompatible with the common scale $U\ge3$. If $U$ is bounded in terms of $M,h$, the companion belongs to another finite polynomial list after restricting the parameter to classes modulo $U$. Each of its parts is a positive rational multiple of the weight. In a counterexample at least one of these two rows has odd coefficient. On these finite lists use (1) with three parts. Restricting an odd-step progression to a fixed residue class and reindexing that class preserves oddness of the induced step. When the discriminant is not a nonnegative square these lists cover all candidates, proving the progression bound.

For square discriminant it remains to count an odd companion in the unbounded-scale case. The [companion cut](quadratic_prime_power_density.md#the-smallest-group-of-predecessor-powers) uses $d\ll_{M,h}X^{1/3}$ and $1\le m\ll_{M,h}d$ (the multiplier called $n$ there). For fixed $d,m$ there are $O_{M,h}(1)$ classes with modulus comparable to $d^2$. After dividing out a power of two from the nonconstant coefficients and rounding the constants down, the two parts of the cut are integral quadratics with nonzero derivatives. Carry freedom of the original cut implies carry freedom of the divided pair. The [derivative content calculation](quadratic_prime_power_density.md#the-derivative-bound-on-a-class) proves $\max e_i\le\max(v_2(d),v_2(m))+O_{M,h}(1)$ for the pair and its sum. Formula (1) therefore bounds the count per class by
$$
 O_{M,h}\!\left((\varphi^{v_2(d)}+\varphi^{v_2(m)})
        (1+K/d^2)^\beta\log(2+K/d^2)\right).
$$
Since $\sum_{j\le Z}\varphi^{v_2(j)}=O(Z)$, summing over $m$ and $d\asymp J$ costs $O_{M,h}(J^2(1+K/J^2)^\beta\log(2+K/J^2))$ per dyadic block. On halving $J$, the factor before the logarithm decreases by at least a factor $2^{2-2\beta}>1$, while the logarithm increases by at most $\log4$. The sum is bounded by a constant times the largest-scale expression, at $J\asymp_{M,h}X^{1/3}$. This absorbs the finite-list contribution $O_{M,h}(K^\alpha\log(2K))$ for $K\le X+1$, proving the interval result. Setting $K\asymp X$ and summing over dyadic intervals gives the stated global exponent.

Run `python3 check_quadratic_prime_power_density.py` for the underlying allocations and quadratic identities. The normalized derivative test compares the determined bits of 26,002 input prefixes directly across 1,212 stopped classes, including signed coefficients and different powers of two in the two derivative contents.
