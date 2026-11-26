# Common endpoint divisors for binomial cuts

Put $M=N-1$. A cut of a multinomial row is the sum $a$ of a nonempty proper selection of its parts. The coefficient of that row is a multiple of $\binom Na$ by grouping the selected and unselected parts. We first work with two binomial coefficients, writing
$$G=\gcd\left(NM,\binom Na,\binom Nb\right).$$

* If each cut is a proper divisor of $N$, then $G>1$. If also $a,b>1$, then **$N-1<G^4$**. This holds for arbitrarily large denominators $N/a,N/b$.
* For arbitrary cuts, if $\gcd(a,M)=\gcd(b,M)=1$ then $G>1$. Thus in a [trinomial counterexample](all_weights_endpoint_reduction.md), one of the two nonendpoint rows must have $\gcd(c,M)>1$ and $\gcd(c-1,M)>1$ at **every** part $c$. Use both $c$ and $N-c$ as possible cuts.
* Suppose $a\ge b>1$, and their reduced fractions of $N$ have respective numerators $r,t$. Then $M<4G^4r^6t^5$. In a counterexample, one of the two rows must have **every proper cut** of reduced numerator greater than $(M/4)^{1/11}$. All such cuts exceed one by the [nonendpoint part bound](million_part_bound.md).
* At a positive limit of $M/G^4$ for divisor cuts with $a,b>1$ and unbounded $N$, the limiting value must be $1/j$ for an integer $j\ge2$. A family attains $1/2$ for the explicit divisor furnished by the factorial identities below.

Bergman’s Theorem 2 bounds the full gcd of two binomial coefficients, with no prescribed set of prime factors. Cui and Yang’s November 2026 preprint counts a single coefficient with specified gcd with its upper index. The factor $N(N-1)$ in the assertions above selects the primes that can divide the endpoint trinomial. [Bergman](https://arxiv.org/abs/0806.0607) · [Cui and Yang](https://arxiv.org/abs/2611.25794)

## Factorial identities for a cut

Write $a/N=r/s$ and $b/N=t/u$ in lowest terms, put $g=\gcd(s,u)$, $m=N/\operatorname{lcm}(s,u)=\gcd(N,a,b)$, and $F=r(s-r)$, $H=t(u-t)$. Then
$$g\,\frac{M}{\gcd(M,\operatorname{lcm}(F,H))}\mid G,\qquad
 G r t>m\gcd(F,H). \tag{1}$$

Indeed $\binom Na\,a/N$ and $\binom Na\,a(N-a)/(NM)$ are integers. The first gives the divisor $s$. In the second substitute $a=Nr/s$ and use $\gcd(sN,M)=1$ to obtain $M/\gcd(M,F)\mid\binom Na$. Applying this to the two cuts and taking gcds proves the divisibility. Bounding the denominator by $\operatorname{lcm}(F,H)$ gives the inequality since $Msu>N(s-r)(u-t)$. This last strict inequality follows from $(s-r)/s\le M/N$ and $(u-t)/u<1$.

## Bounding the weight by the common divisor

Let $R=\gcd(M,\operatorname{lcm}(r,t))$ and $W=GR$. Since $\operatorname{lcm}(F,H)$ divides $\operatorname{lcm}(r,t)(s-r)(u-t)$, (1) gives the positive integer
$$h=\frac{(W/g)(s-r)(u-t)}{M}.$$
Using $m su=gN$, subtract $mh$ from $W$ to obtain
$$k=W-mh=\frac{W}{gM}\big(m(ts+ru-rt)-g\big)>0.$$
Indeed $ts+ru-rt=t(s-r)+ru>g$. In particular $W>m$. The prime divisors of $R$ divide at least one of $a,b$, because the removed denominators divide $N$. If both cuts are coprime to $M$ then $R=1$ and $G>m\ge1$.

Assume $a\ge b>1$, so $ts\le ru$ and $mts/g=b>1$. Multiplying the expression for $k$ by $s$ gives
$$(ks-Wr)M=(W/g)(s-r)(mts-g)>0. \tag{2}$$
Also $ks<2Wr$. To see this, the difference of the numerators after clearing $W/(gM)$ is $ms(ru-ts)+mrts+g(s-2r)>0$, using $g\le s$ and $s>r$. Therefore $s<2Wr/k$. In (2) use $ks-Wr\ge1$, $m<W$, and $R\le rt$. This gives $M<4W^4r^2t\le4G^4r^6t^5$, as asserted.

For the divisor cuts $r=t=1$ we have $R=1$. Set $x=ks-G$ and $y=ku-G$. Both are positive by (2) and its symmetric version. Also $x<G$, $m\le G-1$ and $s\le G+x$. Thus
$$M< G(G-1)\frac{(G+x-1)(G+x)}x\le G^2(G-1)(G+1)<G^4.$$
The middle inequality follows by expanding the quotient as $G(G-1)/x+2G-1+x$ and using $(x-1)(G(G-1)-x)\ge0$. If one of the divisor cuts is 1 the strict inequality $G>1$ still follows from $W>m$, while the weight bound can fail. For instance $N=2p$, $a=p$, $b=1$ at an odd prime $p$ gives $G=2$ since $\binom{2p}{p}\equiv2\bmod p$.

## Recovering the divisor cuts at every weight

In the divisor case the equations give
$$xy=h(mG-kg),\qquad s=(G+x)/k,\qquad u=(G+y)/k,\qquad N=msu/g. \tag{3}$$
This follows by expanding $(ks-G)(ku-G)$ and substituting $k=W-mh$ with $W=G$. For a fixed $G$, loop over $1\le k<G$, $mh=G-k$, $g\mid G$, and all positive factor pairs $xy=h(mG-kg)$. Retain integral $s,u\ge2$ with $\gcd(s,u)=g$ and $N/u\ge2$ for $s\le u$. The actual gcd is checked by counting the borrowing layers of the two binomial coefficients at the prime powers of $N$ and $M$.

In particular, along a family with $M\ge cG^4$ for fixed $c>0$, (2) and $s<2G/k$, $m<G/h$ give $g h k^2 x<4/c$. These four positive integers are bounded. On a subsequence fix all four and let $G$ tend to infinity. Then $m=(G-k)/h$, $s=(G+x)/k$, and (2) gives $M/G^4\longrightarrow 1/(g h k^2 x)$. If the denominator were one then $g=h=k=x=1$, yielding $s=G+1$ and $u=G^2-1$, contrary to $\gcd(s,u)=1$.

To see the scope of the factorial bound in (1), take $w\ge5$ coprime to 6, $s=w+1$, $u=w^2-2$ and $m=(w-1)/2$. Then $g=\gcd(s,u)=1$, $N=msu$ and $M=w^2(w^2-3)/2$. Here $\gcd(s-1,u-1)=1$ and $\operatorname{lcm}(s-1,u-1)=w(w^2-3)$. The divisor on the left of (1) is exactly $w$, and $M/w^4\to1/2$. The actual binomial gcd may include additional factors.

## How many weights can have a small gcd?

Across all weights, the number of pairs of cuts $a,b>1$ with reduced numerators at most $A$ and endpoint gcd at most $B$ is
$$O_\epsilon\big(A^4 B^2 (AB)^\epsilon\big) .$$
In particular, for each $0<\delta<1/2$ and each $\epsilon>0$, at most $O_{\delta,\epsilon}(X^{1-2\delta+\epsilon})$ weights through $X$ admit two nontrivial proper divisor cuts with endpoint gcd at most $N^{1/2-\delta}$. This conclusion applies to arbitrary multinomial rows possessing such cuts.

For the count, retain $r,t,R,W,k,g,m,h$ as above. With $x=ks-Wr$ and $y=ku-Wt$ both positive, expanding their product gives
$$xy=h(mWrt-kg).$$
For fixed $r,t,G$ there are at most $O_\epsilon((AB)^\epsilon)$ choices of $R\mid\operatorname{lcm}(r,t)$ and $g\mid G$. After choosing $1\le k<W\le Grt$, factor the positive integer $W-k=mh$, and then factor $h(mWrt-kg)=xy$. All factored integers are polynomially bounded in $AB$. The divisor bound shows that the choices other than $k$ cost $O_\epsilon((AB)^\epsilon)$ altogether. Each choice specifies $s=(Wr+x)/k$, $u=(Wt+y)/k$, $N=msu/g$, and the two cuts. Summing $Grt$ over $G\le B$ and $r,t\le A$ proves the bound. For the divisor-cut consequence take $A=1$ and $B=X^{1/2-\delta}$.

## Reproduction

Run `python3 check_divisor_cuts.py 512 --direct-through 4096` with python-flint installed. Both cuts are proper divisors greater than one in the table. Repetitions are allowed and exchanging the two cuts counts once. Each bound covers every weight. The enumeration reconstructs 140,051 candidates through gcd 512. The program also checks (1) and the weight bounds on all 1,333,300 unordered pairs of unrestricted cut indices at weights 2 through 200, using exact binomial coefficients. The independent direct-weight search uses exact binomial coefficients through 4,096, covering the entire possible range for gcd at most eight. The largest-weight example in the full table was also checked with exact FLINT binomial coefficients.

| Maximum gcd | Pairs | Largest weight | Cuts at that weight | Actual gcd |
|---:|---:|---:|---:|---:|
| 8 | 4 | 176 | 16, 22 | 5 |
| 16 | 23 | 232 | 8, 58 | 11 |
| 32 | 84 | 17,120 | 160, 856 | 17 |
| 64 | 270 | 17,120 | 160, 856 | 17 |
| 128 | 719 | 3,063,736 | 4,628, 29,459 | 90 |
| 256 | 1,758 | 29,176,686 | 6,237, 147,357 | 190 |
| 512 | 4,156 | 115,980,239 | 5,899, 334,237 | 341 |

Beyond the last displayed weight, any two multinomial coefficients with such cuts have common endpoint gcd greater than 512. The maximum weight is sharp for the binomial formulation. For unrestricted binomial cuts there can be no endpoint prime at all. For example the endpoint gcds of the individual cuts 10 and 23 at weight 46 are 23 and 90, respectively.
