# Endpoint primes for every pair of binomial coefficients at almost every weight

Let $N\ge3$, $M=N-1$, and $C_i=\binom Ni$ for $0<i<N$.

**Theorem.** Outside $O(X^{0.993})$ weights $N\le X$, **every pair** $C_i,C_j$ has a common prime dividing $N(N-1)$. The quantifier over all pairs is inside the assertion at each weight. Thus at almost every weight every single cut of every positive trinomial row already supplies a binomial divisor that shares an endpoint prime with every positive trinomial coefficient.

There is also a growing-divisor result. Outside $O(X^{0.997})$ weights $N\le X$, for all $2\le i,j\le N-2$ the integer
$$P_{ij}=\prod_{p\mid NM}p^{\min(v_p(C_i),v_p(C_j))}$$
exceeds $N^{1/1000}$. This definition retains the full exponents in the two coefficients at the indicated primes.

Erdős and Szekeres proved that two interior binomial coefficients at every weight have a common prime. Bergman gave quantitative lower bounds for their gcd. Here the prime divisors in the conclusions must divide the weight or its predecessor. At $N=46$, the cuts 10 and 23 show that the endpoint restriction can fail. [Bergman, *On common divisors of multinomial coefficients*](https://arxiv.org/abs/0806.0607)

## Count when the modulus specified by a cut is small

Work on $Q/2\le N<Q=2^L$. In the first assertion consider a pair with endpoint gcd one. A cut of size one or its complement already shares a prime of $N$ with every interior coefficient: $N\mid i C_i$ and $N\nmid i$. Thus all four parts of a pair under consideration are at least two in both assertions.

For the general count write
$$D=\gcd(N,C_i,C_j),\qquad J=\gcd(M,C_i,C_j).$$
Use the [retained prime powers](endpoint_pair_bound.md) for the two two-part rows, choosing as first row one with smaller binary valuation. The assigned factors satisfy $DUV=N$, $JST=M$. The first row consists of multiples of $U$ and each of its parts is zero or one modulo every power in $S$. The second has the analogous properties at $V,T$. Write $K=US$. The following estimate adapts the [trinomial modulus count](endpoint_count.md) to two complementary parts:
$$\#\{N\text{ with such a pair},\ K\le Y,\ DJ\le Z\}
 \ll Q^{-1/3+o(1)}Y^{5/3}Z^{8/3}. \tag{1}$$
We use it for $1\le Y,Z\le Q$. Throughout, $Q^{o(1)}$ can be bounded by $O_\epsilon(Q^\epsilon)$, using the divisor bound and logarithmically many dyadic boxes.

Fix $D,J$ and a dyadic box $U\asymp U_0$, $S\asymp S_0$. Split $T=de$ so that the second row is $(Vdc,Vef)$ with positive $c,f$. Each assigned power divides one of the two parts and is coprime to $V$. Swap sides so that $d\le e$. Then
$$V(dc+ef)=JSde+1,\qquad dc+ef=DU,\qquad cf\le(DU)^2/(4T).$$
Here $T\asymp Q/(JS_0)$. The number of multiplier pairs $(c,f)$ is $O(Q^{-1+o(1)}D^2J U_0^2 S_0)$ by counting their possible products. A nonempty box has a product bound at least one.

* Choose $S,V,c,f$. The remaining factors satisfy
  $$(JSd-Vf)(JSe-Vc)=V^2cf-JS.$$
  Both factors are positive since multiplying them by $e,d$ gives respectively $Vdc-1,Vef-1$, and each actual part is at least two. The right side therefore has $Q^{o(1)}$ factor pairs, each determining $d,e$. The resulting box count is $Q^{o(1)}DJ U_0 S_0^2$.
* Choose $S,d,c,f$, using $d\le\sqrt T$. The integer $dc+ef$ divides $JSd^2c-f$, obtained by eliminating $e$ from the two linear expressions. This integer is positive since $JSd>Vf$ and $Vdc\ge2$. A divisor determines $e$, then $V$. This gives a box count $Q^{-1/2+o(1)}D^2J^{1/2}U_0^2 S_0^{3/2}$.

Use the geometric mean of the two upper bounds with exponents $1/3,2/3$. It is $Q^{-1/3+o(1)}D^{5/3}J^{2/3}(U_0 S_0)^{5/3}$. There are at most $Z(1+\log Z)$ pairs of positive $D,J$ with $DJ\le Z$, and $D^{5/3}J^{2/3}\le Z^{5/3}$. This proves (1).

## Recovery from binary digits at large modulus

Put $\kappa=159/200$, $l=\lfloor\kappa L\rfloor$, and suppose $K\ge2^l$. Each first-row part $z$ specifies $U,S$ in at most $Q^{o(1)}$ ways because $U\mid z$, $S\mid z(z-1)$, and $z\ge2$. For each pair, $N\equiv0\pmod U$, $N\equiv1\pmod S$. Once its high $L-l$ bits are fixed, the resulting class modulo $K$ contains at most one weight.

First suppose the binary valuation of the first row is zero. Fix the high binary string of $N$, with $j$ ones, and the number $k$ of its low ones. One of the two first-row parts uses at most $k/2$ low ones. Its high ones are a subset of the fixed $j$ positions. The number of possible parts is at most
$$2^j(4/3)^l3^{k/2},$$
by weighting a choice of $b$ low ones by $3^{-b}$ and summing the binomial expansion. Together with recovery and the unrestricted count $\binom lk$, this gives two bounds in this stratum.

Use the powers $3/4,1/4$ of the unrestricted and recovery bounds. Summing over high strings gives $(1+2^{1/4})^{L-l}$. For the low strings, each $[\binom lk 3^{k/6}]^{3/4}$ is at most $(1+3^{1/6})^{3l/4}$. Absorbing the $l+1$ choices of $k$ gives $Q^{E+o(1)}$, where
$$E=(1-\kappa)\log_2(1+2^{1/4})+
 \kappa\left(\tfrac34\log_2(1+3^{1/6})+\tfrac14\log_2(4/3)\right)<0.993.$$
The verifier proves the strict inequality by raising rational upper bounds to integer powers.

For endpoint gcd one the first row is odd in coefficient, since $2\mid NM$. At $K<2^l$ use (1) with $Z=1$; its exponent is $-1/3+5\kappa/3=119/120$. Both exponents are strictly below $0.993$. A geometric sum over dyadic blocks proves the first assertion.

## Allowing carries

Suppose $P_{ij}\le N^{1/1000}$, and put $\delta=1/1000$. We have $DJ\le Q^\delta$. The smaller of the two binary valuations is at most $\delta L$, since the common power of two divides $P_{ij}$. By Kummer's criterion the first binary addition thus has at most $\delta L$ outgoing carries.

Specify the set of outgoing carry positions. Outside these positions and their successors the two parts still partition the ones of the weight. On the remaining low positions one part uses at most $k/2$ ones; give the exceptional positions of that chosen part as additional binary data. The recovery count acquires a factor at most
$$\left(\sum_{b\le\delta L}\binom Lb\right) 2^{2\delta L}
\le (1025/1024)^L 2^{12\delta L}.$$
The last inequality weights each position by $1/1024$. Only the recovery estimate in the geometric mean acquires this factor, so its exponent increases by less than $0.0035$, using $(1025/1024)^{1000}<4$. At small $K$, (1) gives exponent $119/120+8\delta/3<0.997$. At large $K$ the exponent is below $0.993+0.0035<0.997$. This proves the growing-divisor assertion.

## Reproduction

The verifier checks the algebraic eliminations as polynomial identities, verifies the strict exponent bounds with integers, and independently computes both the binomial gcd classes and the first finite examples. Through weight 5,000 there are exactly ten weights admitting a pair with endpoint gcd one:
$$46,370,476,477,609,974,989,2509,2988,4324.$$
For comparison, replacing one binomial with a positive trinomial first allows failure at [weight 5,752](cut_relaxation.md).

~~~sh
python3 check_binomial_endpoint_density.py  # uses SymPy for the symbolic identities
~~~
