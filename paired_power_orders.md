# Carry tests and multiplicative orders at two predecessor powers

Write $N=wP Q+1$, where $P=p^e$ and $Q=q^f$ are powers of distinct primes, $e,f\ge1$ and $\gcd(w,pq)=1$. Suppose $w$ is a prime power, twice an odd prime power, 12, 15, or 1. The same conclusions hold at any other $w$ for which the two bases are forced onto one row.

In a counterexample there is an integer $a$ with $1\le a<w$ such that **$aQ,(w-a)Q$ add without a carry in base $p$, and $aP,(w-a)P$ add without a carry in base $q$.** In fact the [four-factor inequality](two_prime_powers.md) forces one nonendpoint row to avoid both primes, with some part equal to $aP Q$. Its complement is $(w-a)P Q+1$. Group the other two parts together. In base $p$ the units digit 1 lies below the $P$ shift, and in base $q$ it lies below the $Q$ shift. Kummer's criterion gives the two assertions.

Consequently, if $dQ=p^h+c$ with $h\ge1$, $1\le wc\le d$, and $\gcd(p,d)=1$ with $d>1$, any counterexample requires
$$
 \operatorname{ord}_d(p)>h,\qquad h\le 40 w^{3/2}\sqrt d .
$$
These bounds are independent of the exponent in $P$. They also apply after exchanging the two prime powers and their bases. At $c=1$ we have the stronger asymptotic bound
$$h\ll_{p,w,\epsilon} d^\epsilon\qquad(\epsilon>0).$$
In particular, for fixed $p,w,A>0$ with $p$ prime and $w$ as above, all these weights with $w\le d\le h^A$ and $dQ=p^h+1$ satisfy the full conjecture for sufficiently large $h$. The threshold is uniform over $e$, $q$ and $f$. For prime $d$, the last section gives this conclusion throughout $1\le wc\le d$, even with $p$ varying.

## Counting in two intervals of exponents

Fix distinct primes $p,q$ and a cofactor $w\ge2$ as in the first paragraph. If $e$ varies over $E$ consecutive positive integers and $f$ over $F$ consecutive positive integers, the number of possible exceptional pairs for $wp^e q^f+1$ is
$$O_{p,q,w}\big(E^{\theta_q} F^{\theta_p}\big),\qquad
 \theta_r=\log_{B_r}\!\left(B_r-\lfloor B_r/w\rfloor\right)<1,\qquad
 B_r=r^{\lceil\log_r w\rceil}.$$
The estimate is uniform in the locations of both intervals.

To see this, consider $ax,bx$ adding without carries in base $r$, where $a+b=w$. Group input digits into blocks of radix $B=B_r$. No block of $x$ can have value $B-s$ with $1\le s\le\lfloor B/w\rfloor$. In multiplying $x$ by $a$, the incoming carry $i$ lies in $[0,a-1]$. At that digit the output is $B-as+i$, since $1\le as-i<B$. For $b$ the output is $B-bs+j$ for some $j\in[0,b-1]$. The sum is at least $2B-ws\ge B$. A block overflow contradicts the carry-free addition in base $r$. This argument excludes the same blocks for every split $a+b=w$.

For $s_0$ a fixed prime different from $r$, the order of $s_0$ modulo $B^j$ is at least a positive constant times $B^j$, by lifting. For $r=2$ start at modulus 4. Thus over any interval of $L$ consecutive exponents the residues $s_0^i\bmod B^j$ have bounded multiplicity if $j=\lfloor\log_B L\rfloor$. There are at most $(B-\lfloor B/w\rfloor)^j$ allowed strings. Apply this first to $q^f$ in base $p$ and independently to $p^e$ in base $q$, using the two necessary additions at the start of this note. [Prime-power order calculation](prime_power_ray_density.md#counting-exponents)

**At factor five there is an absolute bound uniform in the prime.** For every odd prime $\ell$ and every pair of intervals of $E$ positive exponents $e$ and $K$ positive exponents $k$, the number of possible exceptions at $5\cdot2^k\ell^e+1$ is
$$O(E^{0.598}K^{0.504}),$$
with an absolute implied constant. The factor-five argument leaves only $\ell>5$ with $\ell\equiv3\bmod8$. Put $Q=\ell^e$ and $T=Q2^k$. The binary test for $3Q,2Q$ depends only on $e$. Since $v_2(\ell^2-1)=3$, the order of $\ell$ modulo $2^j$ is exactly $2^{j-2}$ for $j\ge3$. This makes the bounded-multiplicity estimate absolute. For each remaining $e$ the quinary test for $3T,2T+1$ restricts $k$. Multiplication by the unit $Q$ preserves the order-based multiplicity bound, and $2$ has order $4\cdot5^{j-1}$ modulo $5^j$.

For the binary count use the [four-state matrix](factor_five_order.md#distinct-blocks-of-the-binary-expansion), whose characteristic polynomial is $x^4-x^3-x^2+x-1$. For the quinary count keep the two multiplication carries, starting from $(0,1)$ to include the constant term. In state order $(0,1),(0,0),(1,1),(2,1),(1,0)$, the successor lists, with indices starting from one, are $(2,3),(2),(2,3,4),(2,3,4,5),(2,4)$. Indeed an input digit $t$ at $(u,v)$ produces $3t+u$ and $2t+v$. Admit the digit if their residues sum to less than 5 and pass both quotients to the next state. The characteristic polynomial is $x(x-1)(x^3-2x^2-x+1)$.

The spectral radii of these nonnegative matrices are less than $1513/1000$ and $2247/1000$, respectively. The degree-four and degree-three factors are positive at the indicated rationals and strictly increasing thereafter, by differentiation. Thus the numbers of length-$j$ words are $O((1513/1000)^j)$ and $O((2247/1000)^j)$. The exact power comparisons $(1513/1000)^{1000}<2^{598}$ and $(2247/1000)^{1000}<5^{504}$ yield the stated exponents. The two digit tests and constants are independent of $\ell$.

## A gap in the orbit

The [order-layer lemma](factor_five_order.md#a-carry-at-the-order-modulus) with $a+b=w$ gives the order bound. For $0\le j<h$, put $x=(p^j\bmod d)/d$. The fractional part of $Q/p^{h-j}$ is obtained from $x$ by adding $c/(dp^{h-j})$ modulo one. The increment is at most $1/(2w)$. If $1-1/w<x<1-1/(2w)$, the sum is strictly between $1-1/w$ and 1. At this layer the two fractional parts for $aQ,bQ$ thus sum to $2-w(1-x-c/(dp^{h-j}))>1$. This forces a carry. The orbit segment therefore avoids the displayed interval. In integer form take $m=p^{h-j}$, $v=p^j$ and $t=d-(v\bmod d)$. In the gap, $(a+b)t<d$ and $c<tm$, while $v+t$ is a multiple of $d$. The Lean theorem applies these inequalities directly to $dQ=vm+c$. [Formal gap and order lemmas](CarryOrder.lean)

Take the set $B=\{p^j\bmod d:0\le j\le\lfloor(h-1)/2\rfloor\}$. The elements are distinct, $|B|\ge h/2$, and every product of two elements still belongs to the segment. For $n\ne0$ integer write $g=\gcd(n,d)$. Cauchy–Schwarz followed by additive-character orthogonality gives
$$
 \left|\sum_{x,y\in B}\exp(2\pi i nxy/d)\right|^2
 \le |B|d\,\#\{y,z\in B:d\mid n(y-z)\}
 \le d g |B|^2.
$$
The last inequality holds since there are at most $g$ residues in any fiber modulo $d/g$.

Use the periodic tent of height one supported on the forbidden interval, with half-width $\delta=1/(4w)$. Its mean is $\delta$ and the absolute values of its other Fourier coefficients are at most $\min(\delta,1/(\pi^2\delta n^2))$, by integrating the triangle or convolving two intervals. Its Fourier series converges absolutely. The tent vanishes on every product, so
$$
\delta |B|\le 2\sqrt d\sum_{n\ge1}\sqrt n\,\min(\delta,1/(\pi^2\delta n^2)).
$$
Split at $1/\delta=4w$. The terms through that point sum to at most $\delta^{-1/2}$ and the tail to at most $2\delta^{-1/2}/\pi^2$ by integration. Thus $|B|< (5/2)\delta^{-3/2}\sqrt d=20 w^{3/2}\sqrt d$.

For the factor-five rows of the [binary split](factor_five_order.md#the-necessary-binary-split), a sharper explicit estimate is **$h<164\sqrt d$**. Here $3Q$ has no bit in common with $2Q$. Thus $Q$ cannot contain three consecutive ones. At the middle and final one of such a run the preceding bit is 1. An incoming multiplication carry at the middle gives an overlap there, while no incoming carry creates one at the final place. Since $\lfloor 2^h/d\rfloor=Q-1$ and $Q$ is odd, the first $h-1$ digits of $1/d$ coincide with the padded high digits of $Q$. Thus $(2^j\bmod d)/d$ cannot be in $[7/8,1)$ for $0\le j\le h-4$. For $h\ge4$ take $B$ with indices through $\lfloor(h-4)/2\rfloor$, giving $h\le2|B|+3$. In the tent bound use $\delta=1/16$ and split at $n=5$. The sum on the right, including the factor 2 but omitting $\sqrt d$, is bounded by $2(5\sqrt5/16+32/(\pi^2\sqrt5))<5$. Hence $|B|<80\sqrt d$. Smaller $h$ also satisfy the claim.

## Uniform exclusion for power-sized orders

Fix $p,w$ and let $L=2+\lceil\log_p w\rceil$. Adding $aQ,bQ$ without a carry forbids a run of $L$ maximal digits in $Q$. Indeed in multiplying $Q$ by $a$, the incoming carry lies between 0 and $a-1$. After $L-1$ maximal input digits starting with carry $u$, the carry is $\lfloor(a(p^{L-1}-1)+u)/p^{L-1}\rfloor=a-1$, since $0<a-u<p^{L-1}$. The next output digit is $p-1$. The same calculation holds for $b$. This applies also if the run starts at the lowest input digit.

We give the argument slightly more generally when $\operatorname{ord}_d(p)\le Rh$ for fixed $R$. Assume $p\nmid Q$ as in the two-prime-power setting. The first $h$ base-$p$ digits of $1/d$ form $Q-1$, padded on the left. Its first $h-1$ digits agree with $Q$ since subtracting 1 does not borrow. Choose a fixed positive integer $v$ so that $(1-p^{-L})^v<1/R$, and put $D=vL$. The length-$D$ strings that differ from the all-maximal word in each of $v$ disjoint blocks occupy at most this proportion of the base-$p$ intervals of length $p^{-D}$.

Let $H=\langle p\bmod d\rangle$. At least $h-D$ elements of $H$, divided by $d$ as least residues, belong to that finite union of intervals, by taking successive blocks among the first $h-1$ digits. Suppose that along a sequence with $d\to\infty$ we had $h\ge d^\epsilon$. The normalized residues of $H$ are equidistributed on this sequence by Bourgain's subgroup bound for arbitrary moduli. To check its primitive-frequency hypothesis, fix $n\ne0$, set $g=\gcd(n,d)\le |n|$, and project $H$ to a subgroup $H'$ modulo $d/g$. The image has size at least $|H|/g\ge(d/g)^{\epsilon/2}$ for large $d$. Its normalized exponential sum at the unit frequency $n/g$ equals that of $H$ at $n$ because all fibers have the same size. Bourgain's bound makes these normalized sums tend to zero. Weyl's criterion gives equidistribution, contradicting the lower limiting proportion at least $1/R$ in our finite union. [Bourgain's original subgroup estimate](https://doi.org/10.1007/BF02807410) · [The author's formulation for arbitrary modulus, remark (2) after Corollary 3](https://comptes-rendus.academie-sciences.fr/mathematique/item/10.1016/j.crma.2007.01.019.pdf)

This proves $h\ll_{p,w,R,\epsilon}d^\epsilon$ whenever the order is at most $Rh$. For $c=1$ and $d>2$, the order divides $2h$ and exceeds $h$, so it equals $2h$. Bounded $d$ are already covered by the elementary order bound.

## Prime and fixed-support denominators at every offset

Fix $w$ and $\eta>0$. There are $\delta=\delta(\eta)>0$ and a threshold depending only on $w,\eta$ such that, above that threshold for $d$, a counterexample requires $h<d^\eta$ under either of the following hypotheses.

* $d$ has a prime factor $r>d^{1-\delta}$.
* $C_p(d)\le d^\delta$, where $C_p(d)=\prod_{r\mid d}r^{\min(v_r(d),E_r)}$. The product is over primes, $E_r=v_r(p^{r-1}-1)$ for odd $r$, and $E_2=v_2(p^2-1)$ if $2\mid d$.

Here $1\le wc\le d$ is the only size restriction on the offset in $dQ=p^h+c$, and $p$ may vary. In particular the bound holds for **every prime $d$**, uniformly in $p,c$. It also holds for $d$ supported on a prescribed finite set of primes with $p$ fixed, since then $C_p(d)$ is bounded. Thus, in either family, all $d\le h^A$ satisfy the conjecture for sufficiently large $h$, for any fixed $A>0$.

It suffices to consider $0<\eta<1$, since $h<\operatorname{ord}_d(p)\le d$. For any divisor $u>1$ of $d$, write $T_u=\operatorname{ord}_u(p)$. We have $T_u\ge u/C_p(d)$. Indeed for odd $r^v\mid u$, the valuation of $(p^{r-1})^{T_u}-1$ is $E_r+v_r(T_u)$ by lifting. If $r=2$, $p$ is odd, $p^2\equiv1\bmod8$, and the valuation of $(p^2)^{T_u}-1$ is $E_2+v_2(T_u)$. In either case $r^{\max(0,v-E_r)}$ divides $T_u$. Multiply these divisors over the distinct primes of $u$.

Take $\gamma=\eta/4$, and let $\sigma,\tau>0$ and $K\ge1$ be the projection threshold, cancellation exponent and number of factors in Bourgain's multilinear estimate for $\gamma$. Set $\delta=\min(\eta/4,\sigma/8,1/8)$. Choose a fixed integer $J$ so that truncating the Fourier series of the gap tent to $|n|\le J$ gives uniform error smaller than half its mean. This is possible by absolute convergence and only depends on $w$. [Bourgain, Theorem 2](https://www.numdam.org/item/10.1016/j.crma.2007.01.019.pdf)

Suppose $h\ge d^\eta$. Form $B$ from the consecutive powers $p^j$ for $0\le j\le\lfloor(h-1)/(K(J+1))\rfloor$, reduced modulo $d$. Its size $L$ is at least a constant times $d^\eta$. All products of $K$ elements belong to the segment avoiding the gap. For each $1\le |n|\le J$ set $g=\gcd(n,d)$, $m=d/g$. Reduction of the cyclic group generated by $p$ from modulus $d$ to $m$ has kernel of size at most $g$, so $T_m\ge T_d/g>h/g$. The $L$ terms of $B$ remain distinct modulo $m$.

For every $u\mid m$ larger than $m^\sigma$, check that the fraction of $B$ in any one residue modulo $u$ is less than $u^{-\gamma}$ for sufficiently large $d$. In the first case $d/r<d^\delta<m^\sigma$, so $r\mid u$ and the fiber contains at most $d/u<d^\delta$ elements. Dividing by $L$ gives the assertion. In the second case $C_p(d)\le d^\delta<u^{1/2}$, so $T_u>u^{1/2}$. Consecutive exponents put at most a fraction $1/T_u+1/L$ in any fiber. This is again smaller than $u^{-\gamma}$ for large $d$.

Apply Bourgain's theorem modulo $m$ with all $K$ sets equal to the image of $B$, at the unit frequency $n/g$. The normalized sum of $\exp(2\pi i n x_1\cdots x_K/d)$ over $B^K$ is at most $m^{-\tau}$ in absolute value. Thus the average of the truncated Fourier polynomial over these products tends to its positive constant term. Its difference from the tent average is less than half that constant, yet the tent vanishes on every product. This contradiction proves the assertion with uniformity in the stated parameters.

Run the finite checks with `python3 check_paired_power_orders.py` and python-flint. The program computes both characteristic polynomials, checks the rational power bounds, and compares the order, gap and block tests with exact binomial coefficients. Lean 4.12.0 checks the order and gap lemmas with `lean CarryOrder.lean`.
