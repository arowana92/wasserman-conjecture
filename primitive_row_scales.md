# At most one possible weight for a pair of primitive rows

Fix two positive integer vectors $x,y$, each of length at least two and each with coordinate gcd one. Write $u=\sum x_i$, $v=\sum y_j$, $n=\operatorname{lcm}(u,v)$ and $c=\gcd(u,v)$. At the common weight $N=gn$, $g\ge1$, the actual rows are $A=(N/u)x$ and $C=(N/v)y$. The number $g$ is the gcd of all entries of the two actual rows together. Let $K_A,K_C$ denote their multinomial coefficients.

**All but at most one common weight have a prime dividing $N(N-1),K_A,K_C$.** There is a factorization-free formula for the only candidate. If $c>1$, every common weight works. For all pairs form
$$
 H_x=\gcd_{i<j}(x_i x_j),\qquad H_y=\gcd_{i<j}(y_i y_j).
$$
Let $H$ be the largest divisor of $\operatorname{lcm}(H_x,H_y)$ coprime to $n$. Repeatedly dividing out its gcd with $n$ computes it. For $c=1$, set $k$ equal to the least positive residue of $-H$ modulo $n$.

- If $k\nmid H$, every scale works.
- If $k\mid H$, the only possible weight is **$N=1+H/k$**. It passes the full set of lowest-layer assignment conditions at the weight and predecessor. Kummer tests at this one weight decide the pair for every scale.

For trinomial directions, the [endpoint theorem](all_weights_endpoint_reduction.md) therefore gives an equivalent formulation of the unrestricted conjecture. Range over pairs of primitive positive triples with coprime sums, keep the pairs with $k\mid H$, and test the displayed weight. For binomial directions a genuine exception occurs at $x=(5,18),\ y=(1,1)$. Here $n=46$, $H=45$ and the endpoint gcd at $g=1$ is one. All larger scales have a common endpoint prime.

### Recovering the exceptional weight

The identity $N\mid A_i K_A$ implies $u\mid K_A$. Likewise $v\mid K_C$, hence $c$ divides both coefficients and $N$. At a prime $p\mid N$, the lowest-layer condition requires one whole row to have all entries divisible by $p^{v_p(N)}$. This holds at every weight prime exactly when $c=1$.

At $p^e\parallel M=N-1$, the primitive coordinates and the actual ones differ by a unit modulo $p^e$. Each primitive vector has at least one $p$-unit. Thus $p^e\mid H_x$ holds exactly when all but one of its coordinates are multiples of $p^e$. The exceptional actual coordinate then has residue one because the row sums to $N$. The analogous alternative for $y$ proves that the predecessor conditions hold at every prime exactly when $M\mid H$.

If $c=1$ then $n=uv$ and $H\le u^2v^2/16=n^2/16$, by bounding the product of two positive parts in each primitive row. For any scale passing the conditions the complementary divisor $H/M$ is positive and less than $n$, since $M\ge n-1$. From $H=(gn-1)(H/M)$ it must equal $k$. Conversely, if $k\mid H$, then $n\mid H+k$ and $\gcd(n,k)=1$, so $N=1+H/k$ is a positive multiple of $n$. Its predecessor divides $H$, proving the claim. This uniqueness uses only the lowest layers and positivity.

Every prescribed $g\ge1$ can actually occur in this formula. Take any $u\ge3$, $h=gu$, $d=h+1$, $e=hd+1$ and $f=hde-1$. Use $x=(1,1,u-2)$ and $y=(ef,df,de)$. The numbers $d,e,f$ are pairwise coprime. Writing $v=ef+df+de$ gives $def=hv-1$, because $de-h(d+e)=1$. Also $v\equiv-1\pmod h$, so $u,v$ are coprime. Here $H_x=1$, $H_y=def=N-1$, and the sole scale passing the lowest layers is the prescribed $g$. Both coefficients are even there. In $y$ either all three entries are odd or two have the same positive 2-adic valuation, which forces a binary carry after every common scaling.

### Divisors at two different scales

Let $G_g=\gcd(N(N-1),K_A,K_C)$ at scale $g$, and write $M_g=gn-1$. For any two distinct scales $g<h$,
$$
                 G_g G_h\ge \frac{c^2 M_g M_h}{H(h-g)}. \tag{1}
$$
Indeed at any predecessor prime the [retained-power lemma](endpoint_pair_bound.md#retaining-partial-prime-powers) gives $M_g/\gcd(M_g,K_A,K_C)\mid H$. In particular the predecessor gcd is divisible by $M_g/\gcd(M_g,H)$. For the two removed divisors of $H$, their lcm is at most $H$ and their gcd divides $h-g$ since $M_g$ is coprime to $n$. Each endpoint gcd also contains $c$, coprime to its predecessor. This proves (1).

For triples, $H_x^2\mid x_1x_2x_3$. Prime by prime, the valuation of $H_x$ is the second smallest coordinate valuation and the smallest is zero. Hence $H\le(uv)^{3/2}/27$ by AM–GM. Since $M_h/(h-g)>n$, equation (1) gives the uniform bound
$$
                        G_gG_h>27\sqrt{uv}\,(g-1/n).
$$
For any $Y\ge1$, at most $\tau(H)(1+\lfloor Y/n\rfloor)$ positive scales have predecessor gcd at most $Y$. Indeed for $d=\gcd(M_g,H)$, the positive integer $M_g/d\le Y$ is fixed modulo $n$ by $d(M_g/d)\equiv-1\pmod n$. Each $d$ and quotient give just one scale.

### Removing the common scale from the digit bound

In every possible trinomial counterexample, let $g$ be the common gcd of the six nonendpoint entries. Then
$$
 (12g)^3\le N,\qquad 4(12g)^3\le N\quad(2\mid N),\qquad
 (24g)^3\le N\quad(2\mid N/g).
$$
More generally use the [digit-unit factor $\Lambda(N)$ and the maximum row product $H_p(N)$](digit_unit_products.md#a-test-using-all-primes-of-the-endpoint) defined there. If a prime base has digit sum below three, it is already universal. Otherwise
$$
                 27g^3 (N-1)^2\Lambda(N)
                       \le \min_{p\mid N(N-1)} H_p(N). \tag{2}
$$
The left side may be doubled when $2\mid N/g$.

To see this, $u,v$ are coprime in a counterexample. For each $p\mid N$, choose a row avoiding $p$, say $A$. Then $p\nmid u$ and $N/u$ contains the full $p$-power in $N$. The lowest occupied position in each part selects a distinct occurrence among the digit units of $N$, because the row adds without carries. Their product is divisible by the product of the three smallest digit units of $N$. After subtracting $3v_p(N)$ this gives precisely the factor of $\Lambda(N)$ in the primitive product $\prod x_i$. At $p\mid N-1$ the two row scales are units, so the factor of $(N-1)^2\Lambda(N)$ already occurs in the primitive product of the row avoiding $p$. Altogether $(N-1)^2\Lambda(N)\mid \prod x_i\prod y_j$. If $N/g$ is even, the binary avoiding row has odd primitive sum because it has a unique coordinate of smallest binary valuation. The opposite primitive row has even sum, hence has an even entry. This supplies one extra factor two beyond the contribution from the avoiding row.

The actual product of all six parts is $(gN)^3\prod x_i\prod y_j$. For each tested prime one row has product at most $H_p(N)$ and the other at most $N^3/27$. This proves (2). In particular $H_2(N)<N^3/32$ because the largest occupied bit exceeds $N/2$. At even $N$ the three lowest distinct binary positions in the avoiding primitive row give $\Lambda(N)\ge8$. At odd $N$ the two lowest units of the odd part of $N-1$ give $\Lambda(N)\ge2$. For $N\ge5$, $N^3/(N-1)^2<N+3$. These estimates give the three integer bounds up to an additive two on the right. For $g\ge3$ both sides are multiples of $g$, which removes that slack. For $g\le2$ the [verified initial range](full_range_1b.md) is more than sufficient.

[Lean checks the complementary-divisor congruence](PrimitiveRowScales.lean). Run `lean PrimitiveRowScales.lean` and `python3 check_primitive_row_scales.py` for direct prime-power assignments, exact coefficient and two-scale comparisons, and controls with prescribed scales.
