# All scales of two polynomial families of trinomial rows

Let each of two triples of polynomials in $\mathbb Z[T]$ have positive values at all sufficiently large positive integers. At $T=t$, divide each evaluated triple by its coordinate gcd to obtain primitive directions $x(t),y(t)$ with sums $u(t),v(t)$. Assume at least one triple has nonconstant direction over $\mathbb Q(T)$, so at least one ratio of two of its coordinates varies. At each $t$, allow **every** common weight, with arbitrary positive integer scaling of either direction to that weight.

There are an effectively computable modulus and threshold with the following properties on each residue class beyond the threshold.

1. Either every common weight has a prime dividing both coefficients and $N(N-1)$, or there is just one scale to test. In the second case the scale $g=N/\operatorname{lcm}(u(t),v(t))$ is the value of an explicit rational polynomial in $t$. Its weight and all six actual entries are rational polynomials too. The alternatives classify exactly which scales pass all lowest-layer conditions, including when the candidate polynomial is nonconstant.
2. For some $\delta>0$ depending only on the two polynomial triples, **all scales work simultaneously for all but $O(X^{1-\delta})$ parameters $1\le t\le X$**. Here “work” means that the two coefficients and $N(N-1)$ have a common prime. This addresses the pair of specified directions at every allowed weight.

## Computing the polynomial scale

The gcd of any finite list of nonzero rational polynomials which take integer values on a progression is itself given, after evaluation and taking the positive integer gcd, by a rational polynomial on each of finitely many smaller progressions, for all sufficiently large arguments. To see this, clear a common denominator $D$, extract the positive-leading primitive polynomial gcd $F$ from the resulting integral polynomials and divide by $F$. The quotients have polynomial gcd one over $\mathbb Q$, so a polynomial Bézout combination with integral coefficients gives a nonzero integer $b$. Their evaluated gcd equals their gcd with $|b|$, which depends only on $t\bmod |b|$. The original gcd is this constant times $F(t)/D$ when $F(t)>0$. Products and lcms also have this property, by the gcd formula for an lcm. In that formula the polynomial gcd of the arguments divides their polynomial product, so the quotient on each progression is a polynomial. All steps are effective.

Apply this observation to produce $x,y,u,v$ on finitely many progressions. The condition $\gcd(u,v)=1$ is eventually constant after further subdivision. Indeed its value function is then a positive polynomial, either constant or tending to infinity. When the condition fails there is a common weight prime at every scale. On the other progressions put
$$
 n=uv,\quad H=\operatorname{lcm}\big(\gcd_{i<j}(x_ix_j),\gcd_{i<j}(y_iy_j)\big).
$$
Both are rational polynomials on the refined classes. At least one primitive sum tends to infinity because, after removing the common polynomial factor, its positive coordinates have at least one nonconstant polynomial and their common integer divisor is bounded. Thus $n\to\infty$. The [primitive row criterion](primitive_row_scales.md) gives the exact remaining condition $gn-1\mid H$. Using the unstripped lcm here is valid since $gn-1$ is coprime to $n$. For triples, $H\le n^{3/2}/27$ pointwise.

Divide polynomials over $\mathbb Q$, writing $H=Qn+R$ with $\deg R<\deg n$ or $R=0$. If a scale passes the test at $t$, put $k=H(t)/(g n(t)-1)$. Then $0<k\le H(t)/(n(t)-1)=O(\sqrt{n(t)})$, and
$$
                  gk-Q(t)=\frac{R(t)+k}{n(t)}.
$$
The right side tends to zero uniformly over all possible scales. The left side lies in a fixed discrete lattice since the coefficients of $Q$ have a common denominator. It must vanish beyond an effective threshold. Thus $k=-R(t)$, $R(t)<0$, and $g=Q(t)/(-R(t))$.

If $R$ is eventually negative, divide $Q$ by $-R$ over $\mathbb Q[T]$. For a nonzero remainder and nonconstant divisor, the proper rational fraction tends to zero and is eventually nonzero. The same discrete-lattice argument excludes all sufficiently large integral values of the quotient. A constant divisor already divides exactly. The only possible branch consequently has a polynomial $G$ with $Q=(-R)G$. Conversely, every large $t$ with $-R(t)$ and $G(t)$ positive integers gives
$$
                  H(t)=(-R(t))(G(t)n(t)-1).
$$
This proves sufficiency for the lowest layers as well. Signs stabilize and integrality is a congruence test using fixed denominators. All thresholds above can be bounded from polynomial coefficients. On a surviving class the two actual triples are $Gv x$ and $Gu y$.

In particular $\deg H<\deg n$ rules out **all** sufficiently large parameters and all scales. This includes polynomial directions in which, after removing their common polynomial factors, the three coordinates within each triple are pairwise coprime over $\mathbb Q[T]$. Their numerical pairwise gcds are bounded, so $H$ is bounded. Another case consists of two nonconstant primitive directions given by affine coordinates, with all three coordinate pairs nonproportional in at least one direction. The gcd products have degree zero for that direction and at most one for the other. Indeed at most one proportional pair is possible in a nonconstant affine direction. Here $n$ has degree two.

## Removing the remaining scales for almost every parameter

Restrict to a surviving progression and replace its parameter by its progression index. Split further to give the six actual part polynomials integral coefficients. For a rational polynomial that takes integer values, restriction modulo a common coefficient denominator achieves this because all of its new nonconstant coefficients are integral and its constant term is one of its integer values. The common weight is a nonconstant polynomial.

If either actual triple has two constant parts, the [two-part bound](endpoint_pair_bound.md) forces its coefficient and the other coefficient to share a divisor of $N(N-1)$ tending to infinity. Otherwise each triple has at least two nonconstant parts. Apply the [low-digit estimate for polynomial rows](polynomial_row_divisors.md#divisibility-with-a-power-saving-exceptional-set) at the fixed prime 2 to each. For $k=\lfloor\log_2 X\rfloor$ it gives $v_2(K)\ge k/4-2$ outside $O(X^{1-\delta})$ indices through $X$, for a fixed $\delta>0$. Both coefficients are therefore even away from the union of these two sets. Since $N(N-1)$ is even at every weight, this covers the surviving scale. The other scales already have a common endpoint prime. Sum over the finite progression classes.

### Parameters that are powers

For every fixed odd integer $b\ge3$, all common scales work at $t=b^a$ for all but $O_b(Y^{1-\delta})$ exponents $a\le Y$. Uniformly over odd bases $3\le b\le X$ and exponents $a\le Y$, at most $O(YX^{1-\delta})$ pairs can fail at even one scale. In particular, for all but $O(X^{1-\delta/2})$ odd bases up to $X$, **all scales work at all exponents $1\le a\le X^{\delta/2}$**. Here $\delta$ may be smaller than before, and constants may still depend on the initial polynomials. These estimates apply in particular when the bases are restricted to odd primes.

Indeed on a surviving progression $t=r+Lj$, the bad binary residues of the index $j$ occupy $O(2^{(1-\gamma)k})$ classes modulo $2^k$. The digit argument uses only the low $k$ positions, so the restrictions hold at arbitrarily large indices. Bad $t$ on that progression occupy at most as many classes modulo $2^{k+v_2(L)}$. For a fixed odd $b>1$ its multiplicative order modulo this power is $\gg_b 2^k$. This follows by successively factoring $b^{2^h}-1$, since $b^2\equiv1\pmod8$. Take $2^k$ comparable to $Y$ for the fixed-base estimate. For the uniform estimate take $2^k$ comparable to $X$. At exponent $a$, raising an odd residue modulo $2^m$ to the $a$th power has kernel of size at most $2\gcd(a,2^m)$ by the structure of the unit group. Also $\sum_{a\le Y}\gcd(a,2^m)\le(m+1)Y$ by counting multiples of each $2^i$. Summing the preimage bound, absorbing the logarithm by reducing $\delta$, and taking the finite union proves the assertions. The cases with two constant parts require only the same finite initial exceptions as before.

## An unbounded polynomial scale that passes the lowest layers

Take $x=(1,2,4)$ and, for $h=7t$, set $d=h+1$, $e=hd+1$, $f=hde-1$ and $y=(ef,df,de)$. The three factors $d,e,f$ are pairwise coprime. With $v=ef+df+de$, one has $def=hv-1$ and $v\equiv-1\pmod h$. Thus $n=7v$, $H_y=def=tn-1$ and $H_x=2$. The lcm is $H=tn-1$ for odd $t$ and $H=2(tn-1)$ for even $t$. Both branches recover $g=t$.

A single prime, 7, works at **every** $t$ and **every** scale. The first coefficient contains 7 because its primitive sum is 7. In the second actual row $7g y$, remove the common power of 7 dividing $7g$ and call the remaining unit $s$. The three low residues are $-s,-s,s$ modulo 7 since $d,e,f\equiv1,1,-1$. Their least nonnegative representatives sum to $14-(s\bmod 7)>7$, forcing a carry. Hence 7 divides both coefficients and their weight at every scale.

The [verifier](check_polynomial_row_scales.py) checks the polynomial divisions, including fractional and rejected branches, against direct integer recovery and tests the unbounded example by an independent factorial-floor valuation. Run `python3 check_polynomial_row_scales.py` with SymPy installed. The checks cover 4,708 scalar evaluations, 270 evaluated direction profiles and 540 direct divisibility tests.
