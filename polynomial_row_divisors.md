# Divisors of polynomial trinomial rows and exceptional prime powers

Let $H(X)$ be a nonconstant polynomial with nonnegative integer coefficients and $H(0)=1$. Suppose $H=J^e$ with $J$ irreducible over the rationals, integral with positive leading coefficient, and $e\ge1$. Write $P=p^a$, where $p$ is prime and $a\ge1$. There is an effective $\delta_H>0$ with the following properties.

- For **every fixed prime** $p$, at most $O_{H,p}(K^{1-\delta_H})$ exponents $a\le K$ can fail Wasserman's full three-row conjecture at $H(p^a)$.
- Uniformly for $X\ge2$ and integers $Y\ge1$, the number of exceptional pairs $(p,a)$ with $p\le X$, $a\le Y$ is $O_H(YX^{1-\delta_H})$. In particular, for all but $O_H(X^{1-\delta_H/2})$ prime bases up to $X$, **every** exponent $1\le a\le X^{\delta_H/2}$ works. The exceptional set has relative density zero among primes.
- In each estimate, outside the same exceptional set **every row whose coefficient is coprime to $p$ has coefficient divisible by $J(P)$**. When $e=1$, this divisor is the entire weight.
- If the conjecture holds at $H(1)\ge3$, there is an effectively computable modulus $L$ such that **every sufficiently large prime $p\equiv1\pmod L$ works at every exponent**. This has positive relative prime density. When $H(1)<3$ no carry-free positive row can split its coefficient units into three parts, and every sufficiently large prime power works directly.

In particular, for any fixed positive $d$ and any distinct shifts $0=h_0<h_1<\cdots<h_{n-1}$ with $n\ge2$, all the density bounds apply to $H(P)=1+d\prod_j(P+h_j)$. The factorization lemma for [translated products](translated_products.md) gives $e=1$ or $2$. There is no size or factorization condition on the positive translates.

## A fixed multiplier for arbitrary polynomial rows

Let $B_1,B_2,B_3\in\mathbb Z[X]$ be positive at all sufficiently large integers and put $H=B_1+B_2+B_3$. Here we allow any nonconstant $H$ and signed coefficients. Let $K(t)$ be the trinomial coefficient with parts $B_i(t)$ whenever they are positive. For $j=0,1,2,3$ take $G_j$ to be the primitive positive-leading gcd over $\mathbb Q[X]$ of the three polynomials $B_i-[i=j]$ for $1\le i\le3$. In particular $G_0$ is the gcd of the parts. Put
$$D_0=H/G_0,\qquad D_1=(H-1)/(G_1G_2G_3).$$
These polynomials are integral by Gauss's lemma. Each $G_j$ for $j>0$ divides $H-1$. Two distinct such polynomials are coprime, as their common divisor divides the difference 1 at one coordinate. Their values at the same integer are coprime for the same reason.

There is an effectively computable nonzero integer $C$ such that **$D_0(t)D_1(t)\mid C K(t)$ at every positive row**. For each $j$, the three quotients by $G_j$ have polynomial gcd 1. Clear denominators in a polynomial Bézout identity to express a nonzero integer $C_j$ as an integral polynomial combination of the quotients. Take $C=\prod_j|C_j|$.

Write $b=v_r(K(t))$. If $r^u\parallel H(t)$ and $b<u$, then $r^{u-b}$ divides each $B_i(t)$, since $B_i(t) K(t)/H(t)$ is integral. Thus $u-b\le v_r(C_0G_0(t))$. If $r^u\parallel H(t)-1$ and $b<u$, use the two-slot factorial identities. Both
$$\frac{B_i(B_i-1)K}{H(H-1)},\qquad \frac{B_iB_j K}{H(H-1)}\quad(i\ne j)$$
are integers at this row. The first forces each part to be 0 or 1 modulo $r^{u-b}$. The second forces at most one unit, and the sum of all parts is 1 modulo this power. Thus $r^{u-b}\mid C_jG_j(t)$ for the unique unit slot $j$. The coprimality of $H(t)$ and $H(t)-1$ finishes the divisibility proof.

## Low digits of a polynomial

Fix a prime $r$ and a nonconstant integral polynomial $f$. There is an effective $\gamma>0$ such that, among residues $z\bmod r^k$, at most $O_{f,r}(r^{(1-\gamma)k})$ have
$$\left|s_r(f(z)\bmod r^k)-(r-1)k/2\right|>(r-1)k/12.$$
Here $s_r$ sums the digits of the least nonnegative representative.

First consider a nonzero polynomial $g$ of degree $m\ge1$ and leading coefficient $c$. The solutions of $g(z)\equiv0\pmod{r^j}$ occupy at most $m$ classes modulo $r^s$, where $s=\lceil(j-v_r(c))/m\rceil>0$. Otherwise choose $m+1$ solutions pairwise distinct modulo $r^s$. Lagrange interpolation writes the leading coefficient as the sum of $g(z_i)/\prod_{l\ne i}(z_i-z_l)$. Every summand has valuation at least $j-m(s-1)>v_r(c)$, a contradiction.

For a uniformly chosen output modulo $r^k$, the $k$ digits are independent uniform variables. Hoeffding's inequality applied after division by $r-1$ bounds the exceptional fraction by $2\exp(-k/72)$. To pull this estimate back through $f$, let $u=\lceil k/(1000\log r)\rceil$. By the preceding root bound applied to $f'$, the inputs with $v_r(f'(z))\ge u$ occupy an exponentially small fraction. If $f'$ is a nonzero constant this set is empty for large $k$. On each of the remaining classes modulo $r^u$, for different $x,y$ in the class,
$$v_r(f(x)-f(y))=v_r(x-y)+v_r(f'(x)).$$
This follows by factoring $y-x$ in Taylor's formula, since $v_r(f'(x))<u$. Each output has at most $r^u$ preimages per such class modulo $r^k$, hence at most $r^{2u}$ in all those classes. The factor $r^{2u}$ costs at most a constant times $\exp(k/500)$. Both error fractions still decrease exponentially, proving the bound.

## Divisibility with a power-saving exceptional set

Suppose at least two of the three $B_i$ are nonconstant. At every fixed $r$, outside $O(r^{(1-\gamma)k})$ residues of $t\bmod r^k$, apply the digit bound to all nonconstant parts and to $H$. The sum of the three low digit sums minus the low digit sum of $H$ is at least $(r-1)k/4$, since constant parts contribute nonnegatively. This difference equals $(r-1)$ times the number of carries in the low $k$ places, counted with multiplicity, plus the outgoing carry at the last place. That carry is at most 2. Kummer therefore gives
$$v_r(K(t))\ge k/4-2/(r-1).$$

The root bound for $g=H(H-1)$ shows that $v_r(g(t))<\lfloor k/8\rfloor$ outside another $O(r^{(1-\gamma')k})$ classes. Thus for every fixed $r$, outside a power-saving set of residues, **the entire $r$-part of $H(t)(H(t)-1)$ divides $K(t)$**. The bad residue sets can be obtained by testing low digits and divisibility modulo $r^k$ alone. Their definitions and bounds cover all sufficiently large $k$.

Apply this at the finitely many $r\mid C$, taking $k=\lfloor\log_r X\rfloor$. Each class occurs $O_r(1)$ times among $1\le t\le X$. Combined with the fixed multiplier, this proves
$$D_0(t)D_1(t)\mid K(t)$$
outside $O(X^{1-\delta})$ integer inputs, including at most that many prime inputs. All constants may depend on the row. Polynomial digit sums along integers and primes also have a classical central limit theorem. The estimates above use only low digits and elementary counting to obtain the power saving. [Bassily and Kátai, *Distribution of the values of q-additive functions on polynomial sequences*](https://doi.org/10.1007/BF01874349)

## Passing to prime powers

Return to $H=J^e$ with nonnegative coefficients and $H(0)=1$. For $P$ greater than every coefficient of $H$, a row avoiding $p$ splits the coefficients of $H$ into three nonnegative integral polynomials $B_i$, by grouping the base-$p$ digits into blocks of length $a$. At most one positive part is constant, so the digit estimate applies to every one of the finitely many shapes. The polynomial $G_0$ is $J^v$ for some $v<e$. Indeed $G_0\mid H$ and $J^e\nmid B_i$ for every part, by the degree bound, primitivity, and strict positivity of all three parts. Hence $J(P)\mid D_0(P)$.

For a fixed base $p$, omit $r=p$ from the test on $D_0$, since $p\nmid H(P)$. At any other $r\mid C$ the order of $p$ modulo $r^k$ is $\gg_{r,p} r^k$ by the valuation formula for powers. With $k=\lfloor\log_r K\rfloor$, each bad residue occurs $O_{r,p}(1)$ times among the first $K$ exponents. Sum the bounds over the finitely many shapes and primes. [Order calculation](prime_power_ray_density.md#counting-exponents)

For the uniform estimate take $k=\lfloor\log_r X\rfloor$. Among units modulo $r^k$ the power map $z\mapsto z^a$ has kernel of size at most $2(r-1)\gcd(a,r^k)$, by the cyclic unit group at odd primes and the decomposition into a sign and a cyclic group at 2. Moreover
$$\sum_{a\le Y}\gcd(a,r^k)\le\sum_{j=0}^k r^j\lfloor Y/r^j\rfloor\le (k+1)Y.$$
Count possible prime bases by all unit inputs up to $X$. Their exceptional pairs for this $r$ are $O_{H,r}(Y(k+1)X^{1-\gamma_r})$. The single base $p=r$ adds at most $Y$ and the bounded initial prime powers add $O_H(1)$. Decrease the exponent saving to absorb the logarithm and sum over the finite list. Taking the minimum of the resulting positive savings gives one $\delta_H$ for both estimates.

Finally, in any hypothetical [endpoint counterexample](all_weights_endpoint_reduction.md), at least one of the two nonendpoint rows avoids $p$, since $p\mid H(P)-1$. Its coefficient contains $J(P)$ on the good inputs. Every other positive row has a coefficient divisible by some prime of $H(P)$, since otherwise each full power in $H(P)$ divides all its parts. Every such prime divides $J(P)$ and the endpoint coefficient. This gives the asserted full three-row conclusion.

## An arithmetic progression valid at every exponent

For a possible counterexample choose the row $A=(B_i(P))$ avoiding $p$ and assign the prime powers of $H(P)$ to the two nonendpoint rows. Let $U$ be the product assigned to $A$. Its coefficient is coprime to $U$. The fixed-multiplier result for $D_0=J^{e-v}$ shows that at every $r\mid U$, $(e-v)v_r(J(P))\le v_r(C)$. In particular $U\mid C^e$. The other row is $H(P)(f_1,f_2,f_3)/U$, with the $f_i$ positive integers of sum $U$. These expressions give finitely many rational polynomial pairs, independent of $P$.

Fix an integer $t\ge1$ where $H(t)\ge3$ satisfies the conjecture. There is a modulus $L_t$ such that all sufficiently large prime powers $P\equiv t\pmod{L_t}$ satisfy it as well. First include every possible $U$ in $L_t$ so integrality at $P$ implies integrality at $t$. For each pair with $U\mid H(t)$, the two evaluated rows at $t$ are positive, so the conjecture at $H(t)$ supplies a common prime $r$ with the endpoint. Choose $k$ with $r^k>H(t)$ and include $Ur^k$ in $L_t$. At $P\equiv t\pmod{L_t}$ the parts of both rows agree with their values at $t$ modulo $r^k$. Each carries in the same low position as at $t$, and $r\mid H(P)(H(P)-1)$. This contradicts the hypothesized counterexample. Take $t=1$ and apply Dirichlet's theorem and the prime number theorem in arithmetic progressions to the unit class 1. The [verified initial interval](full_range_1b.md) supplies this anchor whenever $3\le H(1)\le10^9$.

## Reproduction

The verifier checks 640 trinomial coefficients by exact integer binomial products. In 56 cases the fixed multiplier is necessary at the tested input. It also tests 67 polynomial root bounds, 80 nonsingular fibre bounds, 600 power-map bounds and 961 low-digit carry identities. Run it with SymPy installed.

```sh
python3 check_polynomial_row_divisors.py
```
