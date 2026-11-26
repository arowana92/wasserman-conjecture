# A binomial–trinomial gcd up to a subexponential factor

Let $p$ be prime, $d\ge1$, $0\le e\le f$, $H=\max(d,f)$, $b=\gcd(d,e,f)$, and
$$G=\gcd(1+p^d,1+p^e+p^f).$$
For $p\ge5$ define $V=p^{2b}-p^b+1$ if $d/b\equiv3\pmod6$ and $\{e/b,f/b\}\bmod6=\{2,4\}$. In every other case at these primes set $V=1$. Then
$$ V\mid G,\qquad 1\le G/V\le p^{13H^{2/3}}. \tag{1} $$
In particular, for any fixed $\epsilon>0$ and $H>(13/\epsilon)^3$, a gcd above $p^{\epsilon H}$ requires that exact congruence pattern. The assigned scale $U$ in a [seven-unit row whose units-place 1 belongs to a pair](seven_span_and_gaps.md#which-allocations-need-checking) divides precisely this kind of gcd, regardless of the relative starts of its two unshifted blocks.

In bases 2 and 3 the same estimate holds with the explicit divisor below, computed from at most seven polynomial gcds of resulting degree at most two. Extra factors really can occur. For $p=2$, $e=d+1$, $f=2d$ and $d\ge2$, the trinomial is $(1+2^d)^2$. For $p=3$, take odd $k\ge3$ and $(d,e,f)=(3k,(k+1)/2,k)$. With $x=3^k$, $u=3^{(k+1)/2}$, the trinomial $x+1+u$ divides $x^3+1$ since $(x+1-u)(x+1+u)=x^2-x+1$.

## The explicit factor in small bases

Write $E=(d,e,f)$. Start the offset set $\mathcal H_p$ with the origin. At $p=2$ use $\gamma=1/2,1$ and at $p=3$ use $\gamma=1/2$, adjoining
$$ (0,\gamma,0),\qquad (0,0,\gamma),\qquad(0,-\gamma,-\gamma) $$
for each such $\gamma$. At larger primes keep just the origin. For each $t\in\mathcal H_p$, let $v$ be the primitive integer vector in the direction $E-t$, oriented with $v_1>0$. With negative Laurent powers cleared by monomials in $T$, take the primitive integer gcd $P_t(T)$ of
$$1+p^dT^{v_1},\qquad 1+p^eT^{v_2}+p^fT^{v_3}.$$
Set $V=\operatorname{lcm}_{t\in\mathcal H_p}|P_t(1)|$. Each value divides $G$ by Gauss's lemma. In the [verifier](check_binomial_trinomial_gcd.py) the same definition is evaluated after $T\mapsto T/p^n$ for an integer $n$ to keep the coefficients smaller. Under this substitution stripping monomials and passing to primitive parts changes the gcd only by a power of $p$ and a sign. Taking the part prime to $p$ at $T=p^n$ gives the same divisor.

## Locating every common root

Take any primitive integer vector $v$ with $v_1\ge0$, and use the two polynomials above at that direction. If $v_1=0$ there is no common root. Otherwise a common root $r$ satisfies $r^{v_1}=-p^{-d}$, so $r=\zeta p^\lambda$ for a root of unity $\zeta$ and $\lambda=-d/v_1$. Put $t=E+\lambda v$. The three summands at $r$ have absolute values $1,p^{t_2},p^{t_3}$, and their valuations at any place above $p$, normalized by $v_p(p)=1$, are $0,t_2,t_3$. The least two valuations must coincide. If the third is larger by $\gamma$, its squared ratio to the smaller magnitude equals $2+\omega+\omega^{-1}$ for a root of unity $\omega$. All conjugates are real and nonnegative, at most four. A positive rational power $p^{m/n}$ in lowest terms has degree $n$ and is conjugate to its multiples by all $n$-th roots of unity, by Eisenstein and Bézout. Hence $p^{2\gamma}$ is an integer at most four, giving precisely the listed candidates for $t$, including the origin when all valuations agree.

At the unique possible absolute value of $r$, the relative phases of the three terms summing to zero are fixed up to conjugation. The phase of the constant term is fixed and that of $p^d r^{v_1}$ is $-1$. Thus there are at most two triples of phases of $(r^{v_1},r^{v_2},r^{v_3})$. Primitivity makes the map injective. The binomial has simple roots, so the gcd has degree at most two.

At $t=0$ substitute $X=p^b T$. A common root of $1+X^{d/b}$ and $1+X^{e/b}+X^{f/b}$ is a root of unity whose three displayed powers are $-1$ and the two primitive cube roots. Since the three exponents are coprime together, the root itself has order six. Conversely either primitive sixth root is common precisely under the stated congruences. This gives the formula for the contribution from $t=0$ at every prime.

## Bounding the remaining scalar

For $H\ge8$ put $Q=\lfloor H^{1/3}\rfloor$. As in the [two-coordinate approximation](translated_triple_weights.md#a-weighted-determinant-and-an-explicit-bound), there is a primitive direction $v$ and $1\le q\le Q^2$ with $\|v\|_\infty\le q$, $v_1\ge0$, and a real scale $\alpha$ satisfying $\|E-\alpha v\|_\infty\le H/(qQ)$. One obtains $v$ by dividing the integer approximant to $qE/H$ by its content and adjusting $\alpha$.

Use the [coefficient-minor construction](translated_triple_weights.md#a-weighted-determinant-and-an-explicit-bound) for the two polynomials in this direction. Their degrees after clearing powers of $T$ are at most $q,2q$. It supplies a nonzero integer-polynomial combination $R=cP$, with $P$ the primitive gcd and $c$ a nonzero integer. The part of $c$ prime to $p$ is bounded by that of any nonzero coefficient of $R$. Each such coefficient is a single determinant with at most $2q$ binomial columns and $q$ trinomial columns. Its sum of absolute multiplicities is at most $2^{2q}3^q\le p^{4q}$. Every monomial in that coefficient has the same $v$-weight. Its remaining exponent contribution is the sum of at most $3q$ terms from $0$ and the coordinates of $E-\alpha v$. Hence the width of the occurring powers of $p$ is at most $6H/Q$. Removing their smallest power of $p$ yields
$$ |c|_{p'}\le p^{4q+6H/Q}. $$
If the two specialized polynomials have equal degree and are proportional, take the first polynomial itself for $R$. Its nonzero coefficients satisfy the same bound.

If the gcd is constant, $P=1$. Otherwise the root argument identifies this $v$ with the direction from one of the listed offsets, and $|P(1)|$ divides $V$. Evaluating $R$ at 1 gives $G\mid cP(1)$, while $G$ is prime to $p$. It follows that $G/V\le|c|_{p'}$. Now $Q\ge2H^{1/3}/3$ and $q\le H^{2/3}$, proving (1). For $H<8$ the bound follows from $G\le1+p^H$.

Run the following with python-flint installed. It checks the proposed divisors against exact integer gcds, including repeated terms and both exceptional primes.

~~~sh
python3 check_binomial_trinomial_gcd.py
~~~
