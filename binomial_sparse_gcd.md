# The numerical gcd of a binomial and a sparse polynomial

Let $s\ge2$, let $p$ be a prime with $p>s$, and let
$$ F(X)=X^d\pm1,\qquad G(X)=\sum_{j=1}^s X^{a_j},\qquad
d\ge1,\quad 0=a_1\le\cdots\le a_s,\quad H=\max(d,a_s). $$
Repeated positions are allowed. Put $J=\gcd_{\mathbb Q[X]}(F,G)$, monic. For any integer $Q\ge2$,
$$
 |J(p)|\ \mid\ \gcd(F(p),G(p)),\qquad
 \frac{\gcd(F(p),G(p))}{|J(p)|}\ \le\ (2s)^{Q^{s-1}}p^{3H/Q}. \tag{1}
$$
In particular the logarithm to base $p$ of the quotient is $O_s(H^{1-1/s})$, uniformly in the prime and all positions. Every binomial–quadrinomial pair at $p\ge5$ is included. For $s=4$, $H\ge16$ the quotient is at most $p^{7H^{3/4}}$. For $s=3$ the same proof gives $p^{6H^{2/3}}$ for $H\ge8$, for either sign of the binomial. Constant or small heights can always be bounded by $p^H+1$. The [trinomial formula](binomial_trinomial_gcd.md) identifies $J$ explicitly for the plus binomial and also describes additional factors at 2 and 3.

For four terms at the remaining primes, the [small-base classification](quadrinomial_small_bases.md) gives an integer divisor using 68 offsets at 2 or 26 at 3, with a subexponential bound for the quotient.

## Separating the absolute values in a vanishing sum

Suppose $\sum_{j=1}^s p^{b_j}\xi_j=0$ for rational $b_j$ and roots of unity $\xi_j$. Under the prime hypothesis above, the terms at each fixed $b_j$ sum to zero.

First, in any minimal vanishing subsum all differences of the $b_j$ are half-integers. To see this, let $K$ be the union of the cyclotomic fields, choose a common denominator $D$, and put $u=p^{1/D}$. The extension $K(u)/K$ is Galois since all $D$-th roots of unity are in $K$. The automorphisms act on $u$ by multiplication by a cyclic group of order $m\mid D$, so $u^m\in K$ and $1,u,\ldots,u^{m-1}$ are linearly independent over $K$. A real element of $K$ is totally real: in a cyclotomic field containing it every automorphism commutes with complex conjugation. But $u^m=p^{1/(D/m)}$ has a nonreal conjugate if $D/m>2$, by Eisenstein. Group the sum by exponents of $u$ modulo $m$. Minimality puts all terms into one group, proving the assertion.

We need a length bound for cyclotomic integers. Write $\zeta=\exp(2\pi i/p)$, $\lambda=1-\zeta$ and $r=(p-1)/2$. If a nonzero sum $A$ of signed roots of unity has $A/\lambda^r$ integral, it has at least $p-1$ terms.

To prove this, choose a conductor $M$ divisible by $p$ containing all the phases. Descend from $\mathbb Q(\zeta_M)$ to $\mathbb Q(\zeta)$, removing higher prime powers and primes other than $p$. At each step select a nonzero coefficient of $A$ in the relative power basis. It is still divisible by $\lambda^r$ in the lower integer ring and can still be written with at most the original number of signed terms. In detail, when removing one level of $l^a$ for $a\ge2$, group by the exponent of a primitive $l^a$-th root modulo $l$. For removal of $l\ne p$ at level one, collect the $l$ groups $A_j$ and replace the last power by minus the sum of the other $l-1$ powers. Each resulting coefficient $A_j-A_{l-1}$ uses at most the total number of terms. These are integral relative bases: the ring of integers is generated over the lower ring by the adjoined root, whose respective monic minimal polynomial is $T^l-\zeta_{l^{a-1}}$ or $1+T+\cdots+T^{l-1}$ by the cyclotomic degree formula. Applying the same basis expansion to $A/\lambda^r$ proves divisibility of every coefficient. [Cyclotomic integer rings](https://web.math.ucsb.edu/~agboola/teaching/2021/fall/225A/neukirch.pdf)

In the last field we have a nonzero $B=\sum_j \epsilon_j\zeta^{a_j}$ with $0\le a_j<p$, $\epsilon_j=\pm1$. Suppose it uses at most $p-2$ terms. Expand at $\zeta=1-\lambda$. For each $0\le k<r$, $\sum_j\epsilon_j\binom{a_j}{k}=0\bmod p$. Indeed $p$ is a unit times $\lambda^{p-1}$, by $p=\prod_{i=1}^{p-1}(1-\zeta^i)$, and a first coefficient nonzero modulo $p$ at degree $k<r$ would give valuation exactly $k$ in units of $v(\lambda)$. At $k=0$ the positive and negative lists must have equal lengths $t\le (p-3)/2=r-1$. The binomial moments through degree $t$ give equality of the first $t$ power sums over $\mathbb F_p$. Newton's identities recover all $t$ elementary symmetric functions since $t<p$. The two multisets of exponents are equal modulo $p$, hence $B=0$, a contradiction. This proves the length bound.

In a minimal vanishing subsum with unequal heights, divide by its least power of $p$. Write $A$ for the sum of the terms at the lowest height. It is nonzero by minimality. Every other height is at least $1/2$ and is a half-integer. Therefore $A/\sqrt p$ is an algebraic integer. Since $p/\lambda^{p-1}$ is a unit, $\sqrt p/\lambda^r$ is also an algebraic unit in the field with $\sqrt p$ adjoined. Thus $A/\lambda^r$ is integral. But the lowest group has at most $s-1\le p-2$ terms, contradicting the length bound.

The threshold $s<p$ for height separation is sharp at each odd prime. For the quadratic character $\chi$ modulo $p$, the Gauss sum $\tau=\sum_{j=1}^{p-1}\chi(j)\zeta^j$ satisfies $\overline\tau=\chi(-1)\tau$ and $\tau\overline\tau=p$. For the latter identity write the two nonzero indices as $tb,b$ and sum over $b$ first. Thus $\tau$ is a root of unity times $\sqrt p$. Moving this last term to the other side gives a vanishing sum of exactly $p$ terms with two different heights, neither height group vanishing.

## Which factors can occur on an exponent line

Take a nonnegative integer direction $(v_0,v_2,\ldots,v_s)$, put $v_1=0$, and specialize to
$$ f(T)=p^d T^{v_0}\pm1,\qquad h(T)=\sum_j p^{a_j}T^{v_j}.$$
If $v_0=0$ the binomial is constant. Otherwise a common root has the form $T=\xi p^{-d/v_0}$ with $\xi$ a root of unity. Partition the terms of $h$ by the values $a_j-dv_j/v_0$. The preceding lemma forces each group to vanish at this root. The partition is independent of which root is considered.

Let $g$ be the gcd of $d$ and all within-group differences $a_j-a_i$, and $m$ the gcd of $v_0$ and the corresponding $v_j-v_i$. Then $m/g=v_0/d$. In each group choose a least $a_i$ and form $G_b(Z)=\sum_j Z^{(a_j-a_i)/g}$ there. Write $J_0=\gcd(Z^{d/g}\pm1,\{G_b\}_b)$, monic. In the coordinate $Z=p^g T^m$ the group sums are monomial multiples of $G_b(Z)$. Thus every common root is a root of $J_0(p^g T^m)$. The binomial has simple roots and nonzero constant term, so the primitive integer gcd $P(T)$ of $f,h$ divides that integer polynomial by Gauss's lemma. Also $J_0(X^g)$ divides both $F(X)$ and $G(X)$, by summing the original groups. Consequently $|P(1)|\mid |J(p)|$. This includes partitions with individually constant $G_b$.

## The remaining integer factor

Approximate the $s$ coordinates $(d,a_2,\ldots,a_s)$ as in the [weighted determinant bound](translated_triple_weights.md#a-weighted-determinant-and-an-explicit-bound). One coordinate equals $H$. Pigeonholing the other $s-1$ fractional parts gives $1\le q\le Q^{s-1}$ and integer slopes between 0 and $q$ whose products by $H/q$ differ from the coordinates by at most $H/(qQ)$. Use them for the exponent line above.

Both specialized polynomials have degree at most $q$. The coefficient-minor construction in that bound supplies an integer-polynomial combination $cP$, where $c$ is a nonzero integer. A nonzero coefficient is an integer determinant of at most $q$ binomial columns and $q$ sparse columns. Each monomial in a fixed coefficient has the same slope weight. The sum of absolute multiplicities is at most $(2s)^q$. The residual spread per binomial column is at most $H/(qQ)$ and per sparse column at most $2H/(qQ)$. Dividing out the smallest power of $p$ bounds the part of $c$ coprime to $p$ by $(2s)^q p^{3H/Q}$. If both polynomials are proportional with equal degree, use the binomial itself for the combination. Constant cases satisfy the same bound. At $T=1$ the numerical gcd divides $cP(1)$ and is coprime to $p$. This proves (1). For the displayed small values of $s$ take $Q=\lfloor H^{1/s}\rfloor\ge 2H^{1/s}/3$. At $s=3$ use the sharper residual width $2H/Q$ from the exact maximal coordinate, as in the trinomial proof.

Run the [verifier](check_binomial_sparse_gcd.py) with python-flint. It checks that the specialized primitive gcd at 1 divides the original polynomial gcd value, including lines with only some exponent differences proportional.
