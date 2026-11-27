# One cofactor for an entire range of power indices

For every integer $H\ge2$ there is **one cofactor $d\ll H^{24}$** for which the same endpoint row shapes survive every finite digit depth separately at **each exponent $1\le e\le H$**, and at a set of exponents of natural density $1-O(1/H)$. The cofactor is $247s$ with $s$ having at most nine prime factors counted with multiplicity. For each indicated exponent, each prime cutoff and each depth, actual prime bases can be chosen with the full conjecture true throughout the constructed progression. The three rows have 2 as a common divisor, whose first carry on the companion can be placed above the prescribed depth. [Local tests and passage to actual weights](cofactor_index_principle.md)

There are also many cofactors without the size bound. A set of integers $s$ of natural density $\gg H^{-4}/\log^2(2H)$ works for every $e\le H$. Among primes $s$, a set of relative density $\gg H^{-4}/\log(2H)$ works for every $e\le H$ and for a set of exponents of density at least $1-5/H$ for each such prime. The prime bases producing the delayed carries can depend on $e$. For a different theorem with the same two bases for **all** exponents at fixed precision see the [independent prime progressions](all_exponent_cofactors.md).

## Keeping the old four local scalars

Use $A=(114sP,117sP,16sP+1)$ and $C=(V,V,V)$ with $P=q^e$ and $3V=1+247sP$. Write $T=\{2,3,13,19\}$ and use the rational targets $t_r$ in the [individual-index construction](all_exponent_cofactors.md). Fix
$$
 m=\prod_{r\in T}r^{1+{\bf1}_{r=2}+\lfloor\log_r H\rfloor},\qquad
 s\equiv t_r\pmod{r^{1+{\bf1}_{r=2}+\lfloor\log_r H\rfloor}}.
$$
This is one reduced class $s_0\bmod m$, with $m\le2964H^4$. Each $t_r/s$ is a local $e$th power provided $r^{1+\lfloor\log_r H\rfloor}\nmid e$ for each $r\in T$. All $e\le H$ qualify. The union of the four omitted divisibility classes has upper density at most $4/H$.

Let $z=\lceil(1000H)^{5/2}\rceil$. For primes $r\le z$ outside $T$ require $r\nmid s(1+247s)$. There are two distinct forbidden classes at each such prime. CRT gives density
$$
 {1\over m}\prod_{\substack{r\le z\\r\notin T}}(1-2/r)
 \gg {1\over H^4\log^2(2H)}
$$
by Mertens' product theorem, since the ratio of the local factor to $(1-1/r)^2$ is $1-1/(r-1)^2$. Among primes $s>z$ only one reduced residue is forbidden at each such $r$. The prime number theorem in arithmetic progressions gives relative density $\varphi(m)^{-1}\prod(1-1/(r-1))$, with the asserted lower bound by the same product comparison. These are fixed-modulus density statements for each fixed $H$.

Every further prime $r\mid1+247s$ outside $T$ is greater than $z$. It can be handled unless $r(r-1)\mid e$. If $r-1\nmid e$ choose a nonidentity $e$th power modulo $r$, which keeps $r$ out of the weight. If $r\nmid e$ choose the unit $P=-1/(247s)\equiv1\bmod r$, an $e$th power in the principal units. This makes $V=0$ locally. Summing $1/(r(r-1))$ over all integers $r>z$ bounds the additional omitted density by $1/z$. At other primes outside $T$ not dividing $s$ or $1+247s$, $P=1$ suffices. For prime $s>z$ the choice $P=1$ at $s$ also works at every exponent, by shifting the three slopes. This proves the prime-density assertion.

## A small positive representative at a divisor of a composite multiplier

Let $r\mid s$ with $r>z$, and write $s=r^v u$ with $r\nmid u$. If $r\nmid e$ and $g=\gcd(e,r-1)\le\sqrt r/1000$, there is a positive integer $w$ in the residue coset $u(\mathbb F_r^\times)^e$ with $247w<r$. The scalar $P=w/u$ is then an $e$th power locally because the residue has a root and $r\nmid e$. At this scalar $A=(114r^vw,117r^vw,16r^vw+1)$ adds without any carry in base $r$. This handles every $e\le H$ at every prime divisor of $s$, regardless of the size of $s$.

For the coset assertion take $J=\lfloor(r-1)/494\rfloor$ and $I=\{1,\ldots,J\}$. Every sum of two members is positive and at most $(r-1)/247$. For a nonprincipal multiplicative character $\chi$ extended by zero,
$$
  \sum_{x\bmod r}\left|\sum_{y\in I}\chi(x+y)\right|^2=J(r-J).
$$
The diagonal correlations are $r-1$. For distinct shifts the quotient of the two shifted arguments runs through every field element other than 1 after dropping the zero denominator, so the correlation is $-1$. Cauchy–Schwarz bounds the absolute double character sum over $I^2$ by $J\sqrt r$. Use the $g$ characters trivial on the subgroup of $e$th powers to count pairs with sum in the given coset. The count is at least $(J^2-(g-1)J\sqrt r)/g$, positive since $J\ge r/500$ for these primes. This proves the assertion without a hypothesis on the other factors of $s$.

## A polynomial bound and nine prime factors

Write $s=s_0+mt$ with $1\le t\le X=\lfloor z^8\rfloor$. For each prime through $z$ outside $T$, the polynomial $s(1+247s)$ has exactly two roots in $t$. For a squarefree product $a$ of these primes the bad intersection count differs from $X2^{\omega(a)}/a$ by at most $2^{\omega(a)}$, directly by CRT. Use the dimension-two lower sieve, with local densities $2/r$ off $T$ and zero on $T$. Its product condition follows from the preceding Mertens comparison (equivalently, summing $2\log r/r$ gives $2\log z+O(1)$). All sieve constants are uniform in $H$ since the four omitted primes are fixed. Theorem 3.3(ii) of Ford uses level $z^{7.8}$. The sum of absolute errors is $O(z^{7.8}\log z)$ since $2^{\omega(a)}\le\tau(a)$ and $\sum_{a\le y}\tau(a)\le y(1+\log y)$. This is smaller than the lower main term of order $X/\log^2 z$. [Ford, *Sieve Methods*, Theorem 3.3](https://ford126.web.illinois.edu/sieve2023.pdf)

Consequently there is a surviving $s\le m(X+1)\ll H^{24}$. All its prime divisors exceed $z$, including the fact that it is coprime to $T$ by its starting class. For all sufficiently large $H$, $s<z^{10}$ because $m\le2964H^4\ll z^{8/5}$. Hence it has at most nine prime factors counted with multiplicity. For the bounded initial range of $H$, use the choice at a sufficiently large fixed $H$ and increase the absolute constant.

For each $r\mid s$, the extra exponents missed by the coset proof lie in the multiples of $r$ or of some divisor $g$ of $r-1$ with $g>\sqrt r/1000$. Their upper density is at most
$$
                 1/r+1000\,\tau(r-1)/\sqrt r\ll r^{-9/20}.
$$
Here the elementary divisor bound $\tau(n)\ll n^{1/20}$ suffices. Since $r>z\ge H^{5/2}$ and there are at most nine such primes, their total contribution is $O(1/H)$. All the omitted conditions are finite unions of divisibility classes for this fixed $s$, so the resulting sufficient set has natural density. Combining the local choices with CRT and refining the binary prefix for the full conjecture works as in the individual-index proof.

## A single prime cofactor multiplier for the first trillion exponents

The following certificate uses $d=247s$ with
~~~text
s = 23128931932329398448536814804548764654513479119879993
L = 865967223229951653409533606573487133219537
~~~
Both numbers are prime and $1+247s=3\cdot2^{41}L$. The four valuations of $t_r/s-1$, in order at $2,3,13,19$, are $41,26,11,10$. The argument therefore applies unless the exponent is a multiple of at least one of
$$2^{40},\quad 3^{26},\quad 13^{11},\quad 19^{10},\quad L(L-1).$$
The smallest exceeds $10^{12}$. Even over all positive exponents the sum of the five reciprocal bounds is less than $2.024\cdot10^{-12}$. These exclusions specify the guarantee of the witnesses. Multiply $s$ by 247 to obtain the actual cofactor.

The [construction program](check_all_exponent_cofactors.py) searches a CRT class for $s$ and $L$. The [separate verifier](check_all_exponent_cofactors_cross.py) proves their primality recursively using orders modulo the candidate primes, checks the factorization and all four valuations, and computes the density bound with rational arithmetic. It also checks the finite-field coset test directly.
