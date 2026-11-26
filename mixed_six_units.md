# Six digit units on both sides of the weight

**Theorem.** Let $N\ge3$. If primes $p\mid N-1$ and $q\mid N$ satisfy $s_p(N)\le6$ and $s_q(N)\le6$, then every three positive trinomial coefficients of weight $N$ have a common prime divisor. Here $s_b$ counts digit units with multiplicity. Neither the exponents nor the primes are bounded.

The [equal-gap theorem](equal_predecessor_gaps.md#two-sparse-primes) proves this except when both primes are odd, both digit sums are six, and one row $A$ must avoid both primes. It shows that in the remaining case the row is
$$
(1+p^u,\ p^i F,\ p^i F),\qquad F=1+p^d,\qquad u=td,\qquad t\text{ positive and odd}.
$$
Put $v=\min(u,i)$, $m=(N-1)/p^v$. The companion row must avoid every prime of $m$. Its three normalized parts sum to an odd divisor $U$ of $F$, and their product is divisible by $m^2$. In particular
$$
27m^2\le U^3\le(F/2)^3. \tag{1}
$$

## Comparing the two leading powers

All three parts of $A$ are even, so at an odd base each uses at least two digit units. Since $A$ avoids $q$, each uses exactly two. As $q\mid N$, avoiding it makes each part divisible by $q$. Write
$$
1+p^u=q^b(1+q^r),\qquad p^iF=q^a(1+q^f)
$$
with $a,b\ge1$, $r\ge0$ and $f>0$. Indeed $f=0$ would force $p\mid2q^a$. We have $a=v_q(F)$ and $b=a+v_q(t)$ by odd-prime lifting at $p^d\equiv-1\bmod q$. In particular $q^b\le t q^a\le t F/2$.

Set $R=p^{u-v}$, $B=p^{i-v}F$, so $m=R+2B$. Eliminating $p^v$ gives
$$
Bq^{b+r}-Rq^{a+f}=Rq^a-B(q^b-1).
$$
The right side is odd because $R,q$ are odd and $B$ is even. Its absolute value is at most $t(F/2)\max(R,B)$. It is a nonzero multiple of $q^{\min(b+r,a+f)}$. The leading power in each two-unit part is at least half that part, even if the units coincide. Hence
$$
p^v\min(R,B)\le t F\max(R,B).
$$
If $R\le B$, this gives $p^u\le tFB\le t m^2$. If $R>B$, then $p^u\le tFR^2/B\le t m^2$, since $F\le B$. By (1), for $x=p^d\ge3$ we have
$$x^t\le t(x+1)^3/216<t x^3/3.$$
For every integer $t\ge3$, $x^{t-3}\ge t/3$ by induction. Thus $t=1$ and $u=d$.

## The two units in the first part

Now $b=a$ and $F=q^a(1+q^r)$, so
$$p^i(1+q^r)=1+q^f.$$
If $r>0$, then $f=rL$ for an odd $L\ge3$. In fact the order of $q$ modulo $q^r+1$ is $2r$ since it divides $2r$ and no exponent through $r$ gives one. Also $p^i\equiv1\bmod q$, while $p^d\equiv-1\bmod q$, so $i$ is even. This would give a square
$$
(p^{i/2})^2=\frac{(-q^r)^L-1}{-q^r-1}.
$$
Ljunggren's square-quotient theorem for integer bases of absolute value greater than one has only the positive bases 3 and 7 at lengths at least three. His proof treats the negative base explicitly. [Ljunggren, Sats 1, pp. 17–18](https://msleziak.com/temp/Ljunggren.pdf) · [Terai, Lemma 2.1](https://doi.org/10.1017/S0004972713000981)

If $r=0$, we have the reciprocal equations
$$p^d+1=2q^a,\qquad q^f+1=2p^i.$$
The lemma below gives $F\le10$. But $m=R+2B>2F$, and (1) implies $F>864$. This excludes the last case.

## A reciprocal prime-power calculation

For distinct odd primes $p,q$ the only paired identities of this form, in either order, are $3^2+1=2\cdot5$ and $5+1=2\cdot3$.

By interchanging the two equations if needed, assume $i\le d$. Indeed $i>d$ would give $q^f=2p^i-1>q^a$, hence $a<f$. Modulo $p^i$ we have $2q^a\equiv1$ and $q^f\equiv-1$. Raising the first congruence to $f$ gives $p^i\mid 2^f-(-1)^a$. Thus $q^f=2p^i-1\le2^{f+1}+1$. For $f\ge3$, even $3^f$ exceeds this upper bound, by induction. For $f=1,2$ the only possibilities with $p$ odd are $q^f=5,9$.

At $q^f=5$ we have $p^i=3$ and $d=2h$ with $h$ odd, by reduction of $3^d+1$ modulo 5. Lifting gives $a=1+v_5(h)$, hence $9^h+1\le10h$, forcing $h=1$. At $q^f=9$ we have $p^i=5$, $d$ odd, and $a=1+v_3(d)$. Now $5^d+1\le6d$ forces $d=1$. Both identities follow.

## A consequence with one sparse prime

Even without a hypothesis on a weight prime, the conjecture holds if $N-1$ has at most three distinct prime factors and $s_p(N)\le6$ for some prime $p\mid N-1$. For six units the general equal-gap row has two unshifted parts $p^i F,p^j F$. Put $v=\min(u,i,j)$. Thus $m=(N-1)/p^v>2F$ and the normalized companion parts sum to $U\le F$. Each full prime power of $m$ selects one companion part with residue $U$ and the other two with residue zero. Each of the three positive parts is less than $m$, so each must be selected by some power. This needs three distinct primes in $m$ alone. The [five-unit theorem](six_digit_gaps.md) covers smaller digit sums.

The verifier makes 79,941 direct digit comparisons, checks reciprocal prime-power identities with both primes through 251 and all four exponents through 64, and computes 210,608 trinomial coefficients exactly at the 191 weights meeting the theorem through 200. Their endpoint gcd signatures pass all 4,173 pair tests. Run `python3 check_mixed_six_units.py`.
