# Spacing of exceptional exponents for two fixed digit blocks

Fix a prime $p$ and positive integers $B,c,h$. Put $Q=p^h>c$, $D=BQ+c$, and consider $N_a=1+p^aD$ for $a\ge1$. If $B$ has at most two base-$p$ digit units, then in **any interval of $K$ consecutive outer exponents** the number of $a$ at which the full trinomial conjecture could fail is
$$ O_{B,c}(1+K/h). $$
The same bound holds, without the digit restriction, if $\gcd(B,c)$ has no divisor that can occur as a [shifted denominator](shifted_directions.md) in a counterexample. In particular $\gcd(B,c)<364$ or a prime-power gcd suffices. The implied constant is independent of $p$, $h$, and the starting exponent.

**Consequence for three free powers.** Fix *any* positive offset $c$. Among all primes $p$ and all exponents $1\le h<n\le K$ and $1\le a\le K$, the number of parameter choices that can give an exception at
$$1+p^a(p^n+p^h+c)$$
is $O_c(K\log(2K))$. Consequently, at most $O_c(\log X\log\log X)$ possible exceptional weights in this family are at most $X$ as $X\to\infty$.

## Comparing two outer exponents

Use the [endpoint reduction](all_weights_endpoint_reduction.md). In a counterexample a nonendpoint row $A$ avoids $p$. Its parts are $(p^a x,p^a y,p^a z+1)$ with $x,y,z>0$. The normalized entries split the digits of $D$ without a base-$p$ carry. Since $Q>c$, each is a linear form in $Q$, for example $x=b_x Q+c_x$, where the three upper coefficients split the digits of $B$, and the three constants split those of $c$. Fix such an ordered split. There are at most $\binom{B+2}{2}\binom{c+2}{2}$ choices, regardless of $h$ and $a$.

Assign the full powers as in the [three-power argument](two_high_powers.md). Write $UV=N_a$ and $ST=D$, with $U,p^a S$ assigned to $A$. Then $U\ge3$, $U\mid x,y$, $S\mid\gcd(D,xy)$, and $27D^2\le U^3 S^2$. For any entry $x$ whose linear form is not proportional to $BQ+c$, put $E=Bc$. Its gcd with $D$ divides the nonzero determinant $B c_x-c b_x$, of absolute value at most $E$.

**Neither unshifted form is proportional.** Now $S\le E^2$ and $U\ge3D^{2/3}/E^{4/3}$. At two distinct exceptional exponents $a<a'$ for this split, both assigned factors $U,U'$ divide the fixed integer $x<D$. Their gcd is at least $UU'/x>9D^{1/3}/E^{8/3}$. It divides both $1+p^aD$ and $1+p^{a'}D$, so it divides $p^{a'-a}-1$. In particular
$$p^{a'-a}>9D^{1/3}/E^{8/3}.$$

**Exactly one unshifted form is proportional.** Name it $y$, so $By=b_yD$. Coprimality with $D$ gives $U\mid b_y$ and hence $U\le B$. The complementary product $T$ is at most $B^{3/2}$ by $27T^2\le U^3$. Inside $S$ take the product $F$ of powers at which the actual part $p^a x$ has the unique residue 1. For those powers $p^a x\equiv1\pmod F$. Every other power in $S$ divides $x$, including the extra $p$-power in $D$ if present. Thus $S/F\le\gcd(D,x)\le E$ and $F\ge D/(E B^{3/2})$. At two exceptional exponents, $F,F'$ both divide $D$, and $\gcd(F,F')\mid p^{a'-a}-1$. This gives
$$p^{a'-a}>D/(E^2 B^3).$$

**Both unshifted forms are proportional.** Their sum with the third is $D$, so all three are positive multiples of the primitive form $D/g$, $g=\gcd(B,c)$. Write the multipliers $m_1,m_2,m_3$ with sum $g$. Each is a positive integer by Bézout applied to the coefficients of the primitive form. Thus each of the three upper coefficients must receive a digit unit of $B$. Alternatively the actual row has shifted denominator $g/\gcd(m_1,m_2,m_3)$ dividing $g$. Either hypothesis of the theorem excludes this case. A prime-power denominator is impossible since each of its full prime powers $w$ must satisfy $3w<m$ at denominator $m$. Denominator one cannot have three positive slopes.

In the first two cases $D\ge p^h$. For all $h$ above a constant depending only on $B,c$, the gaps exceed $h/4$, since $p\ge2$. For bounded $h$ the trivial interval count $K$ suffices. Summing over the fixed splits proves the interval estimate.

For the three-power consequence the [uniform high-power theorem](two_high_powers.md) excludes every exception with $p^n>(2c)^{64}$ and $p^{n-h}>c$. The complementary range $p^n\le(2c)^{64}$ has only finitely many triples $(p,n,h)$ depending on $c$. In the remaining range $p^{n-h}\le c$ there are only finitely many $(p,n-h)$, and $B=1+p^{n-h}$ has exactly two digit units. For $p^h\le c$ again only finitely many inner pairs occur. Otherwise sum the interval bound over $h\le K$, using $\sum_{h\le K}1/h\le1+\log K$. Finally $N_a\le X$ bounds each exponent by $\log_2 X$.

Run `python3 check_two_block_spacing.py` for exact comparisons of the assignment moduli at distinct exponents. It exercises 73,493 pairs with a nontrivial shared weight factor, 8,626 with a nontrivial shared predecessor factor, and 38 fully proportional splits.
