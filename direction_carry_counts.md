# Counting primitive trinomial pairs using the higher carry layers

Choose two primitive positive triples of sums at most $X$, and allow every positive common weight. Write $\varphi=(1+\sqrt5)/2$. Except for at most
$$
 O_\epsilon\!\left(X^{\,2+\frac12\log_2\varphi+\epsilon}\right)
$$
pairs, the two coefficients share a prime divisor of $N(N-1)$ at **every** common weight $N$. The exponent before $\epsilon$ is $2.3471209568\ldots$. The [count using just the lowest layers](counting_primitive_directions.md) gives exponent $5/2$. By the [endpoint theorem](all_weights_endpoint_reduction.md), the pair of nonendpoint rows in any counterexample to the three-row conjecture must belong to the exceptional set here.

If 2 divides one primitive sum and 3 divides the other, the exponent improves to **less than $2.1626$**. These are counts of pairs of directions. The size of the common weight is unrestricted.

## A prime-digit bound for multiples in a progression

For a fixed prime $p$ put
$$
 \rho_p=\frac1{2\sin(\pi/(4p+2))},\qquad
 b_p=\log_p\rho_p,\qquad c_p=\log_p\binom{p+1}{2}.
$$
For positive integers $a,b$, let $t=|v_p(a)-v_p(b)|$. In any progression of $H$ positive integers $w$, with arbitrary positive common difference, the number for which $aw$ and $bw$ add without a base-$p$ carry is
$$
                   O_p\big(\rho_p^t H^{b_p}\big).                 \tag{1}
$$
The exponent is sharp already for $a=1,b=p$. The delay $t$ in the prefactor involves only the two slopes, even when the step has high valuation. The same matrix occurs in the [digit sieve for three predecessor primes](inverse_digit_sieve.md).

Write $w=w_0+Pj$, $0\le j<H$, and extend the range of $j$ to all words of $h=\lceil\log_p H\rceil$ digits. Suppose $v_p(b)>v_p(a)$. Start the digit window at $v_p(aP)$. On reading digit $i$ of $j$, the current digit of $aw$ runs through all $p$ values, as does the digit of $bw$ lying $t$ places later. Each is a permutation because its slope after removal of the appropriate power of $p$ is a unit. The next $t$ digits of $bw$, starting at the current position, are already known from the prefix, including at the initial empty prefix.

Let $A_{ij}=1$ for $i+j<p$ and zero otherwise, with $0\le i,j<p$. Write $\theta=\pi/(4p+2)$ and $v_i=\cos((2i+1)\theta)$. Summing the cosines gives $Av=\rho_p v$. The entries decrease and $v_0/v_{p-1}=\rho_p$ by the last row. Give each prefix the product of $v_i/v_0$ over the digits in its queue of $t$ known digits. If its first digit is $i$, there are at most $p-i$ choices for the new input digit. Their newly revealed future digits are distinct. The total child weight relative to the prefix weight is therefore at most $\sum_{j<p-i}v_j/v_i=\rho_p$. The initial weight is at most one and every final weight is at least $\rho_p^{-t}$. This proves (1) for different valuations. With $a=1,b=p$ and $0\le w<p^h$, the digit words are exactly the paths with consecutive digits summing to at most $p-1$. Since $v$ is positive and $Av=\rho_p v$, their count grows as $\Theta_p(\rho_p^h)$. Omitting zero preserves that rate.

For equal valuations, remove their common power of $p$ from $a,b$. If $p\mid a+b$, the first occupied digits of $aw,bw$ already carry for every positive $w$. Otherwise the two current digits and their sum reduced modulo $p$ each permute $0,\ldots,p-1$ as the next input digit varies, starting the window at $v_p(P)$. Subtracting the third digit from the sum of the first two and summing over all choices shows that exactly $(p-1)/2$ choices carry. This case only arises for odd $p$. The remaining $(p+1)/2$ choices are at most $\rho_p$, since the Rayleigh quotient of the all-ones vector for $A$ is $(p+1)/2$ and its positive eigenvector gives its spectral radius. This finishes (1).

## Counting assignments on the two rows

Use the [lowest-layer parameterization](counting_primitive_directions.md). The directions $x,y$ have coprime sums $u,v$ and actual rows $mx,(N/v)y$, with $m=N/u=gv$ for an integer $g\ge1$. Assign full prime powers of $N-1$ to rows passing their lowest layers. Their respective products satisfy $ST=N-1$. At least one row has $S\ge u$ after relabeling. If both assigned products reach their sums, the earlier count already gives $O_\epsilon(X^{2+\epsilon})$. Assume $T<v$ in what follows. Work in dyadic boxes with $u,v,S,T$ of orders $U,V,D,E$. We have $DE\gg UV$, $m\asymp DE/U$ and $D\ll U^{3/2}$, the last because $S^2$ divides the product of the three entries of $x$.

**When the row with $S\ge u$ avoids $p$.** Group the prime powers of $S$ by the unique nonvanishing coordinate. With the group products ordered $d\le e\le f$, write
$$
       S=def,\qquad x=(z_1ef,z_2df,z_3de),\qquad P=de,\quad W=mf.
$$
All three groups exceed one since $x_i<u\le S$, and they are pairwise coprime. The actual parts must have residue one at their designated groups. Thus $W z_1 e\equiv1\bmod d$ and $W z_2 d\equiv1\bmod e$. These specify one class of $W$ modulo $P$ whenever feasible.

Fix $d,e,z_1,z_2$, further stratifying $d,e,f$ by dyadic size. There are $O(1+mD/P^2)$ values of $W$ to test since $W=mf$ and $f\asymp D/P$, where the scale $m$ on the right denotes its box size. Apply (1) with slopes $a=z_1 e,\ b=z_2 d$ to this one progression. The first two parts of any $p$-avoiding row must themselves add without a carry. For each resulting $W$, there are $O_\epsilon(X^\epsilon)$ factorizations as $mf$. The third congruence $m z_3 de\equiv1\bmod f$ fixes $z_3$ because $0<z_3<f$. This determines $u,N,T$. The choice $v\mid m$ costs only another divisor bound. All integers being factored have polynomial size in $X$, since $N-1=ST\ll X^3$.

In the sum over the fixed parameters the factor $\rho_p^t$ in (1) costs only a constant on average. Indeed it is at most $\rho_p^{v_p(d)+v_p(e)+v_p(z_1)+v_p(z_2)}$, and $\sum_{j\le B}\rho_p^{v_p(j)}\ll_p B$ since $\rho_p<p$. The four positive variables range up to constant multiples of $d,e,Ud/D,Ue/D$ respectively, with $d,e$ here denoting their box sizes. This bounds the number of row assignments and choices of $v$ per box by
$$
 \ll_{p,\epsilon} X^\epsilon (U/D)^2 P^2
                  (1+mD/P^2)^{b_p}.                             \tag{2}
$$

For fixed $v,T$ and an assignment of the prime powers of $T$ to positions, each companion coordinate is fixed modulo $T$. There are $O_\epsilon(X^\epsilon)$ such assignments. For an unrestricted companion there are $O((V/E)^2)$ rows since the first two coordinates determine the third. If the companion must avoid a prime $q$, the count improves to $O_q((V/E)^{c_q})$. To see this, write $y_i=r_i+Tj_i$ using the least positive representative of each class. With $0\le j_i<q^h$ and $q^h\asymp_q 1+v/T$, consider the $h$ actual part digits starting at $v_q((N/v)T)$. For each coordinate this block determines $j_i$ by invertibility after removing the power of $q$. At every position the three digits add to the corresponding digit $n_j$ of $N$, so there are exactly $\binom{n_j+2}{2}\le\binom{q+1}{2}$ possible digit triples there. Multiplying over the positions gives the claim.

Let $c=2$ for an unrestricted companion and $c=c_q$ otherwise. Multiply (2) by $(V/E)^c$. Substituting $P\ll D^{2/3}$, $m\asymp DE/U$, then $E\gg UV/D$ gives
$$
 \ll_{p,q,\epsilon} X^\epsilon U^{2-c}
       \left(D^{c-2/3}+D^{c-(2+b_p)/3}V^{b_p}\right)
 \ll_{p,q,\epsilon} X^{1+(c+b_p)/2+\epsilon}.                   \tag{3}
$$
The substitution for $E$ is valid for the second term because $c>b_p$. The final inequality uses $D\ll U^{3/2}$ and $U,V\ll X$.

**When only the row with $T<v$ is required to avoid $p$.** The [first counting argument](counting_primitive_directions.md) bounds the assignments up to the choice of the companion coordinates by $O_\epsilon(X^\epsilon U^3D^{-1}(1+E/U))$. Insert the companion bound $(V/E)^{c_p}$. If $T\ge u$ the earlier count even with exponent two is $O_\epsilon(X^{2+\epsilon})$. For $T<u$, use $E\gg UV/D$ to get
$$
     O_{p,\epsilon}(X^\epsilon U^{3-c_p}D^{c_p-1})
       = O_{p,\epsilon}(X^{(3+c_p)/2+\epsilon}).                \tag{4}
$$
All dyadic and divisor losses can be absorbed into the arbitrary $\epsilon$.

## Consequences and finite checks

For any fixed $p$, the count of lowest-layer pairs with at least one $p$-avoiding row is $O_{p,\epsilon}(X^{2+b_p/2+\epsilon})$, since $c_p-1\le b_p$. At $p=2$ an actual failure must have an avoiding row because $2\mid N(N-1)$. This proves the first assertion. If one row avoids $p$ and the other avoids $q$, the exponent from (3) and the case where both products reach their sums is
$$
       \max\big(2,\ 1+(b_p+c_q)/2,\ 1+(b_q+c_p)/2\big).
$$
It is less than $2.1626$ at $(p,q)=(2,3)$ and less than $2.1397$ when both rows avoid 2. If one primitive sum is even and the other is divisible by 3, each scaled coefficient is divisible by its own primitive sum: $N\mid (N/u)x_i K$ for each $i$ and $\gcd_i x_i=1$. Since $6\mid N$, any failure must then avoid 2 and 3 on opposite rows.

Through primitive sum 6,000, the [enumeration of every lowest-layer pair](counting_primitive_directions.md) gives the following numbers with at least one avoiding row at their unique candidate weight. The test at $p$ is made even when $p$ is outside the endpoint.

| Sum bound | Lowest-layer pairs | Avoids 2 | Avoids 3 | Avoids 5 | Avoids 7 |
|---:|---:|---:|---:|---:|---:|
| 3,000 | 75,278 | 1 | 7 | 12 | 31 |
| 6,000 | 551,062 | 5 | 11 | 41 | 109 |

On this entire 6,000 list, $\{2,3\}$ suffices as a common-prime certificate at the listed weights, as does $\{2,5\}$. Each of the five pairs with an odd coefficient has both 3 and 5 as common endpoint primes. These certificates have minimum size two. The list already contains weights 26,599 and 96,521 whose endpoint coefficients have gcd 2, and the five odd coefficients exclude that last singleton. Through sum 2,768 the prime 2 alone suffices for the entire list. The first parity survivor occurs at primitive sum bound 2,769, with $N=3,126,201$ and directions $(49,296,784),(275,319,2175)$. The first actual row has odd coefficient and the second even. Their exact gcd with $N(N-1)$ is $1,221,641,195,775=3\cdot5^2\cdot7^2\cdot11\cdot13\cdot29\cdot71\cdot1129$.

The [checker](check_direction_carry_count.py) tests the digit estimates including their sharp examples, compares the two-congruence reconstruction with direct generation of the primitive entries on a bounded box, and checks digit windows for each companion residue class on a small range. The optional exact check uses two binomial products with python-flint.

~~~sh
python3 check_direction_carry_count.py --limit 6000 --exact
~~~
