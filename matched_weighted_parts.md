# Equal ratios of weighted parts and a bound on two plane coefficients

Let $A=(x,y,a)$ and $B=(z,w,b)$ be positive rows of the same weight $N$. Write $K_A,K_B$ for their trinomial coefficients, $M=N-1$, and $G=\gcd(NM,K_A,K_B)$. The parts of either row can be reordered independently.

* If $(x+2y)/a=(z+2w)/b$, then $G>1$ at every weight. The common rational value can have arbitrarily large numerator and denominator.
* More generally, suppose $(r x+s y)/a=(u z+v w)/b>0$ for integers $r,s,u,v$. Write $P=\max(r,s)$ and $Q=\max(u,v)$, exchanging the two rows so that $P\le Q$. Then $G>1$ whenever $P=1$, or $(P,Q)$ is coordinatewise at most one of $(2,11),(3,5),(4,4)$. Negative coefficients may have arbitrarily large magnitude. In the case $P=1$ there is no bound on the second pair of coefficients.
* Suppose the two rows satisfy the same nonzero homogeneous integer relation. If its primitive normal has at least two coefficients of magnitude at most **211**, then $G>1$ at every weight. The third coefficient is unrestricted.

By the [endpoint reduction](all_weights_endpoint_reduction.md), each result excludes the stated pair from occupying the two other positions in a counterexample.

## A divisor from one position in each row

Suppose $r x+s y=t a$ and $u z+v w=h b$ for arbitrary integer coefficients with $th\ne0$. Then
$$
 NM\mid G\,\operatorname{lcm}(a,b,|t|,|h|).                         \tag{1}
$$

To see this, retain the prime powers of $NM/G$ as in the [factorial-floor lemma](endpoint_pair_bound.md#retaining-partial-prime-powers). At a retained weight power the assigned row has all entries zero modulo that power. At a retained predecessor power it has residues $0,0,1$ in some order. If its selected entry is zero, the power divides $a$ or $b$. If that entry is one, the relation says that the power divides $t$ or $h$. Each retained power thus divides the lcm.

Suppose the two weighted ratios have the same positive value $\lambda=T/d$ in lowest terms. Multiplying both relations by $d$ shows $d\mid a,b$. Use $t=h=T$ in (1). Since $\operatorname{lcm}(a,b,T)$ divides $Tab/d$, this gives the integer divisibility
$$
 NM\mid G\,\lambda ab.                                                \tag{2}
$$
In particular, **any** matched value $0<\lambda\le1$ forces $G>1$ since $ab\le(N-2)^2<NM$. The smallest common positive value obtainable with any integer numerator coefficients is
$$
 \lambda_0=\frac{\operatorname{lcm}\big(b\gcd(x,y),\,a\gcd(z,w)\big)}{ab}.
$$
Indeed, after multiplication by $ab$, the two sets of attainable values are the multiples of the two lcm arguments respectively.

For the coefficient bounds in the theorem, $P,Q$ are positive and $\lambda a\le P(N-a)$, $\lambda b\le Q(N-b)$. When $P=1$, (2) already gives $G>1$ since $\lambda ab\le (N-a)b<NM$. In general
$$
 G\ \ge\ \frac{M}{N}\frac{(\lambda+P)(\lambda+Q)}{\lambda P Q}
 \ \ge\ \frac{M}{N}\left(\frac1{\sqrt P}+\frac1{\sqrt Q}\right)^2.    \tag{3}
$$
The last step applies the arithmetic-geometric mean inequality to $1/\lambda+\lambda/(PQ)$. At $(2,11)$ and $(3,5)$ the square exceeds $101/100$. This proves $G>1$ for $N\ge101$ there and for every smaller coefficient pair. The [finite certificates](RESEARCH.md) cover the smaller weights.

## The equality case at four

Assume $G=1$ and $P,Q\le4$. Assign the full weight powers to obtain coprime $UV=N$ dividing all entries of $A,B$, respectively. Apply the rational relation after dividing each row by its weight factor. It gives $a=Udi$, $b=Vdj$ for positive integers $i,j$. Thus (2) implies $M\mid dTij$. Positivity and the coefficient bounds give
$$
 4V\ge(4d+T)i,\qquad 4U\ge(4d+T)j .
$$

If both are equalities then $\lambda a=4(N-a)$ and $\lambda b=4(N-b)$, so $a=b$ is divisible by $UV=N$, which is impossible. At least one integer inequality therefore has slack at least one. Multiplying and using $N\le dTij+1$ gives
$$
 (4d-T)^2ij+(4d+T)\min(i,j)\le16 .
$$
If $T\ne4d$, this forces $ij\le16$ and $4d+T\le16$, hence $dT\le16$ and $N\le257$. These weights have finite certificates. If $T=4d$, then $V\ge2di$ and $U\ge2dj$ with at least one strict inequality, giving $N\ge dTij+2$. This exhausts the cases.

The weighted-ratio conditions also exclude sufficiently close *unequal* values. If $G=1$ then $UV=N$ divides all cross products of parts. In particular $(rx+sy)b-(uz+vw)a$ is a multiple of $N$. For coefficients covered by the theorem, any two positive values must therefore differ by at least $N/(ab)$.

## Two bounded coefficients in a common plane

A normal with a zero coefficient forces the same ratio of two positive parts in the two rows, already [excluded at every weight](row_line_obstructions.md). Otherwise permute positions and orient the normal as $(r,s,-T)$, where $T>0$ has largest magnitude. Put $C=\max(|r|,|s|)$, $P=\max(r,s)>0$, and $m=|rs|$. Assume $G=1$. At each full power of $M$ the exceptional coordinate of the assigned row has its normal coefficient divisible by that power. Hence $M\mid\operatorname{lcm}(|r|,|s|,T)$. Equation (1) gives $NM\mid\operatorname{lcm}(a,b,T)$. Consequently
$$
 k=\frac{mT}{M},\qquad h=\frac{Tab}{NM}
$$
are positive integers. The two plane equations give $a,b\le PN/(P+T)$. Thus
$$
 kh\le \frac{m P^2 N}{M^2}\frac{T^2}{(P+T)^2}
       < \frac{C^4 N}{M^2}.                                          \tag{4}
$$
In particular $C^4\ge N-1$ at every endpoint failure, since $M^2/N=N-2+1/N$. This bound applies regardless of how large the third coefficient is.

Take $C\le211$. The [verified range](full_range_1b.md) covers $N\le10^9$. At larger weights the last bound in (4) is less than two, so $k=h=1$. The two lcms above now equal their respective products. It follows that
$$
 N=1+mT,\quad ab=mN,\quad
 \gcd(|r|,|s|)=\gcd(m,T)=\gcd(a,b)=1 .
$$
Also $1=h\le TP^2N/(M(P+T)^2)$ gives $T<P^2N/M<P^2+1$, hence $T\le P^2$. Up to permuting the first two positions there are only three sign choices for $r,s$ at given absolute values. Equal magnitudes force $m=1$ and $T\le1$ here.

For a first enumeration write $a=U i$, $b=V j$ with $ij=m$. The two normalized plane equations yield
$$
 (T+s)i\equiv sV,\qquad (T+s)j\equiv sU\pmod{|r-s|}.
$$
Thus $(T+s)^2m\equiv s^2 N\pmod{|r-s|}$. At each full prime power $q\parallel T$ the assigned row must have exceptional third position since $\gcd(m,T)=1$. Its selected part is one modulo $q$, and the other is $m$ since $ab=mN\equiv m\pmod q$. Hence $a+b\equiv1+m\pmod T$.

Set $L=\lfloor PN/(P+T)\rfloor$. The possible sums lie between $\lceil2\sqrt{mN}\rceil$ and $L+\lfloor mN/L\rfloor$ when $L^2\ge mN$. The upper bound follows by maximizing $v+mN/v$ on $[\sqrt{mN},L]$. An integer square discriminant gives both cuts, and the plane recovers the first coordinate of each row as $((T+s)a-sN)/(r-s)$ or its counterpart for $b$.

The Python search visits 24,100,453 signed coefficient triples in the size interval. Of these, 821,417 pass coprimality and the plane congruence. The sum search gives 180 integral cut pairs before coprimality and positivity of the recovered rows. Two positive pairs remain. Both have an endpoint common prime.

| $N$ | $(r,s,-T)$ | Selected parts $a,b$ | $v_2(K_A),v_2(K_B)$ | $v_5(K_A),v_5(K_B)$ |
|---:|---|---|---|---|
| 1,006,297,825 | $(151,211,-31584)$ | $5384509,5954425$ | $28,25$ | $14,10$ |
| 1,694,186,605 | $(197,209,-41148)$ | $8215885,8490229$ | $28,29$ | $5,9$ |

A second enumeration assigns $a\equiv1$ or $m$ separately at every full prime power of $T$ and combines the residues by the Chinese remainder theorem. It visits possible divisors $a$ of $mN$ in those progressions and reconstructs the rows directly. It uses neither the plane congruence nor the sum discriminant. No positive row pair survives. The two weaker survivors above have residues $(29,25)$ modulo 32 and $(55,52)$ modulo 81 respectively, so neither cut is one at the indicated power. Direct divisor scans at the two listed planes recover both pairs and confirm the binary valuations by factorial floors.

## Reproduction

The Python program also checks (1) and (2) against exact trinomial coefficients at every weight through 32, allowing distinct relations. The minimal possible selected coefficient on a single row is $\gcd(x,y)/\gcd(x,y,a)$. Every other nonzero selected coefficient is its multiple. Likewise the reduced numerator of any positive multiple of $\lambda_0$ is a multiple of its reduced numerator.

~~~sh
python3 check_common_plane_median.py
c++ -std=c++17 -O3 check_common_plane_median_cross.cpp -o .work/check_common_plane_median_cross
.work/check_common_plane_median_cross
~~~
