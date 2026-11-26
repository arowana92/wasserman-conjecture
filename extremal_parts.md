# A growing gap above the lower bound for an individual part

Let $a$ be a part of a nonendpoint row in a counterexample, $M=N-1$, and $\rho=27M^2/a^5$. The [individual-part inequality](endpoint_factors.md) gives $\rho<1$. The following stronger bounds hold.

* At every weight, $a^5-27M^2>a^{9/2}/5$. In particular $a>3^{3/5}M^{2/5}+\sqrt a/25$. Along any unbounded sequence the lower limit of $(a^5-27M^2)/a^{9/2}$ is at least $3\sqrt3/2$.
* If $a\not\equiv3\pmod8$, then $a^5>54M^2$. This also holds for every $a\le1{,}811{,}939{,}328$. The constant improves to 81 for odd $a\not\equiv3\pmod8$, and to 108 for $a\equiv0,7\pmod8$. For every prime-power part, $243M^2\le(a-1)^2U^3<a^5$, where $U$ is the full power of its prime in $N$.
* In $X\le N\le2X$, the number of possible exceptional weights with some $a\le Y$ is $O_\epsilon(X^\epsilon Y^{13/3}/X^{4/3})$ for $1\le Y\le2X$. In particular, for any fixed $C$, a part $a\le C N^{2/5}$ allows only $O_{C,\epsilon}(X^{2/5+\epsilon})$ weights up to $X$.
* For $0<\eta\le1/2$, a part $a\le Y$ with $1-\rho\le\eta$ allows $O_\epsilon((1+\eta Y)Y^\epsilon)$ weights, independent of any upper cutoff on $N$. At the additive scale $a-3^{3/5}M^{2/5}\le C\sqrt a$, this gives $O_{C,\epsilon}(X^{1/5+\epsilon})$ possible weights through $X$.

## The equality shape

Use the [fixed-part parametrization](million_part_bound.md): $a=Ux$, $UV=N$, $ST=M$, $S=FB$, $F\mid x$, $kB=a-1$. Write the other parts of its row $A$ as $UBj_1,UBj_2$. Thus $FT=k+U(j_1+j_2)$. Partition the full powers of $T$ into three products $d_i$ according to the unique nonzero coordinate of row $C$ at each power. For positive integers $c_i$,
$$C_i=V y_i,\qquad y_i=c_iT/d_i,\qquad \sum_i y_i=U.$$
Set $m_0=c_1c_2c_3$. Then $27m_0T^2\le U^3$ by AM–GM. Substituting $M=F(a-1)T/k$ gives
$$\rho\le(1-1/a)^2\frac{F^2}{x^3k^2 m_0}<\frac{1}{x k^2 m_0}.$$

If $\rho\ge1/2$, then $x=F=k=m_0=1$. Write the bins as $d,e,f$, so that
$$a=U=de+df+ef,\quad S=U-1,\quad T=def=1+U\ell,\quad \ell=j_1+j_2\ge2.$$
The bins are pairwise coprime. Each exceeds one since $T>U$ and $T/d_i<U$. At most one is even, so $U$ is odd. Then $S$ is even and $\gcd(S,T)=1$, forcing all bins odd. Since $U+d^2=(d+e)(d+f)$, we have $U\equiv3\pmod4$. The lattice argument below strengthens this congruence to one modulo 8.

Put $H=(d+e+f)U-9T=\sum_{\rm cyc}d(e-f)^2$. The bins are distinct, so $H>0$. The integer
$$m=(H+9)/U=d+e+f-9\ell$$
is positive. Since $4\mid H$ and $U\equiv3\pmod4$, we have $m\equiv3\pmod4$. In particular $H\ge3U-9$.

## An integral lattice at equality

The matrix $L=\operatorname{diag}(d,e,f)-\ell\mathbf1\mathbf1^{\mathsf T}$ is positive definite and has determinant one. Indeed $\ell(1/d+1/e+1/f)=1-1/T<1$, so Cauchy–Schwarz gives positivity, and expansion gives $\det L=T-\ell U=1$.

Such an integral ternary Gram matrix equals $B^{\mathsf T}B$ for an integral basis matrix $B$ with determinant $\pm1$. To see this, realize the lattice in Euclidean space with covolume one. In dimensions two and three the open ball of radius $1/\sqrt2$ has volume greater than one. Two of its points differ by a nonzero lattice vector, by summing the volumes of its translates in a fundamental domain. The squared length of this vector is a positive integer less than two. Split off the unit vector orthogonally: its inner product with every lattice vector is integral, so the orthogonal complement is again integral of covolume one. Induction gives an orthonormal integral basis.

Since $\ell$ is even and the bins are odd, $B^{\mathsf T}B\equiv I\pmod2$. Each column has odd parity and cannot be the all-ones column modulo 2, since its dot product with either other odd column would be one. Thus $B$ modulo 2 is a permutation matrix. Its column norms have squared lengths $1\pmod4$, giving $d,e,f\equiv\ell+1\pmod4$.

In fact $b=B\mathbf1$ and $w=B^{-\mathsf T}\mathbf1$ are both integral vectors with all three entries odd. Their dot product is 3, and
$$\|b\|^2=m,\qquad \|w\|^2=U.$$
For the second identity, $\det(L+t\mathbf1\mathbf1^{\mathsf T})=T-(\ell-t)U=1+tU$, which also equals $1+t\|w\|^2$ on using $B$. Both $m$ and $U=a$ are sums of three odd squares and are $3\pmod8$. For instance, if $m=3$, changing coordinate signs gives $b=(1,1,1)$ and $w=(1+2i,1+2j,1-2i-2j)$, so $(U-3)/8=i^2+ij+j^2$.

## Congruences at the next two thresholds

Suppose $\rho\ge1/4$. Then $xk^2m_0<4$, so $k=1$. If $a$ is odd, $S$ is even since $B=a-1$, and all bins are odd. Here $x,U$ are odd. A product $m_0=2$ would give $U$ even, since exactly one multiplier is even. Thus, if $\rho\ge1/3$, the only possible case for an odd $a$ is the equality shape, giving $a\equiv3\pmod8$. At $\rho\ge1/4$ the extra possibilities are $x=3,m_0=1$ and $x=1,m_0=3$. In the former $U\equiv3\pmod4$ and $a=3U\equiv1\pmod4$. In the latter the sum for $U$ adds twice an odd bin product to $de+df+ef\equiv3\pmod4$, again giving $a\equiv1\pmod4$. This proves the bounds for odd residues.

If $8\mid a$, then $x\le3$ gives even $U$ and hence odd $T$. All bins are odd. Thus $m_0=1$ and $m_0=3$ are impossible since all three multipliers are odd. The remaining case is $x=1,m_0=2$, with $c=(2,1,1)^{\mathsf T}$ after reordering. Here $T-\ell U=1$, $U=2ef+df+de$. The integral positive definite matrix $L=\operatorname{diag}(2d,e,f)-\ell c c^{\mathsf T}$ has determinant two by the rank-one determinant formula. Its lattice of covolume $\sqrt2$ has a unit vector by the same ball argument, since the ball of radius $1/\sqrt2$ has volume greater than $\sqrt2$ in dimensions two and three. Splitting twice gives $L=B^{\mathsf T}\operatorname{diag}(1,1,2)B$ with $B$ unimodular. If $\ell$ were odd all diagonal entries of $L$ would be even, excluding a unit vector. Thus $\ell$ is even. Modulo 2 the first column of $B$ has first two entries zero, since it pairs trivially with a spanning basis at the form $\operatorname{diag}(1,1,0)$. Its third entry is odd. The other two columns project to the two different unit coordinate vectors in the first two positions. Hence $w=B^{-\mathsf T}c$ has first two entries odd and third entry even. Expanding $\det(L+tcc^{\mathsf T})=2(1+tU)$ yields $U=w_1^2+w_2^2+w_3^2/2\equiv2$ or $4\pmod8$, contradicting $8\mid U=a$.

## The additive improvement

Let $y_1=ef,y_2=df,y_3=de$. Expanding the symmetric polynomials gives
$$D=U^3-27T^2
=\frac U2\sum_{i<j}(y_i-y_j)^2+3\sum_{\rm cyc}y_i(y_j-y_k)^2.$$
In the first sum each difference is $d_i(d_j-d_k)$ up to sign, so $D\ge (U/2)\min(d,e,f)H$. Also $\min(d,e,f)>T/U>\sqrt{U/54}$ since $\rho\ge1/2$. Consequently $D>U^{5/2}/5$ when $U>10^6$, which holds by the [absolute bound](million_part_bound.md). We get $a^5-27M^2 \ge a^2D$. For $\rho<1/2$ the claimed inequality is immediate. Taking the difference of fifth powers gives the displayed margin for $a$.

For the lower limit it suffices to consider sequences with $\rho\to1$ since any bounded deficit on the $a^{9/2}$ scale forces that convergence. The normalized product $27y_1y_2y_3/U^3$ tends to one. Uniqueness of equality in AM–GM on the compact simplex gives $y_i\sim U/3$ and then $d_i=T/y_i\sim\sqrt{U/3}$. In the exact expression for $D$ the coefficients of every squared difference are asymptotic to $3U/2$. Moreover $\sum_{i<j}(y_i-y_j)^2\sim\sqrt{U/3}\,H$. Hence $D\sim(\sqrt3/2)U^{3/2}H$, giving the stated lower limit.

## A further restriction for prime-power parts

For $a=p^i>1$, the divisor $U\ge3$ of $a$ must equal $p^{v_p(N)}$ since $U,V$ are coprime. Now $F=1$ and $B=S$. Write $l=j_1+j_2=(T-k)/U<\sqrt{U/27}<U$. The row $A$ must avoid $p$ even after stripping the full factor $U$ from its entries. Thus addition of $(x,Bj_1,Bj_2)$ cannot overflow a block of length $v_p(U)$ in base $p$, in particular at the lowest block modulo $U$.

If $k=1$ then $B\equiv-1\pmod U$. The last two residues sum to $2U-l>U$. If $k=2$ then $U$ is odd and both $j_i$ must be odd. Indeed an even $j_i$ gives residue $U-j_i/2$, and an odd one gives $(U-j_i)/2$, since $2B\equiv-1\pmod U$. With any even $j_i$ the sum exceeds $U$ because $l<U$. Consequently $T=2+Ul$ is even, $S=(a-1)/2$ is odd and $a\equiv3\pmod4$.

In this last case $C$ avoids 2. A part of $C$ receives the leading binary unit of $N$ and exceeds $N/2$. Hence $y_1y_2y_3<U^3/32$, by maximizing the product of the other two with their sum fixed. The two normalized entries divisible by the full 2-power in $T$ must have different 2-adic valuations in the binary carry-free row. Before their multipliers $c_i$ are included, both have the same 2-adic valuation, namely $v_2(T)$. Thus $2\mid m_0$ and $2T^2\le y_1y_2y_3<U^3/32$. This gives $256M^2<(a-1)^2U^3$ when $k=2$. For $k\ge3$, $243M^2\le(a-1)^2U^3$ follows directly from the size bound on $T$.

More generally let $q$ be the least prime divisor of $a-1$. Since $k>1$ divides $a-1$, we obtain $27q^2M^2\le(a-1)^2U^3$. For $a=2^i$ with odd prime $i$, each such divisor is $1\pmod{2i}$ by the order of 2, so $q\ge2i+1$. If $2^i-1$ is prime then $S=1$ and $27M^2\le U^3$.

One version of the block obstruction applies to composite scales as well. If $k=1$ and $Q$ is any full prime power in $U$, put $b_0=U/Q$. Stripping $Q$ from the entries leaves the last two equal to $b_0Bj_i$. Modulo $Q$ these equal $-b_0j_i$. If $b_0(j_1+j_2)\le Q$ their least residues sum to at least $Q$, forcing a carry. Therefore $Q^2< U(j_1+j_2)=FT-1$. In the equality shape every full power in $a$ is consequently less than $a^{3/4}/27^{1/4}$.

## Counting a range of small parts

Group $U,T$ in dyadic boxes of sizes $u,t$ at weights $X\le N\le2X$. Use $X^{o(1)}$ for a quantity bounded by $O_\epsilon(X^\epsilon)$ for each $\epsilon>0$. There are $O(Y/u)$ choices for $x$ in a nonempty box and $X^{o(1)}$ for $F\mid x$. Since $k=F(a-1)T/M$, there are $O(Y^2t/(uX))$ choices of $k$ there. Positive integer counts justify this estimate even if the upper bound is less than one, in which case the box is empty. The product of the $c_i$ is at most $O(u^3/t^2)$, so there are at most $O(u^3/t^2)X^{o(1)}$ ordered triples. Order the bins $d_1\le d_2\le d_3$. There are at most $t^{2/3} X^{o(1)}$ choices for $d_1,d_2$ since their product is at most $T^{2/3}$.

Put $L=c_1d_2+c_2d_1$, $B'=c_3d_1d_2$. Then $U=Ld_3+B'$ and $U\mid FT-k$. Eliminating $d_3$ gives
$$U\mid Fd_1d_2B'+kL.$$
The right side is positive of polynomial size in $X$. Each divisor $U$ fixes $d_3,a,N$, with at most $X^{o(1)}$ choices. The box count is therefore $X^{o(1)}Y^3 u/(X t^{1/3})$. Since $S\le x(a-1)\ll Y^2/u$ we have $t\gg Xu/Y^2$. Substitution, $u\ll Y$ and summation over boxes prove the stated estimate. It improves the [$O_\epsilon(Y^{3/2+\epsilon})$ estimate](singleton_scales.md) when $Y=X^\gamma$ and $2/5\le\gamma<8/17$.

## The cubic discriminant and the narrowest range

With $\rho\ge1/2$ order the bins as $d,d+h,d+j$, so $0<h<j$ are even. For the monic cubic
$$P(z)=z(z+h)(z+j)-1$$
we have $P'(d)=U$ and $P(d)=T-1\equiv0\pmod U$. Its discriminant is translation invariant. In $P(d+z)$ both the linear and constant coefficients are multiples of $U$, so the cubic discriminant formula gives
$$U\mid \Delta=(hj(j-h))^2+4(h+j)^3-18hj(h+j)-27.$$
The discriminant is nonzero. Indeed, a repeated root of a rational cubic must be rational, hence an integer dividing its constant term. Neither $1$ nor $-1$ is a root here. Since $|\Delta|\ll (1+j)^6$, fixed gaps allow $O_\epsilon((1+j)^\epsilon)$ choices for $U$, and $U=3d^2+2(h+j)d+hj$ determines $d$ uniquely.

If $1-\rho\le\eta\le1/2$, then $D\le\eta U^3$. One squared difference of the $y_i$ is $e^2j^2$ with $e=d+h>\sqrt{U/54}$. Thus $j^2\le108\eta U$. There are $O(1+\eta Y)$ ordered gap pairs for $a=U\le Y$, giving the asserted count. For $a-3^{3/5}M^{2/5}\le C\sqrt a$ one has $1-\rho\le5C/\sqrt a$ by the difference of fifth powers, and $a\ll_C X^{2/5}$ when $N\le X$. Sum dyadically over $a$ to obtain $O_{C,\epsilon}(X^{1/5+\epsilon})$.

For a separate enumeration set $s=h+j$, $b=hj$. The identity $H+9=mU$ becomes
$$3m d^2+[2ms-2(h^2-hj+j^2)]d+(m-s)b-9=0.$$
We have $m\le s+1$ since $H=2(h^2-hj+j^2)d+sb\le sU$ and $U>9$. The second verifier loops over $m\equiv3\pmod4$ in this interval and uses exact integer roots. Both verifiers impose pairwise coprimality, $T\equiv1\pmod U$, $T\ge2U+1$ and $\gcd(T,U-1)=1$. Every surviving triple determines $N=(U-1)T+1$.

For $j\le4096$ both algorithms produce exactly the same 74 triples from 2,096,128 pairs of even gaps. Of these 59 also have $\gcd(U,N/U)=1$. In all 74 triples, $\binom NU$ and the forced coefficient of $C$ are even, by exact factorial floors and by disjoint-bit tests. At the largest weight $268{,}597{,}369{,}582{,}527{,}663{,}236{,}965{,}715$ the part is $72{,}098{,}004{,}035$ and these two valuations at 2 are 13 and 69.

## Bounding the multiplier without bounding the bin gaps

Independently fix even $\ell$ and order $d<e<f$. The identity $\ell(1/d+1/e+1/f)+1/T=1$ implies $\ell<d\le3\ell$. For the upper bound, $d\ge3\ell+1$ would give $1\le3\ell/d+1/d^3<1$. Let $z=d-\ell\equiv1\pmod4$. Rearranging the identity gives
$$(ze-\ell d)(zf-\ell d)=(\ell d)^2+z.$$
Both factors are positive since $(ze-\ell d)f=1+\ell de$. The smaller is at most $\ell d$ since $z\le2\ell$. Also $\gcd(T,U-1)=\gcd(T,\ell+1)$, using $\ell U\equiv-1\pmod T$. This excludes $z=1$. Enumerating the divisors of the right side therefore recovers all choices for $e,f$. This search imposes no upper limit on either gap.

For $2\le\ell\le8192$ the program factors 5,411,799 integers of the displayed form. After unitarity it recovers 19,736 triples. In each case the required binomial $\binom NU$ and the coefficient of $C$ are even. The largest reconstructed weight was $1195508910453436934561083704833025954151023757651$. Through $\ell=128$ a direct scan of the smaller factor reproduces all 71 triples. Hence $\rho\ge1/2$ requires $\ell>8192$ and $a=U>27\ell^2>1{,}811{,}939{,}328$ by $27T^2\le U^3$ and $T=1+\ell U$.

The Python programs require python-flint. The separate prime-power block test uses 76,146 rows and includes 1,378 carry-free examples. The fixed-two-bin test also checks the lattice congruences on 136 and 214 positive examples of determinant one and two.

~~~sh
mkdir -p .work
c++ -std=c++17 -O3 check_extremal_gap_roots.cpp -o .work/check_extremal_gap_roots
.work/check_extremal_gap_roots 4096 > .work/extreme_roots_4096.log
python3 check_extremal_parts.py 4096 > .work/extreme_4096.log
diff <(rg '^BINS ' .work/extreme_4096.log | sort) <(sort .work/extreme_roots_4096.log)
python3 check_extremal_parts.py --audit .work/extreme_roots_4096.log --blocks
python3 check_extremal_multipliers.py 2 8192
python3 check_extremal_multipliers.py 2 128 --direct
~~~
