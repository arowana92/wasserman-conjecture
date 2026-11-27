# Perturbing one part of a scaled prime-power triple

Let $p$ be prime, $0\le i\le j$, $M\ge1$, and $h$ a nonzero integer. Form a positive row $A$ by adding $h$ to any one part of $M(1,p^i,p^j)$. Put $s=1+p^i+p^j$, $N=Ms+h$ and $K_A=N!/\prod A_i!$. For any positive companion row $B'$ of weight $N$, write $G=\gcd(N(N-1),K_A,K_{B'})$.

**Theorem.** For every $1\le|h|\le512$, the two rows share a prime of $N(N-1)$, at all primes, exponents and scales. Thus no counterexample to the full conjecture can contain such a nonendpoint row.

There is also a uniform result at arbitrary shifts. Write $H=|h|$, $B=|h-1|$ and suppose $h\ne1$. If $G\le2$, or $h\ge2$, or $M>B$, then
$$ N\le5G^6(H+1)^7. \tag{1}$$
In particular for **each fixed positive shift $h\ge2$**, the common gcd with the endpoint grows without bound with the weight, uniformly in the companion and in the three powers. The shift $h=1$ is treated by the separate [carry-layer theorem](shifted_power_rows.md).

## Two bounds uniform in the exponents

Let $D=\gcd(N,K_A,K_{B'})$, $J=\gcd(N-1,K_A,K_{B'})$, so $G=DJ$. Assign retained prime powers to the two rows as in the [endpoint product bound](projected_rows.md). Write $U,V$ for the weight products assigned to $A,B'$, and $S,T$ for their predecessor products. Then $N=DUV$, $N-1=JST$, $DU\ge3$, each part of $A$ is divisible by $U$, and
$$ T^2\le(DU/3)^3,\qquad 27(N-1)^2\le D^3J^2U^3S^2. \tag{2}$$
Each power in $S$ divides at least one of the two unchanged entries, so $S\mid Mp^j$. We can assume $N\ge4H$ in proving (1).

If $p\mid N$, then $S\mid \gcd(M,B)$. We claim $U^2\le HN$. If any unchanged entry is $M$, then $U\mid M,h$, and $U\le H$. Otherwise the shift is at the unique minimal power. Now $U\mid M+h$ and $U\mid Hp^i$ by subtracting a multiple of $M+h$ from $Mp^i$. If $M<H$, positivity gives $U\le M+h<2H$. If $M\ge H$ and $h>0$, then $(M+h)p^i\le2Mp^i<N$. For $h<0$, positivity gives $(M+h)p^i<N$ as well. This proves the claim. From (2),
$$729N(1-1/N)^4\le G^6 H^3 B^4. \tag{3}$$
Since $N\ge3$, this implies (1). At $G=1$, $N>10^9$, the left side is greater than $728N$.

If $p\nmid N$, then $U\mid\gcd(M,h)$. Put $g=\gcd(M,B)$, $m=M/g$, $\beta=(h-1)/g$. We have $\gcd(m,\beta)=1$ and $U\mid m$ since $\gcd(U,g)=1$. As $S\mid \gcd(N-1,Mp^j)$, it divides $gP$ where $P=\gcd(ms+\beta,p^j)$. Thus
$$ ms+\beta\le K P,\qquad K=J(DU/3)^{3/2}\le (GH)^{3/2}. \tag{4}$$
Under any of the three hypotheses of (1), the integers $x=m+\beta$ and $z=m(1+p^i)+\beta$ are nonzero. For $h\ge2$ or $M>B$ they are positive. For $G\le2$, $DU\ge3$ gives $U\ge2$, so $m\ge2$ and $\gcd(m,\beta)=1$ excludes both zeros.

If $z<0$, then $P\mid z$, and $N-1\le gK|\beta|=BK$. If $z>0$, then $ms+\beta>mp^j$, so $m<K$. For $p^i>|x|$, the valuation of $x$ is less than $i$, hence $P\le |x|$. Otherwise $P\mid z$ gives $P\le(m+1)(m+|\beta|)$. This last bound consequently holds in both cases. It follows that $N-1\le B K(K+1)(K+B)$. Taking $L=(G(H+1))^{3/2}$, with $K,B,1\le L$, gives $N-1\le4B L^3$, proving (1).

For example, take any $Y\ge2$, fix $h\ge2$, and let $M$ be divisible by $h(h-1)\prod_{q\le Y,\ q\text{ prime}}q$. At every prime through $Y$ the valuations of $N(N-1)$ and $h(h-1)$ agree. Once $N>5[h(h-1)]^6(h+1)^7$, (1) forces the endpoint, the shifted row and every companion to have a common prime **greater than $Y$**. This works simultaneously for all prime-power directions meeting the weight bound.

## Exhausting the shifts through 512

It suffices to exclude $G=1$, so $D=J=1$, $U\ge3$, and $U$ is a product of full powers of $N$ at which $A$ adds without a carry. The [verified initial interval](full_range_1b.md) lets us assume $N>10^9$. We describe both finite searches, including their stopping bounds.

When $p\mid N$, (3) first excludes $H^3B^4<728\cdot10^9$. Equation (2) gives $U^3 B^2\ge27\cdot10^{18}$, which forces $U>16H$ over this shift range. Hence the shift must be at the unique minimum, $U\mid Hp^i$ and $M\ge U-H>15U/16$. Put $q=p^i$, $r=p^j$. Then
$$ 1+2q\le s\le\frac{16}{15}B\sqrt{U/27}+\frac{16B}{15U}
       <1+\frac{16}{15}B\sqrt{Hq/27}.$$
Here $16B/(15U)<1$ follows from $U>16H$. In particular $6075q<64B^2H$ and $r\le1+\sqrt{256B^2Hq/6075}$. Enumerate these prime powers and the divisors $U$ of $Hq$ with a nontrivial $p$-part. Such a part is necessary since $U>H$. Write $M+h=Uk$. The bounds
$$10^9<N=sUk-(s-1)h\le1+\sqrt{U^3 B^2/27}$$
give a finite interval for $k$. Test (2) again with $S\le\gcd(M,B)$, the unitary condition $\gcd(U,N/U)=1$ and actual carries in $A$ at all primes dividing $U$.

When $p\nmid N$, enumerate $U\mid H$, $U\ge3$, and let $T_{\max}=\lfloor\sqrt{U^3/27}\rfloor$. Necessarily $m<T\le T_{\max}$. Indeed if $T\le m$, (4) implies $z<0$, giving $m\le B$, $T\le B$ and $N-1\le BT\le B^2$. This is in the initial interval. Enumerate multiples $m$ of $U$ below $T_{\max}$ and divisors $g$ of $B$ with $\gcd(m,\beta)=1$. The proof of (4) bounds $P\le(m+1)(m+|\beta|)$ even if $z<0$. Moreover $P>1$, since $gT_{\max}<10^9$. Thus $p$ divides one of the three nonzero integers $km+\beta$ for $k=1,2,3$, by reducing $ms+\beta$ modulo $p$. Finally $r=p^j\le(T_{\max}(m+1)(m+|\beta|)+|\beta|)/m$. Enumerate $q,r$ through this bound and all shift positions, and test size, unitary divisibility and carries.

For $1\le |h|\le512$, excluding the already proved $h=1$, there are 249,897 size candidates in the weight-prime search and 356 in the other search. After the carry tests, 1,725 distinct sorted rows remain at 1,719 weights, the largest $1{,}451{,}078{,}059{,}985{,}101$. The [full endpoint solver](check_binary_gap_finite.py) tests all 146,045 size-feasible assignments at these weights. No endpoint pair survives.

The [independent companion check](check_signed_power_companions.py) instead computes the factorial valuations of each designated row at every endpoint prime. All full powers at primes dividing its coefficient are compulsory for its companion. Let $v$ be the product from $N$, $u=N/v$ and $t$ the product from $N-1$. A normalized companion of sum $u$ requires $27t^2\le u^3$ and $2w<u$ for every full power $w\mid t$ belonging to a single prime. These conditions already exclude 1,721 rows. For the other four, assigning each predecessor power to its exceptional coordinate fixes the three residues modulo $t$ by the Chinese remainder theorem. Enumerating the positive normalized parts in those residues excludes all four.

~~~sh
# Requires python-flint for certified prime factorization in the final checks.
python3 check_signed_power_rows.py --limit 512
~~~
