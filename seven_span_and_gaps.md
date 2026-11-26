# Seven digit units through normalized span 100 and across large gaps

Let $p$ be prime, $a,b\ge1$, $B=p^b$, and
$$N=1+p^a\sum_{i=1}^6 B^{h_i},\qquad 0=h_1\le\cdots\le h_6=H.$$
Allow repetitions at odd primes and require all six positions to be distinct at $p=2$.

**Span theorem.** Wasserman's three-row conjecture holds at every such weight when $H\le100$. Both the outer shift $a$ and the dilation $b$ are unrestricted. There are $\binom{105}{5}=96{,}560{,}646$ normalized multisets in the odd-prime statement. The [earlier theorem](seven_progressions.md) covered normalized span 20 at odd primes.

**Gap theorem.** At any odd prime, a counterexample requires, for each successive gap $G=h_{i+1}-h_i>0$,
$$ B^{7G-5H}\le 41{,}472. \tag{1}$$
The conclusion also holds in binary when the unit at 1 belongs to a pair in the row avoiding 2. In particular it holds for all binary dilations $b\ge2$. At an odd prime two clusters with total width at most 26 certify the conjecture for **every separation**. Here a cluster's width is the difference between its largest and smallest positions.

One consequence for products, valid at all primes under the same binary distinctness assumption, is
$$N=1+p^a(1+B^i+B^j)(1+B^k),\qquad 0\le i\le j\le13,\quad k\ge0.$$
In this formula $k,a,b$ have no upper bound.

## Which allocations need checking

At an odd prime, any carry among the six specified units lowers the digit sum to at most five, covered by the [five-unit theorem](six_digit_gaps.md). Assume henceforth that the digits add freely. By the [endpoint reduction](all_weights_endpoint_reduction.md) take a nonendpoint row avoiding $p$ and assign $p$ to it. The [seven-unit allocation theorem](seven_large_bases.md) and its [binary counterpart](seven_binary_alignment.md) give two possible forms.

When 1 belongs to a pair the row is $(p^a x,p^a y,1+p^a z)$, where
$$x=B^i F,\quad F=1+B^d,\qquad y=B^j J,\quad J=1+B^e+B^f,\qquad z=B^k.$$
Here $d>0$, $0\le e\le f$, $f>0$. A pure power times a coefficient of at most three cannot supply the assigned weight factor $U\ge7$ at an odd base. At $p=2$ all six units have distinct positions and $U\ge3$. The factor $U$ is odd in both cases. Put $m=x+y+z$. The full powers in $m$ assigned to this row give factors $A,C$, according as the other zero slot alongside $1+p^a z$ is $x$ or $y$. Two unshifted zeros would force $z=0$ modulo the assigned power. Thus
$$U\mid F,J,\quad A\mid F,y+z,\quad C\mid J,x+z,\qquad
27m^2\le U^3 A^2 C^2. \tag{2}$$
The bound follows from the product of the three normalized companion parts after removing its assigned weight factor. The factors $A,C$ are odd. At odd $p$ the two even pairs in the first row carry at 2, as proved in the [seven-unit reduction](seven_digit_units.md). At $p=2$ the integer $m$ is odd.

The other possibility occurs only at $p=2,b=1$. The positive bits split into unshifted pairs at $(i,i+n),(j,j+n)$ with $i<j$, and the positive terms of the shifted triple at $(k,k+h)$. Use $W=2^a$, $F=2^n+1$ and $z=2^k(1+2^h)$. The binary alignment theorem requires either $h=n+1$ with $a+k\equiv0\bmod 2n$, or $h=n-1$ with $a+k\equiv n+1\bmod 2n$. In both cases $1+Wz\equiv0\bmod F$ and $\gcd(z,F)=1$. Every prime dividing $m=2^iF+2^jF+z$ must be avoided by the companion, and $27m^2\le U^3$ with $U\mid F$. When $b\ge2$ both $h,n$ in actual bit units would be multiples of $b$, excluding $h=n\pm1$.

## The determinant at a gap

Use a cut at the gap in (1), with low positions in $[0,L]$ and high positions in $[L+G,H]$. Write $R=H-L-G$. Any inner block sum $v$ has the form $v_0+B^{L+G}v_1$. For two sums $v,w$ of $r,t$ units, respectively, their common divisor divides $v_0w_1-v_1w_0$, of absolute value at most $rt B^{L+R}$. If both sums occupy only the same side, their common divisor away from $p$ is also at most this bound by removing the common base power. If this determinant is zero and the two sums together occupy both sides, they each occupy both, and the ratios of the high portions to the low portions agree.

Apply this to $x,y$ and to $x,y+z$ and $y,x+z$ in (2). In each case a proportional split across the cut is impossible.

- For $x,y$, the two parts of $x$ across the cut are monomials with ratio $B^d$. The high part of $y$ would be $B^d$ times its low part. The number of base-$p$ digit units in $y$ would be even, whereas $y$ has three.
- For $x,y+z$ the same ratio makes $y+z$ divisible by $1+B^d$. Since $U$ divides both $1+B^d$ and $y$, it would divide the monomial $z$, contradicting $U>1$ and $p\nmid U$.
- For $y,x+z$, both sums have three terms. If their singleton portions are on the same side, proportionality of the full sums gives $y=B^t(x+z)$ for some integer $t$. Working away from $p$ gives the same contradiction. If the singletons are on opposite sides, their product equals the product of the two two-unit portions. The former is a power of $p$. Each two-unit portion has base-$p$ digit sum two, so it cannot be a power of $p$.

Therefore $U,A,C$ are bounded by $6B^{L+R},8B^{L+R},9B^{L+R}$, respectively. With $m\ge B^H$, (2) gives (1). For $H\ge101$ and $L+R\le26$ we have $7G-5H=2H-7(L+R)\ge20$. Smaller $H$ is covered by the span theorem.

There is a complementary bound in terms of the unshifted pair itself. Every successive gap satisfies **$G<3d/2$**, and hence **$H<7d$**. Indeed $U A\mid F$ by coprimality. Thus $m/C\le F^{3/2}/\sqrt{27}$. The coprime positive integers $v,w$ with $v y=w(x+z)$ have sum $m/\gcd(y,x+z)\le m/C$. If $B^G\ge 3m/C$, the low portions on each side of this relation after multiplication are strictly smaller than $B^{L+G}$. Equality modulo that power forces separate equality on both clusters. This is exactly the impossible proportional split for $y,x+z$ above. Hence $B^G< F^{3/2}/\sqrt{3}$. For $B^d\ge3$ this is less than $B^{3d/2}$, since $(4/3)^{3/2}<\sqrt3$. For $B^d=2$, $F=3$ and $m/C\le1$ is impossible. The interval between the two positions of $x$ uses at least one of the five consecutive gaps, whose aggregate there is $d$. The other gaps sum to less than $6d$ (or vanish), proving the span bound.

For the product example with $k\le j$ use the span theorem directly. For $k>j$ the clusters both have width $j$. In the outstanding binary aligned case assume $H=k+j>100$. The gap $k-j$ exceeds $j+1$. If some paired bits lie on the same side, at least one pair still crosses since each cluster has three bits. Their lengths cannot differ by only zero or one. If all pairs cross, the sum of their lengths is $3k$, by subtracting the sum of the three low positions from the three high positions. But the lengths $n,n,n\pm1$ have sum $3n\pm1$. This proves the example also in binary.

## Exhausting the pair containing 1

The [effective-base theorem](seven_refined_cutoffs.md) proves the conjecture at all $B\ge64$. There are 26 prime powers $B<64$. The verifier uses
$$g=\operatorname{oddpart}(\gcd(F,J)),\qquad U\mid g,\qquad \gcd(g,m)=1.$$
The last assertion follows from $m\equiv z\bmod g$. A uniform upper bound for $A$ is the odd part of $F$ after removing all primes dividing $g$. For $C$ take the odd part of $J$ after removing the primes of $pg$. Denote these bounds by $A_0,C_0$. The bound for $A$ can be reduced further to $\gcd(A_0,1+B^{j-k}J)$, with the gcd interpreted after clearing powers of $B$. The analogous bound for $C$ is $\gcd(C_0,1+B^{i-k}F)$. These depend on different relative positions.

For each $d,e,f$ through 100, enumerate $r=i-k,t=j-k$ between $-100$ and 100. The unique normalized placement has $k=-\min(0,r,t)$. Maxima of the two independent gcd bounds prune whole intervals using $m\ge B^{\max(d,f)}$. The verifier discards positions above 100 and multiplicities reaching $p$, then uses the exact $m$ and the two bounds in (2). Only 654 allocations survive, with effective bases $2,3,4$.

## Testing the companion at every outer shift

Let $g$ be a known multiple of $U$ coprime to $m$, and let $S_*$ be a multiple of the predecessor powers the first row can take. In the preceding branch use the product of the two gcd bounds for $S_*$. Each full power $w=q^e\parallel m$ that does not divide $S_*$ must go to the companion. Let $T$ be their product. Write the companion as $(N/U)(v_1,v_2,v_3)$ with $\sum v_i=U$. The integers $r_i=(g/U)v_i$ sum to $g$. At each $w$ their residues are a permutation of $(0,0,g)$, since $g$ is a unit there and $N\equiv1\bmod w$. In particular $27T^2\le g^3$ and $2w+1\le g$. For each surviving assignment the verifier enumerates the possible individual $r_i$ by the Chinese remainder theorem and matches their integer sums. In the shifted-pair branch, 404 allocations reach this step. None admits a triple.

In the aligned binary branch start with $g=F$ and $S_*=1$. For $q^e=w\parallel F$ the actual two unshifted parts divided by $w$ reduce modulo $w$ to
$$-2^i(F/w)z^{-1},\qquad -2^j(F/w)z^{-1},$$
since $Wz=-1\bmod F$. If these two residues already add with a base-$q$ carry, the first row cannot avoid $q$ and $w$ may be removed from $g$. At span 100, 83,230 allocations pass the size test before this reduction and 12,123 pass afterwards. There are 8,330 CRT searches and eight surviving normalized triples. Two exact tests exclude them.

- For any two companion parts the integer $(r_i+r_j)/\gcd(r_i,r_j)$ divides their pair binomial coefficient, by multiplying it by $r_i/(r_i+r_j)$. That coefficient divides the trinomial. Three triples have a forced prime in such a divisor.
- For $q^e=w\mid T$ the two zero slots of the actual companion divided by $w$ have residues $(r_i/w)g^{-1}\bmod w$. This uses $g C_i=N r_i$ for the actual parts $C_i$. In the other five triples these digit blocks already carry in base $q$.

All primality and divisibility checks use integer arithmetic. With python-flint installed run the following from the repository root. The second verifier directly allocates small profiles, checks positive examples of the size sieve at span 100, and checks the aligned carry at a gap of 85 using factorial-layer residues.

~~~sh
python3 check_seven_span.py
python3 check_seven_span_cross.py
~~~
