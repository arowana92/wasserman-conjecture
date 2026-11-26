# Three distinct parts and their two successive gaps

In either nonendpoint row of a counterexample, suppose the sorted parts have successive positive differences $d,f$. Write $w=d+f$ and $M=N-1$. Then
$$w(df)^4\ge27M^4,\qquad w^9\ge6912M^4. \tag{1}$$
The first inequality bounds the smaller gap even when the two gaps are very unequal. Also, the number of possible exceptional weights through $X$ admitting such a row with $w\le L M^{4/9}$ is $O_L(X^{1/9})$ for each fixed $L$. In a dyadic interval only $O_L(1)$ gap pairs can occur. Their limiting ratios belong to a finite list of quadratic irrationals.

The [finite search and prime-power families](spread_families.md) show that every nonconstant nonendpoint row has spread greater than one million, and use (1) to prove the full conjecture at explicit infinite families of weights.

## Two integer quotients

For any two positive rows $A,C$ of weight $N$, let $D=\gcd(N,K_A,K_C)$, $J=\gcd(M,K_A,K_C)$, where $K_A,K_C$ denote their coefficients. Use the [retained powers](endpoint_pair_bound.md#retaining-partial-prime-powers) with $DUV=N$, $JST=M$, assigning $U,S$ to $A$ and $V,T$ to $C$. The entries of $A$ are multiples of $U$. The product $S=S_1S_2S_3$ is split into pairwise coprime bins by the position of residue 1. The other two actual parts vanish at each assigned power. In particular $27T^2\le(DU)^3$ since the three parts of $C/V$ sum to $DU$ and their product is divisible by $T^2$.

Sort $A$ as $(a_*,a_*+d,a_*+d+f)$ with $d,f>0$. Set $g=\gcd(d,f)$, $r=d/g$, $s=f/g$, $h=r+s$ and $v=g/U$. Each prime in $S$ divides a difference of two parts while another difference is a unit at that prime. Thus $\gcd(g,S)=1$ and
$$a=r/S_3,\qquad b=s/S_1,\qquad c=h/S_2,\qquad P=abc=rsh/S$$
are positive integers. The integers $M-3A_1=g(r+h)-1$ and $3A_3-M=g(s+h)+1$ are positive multiples of $S_2S_3$ and $S_1S_2$, respectively. Write their quotients as $\ell,j$, and put $B=b\ell$, $C_0=aj$. Their defining relations are
$$Brh=P(g(r+h)-1),\qquad C_0sh=P(g(s+h)+1).$$
Adding them and dividing by $h$ gives the first of the following identities. Multiplying the first by $s+h$ and the second by $r+h$, subtracting and dividing by $h$ gives the second.
$$Br+C_0s=3Pg,\qquad Br^2+2(B-C_0)rs-C_0s^2=-3P. \tag{2}$$
The derivation of (2) from the two displayed relations and $h=r+s\ne0$ is checked in [Lean](SpreadIdentities.lean).

As $B,C_0\ge1$, we have $3Pg\ge h$ and hence $3Pg^2\ge w$. Using $S=rsh/P$ and the companion bound gives
$$M^4\le\frac{J^4D^6(rsh)^4g^6}{729 P^4v^6},\qquad
D^6J^4\,w(df)^4\ge27 P v^6 M^4. \tag{3}$$
Thus (3) applies even when the endpoint and the two given coefficients already have a common factor. In a counterexample the [endpoint theorem](all_weights_endpoint_reduction.md) gives $D=J=1$ for its two nonendpoint rows. Now $df\le w^2/4$ proves (1).

## Another obstruction from the other two parts

In the same assignment choose any part $a_*$ of $A$ and denote its two companions by $b_*,c_*$. The bin at $a_*$ divides $\gcd(b_*,c_*)/U$ and the other two bins together divide $a_*/U$. Consequently
$$U^2S\mid a_*\gcd(b_*,c_*),\qquad
27M^2 U\le J^2D^3\big(a_*\gcd(b_*,c_*)\big)^2. \tag{4}$$
In a counterexample $U\ge3$ since the positive parts of $C/V$ sum to $U$. Thus $a_*\gcd(b_*,c_*)\ge9M$. In particular $a_*|i b_*-j c_*|\ge9M$ for every nonzero integer combination of that form, and $a_*b_*\ge9M$. If $b_*:c_*=i:j$ is in lowest terms, then $a_*>9(i+j)$, since $a_*\ge U\ge3$ and $b_*+c_*=N-a_*<M$.

## Weights near the spread bound

Assume $D=J=1$. For $1\le L\le X$ consider weights $X\le N\le2X$ with some row as above of width $w\le L M^{4/9}$. An interval of length $H\le X$ contains at most
$$O_\epsilon\big(X^\epsilon(L^6+H L^7 X^{-8/9})\big) \tag{5}$$
such weights. Across $N\le X$ the bound is $O_\epsilon(L^7 X^{1/9+\epsilon})$. With $L$ fixed, the local bound improves to $O_L(1+H X^{-8/9})$. In particular, at the larger cutoff $w\le\sqrt N$ the global count is $O_\epsilon(X^{1/2+\epsilon})$.

To prove the parameter bounds put $z=r/h$ and $K=\max(B,C_0)$. Substituting the first identity in (2) into the upper bound for $M^4$ in (3) gives
$$L^9\ge\frac{w^9}{M^4}\ge
27P v^6\,\frac{(Bz+C_0(1-z))^3}{z^4(1-z)^4}
\ge27P v^6 K^3. \tag{6}$$
Furthermore $a\mid C_0$, $b\mid B$ and $c\mid B-C_0$. The last divisibility follows by reducing the quadratic in (2) modulo $c\mid r+s$ and using $\gcd(c,r)=1$. Here $B\ne C_0$. Indeed $B-C_0=P(g(s-r)-1)/(rs)$ by the defining relations, and $g\ge U\ge3$.

For each $v$ we have $K\ll L^3/v^2$. Specifying the two positive integers $B,C_0$ and the three divisors $a,b,c$ therefore costs $O_\epsilon(L^6 X^\epsilon)$ tuples summed over $v$. The quadratic in (2) becomes
$$(Br+(B-C_0)s)^2-(B^2-BC_0+C_0^2)s^2=-3PB. \tag{7}$$
When the discriminant is square, factor the nonzero right-hand side. Otherwise the [quadratic norm count](row_gap_products.md#counting-the-possible-weights) gives $O_\epsilon(X^\epsilon)$ solutions in the required polynomial range, uniformly in these coefficients. That argument bounds the number of ideals by the divisor count of $3PB$ and the number of generators in the range by $O(\log(2X))$. Equation (2) determines $g$, then $U=g/v$ and $S=rsh/P$ are fixed.

The weights from one solution satisfy $N\equiv0\bmod U$, $N\equiv1\bmod S$. Their spacing is at least $US\ge\sqrt{27}M/\sqrt U$. Moreover,
$$M\le \frac{w^3}{4\sqrt{27}P v^3 U^{3/2}},\qquad
\sqrt U\ll \frac{L X^{1/9}}{P^{1/3}v}. \tag{8}$$
The first inequality uses $rs\le h^2/4$. Each tuple and solution thus supplies at most $1+O(HL X^{-8/9})$ weights, giving (5). For the global count use the fuller range $1\le T\le U^{3/2}/\sqrt{27}$, with $ST\equiv-1\bmod U$. It gives $O(\sqrt U)$ weights per solution even when the dyadic interval is near the start of the count. Applying the tuple estimate on each dyadic interval and summing proves the global bound. For fixed $L$ the small initial intervals are harmless. When $L$ varies, use $X^\epsilon$ throughout, since all parameters are polynomially bounded in the global $X$.

For fixed $L$ there are only finitely many parameter tuples and (6) bounds $z$ away from 0 and 1. In a dyadic weight interval, $w\asymp_L X^{4/9}$ by (1) and $g/h=(Bz+C_0(1-z))/(3P)\asymp_L1$. Thus $r,s\asymp_L X^{2/9}$. For a nonsquare discriminant $\Delta$ in (7), the ratio $(Br+(B-C_0)s)/s$ is within $O_L(s^{-2})$ of one of $\pm\sqrt\Delta$. Distinct rational numbers with denominators at most $O_L(X^{2/9})$ are spaced by $\gg_L X^{-4/9}$. An equal ratio specifies $s$ by (7). Hence only $O_L(1)$ solutions occur in the dyadic interval. This proves the stronger local bound and the asserted restriction on limiting ratios. A square discriminant permits only finitely many solutions in total. Each solution fixes $g$ and bounds $N$ by (8) or the companion inequality.

## When exactly two parts agree

A nonconstant row with exactly two equal parts has just one possible nonempty predecessor bin, at the remaining part. Put $w$ for its spread and $\sigma=1$ if the singleton is largest, $\sigma=-1$ if it is smallest. Then $S\mid w-\sigma$ and $U\mid w$. In a counterexample $w^5>27M^2$. Indeed a binary-free row has spread greater than $N/4$: the part containing the highest bit $Q$ is at least $Q$, the other two parts sum to less than $Q$, and $2Q>N$. For $w<N/4$ the companion is therefore assigned 2, making $U,S$ odd. If $w$ is even then $U\le w/2$. If $w$ is odd then $S\le(w+1)/2$. In both cases $27M^2\le U^3S^2<w^5$ since $w\ge U\ge3$. For $w\ge N/4$ the inequality holds directly at $N>10^9$, and smaller weights are [already covered](full_range_1b.md). For a distinct row (1) also implies $w^5>27M^2$.

At fixed $L$, at most $O_{L,\epsilon}(X^{2/5+\epsilon})$ weights through $X$ can have a counterexample with two equal parts and $0<w\le L M^{2/5}$. To see this put $k=w/U$, $l=(w-\sigma)/S$. The inequality $27 M^2\le U^3 S^2$ bounds $k^{3/2}l$ in terms of $L$. At fixed $k,l$ we have $S\asymp_L U\asymp_L M^{2/5}$ and $T\asymp_L U^{3/2}$. Since $lUV=(kU-\sigma)T+l$, write $T=zU+\sigma l$ with integral $z\asymp_L\sqrt U$ outside a bounded initial range.

Split $T=q_1q_2q_3$ by the companion's exceptional positions. Its normalized parts are $f_i T/q_i$, and $\prod f_i\le U^3/(27T^2)$ is bounded. Each $q_i\gg_L\sqrt U$ because $T/q_i\le U$. Their product then gives $q_i\ll_L\sqrt U$. The two equations for the companion yield
$$q_1q_2q_3-z(f_1q_2q_3+f_2q_1q_3+f_3q_1q_2)=\sigma l.$$
For each choice of $z,q_1\ll_L X^{1/5}$ let $A=q_1-zf_1$. Then
$$(Aq_2-zf_2q_1)(Aq_3-zf_3q_1)=z^2 f_2 f_3 q_1^2+\sigma l A.$$
The preceding equation gives $A>0$ for large $U$, and the right-hand side is then positive for large $U$ by the size bounds. Factoring it recovers $q_2,q_3$ in $O_{L,\epsilon}(X^\epsilon)$ choices and determines the weight. This proves the count.

With python-flint installed, run `python3 check_spread_bounds.py` for the exact coefficient checks. At weights through 40 it tests (3) on 124,306 distinct-row pairs and (4) on 452,655 choices of a part. The width bound alone forces a common factor in 59,148 of the distinct-row comparisons. The proof of the two identities in (2) also compiles under Lean 4.12 with `lean SpreadIdentities.lean`.
