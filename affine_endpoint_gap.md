# A small linear form on one row excludes every companion

Let $A,C$ be positive trinomial rows of the same weight $N$. Write $K_A,K_C$ for their coefficients and $G=\gcd(N(N-1),K_A,K_C)$. For any integer vector $v$ put $e=v\cdot A$.

**At every weight, if $|v_i|\le2$ and $1\le |e|\le2048$, then $G>1$.** The companion $C$ is arbitrary. For example, the condition $0<|A_1+2A_2-2A_3|\le2048$ certifies a common prime, including offsets on which the [gcd over the entire line](translated_line_minima.md) is one. The prime used can depend on the rows and the weight.

By the [endpoint reduction](all_weights_endpoint_reduction.md), the following separation is necessary in either nonendpoint row of a counterexample. Form the 27 sums $\sum_i d_i A_i$ with $d_i\in\{0,1,2\}$. Any two distinct values differ by at least 2,049.

More generally, for arbitrary $v$ and $e\notin\{0,v_1,v_2,v_3\}$ let $L_e=\operatorname{lcm}_i|e-v_i|$. Then
$$
                       27(N-1)^2\le G^3 |e|^3 L_e^2.                       \tag{1}
$$
For each fixed nonexceptional offset this gives a lower bound proportional to $N^{2/3}$ on the common gcd with **every** companion. With three distinct nonzero coefficients, at most $O_{v,B,\epsilon}(X^{2/9+\epsilon})$ weights $N\le X$ can have $G=1$ and $|e|\le B N^{2/9}$ outside the four offsets, for any fixed $B>0$.

## Divisibility bounds for the form

Put $D=\gcd(N,K_A,K_C)$ and $J=\gcd(N-1,K_A,K_C)$, so $G=DJ$. Use the [retained prime powers](projected_rows.md) with $N=DUV$, $N-1=JST$. The powers assigned to $A$ from $N$ have product $U$ dividing all its entries. At each power assigned to $A$ from $N-1$ the residues of its entries are $0,0,1$ in some order. Thus $U\mid e$, and each such predecessor power divides some $e-v_i$. Their product $S$ divides $L_e$. On the other side the three positive entries of $C/V$ sum to $DU$, and two vanish at each retained power of $T$. Their product is at least $T^2$, so $27(N-1)^2\le J^2 D^3 U^3 S^2$. This proves (1).

Because $U,S$ are coprime, their product even divides $\operatorname{lcm}(|e|,L_e)$. With $m$ distinct values among $0,v_1,v_2,v_3$, this yields $G^3(1+|e|)^{2m+1}\gg_v (N-1)^2$. At a failure $G=1$, $U,V\ge3$ since the normalized rows are positive triples. Thus $|e|=1,2$ are impossible even when $e$ equals one of the $v_i$.

## The finite calculation at all weights

For the height-two theorem it suffices to take primitive normals, permute their entries and change overall sign. Multiplying a primitive form by two can only increase the magnitude of a nonzero value. Constant forms with a nonzero value in the stated range have $N\le2048$. There are twelve remaining normals. Order the coordinates as $(a,b,c)$ with $b>c$. The [full range through one billion](full_range_1b.md) handles the initial weights.

At a failure, with $3\le|e|\le2048$, enumerate $U\mid |e|$, $U\ge3$, and $S\mid L_e$, coprime to $U$. The possible $T$ satisfy
$$
 27T^2\le U^3,\qquad ST\equiv-1\pmod U,\qquad N=1+ST,\qquad V=N/U .
$$
Require $V\ge3$ and $\gcd(U,V)=\gcd(S,T)=1$. Every full power $w$ in $T$ must satisfy $2w+1\le U$, since two normalized companion parts are positive multiples of $w$.

Write $A=U(x,y,z)$. For each full power $w$ of $S$, a position of residue one in $A$ can occur only at an index $i$ with $v_i\equiv e\bmod w$. At that position the normalized coordinate is $V\bmod w$, while both others vanish. The [first verifier](check_affine_endpoint.py) constructs all possible residues of $x\bmod S$, lifts them in the positive range, and solves
$$
               (a-c)x+(b-c)y=e/U-cV,\qquad z=V-x-y .
$$
It then tests carry-free addition at every prime of $US$, as a failure requires. The [independent congruence test](check_affine_endpoint_cross.py) instead solves three linear congruences, using the product of full powers $w\parallel S$ with $w\nmid e-v_i$ as a divisor of each normalized part. If two coefficients agree, it also tests the binomial for the uniquely determined cut. It checks the resulting progression using factorial-floor valuations.

There are 1,085,290 assignments above one billion before the individual $T$-power test. The largest weight among them is 144,214,167,138,805. Both row tests exclude the 524,910 assignments left by the individual power test. These are assignments for the indicated forms and offsets. The largest weight is not a verification of every intervening weight for unrestricted triples.

Run the following commands.

```sh
python3 check_affine_endpoint.py --small-check
python3 check_affine_endpoint_cross.py
```

The exact small check compares (1) against 611,074 form/coefficient pairs. The two row tests also agree with direct factorials on 8,273 positive and 41,607 negative decisions through weight 26.

## Two exact relations of height at most four

If $A,C$ each satisfy a nontrivial homogeneous relation with integer coefficients of magnitude at most four, then $G>1$ at every weight. Thus every three positive rows share a prime when **each** admits either such an exact relation or a nonzero form value in the height-two range above.

Suppose $G=1$. For a relation with all three coefficients nonzero, every full predecessor power assigned to that row divides one of the coefficients. If both relations have all coefficients nonzero, $N-1$ divides $\operatorname{lcm}(1,2,3,4)=12$. These weights are already certified. More generally, with unrestricted nonzero coefficients in both normals, $N-1$ must divide the lcm of the six coefficient magnitudes.

If a coefficient is zero, positivity gives two proportional parts $ag,bg$ with $a,b$ positive, coprime and at most four. The assigned weight factor $U$ divides $g$. Every predecessor power assigned to this row divides $g$, $a$ or $b$. Indeed, if the exceptional position is outside the pair, both parts vanish modulo the power and it divides $g$. Otherwise $g$ is a unit at that prime and the coefficient of the vanishing part is divisible by the power.

First suppose only one row has such a pair. Put $d=\gcd(12,N-1)$. The preceding alternatives on both rows give $(N-1)/d\mid g$. Since $U$ divides $N$, we have $U(N-1)/d\mid g$. The third part is at least $U\ge3$, so
$$
                         (a+b)U<d .
$$
Here $d\mid12$ and $\gcd(d,U)=1$. These conditions force $d=12$, $U=5$ and $a=b=1$. The first row has two equal parts. The normalized companion has sum five, hence is a permutation of $(1,1,3)$ or $(1,2,2)$. Both actual rows therefore have even coefficients by a binary carry, contradicting $G=1$.

If both rows have proportional pairs, write the second as $ch,fh$ in lowest terms and let $L=\operatorname{lcm}(a,b,c,f)$. Then $N\mid gh$ and $N-1\mid ghL$, so $N(N-1)\mid ghL$. Positivity gives $g(a+b),h(c+f)\le N-1$, implying $(a+b)(c+f)<L$. This is impossible with the terms at most four. Indeed $L\mid12$. If the product is below 12, either both pairs sum to three (then $L=2$), or one pair is $(1,1)$ and the other sums to at most five. In the latter case $L$ is at most the product of the remaining two terms, which is less than twice their sum in this range. This proves the assertion.

## Counting weights near the lower bound

Assume $v_1v_2v_3\ne0$, all three coefficients are distinct, $G=1$ and $X\le N<2X$ with $|e|\le B N^{2/9}$ outside the four exact offsets. Write $P(e)=\prod_i(e-v_i)$. From the bounds above,
$$
 U\asymp X^{2/9},\qquad S\asymp X^{2/3},\qquad T\asymp X^{1/3}.
$$
Indeed $U\le|e|\ll X^{2/9}$, $S\le|P(e)|\ll X^{2/3}$, and $U^3S^2\ge27(N-1)^2$. Thus the nonzero integers $e/U$ and $b_0=P(e)/S$ lie in finite sets depending only on $v,B$.

Reducing $P(e)T=b_0(N-1)$ modulo $U$ gives $A_0 T=wU+c_0$, where $A_0=|P(0)|>0$, $c_0=-\operatorname{sgn}(P(0))b_0$, and $w\asymp X^{1/9}$ for large $X$. Split $T=l_1l_2l_3$ according to the exceptional positions in the companion; the normalized parts are $h_i T/l_i$, with each $h_i$ positive. Their sum $U$ gives $l_i>T/U\gg X^{1/9}$, whence all $l_i\asymp X^{1/9}$. The bound $h_1h_2h_3\le U^3/(27T^2)$ places the $h_i$ in a fixed finite set.

Fix these bounded data, $w$, and $l_1$. Put $F=A_0 l_1-w h_1$. Substituting the parts into $A_0T=wU+c_0$ gives
$$
 (F l_2-wl_1h_2)(F l_3-wl_1h_3)
       =F c_0+w^2l_1^2 h_2h_3 .
$$
Here $F\ne0$ for large $X$ since $F l_2l_3=w l_1(h_2l_3+h_3l_2)+c_0>0$. The integer on the right is positive for large $X$ since $|F c_0|\ll X^{1/9}$ and the second term has order $X^{4/9}$. Its signed divisor pairs give $O_\epsilon(X^\epsilon)$ possibilities for $l_2,l_3$. Each determines $T,U,e,S,N$. There are $O(X^{2/9})$ pairs $w,l_1$. Summation over dyadic intervals proves the claimed count. [Lean check of the valuation-bound arithmetic and the signed factorization](AffineEndpointGap.lean)
