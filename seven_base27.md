# Seven digit units at effective base 27

The [general seven-unit theorem at effective base 31](seven_base31.md) holds already at **$B=p^b\ge27$**. Thus it covers the ordinary predecessor prime 29 and the dilation $B=3^3$. In particular, for every $N$ with $p\mid N-1$ and $s_p(N)=7$ at a prime $p\ge29$, every three positive trinomial coefficients of weight $N$ have a common prime divisor. The six positive units may repeat at odd primes and must be distinct in the binary version. The outer shift is unrestricted.

Retain $H,d,F,y,z,U,A,C,k,h_U,L$ from the previous proof, with $H\ge101$. We refine the bounds for an opposite pairing of the four terms of $y+z$ and for two trinomials with repeated positions. We always need $27 B^{2H}\le U^3A^2C^2$.

## A matching with index three

Suppose $h_U=1$, the three residues of $y$ modulo $d$ are distinct, $L\le2$, and a matching gives $X^w+1$ as its common factor with $w=d/3$. The pair of this matching lying completely in $y$ has gap $mw$ with $m$ odd. Write $y=X^a(1+X^{mw}+X^c)$ on the torus, allowing a signed $c$. The position of $z$ differs from $a+c$ by an odd multiple of $w$. For $g=\gcd(c,w)$ all four positions are in one class modulo $g$.

If $g\ge2$, use base $B^g$ in the [mean-square and signed shift estimates](seven_base31.md#three-distinct-residues-in-the-unshifted-trinomial). Writing $r=k/d\le7/15$, they give
$$U^3 A^2\le [4(19/4)^{3/2}]^{(d-k)/g}(B+1)^{2k}.$$
The same reduced norm bound $L\le2$ applies after removing the common monomial because the signed residue patterns are identical. With $C\le3^H$ the upper bound on $(U^3 A^2 C^2)^{1/H}$ is $9[4(19/4)^{3/2}]^{4/15}(B+1)^{14/15}$.

Suppose now $g=1$. At the roots of $X^w+1$, $y$ has absolute value one. At every other root of $F$, $X^w$ is one of the two primitive sixth roots $\xi$. On the $w$ roots in each such group, $X^c$ runs over a rotation of all $w$-th roots of unity since $\gcd(w,c)=1$. The product of $|1+\xi^m+X^c|$ there is at most $|1+\xi^m|^w+1$, by the polynomial identity for the product over a full set of $w$-th roots. Here $|1+\xi^m|$ is zero or $\sqrt3$. Taking the resultant over all roots gives
$$U\le (3^{w/2}+1)^2\le4\cdot3^w.$$
Also $k\le w+2$. If there is a second matching increasing $k$, the previous proof gives $k=w+2g'$ with all four exponents in one class modulo $g'=d/\operatorname{lcm}(3,n)$. Since $g'\mid w,c$, it equals one. Otherwise $k=w$. Use $A^2\le4^{d-k}(B+1)^{2k}$ at $k\le w+2$, $d=3w\le H$. This gives the upper bound
$$27\cdot4^{2/3}(B+1)^{2/3}\big[4(B+1)^4\big]^{1/101}$$
on $(U^3 A^2 C^2)^{1/H}$. Both bounds in this section are strictly below $B^2$ at 27.

## Other matching indices

With the same hypotheses on $h_U,y,L$ but no index three, first take $k\le d/4$. The mean-square bound gives $U\le(11/3)^{(d-k)/2}$. Therefore $9[4(11/3)^{3/2}]^{3/4}(B+1)^{1/2}$ suffices. For $k>d/4$ two distinct matchings are necessary since one index is at least five. Their smaller index must be five: at $7,9$ the union fraction is $15/63$, and for $n_1\ge7,n_2\ge11$ even $1/n_1+1/n_2<1/4$. Thus the indices are $5,n$ with $7\le n<20$, $n$ odd and not divisible by five. Write $g=d/(5n)$. All four terms are congruent modulo $g$ and $k/d\le1/3$. If $g\ge2$ the normalized bound is $9\cdot32^{1/3}(B+1)^{2/3}$. This uses $((3-r)/(1-r))^{3/2}4\le32$ at $r\le1/3$. If $g=1$, then $d=5n$ and $k=n+4$. For each of $n=7,9,11,13,17,19$, $9^{101}32^{d-k}(B+1)^{2k}<B^{202}$ at 27. The ratios decrease for larger $B$ and $H$.

## Two repeated trinomials and the other cases

When $C>3^H$, both trinomials $y,x+z$ repeat and their gaps $l,d$ exceed $H/2$. The six relative positions then occupy just the endpoints of two overlapping intervals of lengths $d,l$. Their union has length $H$, so the overlap length $G=d+l-H$ is a successive gap between occupied positions. The [gap theorem](seven_span_and_gaps.md) requires $7G-5H\le\log_B41472<4$. It follows that $C\le2^{d+l}+1\le2^{(12H+4)/7+1}$. Also $U\le3\cdot2^{H-1}$, $k\le1$ and $A^2\le6^H(B+1)^2$ in this case. The resulting upper bound on $(U^3A^2 C^2)^{1/H}$ is
$$48\cdot2^{24/7}\big[(27/8)2^{22/7}(B+1)^2\big]^{1/101}.$$

For $h_U\ne1$ the earlier $d=3t$ proof applies at 27 as well. For $d\ge9t$ and $t\le3$ the [normalized binomial–trinomial estimate](seven_base31.md#residue-collisions-and-a-common-polynomial-factor) gives the bound $54\cdot3^{3/2}B^{18/101}$. For $t\ge4$ use $3(d/t-2)/2\le3(d-2t)/8$ and maximize at $t/d=1/9$. The bound is $54\cdot3^{7/24}B^{2/3}$.

For distinct residues with $k=0$ or $L>2$, and for two residues with $C\le3^H$, the earlier bounds apply directly at 27. The closest comparison from the latter cases is $720(27/8)^{1/101}<729$, verified by exact powering. The previous exclusions of cancellation or a single residue continue to apply. Across all cases the verifier checks 13 rational-power bounds strictly below $B^2$ at 27 and the six short comparisons. Every quotient decreases with $B$. It also checks 1,152 integer resultants for the index-three argument.

~~~sh
python3 check_seven_base27.py
~~~
