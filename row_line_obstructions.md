# Common divisors from rational lines through two rows

Write $K_A=N!/\prod_i A_i!$ for a positive trinomial row $A$ of weight $N$. For two such rows put $D=\gcd(N,K_A,K_C)$, $J=\gcd(N-1,K_A,K_C)$, and $G=DJ$. The [endpoint reduction](all_weights_endpoint_reduction.md) identifies $G=1$ as the remaining obstruction.

**Two rows sharing a part ratio have $G\ge5$.** The parts can be in different positions. More generally, if each row contains one part dividing another, the two coefficients share a prime with $N(N-1)$, even when the ratios differ.

**The rows in any counterexample must also differ substantially in every ordering.** Align the two nonendpoint rows in any way and write $C-A=g\delta$, with $g>0$ and $\delta$ primitive. Then
$$4\|\delta\|_\infty^3>27(N-1). \tag{1}$$
No pair of parts can be proportional across the rows as above. For each fixed primitive displacement $\delta$, the proof gives an explicit divisor list containing all possible exceptional weights, independent of the multiplier $g$.

## A divisor common to all rows on a line

Let $v=(v_1,v_2,v_3)$ be a primitive integer vector and consider all positive rows of weight $N$ satisfying $v\cdot A=0$. Assume this set is nonempty and $v\ne0$. Put
$$b=\gcd(v_1-v_3,v_2-v_3).$$
Then $b>0$ and $b$ divides $N$ and every $K_A$. Indeed, if $m=N/\gcd(A_1,A_2,A_3)$, the relation on the primitive row gives $b\mid m v_3$. Since $\gcd(b,v_3)=1$, $b\mid m$. Also $m\mid K_A$ because $K_A A_i/N$ is an integer at every position, by cancelling the distinguished part in the factorial formula.

If all three $v_i$ are nonzero, put $L=\operatorname{lcm}(|v_1|,|v_2|,|v_3|)$. Every row on the line has coefficient divisible by
$$b\,\frac{N-1}{\gcd(N-1,L)}. \tag{2}$$
At a prime power $p^e\parallel N-1$ let $l=v_p(K_A)$. If $l<e$, the residues of $A$ at $p^{e-l}$ are $0,0,1$ in some order by the [factorial-floor lemma](endpoint_pair_bound.md#retaining-partial-prime-powers). Reducing $v\cdot A=0$ shows that this power divides one of the $v_i$. This proves the predecessor factor, which is coprime to $b\mid N$. For instance, every arithmetic progression of positive parts at weight $3m$ has coefficient divisible by $3(3m-1)/\gcd(3m-1,2)$, using $v=(1,-2,1)$.

## Lines through a vertex and the ratio bound

Suppose two parts of $A$ are $ux,uy$ and two parts of $C$ are $vs,vt$, with $x,y,s,t>0$, $\gcd(x,y)=\gcd(s,t)=1$. Set $L_0=\operatorname{lcm}(x,y,s,t)$. Then
$$N(N-1)\mid G u v L_0,\qquad
G\ge \frac{N}{N-1}\,\frac{(x+y)(s+t)}{L_0}. \tag{3}$$
Use the [retained powers](endpoint_pair_bound.md) with $DUV=N$ and $JST=N-1$. Since the parts of $A,C$ are multiples of $U,V$, respectively, $U\mid u$ and $V\mid v$. At any retained predecessor power at least one of the two selected entries of its assigned row vanishes modulo the power. That power therefore divides $uvL_0$. The product $ST$ is coprime to $UV$, so $UVST\mid uvL_0$. Finally $u(x+y),v(s+t)\le N-1$.

For equal reduced ratios, $L_0=xy$ and $(x+y)^2/(xy)\ge4$. If each pair has an integer ratio, take $x=s=1$, giving $L_0\le yt<(1+y)(1+t)$. This proves both exclusions. For an identical pair of rows one can also bound $G\ge3(N-1)^{2/3}$. Assign all retained powers to that single row. Its three parts divided by $N/D$ sum to $D$, and their product is divisible by $((N-1)/J)^2$. Thus $27(N-1)^2\le J^2D^3\le G^3$.

## The displacement

For distinct aligned rows take $v$ to be either primitive integer normal to $A$ and $C$, obtained by dividing their cross product by the gcd of its entries. A zero coordinate of $v$ means that the other two positions have the same ratio across the rows, already covered by (3).

When all coordinates are nonzero, write $W=\max v_i-\min v_i$. The primitive direction of the intersection of $v\cdot A=0$ with the plane of fixed weight is, up to sign,
$$\delta=(v_2-v_3,v_3-v_1,v_1-v_2)/b.$$
Hence $W=b\|\delta\|_\infty$. Positivity forces both signs among the $v_i$. If the two entries of the same sign have larger absolute value $z$, the third has absolute value $W-z$. Thus
$$L\le |v_1v_2v_3|\le z^2(W-z)\le4W^3/27.$$
Formula (2) gives $N-1\le JL$ and $b\mid D$. Therefore $27(N-1)\le4JD^3\|\delta\|_\infty^3$. In a counterexample even one equal aligned coordinate would be impossible because $U,V$ would both divide that coordinate, making it divisible by $N$. Thus the $v_i$ are distinct and the product bound is strict. This proves (1).

For fixed $\delta$ in a counterexample, $b=1$. The pairwise differences of $v$ are specified by $\delta$ up to sign. Only $O(\|\delta\|_\infty)$ integer translates have both signs. Omit any with a zero coordinate. For each remaining vector (2) requires $N-1\mid L$. In particular there are $O_\epsilon(\|\delta\|_\infty^{1+\epsilon})$ candidate weights, using the divisor bound and $L\le4\|\delta\|_\infty^3/27$.

Run `python3 check_row_lines.py`. It uses exact trinomial coefficients through weight 40. It checks 687,267 reduced-ratio comparisons, including 12,648 equal ratios, and 403,925 aligned pairs. Of these, 398,179 have a normal with all coordinates nonzero. It also tests the single-row divisor on separately chosen lines.
