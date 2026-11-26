# Progressions above two powers with an arbitrary leading coefficient

Let $p$ be any prime and $m\ge1$. Take an odd positive integer $i$ and an even integer $j>i$, and set $P=p^{im}$ and $Q=p^{jm}$. Thus $Q^i=P^j$. For any positive integer $t$, the full trinomial conjecture holds at
$$N=tQ+P+D,\qquad D\ge2,\quad p\mid D(D-1),\quad t^{5j+2}D^{5j}\le Q^2. \tag{1}$$
For **$t=3$** one can take any positive $i<j$, with the sharper condition
$$1\le D,\qquad p\mid D(D-1),\qquad D^{5j}\le Q^2. \tag{2}$$
The lower digits of $D$ are unrestricted. At $p=2$ every offset in the interval qualifies. For example, at every $m\ge1$ the conjecture holds at every integer in both intervals
$$[3\cdot 2^{2m}+2^m+1,\ 3\cdot 2^{2m}+2^m+\lfloor 2^{2m/5}\rfloor],$$
$$[3\cdot 2^{4m}+2^{3m}+1,\ 3\cdot 2^{4m}+2^{3m}+\lfloor 2^{2m/5}\rfloor].$$
The first interval has length of order $N^{1/5}$, the second of order $N^{1/10}$. Replacing 2 by $p$ gives the corresponding offset ranges in the two indicated residue classes.

For independently chosen exponents $k>l\ge1$, put $Q=p^k$, $P=p^l$. The conjecture holds at $3Q+P+D$ whenever
$$D\ge1,\qquad p\mid D(D-1),\qquad D^{25}\le Q,\qquad D<\min(P,Q/P). \tag{4}$$
At $p=2$ it suffices to require $D<Q/P$ in the final condition. In particular at every binary exponent $k$ and $1\le l<24k/25$, all offsets from 1 through $\lfloor 2^{k/25}\rfloor$ work.

## Reducing to the two parts missing the middle power

Put $M=N-1$. In both (1) and (2), $D\le p^{2m/5}<P$ and $P+D<Q$. The [endpoint theorem](all_weights_endpoint_reduction.md) forces a row other than the endpoint to avoid $p$, since $p\mid N(N-1)$. By Kummer's criterion each upper digit unit in that row can be split into copies of $Q$ and the single unit $P$ occurs in one part. This applies even when the leading coefficient $t$ exceeds $p$. Write the parts as
$$x=aQ+b,\qquad y=cQ+d,\qquad z=eQ+P+f,\qquad a+c+e=t,\qquad b+d+f=D.$$
All these parameters are nonnegative integers. Give this row scale $U$ from the weight and predecessor bins $S_x,S_y,S_z$ of product $S$ by the [endpoint assignment](endpoint_reduction.md). Then $U\ge3$ divides every part, and $27M^2\le U^3S^2$. If $a$ or $c$ vanishes, a positive part is at most $D$ and contradicts the [individual-part bound](endpoint_factors.md#lower-bound-on-every-part), since $D^5\le Q^2$. Hence $a,c\ge1$.

At $W=U$ take all the shifts $(u_x,u_y,u_z)$ to be zero. At $W=S_k$ take a single 1 at position $k$. Write $b'=b-u_x$, $d'=d-u_y$, $f'=f-u_z$. Then $W$ divides $aQ+b'$, $cQ+d'$ and $eQ+P+f'$. In particular it divides both
$$c b'-a d',\qquad R=(e b'-a f')^j-a^{j-i}(-b')^i. \tag{3}$$
Indeed $aP\equiv e b'-a f'$ and $aQ\equiv-b'$ modulo $W$, so the second divisibility follows from $P^j=Q^i$.

Suppose $D\ge2$ and the first integer in (3) vanishes. Then $b',d'\ge0$. A negative value for either would require both to be negative, whereas at most one shift is 1. If $b'>0$, then $R>0$ by the parities of $i,j$. If $b'=0$, also $d'=0$ and $f'\ne0$ since $D\ge2$. Again $R>0$. As $|b'|,|f'|\le D$ and $a+e\le t$, in this case $W\le2(tD)^j$.

Put $L=cb-ad= cx-ay$. The first integers in (3) for $U,S_x,S_y,S_z$ are respectively $L,L-c,L+a,L$. Each has absolute value at most $tD$, using the shifted expression. If all are nonzero, $U^3S^2\le(tD)^9$. If $L=0$, then $S_x\mid c$, $S_y\mid a$, while both $U$ and $S_z$ are at most $2(tD)^j$. Since $ac\le t^2/4$ this gives $U^3S^2\le2t^4(tD)^{5j}$. Finally if $L=c$ or $L=-a$, then $U$ and the two bins with nonzero first integer are at most $t$, and the remaining bin is at most $2(tD)^j$. The resulting bound $4t^7(tD)^{2j}$ is at most $2t^4(tD)^{5j}$ for $t,D,j\ge2$. This exhausts the cases, including the first since $5j\ge10$. Under (1) we obtain $27M^2\le 2t^4(tD)^{5j}\le2t^2 Q^2$, contradicting $M\ge tQ$.

## The coefficient three

First take $D\ge2$. After interchanging $x,y$ we have $a=1$ and $c\le2$. If $L\ne0$ then $U\mid L$ and $|L|\ge3$. The residues for the three bins give
$$S\mid (L/U)(L-c)(L+1).$$
Here $S_z$ divides $L/U$ because $\gcd(U,S_z)=1$. Neither of the other factors vanishes. Since $|L|\le2D$, the complementary bound gives
$$27M^2\le |L|^3(|L|+2)^2(|L|+1)^2\le648D^7<243D^{5j}.$$
But $D^{5j}\le Q^2$ and $M\ge3Q$.

If $L=0$, then $S_xS_y\le2$. For the two possibly large factors $U,S_z$ we have $b'=b$, $a=1$, and $e=0$ or 1. Thus $|e b'-f'|\le D$ for both. Indeed if $e=1$ then $c=1$, $b=d$, and $2b+f=D$ gives $|b-f+1|\le D$ for $D\ge2$. Under the original parity hypotheses (3) gives $U,S_z\le2D^j$. Hence $27M^2\le U^3S^2\le128D^{5j}\le128Q^2$, again impossible.

At $D=1$ we have $p\mid M$, hence $\gcd(U,p)=1$. If $a,c\ge1$, at least one of $x,y$ is either $Q$ or $2Q$ because there is just one low unit. Its divisibility by $U$ contradicts $U\ge3$. The small-part argument already covers $ac=0$.

## A gcd bound for two shifted powers

Return to $L=0$, write $Q=p^k$, $P=p^l$, and let $v=f-eb$. The two moduli $U,S_z$ divide $\gcd(Q+b,P+v)$ and $\gcd(Q+b,P+v-1)$ respectively, with $|v|,|v-1|\le D$. If $b=0$ and $D<P$, then $U$ divides the pure power $Q$ and also $N=tQ+P+D$. Thus $U\le D$ by the $p$-adic valuation of $P+D$. If $p\nmid N$ then $U=1$, and otherwise $S_z=1$ since $S_z\mid Q$ and $\gcd(M,p)=1$. This contradicts $U\ge3$ or $27M^2\le4D^3$ in all ranges considered. We can assume $1\le b\le D$.

For $s=v$ or $v-1$ put $G_s=\gcd(Q+b,P+s)$. For any coprime $1\le h\le n$ the gcd divides
$$p^{(nl-hk)_+}(-b)^h-p^{(hk-nl)_+}(-s)^n. \tag{5}$$
If the expression vanishes, then $s\ne0$ and $n$ is odd. Indeed for $n$ even, $h$ is odd and the two terms have opposite signs. Also $(Q/b)^h=(P/|s|)^n$, so there are coprime positive integers $u>w$ such that $Q/b=(u/w)^n$ and $P/|s|=(u/w)^h$, by comparison of prime valuations. Write $Q=A u^n$, $b=A w^n$, $P=C u^h$, $|s|=C w^h$. Then $A,C$ are powers of $p$ including 1, and $u$ is a positive power of $p$. The factors $u^n+w^n$ and $u^h+\operatorname{sgn}(s)w^h$ are both coprime to $p$. Their gcd divides $u+w$: modulo the gcd the unit $u/w$ has square 1 by Bézout applied to $2h,2n$, and its odd $n$-th power is $-1$. It follows that
$$G_s\le A(u+w)<2A^{1-1/n}Q^{1/n}\le2D^{1-1/n}Q^{1/n}. \tag{6}$$
This also covers a negative shift. If $s=0$ then $G_s=\gcd(b,P)\le b$ directly.

## All ratios of exponents for the coefficient three

For (2) we can suppose $i,j$ coprime, after dividing the indices by their gcd and multiplying $m$ by it. The reduced condition is at least as permissive. Write $B=p^m$, so $D\le B^{2/5}$. Only odd $j\ge3$ remains. In (5) use $(h,n)=(i,j)$, when both prefactors are 1. Each nonzero difference bounds its gcd by $2D^j$. At most one difference vanishes since $s\mapsto s^j$ is injective. Its gcd is less than $2BD$ by (6).

For $j\ge5$ this gives
$$U^3 S^2\le128\max\big(D^{5j},B^3D^{2j+3}\big)\le128B^{2j}.$$
Indeed $S_xS_y\le2$, there is at most one gcd larger than $2D^j$, and $3+(4j+6)/5\le2j$. If $j=3$ and a difference vanishes, then $b=z^3$, $|s|=z^i$ for a positive integer $z$ and $i\le2$. The adjacent shift's determinant has absolute value at most $3z^{2i}+3z^i+1\le7z^4\le2D^2$ since $b+d\ge2b$. Now $BD\ge D^2$, so $U^3S^2\le128 B^3D^7\le128B^6$. If both determinants are nonzero the original $128D^{5j}$ suffices. Each bound contradicts $27M^2\ge243B^{2j}$, proving (2).

## Independently varying exponents

Assume (4). The preceding small-part and nonzero-$L$ exclusions use only $D^{10}\le Q^2$, and the $D=1$ argument also applies. Suppose $L=0$, $b>0$. Among the six fractional parts of $zl/k$, $0\le z\le5$, two belong to the same one of five equal subintervals of $[0,1)$. Subtract and reduce the resulting pair to obtain coprime $0\le h\le n\le5$ with $n\ge1$, $|hk-nl|\le k/5$. If $h=0$ then $P\le Q^{1/5}$ and both gcds are at most $P+D<2Q^{1/5}$. For $h>0$ a nonzero difference in (5) gives
$$G_s\le2Q^{1/5}D^5\le2Q^{2/5}.$$
Vanishing with $n=1$ would give $h=1$ and $Q/P=b/|s|\le D$, excluded by hypothesis. Otherwise $n\ge3$ is odd, and (6) gives $G_s\le2D Q^{1/3}\le2Q^{2/5}$. Consequently $U^3 S^2\le4(2Q^{2/5})^5=128Q^2$, proving (4).

Finally take $p=2$ and $P\le D$ with $D^{25}\le Q$. We have $P+D\le2D<Q$. The two leading binary units of $N$ are $2Q$ and $Q$, so a positive binary carry-free row has a part at most $2D$. This violates the individual-part bound since $(2D)^5<27M^2$. This proves the wider binary interval.

## Computational checks

Run `python3 check_dense_offsets.py` to reproduce the exact row checks. The verifier assigns ordinary base-$p$ digits directly, allowing only positive sorted rows. It uses their actual gcd $G$ and the product $B$ of the three pairwise gcds with $M$. Every admissible $U$ divides $G$, and $S$ divides $B$.

At 7,673 weight and prime pairs the verifier checked all 500,023 positive unordered carry-free rows. Each failed either the scale or the complementary size inequality. Direct factorial coefficients check the digit generator at all weights through 39 and the four test primes. The signed divisibility in (5) and the carry at a shifted binary row are formalized in [Lean](OffsetPowers.lean) and compile with Lean 4.12. The [binary offset classification](binary_offset_shapes.md) treats the exponent gaps outside (4).
