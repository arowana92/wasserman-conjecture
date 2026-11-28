# The gcd across every integer translate of a rational line

Fix positive integers $r\le s$, $t$ with $\gcd(r,s,t)=1$. Let $G_e$ be the gcd of all positive trinomial coefficients on $rX+sY-tZ=e$, over all weights. Every integer offset $e$ has a positive row, since an integer solution exists by primitivity and a sufficiently large multiple of $(t,t,r+s)$ makes all parts positive.

* Every prime dividing any $G_e$ is at most $s+t$. This is sharp, for example when $r=s=1$, $e=0$ and $t+1$ is prime.
* The least positive integer divisible by every $G_e$ is exactly $\operatorname{lcm}_{-s\le c\le t}G_c$. The [affine carry algorithm](check_line_divisors.py) computes each of these finite gcds.
* For some effective $\delta>0$ depending on $r,s,t$, the number of $|e|\le B$ with $G_e>1$ is $O_{r,s,t}(B^{1-\delta})$.
* At $r=s=1$ there is a closed formula for $G_e$ at every $t,e$, given below.

These statements concern a single prime working across every positive row of a translate. At a fixed weight, the endpoint and a particular pair of rows can still have a common prime even when $G_e=1$.

## A bounded offset suffices for the worst valuation

At a prime $p$, choose any coefficient $a$ of the primitive normal $(r,s,-t)$ with $p\nmid a$. For $Q=p^j>|e|$ choose $0\le d<Q$ with $ad\equiv e\pmod Q$. Then $c=(e-ad)/Q$ lies between $0$ and $-a$, inclusive. A row of valuation $v_p(G_c)$ at offset $c$, multiplied by $Q$, can be combined with the low block of value $d$ in the single coordinate having coefficient $a$. This gives a positive row at $e$ with the same valuation. Thus the supremum of the minima over all integer offsets is attained in the indicated interval. Applying this at each prime gives the lcm assertion.

For $p>s+t$, take $a=r$. Each $-r\le c\le0$ has a nonnegative solution of weight at most $s+t-1$. Indeed $r,s$ generate the additive group modulo $t$. A shortest path from 0 to $c\bmod t$ using additions of $r,s$ has length at most $t-1$, so it supplies $x,y\ge0$ with $x+y\le t-1$. Set $z=(rx+sy-c)/t$. Then $0\le z\le(s(t-1)+r)/t\le s$. The resulting one-digit row has no carries. Both homogeneous edges have primitive weight at most $s+t<p$, so one can append a positive homogeneous carry-free row on high digits by the [two-edge construction](rational_row_line_divisors.md). This proves the prime cutoff.

## Almost every translate has gcd one

It suffices to prove a power saving for each of the finitely many candidate primes $p$, on $B\le |e|<2B$ of either fixed sign. Select a coefficient $a$ of that sign, and any coefficient $b$ coprime to $p$. They may refer to the same coordinate. Put $m=|a|/p^\alpha$ with $p\nmid m$.

For each class $c\bmod m$ choose a positive carry-free row $R_c$ whose offset $f_c$ is congruent to $c$. Such rows exist. Primitivity gives residues for the three parts solving the linear congruence. Lift them to positive integers, and multiply each by a separately chosen power of $p$ congruent to 1 modulo $m$ so that their digit blocks are disjoint. Choose $L\ge1$ with every template weight less than $p^L$, and put $P=p^{L+\alpha}$.

Fix $e\bmod m$. Read the digits of $e/b$ from the low end in blocks of length $L+\alpha$, meaning the compatible residues modulo powers of $p$ obtained using the inverse of $b$. After $j$ blocks let $Q=P^j$, $d=(e/b)\bmod Q$ in $[0,Q)$, and $U=(e-bd)/Q$. The class $c=U\bmod m$ is already determined by the fixed residue and $d$, since $Q$ is a unit modulo $m$. Among the $P$ possibilities for the next block there is exactly one yielding $U\equiv f_c\pmod P$. Indeed the block is $U/b\bmod P$.

When that block occurs, $n=(U-f_c)/a$ is an integer divisible by $p^L$. It is nonnegative whenever $Q\le B/C$ with $C>|b|+\max_c|f_c|+1$, since then $U-f_c$ has the sign of $e$. Use the single coordinate $b$ with value $d$ below position $\log_p Q$. Above it append $R_c$ with $n$ added to the coordinate $a$. The template fits in the low $L$ positions of this higher block and $p^L\mid n$. All parts are positive, there are no carries, and the offset is $bd+Q(f_c+an)=e$.

Taking $j=\lfloor\log_P(B/C)\rfloor$ for sufficiently large $B$, at most $(P-1)^j$ low digit prefixes of length $j$ blocks can avoid all these choices. Each surviving prefix and the fixed class modulo $m$ specify one residue modulo $mP^j$, and hence only $O_{r,s,t,p}(1)$ integers in the given interval. This gives exponent $\log_P(P-1)<1$. Summing over signs, primes and dyadic intervals yields the stated power saving (one may enlarge the exponent slightly if it is zero).

## The exact gcd when two normal coefficients equal one

On $X+Y=tZ+e$ the answer is
$$
G_e=
\begin{cases}
t & e=-1,\\
t+1 & e=0,\ t\ge2,\\
4 & e=0,\ t=1,\\
p & e>0,\quad t+1=p^k,\ e=p^j\quad(k\ge1,\ j\ge0),\\
1 & \text{otherwise}.
\end{cases}
$$
In particular, once $t+1$ is a prime power, one prime certifies infinitely many parallel lines with positive offsets, and these are exactly the positive offsets admitting any fixed prime. If two rows at the same weight each have a part $Z$ with $N=p^kZ+p^j$, the exponents can differ between the rows. Both coefficients share $p$ with the endpoint since $N$ is zero or one modulo $p$.

Put $A=X+Y=tZ+e$ and $N=A+Z$. At a prime $p$, if $A+Z$ is carry-free and $s_p(A)\ge2$, partition the digit tokens of $A$ into two nonempty groups to obtain a positive carry-free trinomial. When $p\nmid t+1$, a nonnegative carry-free seed with $X+Y=A$ suffices too. Append a positive carry-free homogeneous row made from the two edges of primitive weight $t+1$.

**Positive offsets.** If $p\nmid t+1$, set $Z=0$, $X=e$, $Y=0$ and append the homogeneous row. If $p\mid t+1$, choose $Z=p^j$ at the highest occupied position of $e$. Then $A=(e-Z)+(t+1)Z$, whose two blocks are disjoint, and adding $Z$ restores the removed digit without a carry. Thus $s_p(A)=s_p(e)-1+s_p(t+1)$. Only when both digit sums are one can $p$ be forced. In that case $N=p^kZ+p^j$ for any $Z$ and the factorial digit-sum formula gives
$$
 v_p(K)=\frac{s_p(X)+s_p(Y)-1}{p-1}+v_p\binom{p^kZ+p^j}{p^j}.
$$
The first term is a positive integer. With $Z=p^j$, split $A=p^{j+k}$ into $p^{j+k-1}$ and $(p-1)p^{j+k-1}$ to attain valuation one.

**Zero offset.** Since $N=(t+1)Z$, factorial decrement shows $t+1\mid K$. At $p\mid t+1$ and $t\ge2$, the integer $t$ has at least two base-$p$ digit tokens (its low digit is $p-1$, and for $p=2$ it is odd and exceeds one). At $Z=1$, splitting $t$ without carries leaves precisely the valuation of $t+1$ in $K=(t+1)\binom{t}{X}$. For $t=1$, $p=2$, the valuation is $s_2(X)+s_2(Y)\ge2$ since $N=2Z$. Equality holds at $(1,1,2)$. Primes outside $t+1$ are avoided by the two edges.

**Negative offsets.** Write $e=-b$, $b>0$. If $p\nmid t$, let $Z=p^J>b$ with $J$ sufficiently large. The digit of $A=tp^J-b$ at $J$ is $(t-1)\bmod p\le p-2$. Thus $A+Z$ has no carry and the long low string of digits $p-1$ allows a positive split.

If $p\mid t$ and $b\ge2$, choose $Q=p^J\le b-1<pQ$ and set $Z=Q$. Then $A=tQ-b\ge0$, and its digit at $J$ is $-\lceil b/Q\rceil\bmod p$, at most $p-2$ because $2\le\lceil b/Q\rceil\le p$. Use $A$ as the first coordinate and append a positive homogeneous row, available since $p\nmid t+1$.

Finally, for $p\mid t$ and $b=1$ the factor $\binom{(t+1)Z-1}{Z}$ equals $\tfrac{t}{t+1}\binom{(t+1)Z}{Z}$. Its valuation is at least $v_p(t)$ and is exactly that at $Z=p^J$: the summands $tZ,Z$ have disjoint occupied positions at $J$. Taking $J$ large allows a carry-free positive split of $tZ-1$. This proves the formula.

Run `python3 check_line_divisors.py 1 1 8 --offset 9` for a concrete example. The [independent checker](check_translate_spectra.py) compares the closed formula with both the carry graph and gcds of exact multinomials. It verifies 1,296 cut-line gcds with 276,517 exact coefficients, 402 prime/offset windows with 2,412 further valuations, and 1,632 digit-block witnesses for the density argument.
