# Products of minors from opposite rows

Let $A,B$ be positive trinomial rows of weight $N$, let $M=N-1$, and write $K_A,K_B$ for their coefficients. Fix nonnegative integer vectors $r,u$ of positive sums $s,t$. Write
$$X=\max_i|sA_i-r_i N|,\qquad Y=\max_j|tB_j-u_j N|.$$
For each position $i$ let $L_i=r_k A_l-r_l A_k$, with $(i,k,l)$ in cyclic order. Let $H_j$ be the analogous minor of $u,B$. A row is exact if it equals its target $Nr/s$ or $Nu/t$.

**Theorem.** Put $D=\gcd(N,K_A,K_B)$ and $J=\gcd(M,K_A,K_B)$. For all nine choices of $i,j$,
$$N^2 M\mid D^2 J L_i H_j\big(st A_i B_j-NM r_i u_j\big). \tag{1}$$
Consequently, if $D=J=1$ and
$$XY\big[N(st+sY+tX)+XY\big]<N^2 M, \tag{2}$$
then every pair of nonzero minors gives an exact equality
$$st A_i B_j=NM r_i u_j. \tag{3}$$
For fixed targets, (2) holds when both widths are $o(N^{2/3})$. It also allows unequal widths.

**Proof.** In the [retained prime-power assignment](endpoint_pair_bound.md), put $n=N/D$ and $m=M/J$. The weight factors imply $n\mid L_iH_j$ and $n\mid A_iB_j$. At each retained predecessor power assigned to $A$, either $A_i$ vanishes or both other positions vanish. Thus the power divides $A_i L_i$. The corresponding assertion for $B$ shows $m\mid A_iB_jL_iH_j$. Since $n,m$ are coprime, $n^2m$ divides this last product. It also divides $NM L_iH_j$. Subtraction gives (1). We have $|L_i|\le X$ by expressing it in the discrepancies $sA_l-r_l N$ and $sA_k-r_k N$, and $|H_j|\le Y$ likewise. Expanding the last factor in (1) gives the bound in (2).

## The shapes forced in two nonexact rows

Suppose (2) holds, $D=J=1$, and both rows are nonexact. Each row then has exactly one vanishing minor. At least two minors in each row are nonzero since a single nonzero coordinate of the cross product cannot be perpendicular to the positive actual row. Fixing a nonzero minor of the other row in (3) shows that all the parts indexed by nonzero minors have the same ratio to their target entries, which must be positive. Three such positions would make the row exact. This proves the assertion.

Let $i,j$ be the positions of the vanishing minors. There are nonzero rationals $d,f$ such that
$$A=(N-d)r/s+d e_i,\qquad B=(N-f)u/t+f e_j,\qquad (N-d)(N-f)=NM. \tag{4}$$
Here $e_i$ is the corresponding coordinate unit vector. The two aligned positions give the formulas away from $i,j$ and the row sums give the remaining entries. Substituting at any of the four pairs of aligned positions gives the product in (4).

**Counting weights.** Fix the two targets. In $Z\le N\le2Z$, prescribe upper bounds $X_0,Y_0$ on $X,Y$ for which (2) holds throughout. The number of possible weights with two nonexact rows and $D=J=1$ is
$$O_{r,u,\epsilon}\big((1+\min(X_0,Y_0))Z^\epsilon\big).$$
Indeed a common fixed denominator $b$ for $d,f$ is the least common multiple of all positive target entries. At any aligned position $r_k d=r_k N-sA_k$, so $|d|\le X$, and likewise $|f|\le Y$. With $v=bd$, $w=bf$, (4) gives $bN(v+w-b)=vw$. Thus $v+w-b$ divides $v(b-v)$. The divisor and the dividend are both nonzero. Otherwise (4) would give an exact row or $N=1$. For each $v$, enumerate the signed divisors, recover $w$, then $N$. The divisor bound applies since the actual discrepancies are $O_{r,u}(Z)$. Interchange the two rows for the minimum.

**Example.** For two targets $(1,2,4)/7$ and both widths $o(N^{2/3})$, the two coefficients always share an endpoint prime for sufficiently large $N$. In the nonexact case each aligned pair has one part dividing the other, excluded by the [ratio theorem](row_line_obstructions.md). For an exact row, $7\mid N$, and the exact-row argument below forces the other row to have an entry $u_i M/7$, requiring $7\mid M$. This raises the [previous exponent](proportional_part_pairs.md) from $3/5$ to $2/3$.

## One exact row

Suppose $B$ is exact, with primitive direction $u$ and $\Gamma=\prod_{i<j}\gcd(u_i,u_j)$, where its entries are positive. With $D=J=1$, let $S,T$ be the predecessor products assigned to $A,B$. At each power in $T$ two entries of $B=(N/t)u$ vanish. Since $N/t$ is coprime to $M$, this implies $T\mid\Gamma$. Consequently $S\ge M/\Gamma$.

For **any** integers $c_i,d_i$, including different denominators at different positions, two of $c_i A_i-d_iM$ vanish at each power in $S$. Hence
$$S^2\mid\prod_i(c_i A_i-d_i M). \tag{5}$$
An absolute product below $(M/\Gamma)^2$ forces an exact prescribed predecessor proportion in at least one position. In particular this applies to $c_i=s,d_i=r_i$ at widths $o(N^{2/3})$. No positive part of either row can divide $M$, since each is a multiple of its assigned weight scale, which divides $N$ and is at least three.

## Two and up to nine leading blocks

Let $p,q$ be distinct primes dividing $NM$. For powers $P=p^a,Q=q^b$, $a,b\ge1$, suppose
$$N=2P+h=eQ+k,\qquad 1\le e\le9,\qquad 0\le h<P,\quad 0\le k<Q.$$
Use $s=2,t=e,X=2h,Y=ek$ in (2). Require also
$$(2h+2)^3 e^6<M^2. \tag{6}$$
Suppose one of the following tests holds.

- Powers of $p$ and $q$ dividing $N$ have product greater than $h$.
- $e$ is odd, $2e(h+k)<N$, and either $\min(h,k)^7\le6^6M^4$ or two powers of distinct primes dividing $M$ have product at least $\min(h,k)$.

Then **the full three-row conjecture holds at $N$**.

In a counterexample, use the [all-weight endpoint reduction](all_weights_endpoint_reduction.md). Each designated base must be avoided by a nonendpoint row. In the first test a row avoiding both would have all parts divisible by the specified product, whereas at least one part is at most $h$ by the two upper blocks at $P$. In the second test let $r_i,u_i$ count the upper blocks in that same row. Carry freedom gives $|2A_i-r_iN|\le2h$ and $|eA_i-u_iN|\le ek$. Hence $|(e r_i-2u_i)N|\le2e(h+k)<N$. Coprimality forces all upper blocks at both powers into a single position. The two other parts sum to at most $\min(h,k)$, contrary to the [seventh-power bound](endpoint_pair_bound.md) or the [two-power bound](two_predecessor_pair_sum.md).

The avoiding rows must therefore be different. Take their upper blocks as $r,u$. Their entries sum to $2,e$ and the indicated $X,Y$ bound the discrepancies. The first row cannot be exact since some $r_i=0$. If the second is exact, its primitive target has $\Gamma\le e^3$. By (5) and (6), $2A_i=r_iM$ for some position. The only positive possibilities are divisors of $M$, already excluded.

If both are nonexact, (3) applies. Any zero component of $r$ must index the unique vanishing minor. Thus $r$ is a permutation of $(1,1,0)$ and its two corresponding actual parts are equal. If $u$ has a zero entry, the second row has a part at most $k$. The [opposite-part inequality](cross_parts.md) requires $hk(h+k)>NM$, whereas (2) with the stated bounds implies the reverse inequality. If either remainder is zero a positive small part is already impossible in the corresponding zero position.

Finally suppose $u$ is positive. The aligned parts of $B$ have reduced ratio $v:w$ with $v+w\le e-1\le8$. For every such coprime pair, $vw\le2(v+w)$. When the smaller entry is at least three, the only possibilities are $3:4$ and $3:5$. The [two-ratio bound](row_line_obstructions.md) for $1:1$ and $v:w$ now gives an endpoint common divisor at least $2N(v+w)/(Mvw)>1$. This completes the proof. For fixed $e$, (2) and (6) hold for $h,k=o(N^{2/3})$; the same-row test is a separate hypothesis.

## A single small part against any rational target

Suppose $A$ has a part $a$, $u$ is primitive, and $D=J=1$. Use $s=1,X=a,Y=\max_j|tB_j-u_jN|$ in (2). If the bound holds, then one of the following is necessary.

- $B$ is exact and $t\mid N$.
- $B$ is nonexact, $t\mid M$, and at the two aligned positions $j$ one has $t(N-a)B_j=NM u_j$.

To see this, group $A$ into the two parts $a,N-a$. At every predecessor power assigned to $A$, one of these two parts vanishes because $a$ is zero or one there. The proof of (1) applies with $L=a$, distinguished part $N-a$, and target entry one. For each nonzero $H_j$ it forces $t(N-a)B_j=NM u_j$. Three such equalities would give $N-a=M$ by summing, or $a=1$, impossible by the weight scale. Thus exactly two hold. Since $N\mid (N-a)B_j$, both their $u_j M$ are multiples of $t$. The sum condition gives the same for the third, and primitivity implies $t\mid M$. If two different part values of $A$ satisfy this test, $B$ must be exact, since the displayed equality determines $a$ uniquely.

In particular let $N=dP+h=eQ+k$ at distinct endpoint bases as above, with $s_p(d)\le2$ and **arbitrary $e$**. Require (2) with $s=1,t=e,X=h,Y=ek$, and
$$\max\{\gcd(e,N),\gcd(e,M)\}\le2. \tag{7}$$
The full conjecture holds whenever the two bases cannot be avoided in the same nonendpoint row. Explicit sufficient tests for that last condition are the weight-power product greater than $h$, or $\gcd(d,e)=1$, $de(h+k)<N$ and either pair-sum test above. These follow by the same upper-block argument with $d$ in place of two.

Indeed the first avoiding row has a positive part at most $h$. If an upper block of the second row is zero, it too has a small part, and (2) implies $hk(h+k)\le NM$, contradicting the opposite-part inequality. Otherwise its primitive upper vector has positive entries and sum $t'\ge3$ dividing $e$. The bound still holds on dividing that target by its content and substituting the actual part $a\le h$. The alternatives force $t'\mid N$ or $t'\mid M$, contrary to (7).

## More than 110 quadrillion certified weights

Fix $P=5^{168}$. For $e=3$, take
$$q=42211883872453072075928572096340406768795975134335019729909,\qquad Q=q^2.$$
Let $L=10^{18}$ and take $N=q(3q+m)+1$ for $L\le m<2L$ with $5\mid 3q+m$. There are $K=L/5$ such parameters. **More than $110{,}000{,}000{,}000{,}000{,}000$ of these distinct 118-digit weights satisfy the full conjecture.** At every weight in the interval both remainders $h=N-2P$ and $k=N-3Q$ exceed $N^{13/20}$. This puts both outside the earlier fifth-degree test for these upper blocks. The verifier bounds (2), (6), and $6(h+k)<N$ uniformly by their endpoint values. Both 5 and $q$ divide $M$.

Write $z=(3q+m)/5$. If any prime power $W>2L$ divides $z$, its base differs from $q$ since $0<m<q$. Now $qW\mid M$ and $qW>k$, so the same-row test holds. The values of $z$ form $K$ consecutive integers, each exceeding $2^{193}$. Let $b$ count those with all full prime powers at most $2L$. Their product divides $K!\operatorname{lcm}(1,\ldots,2L)$. For each prime remove a single bad term of highest valuation when counting that prime. At any prime-power level the number of other multiples is at most $\lfloor(K-1)/p^j\rfloor$, since they lie within $K$ consecutive integers. The removed prime power divides the indicated least common multiple.

For $H=2^{61}\ge2L$ this least common multiple is at most $2^{2H}$. Indeed $\operatorname{lcm}(1,\ldots,2n)$ divides $\binom{2n}{n}\operatorname{lcm}(1,\ldots,n)$: for each prime there is at most one new power above $n$, and at that level the binomial factorial floor contributes one. Iterate at powers of two. Since $K!\le2^{58K}$ and $2H<24K$, one gets $193b<82K$ whenever $b>0$. Rounding the strict inequality gives at least $115{,}025{,}906{,}735{,}751{,}296$ qualifying parameters. This counts qualifying parameters without identifying every one. One individually certified member has $m=1108391398016982148$ and $q\,5^{27}\mid M$.

For the arbitrary-quotient result use $e=101$,
$$q=7275028100510113835355475391056051625808232326939586246027,\quad m=4486974951398526273,\quad N=q(101q+m)+1.$$
Here $q\,5^{27}\mid M$, (7) holds with both gcds one, and both remainders exceed $N^{16/25}$. The primality checks for both displayed bases $q$ use recursive multiplicative-order certificates.

## Reproduction

The Python check uses exact trinomial coefficients for 1,803,816 instances of (1) and the same number of grouped-cut instances through weight 32. Of the first group, 1,562,152 require the common endpoint factor in (1). It checks the large certificates with integers. It also exercises the exact equalities and divisor enumeration on a family satisfying the assigned lowest layers. The Lean file checks the algebraic divisibility implication from the retained factors and the short-product implication.

~~~sh
python3 check_mixed_row_minors.py  # uses python-flint
lean MixedRowMinors.lean
~~~
