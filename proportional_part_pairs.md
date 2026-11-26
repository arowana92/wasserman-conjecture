# Exact ratios forced between parts of two rows

Fix primitive positive integer triples $r,u$ with sums $s,t$, giving target proportions $\rho=r/s$ and $\sigma=u/t$. For positive rows $A,B$ of weight $N$, write $M=N-1$, $H_A=\max_i|A_i-\rho_i N|$, and $H_B=\max_i|B_i-\sigma_i N|$. Call a pair of positions $j,k$ aligned in $A$ if $r_j A_k=r_k A_j$, and define alignment in $B$ using $u$. A row is exact when all three target proportions hold.

**For every weight,** if the two coefficients have no common prime dividing $NM$, the following implications hold.

- If neither row has an aligned pair, then $(s t H_A H_B)^3\ge N^3 M$.
- If $A$ is not exact and $B$ has no aligned pair, then $(sH_A)^2(tH_B)^3>2N^2 M$. The roles can be exchanged.
- If $A$ is exact and $H_B=o(N^{2/3})$, then at all sufficiently large weights at least one component of $B$ is **exactly** $u_i M/t$. An explicit test below gives a finite threshold.

Consequently if both rows lie within $o(N^{2/3})$ of their respective targets, at least one row has an aligned pair. If neither row is exact and both lie within $o(N^{3/5})$, **each** has an aligned pair. These restrictions concern the two nonendpoint rows of every counterexample by the [all-weight reduction](all_weights_endpoint_reduction.md).

In particular, consider two rows each within $o(N^{3/5})$ of **$(N/7,2N/7,4N/7)$**, with coordinates permuted independently if desired. For all sufficiently large weights these two trinomial coefficients share a prime of $NM$. Both rows can be quite far from the central or balanced proportions.

## The three minors in each row

Use the [endpoint assignment](endpoint_reduction.md) with $UV=N$ and $ST=M$. The parts of $A,B$ are multiples of $U,V$ respectively. Write $S=S_1S_2S_3$ with the exceptional residue 1 in position $i$ at every power in $S_i$. Use $T_i$ similarly for $B$. For complementary indices $j,k$ we have
$$US_i\mid r_j A_k-r_k A_j,\qquad |r_j A_k-r_k A_j|\le (r_j+r_k)H_A\le sH_A.$$
Also $2S_i<V$ since the normalized positive parts at $j,k$ are multiples of $S_i$ and the third part is positive. At most one minor vanishes unless $A$ is exact. With three nonzero minors this gives $S\le(sH_A/U)^3$. With at least two nonzero minors use $S<(V/2)(sH_A/U)^2$. Multiplication by the corresponding bounds for $T$ proves the first two inequalities, using $UV=N$ without any comparison between the two scales.

The same proof yields quantitative gcd bounds for arbitrary rows. Let $D=\gcd(N,K_A,K_B)$, $J=\gcd(M,K_A,K_B)$ and $L_A,L_B$ be the largest absolute minors in the two rows. Using the [retained powers](endpoint_pair_bound.md) gives $DUV=N$, $JST=M$ and $2S_i<DV$. If all six minors are nonzero then $N^3 M\le D^3J L_A^3L_B^3$. If $A$ is nonexact and $B$ has three nonzero minors then $2N^2M<D^3J L_A^2L_B^3$.

## When one row is exact

Let $A=(N/s)r$ and $\Gamma=\prod_{j<k}\gcd(r_j,r_k)$. We have $s\mid N$ by primitivity. At each factor $S_i$, the two zero positions imply $S_i\mid\gcd(r_j,r_k)$ since $N/s$ is coprime to $M$. Thus $T\ge M/\Gamma$. For $B$ put
$$z_i=2tB_i-2u_i M-t,\qquad Z=\max_i|z_i|\le t(2H_B+1).$$
If
$$T\ge16t^2,\qquad 27tZ^3<8T^2, \tag{1}$$
then some $B_i=u_iM/t$. The stronger conditions obtained by replacing $T$ with $M/\Gamma$ also suffice.

At each assigned power of $T$ every $B_i$ is zero or one. Therefore $z_i^2=t^2+4b_iT$ for integers $b_i$, nonnegative under (1). We have $\sum z_i=-t$. Suppose first that every $b_i>0$ and put $F=\sum_i\operatorname{sgn}(z_i)\sqrt{b_i}$. Each $|z_i|-2\sqrt{b_iT}$ lies between zero and $t^2/(4\sqrt T)$. Consequently $0<|F|<t/\sqrt T<1$. If one signed sum of the three square roots vanishes then one root equals the sum of the other two, and every nonzero signed sum has magnitude at least two. Thus all four signed sums, up to global sign, are nonzero. Their product up to sign is the nonzero integer $2\sum_{i<j}b_i b_j-\sum_i b_i^2$. Each of the other three factors is bounded by $3Z/(2\sqrt T)$. This gives $1<(t/\sqrt T)(3Z/(2\sqrt T))^3$, contrary to (1).

Now some $z_i$ is $-t$ or $+t$. The former gives the claim. In the latter case the other two sum to $-2t$, so the absolute difference of their squares is at most $4tZ<4T$ by (1). It is a multiple of $4T$, hence is zero. Both of these other entries must be $-t$. This proves the assertion. In particular (1) holds eventually for $H_B=o(N^{2/3})$.

## Proportions forming a divisibility chain

In the row with $S\ge V$ or $T\ge U$, every exceptional bin must be nonempty, since a positive normalized part is smaller than the whole assigned predecessor product. In particular **no part of this row divides another**. Indeed if one part divides another then a power with residue 1 at the latter and 0 at the former is impossible.

Suppose the entries in each target triple form a chain under divisibility after sorting. Any aligned pair would then have one actual part dividing the other. With both widths $o(N^{3/5})$, the two-nonexact case is impossible by the two-minor bound. The exact case is also impossible if, for every coordinate of each target, its reduced denominator shares a prime with the full denominator of the other target. The exact row needs its denominator to divide $N$, while at least one coordinate denominator of the opposite row must divide $M$. These hypotheses hold, in particular, for any two chains with the **same full denominator**, including two copies of $1:2:4$.

For that example a fully numerical sufficient test uses $H=\max(H_A,H_B)$. The conclusion holds whenever
$$M\ge1568,\qquad (7H)^5\le2N^2M,\qquad 189(14H+7)^3<2M^2.$$
Here $\Gamma=2$ and all component denominators are seven. The coefficient at the endpoint and the coefficients of these two rows therefore have a common prime.

## Entire weights near two and three prime-power blocks

Let $P=p^a,Q=q^b$ with $a,b\ge1$ and distinct primes $p,q\mid NM$. Write
$$N=2P+h=3Q+k,\qquad 0\le h<P,\quad 0\le k<Q.$$
For all sufficiently large $N$, the **full three-row conjecture holds whenever $h,k\le N^{4/7}$**. In particular the two remainders can both exceed the square root scale. The following entirely explicit conditions suffice:
$$\min(h,k)^7\le6^6M^4,\qquad 6(h+k)<N,$$
$$\max\{(2h)^2(3k)^3,(2h)^3(3k)^2\}\le2N^2M,\qquad
M\ge64,\qquad 54(4h+2)^3<8M^2.$$

Take the two putative nonendpoint rows. A row avoiding the base of $P$ has upper blocks summing to two and lower blocks summing to $h$. For the base of $Q$ the analogous sums are three and $k$. No row can avoid both. If a single part $a$ had upper blocks $j,l$ at the two powers, then $|2a-jN|\le2h$, $|3a-lN|\le3k$ give $|(3j-2l)N|\le6(h+k)<N$. Applying this at every part forces all upper blocks into one position by coprimality. The remaining two positive parts would have sum at most $\min(h,k)$, excluded by the [pair-sum bound](endpoint_pair_bound.md). For this step one may also replace the seventh-power inequality by $\min(h,k)\le W_1W_2$ for powers $W_1,W_2>1$ of distinct primes with $W_1W_2\mid M$, by the [two-power pair-sum theorem](two_predecessor_pair_sum.md).

Use the upper vectors of the two different rows as targets of sums two and three. The widths are at most $h,k$. If a target has two zeros then its two small actual parts have sum at most the corresponding remainder, say $h$. At each of its three bins at least one of those parts is a positive multiple of $US_i$, so $S\le(h/U)^3$. Such a row therefore supplies the cubic bound in the mixed fifth-degree inequalities. For any nonexact companion the quadratic bound still holds since at most one minor vanishes. If neither target has two zeros, the possible upper vectors are $(1,1,0)$ for the first row and $(1,2,0)$ or $(1,1,1)$ for the other, up to permutation. A minor involving a zero and a positive target entry is nonzero since the actual parts are positive. The nonzero entries in each target form a divisibility chain. An aligned pair would therefore be an actual divisible pair, impossible in the row whose predecessor product is at least its normalized sum. This row supplies the cubic bound. In either case the fifth-degree inequalities exclude two nonexact rows.

The only possible exact row has three equal parts $N/3$. It is even in coefficient. Its three pairwise gcds in the primitive direction are one, so the opposite row takes the entire predecessor $M$. Apply (1) at the target $(1,1,0)/2$ or $(2,0,0)/2$. Its proof uses only the sum of its target entries, and $Z\le4h+2$. It forces an actual part $(N-1)/2$, since zero and $N-1$ are both impossible in a positive triple. This also forces an even coefficient. Indeed if the singleton $z=(N-1)/2$ belongs to a binary carry-free row, its ones must be present in $N=2z+1$. Every one of $z$ above position zero then requires a one immediately below it, so $z=2^v-1$. The other two positive parts cannot sum without carries to $z+1=2^v$. This completes the exact case, including weights divisible by three. The asymptotic assertion follows since $5(4/7)<3$ and $3(4/7)<2$.

An example uses $P=5^{168}$, $Q=q^2$ with the prime
$$q=42211883872453072075928572096340406768795975134335019729909,\qquad
N=q(3q+256640162).$$
This weight has 118 digits. Both remainders satisfy $N^4/1000<h^7,k^7<N^4$, and $5\cdot7\cdot11\cdot13\mid M$. The explicit conditions certify all choices of rows without factoring $N$ or $M$. In fact the same $P,Q$ certify all 40,000,000 weights $q(3q+m)$ with $200000000\le m\le300000000$ and $m\equiv2,3\bmod5$. These are the residue classes where 5 divides $N(N-1)$, while $q$ divides $N$ throughout. Over this interval $N,h,k$ increase with $m$. The verifier checks the upper sides of all the sufficient bounds at their largest remainders against their smallest right sides, and also verifies $N^4/1000<\min(h,k)^7$ throughout.

The predecessor-product alternative also works beyond the seventh-power bound on *both* remainders. With the same $P,q,Q$, set $W=5^{16}$, $m=132147841523$ and $N=q(3q+m)+1$. Here $0<m<W$ and $W\mid3q+m$, so $qW\mid M$ and $k=qm+1\le qW$. The fifth-degree and exact-row tests still pass even though $\min(h,k)^7>6^6M^4$. In fact $\min(h,k)^5>N^3/28282$.

## Both nonexact rows in square-root windows

For any fixed $C>0$, suppose both nonexact rows have $H_A,H_B\le C\sqrt N$. At all sufficiently large hypothetical exceptions their normalized shapes obey
$$A/U=(V-\alpha)\rho+\alpha e_i,\qquad B/V=(U-\beta)\sigma+\beta e_j,$$
$$\alpha U+\beta V=1+\alpha\beta,\qquad \alpha\beta<0.\tag{2}$$
Here $e_i,e_j$ are coordinate unit vectors, and the nonzero rationals $\alpha,\beta$ belong to finite sets depending only on $r,u,C$. In particular **all** these weights belong to finitely many explicitly enumerable quadratic polynomials of one integral parameter. This conclusion does not impose divisibility-chain assumptions.

Indeed $U$ divides every $sA_i-r_iN$ and $V$ divides every $tB_i-u_iN$. As neither row is exact, $U,V\ll\sqrt N$, hence both are also $\gg\sqrt N$. All normalized discrepancies $sA_i/U-r_i V$ and $tB_i/V-u_i U$ belong to fixed finite integer sets. By the minor inequalities each row has exactly one aligned pair, yielding the displayed vertex directions for its normalized discrepancies. In each row the two bins for nonzero minors are bounded, since the normalized minors are bounded. The remaining bin $S_i$ is less than $V/2$ and $T_j<U/2$. Since $ST=M\asymp UV$, this gives $S_i\asymp V$ and $T_j\asymp U$.

For any $k\ne i$, $A_k/U=(V-\alpha)\rho_k$ is a positive multiple of $S_i$ and its quotient is bounded. Thus $(V-\alpha)/S_i$ takes finitely many positive rational values. Likewise $(U-\beta)/T_j$ and the products $S/S_i,T/T_j$ take finitely many values. Their combination gives $M=c_0(V-\alpha)(U-\beta)$ with $c_0$ in a fixed finite rational set. As $U,V\to\infty$ this ratio $c_0$ tends to one, so $c_0=1$ for every sufficiently large weight. Expanding proves (2). Equal signs are impossible since $U,V$ grow and $\alpha,\beta$ range over a finite set bounded away from zero. Substituting the resulting affine expression for $V$ into $N=UV$ proves the quadratic assertion.

All finite choices here can be obtained from $s,t,C$. For example $U\le sC\sqrt N$, $V\le tC\sqrt N$ and their product gives the reciprocal lower bounds. These bound each normalized minor and each of the four small bins explicitly. Using $M\ge N/2$ gives corresponding lower bounds for $S_i/V,T_j/U$. Thus neither the finite lists nor the existence of thresholds in this argument require an ineffective integral-point theorem.

## An infinite family passing all lowest endpoint layers

The line $\alpha=1,\beta=-2$ in (2) has $U=2V-1$, $N=UV$ and $M=(V-1)(2V+1)$. Take
$$V=11972031+75150075h,\quad h\ge0,\qquad m=(V-1)/35,\quad l=(2V+1)/143.$$
The two positive rows
$$A=U(23m+1,5m,7m),\qquad B=V(119l-2,11l,13l)$$
are in square-root windows at the fixed target proportions $(23,5,7)/35$ and $(119,11,13)/143$. Every full prime power of $N$ or $M$ has a row whose residues at that layer sum to the weight residue. Nevertheless **both coefficients are even** for every $h$.

To check the assignments take $S=V-1$, $T=2V+1$. They are coprime since $3\mid V$, and $\gcd(2V-1,V)=1$. The parameters obey $m\equiv3\bmod35$, $l\equiv131\bmod143$. The powers in $U,V$ have all residues zero in the corresponding rows. At every full power dividing $m$, the second and third parts of $A$ are zero and its first part is $1\bmod m$ since $U\equiv1\bmod m$. The factors 5 and 7 occur once in $S$, with exceptional positions three and two respectively. At every full power dividing $l$, the first part of $B$ is $1\bmod l$ since $-2V\equiv1\bmod l$. At 11 and 13, which occur once in $T$, the exceptional positions are respectively three and two by the residues for $l$. This verifies every assignment. The two unshifted parts in each row have the same 2-adic valuation and hence a common binary one. The higher binary digits thus prevent these rows from being a counterexample even though every assigned lowest layer passes. The first weight is $286659040557891$.

## Reproduction

Use `python3 check_proportional_part_pairs.py` with python-flint installed. The verifier uses exact trinomial coefficients for the quantitative minor bounds, independently tests the square-root sign argument with modular triples, and factors the endpoints in samples from the explicit progression. It verifies the tables of low residues and checks the actual parity by factorial valuations.

A further explicit progression passes every lowest-layer check and has an odd coefficient in one row for every permitted exponent. In the other row the first ternary carry occurs at 9. [Full binary carry example](binary_admissible_endpoint_rows.md)
