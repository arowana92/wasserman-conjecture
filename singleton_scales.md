# A part coprime to the predecessor

In a counterexample take the [mandatory endpoint](all_weights_endpoint_reduction.md) and the other two rows $A,C$. Assign the full powers as $UV=N$ and $ST=N-1$, with $U,S$ assigned to $A$. In particular $A$ has all entries divisible by $U$, and $C=Vy$ where the three positive entries of $y$ sum to $U$.

**If a part of $A$ is coprime to $N-1$, then $U\ge1438$. If $U\le5000$, the weight must exceed $10^{1000}$, and the number of possible exceptional weights up to $X$ in this case is $O_\epsilon(X^\epsilon)$ for each $\epsilon>0$.** These assertions also hold if all the full powers in $S$ have their exceptional residue 1 at a single position, including $S=1$.

**If either nonendpoint row has a part equal to a prime power $p^i>1$, then $p\mid N$ and $p^{v_p(N)}>5000$.** This bound has no restriction on $i$, the base-$p$ digit sum of $N$, or the factorization of the cofactor. For a different result under a [five-unit hypothesis](sparse_weight_factors.md) the bound on the full power is 200 million.

## Reducing to finitely many row shapes

At a predecessor power assigned to $A$, the actual residues are $1,0,0$. A part coprime to $N-1$ must be the exceptional one at every such power. Put it first. There are positive integers $k,f,g$ with
$$A=(kS+1,UfS,UgS),\qquad T=k+U(f+g)>2U.$$
Here $k>0$ because the part is a positive multiple of $U\ge3$. Assume $f\le g$. Write $d,e,h$ for the products of powers in $T$ exceptional at the three positions of $C$. They are pairwise coprime, and the normalized entries are multiples of $T/d,T/e,T/h$, respectively. Every bin is greater than two because each such entry is positive and less than $U<T/2$. Order the bins as $3\le d<e<h$. This gives
$$de+(d+e)h\le U,\qquad 27T^2\le U^3.$$
The latter follows from $T^2\mid y_1y_2y_3$. Thus a bound on $U$ yields a finite enumeration of bins, normalized entries, and slopes $k,Uf,Ug$, without bounding $S$.

At each $q^e=w\parallel U$, $S\equiv-T^{-1}\pmod w$ and $TS+1$ must have exact $q$-valuation $e$. This excludes one of the $q$ lifts modulo $qw$. At every $q^e=w\parallel T$, the analogous condition is $V\equiv U^{-1}\pmod w$ with one forbidden lift, since $\gcd(S,T)=1$. The assigned rows must add without carries in their bases. The [finite carry graph](carry_intersections.md#exhaustive-search-by-intervals) reads the fixed low digits and then tests whether any continuation can terminate without a carry. Its states are the remaining multiplication carries of the three affine expressions. Every accepted integer extension is covered. The local tests allow the scalar to vary independently at each prime.

If $U,T$ are odd, either $S$ is even or $V$ is even. In the first case $A$ must avoid 2, and its new units bit at the shifted part is separate. After stripping the power of 2 from $S$, some odd scalar times $(k,Uf,Ug)$ is carry-free. In the second case some odd scalar times $y$ is carry-free. Require at least one possibility.

For $U\le5000$ the bin enumeration gives 345,560 systems $(U,T,y)$, of which 1,950 pass the companion local tests. The first-row and parity tests leave **61 shapes**, the smallest at $U=1438$. An independent enumeration factors each possible $T$ first and partitions its prime powers. It reproduces the entire set of companion systems. Applying factorial floors at the prescribed low digits and then a separate digit-extension implementation gives exactly the same 61 shapes.

## The prime-power part

For a part $p^i$, $U\mid p^i$ forces $U$ to be the full power of $p$ in $N$. The part is coprime to $N-1$ and the parametrization requires $kS+1=p^i$. Among the 61 shapes there are just 28 with $U$ a prime power. When $i\ge16$ the low digits of $S$ satisfy
$$S\equiv-k^{-1}\pmod {p^{16}}.$$
The inverse exists since $k\equiv T\pmod U$. A carry between the low digits of $fS,gS$ contradicts avoidance of $p$ on $A$, after dividing its parts by $U$. In two shapes $p,k$ are odd and $T$ is even, which would give $S$ even. The digit test excludes the other 26. In the digit exclusions a layer $p^j$ with $j\le8$ already suffices. For $v_p(N)\le i<16$ check divisibility by $k$, unitarity, and the forced carries directly. No shape survives. The audit also checks the required base-$p$ carry using the sum of two residues at one layer.

## The remaining unrestricted scales

For each of the 61 shapes, $V\equiv U^{-1}\pmod T$ and all primes dividing $T$ must be avoided by $Vy$. Two independent exhaustive algorithms find no such $2\le V\le10^{1000}$. Intersecting high-digit prefixes in the different bases visits 4,428,789 nodes in total. Covering the same progressions by intervals on which a factorial floor forces a carry uses 590,192 intervals. Since $N=UV$, this gives the stated lower bound.

For each $q\mid T$ count the paths of length 128 in the carry graph for $Vy$ that can reach an accepting state. Every allowable block in base $q^{128}$ occurs on such a path, including any padding by zeros. If the count is $H_q$, the resulting integer missing-digit set has dimension at most $\log H_q/(128\log q)$. The verifier rounds up to hundredths using powers of integers. For 59 shapes two distinct forced bases have bounds summing to at most one. The intersection of the two integer sets then contains $O_\epsilon(X^\epsilon)$ inputs up to $X$. [Glasscock, Moreira and Richter, Theorem B](https://arxiv.org/abs/2007.05480)

Both other shapes have $U=4762$, $T=36049$. The companion bound at 13 is $0.43$, and the first row avoids 2. The [binary affine count](binary_affine_counts.md) gives a block dimension at most $\alpha=\log_2\lambda<0.56$, where $\lambda^3=\lambda^2+1$. For the block-dimension deduction, a block beginning at a reachable carry state is itself an extendible input prefix for the same slopes with those state coordinates as its offsets. The binary bound holds at every such nonnegative offset. The number of states is fixed, so the union of their length-$r$ block languages has size $O(2^{\alpha r})$. Taking $r$ large gives dimension below $0.56$ for a missing-digit superset. Their dimension sum is below one, so the affine form of Theorem B applies to $TS+1$ and $UV$. Summing over the finite list proves the count.

## Bounds without a cutoff on the scale

For any chosen part $a$ in a nonendpoint row, let $F$ be the product of the powers in $S$ at which it has residue zero. Thus $F\mid\gcd(a,N-1)$. Put $B=S/F$. The same row equations become
$$A=(kB+1,UfB,UgB),\qquad k+U(f+g)=FT,\qquad N=FTB+1.$$
Each power of $F$ divides $f$ or $g$. For $U\le K$, $F\le Y$, the number of possible exceptional weights through $X$ is
$$O_\epsilon\bigl(X^\alpha K^{5/2+\epsilon}Y^{2+\epsilon}\bigr),\qquad 1\le K,Y\le X. \tag{1}$$
In particular the bound applies with $F$ replaced by $\gcd(a,N-1)$ in the hypotheses.

To count shapes, put $U$ in a dyadic range $Z/2<U\le Z$. At each $T\le (Z/3)^{3/2}$ partition the full powers of $T$ among three exceptional positions in $3^{\omega(T)}$ ways, allowing empty bins. Their products $d_i$ require $T/d_i\mid y_i$. There are at most $Z^3/T^2$ ordered choices for the positive $y_i$, since $\prod d_i=T$. Each fixes $U$. Split the full powers of $F$ among $f,g$ in at most $2^{\omega(F)}$ ways. Both positive slopes are less than $FT/U$, so after imposing the two divisibilities there are at most $(FT/U)^2/F$ choices per split. These fix $k$. Multiplication cancels $T^2$, giving $O_\epsilon((FT)^\epsilon FZ)$ choices per $F,T$, with the divisor-bound epsilon adjustable. Summing proves the shape count in (1).

For each shape one of the two rows avoids 2. Its input is $B\le X/(FT)$ or $V\le X/U$. The respective slope totals are $FT$ and $U$. In the binary affine count the extra digit delay is at most $\log_2$ of this slope total. It cancels the denominator in the input bound, giving $O(X^\alpha)$ weights per shape with an absolute constant.

There is also a direct bound independent of any cutoff on $N$. The number of possible exceptional weights having a nonendpoint part $a\le Y$ is $O_\epsilon(Y^{3/2+\epsilon})$. For a fixed $a$, the choices $U\mid a$, $S\mid a(a-1)$ number $O_\epsilon(a^\epsilon)$, and $a\ge U\ge3$. The companion requires $T\le (U/3)^{3/2}$ and $ST\equiv-1\pmod U$, so it leaves $O(1+\sqrt U)$ choices for $T$. Each gives a single $N=ST+1$. For parts $a=p^i$ with $P=p^{v_p(N)}\le K$, this argument counts $O_\epsilon(K^{3/2}X^\epsilon)$ possible weights $N\le X$ by summing over $P$ and $i\le\log_2 X$ with $U=P$.

~~~sh
python3 check_singleton_scales.py
python3 check_singleton_scale_audit.py
python3 check_singleton_intersections.py --jobs 8 --floor-digits 1000
~~~
