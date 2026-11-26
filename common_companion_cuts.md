# Two exceptional cuts against the same trinomial

Write $M=N-1$, $C_a=\binom Na$, and let $K_B$ be the coefficient of a positive trinomial row $B$ of weight $N$. Call $a$ exceptional against $B$ when $\gcd(NM,C_a,K_B)=1$. Identify $a$ with $N-a$ when counting unoriented cuts. Call a cut universal if it is exceptional against no positive trinomial of that weight.

**At $N=4{,}218{,}421{,}148$, the cuts $141{,}136{,}768$ and $1{,}393{,}201{,}304$ are both exceptional against the same trinomial.** Among the 2,099,636,978 weights for which $M$ has three distinct prime factors and each full prime power is at most 20,000, this is the unique weight with that property. The companion and the cut pair are also unique at this weight, up to permutation and reflection.

**At every weight in this search, for every pair of positive rows $A,B$, at least one part $a$ of $A$ has $\gcd(NM,C_a,K_B)>1$.** The choice of part may depend on $B$. This strengthens the [three-factor computation for full trinomial coefficients](endpoint_factors.md).

**At every weight $N$, the total number of cuts that are exceptional against any trinomial is $O_\varepsilon(N^{1/3+\varepsilon})$ for each $\varepsilon>0$.** In particular, all but $O_\varepsilon(N^{2/3+\varepsilon})$ positive rows at that weight have a universal cut. The estimates are uniform in $N$, including weights outside the search.

**The minimum prime certificate at the example weight has size five.** One of its primes lies outside $N(N-1)$, as is necessary for any certificate at this weight.

## The repeated failure

Put $q=1{,}054{,}605{,}287$. The factors are $N=4q$ and $M=269\cdot3137\cdot4999$, with all four odd factors prime. Take $B=(q,q,2q)$. The exact factorial valuations at every endpoint prime are below.

| Coefficient | 2 | $q$ | 269 | 3137 | 4999 |
|---|---:|---:|---:|---:|---:|
| $K_B$ | 50 | 0 | 2 | 2 | 1 |
| $C_{141136768}$ | 0 | 1 | 0 | 0 | 0 |
| $C_{1393201304}$ | 0 | 1 | 0 | 0 | 0 |

The four exceptional oriented cuts are $141136768,1393201304,2825219844,4077284380$. No three, allowing repetitions, sum to $N$. They also form an antichain for binary digit containment. For the two smaller cuts $a,b$, the four intersections of the binary one-positions of $a$ or $N-a$ with those of $b$ or $N-b$ all contain a bit. Their smallest bit values in the order $(a,b),(a,N-b),(N-a,b),(N-a,N-b)$ are $128,256,8,4$. The two orientations of a single cut have disjoint nonempty supports. Thus no two cuts can be the partial sums of a binary carry-free positive row. Since $2\mid K_B$, this rules out the [comparability criterion for a counterexample](cut_relaxation.md#a-criterion-equivalent-to-the-original-problem) for this $B$.

## Five primes are necessary and sufficient at the example

A minimum list certifying every choice of three rows at $4{,}218{,}421{,}148$ is
$$\{269,3137,4999,1{,}054{,}605{,}287,4{,}218{,}421{,}093\}.$$

Here is a sufficient condition yielding this list. Suppose $N=hQ$, with $h\in\{3,4\}$ and $Q$ a prime power of base $q$. Let $M=w_1w_2w_3$ consist of three full powers on distinct prime bases $p_i$. Suppose the three CRT singleton residues $e_i$ have sum different from $N$. Suppose also that a prime power $R$ of base $r$ satisfies $N/2<R\le N$ and $N-R\le\min w_i$. Then $q,r,p_1,p_2,p_3$ form a certificate.

Indeed, any row avoiding $q$ consists of multiples of $Q$. After dividing by $Q$, the only positive partition is $(1,1,1)$ or $(1,1,2)$. Call the corresponding actual row $B$. Two entries of $B$ equal $Q$ and are coprime to $M$, so $B$ carries at each $p_i$. It also carries at $r$, since a row avoiding $r$ must have an entry at least $R$ while $\max B\le N/2$. Any row avoiding $r$ has two positive small parts with sum at most $N-R$, so neither is divisible by any $w_i$; that row carries at every $p_i$. Lastly every row carries at some $p_i$: a row avoiding them all must have each of its three positive parts less than $M$ and receiving an exceptional residue, forcing the parts to be the three $e_i$. In a triple missing $q$ and $r$ as common divisors, the respective avoiding rows are distinct, both carry at every $p_i$, and the third row carries at one. This proves the certificate.

In the example $R=r=N-55$ is prime, and the three residues have sum $2M+1$. For the lower bound on certificate size, every prime is individually avoidable. The endpoint avoids primes outside $NM$, the row $B$ avoids $q$, and the following witnesses cover the remaining four primes. Each displayed row avoids the two indicated primes, as can be checked by factorial floors.

| Pair avoided | First two parts (the third is their complement in $N$) |
|---|---|
| $2,269$ | $5380,10760$ |
| $2,3137$ | $12548,50192$ |
| $2,4999$ | $459908,559888$ |
| $269,3137$ | $1,843853$ |
| $269,4999$ | $1,1344731$ |
| $3137,4999$ | $1,15681863$ |

If one row avoids two primes of a four-prime list, choose one avoiding row for each remaining prime to defeat the list. Thus a four-prime certificate would contain at least three endpoint primes, since the endpoint row avoids every prime outside $NM$. At least two of those three would belong to $\{2,269,3137,4999\}$, and the table excludes every such pair. Smaller lists can be enlarged to four. In fact a prime outside the endpoint is necessary in a certificate of **any** size at this weight: the row $B$, the witness for $269,3137$, and the witness for $2,4999$ have no common endpoint prime. The displayed list supplies $r$ for this triple, with valuation exactly one in each coefficient.

## All-weight count in any interval

Fix $B$. Assign to it every prime of $NM$ that its coefficient avoids, and assign the remaining primes to the cut. Let $V,T$ be the respective products of full powers from $N,M$ assigned to $B$, and let $U=N/V$, $S=M/T$. By Kummer's criterion, $B=V(y_1,y_2,y_3)$ with $\sum y_i=U$. At a full power $w$ in $T$ two of these integers are divisible by $w$, because the actual parts have residues $0,0,1$ and $V$ is invertible there. Therefore
$$T^2\mid y_1y_2y_3,\qquad 27T^2\le U^3,\qquad US\ge3 M^{2/3}S^{1/3}.$$
Every exceptional index satisfies $a\equiv0\pmod U$ and $a\equiv0$ or $1\pmod w$ for each full power $w$ in $S$. These are at most $2^{\omega(S)}$ residue classes modulo $US$. An interval of length $H$ meets each class in at most $1+H/(US)$ integers. Summing over all unitary divisors $U$ of $N$ and $S$ of $M$ bounds the number of oriented indices in the interval exceptional against **any** companion by
$$2^{\omega(N)}3^{\omega(M)}\left(1+\frac{H}{3M^{2/3}}\right).$$
Here $\omega$ counts distinct primes and the sum of $2^{\omega(S)}$ over the unitary divisors of $M$ is $3^{\omega(M)}$. The prefactor is $O_\varepsilon(N^\varepsilon)$ by the divisor bound. Taking $H=N$ proves the cut estimate. If a row has no universal cut, both its first two parts belong to this set, and they determine the third. Applying the cut estimate with $\varepsilon/2$ proves the row estimate.

More generally, for a companion with $k\ge2$ positive parts, $T^{k-1}$ divides the product of its normalized parts, giving $US\ge k M^{1-1/k}S^{1/k}$ by the arithmetic-geometric mean inequality. The same proof bounds the total number of cuts exceptional against any such companion by $O_{k,\varepsilon}(N^{1/k+\varepsilon})$ for each fixed $k$. For $k=2$, this gives a square-root bound at every weight for the cuts that can fail against another binomial.

Two distinct exceptional unoriented cuts against the same $B$ also satisfy a spacing restriction. Orient them so that $0<a<b\le N/2$. Put $d=b-a$, $h=N-a-b$, and $g=\gcd(d,h)$. At each full power $w$ in $S$, exactly one of $d,h$ is zero and the other is $1$ or $-1$ modulo $w$. Thus $\gcd(g,S)=1$ and $S\mid dh/g^2$. Since $U\mid g$, the bound on $T$ gives $$(dh)^2\ge27g M^2.$$ This applies even if the two cuts are incomparable in some required base.

## Why the search covers every common companion in its domain

Suppose first that $T$ contains one or two of the three full powers of $M$. At most two positions of $B/V$ receive exceptional residues, so one positive entry is divisible by $T$. This gives $T<U$ and $S>V$ whenever there is an exceptional cut. Indeed $V\ge2$, since $U=N$ would leave no interior multiple of $U$. Now $x=a/U$ lies in $(0,V)$ and is zero or $V$ at each full power in $S$. The empty and full choices are impossible since $S>V$. With at most two powers in $S$, there is at most one remaining choice up to reflection.

If $T=M$, a companion allowing two distinct unoriented cuts must have $V\ge3$, since for $V=2$ only $a=U$ is possible. Such a companion would avoid all three predecessor primes and full powers of $N$ of product at least three. Its three parts must be distinct CRT singleton residues, one for each predecessor power. Each part is a multiple of $V$ lying strictly between 1 and $M$, so having no exceptional residue would make it a multiple of $M$, which is impossible. The [three-factor row enumeration](endpoint_factors.md) excludes this case throughout the reported domain.

It remains to enumerate cut pairs that avoid every predecessor prime. For each full power $w_i$ of $M$ form
$$e_i=(M/w_i)((M/w_i)^{-1}\bmod w_i),\qquad 0<e_i<M.$$
A cut exceptional against a companion has $U\ge3$, so the boundary choices $1,N-1$ are impossible. The six remaining predecessor residue patterns are represented up to reflection by the three $e_i$. Test each binomial $C_{e_i}$ at all three predecessor primes. For each passing pair $x,y$, compute the product $u$ of the full powers in $N$ whose bases **both** binomials avoid. Every companion must avoid every prime of $v=N/u$, so it is $v$ times a positive partition of $u$ into three parts. Conversely, every such row avoiding all primes of $v$ witnesses the pair, since both cuts already avoid every prime of $M$. This proves the sufficiency of the last test as well as completeness.

| Test | Count | Unit |
|---|---:|---|
| Three full powers on distinct bases, each at most 20,000 | 2,099,636,978 | weights |
| At least two binomials avoid all predecessor primes | 43,080,750 | weights |
| Two such cuts have $u\ge3$ | 75 | pairs at 73 weights |
| Unordered positive partitions examined across those pairs | 54,435 | rows counted per pair |
| Both cuts are exceptional against the same row | 1 | pair and row |

For 73 pairs the prime 2 alone eliminates all possible rows. For one further pair the prime 3 does so. The remaining pair is the displayed example.

For the assertion about every pair $A,B$, a single exceptional cut up to reflection can supply all three parts of $A$ only if it equals $N/3$. In that case $3\mid V$, because $a=N/3$ is a multiple of $U$. Also $\gcd(a,M)=1$, so $a\equiv1$ at every power in $S$, forcing $S\mid 3a-3=N-3$ and $S\mid2$. If $T$ contains at most two powers, $T<U$ and $S>V\ge3$, a contradiction (the same argument includes $T=1$). The remaining case $T=M$, $V\ge3$ was excluded by the trinomial enumeration. At the one weight and companion with two cuts, the explicit four orientations cannot sum to $N$ in triples. This proves the stated finite strengthening.

An intermediate test can fail even when all three CRT residues sum to $N$ and all their binomials pass. At $N=796{,}691{,}358{,}564$ the residues are $44{,}733{,}620{,}930$, $372{,}887{,}368{,}304$, and $379{,}070{,}369{,}330$. All three binomials avoid $31,5521,11657,12379$, with the latter three the predecessor primes. The combined trinomial has valuations $2,0,0,1$ at those primes. A carry can therefore be present in a full row while absent from each of its three binomials.

## Reproduction

The C++ program computes cut residues using a table of pairwise inverses and checks Kummer digits. The Python program independently reconstructs all three inverses and uses factorial floors. It checks the entire initial enumeration at bound 500, every retained pair at the full bound, every possible companion of those pairs, and the successful example by direct Lucas residues as well. It uses only the Python standard library. The original trinomial mode supplies the separate $T=M$ exclusion.

~~~sh
mkdir -p .work
c++ -std=c++17 -O3 check_three_factors.cpp -o .work/three_factors
.work/three_factors 20000
.work/three_factors 500 --cut-pairs
.work/three_factors 20000 --cut-pairs > .work/cut_pairs
python3 check_common_companion_cuts.py .work/cut_pairs
~~~
