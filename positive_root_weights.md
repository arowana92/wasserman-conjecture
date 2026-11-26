# A positive root of the predecessor polynomial

Let $F\in\mathbb Z[X]$ be nonconstant with nonnegative coefficients and $F(0)=0$. The equation $F(x)=1$ has a unique positive real solution $\xi$. Let $H$ be its primitive irreducible polynomial over $\mathbb Q$, with positive leading coefficient. It divides $F-1$ in $\mathbb Z[X]$.

**Common divisor of rows avoiding the base prime.** There is an effectively computable positive integer $c_F$ such that, for every prime power $P=p^a$ at which $H(P)\ne0$ and $N=F(P)\ge3$, the integer
$$\frac{|H(P)|}{\gcd(H(P),c_F)}$$
divides **every** positive trinomial coefficient of weight $N$ that is coprime to $p$. In particular, for every such polynomial $F$, the exact central cut at $N=F(2^a)$ is universal for all sufficiently large $a$. Together with the [displacement bound](universal_cut_offsets.md), this includes the center itself in a whole interval of universal cuts of width $\gg_F N^{2/7}$.

**Full conjecture on polynomial families.** Let $G$ be monic, irreducible over $\mathbb Q$, of positive degree, with nonnegative integer coefficients and $G(0)>0$. For $m=1$ or $2$ the full conjecture holds at every sufficiently large $N=mP G(P)$, uniformly over all prime powers $P=p^a$. The threshold is effective in $G,m$. This covers reducible predecessors $mXG-1$. Monicity can also be replaced by primitivity and $G(0)\ge\operatorname{lc}(G)$, where $\operatorname{lc}$ is the leading coefficient. For general fixed $m\ge1$ the conclusion holds on the bases $p$ for which either $m/p^{v_p(m)}\le2$, or some nonzero coefficient of $mG$ has base-$p$ digit sum at most two.

**An explicit two-exponent family.** For $a,j\ge1$ and $p$ prime, write $P=p^a$, $Z=p^j$ and put
$$ N=(P^2+2P+2)Z+2P(P+1). \tag{1} $$
The full conjecture holds at every weight in (1) in each of the following cases:

- $p=2$, with **all** $a,j$.
- $j\ge8a$ or $a\ge24j$, with independent exponents at every prime.
- $j/a$ is an integer at least four, at every prime.

For example, setting $j=ka$ for any even $k\ge4$ gives a predecessor polynomial divisible by $X+1$. The general family theorem applies since $X^{k+1}+2X^k+2X^{k-1}+2X+2$ is Eisenstein at 2. The explicit result covers every evaluation and includes nonintegral exponent ratios.

## Why a positive root gives a common divisor

A row $A$ avoiding $p$ distributes the base-$p$ digit units of $N$ among its three parts. Therefore $A_i=F_i(P)$ for three nonzero polynomials with nonnegative integer coefficients summing coefficientwise to $F$. This remains true for overlapping digit blocks: start with the monomials counted with multiplicity, and merge $p$ equal powers of $p$ until the base-$p$ expansion of $N$ is reached. Assign every original unit to the part receiving its merged unit. There are only finitely many polynomial splits, independently of the base and exponent.

Every possible $F_i$ has $F_i(\xi)>0$, so its resultant with $H$ is nonzero. Take $c_F$ to be a common positive multiple of the absolute values of these resultants, or take 1 if there are no splits. By the resultant identity, $\gcd(H(P),F_i(P))\mid c_F$. Put $M=N-1$. At $q^e\parallel M$, if $l=v_q(K_A)<e$, at least two of the parts are zero modulo $q^{e-l}$, by the [retained-power lemma](endpoint_pair_bound.md). Since $H(P)\mid M$, this gives $l\ge v_q(H(P))-v_q(c_F)$ and proves the common-divisor assertion.

More generally, at any fixed reduced fraction $r/t$ with $0<r<t$, $t\mid N$ and $p\mid t$, the cut $b=rN/t$ is universal at all sufficiently large $P$ in the family. Indeed $t\mid\binom Nb$ and $M/\gcd(M,r(t-r))\mid\binom Nb$ by the [cut factorial identities](divisor_cuts.md#factorial-identities-for-a-cut). If $p\mid K_A$ it supplies the prime. Otherwise $\binom Nb$ and $K_A$ share the integer $|H(P)|/\gcd(H(P),c_F r(t-r))$, which grows without bound. Taking $r/t=1/2$ gives the asserted binary result.

## Removing the irreducibility assumption on the predecessor

Take $F=mXG$ and write $b=\operatorname{lc}(G)$. Since its predecessor has constant $-1$ and leading coefficient $mb$, any positive rational root would be at least $1/(mb)$ by the rational root test. But $F(1/(mb))=G(1/(mb))/b>1$ since $G(0)\ge b$. Thus $e=\deg H\ge2$ without any further factorization hypotheses.

Write the polynomial split of a $p$-avoiding row as $A_i=P L_i(P)$, where the $L_i$ are nonzero and sum coefficientwise to $mG$. We can obtain a split which is also carry-free at each coefficient of $mG$ in base $p$: first form the base-$p$ units of each coefficient, then merge across powers of $P$ if needed. If every $L_i$ is a multiple of $G$, the degree bound and primitivity give $L_i=c_i G$ with positive integers $c_i$ summing to $m$. This is impossible for $m\le2$. More generally, write $v=v_p(m)$ and use any coefficient of $G$ coprime to $p$, which exists by primitivity. The low $v$ digits at this coefficient of $mG$ vanish, so each $c_i$ is divisible by $p^v$, so $m/p^v\ge3$. At every nonzero coefficient of $G$ the three positive summands also need at least three base-$p$ units. This excludes the multiples case under each of the stated hypotheses.

Therefore one $L_i$ has a nonzero resultant with $G$. Uniformly over the finite list of splits,
$$\gcd(A_1,A_2,A_3)\le P\gcd(mG(P),L_i(P))\ll_{G,m} P.$$
For any second trinomial row $B$ put $D=\gcd(N,K_A,K_B)$ and $J=\gcd(M,K_A,K_B)$. With the [retained factors](endpoint_pair_bound.md), $U$ divides the three parts of $A$, $JST=M$, and $27T^2\le (DU)^3$. At each power in $S$ two parts of $A$ vanish, giving $\gcd(S,H(P))\mid c_F$. Thus $|H(P)|\le c_F JT$. It follows that
$$DJ\gg_{G,m} P^{(2e-3)/3}.$$
This grows for $e\ge2$. In an endpoint counterexample some nonendpoint row must avoid $p$, since $p\mid N$. The [all-weight endpoint theorem](all_weights_endpoint_reduction.md) now proves the full assertion for large $P$. The finite split lists, integer resultants and leading term of $H$ make the threshold effective. The earlier [polynomial content theorem](polynomial_content.md) gives further results, including an explicit norm bound under its irreducible-predecessor hypothesis.

## The binary proof at all exponent pairs

Write $B=P^2+2P+2$ and $C=2P(P+1)$. At $p=2$ the weight $BZ+C$ is the sum of five powers of two with exponents
$$ a+1,\quad 2a+1,\quad j+1,\quad j+a+1,\quad j+2a. $$
In a counterexample some nonendpoint row distributes these five units without a binary carry, including any mergers when exponents coincide. At least one part is a pure power of two. Assigning the prime 2 to this row gives $U=2^{v_2(N)}$, since $U$ divides that part. Each predecessor power in $S$ must take its residue 1 at the pure part, because this part is coprime to $M$. If a second part is pure then $S=1$. Otherwise the other two parts are each sums of two units with gaps $d_1,d_2$, so
$$S\mid\gcd(1+2^{d_1},1+2^{d_2}).$$
For five ordered exponents $e_1\le\cdots\le e_5$, the smaller of the two gaps is at most $e_4-e_1$. Indeed at least one of the two pairs omits $e_5$.

For $j>a$ this yields $S\le Z+1$ and $U=2P$. For $j<a$ it yields $S\le P^2/Z+1$ and $U=2Z$. In both cases $M>P^2 Z$, and $27M^2>U^3 S^2$: use $Z+1\le2Z$ in the first case and $P^2/Z+1\le2P^2/Z$ in the second, with $27P,27Z>32$. The same bounds cover $S=1$. When $j=a$ the weight is $P(P+2)^2=P^3+4P^2+4P$. It has at most three binary ones, so any odd trinomial must consist of three pure powers, giving $S=1$. For $P\ge4$ we have $U=4P$ and $27M^2>U^3$. At $P=2$ the weight is 32, a power of two. This proves the family at every binary exponent pair.

## Two independent exponents at every prime

Suppose $j\ge8a$. Allocate the coefficients of the two blocks among three parts, writing $A_i=b_i(P)Z+c_i(P)$ with $\sum b_i=B$, $\sum c_i=C$ as polynomials. The merging argument again applies. For each of the four choices $\delta=0,e_1,e_2,e_3$ put $c'_i=c_i-\delta_i$. The modulus assigned to this vertex is respectively $U$ or the bin $S_i$ of predecessor powers with exceptional position $i$. It divides all three values $A_i-\delta_i$, hence all evaluations of
$$ b_i c'_l-b_l c'_i. $$
Its part coprime to $p$ also divides every $b_i(P)$ where $c'_i$ is the zero polynomial, and every $c'_i(P)$ where $b_i$ is zero. In all these bounds one may strip factors of $X$. The [polynomial gcd and resultant procedure](polynomial_prime_powers.md#an-effective-test-including-reducible-polynomials) gives a polynomial bound divisible by this common numerical divisor. All predecessor bins are coprime to $p$. For $U$, restore a factor $2P$, since the exact $p$-power in $N$ is $P$ at odd $p$ and $2P$ at $p=2$. Taking the integer-polynomial lcm of the three bin bounds gives a bound for $S$, by their coprimality. Call the resulting signed bounding polynomials $u,s$.

A zero at an actual $P$ needs separate care. If a common polynomial bound vanishes, all the corresponding determinants vanish there. For $\delta=0$, proportionality of the two numerical vectors gives $A_i=(b_i/B)N$. Each $b_i$ is then positive and a multiple of $B/\gcd(B,C)$. But $\gcd(B,C)\le2$ because $B\equiv1\bmod(P+1)$ and $\gcd(B,2P)\le2$. There cannot be three such terms. For $\delta=e_l$, every $A_i-\delta_i$ is positive in a counterexample and equals $(b_i/B)M$. Now $g=\gcd(B,C-1)$ divides 13, since $2B-(C-1)=2P+5$ and $4B-(2P+5)(2P-1)=13$. The three positive parts after shifting require $g=13$. Write them as $h_i M/13$ with positive $h_i<13$ summing to 13. At the full 13-power in $M$, neither of the two unshifted parts vanishes because its multiplier is less than 13. Thus the first row carries at 13. The companion would have to avoid 13, requiring two normalized parts divisible by 13 and hence $U\ge27$. Yet $U$ divides the two unshifted $h_i$, since $M/13$ is coprime to $N$. This excludes the zero case.

There are 575 coefficient splits up to swapping the two positions without the leading monomial. For each, put $M_0=B X^8+C-1$. The verifier proves
$$729 M_0(P)^4-u(P)^6 s(P)^4>0\qquad(P\ge8)$$
by expanding at $P=8+Y$ and checking all coefficients and the strict constant term. Since $Z\ge P^8$, this contradicts the necessary size bound. For prime powers $P<8$, the code directly tests the 391 residual shape and exponent cases where the polynomial inequality alone is insufficient. With the actual row it uses
$$ U\le\gcd_i A_i,\qquad S\le\gcd\big(M,\operatorname{lcm}_{i<l}\gcd(A_i,A_l)\big). $$
All 391 contradict $27M^2\le U^3 S^2$.

For the other order of scales, suppose $a\ge24j$. As a polynomial in $P$ the weight is $(Z+2)P^2+2(Z+1)P+2Z$. Both this quadratic and its predecessor are irreducible over $\mathbb Q$ since their discriminants divided by four are $2-(Z+1)^2$ and $3-Z-Z^2$. The weight has content at most two, while its predecessor has content one. A $p$-avoiding row has a polynomial split in $P$ as before. At least one part is coprime over $\mathbb Q$ to the weight, since three positive multiples of its primitive part cannot sum to it. Every part is coprime to the predecessor: a multiple would equal the entire predecessor by the leading coefficient bound and primitivity, leaving two positive parts summing to one. The quadratic coefficient vectors, including the parts, have Euclidean norms below $R=4(Z+2)$. The Sylvester determinant bound gives $U\le R^4$ and $S^2\le R^{12}$, using two zero slots at each power of $S$. But $\sqrt{27}M>R^{12}$, since $P\ge Z^{24}$, $Z+2\le2Z$, and $Z\ge2$. This contradicts the size bound.

## The four smaller integral ratios

For $k=j/a=4,5,6,7$, put $G_k=(B X^k+C)/X$. The verifier takes all 575 coefficient splits of $G_k$ and computes the existing [single-variable bounds](polynomial_prime_powers.md#an-effective-test-including-reducible-polynomials) for $U,S$. In every split all coefficients of $729(XG_k-1)^4-u^6s^4$ are nonnegative after translating by 34 and the constant is positive. The 73,600 shape evaluations for integers $2\le P\le33$ all fail the necessary size inequality by numerical gcds. Larger integral ratios are already covered by $j\ge8a$.

The verifier also checks the polynomial divisibility bounds against 10,350 direct numerical gcds, the quadratic identities, and 38,340 binary pairings at small exponents. Run it with python-flint installed.

~~~sh
python3 check_positive_root_weights.py
~~~
