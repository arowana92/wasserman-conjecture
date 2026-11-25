# Products of translates of a prime power

Fix $n\ge2$, distinct integers $0=h_0<h_1<\cdots<h_{n-1}$, and a positive integer $d$. Write $F(X)=\prod_j(X+h_j)$ and $P=p^a$, with $p$ prime and $a\ge1$. For the two asymptotic assertions below, one may take $P>(2dF(1))^{\max(83,18n)}$. In this explicit form the multiplier and the shifts can grow with $P$.

**Two consecutive weights.** Suppose that for every $j$ there is a prime power $w_j\mid P+h_j$ with
$$ (P+h_j)/w_j\le P^{1/(2n)}. \tag{1}$$
For every sufficiently large such $P$, Wasserman's full three-row conjecture holds at **both** $dF(P)$ and $dF(P)+1$. In particular this holds whenever all the $P+h_j$ are prime powers. At the successor alone, (1) is needed only at the largest shift.

**Successors with unrestricted factors.** Suppose $n\ge4$, or $n=3$ and the shifts are different from $(0,1,2)$, or $n=2$ and $h_1\ge2$. At every sufficiently large $P$, without (1), any counterexample at $dF(P)+1$ must have a row
$$(l_1 F(P),\ l_2 F(P),\ l_3 F(P)+1),\qquad l_i\ge1,\qquad \sum l_i=d. \tag{2}$$
Consequently the conjecture holds throughout these families for $1\le d\le363$. More generally it holds if every divisor of $d$ is excluded by the [shifted-denominator theorem](shifted_directions.md). It also holds for every $d\le1000$ when the shifts include 19 consecutive integers. The 21 surviving shapes in that theorem all require a prime from $7,11,13,17,19$ dividing $N$. All five primes divide $F(P)$ when the shifts contain 19 consecutive integers.

**All prime powers in explicit families.** For every $1\le h<k\le5$, both $P(P+h)(P+k)$ and $P(P+h)(P+k)+1$ satisfy the conjecture for **every** prime power $P$. The same unrestricted statement holds for the two successors $P(P+1)(P+2)(P+3)+1$ and $P(P+1)(P+2)(P+4)+1$. No factor hypothesis on the translates is needed in these computations.

## Polynomial rows and their assignments

Use the [all-weight endpoint reduction](all_weights_endpoint_reduction.md). A counterexample consists of the endpoint and two further rows $A,C$. Assign every full prime power of $N(N-1)$ to one of the latter rows avoiding its base. Let $U,V$ be the assigned products from $N$ and $S,T$ those from $N-1$, in row order $A,C$. Thus $UV=N$, $ST=N-1$, and $U,V\ge3$. The row $C/V$ consists of three positive integers of sum $U$. Each power in $T$ divides two of these integers, so $27T^2\le U^3$.

Choose $A$ to receive $p$. For $P$ larger than the coefficients of $N(X)=dF(X)+\epsilon$, $\epsilon\in\{0,1\}$, the parts of $A$ have the form $B_i(P)$ for fixed choices among finitely many integer polynomials. They have nonnegative coefficients summing coefficientwise to $N(X)$. This follows by grouping the carry-free base-$p$ addition into blocks of $a$ digits. Every $B_i$ is nonconstant since $U$ divides its positive value and the only possible positive constant part is 1. At $X=0$ their values are all zero for $\epsilon=0$, or a coordinate unit vector for $\epsilon=1$.

Two factorization facts keep the other assignments under control. If $dF-1=HK$ with nonconstant integer factors, $H(-h_j)=-K(-h_j)=\pm1$ at all $n$ roots. The polynomial $H+K$ of degree less than $n$ would vanish identically, contradicting the positive leading coefficient of $dF-1$. Thus $dF-1$ is irreducible over $\mathbb Q$, by primitivity. For $dF+1$, the same argument gives $H=K$. Thus $dF+1$ is an irreducible polynomial or the square of one. If a square root factored further, some factorization of $dF+1$ would have unequal factors. Normalize the irreducible factor to positive leading coefficient.

At the unshifted weight, $S$ is bounded independently of $P$. For each $i$, $dF-1$ is coprime to $B_i$ as polynomials. Their degrees are at most $n$, and their constant terms rule out a nonzero proportionality. Each full power in $S$ divides at least two of the $B_i(P)$ and also divides $dF(P)-1$. Hence $S^2$ divides the product of the three nonzero absolute resultants.

At the shifted weight, $U$ is bounded independently of $P$. Write $dF+1=H^e$ with $e=1$ or 2. Choose $B_i=H^j L$ with $j<e$ and $H\nmid L$. Such a part exists, since $N(X)$ is primitive and no part equals the whole polynomial. If $r^v$ is the full $r$-power in $H(P)$ and $r\mid U$, then $r^{ev}\mid B_i(P)$. Thus $r^v\mid L(P)$, which implies $U\mid |\operatorname{Res}(H,L)|^e$. Now $T$ is bounded by $27T^2\le U^3$. All these bounds are effective over the finite collection of polynomial rows.

## The values at the roots of the predecessor

At the successor, fix a shift $h$ and factor $Y=P+h$ among $T$ and the three groups in $S$ according to the position of the residue 1. This uses the full powers in $N-1$ even if several translates share a base. In group $k$, the portion of $Y$ divides each $B_i(-h)-[i=k]$ since $P\equiv-h\pmod Y$. If none of the three unit vectors equals $(B_i(-h))_i$, each group is bounded by a nonzero absolute value from its three tests. The remaining portion divides $T$. This bounds $Y$, a contradiction for large $P$. Thus **at each of the $n$ roots exactly one polynomial has value 1 and the others vanish**.

For two integral nodes differing by at least 2, these unit vectors agree. The difference of any two evaluations of the same integer polynomial is divisible by the node difference. The graph on the shifts with edges for differences at least 2 is connected if the extreme shifts differ by at least 3. The extreme shifts then join each other and every intermediate shift joins at least one of them. For two shifts a single edge suffices. Consequently all root vectors agree under the hypotheses of the unrestricted-factor assertion. Subtract this vector and divide by the monic polynomial $F$. The quotients are integer constants since each part has degree at most $n$. They are positive, including in the shifted slot because a part of 1 is impossible. This proves (2), whose shifted denominator is $d/\gcd(l_1,l_2,l_3)$.

## Carries after a change of origin

Fix a shift $h$ and write $P+h=tw$ with $w$ a power of a prime assigned to $A$. Express each $B_i(P)$ as $B_i(tw-h)$. Its coefficient of $w^j$ is $t^j$ times the coefficient of $Z^j$ in $B_i(Z-h)$. The analogous statement holds for the total. For every $j<n$ these signed coefficients have absolute value at most $C t^{n-1}$, with an effective constant $C$ for all rows and nodes. Under (1), $w\ge P^{1-1/(2n)}$, so these absolute values are less than $w/4-1$ for large $P$.

Convert each polynomial evaluation to ordinary base-$w$ digits from the bottom. At indices below $n$ the incoming borrow is 0 or 1. Each strictly negative coefficient produces a digit greater than $3w/4-1$, and each strictly positive coefficient produces a digit less than $w/4$. Two negative coefficients from different rows at the same index would force a carry. A negative coefficient also cannot occur where the total has a strictly positive coefficient, as the sum of the three digits must equal the total digit. This use of base $w$ is valid since $w$ is a power of a prime at which $A$ is carry-free. No bound on the top coefficient after substitution is needed. Constants at the nodes used below are nonnegative, and leading coefficients are positive. Every strict negative counted therefore occurs at an index from 1 to $n-1$.

At the shifted weight choose the largest shift $h$. By (1) at that node, $w>T$ for large $P$, so its base is assigned to $A$. In the variable $Z=X+h$ there are $n-1$ earlier nodes at **positive** real numbers. At each, two part polynomials vanish. By Descartes' rule, a polynomial with positive leading coefficient and $b$ distinct positive roots has at least $\lceil b/2\rceil$ strictly negative coefficients. Hence there are at least $n-1$ negative coefficients across the three parts. The nonconstant coefficients of $dF(Z-h)+1$ alternate strictly in sign, with only $\lfloor n/2\rfloor$ negative positions. All other positions below the top are positive or are the constant position. This is impossible for $n\ge3$. For $n=2$, two different parts vanish at the earlier node, requiring two negative coefficients for one available position. This proves the one-factor version for the successor.

At the unshifted weight, consider the indices $j$ at which the base of $w_j$ is assigned to $A$. Let $b$ be their number. At each such node every $B_i(-h_j)$ vanishes, since it is a bounded integer divisible by $w_j$. Choose the largest such node if $b>0$. Now each part has $b-1$ distinct positive roots in the translated variable. This requires at least $3\lceil(b-1)/2\rceil$ negative coefficients. The total also has at least $b-1$ positive roots and positive leading coefficient, hence at least $1+\lfloor(b-1)/2\rfloor$ strictly positive coefficients, all at positive degrees. Counting disjoint positions in degrees 1 through $n$ gives
$$3\lceil(b-1)/2\rceil+1+\lfloor(b-1)/2\rfloor\le n.$$
Equivalently $b\le B=\lfloor n/2\rfloor+[n\equiv1\pmod4]$. This bound also holds when $b=0$.

Because $U$ is unitary in $N$, its part from each of the remaining $n-b$ factors is at most $t_j\le P^{1/(2n)}$. For $P$ above the largest shift,
$$ U\le d(2P)^b P^{(n-b)/(2n)}.$$
But $27(N-1)^2\le U^3 S^2$ and $S$ is bounded. We need $b+(n-b)/(2n)\ge2n/3$ as $P$ grows. For $n=2,3,4,5$ substitution of $B=1,1,2,3$ gives strict failure. For $n\ge6$, $B\le(n+1)/2$ and $B+(n-B)/(2n)<n/2+1\le 2n/3$. This finishes the proof. For an explicit threshold, take the maxima of the finitely many resultant and node bounds, impose the digit bound above, and solve the last strict power inequality.

## An explicit bound for varying shifts and multiplier

Put $R=2dF(1)$. We have $R\ge d2^n$, $R\ge1+h_{n-1}$ and $R\ge n$. The coefficient $1$-norms of the parts, the total and $dF\pm1$ are at most $R$. Hadamard's inequality on the Sylvester matrix bounds a nonzero resultant of two such polynomials of degrees at most $n$ by $R^{2n}$. At the unshifted weight this gives $S\le R^{3n}$.

At the shifted weight write $H=dF+1=J^e$ as before and take a part $B$ not divisible by $H$. If $J\nmid B$, then $U\mid\operatorname{Res}(H,B)$. Otherwise $e=2$ and $B=JL$ with $J\nmid L$. For every $r^v\parallel J(P)$ with $r\mid U$, we have $r^v\mid L(P)$ since $r^{2v}\mid B(P)$. Therefore $r^v\mid B'(P)$ by the product rule. The resultant $\operatorname{Res}(H,B')$ is nonzero since $J\nmid J'L$ in characteristic zero. Thus $U\mid\operatorname{Res}(H,B')^2$. The derivative has norm at most $nR$, so in all cases $U\le(nR)^{4n}\le R^{8n}$ and $T\le R^{12n}$. Each root evaluation used in the three CRT groups is bounded in absolute value by $R(1+h_{n-1})^n+1\le R^{n+2}$. Failure to obtain a unit vector would force $P\le R^{15n+6}\le R^{18n}$.

Take $P>R^{\max(83,18n)}$. Every $w_j$ satisfying (1) is larger than $R^{12n}$. The absolute signed coefficients at degrees below $n$ in the change of origin are at most $R^{n+1}P^{(n-1)/(2n)}<w_j/8$, as $P^{1/2}>8R^{n+1}$. The integer node and digit arguments therefore apply with the stated threshold, including the successor result with just the last factor condition.

At the unshifted weight set $m=2n-3B-3(n-B)/(2n)$, where $B$ is the upper bound on the number of assigned nodes. The inequalities $N-1\ge P^n/2$, $S\le R^{3n}$ and $U\le R P^{B+(n-B)/(2n)}$ would give $P^m\le(4/27)R^{6n+3}$. For $n=2,3,4,5$, the values of $m$ are $1/4,2,5/4,2/5$, and $(6n+3)/m\le83$. For $n\ge6$, using $B\le(n+1)/2$ gives $m\ge n/2-9/4\ge n/8$ and $(6n+3)/m\le52$. Thus $P>R^{83}$ contradicts the size inequality. This completes the uniform bound.

For example, take any $n\ge3$ distinct full prime powers $w_0<\cdots<w_{n-1}$ of distinct bases. Put $P=w_0$ and $D=w_{n-1}-P$. Both their product and its successor satisfy the conjecture whenever
$$P>\big(2(1+D)^{n-1}\big)^{\max(83,18n)}.$$
The exponents are explicit sufficient bounds. For three factors the direct gap test below gives a stronger bound at the successor.

## Products from a growing block of consecutive primes

Fix $n\ge2$ and $d\ge1$, and write $E_n=\max(83,18n)$. Among starting primes in $[X,2X]$, at most
$$O_{n,d}\big(X^{1-1/((n-1)E_n)}\big)$$
can fail the following assertion. For *any* subset of at least two of the $n$ consecutive primes beginning at the starting prime, both $d$ times its product and the successor of that multiple satisfy the conjecture. In particular this proves infinitude without any hypothesis on prescribed prime tuples.

There is also a simultaneous result with growing parameters. Put $m=\lfloor\tfrac17\sqrt{\log X/\log\log X}\rfloor$. At all but $O(mX/(\log X)^2)$ starting primes in $[X,2X]$, take **any subset of at least two of the next $m$ primes, counting the starting prime**, and **any integer $1\le d\le X^{1/(100m)}$**. Both $d$ times the subset product and one more than that multiple satisfy the full conjecture. The proportion of excluded starts among primes tends to zero.

To prove both claims, the prime number theorem puts the last prime of each block below $3X$ for sufficiently large $X$. Each prime gap contributes to the sum of the block diameters at most $n-1$ times in the fixed case, and at most $m-1$ times in the growing case. The respective diameter sums are at most $2nX$ and $2mX$. For fixed $n,d$ the explicit bound covers every diameter $D$ with $(2d(1+D)^{n-1})^{E_n}<X$. Counting larger diameters gives the first assertion. For growing $m$ discard all diameters exceeding $(\log X)^2$. For any remaining subset, at its smallest prime $P\ge X$ the factors are translates by at most $D$ and each is prime. Its threshold is at most $(2d(1+D)^{m-1})^{18m}$ for large $X$. The logarithm of this bound is at most
$$\big(18/100+36/49+o(1)\big)\log X<\log X.$$
All the indicated choices are consequently covered at once.

## A gap bound with no asymptotic threshold

Let the predecessor have exactly three distinct prime bases, with its full prime powers in increasing order $x$, $x+h$, $x+h+k$. A counterexample requires
$$x\le hk(2h+k).$$
Thus the theorem applies even when the gaps grow with the powers. Together with the [finite three-factor verification](endpoint_factors.md), it gives the full conjecture whenever both adjacent gaps are at most 18, or whenever $h=1$ and $k\le139$, without an upper bound on the weight.

Write $M=N-1$ and use the [three-factor CRT row](endpoint_factors.md). That row avoids all three bases and has parts $A_i=(M/w_i)c_i$ with $3\le U\mid c_i<w_i$, where the $w_i$ are the three powers. Let $D_i=\prod_{j\ne i}(w_j-w_i)$. Since $c_iD_i\equiv1\pmod{w_i}$ and $c_i\ge3$, we have $|c_i D_i|\ge w_i-1$. Hence $A_i\ge M(1-1/x)/|D_i|$. The sum of the two smallest among the reciprocals $1/|D_i|$ is $1/(hk)$.

The trinomial sum is carry-free in the base of $x$, and $N>x^3$. At least one part therefore reaches $x^3$. Combining the other two parts gives $M(1-1/x)/(hk)\le M-x^3+1$. Direct expansion gives $(M-x^3+1)/M<(2h+k)/x$. It follows that $x-1<hk(2h+k)$, proving the integer bound. In each stated fixed gap range the possible $x+h+k$ left by the bound is below 20,000.

## Finite certificates

The explicit cubic and quartic families use the [polynomial specialization verifier](polynomial_prime_powers.md#an-effective-test-including-reducible-polynomials). One extra simplification for the unshifted weight removes the common polynomial divisor of the three normalized parts. Its value divides $N$ and is coprime to the predecessor. If two reduced parts are positive constants $a,b$, every power of $S$ divides at least one of them, so $S\mid\operatorname{lcm}(a,b)$. Delete any base dividing $(a+b)/\gcd(a,b)$, a fixed divisor of the paired binomial coefficient, and any base at which $N(X)-1$ has no root.

The complete run handles 430,530 polynomial rows and 92,088 integer evaluations above the previously verified initial interval. Every such evaluation already fails the necessary numerical size or part-gcd bound. The weight at every smaller input where coefficients could overlap lies within the initial interval. Run the following with python-flint and SymPy installed.

```sh
python3 check_translated_products.py
python3 check_translated_products_cross.py
```
