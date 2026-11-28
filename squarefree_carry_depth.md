# Every carry depth from three onward occurs at squarefree endpoints

**Theorem.** Fix any integer $H\ge2$. For $\gg_H X/(\log X)^2$ weights $N\le X$, the integer $N(N-1)$ is squarefree and the **least depth sufficient for every triple of rows is exactly $H+1$**, allowing all primes. At each such weight every pair of rows has a common endpoint prime by that same depth, drawn from an explicit list of at most six primes greater than ten. Here an endpoint prime divides $N(N-1)$. There is also an absolute constant $K_*$, independent of $H$, such that a list of at most $K_*$ primes attains this sharp depth for all triples at $\gg_H X/(\log X)^2$ weights. The six-prime bound applies to the pair list. The prescribed hard triple includes $(1,1,N-2)$.

The units digit is digit one. Each endpoint prime occurs just once in $N(N-1)$, so even inspecting every full endpoint power and any fixed number of further digits fails at some weights. This strengthens the [unbounded absolute-depth theorem](unbounded_absolute_carry_depth.md). The construction also gives a criterion proving Wasserman's conjecture when a cofactor of $N$ has arbitrarily many prime divisors.

## A row forced by three exceptional residues

Choose a prime $u>10$, to be taken sufficiently large below. Dirichlet's theorem gives primes $P>u$ and then $Q>2P+1$ satisfying
$$
 \begin{aligned}
 P&\equiv1\pmod{u^H}, &P&\equiv 1+2^{H-1}\pmod{2^H},\\
 Q&\equiv-1\pmod{u^H P^H}, &Q&\equiv 2^H-1\pmod{2^{H+1}}.
 \end{aligned}
$$
Write $a=2PQ$, $N=1+aq=uV$, and require
$$
 q\equiv1\pmod2,\qquad 2q\equiv1\pmod{P^H},\qquad
 q\equiv1+Q^H\pmod{Q^{H+1}},\qquad 1+aq\equiv u\pmod{u^2}.       \tag{1}
$$
This specifies one class $q_0\bmod M$, $M=2P^H Q^{H+1}u^2$. Suppose $q>a$ is prime, $V>1$, and all prime divisors of $V$ exceed $a$. Put
$$
 A=((Q+1)q,\ 2(P-1)Qq,\ (Q-1)q+1).
$$
Its parts sum to $N$. At 2 and $u$ the first two parts vanish through depth $H$, so the addition passes there. At $P^H$ the first and third parts vanish, by $Q\equiv-1$ and $2q\equiv1$. These tests at a power imply the residue identities at every smaller power in the same base. At $Q$ the first $H$ digits agree with the seed $q=1$, where the parts are $Q+1,2(P-1)Q,Q$. These add without a carry, with sum less than $Q^2$. The added $Q^H$ in (1) makes the next column of $A$ have digits $1,0,Q-1$. Its first carry at $Q$ occurs at digit $H+1$. At $q$ the upper counts add to $a<q$, so $A$ avoids $q$ completely.

Consider the list
$$
                 \{q,u,P,Q\}\ \cup\ \{s:s\text{ prime},\ s\mid V\}.       \tag{2}
$$
For two rows, test the units digit at $q,u,P$, the first $H+1$ digits at $Q$, and the full weight power at each base dividing $V$. If a pair could pass these tests between them, a row passing $q$ would have two positive parts $i q,j q$ with $i,j\le a$. It fails the units digit at each $s\mid V$, since $N\equiv0\bmod s$ and $0<i<s$. The companion must therefore have all parts divisible by $V$. Its three positive multipliers sum to $u$. It fails at $u$ since $u\nmid V$ and each multiplier is less than $u$. It also fails at $P,Q$, since all its residues are nonzero and their sum would have to be one.

The $q$-row would thus pass the units digits at $q,P,Q,u$. Each of $q,P,Q$ gives one exceptional position with residue one and two positions with residue zero. All parts are multiples of $u$. Every position must receive an exception, since a position receiving none would be a positive multiple of $uPQq=u(N-1)/2>N$. With three exceptions on three positions, their assignment is unique up to permutation. CRT determines each part modulo $uPQq>N$. The row $A$ has these residues, so it is the unique candidate. Its carry at $Q$ contradicts the test. This proves the assertion for pairs and proves the ordinary conjecture at every weight satisfying the hypotheses by the [all-weight endpoint reduction](all_weights_endpoint_reduction.md).

There are infinitely many such weights with $V$ having arbitrarily many distinct prime factors. For each prime $\ell\le a$ outside $M$, avoid the at most two forbidden classes where $q$ or $V$ vanishes modulo $\ell$. At the primes of $M$ both are already units. Choose any finite set of additional primes $\ell>a$ and require $1+aq\equiv0\bmod\ell$ there. CRT and Dirichlet give a reduced progression containing infinitely many primes $q$. The depth assertion for (2) uses just the units digit at the factors of $V$ when $V$ is squarefree.

## A hard triple at these weights

Let $V$ be prime or $V=st$ with large distinct primes $s<t$. Set $E=(1,1,N-2)$ and
$$
                C=(kV,kV,(u-2k)V),\qquad 1\le k<u/2 .
$$
When $C$ avoids every prime dividing $V$, no prime witnesses $E,A,C$ through digit $H$. In fact the endpoint restricts the primes to $2,u,P,Q,q$ and the factors of $V$. The tests above cover every one. At $Q$ both $E$ and $C$ carry in their units digits, since the three positive multipliers of $C$ are less than $Q$ and $N\equiv1\bmod Q$. Thus the **first** common depth is $H+1$. The argument with (2) proves the matching upper bound for every pair with the endpoint. The factors of $N(N-1)$ are all distinct.

## Counting enough prime and almost-prime values

For clarity, we use the following unweighted consequence of the two-linear-form Chen sieve. For distinct nonproportional admissible fixed positive linear forms $L_1,L_2$ there are $\gg \mathfrak S(L_1,L_2)x/(\log x)^2$ indices $x\le n<2x$ with $L_1(n)$ prime and $L_2(n)$ having at most two prime factors, all at least $x^{1/10}$. The main-term constant is absolute, and the threshold may depend on the coefficients.

Here is the specialization of the quantitative sieve needed for this assertion in our application. Take
$$
 L_1(n)=q_0+Mn,\qquad L_2(n)=V_0+Kn,\qquad
 V_0=(1+aq_0)/u,\quad K=aM/u ,
$$
with $0<q_0<M$. The two slopes and their determinant $MV_0-Kq_0=M/u$ have prime support $\{2,u,P,Q\}$. Both intercepts are units at those primes. At each other prime there are two distinct forbidden residues. Thus the singular series is bounded below by an absolute positive constant (compare the convergent twin-prime product over primes above two).

Apply Matomäki and Shao's Theorem 6.4 with weight identically one. In their Hypothesis 6.3(1) the allowed moduli $d$ are coprime to the two slopes and determinant. The restrictions on $L_1$ give a single reduced class modulo $Md$. The usual Bombieri–Vinogradov theorem gives the required averaged error, since $M$ is fixed and $d\le x^{1/2-\eta}$ for a fixed $\eta>0$. The extra convolution weights with prime factors at least $x^{1/10}$ are bounded absolutely, since at most five such primes divide $d$. [Matomäki and Shao, Hypothesis 6.3 and Theorem 6.4](https://arxiv.org/pdf/1512.03213) · [Bombieri, *On the large sieve*](https://doi.org/10.1112/S0025579300005313)

For Hypothesis 6.3(2) the restrictions on their three-prime products $L_2(n)\in B_j$ specify a single reduced class modulo $Kd$. Use the unweighted type II Bombieri–Vinogradov estimate for the $B_j$ indicators, as in Bienvenu, Shao and Teräväinen's proof of (8-8), which treats the constant term of the weight separately. The product interval here is $[L_2(x),L_2(2x))$ with fixed slope, obtained by using the two endpoints. The decomposition into restricted prime factors uses the same size conditions, independently of the reduced residue being counted. Multiplying the modulus by the fixed $K$ stays below the level of distribution. The target main term is the count in the class modulo $K$ divided by $\varphi(d)$. Removing products sharing a prime with $d$ from that main term costs $x^{9/10+o(1)}$ in the sum over $d$ because each product has at most three prime divisors $p\ge x^{1/10}$, and $\sum_{p\mid d,\ d\le\sqrt x}1/\varphi(d)\ll (\log x)^2/p$. The fixed-modulus estimate at $K$ also gives Hypothesis 6.3(3), since all these products are coprime to $K$ for large $x$. This verifies the hypotheses for our forms. [Bienvenu, Shao and Teräväinen, §8B, unweighted estimates for (8-7), (8-8) and (8-11)](https://pure.mpg.de/rest/items/item_3509190_4/component/file_3509192/content)

Put $Y=V_0+2Kx$. We have supplied at least $c_0Y/(K\log^2 Y)$ possibilities with $Y/2\le V\le Y$, where $c_0>0$ is absolute once $Y$ is sufficiently large for the fixed forms. The choice $k=1$ works at a prime $V>u$. In the composite case we adapt the [digit-box and upper sieve proof](unbounded_absolute_carry_depth.md) as follows.

* For $s\le t$ we have $Y^{19/200}<s\le\sqrt Y$ for large $Y$. Remove the boundary ranges $Y^{1-\delta}\le s^i\le Y^{1+\delta}$ for $2\le i\le10$, with $0<\delta<1/200$. Their reciprocal-prime sum is $O(\delta)+o(1)$. Else $s^{J+1}<Y^{1-\delta}$ and $s^{J+2}>Y^{1+\delta}$ for some $1\le J\le9$, since $s^{11}>Y^{1.045}$. We have $t>us$ and $ut<s^{J+1}$. Thus $C$ avoids $t$ for any allowed $k$. It avoids $s$ precisely when $2\{kt/s^j\}\le\{ut/s^j\}$ for $1\le j\le J$.
* Choose ternary blocks of width $w$ and number $L$ as in the linked proof, now allowing $J\le9$. Choose the prime $u>2\cdot3^{wL}$, so all block multipliers $k=3^{wi}$, $0\le i<L$, are available. The exceptional volume is at most $(1-3^{-9w})^L+18\cdot3^{-w}$. The final marginal condition $\{u x_j\}\ge2\cdot3^{-w}$ uses only that $u$ is a positive integer. A union bound suffices even if this condition depends on the digit blocks. All exceptional regions are finite unions of boxes for fixed parameters.
* For each $s$, integrality and the progression for $q$ specify $t\equiv V_0/s\bmod K$, a single reduced class. Sieve $t$ and $ust-1$ at primes through $R$ omitting the prime factors of $6Ku$. At each retained prime below $s$ there are two roots. At most five primes were omitted, uniformly in the construction. Thus (5) of the linked proof applies with $a$ there replaced by $K$, with absolute main-term constants. The mass is $\rho Y/(Ks)$, with $\rho=1$ on boundary primes and $\rho\le\epsilon+o(1)$ on the exceptional boxes if the displayed volume is at most $\epsilon$. The same levels $R=Y^{1/20},Y^{\delta/8}$ on the respective ranges give a discarded count at most $O(\delta+\epsilon\delta^{-2}+o(1))Y/(K\log^2Y)$. In the first range $R<s$ also holds since $s>Y^{19/200}$.

First fix $\delta$ and $\epsilon$ small compared with the absolute lower-bound and sieve constants. Then choose $w,L,u,P,Q$. At least half the lower bound survives for sufficiently large $Y$. Squares $V=s^2$ occur only in the discarded boundary range since $V\ge Y/2$. We have $\gg_H X/(\log X)^2$ suitable weights through $X$ by taking $uY\le X$ with $uY\sim X$.

## Transferring the sharp bound to all triples

Outside $o(X/(\log X)^2)$ weights in a dyadic block $[X,2X)$, any three rows each having its two smallest parts summing to more than ten have a common **units-digit** carry. We include the error estimate because a density-one assertion alone would not suffice on the preceding set.

For $n\asymp X$ put $R(n)=n/\prod_{p\mid n,\ p>10}p$. Its square-root moment has bounded mean. Indeed $f(n)=\sqrt{R(n)}$ is multiplicative, with $f(p^e)=p^{e/2}$ for $p\le10$ and $p^{(e-1)/2}$ otherwise, $e\ge1$. Expand $f(n)=\sum_{d\mid n}(f*\mu)(d)$, whose terms are nonnegative by monotonicity on prime powers. Summing multiples bounds the mean by the Euler product with local sum $1+\sum_{e\ge1}(f(p^e)-f(p^{e-1}))/p^e$. The factors at $p>10$ are $1+O(p^{-3/2})$. Thus the six values $R(N-j)$, $0\le j\le5$, are each at most $X^{1/200}$ outside $O(X^{1-1/400})$ weights. The largest prime below or equal to $N$ has gap $g\le X^{0.18}$ outside $o(X/(\log X)^2)$ weights by the short-interval estimate used in the [bounded-depth proof](bounded_carry_depth.md#three-arbitrary-rows). [Guth and Maynard, Corollary 1.4](https://arxiv.org/abs/2405.20552)

With every two-smallest sum above ten, the second part in increasing order is at least six. The same height-six units-digit test in that proof gives a lower bound $\gg X^{21(1-1/200)}$ on the product of row costs if the three rows pass all the tested primes between them. Its upper bound is $O(X^{18}g^{15})=O(X^{20.7})$. Since $21(199/200)>20.7$, this is impossible for large $X$.

For any triple containing a row with two parts summing to at most ten, apply the pair theorem to the other two rows. Every prime on (2) exceeds ten. The small parts have positive residues of sum greater than one, whereas $N$ has units residue zero or one. The small row therefore carries at the units digit of that same prime. Intersecting with the preceding sieve family now proves the theorem. Dyadic summation preserves the little-oh exceptional count through $X$.

## A bounded number of primes also suffices for the full triple test

We can keep a positive fraction of the preceding family and use an absolute number of primes for all triples. Work in the sieve interval $x\le n<2x$ above, with $N=N_0+uK n\asymp_H x$. At this stage at least $c_1x/\log^2 x$ indices survive for an **absolute** $c_1>0$. Call these the retained indices. Include the pair list (2), the greatest prime at most $N$, and the primes above $x^\alpha$ dividing $N,N-1,\ldots,N-5$ for a sufficiently small absolute $\alpha>0$.

For $0\le j\le5$ we bound the product of the discarded small prime powers of $N-j$ on the retained indices. The at most five primes dividing $6uK$ contribute product at most $x^{1/1000}$, after discarding $o_H(x/\log^2 x)$ indices. Indeed for each fixed prime discard multiples of the first power above $x^{1/5000}$, counting all integers $N$ in the weight interval. We may also discard all $N$ with any $p^e\mid N-j$, $e\ge2$, $p^e>x^{1/4}$. There are $O(\sqrt T\log T)$ powers with exponent at least two through $T$, counting repetitions. Summing $O_H(x/p^e+1)$ in dyadic blocks gives $O_H(x^{7/8}\log^2x)$ discarded weights. Finally the weights with $p^2\mid N-j$ at any $p>x^\alpha$ cost $O_H(x^{1-\alpha}+\sqrt x)$ by the same multiple count.

Fix $d=p^e\le x^{1/4}$ with $p\nmid6uK$. Imposing $d\mid N-j$ specifies one class of $n\bmod d$. In that class sieve $L_1(n)L_2(n)$ at the primes through $x^{1/30}$ omitting primes dividing $6uKp$. The retained indices have no prime factor of either form in the sieving range. At every tested prime there are exactly two forbidden residues by the determinant computation. The Selberg estimate used above now has main mass $x/d$ and unit period for the digit condition. It bounds the count of retained indices in the class by
$$
               C\frac{x}{d\log^2 x}+O(x^{1/15+o(1)})
               \ll \frac{x}{d\log^2 x}.
$$
The constants are absolute since at most six primes are omitted. The error is absorbed uniformly for $d\le x^{1/4}$.

On summing with weight $\log p$ for $p\le x^\alpha$, the total log contribution at these primes is $O((\alpha\log x+1)x/\log^2 x)$ per shift, by $\sum_{p\le z}\sum_{e\ge1}(\log p)/p^e=\log z+O(1)$. Markov's inequality at $(\log x)/300$ for each shift costs at most $O(\alpha+o(1))x/\log^2 x$ in total. Choose $\alpha<1/100$ small enough compared with $c_1$ and the absolute main-term constant. On at least half the retained indices, the omitted products are at most $x^{1/200}$, since $1/1000+1/300<1/200$. Every retained large prime occurs to exponent one in its shift after the square divisibility exclusion. The same height-six inequality and the same prime-gap estimate now prove that this list suffices for all triples at the sharp depth.

Each shift has at most $2/\alpha$ distinct primes above $x^\alpha$ once $x$ is sufficiently large for $H$. Thus one can take $K_*=7+\lceil12/\alpha\rceil$. This proof supplies a uniform bound without optimizing its numerical value.

## Exact example

With $H=2$, $u=13$, $P=3719$, $Q=28049188907$ and $q=q_0+168324M$, the weight is $3622879115469537062805584912986808464521482502471662496288083$. The checker constructs this 61-digit example, certifies its prime factors and the nearby prime $N-6$, independently reconstructs the unique row by CRT, and checks the first carry at each possible common prime. At this weight a six-prime list certifies all triples by digit three. Any row passing the nearby prime has two parts with sum at most six and therefore carries at every prime on (2) in its units digit. The pair theorem applies to the other two rows. [Checker](check_squarefree_carry_depth.py) · [Lucas primality certificates](squarefree_depth_primes.json)

~~~sh
python3 check_squarefree_carry_depth.py
~~~
