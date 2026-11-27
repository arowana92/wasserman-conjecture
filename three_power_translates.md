# Adding three powers at unrestricted positions

**Theorem.** Outside $O(X^{0.9})$ positive integers $C\le X$, every weight $N\ge3$ obtained from $C$ by adding at most three powers of two satisfies Wasserman's three-row conjecture. All three exponents are independently arbitrary nonnegative integers. Repetitions and carries into the digits of $C$ are allowed.

Consequently, outside $O(2^{0.982n})$ binary words of length $n$, one may first change any $\lfloor n/100\rfloor$ or fewer of their digits and then add any three powers of two, at positions of any size. All the resulting weights $N\ge3$ work. Indeed each exceptional edited word has at most $2^{H(1/100)n}$ preimages by the binomial theorem, where $H(1/100)<0.081$. Include the edited value zero among the exceptions.

**Prime-base version.** There are absolute $c_0,c_1>0$ such that for $2\le P\le(\log X)^{c_0}$, at most
$$ O\big(X\exp(-c_1\log X/\log(2P))\big) $$
starting values can fail at even one prime $p\le P$, any choice of at most three powers of $p$, and any resulting $N\ge3$ with $p\mid N(N-1)$. This endpoint divisibility condition holds in particular if all added powers are multiples of $p$ and $C\equiv0,1\bmod p$.

At a larger explicit height the [single large divisor argument](three_power_diagonals.md) forces all three added powers to have one of finitely many pairs of exponent differences, for every individual offset.

## Bounded heights

Discard $C=1$. At a counterexample the [endpoint reduction](all_weights_endpoint_reduction.md) supplies a pair of rows $A,B$ and prime-power assignments from $N,N-1$. Choose $A$ to receive $p$. Its products $U\mid N,S\mid N-1$ satisfy
$$
 27(N-1)^2\le U^3 S^2,\qquad U\mid y,\quad S\mid y(y-1),\quad y\ge U\ge3 \tag{1}
$$
for each part $y$ of $A$, by the [complementary row bound](projected_rows.md). In particular $y^5>27(N-1)^2$, since $S\mid(y/U)(y-1)$. Full carry freedom at $p$ permits allocating the added singleton units and the digits of $C$ to the three parts of $A$. If original units overlap, first merge $p$ equal units and then give all constituents of a merged unit to its eventual recipient in the carry-free row.

Put $n_p=\lceil\log_p(2X+2)\rceil$. Fix the list of added exponents of maximum $a\le90n_p+9$, allowing an empty list. There are $O((1+\log X)^3)$ choices for each prime. Every relevant weight has $O(\log X)$ bits uniformly for $p\le X$. Below, a factor $E=\exp(O(\log X/\log\log X))$ covers polynomial logarithmic losses and any fixed number of divisor counts on integers of $O(\log X)$ bits. This follows, for example, by splitting prime factors at $\sqrt{\log X}$.

For $K=US\le Y\le X$ the [complementary-row reconstruction](two_predecessor_congruences.md#recovering-the-row-from-a-product) gives $O(YE)$ possible weights, by summing its two dyadic bounds and counting weights below $Y$ directly. Each weight fixes $C$ at this list of exponents. For $K>Y$, the low contribution to one part uses no more than a third of the base-$p$ digit sum of $C$. It has at most eight choices of high contribution. The part $y$ allows at most $E$ labels $(U,S)$ by (1), and the label fixes $N\bmod K$, hence admits at most $1+X/Y$ constants.

At $p=2$ discard the $O(X^{H(7/10)})$ inputs with more than $7n_2/10$ ones. In the other cases the low contribution has at most $7n_2/30$ ones and $O(X^{H(7/30)})$ values. Since $H(7/10)<0.89$ and $H(7/30)<0.79$, take $Y=X^{0.895}$. This gives the claimed bound after the $E$ loss.

For all primes put $\rho=3\sqrt[3]2/4<1$. The number of $n_p$-digit strings of digit sum at most $(p-1)n_p/3$ is at most $p^{n_p}\rho^{n_p}$. For a uniform digit $d$ the convex chord bound gives $\mathbb E[2^{-d/(p-1)}]\le3/4$. Multiply over the digits and apply Markov's inequality. Thus take
$$ Y=X\exp\big((\log\rho)\log X/(2\log p)\big). $$
Since $p^{n_p}\le p(2X+2)$, both ranges cost at most $O(pYE)$ for this prime. Summing over $p\le P$ and choosing $c_0$ sufficiently small absorbs the divisor losses into the displayed saving of the prime-base theorem.

## What three very high powers force

Now $a>90n_p+9$. With at most two added powers, some part is at most $C$, contradicting (1). With three, write their exponents as $a_1=a\ge a_2\ge a_3$. Each of the three parts must now receive exactly one of them, so
$$ A_i=p^{a_i}+c_i,\qquad c_1+c_2+c_3=C $$
with the three low terms adding freely at $p$. Moreover $a_3>4n_p$, since otherwise one part is at most $p^{4n_p}+X$, whose fifth power is less than $27(N-1)^2$.

Split $S=S_1S_2S_3$ according to the part having residue one. By (1), some factor among $U,S_1,S_2,S_3$ is at least $(27(N-1)^2)^{1/9}$. Its part $W$ coprime to $p$ satisfies $W\ge p^{2a/9-n_p}$, since $N\equiv C\bmod p^{a_3}$ and $2\le C<p^{a_3}$. Let $e$ be zero for $U$ and the corresponding coordinate unit for $S_i$. Put $x_i=e_i-c_i$. Then
$$ W\mid p^{a_i}-x_i,\qquad 1\le |x_i|\le X.$$
A zero $x_i$ is impossible because $W>1$ is coprime to $p$.

For each $i=2,3$ Dirichlet approximation gives $1\le u_i\le9$, $0\le v_i\le u_i$ with $|r_i|=|v_i a-u_i a_i|\le a/9$. Compare $x_1^{v_i}$ and $p^{r_i}x_i^{u_i}$ after clearing the power of $p$ if $r_i<0$. The difference is divisible by $W$ and of absolute value at most $2p^{9n_p+a/9}<W$. Therefore
$$ |x_1|^{v_i}=p^{r_i}|x_i|^{u_i}. \tag{2}$$
If $v_i=0$ then $|x_i|=p^{a_i}>X$, so both ratios $v_i/u_i$ are positive. Let $d$ be the lcm of their reduced denominators. At primes different from $p$, (2) shows that every valuation of $|x_1|$ is a multiple of $d$. Hence there exists a positive integer $t$ prime to $p$ such that
$$
 |x_1|=p^{j_1}t^d,\qquad |x_i|=p^{j_i}t^{d v_i/u_i}\quad(i=2,3),\qquad 0\le j_i\le n_p. \tag{3}
$$
In particular, the three parts of $|c_i-e_i|$ after removing their $p$-powers are all positive powers of the same integer. This conclusion is an effective necessary condition at each individual $C$.

If $d\ge2$, then $t\le\sqrt X$. The bounded indices, signs and $O((1+\log X)^3)$ valuation choices show that only $O(\sqrt X(1+\log X)^3)$ constants arise. If $t=1$, the same count holds directly. Suppose $d=1$ and $t>1$. Now both ratios equal one. Every $x_i$ is negative, as $x_i\le1$ and $|x_i|\ge t$. By (2), the three low parts have the form
$$ (c_1,c_2,c_3)=l(p^{a-a_3},p^{a_2-a_3},1)+e,\qquad l>0. \tag{4}$$
Both exponent differences are at most $n_p$. Put $s=p^{a-a_3}+p^{a_2-a_3}+1$ so $l\le X/s$.

At $p=2$ the [affine-triple prefix count](binary_affine_counts.md) gives $O(X^{0.551464})$ possible $l$ at each pair of differences: the input needs $\log_2(1+X/s)+O(1)$ bits and its slope delay spread is at most $\log_2 s$. Thus for heights $a>90n_2+9$ alone the total is $O(X^{0.551464}(1+\log X)^3)$.

For a bound uniform in $p$, at $s\ge X^{1/4}$ count all $l\le X^{3/4}$. Otherwise take two unshifted coordinates in (4). Their exponents differ by some $g\le\log_p s$. Among the $k=\lceil\log_p(X+1)\rceil$ possible digits of $l$, two at distance $g$ must sum to at most $p-1$. For $g>0$ the $k-g$ pairs form a graph of maximum degree two, so a greedy matching has at least $(k-g)/3\ge k/4$ edges. Each edge allows the fraction $(p+1)/(2p)$ of ordered digit pairs. At $g=0$ each individual digit must be at most $(p-1)/2$, which gives at least as much saving. Thus there are at most $p(X+1)(3/4)^{k/4}$ scales. The valuation choices cost only powers of $1+\log X$. This gives an absolute exponential saving in $\log X/\log p$, as required.

## Checking individual binary offsets

For each $C$ the verifier checks the necessary condition after (3). It enumerates numbers of the form $2^j t^e$ with $e\ge1$, including $t=1$, and tests whether three such numbers, with the permitted signs and one optional unit adjustment, give a carry-free decomposition of $C$. A negative unadjusted low term can only be $-1$ and must receive the unit. For $t>1$ all three magnitudes are at least $t$, so $t\le C/3$. The independent check partitions the actual ones of $C$ among the three parts and compares the primitive perfect-power bases of their odd parts.

~~~sh
python3 check_three_power_translates.py
~~~

Among $2\le C\le1,048,576$, exactly 82,740 offsets pass the necessary test and 965,835 are excluded. Through $C=4096$ both generation methods agree on all 1,899 passing and 2,196 excluded offsets. The cutoff and the decomposition criterion use each individual $C$, independently of the upper limit used for tabulation.

Passing this finite test is only a necessary condition for a large-height exception. An offset excluded by it works for every added exponent triple whose maximum exceeds $90\lceil\log_2(2C+2)\rceil+9$.
