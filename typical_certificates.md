# Short prime certificates at typical weights

**Theorem.** There is an explicit prime list $\mathcal P(N)$ which certifies the full three-row conjecture for a set of weights of natural density one, and
$$\limsup_{X\to\infty}\frac1X\sum_{3\le N\le X}|\mathcal P(N)|\le16+6\log25<35.32.$$
Here a list certifies $N$ if every selection has a common prime divisor in the list. The primes may depend on $N$. In particular, the set of weights with a certificate of at most **70 primes** has lower natural density greater than $1/2$. For every $\epsilon>0$ a fixed list-size bound suffices at lower density at least $1-\epsilon$.

## Constructing the list

Put $\delta=1/25$. For $n>1$, order its full prime powers by increasing prime base. Discard the longest initial segment with product at most $n^\delta$ and retain the bases of the remaining powers. Use an empty list for $n\le1$. Take the union of the retained lists for $N,N-1,\ldots,N-5$, together with $2,3,5$ and the largest prime at most $N$. Primes exceeding $N$ may be dropped.

Suppose three rows collectively avoid every prime on this list. We count the weights where this could happen in $[X/2,X)$.

## When all second parts are at least 5

Write each row in increasing order. For almost every weight in the block the largest prime $r\le N$ satisfies $g=N-r\le N^{0.14}$. For example, apply Guth–Maynard’s almost-all theorem at interval length $\lfloor X^{0.138}\rfloor$, shifting the interval start to the left of $N$. The lower boundary costs only the interval length. [Guth–Maynard, Corollary 1.4](https://arxiv.org/pdf/2405.20552)

The row avoiding $r$ has its first two parts at most $g$. Use height 6 in the [weighted lemma](weighted_obstruction.md). On each ordinary row choose $(2,2,2)$ if its minimum is at least 2, or $(1,2,3)$ if it is 1. Their combined cost is at most $N^{18}$. On the near-prime row choose $(2,4,0)$ or $(1,5,0)$ respectively. Its cost is at most $g^{15}$, as its second part is at least 5. For each offset the full product of discarded powers is at most $N^\delta$. After removing them and the layers through 6, every remaining prime is on the list. The proof of the weighted lemma therefore gives a lower bound $\gg N^{21(1-\delta)}$ for the three-row product. This is impossible for large $N$ since $21(1-\delta)=20.16>18+15(0.14)=20.1$.

## When a row has two parts at most 4

Call this row $R=(a,b,N-a-b)$. For $1\le a\le b\le4$, at least one $p\in\{2,3,5\}$ divides $\binom{a+b}{a}$, as is seen from the ten values $2,3,4,5,6,10,15,20,35,70$. This $p$ divides the coefficient of $R$, so one other row $B$ must be carry-free in base $p$.

We may assume that for each of $N$ and $N-1$, the product of its powers of $2,3,5$ is at most its $\delta$-th power. Violations have density zero, by counting multiples of powers of one of these three primes exceeding $(X/3)^{\delta/3}$. Under this assumption all retained primes from $N$ and $N-1$ exceed 5. The row $R$ avoids none of them: its residues $a,b$ are positive and already sum to more than 1.

Assign the retained full powers of $N$ and $N-1$ to the other two rows which avoid their primes. Write $U,S$ for the products assigned to $B$ and $V,T$ for the complementary products. With $W,J$ the two discarded products,
$$N=UVW,\qquad N-1=STJ,\qquad W,J\le X^\delta.$$
The third row has parts $Vy_i$ summing to $N$, and $T^2\mid y_1y_2y_3$. Set $K=US$.

If $K<X^{71/100}$, apply the [small-modulus counting lemma](endpoint_count.md) with the parameters $U'=UW$, $S'=SJ$. That lemma uses only $U'V=N$, $S'T=N-1$ and the complementary row’s exceptional bins for $T$. These hold because each assigned power in $T$ divides two normalized parts and has its own exceptional position. The positivity used in its divisor identities holds as well, since $V=N/(UW)\gg X^{1-0.71-\delta}$ in this case. The lemma bounds the count by $X^{-1+o(1)}(X^{71/100+2\delta})^{5/2}=X^{0.975+o(1)}$.

For the remaining weights suppose first that every part of $B$ is at least 2. Write $L_p=\lceil\log_p X\rceil$ and $\mu_p=(p-1)/2+1/100$. For almost every weight, $s_p(N)\le\mu_p L_p$ simultaneously for $p=2,3,5$, where $s_p$ is the base-$p$ digit sum. Indeed, over the complete range of $L_p$-digit strings the uniform digits are independent with mean $(p-1)/2$ and bounded variance. Chebyshev’s inequality already suffices.

Since $B$ is carry-free, it has a part $z$ of digit sum at most $\mu_p L_p/3$. With $F_p(x)=1+x+\cdots+x^{p-1}$, the number of possible $z$ is at most $(F_p(1/3)3^{\mu_p/3})^{L_p}$. For each of the three primes this is $O(X^{c_p})$ for a fixed $c_p<0.71$. The integer-power comparison is $F_p(1/3)^{300}3^{50(p-1)+1}<p^{213}$. Because $z\ge2$, the relations $U\mid z$ and $S\mid z(z-1)$ leave only $X^{o(1)}$ pairs $(U,S)$. Each specifies a residue class of $N$ modulo $K\ge X^{0.71}$. Thus there are $o(X)$ weights of this type.

Finally suppose $B$ contains a 1. Then $U=1$ and the normalized third row sums to $W$. We obtain $T\le(W/3)^{3/2}$ and $S\gg X^{1-5\delta/2}=X^{0.9}$, so both of the other parts of $B$ are positive multiples of $S$. One of them has digit sum at most $\mu_p L_p/2$. There are $O(X^{0.86})$ possibilities, as $F_p(1/2)^{200}2^{50(p-1)+1}<p^{172}$ for each $p$. Each such part has $X^{o(1)}$ possible divisors $S$, and $N\equiv1\pmod S$. This case contributes at most $X^{0.96+o(1)}$ weights. The block bounds establish the asserted density one.

## Mean list size

For each retained prime $p$ of $n$, the product of the powers of primes $q\le p$ dividing $n$ exceeds $n^\delta$. In averaging over $n\le X$, we may drop the initial range $n<X/(\log X)^2$ for an $o(1)$ contribution to the mean, since a number has at most $\log_2 X$ distinct prime factors. Set $Y=(X/(\log X)^2)^\delta$.

The expected count of primes $p>Y$ dividing $n$ is at most $\sum_{Y<p\le X}1/p=\log(1/\delta)+o(1)$. For $p\le Y$, upper-bound the retention indicator on multiples of $p$ by $\sum_{q\le p}v_q(n)\log q/\log Y$. Averaging by counting multiples of $p$ and $q^e$ gives the upper bound
$$\frac1{\log Y}\sum_{p\le Y}\left(\frac1p\sum_{q<p}\frac{\log q}{q-1}+\frac{\log p}{p-1}\right)=2+o(1),$$
where $p,q$ range over primes. Here we used $\sum_{p\le y}\log p/(p-1)=\log y+O(1)$ and partial summation. The former follows from $\log(\lfloor y\rfloor!)=\sum_{d\le y}\Lambda(d)\lfloor y/d\rfloor$, Stirling’s formula and Chebyshev’s bound $\sum_{d\le y}\Lambda(d)=O(y)$, since the higher prime-power tails contribute $O(1)$ after division by $y$.

Six shifts and the four additional primes give $4+6(2+\log25)$ as claimed. Markov’s inequality bounds the upper density of list sizes at least $k+1$ by $(16+6\log25)/(k+1)$. The list fails to certify at only a density-zero set of weights.

## Reproducing the finite check

Use `python3 check_typical_certificates.py N` to generate the proposed list at a single weight. At an individual untested weight it is a candidate list. Without arguments the program checks its lists for $3\le N\le5000$. All 4,998 lists certified their weights, with mean size $56{,}147/4{,}998$ and maximum size 16. For 103 weights the program checked a different set from the one already supplied in the [minimum-size table](certificates_5000.txt). In a separate verification, sublists of at most five primes from these candidate lists passed the C++ factorial-valuation check on all 3,473,263,472 unordered rows through 5,000. The Python program also checks the six exact digit inequalities used above.
