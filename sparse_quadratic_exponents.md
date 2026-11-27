# Prime-power exponents at the twelve remaining quadratic shapes

Let $P$ be a prime power, $N=P(3P+2)$ and $3P-1=DQ$ with $D\le2400$ and $Q$ a prime power. Then the full three-row conjecture holds under **any** of the following conditions.

* $P\le10^{1000}$.
* The prime base of either $P$ or $Q$ is at most **one million**, with no bound on its exponent.
* Either power is a square, as proved in the [quadratic cofactor theorem](sparse_quadratic_powers.md).

For each fixed $k\ge2$, only finitely many inputs in this family can fail with either $P$ or $Q$ a $k$th power, even with unrestricted prime bases. With no further hypotheses on either power, the number of possible exceptions with $P\le X$ is $O((\log X)^{1-\delta})$ for some effective $\delta>0$.

## A finite set of leading binary values

The [cofactor theorem](sparse_quadratic_powers.md) leaves twelve triples $(D,i,j)$. We may assume $Q$ odd and full, $P>200{,}000{,}000$, and the row with its first two parts $it,jt$, where $t=Q(DQ+3)$, must have odd coefficient. Write $n$ for the binary length of $N$, $h=\lfloor\log_2Q\rfloor$, and $s=(3i,3j,D-3(i+j))$. The binary prefix discrepancy bound in the cofactor proof holds for every one of the first $h$ columns, since $2^h\le Q$.

For each shape, construct the entire discrepancy graph with its distinct initial state requiring a nonzero digit. It has at most 304 vertices. Repeatedly delete vertices with no outgoing edge to the remaining set. In the resulting graph, **every vertex on a cycle has exactly one outgoing edge**. The checker enumerates paths until their first repeated vertex and verifies this assertion at each vertex of the cycle. Starting from the root this also checks all reachable cyclic continuations, since any first cycle entered cannot be left in this graph. [Graph and interval verifier](check_sparse_quadratic_prefixes.py)

Thus there are finitely many infinite paths, and their binary total digits are eventually periodic. Let $\lambda$ be the number in $[1/2,1]$ represented by such a path. A finite walk can spend at most 304 final steps outside the retained graph. Its remaining prefix agrees with some infinite path. In particular each possible weight has
$$
          |N/2^n-\lambda|\le 2^{305-h}.
$$
The finite list can include $\lambda=1$, represented by repeating ones. Since $(DQ+2)^2=3N+1$ and $N\ge2^{n-1}$, we have $Q\ge2^{n/2}/(2D)$ at the scales in question. Consequently $N=\lambda2^n+O(2^{n/2})$ and each $n$ leaves only $O(1)$ possible values of $P$ or $Q$.

For $z=P$ set $c=3/\lambda$, and for $z=Q$ set $c=D^2/(3\lambda)$. There are 43 distinct positive rationals $c$ over the twelve shapes. Every reduced numerator and denominator has at most 42 bits. In either case
$$
          |cz^2\,2^{-n}-1|<2^{340-n/2}.                         \tag{1}
$$
Indeed $2^{-h}\le2/Q$, $D\le2400$ and $\lambda\ge1/2$. The differences of $N$ from $3P^2$ and from $D^2Q^2/3$ are $2P$ and $(4DQ+3)/3$, both below $4\sqrt N$.

## Finite intervals and unbounded prime bases at fixed exponents

For the numerical bound, propagate every binary prefix $a$ of the total through the first $h$ columns. Restrict $2^h\le Q<2^{h+1}$ and $200{,}000{,}000<P\le10^{1000}$. For each possible binary length $n$ solve the following integer inequalities by square roots.
$$
       3a2^{n-h}< (DQ+2)^2 \le 3(a+1)2^{n-h}-2 .
$$
There are 156,259 candidate scales counted across the shapes. Of these 25,995 have $Q$ odd and $P$ integral. Every one has a binary carry in the required row. Three direct parity controls outside the residual list are recovered. An independent check splits the input $P$ into binary intervals. It tests each total prefix $a$ by directly enumerating integers $v_\ell$ in $[s_\ell(a+1)/D-1,s_\ell(a+1)/D]$ with disjoint bits and sum $a$. This permits different part-prefix witnesses at different depths. It reconstructs 120,253 integral inputs counted across the shapes and rejects every odd-$Q$ case by exact binary digit sums. [Independent verifier](check_sparse_quadratic_prefixes_cross.py)

For the qualitative fixed-exponent result, squares were already excluded. If $z=u^k$ with fixed $k\ge3$, write $n=2km+r$ with $0\le r<2k$. Taking a $2k$th root in (1) gives
$$
 \left|u/2^m-(2^r/c)^{1/(2k)}\right|=O_k(2^{-km}).
$$
The fraction on the left has odd numerator. If the fixed algebraic number is irrational, Roth's theorem gives finitely many solutions since $k>2$. If it is rational of denominator $b$, every distinct fraction has distance at least $1/(b2^m)$, and equality of the fractions forces $2^m\mid b$. There are finitely many choices of $r,c$. This argument gives no numerical upper bound for the possible exceptions. [Roth, *Rational approximations to algebraic numbers*](https://doi.org/10.1112/S0025579300000644)

## Every exponent at a small prime base

Suppose $z=b^e$ for an odd prime $b\le10^6$. For $P>10^{1000}$ we have $n>6600$, $2e\le n$ and $cb^{2e}2^{-n}\ne1$ because the 2-adic order of $c$ has absolute value below 42. The rational multiplicative form of Matveev's theorem applied to $b,2,c$ gives
$$
          \log|cb^{2e}2^{-n}-1|>-10^{15}(1+\log n).
$$
The degree is one, the exponents are $2e,-n,1$ and admissible logarithmic heights are $14,1,30$. The constant $1.4\cdot30^6\cdot3^{4.5}\cdot14\cdot30$ is below $10^{15}$. Comparing with (1) yields $n<10^{20}$, since $(n/2-340)\log2>10^{15}(1+\log n)$ at and beyond that point. [Matveev's theorem](https://www.mathnet.ru/eng/im314) · [Explicit multiplicative formulation, Theorem 2.1](https://doi.org/10.1016/j.jnt.2020.05.005)

Set $\alpha=\log_2 b$ and $\beta=\log_2 c$. Equation (1) gives $|2e\alpha+\beta-n|<2^{342-n/2}$. For all 78,497 odd prime bases through one million, the checker obtains coprime integers $t,u$ with $2^{83}\le t<2^{117}$ and $|t\alpha-u|<2^{-78}$. For every pair $(b,c)$ outside the four cases below it selects an approximant with $\|t\beta\|>1/16$, where the norm denotes distance to the nearest integer. The triangle inequality and $2e<10^{20}$ give
$$
              t|2e\alpha+\beta-n|>1/32,
$$
incompatible with $n>6600$. Continued fractions propose $t,u$ and 600-bit Arb intervals certify all inequalities.

In the four cases with $c=2^v b^w$, we have $|v|,|w|\le42$. If $2e+w\ne0$, then $0<|2e+w|<t$ and $(2e+w)u-(n-v)t$ is a nonzero integer by coprimality. Replacing $u$ by $t\alpha$ costs less than $1/32$ and again gives a contradiction. If $2e+w=0$, the absolute linear form is the nonzero integer $|v-n|$. The base 2 was already covered by the cofactor theorem.

## Counting the remaining weights

For a fixed residual shape take any odd prime $r\mid D$. The same row must add without carries at $r$. The identity $3N=Dt+DQ+3$ and the binary estimate give
$$
             t=(3\lambda/D)2^n+O(2^{n/2}).
$$
Write $x=t/r^{\lfloor\log_r t\rfloor}\in[1,r)$, and choose $w$ with $r^w>i+j$. The fractional base-$r$ digits of $x$ cannot contain a block of $w$ consecutive digits $r-1$. Such a block starting just after place $k$ puts $\{r^kx\}$ in $[1-r^{-w},1)$. This makes $\{ir^kx\}+\{jr^kx\}>1$, contradicting carry freedom of the first two parts at that layer. Use the terminating expansion at an endpoint.

Put $L=\lceil\log_2(4X^2)\rceil$, so it suffices to count $n\le L$. The rotation by $\theta=\log_r2$ has interval discrepancy $O(L^{-\eta})$ for an effective $\eta>0$. Indeed the two-logarithm bound gives $\|v\theta\|\gg v^{-M}$ for an effective fixed $M$, by applying it to $2^v r^{-a}-1$ with $a$ the nearest integer to $v\theta$. Dirichlet approximation gives a reduced fraction $a/v$ with $v\le\sqrt L$ and $|v\theta-a|\ll1/\sqrt L$. Hence $v\gg L^{1/(2M)}$. Each block of $v$ rotation points differs from a translated equally spaced grid by $O(1/v)$ in each position, and the final incomplete block costs $O(v/L)$. This proves the discrepancy estimate uniformly in a translation. [Matveev's theorem](https://www.mathnet.ru/eng/im314)

Avoidance in $m$ disjoint $w$-digit blocks cuts the length in $[1,r)$ to at most $(r-1)(1-r^{-w})^m$ and leaves at most $r^{wm+1}$ intervals. Taking $\log_r$ changes length by at most a fixed factor. For $n>\sqrt L$, $\log_r t$ differs from $n\theta+\log_r(3\lambda/D)$ by $O(2^{-\sqrt L/2})$. Enlarge the allowed intervals on the circle by this amount and use the discrepancy estimate. Choosing $m=\lfloor\eta\log L/(2w\log r)\rfloor$ bounds the number of admissible $n$ by $O(L^{1-\delta})$ for an effective $\delta>0$. The discarded $\sqrt L$ indices and the bounded number of inputs per index preserve this estimate. Taking the minimum over finitely many shapes and limiting values proves the assertion.

Run the certificate from the repository root with python-flint installed.

~~~sh
python3 check_sparse_quadratic_prefixes.py
python3 check_sparse_quadratic_prefixes_cross.py
~~~
