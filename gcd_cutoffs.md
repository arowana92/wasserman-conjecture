# Pointwise lower bounds on the common divisor

Write $K_R=N!/(a!b!c!)$ for a positive trinomial row $R=(a,b,c)$ and $s_R$ for the sum of its two smallest parts. For **every** selection of three rows of weight $N$, write $G$ for the gcd of their coefficients.

* If every $s_R\ge5$, then $G\gg N^{1/3}$. If every $s_R\ge6$, then for each fixed $\epsilon>0$ one has $G>N^{1/2-\epsilon}$ for all sufficiently large $N$.
* Let $m$ be the smallest of all nine parts. Among the prime divisors of $N(N-1)\cdots(N-m+1)$, the number $r$ dividing all three coefficients satisfies $1+r\gg m\log(N/m)/\log N$. In particular, whenever $m$ tends to infinity, so does the number of distinct common primes. For $m\ge N^\alpha$ with fixed $0<\alpha<1$ there are $\gg_\alpha N^\alpha$ such primes.
* For **any two** positive rows $A,B$, their coefficients have a common divisor with $\binom N5$ of size $\gg N^{1/3}$, and with $\binom N6$ of size $\gg\sqrt N$ ($N\ge7$). In a hypothetical counterexample the two companions of the endpoint row must have at least $(1/10-o(1))N^{2/5}$ distinct common primes larger than $N^{2/5}$ dividing integers $N-j$ with $2\le j<N^{2/5}$.

All statements are pointwise in the weight. The first leaves the four possible small pairs $(1,1),(1,2),(1,3),(2,2)$ when the gcd stays bounded. The [endpoint theorem](all_weights_endpoint_reduction.md) narrows this to $(1,1)$ when the gcd equals one.

## A binomial product that retains the common valuations

For $1\le H\le N$ put $M_H=\operatorname{lcm}(1,\ldots,H)$, $D_H=\operatorname{lcm}(N,N-1,\ldots,N-H+1)$ and $F_t(x)=\prod_{j=1}^t\binom{x}{j}$, with $F_0=1$. For any nonempty collection of rows let $G_H$ be the gcd of $D_H$ and all their coefficients. In each row $R$ choose cutoffs $0\le t_{R,i}\le R_i$ summing to $H$. Then
$$
 F_H(N)\ \mid\ (M_H G_H)^H\prod_{R,i} F_{t_{R,i}}(R_i).               \tag{1}
$$

At a prime power $w=p^e$, each binomial has Legendre layer contribution zero or one. For $w\le H$ the contribution of $F_H(N)$ is at most $H$, supplied by $M_H^H$. For $w>H$ the contribution is $(H-h)_+$ with $h=N\bmod w$. If it is positive then $w\mid D_H$. If some row has zero carry across this layer, its residues $r_i$ sum to $h$ and its right-hand contribution $\sum_i(t_{R,i}-r_i)_+$ covers $(H-h)_+$. Every other relevant layer gives at least one unit of valuation to every row coefficient, and there are at most $v_p(D_H)$ such layers. Their number is at most $v_p(G_H)$, so $G_H^H$ covers them all. This proves (1) and allows $G_H$ to be replaced by $G$.

Here are estimates used below. Write $W(t)=t(t+1)/2$ and $L=\log(N/H)$. Then $F_H(N)\ge (N/H)^{W(H)}$ since $\binom Nj\ge(N/j)^j$. For $q=2,3$ dividing $H$, putting $t=H/q$ on $q$ positions of a row gives
$$
 \log\prod_i F_t(R_i)\le qW(t)L+O(H^2).
$$
Indeed $\binom{x}{j}\le(eN/j)^j$ and $\sum_{j\le t}j\log(t/j)=O(t^2)$. Also $\log M_H=O(H)$ by the elementary Chebyshev bound. If a selected part $m$ is within five units of $H$ and at least $H$, assigning it all $H$ costs at most $2^{mH}=\exp(O(H^2))$. The constants are absolute.

## Counting distinct common primes

For three rows with overall minimum $m\ge3$, let $H=3\lfloor m/3\rfloor$. Assign all $H$ to that minimum and split the other two rows equally into three cutoffs. The estimates in (1) yield
$$
 \log G_H\ge (H/6-1/2)\log(N/H)-O(H),                              \tag{2}
$$
because $W(H)-6W(H/3)=H^2/6-H/2$. For two rows the same bound holds at **every** positive multiple of three with $H/3$ at most every part, by omitting the selected row. Every exact prime power in $D_H$ is at most $N$, so $\omega(G_H)\ge\log G_H/\log N$.

Take $m$ large and $m\le cN$ for a sufficiently small absolute $c>0$. The error in (2) is absorbed by a fixed fraction of the main term, proving the stated bound for $r$. Bounded $m$ are covered by the additive one. For $m>cN$, every part is at most $N-2m$. Every prime in $((1-c/2)N,N]$ divides all coefficients and the product of the last $m$ integers. The prime number theorem gives $\gg_c N/\log N$ such primes for large $N$, which is enough since $m\le N/3$. For the consequences, $m\log(N/m)$ is increasing up to $N/3$. When $m\le\sqrt N$ the lower bound is a constant times $m$, and for larger $m$ it is at least a constant times $\sqrt N$.

At a possible counterexample the [part bound for the other two rows](endpoint_factors.md) gives $m_{A,B}^5>27(N-1)^2$. Apply (2) to that pair with $H=3\lfloor N^{2/5}/3\rfloor$. Its main term after division by $\log N$ is $(1/10+o(1))N^{2/5}$. There are $o(N^{2/5})$ primes at most $N^{2/5}$. A common prime of the two rows cannot divide $N$ or $N-1$, since the endpoint coefficient is $N(N-1)$. The remaining primes thus divide the stated shifts. Each of those integers can contain at most two primes greater than $N^{2/5}$, yielding at least $(1/20-o(1))N^{2/5}$ distinct offsets.

## A small gcd forces a row with second part below 60

Fix $0<d<1/2$ and suppose $G\le N^d$. For every sufficiently large $N$ some row has second-smallest part less than $60$.

First choose $1/2<\theta<\min(3/5,1-d)$. The [four nearby multiples lemma and its existence proof](near_row.md#producing-four-coprime-cofactors) give four primes exceeding $N^{9/10}$ such that, if none divides all three coefficients, one row has $s_R\le N^\theta$. In that construction, each prime is $\asymp N^{25/26-\eta/2}$ for a sufficiently small fixed $\eta>0$, and the four offsets are $O(N^{1/2+\eta})$. All four exceed $G$. The underlying short-interval count is due to [Jia and Liu, Sections 3 and 6](https://doi.org/10.4064/aa-95-1-17-48). For the cube-root conclusion alone one may instead use the largest prime below $N$: it exceeds $G$, its avoiding row has $s_R\le N-p$, and $N-p=O(N^{21/40})$ by [Baker, Harman and Pintz](https://doi.org/10.1112/plms/83.3.532).

For arbitrary parts $a,b,c$ chosen from different rows,
$$
                         N\mid Gabc .
$$
Indeed $N$ divides each row coefficient times any one of its parts, by decrementing that part. At a prime of $N$ use a row of least valuation. In particular, at most two rows have any part below $60$ for large $N$. If two have such parts, all parts of the third are at least $N/(60^2G)>N^\theta$ for large $N$.

Suppose now that all second parts are at least $60$. The row with $s_R\le N^\theta$ supplies a part in $[60,N^\theta]$. Select the least part $m\ge60$ across the three rows and set $H=6\lfloor m/6\rfloor$. When there are two rows with small minima, the selected part belongs to one of those two rows, since the third has all parts above $N^\theta$. In every case the other two rows have at least three and two positions respectively of size at least $H$ after possibly interchanging the rows. Use equal cutoffs in these positions and assign all $H$ to $m$ in its own row. Equation (1) gives
$$
 (1/12-1/(2H))L\le (\log G)/H+O(1),                               \tag{3}
$$
using $W(H)-3W(H/3)-2W(H/2)=H^2/12-H/2$. But $L\ge(2/5)\log N$, $H\ge60$, and $(1/12-1/120)(2/5)>1/120$. This contradicts (3) as $N$ grows. The bound 60 is independent of $d$, although the threshold for $N$ may depend on $d$.

## Height six and the four small pairs

Suppose every $s_R\ge5$. If $G\le N^d$ for a fixed $0<d<1/2$ and large $N$, the preceding result supplies a row with smaller parts $x,y\le59$. Allocate $\min(6,x+y)$ cutoffs to them and any remainder to its largest part. The first two cost at most $2^{354}$ in total. The remaining cost is at most $N$ when $s_R=5$ and is one when $s_R\ge6$. For each ordinary row whose minimum is at least two use $(2,2,2)$. For minimum one use $(1,2,3)$ in increasing part order. These cutoffs fit, and each ordinary row costs at most $N^9$ since $F_1(1)=1$. With $H=6$, (1) now yields
$$
 G^6\ge\frac{N^{2}}{6^{21}60^6 2^{354}},
 \qquad
 G^6\ge\frac{N^{3}}{6^{21}60^6 2^{354}}\quad\text{if every }s_R\ge6.
$$
For the first conclusion take $d=1/3$ and include the case $G>N^d$ directly. For the second, a fixed $d<1/2$ is incompatible at sufficiently large $N$ with both $G\le N^d$ and the square-root lower bound. This proves both pointwise assertions.

## Two rows and a fixed binomial coefficient

For $h=5,6$ and $N\ge7$ take $R=(1,h-1,N-h)$, whose coefficient is $h\binom Nh$. For $A,B$ different from the endpoint, the two ordinary cutoff choices just used fit, including when the second part is two. For $R$ give its first two parts cutoffs $1,h-1$ and its last part cutoff $6-h$. This costs at most $F_{h-1}(h-1)N^{6-h}$. Therefore
$$
 \gcd\left(K_A,K_B,\binom Nh\right)^6\ \ge \frac{N^{h-3}}{6^{21}(60h)^6 F_{h-1}(h-1)} .
$$
Here dividing the selected coefficient by $h$ reduces its gcd by at most $h$.

If either ordinary row is the endpoint, let $K$ be the other coefficient and set $D=\gcd(N,K), J=\gcd(N-1,K)$. The parts divided by $N/D$ are positive integers of sum $D$. At each prime power of $(N-1)/J$ with full exponent there, at least two normalized parts are divisible by it, by the [partial-power carry lemma](endpoint_pair_bound.md#retaining-partial-prime-powers). Consequently $((N-1)/J)^2\le (D/3)^3$ and $DJ\ge3(N-1)^{2/3}$. Since $N(N-1)\mid h!\binom Nh$, the common divisor with the fixed binomial coefficient is at least $DJ/h!$. This covers repeated endpoint rows as well.

The standard-library checker uses exact coefficients to test prime valuations for all allowed cutoffs through height seven, divisibility by the weight for triples, and the fixed-binomial bounds for pairs.

~~~~sh
python3 check_gcd_cutoffs.py
~~~~
