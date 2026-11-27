# Prime powers forced near two equal part proportions

Fix an integer $m\ge1$. For a prime $p$ let $g_p(m)$ be the largest $j$ such that **all** positive trinomial coefficients in some neighborhood of the proportions $(1,1,m)/(m+2)$ are multiples of $p^j$, at every weight. The two smaller target coordinates may be perturbed independently, including when all three target coordinates coincide.

**The exponent $g_p(m)$ is either zero or exactly $v_p(m+1)$.** A finite graph below decides which value occurs, for every $m$ and $p$. The full value is always forced at $p=2$. At every odd prime it is forced whenever
$$
                   m+1=p^u-p^v,\qquad u>v\ge1.
$$
Consequently, for every odd prime $p$ and positive $a,b$, there is an open set of proportions throughout which every coefficient is divisible by $2^a p^b$. Take $v=b$ and $u=b+h$ with $2^a\mid p^h-1$. Any three rows of the same weight in the open set share this divisor.

**Complete classification through $m+1=1{,}000{,}000$.** At odd primes dividing $m+1$, the differences of powers above give all positive values except for
$$
                       m=18197,\qquad p=3,\qquad g_3(m)=3 .
$$
For instance the largest forced divisors at $m=2183$ and $18197$ are 312 and 54. Here $2184=3^7-3=13^3-13^1$ and $18198=2\cdot3^3\cdot337$. The finite classification imposes no bound on the weights or on the scales of the actual rows.

| $m$ | Largest forced divisor | Sufficient absolute error in each part proportion |
|---:|---:|---:|
| 5 | 6 | $1/200$ |
| 17 | 18 | $10^{-7}$ |
| 19 | 20 | $10^{-6}$ |
| 23 | 24 | $10^{-6}$ |

The [general rational-direction algorithm](rational_proportion_divisors.md) provides a computable tolerance at every $m$. The table uses shorter leading-prefix certificates.

## The two ways a limiting multiplier can behave

Put $D=m+2$. The rows $(p^n,p^n-1,mp^n)$ tend to the target direction. Their low $n$ columns have no carries. Above them the addition is $1+0+m$, giving valuation $v_p(m+1)$. Thus $g_p(m)\le v_p(m+1)$. On the exact direction the fixed divisor is $(m+2)(m+1)$, since the primitive sum $m+2$ and the primitive pair sum $m+1$ both divide the coefficient at every scale, and scale one attains their product.

Suppose $v=v_p(m+1)>0$ and there are rows approaching the target direction with valuation below $v$. Normalize by a power of $p$ as in the [proportion graph](rational_proportion_divisors.md). Passing to a subsequence, the limits of the normalized parts are $(t,t,mt)$ with $1/p\le Dt\le1$. Also pass to limiting digit words and carry amounts column by column, using the finite alphabet. In any finite prefix the outgoing carries sum to at most $v-1$, so past some column there are no further carries.

If $t$ has a terminating base-$p$ expansion, let $L$ be least with $p^L t$ integral. Since $t\le1/D$ and $D>p^v$, $L>v$. For $j=L-v,\ldots,L-1$ write $r=\{p^j t\}>0$. Then $(m+1)p^jt$ is integral, and the limiting fractional parts of the three scaled parts and the total are $r,r,1-r$ and $r$, respectively. All four are strictly between zero and one. The factorial floors at each of these $v$ places therefore stabilize on the approximating rows, with deficit one at every place. This contradicts the valuation bound. The witness rows may always be scaled by extra powers of $p$ to put these columns inside their expansions.

Otherwise both copies of $t$ have the same unique infinite expansion. After all carries have stopped, put $r_j=\{p^jt\}$ and let $b_j$ be the value of the unread fractional digits of the third part. The carry-free digit sums give $2r_j+b_j\le1$. Thus $0<r_j\le1/2$, $0\le b_j<1$ and $i_j=mr_j-b_j=\lfloor mr_j\rfloor$ is an integer in $[0,\lfloor m/2\rfloor]$. If $d$ and $e$ are the next digits of $t$ and the third part, then
$$
         m d+i_{j+1}=p i_j+e,\qquad 2d+e<p .
$$

## A graph and a sequence of witnesses

Use the integer states $0,\ldots,\lfloor m/2\rfloor$. From $i$, for a digit $d$, compute $m d+i=p i'+e$ with $0\le e<p$. Keep the edge to $i'$ when $2d+e<p$ and $i'$ is in range. The preceding infinite digit sequence traces this graph backwards and uses positive $d$ infinitely often. A finite graph then has a cycle with a positive digit. For $p\mid m+1$ the transitions have the convenient form
$$
 i\longmapsto \lfloor i/p\rfloor+\frac{m+1}{p}d,
 \qquad 0\le d\le\min(i\bmod p,p-1-(i\bmod p)).
$$

Conversely, take such a cycle of length $h$ starting at $j$. Read its input digits least significant first to obtain $x>0$, and similarly let $y$ be its output word. With $M=p^h-1$ one has $mx=jM+y$, and $x,x,y$ add without a carry, with $2x+y\le M$. Write $R_n=(p^{hn}-1)/M$. At every $n\ge1$ the parts
$$
          xR_n,\quad xR_n,\quad mxR_n+j=j p^{hn}+yR_n
$$
add without a carry. They tend to the target direction. This proves the dichotomy, gives explicit witnesses whenever $g_p(m)=0$ and $p\mid m+1$, and handles the possibility of different digit expansions at a terminating limit. At 2 the condition $2d+e<2$ prevents any positive cycle digit.

## The difference of powers

Suppose $m+1=p^u-p^v$ and a positive cycle exists. Its cyclic rotations have input integers $x_i$ and totals $w_i$, obtained by shifting a block cyclically by $i$ base-$p$ places. In particular $x_{i+a}\equiv p^a x_i\pmod M$ and $w_i\equiv(m+2)x_i\pmod M$. Each has $1\le x_i$ and $2x_i\le w_i\le M$ since addition is carry-free even after rotation.

Choose a largest input $X=x_i$ and put $Y=x_{i+u}$, $Z=x_{i+v}$. Then $w_i\equiv X+Y-Z\pmod M$, but $0<X+Y-Z<2X\le w_i\le M$. This contradiction excludes all cycles and proves the full valuation. [Lean checks the modular inequality](EqualPartCones.lean).

## The finite classification

For odd $p$ with $p\mid m+1$ but $p-1\nmid m+1$, take the single digit $x=1$ at $M=p-1$. The least residue $y=m\bmod M$ is at most $p-3$, and $1,1,y$ add without a carry. Thus only $p(p-1)\mid m+1$ needs examination. For each remaining pair outside the difference-of-powers theorem, search for $M=p^h-1$ and $x>0$ satisfying the block conditions above. The first search tries digit words with $x\le80000$. A complete graph search supplies the remaining positive witnesses, of lengths 54 and 30 at $(m+1,p)=(492048,3),(518238,3)$. Every reported block is checked by an independent factorial digit-sum calculation on two repetitions.

This gives 272,601 nontrivial positive-cycle pairs and 350 difference-of-powers pairs. At the last pair $(18198,3)$ the graph has 9,099 states. The two possible successors are $\lfloor i/3\rfloor$ and, only when $i\equiv1\pmod3$, $6066+\lfloor i/3\rfloor$. Removing zero and then states whose every successor has been removed proves that every path reaches zero within 72 steps. Independently, the Lean kernel checks a rank at every state. It is at most 72 and strictly decreases along every permitted edge from a nonzero state. Its rank table is proposed by a computation and its inequalities are proved by kernel reduction. No positive digit is allowed from zero.

To reproduce the census, the interval bounds, and the formal checks, run these commands in the project directory.

~~~sh
python3 check_equal_part_cones.py
python3 check_rational_cones.py
lean EqualPartCones.lean
~~~

The general proportion graph independently compares 1,055 valuation tests at $1\le m\le250$ against the equal-part classification. Exact binomial products also check the sharp upper-bound construction. The finite graph method here concerns common divisors in real neighborhoods of trinomial proportions. Related multiplication graphs for restricted digits in the $p$-adic integers are studied by [Abram and Lagarias](https://arxiv.org/abs/1308.3133).
