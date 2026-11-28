# A bound on the missing divisor for every primitive sum through 2,000

For a positive primitive triple $y$, write $u=\sum y_i$, $N=uV$ and $K=N!/\prod_i(Vy_i)!$. The full integer quotient
$$
                       R=\frac{N-1}{\gcd(N-1,K)}
$$
satisfies **$R<u$** in the following ranges.

- For every scale $V\ge1$ if $u\le126$.
- For every weight $N\le10^{1000}$ if $u\le2000$.
- For $u\le2000$, outside a common set of $O_\epsilon(X^\epsilon)$ weights $N\le X$, for every $\epsilon>0$.

This includes prime powers partly missing from the coefficient. The [earlier bound](avoided_predecessor_directions.md) counts full powers only at primes entirely absent. Since $u$ divides both $N$ and $K$, the endpoint gcd is at least $u(N-1)/R\ge N$ throughout these ranges. Any two rows enjoying this last bound share an endpoint prime. Hence a counterexample through weight $10^{1000}$ would have to include a nonendpoint row whose parts sum to more than 2,000 after division by their common gcd.

## Listing all possible violations

The [ray constant](sharp_direction_predecessor.md) satisfies $R\mid E_y\mid H$, with $H$ the product of the three pair gcds. Equality $R=u$ is impossible because $u$ divides $N$ and $R$ divides $N-1$. The exact constant computation gives $E_y<u$ for every $u\le126$. At $u=127$ distinct scales can witness the prime factors of $E_{(10,52,65)}=130$.

For a larger failure set $t=R>u$. Each prime power in $t$ divides exactly two primitive coordinates by the retained-power lemma and primitivity. Assign it to the remaining coordinate. Every coordinate must receive some power, since an empty group would make that positive coordinate a multiple of $t$ despite being less than $u$. With three group products $d,e,f>1$, the coordinates are positive multiples of $ef,df,de$, in that order. Conversely we can enumerate such primitive triples whose sum is below $t=def$ and test whether each prime power in $t$ divides $E_y$. Allowing all $t>u$ dividing $E_y$ includes every violation. A second enumeration visits the primitive coordinates directly, computes the maximal local losses using the original row carry graph, and lists every divisor of the resulting constant above $u$.

| Bound on $u$ | Group assignments before local tests | Locally possible $(y,t)$ |
|---:|---:|---:|
| 512 | 3,890 | 47 |
| 1,000 | 22,850 | 347 |
| 2,000 | 133,318 | 2,243 |

The two lists agree entry by entry through 2,000. The last list contains 2,209 distinct sorted directions.

## Testing all scales in the finite range

For each surviving pair use $V\equiv u^{-1}\pmod t$. At $p^j\parallel t$ let $p^e\parallel H$ and decrement the unit coordinate and a coordinate of valuation $e$ by one each, **after** multiplying by $V$. Denote these offsets by $\delta_i$. The loss is at least $j$ exactly when the valuation on the decremented triple is at most $e-j$, by the factorial identity for $E_y$. For $w=p^k$ its layer cost is
$$
       \left\lfloor\frac{uV-2}{w}\right\rfloor-
                  \sum_i\left\lfloor\frac{y_iV-\delta_i}{w}\right\rfloor .
$$
Every cost is nonnegative. Moving forward, any positive layer stays at least as large until a denominator floor increases. Moving backward, it stays at least as large while its numerator floor stays fixed. Intersect these intervals for enough layers of a single prime to exceed its budget. The entire resulting interval is excluded. The two traversals check respectively 23,877,209 and 25,615,292 intervals, with zero simultaneous scales through the stated weight bound.

## Beyond the finite range

For the two smallest prime bases of each pair impose their local tests. At $p^j\parallel t$, write $V=r+p^j n$ with $r=u^{-1}\bmod p^j$, $0<r<p^j$. The low $j$ digits have zero carry after decrement. The remaining digits of $n$ can be read with the multiplication carries of the three affine parts and their incoming addition carry. Track the remaining budget $e-j$ and charge the value of every outgoing addition carry. For block length 96 let $W_p$ be the number of labeled paths from all reachable states, allowing every endpoint. The code checks integer inequalities $W_p^{1000}<p^{96 k_p}$ for both bases. The largest sum of the two $k_p$ is 966.

Allow every block label in the corresponding base-$p^{96}$ digit set. Its dimension is at most $\log_p(W_p)/96$. A locally admissible $n$ has its blocks in this set up to any desired length by zero padding. Zero blocks are also available after flushing the carries. For $V\le X$ choose a power $B$ of $p^{96}$ with $X<B\le p^{96}X$. Then $n/B$ is in the compact digit set and $V/X$ lies in its affine image of slope $p^jB/X$ and translation $r/X$. The slopes are bounded independently of $X$. The two bases are multiplicatively independent and the sum of their dimensions is below one. Corollary 1.17 therefore gives $O_\epsilon(X^\epsilon)$ intervals of length $1/X$ for their intersection. Each contains at most two distinct scaled integers. Taking the finite union over all 2,243 assignments proves the assertion. [Corso and Shmerkin](https://arxiv.org/pdf/2409.04608)

Run the [verifier](check_retained_predecessor_directions.py) with the independent coordinate enumeration enabled.

~~~sh
python3 check_retained_predecessor_directions.py --limit 2000 --digits 1000 --independent
~~~
