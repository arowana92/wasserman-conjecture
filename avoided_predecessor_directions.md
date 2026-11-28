# How much of the predecessor can one primitive direction avoid?

Let $y$ be a primitive positive triple, $u=\sum y_i$, and let $K$ be the coefficient of $Vy$ at weight $N=Vu$. Define
$$
             T(N,Vy)=\prod_{\substack{p^e\parallel N-1\\ p\nmid K}}p^e .
$$
Then $T(N,Vy)<u$ in each of the following ranges.

* At **every weight** if $u\le346$.
* At **every $N\le10^{1000}$** if $u\le2000$.
* For $u\le2000$, outside a single set of $O_\epsilon(X^\epsilon)$ weights $N\le X$, for every $\epsilon>0$. This exceptional set is uniform in the direction.

One consequence applies to a second row with two parts in any reduced ratio $a:b$ satisfying $C(a,b)\le a+b$. All integer ratios and all the [geometric ratios](part_ratio_predecessor.md) qualify, with arbitrarily large terms. In all three ranges **the two coefficients share a prime dividing $N(N-1)$**. Indeed, in a putative failure, at a full power of $N-1$ either the direction avoids the prime, or the other row does. The [companion divisibility proof](ratio_companion_divisor.md) gives the strict inequality $u(a+b)<\operatorname{lcm}(C(a,b),T(N,Vy))$. This would force $T(N,Vy)>u$.

The one-row assertion cannot hold for unrestricted primitive sums at all scales. The [explicit weight 222,788,354](endpoint_factors.md#an-obstruction-to-bounding-a-single-rows-avoided-product) has a row with scale divisible by 283 and avoided predecessor product $29\cdot127\cdot241>222788354/283$.

## Enumeration independent of the weight

Equality $T(N,Vy)=u$ is impossible since $u\mid N$ and $u\ge3$. If the product exceeds $u$, group its full powers by the unique coordinate congruent to one in the actual row at that power. All three groups are nonempty since each positive primitive coordinate is less than $u$. Write their products as $d,e,f$ and put $t=def$. Up to permutation the coordinates are positive multiples of $ef,df,de$, with sum less than $t$. They are primitive and satisfy $V\equiv u^{-1}\bmod t$. These conditions permit a finite enumeration of all pairs $(y,t)$ for any bound on $u$. Conversely, any such pair and multiplier avoiding all primes of $t$ gives $T(N,Vy)\ge t>u$, even if $t$ uses only part of a full power in $N-1$.

For $q=p^e\parallel t$ test whether **some** positive $V\equiv u^{-1}\bmod q$ makes $Vy$ carry-free in base $p$. After reading the fixed suffix, find the reachable multiplication-carry states with arbitrary higher digits. Accept exactly if a finite input completes the addition without carries. Different bases may use different local multipliers in this test.

| Bound on $u$ | Grouped rows with assigned product | Pairs $(y,t)$ passing every local test | Distinct directions |
|---:|---:|---:|---:|
| 512 | 3,890 | 4 | 4 |
| 1,000 | 22,850 | 87 | 87 |
| 2,000 | 133,318 | 713 | 708 |

The smallest locally admissible primitive sum is 347, with $y=(35,130,182)$ and $t=455$. The three bases are 5, 7 and 13. Existence of a separate multiplier for each base leaves the simultaneous question open.

The [independent enumerator](check_dominant_predecessor_cross.cpp) instead lists positive sorted triples by coordinates. For a primitive triple it forms $H=\gcd(y_1,y_2)\gcd(y_1,y_3)\gcd(y_2,y_3)$. For each prime power in $H$ it finds the largest locally admissible power. These powers and their divisors exhaust all choices of $t$, since the usable prime powers form an initial segment at each prime. The two lists of 713 pairs agree entry by entry.

## A common multiplier

For each pair, both finite scans enumerate the progression of $V$ up to $10^{1000}/u$. At power $q$ a violation is $\lfloor uV/q\rfloor>\sum_i\lfloor y_iV/q\rfloor$. Moving forward, this violation persists until some denominator floor increases. Moving backward, it persists as long as the numerator floor stays fixed. The two interval searches make respectively 6,982,872 and 7,452,208 checks and both find zero admissible multipliers.

For the asymptotic assertion it suffices in every one of the 713 cases to impose the first two prime bases in increasing order. At each base $p$ let $W_p$ count length-64 paths in the multiplication-carry graph starting at any reachable state. The code certifies integers $k_p$ with $W_p^{1000}<p^{64 k_p}$. The largest sum of the two numerators is 990. The independent bases thus have digit-block dimensions of sum strictly below one. The [uniform intersection argument](six_coefficient_row_relations.md#the-count-beyond-the-finite-range) gives $O_\epsilon(X^\epsilon)$ joint multipliers. There are finitely many pairs and each multiplier specifies a single weight. [Corso and Shmerkin, Corollary 1.17](https://arxiv.org/pdf/2409.04608)

~~~sh
c++ -O3 -std=c++17 check_dominant_predecessor_cross.cpp -o /tmp/check_dominant_predecessor_cross
python3 check_dominant_predecessor.py --cross /tmp/check_dominant_predecessor_cross
~~~
