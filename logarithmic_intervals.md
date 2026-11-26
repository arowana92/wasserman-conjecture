# Whole orders of magnitude for prime-power cofactors

Write $W=q^e$ for a prime $q$ and $e\ge1$. Call the binary order $j$ exceptional if some $2^j\le W<2^{j+1}$ and some $1\le d\le1000$ give a counterexample at $N=1+dW$.

**In every block of $H$ consecutive nonnegative binary orders, at most $O(H^{1-\delta})$ are exceptional, with $\delta=1/2500000$.** The implied constant is independent of the first order. The same bound holds when orders are taken on the weights $N$ themselves. In particular, every sufficiently long block contains a run of $\gg H^\delta$ consecutive orders in which the conjecture holds throughout this family.

The following finite upper bounds have no threshold and hold for every starting order. The second column takes orders of $W$ or of $N$. The third takes base-seven orders of the cube $W=q^{3e}$. The fourth takes binary orders of $t=p^a q^b$ in the family of **even** weights $7t+1$, with both prime powers unrestricted.

| Length $H$ of block | Binary orders, all $d\le1000$ | Base-seven orders, cubes | Binary orders, even $7t+1$ |
|---:|---:|---:|---:|
| $10^{30}$ | $<1.37\cdot10^{26}$ | $<4.92\cdot10^{23}$ | $<8.36\cdot10^{23}$ |
| $10^{100}$ | $<3.08\cdot10^{80}$ | $<1.66\cdot10^{74}$ | $<3.32\cdot10^{64}$ |

Thus every block of $10^{100}$ binary orders contains at least $10^{19}$ consecutive orders free of exceptions for all the bounded cofactors at once. The [explicit recurring interval](recurring_intervals.md) certifies a particular run and also treats both parities at the cofactors five and seven with two unrestricted prime-power factors.

## Two leading prefixes and a rotation

For bases $p<r$, suppose every admissible positive integer has at most $C_p p^{sn}$ possible leading words of length $n$ in base $p$ and at most $C_r r^{tm}$ of length $m$ in base $r$, with $0<s+t<1$. Write $W=p^i u=r^j v$ with $1\le u<p$ and $1\le v<r$. Set $P=p^{n-1}$ and $R=r^{m-1}$. If $W\ge\max(P,R)$, its two leading words, regarded as integers $a,b$, give
$$
 a/P\le u\le(a+1)/P,\qquad b/R\le v\le(b+1)/R.
$$
Consequently $i\log_r p$ modulo one belongs to the circular image of the interval of possible $\log_r(v/u)$. Its length is at most
$$
 \frac{\log(1+1/a)+\log(1+1/b)}{\log r}
 \le \frac{P^{-1}+R^{-1}}{\log r}.
$$
This covers all compatible phases by one arc per prefix pair.

If $A/B$ is reduced and $|\log_r p-A/B|<B^{-2}$, the $B$ successive phases in any block starting at $i_0$ are within $1/B$ of the translated grid $i_0\log_r p+kA/B$, $0\le k<B$. An arc of length $l$ therefore receives at most $Bl+3$ visits. The same circular grid count applies across the origin. Take $n=1+\lfloor\log_p B\rfloor$ and $m=1+\lfloor\log_r B\rfloor$. Each arc receives $O_{p,r}(1)$ visits, and there are $O(B^{s+t})$ arcs. Orders with $W<\max(P,R)$ contribute at most $O(\log(2B))$. This proves a uniform bound $O(B^{s+t})$ for **every block of $B$ orders**, in particular for all convergent denominators of the irrational number $\log_r p$.

For general $H$ take its largest convergent denominator $B\le H$, with next denominator $B'>H$. A lower bound $|B\log_r p-A|\gg B^{-K}$ gives $H<B'\ll B^K$ by the convergent inequality. At most $\lceil H/B\rceil$ blocks suffice, giving $O(H^{1-(1-s-t)/K})$ orders. Even the elementary inequality $|p^B-r^A|\ge1$ gives $B\gg\log(2H)$ and hence the uniform bound $O(H/(\log(2H))^{1-s-t})$.

For the numerical exponent use Gouillon's Corollary 2.3 on $B\log p-A\log r$, with $D=1$ and $\log A_1=\max(1,\log p)$, $\log A_2=\max(1,\log r)$. Distinct primes are multiplicatively independent. The coefficient multiplying $h$ in the negative log bound is $7200(3.409+1.705)\log A_1\log A_2$. In his notation $b=B/\log A_2+A/\log A_1=O(B)$ and $h=\max(\log b+3.1,1000,768)\le\log B+O(1)$. All selected bases are at most 19, so both heights are less than 3 and $K=332000$ suffices. Finitely many initial convergents, including a zero numerator, are absorbed into the constant. [Nicolas Gouillon, *Explicit lower bounds for linear forms in two logarithms*, Corollary 2.3](https://www.numdam.org/article/JTNB_2006__18_1_125_0.pdf)

## Applying the bounds to trinomial rows

In any failing $1+dq^e$, $d\le1000$, the [primitive denominator reduction](shifted_directions.md) forces a row $(aW,bW,cW+1)$ with $W=q^e$ and $a+b+c=d$. The [two subsequent exclusions](denominator_residue_hulls.md) leave 14 denominators. For every remaining row the checker selects two required bases $p,r$ at which the row must add without carries. It obtains 16 distinct choices of denominator, row and pair of bases.

Let $T_p(L)$ count paths of length $L$ in the [carry graph](carry_intersections.md), starting and ending at arbitrary reachable states. Each leading word is readable along some path. Splitting such paths into blocks of length $L$, with a bounded remainder, gives at most $O(p^{sn})$ leading words whenever $T_p(L)\le p^{sL}$. Length 64 gives $s+t\le 0.86$ throughout. Therefore $(1-s-t)/332000>1/2500000$. Each base-$p$ order meets at most $1+\lceil\log_2 p\rceil$ binary orders. The same is true for $N=1+dW$ at a fixed $d$, since the ratio of its endpoints across one base-$p$ order is less than $p$. A window of $H$ binary orders of $W$ or of $N$ meets at most $H+2$ base-$p$ orders of $W$ per $d$: after subtracting 1 the ratio of positive endpoints for $N\ge2$ is at most $2^{H+1}$. Sum over the 16 cases. At most $M$ occupied orders leave a run of at least $(H-M)/(M+1)$ empty orders.

At cubes only $d=986$ survives. The two possible slopes are $(119,203,664)$ and $(203,714,69)$, both required to avoid 7 and 17. Length-1024 path counts give exponents $(0.287,0.439)$ and $(0.251,0.412)$, respectively. Thus every window of length a convergent denominator $B$ of $\log_{17}7$ contains at most $O(B^{0.726})$ possible base-seven orders. [Cube restrictions](common_power_exponents.md)

For even $7t+1$ with $t=p^a q^b$, a failure forces $6t$ and $t+1$ to add without carries at 2. Also the positive base-seven digits of $t$ are weakly increasing from the most significant to the least significant, and no internal digit is zero. Hence there are at most $\binom{m+5}{5}$ leading words of length $m$ at 7. At 2 the length-1024 path count has 569 bits. The preceding proof gives $O(B^{569/1024}(1+\log B)^5)$ occupied binary orders at convergent denominators of $\log_7 2$, and the power saving for general $H$ follows as above. [Full-factor-seven reduction and digit inequality](predecessor_five_seven.md)

There is also a consequence for **orders of the prime bases $q$**, simultaneously over a growing range of exponents. In any $H$ consecutive binary orders of $q$, all but $O(H^{1-\delta/3})$ orders work at every $d\le1000$ and every $1\le e\le H^{\delta/3}$. For each $e$, the possible $q^e$ occupy a window of $eH$ orders. Every occupied order of $q^e$ determines a unique order of $q$, by integer division by $e$. Summing $O((eH)^{1-\delta})$ over the indicated exponents proves the claim.

## Finite computation

For the table, count distinct prefixes by reading backward from the accepting set, grouping words that give the same set of predecessor states. With $L_n,M_m$ allowable prefixes, a block of $B$ indices has at most
$$
 L_n M_m\big(3+B(P^{-1}+R^{-1})/\log r\big)
$$
large-index visits. The checker counts the small orders separately and tests nearby choices of $n,m$. It encloses each logarithm with outward integer rounding in the series $\log(a/b)=2\sum_{j\ge0}x^{2j+1}/(2j+1)$, $x=(a-b)/(a+b)\le1/3$ after binary range reduction. With $K$ terms the tail is below $9^{-K}$, so the chosen $K$ makes it smaller than one unit of fixed precision. Continued fractions of an interior rational supply candidates. The inequality $|B\log_r p-A|<1/B$ is verified at **both** rational endpoints.

Run from the repository root with Python and python-flint.

~~~~sh
python3 check_logarithmic_intervals.py
~~~~
