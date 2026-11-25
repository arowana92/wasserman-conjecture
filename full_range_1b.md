# The three-row conjecture through one billion

**Finite theorem.** For every integer $3\le N\le 1{,}000{,}000{,}000$, any three positive trinomial coefficients of weight $N$ have a common prime divisor. The triples may repeat. For $5{,}001\le N\le 10^9$ the common prime can always be chosen from the largest prime at most $N$ and the prime factors of the thirty consecutive integers $N-29,\ldots,N$.

## Selections containing the endpoint

The coefficient of $E=(1,1,N-2)$ is $N(N-1)$. If the other two rows avoided all its primes between them, we could assign their full powers as $UV=N$, $ST=N-1$. The first row has parts $Ux_i$ with $\sum x_i=V$. The second has parts $Vy_i$ with $\sum y_i=U$. Choose the names so that $S\ge V$. For every full power $w$ assigned to $S$, two normalized parts $x_i$ are 0 and the third is $V$ modulo $w$. All three positions must receive a nonzero residue from some power since $1\le x_i<V\le S$. Each assignment of positions determines the three $x_i$ uniquely. [Proof of the factor assignment](endpoint_reduction.md)

On each side the normalized sum must be at least 3, at least $2w+1$ for each assigned power of the predecessor, and its cube must be at least 27 times the square of the assigned predecessor product. After these tests the enumeration found **6,151,200,914** assignments of $(N,U,S)$ through one billion. For each, it enumerated every partition of the powers in $S$ into three nonempty groups, computed the normalized residues by the Chinese remainder theorem, tested their sum, and added the actual parts in every base assigned to $U$ or $S$. Zero rows passed. The counts in consecutive blocks of $10^8$ weights, starting at 3 in the first block, are

```
420279528 529308258 575443385 606473600 630129019
649360079 665540509 679592823 691987971 703085742
```

## Two-row divisibility at every prime-power layer

Write $B_t(x)=\prod_{j=0}^{t-1}(x-j)^{t-j}$ and $W(t)=t(t+1)/2$. On a row choose cutoffs $0\le t_i\le x_i$ of total $H$. At a layer $w=p^e$, write $u_i=x_i\bmod w$. The contribution of this layer to the $p$-valuation of $\prod_i B_{t_i}(x_i)$ equals $\sum_{k\ge0}\sum_i(t_i-u_i-kw)_+$. Here $(x)_+=\max(x,0)$. It is at least
$$G_w(H-\textstyle\sum_i u_i),\qquad G_w(d)=\sum_{k\ge0}(d-3kw)_+.$$
Put $h=N\bmod w$ and let $m$ be the largest integer at most $3w-3$ congruent to $h$ modulo $w$. Every row contributes at least $G_w(H-m)$. A row avoiding $p$ contributes at least $G_w(H-h)$ because the low digits add without a carry.

Fix a row $R$. If $p$ divides its coefficient, at least one of the other two rows must avoid $p$ in a counterexample. Thus a valid lower bound for their combined layer contribution is $G_w(H-h)+G_w(H-m)$. Even without this divisibility the two rows contribute at least $2G_w(H-m)$. These inequalities account for the small layers $w\le H$ as well as the large ones. [Earlier weighted lemma](weighted_obstruction.md)

## Exhausting the possible near-prime rows

The [direct prime tables](RESEARCH.md) settle weights through 5,000. At each higher weight the largest prime $r\le N$ exceeds $N/2$. A selection avoiding $r$ has a row $R=(a,b,N-s)$ with $1\le a\le b$, $s=a+b\le g=N-r$. In a counterexample omitting $E$, the two smallest parts of **every** row sum to at least 11 by Bergman’s Proposition 10, which uses only primes dividing $N(N-1)(N-2)$. [Bergman, Proposition 10](https://arxiv.org/pdf/0806.0607)

The checker tries to exclude all $11\le s\le g$ at once, then individual sums, then individual pairs. If $w>g$ and $N\bmod w<s$, the two small parts already exceed the permitted residue sum, so $p$ divides the coefficient of $R$. At a known $s$, it can also use all prime divisors of $\binom Ns$, and at a known $a$ the additional primes dividing $\binom sa$. It detects the divisibility by comparing the low residues at every power of the prime. For $w>H$, the only contributing layers have $w\mid N-h$ with $0\le h<H$, and contribute $H-h$ when the prime is forced. For $w\le H$ it uses the preceding two-row bounds. All tested heights are at most 30. In particular every prime used divides one of the thirty shifts or is $r$.

For the upper bound, order an ordinary row $x\le y\le z$. If all parts are at least $\lceil H/3\rceil$, give them cutoffs as equal as possible, putting the larger cutoffs on the smaller parts. With $d=\sum W(t_i)$ the cost is at most $(N/3)^d$ by the reverse ordering and the arithmetic-geometric mean inequality. More explicitly, opposite ordering gives $\sum_i W(t_i)\log x_i\le(d/3)\sum_i\log x_i$.

For a smaller first part $x$, assign it cutoff $x$. Put $u=\lceil(H-x)/2\rceil$, $v=H-x-u$. When $y\ge u$, the other two cost at most $(N/2)^{W(u)+W(v)}$ by the same argument. When $y<u$, both small parts are bounded, so assign them their actual sizes and give the rest of the cutoff to $z$. This fits since $z\ge N/3\ge H$. The cost is at most $B_x(x) B_y(y)N^{W(H-x-y)}$.

At least one of the two ordinary rows has minimum at least 3. Otherwise the full prime powers of $N$ assigned to those rows together divide 2, leaving the third row with three positive normalized parts of sum at most 2. If $\gcd(a,b,N)\le2$, both ordinary rows have minimum at least 3 by the same argument applied to $R$ and either ordinary row. Knowing $\gcd(N,s)\le2$ already suffices. The checker maximizes the upper bounds over the permitted minima and compares with the prime-layer lower bounds.

All comparisons use integer bounds on eight times the binary logarithm. For $e=\lfloor\log_2 x\rfloor$ and $m=\lfloor128x/2^e\rfloor$, a lower bound is $8e+\lfloor\log_2(m^8)\rfloor-56$, computed from bit lengths. Adding two gives an upper bound because $(129/128)^8<2$. At powers of two the lower bound is exact. The checker rejects a candidate only when the lower bound strictly exceeds the upper bound.

| Stage from 5,001 through one billion | Count |
|---|---:|
| Weights with $g<11$ | 462,506,183 |
| Additional weights excluded without fixing $s$ | 526,982,774 |
| Rows excluded by fixing $s$ at the remaining weights | 13,059,589,355 |
| Additional rows excluded by fixing $a,b$ | 1,438,045 |
| Unresolved rows | 0 |

## Reproduction and cross-checks

Use a C++17 compiler supporting a 128-bit integer type, such as Clang or GCC, and Python 3.9 or later. The executables accept inclusive bounds and can be run on disjoint intervals.

```sh
c++ -std=c++17 -O3 check_endpoint_range.cpp -o .work/endrange
c++ -std=c++17 -O3 check_nonendpoint_range.cpp -o .work/nonendrange
.work/endrange 3 1000000000
.work/nonendrange 5001 1000000000
python3 check_range_arithmetic.py
```

Create the local build directory with `mkdir -p .work` if needed. Both range commands write their counts to standard error and any survivors to standard output. The endpoint range command reports the first row of a surviving assignment, which would still require a complementary-row check.

The Python test uses exact binomial products and weighted factorials for 150,375 valuation checks, including 30,158 carry-free instances. It checks the integer logarithm bounds at the ends of every mantissa interval and reproduces the 6,084 and 5,066 endpoint assignment counts in the windows below by trial division and the earlier Python CRT search. Both implementations also produce the same actual row in a positive control at $N=222{,}788{,}354$, where the complementary size test fails.

```sh
.work/endrange 3 20000
.work/endrange 99999001 100000000
.work/endrange --assignment 222788354 283 887603
python3 check_wasserman_785.py
```
