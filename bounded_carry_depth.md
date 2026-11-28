# A fixed carry depth for almost every weight

**Theorem.** Outside a set of positive integers of natural density zero, every three positive trinomial rows of weight $N$ have a common prime $p$ at which each row has a carry within its first **128 base-$p$ digits**, counting from the units digit. The prime can depend on the three rows. The carries can occur at different positions.

More precisely, put $L=\log X$, $\ell=\log L$, $Z=X^{1/(512\ell)}$. Fix any cutoff $y=y(X)\ge1000$ with $\log y=o(\ell)$. Apart from $O(X\log(2y)/\ell)$ integers in $[X,2X)$ the following tests suffice simultaneously for all triples.

* Test the units digit at the largest prime at most $N$ and at primes $p>y$ dividing any of $N,N-1,\ldots,N-5$.
* At $y<p\le Z$ dividing $N(N-1)$, test the first 128 digits.

For **any two rows**, the same conclusion holds using just the tests at primes dividing $N(N-1)$. All witnesses in these assertions exceed $y$ for sufficiently large $X$. The digit bound is absolute, even at primes appearing to high multiplicity in a shift. The theorem concerns almost every weight. The examples below give a separate, finite result and exhibit the difference between absolute depth and depth measured past a full endpoint power.

Write $M=N-1$. A row passes a prime through depth $h$ when none of its first $h$ columns produces a carry. Equivalently its residues sum to $N\bmod p^j$ for each $1\le j\le h$. Thus if no common witness exists, each tested prime can be assigned to a row passing every column tested there.

## Counting the exceptional pairs

Keep only one factor $p$ at each tested prime $p>y$ of $NM$. Write $UV$ for the retained product from $N$ and $ST$ for that from $M$. The first row receives $U,S$ and the second $V,T$. At assigned weight primes all parts vanish modulo $p$. At assigned predecessor primes the residues are $1,0,0$ in some order.

Put $R_y(n)=n/\prod_{p\mid n,\ p>y}p$. Except at $O(X\log(2y)/\ell)$ weights,
$$
 \prod_{j=0}^5 R_y(N-j)\le L .                                      \tag{1}
$$
Indeed the mean of its logarithm is $O(\log(2y))$. For each shift, the first powers up to $y$ contribute $O(\sum_{p\le y}(\log p)/p)=O(\log(2y))$, and the sum of contributions from all powers $p^j$, $j\ge2$, is $O(1)$. These estimates follow by counting multiples among the integers up to $2X$, then dividing by $X$. Markov's inequality proves (1).

Outside another $O(X/L)$ weights **both** $US$ and $VT$ exceed $X/L^{12}$ in every passing assignment. Here are the details of the uniform divisor bound used for this assertion. Restrict to $\tau(N)\tau(M)\le L^4$, losing $O(X/L)$ by Cauchy–Schwarz and the second divisor moment. For a row receiving $V,T$, the [two-congruence reconstruction](two_predecessor_congruences.md#recovering-the-row-from-a-product) bounds the number of retained weights with $(N/V)(M/T)\le Y\le X$ by $O(Y L^9)$. In that proof the only divisor multiplicity in the first count comes from splitting $W=Vd_3\mid NM$, hence is at most $L^4$. Its second count has no divisor multiplicity. The smaller of the two bounds is therefore $O(Y L^4)$ in each of $O(L^5)$ dyadic boxes. The argument allows partial divisors of the weight and predecessor. Taking $Y=X/L^{11}$ suffices, since $(N/V)(M/T)=US R_y(N)R_y(M)\le US L$. Exchange the rows for the second assertion.

### Many primes with moderate digits in the weight

Put $k=127$ and $r=\lfloor 0.90\ell\rfloor$. For each $h=0,1$ call a prime $y<p\le Z$ eligible if $N\equiv h\bmod p$ and its digits at positions $1,\ldots,k$ have sum at most $3k(p-1)/5$, with the units position numbered zero. Apart from $O(X/\ell)$ weights, there are at least $r$ eligible primes for **each** $h$.

In a full residue system modulo $p^{128}$ the density is at least $(1-e^{-2k/100})/p>0.92/p$. Conditional on the units digit, the other digits are independent and uniform. The upper-tail estimate follows from the bounded-digit exponential moment, or Hoeffding's inequality, at deviation $1/10$ from the mean after division by $p-1$. Mertens's reciprocal-prime sum gives $\sum_{y<p\le Z}1/p=(1+o(1))\ell$ in this range. The variance of each count is $O(\ell)$. For distinct primes the residue-system indicators are independent by CRT. On the interval of length $X$, the total error in first and second moments from incomplete periods is $O(Z^{258}/X)=o(1)$, bounding each joint period by $Z^{256}$ and counting at most $Z^2$ pairs. Chebyshev's inequality gives the claim.

In a passing pair there is a row with **every part at least two** receiving at least $r$ eligible primes. If both rows have this property, assign the at least $2r$ primes and take the larger group. If a row contains a 1, it cannot take any weight prime, so the other row gets all the eligible primes with $h=0$ and all its parts are multiples of them.

Choose $r$ such primes assigned to this row. At each of them its digit sums agree column by column with the weight. Sum over the $k$ positions and the $r$ primes, dividing every digit by its own $p-1$ before summing. At least one part $z\ge2$ has total at most $rk/5$. Its units digit at each chosen prime is 0 or 1.

### Few possible values for that part

For a set of $r$ fixed primes, the product of the moduli $p^{128}$ is at most $X^{1/4}$. Thus the count of $z\le2X$ obeying the two digit conditions is at most $3X$ times their joint density in the product residue system. The independent normalized digits above the units place $D$ have
$$
 {\mathbb E}e^{-5D}=\frac1p\sum_{j=0}^{p-1}e^{-5j/(p-1)}
 \le \frac1p+\frac15 .
$$
This follows by bounding the decreasing sum after its first term by an integral. Exponential Markov at total $rk/5$ gives density at most $(11/20)^{rk}\prod_p(2/p)$, since $e(1/1000+1/5)<11/20$. Sum over unordered sets of $r$ primes. The sum of the products of $2/p$ is at most $(2e\sum_{p\le Z}1/p\,/r)^r\le7^r$ for large $X$. Hence the set $\mathcal Z$ of all possible parts satisfies
$$
       |\mathcal Z|\ll X L^{-65}.
$$
Indeed $r\ge0.89\ell$ eventually and $0.89(127\log(20/11)-\log 7)>65$.

Name the chosen row first. Then $U\mid z$ and $S\mid z(z-1)$. The primes in $U,S$ are disjoint, so there are at most $F(z)=3^{\omega(z)}2^{\omega(z-1)}$ choices. Each coprime pair gives one class for $N$ modulo $US$ and hence $O(L^{12})$ weights. For any fixed $b\ge1$ the mean of $b^{\omega(n)}$ up to $2X$ is $O_b(L^{b-1})$, by expanding over squarefree divisors with weights $(b-1)^{\omega(d)}$ and taking the Euler product. Hölder with exponents $3/2,3$ now gives
$$
 \sum_{2\le z\le2X} F(z)^2
 \le \left(\sum 27^{\omega(z)}\right)^{2/3}
      \left(\sum 64^{\omega(z-1)}\right)^{1/3}
 \ll X L^{115/3}.
$$
Cauchy–Schwarz on $\mathcal Z$ bounds the count of remaining weights by $O(X L^{12+(115/3-65)/2})=O(X L^{-4/3})$. This proves the pair assertion.

## Three arbitrary rows

Write $g=N-q$ for the gap below $N$ to its largest prime $q$. Outside $O(X/\ell)$ weights, $g\le X^{0.18}$ by Guth and Maynard's Corollary 1.4, since $2/15<0.18$. Take the integer interval length $J=\lfloor X^{0.18}\rfloor$ and the starting integer $N-J$. The corollary supplies a prime in $(N-J,N]$ for all but $O(J+X\exp(-cL^{1/4}))$ weights, with the $O(J)$ term covering the lower boundary. [Guth and Maynard, Corollary 1.4](https://ora.ox.ac.uk/objects/uuid%3Aad11b8bf-ad2b-4ebf-a627-647f023c378f/files/s73666708v)

If any row has its two smaller parts summing to at most 10, the other two have an endpoint witness $p>y$. At that prime the small parts already exceed the permitted units residue $N\bmod p\le1$. This witnesses a carry in the third row too.

Otherwise suppose none of the indicated tests gives a common witness. A row passing $q$ has its two smaller parts $a,b$ summing to at most $g$, since $q>N/2$. Every row has second part at least 6. Apply the [height-six weighted product argument](weighted_obstruction.md#a-weighted-divisibility-lemma) using only the retained primes at the units layer. In its notation the product of the three row costs is divisible by
$$
       \prod_{j=0}^5 \left((N-j)/R_y(N-j)\right)^{6-j}.
$$
Indeed any retained prime exceeds 6 and occurs in just one shift. At such a prime the assigned row's residues sum to $j$, so its factors supply at least $6-j$ copies. By (1) the displayed integer is $\gg N^{21}/L^6$.

On each of the two ordinary rows use cutoffs $(2,2,2)$ if its minimum is at least two and $(1,2,3)$ otherwise, in increasing part order. Each costs at most $N^9$. On the row passing $q$, use $(1,5,0)$ when $a=1$ and $(2,4,0)$ otherwise. Its cost is at most $g^{15}$. The product bound would give $N^{21}\ll L^6N^{18+15(0.18)}$, which fails for sufficiently large $X$. This proves the full assertion.

## Sharp finite depths through one hundred million

For the endpoint and two other rows there is an exact calculation measuring depth *relative* to the full power. If $p^e\parallel N(N-1)$ allow layers through $p^{e+j}$. With $j=0$ the single full-power residue test implies the tests at all lower layers, since the target residue there is always 0 or 1.

The table gives the **first weight where the indicated extra depth is insufficient**, along with two actual rows that pass all the tests between them. For every $3\le N\le100{,}000{,}000$, four extra digits suffice for every pair. The third and fourth rows of the table each have at most three common endpoint primes even in the full coefficients.

| Extra digits $j$ | First weight | First row | Second row | Exact endpoint gcd |
|---:|---:|---|---|---:|
| 0 | 771 | $(210,231,330)$ | $(257,257,257)$ | 2,310 |
| 1 | 97,876 | $(14500,27000,56376)$ | $(24469,24469,48938)$ | 391,500 |
| 2 | 2,556,249 | $(576264,958593,1021392)$ | $(852083,852083,852083)$ | 24 |
| 3 | 25,008,076 | $(5594400,7409800,12003876)$ | $(6252019,6252019,12504038)$ | 2,700 |

The C++ checker assigns the **full powers** as $UV=N$ and $ST=M$ to rows of the form $U x,V y$ with normalized sums $V,U$. Name the rows so $S\ge V$, possible since $ST=UV-1>(U-1)(V-1)$. Each full power $q\parallel S$ specifies one position with normalized residue $V\bmod q$ and the others zero. All three positions must be used on this side since its parts lie strictly between zero and $V\le S$. The other row uses the analogous residues with $T,U$, including unused positions. CRT and all positive lifts enumerate each possible row. At each assigned prime only the next $j$ residue sums need testing. Within an assignment any passing rows on the two sides can be combined. The size checks $2q+1\le V$ and $27 S^2\le V^3$ follow from the two positive multiples of each $q$ and the divisibility of the product of normalized entries by $S^2$.

The search increases its tested depth each time it finds a passing pair, retesting that same weight before moving on. This yields the first failures at all lower depths because passing is monotone in the depth. A Python cross-check enumerates **individual cuts** with the two possible residues at each predecessor power and then assembles positive triples. It independently recovers the first two weights through 100,000. It checks all four witnesses, their first common carry depths and the displayed gcds by exact factorial floors. At 771 it also constructs the coefficients by exact binomial products.

~~~sh
mkdir -p .work
c++ -std=c++17 -O3 check_endpoint_depth.cpp -o .work/check_endpoint_depth
.work/check_endpoint_depth 3 100000000 0
python3 check_endpoint_depth_cross.py
~~~

The [fixed-prime construction](fixed_prime_depth_obstruction.md) gives a positive-density obstruction to restricting the additional digits to a fixed finite list, even when every listed prime divides the endpoint exactly once.
