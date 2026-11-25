# Optimal four-prime certificates at consecutive weights

**Theorem.** For all sufficiently large $X$, at least $cX/(\log X)^2$ distinct integers $n\in[X,2X-1]$ have the following property, with an absolute $c>0$. At each of $n$ and $n+1$ there is a list of four primes such that any three positive trinomial coefficients share a prime on the list. The minimum list size at each weight is exactly four. Three of the four primes can be the same at both weights.

For example, at 187 the primes are $11,17,31,181$, and at 188 they are $11,17,47,181$. The existence proof below applies directly to all triples of rows. The [earlier construction](bounded_certificates.md) gave lists of size at most seven at individual odd composite weights.

## A criterion for one weight

Call two distinct primes incompatible at weight $N$ if no positive row of weight $N$ avoids both, where avoidance means that the prime does not divide the coefficient. Four pairwise incompatible primes certify $N$ because three rows avoid at most three of them.

For full prime powers of distinct bases use the following rules from Kummer's criterion. At $P\mid N$, an avoiding row has all three parts divisible by $P$. At $S\mid N-1$, the residues of an avoiding row are $0,0,1$ in some order. Thus:

- Two powers $P,Q$ of $N$ are incompatible if $3PQ>N$.
- Two powers $S,T$ of $N-1$ are incompatible if $ST>N-2$. At least one positive part would be divisible by both, and no part exceeds $N-2$.
- A power $P$ of $N$ and a power $S$ of $N-1$ are incompatible if $P(2S+1)>N$. Two parts would be multiples of $PS$, and the third a positive multiple of $P$.

These rules also give a sufficient condition for the full conjecture at any single weight even without a fourth prime. Form the graph on prime divisors of $N(N-1)$, joining incompatible primes. If it contains an odd cycle, the [all-weight endpoint reduction](all_weights_endpoint_reduction.md) excludes a counterexample. The endpoint coefficient is $N(N-1)$, and assigning every cycle prime to one of the two other rows avoiding it would give a proper two-coloring.

To obtain a four-prime list directly, let $r>N/2$ be a prime at most $N$. A row avoiding $r$ has two parts of sum at most $N-r$. If $N-r<2P$ it cannot avoid the base of $P\mid N$, and if $N-r\le S$ it cannot avoid the base of $S\mid N-1$. At least one of the two small parts in the latter case would be a positive multiple of $S$.

In particular take $n=pq$ for distinct odd primes. Suppose odd primes $s_-\mid n-1$ and $s_+\mid n+1$ exceed $\max(p,q)$. If a prime $r\le n$ satisfies $2r>n+1$ and $n+1-r<\min(p,q)$, then
$$\{p,q,s_-,r\}\text{ certifies }n,\qquad \{p,q,s_+,r\}\text{ certifies }n+1.$$
At the first weight $p,q$ divide the weight; at the second they are the two factors of its predecessor. Each of the two remaining mixed pairs is incompatible by the third rule. The near prime is distinct from the others since the neighbors of $n$ are even and $s_\pm\le(n+1)/2$.

## Large prime factors on both sides of a semiprime

Put $Y=\sqrt X$, $I=(Y,4Y/3]$, and let $p,q$ range independently over primes in $I$ congruent to 1 and 3 modulo 4, respectively. Each ordered pair gives a different product. Write $A,B$ for the two prime counts and $W=AB\asymp Y^2/(\log Y)^2$ by the prime number theorem in arithmetic progressions. Use uniform probability on the $W$ pairs. For each sign put $M=pq\pm1$ and let $B_\pm$ be the event that $M$ has a prime factor exceeding $2Y$.

We prove $\Pr(B_-\cap B_+)\ge1/240-o(1)$. All estimates below hold separately for either sign. Throughout, $M\le2Y^2$ for large $Y$, so at most one prime divisor can exceed $2Y$. The main estimate is
$$\mathbb E\sum_{\ell\le2Y}v_\ell(M)\log\ell\le(1+o(1))\log Y, \tag{1}$$
with $\ell$ prime. The expected contribution from a prime exceeding $2Y$ is then at least $(1-o(1))\log Y$, since $\log M=2\log Y+O(1)$.

**Moduli through $D=Y^{2/5}$.** For any fixed $K>0$,
$$\sum_{d\le D,\ d\text{ odd}}\left|\Pr(d\mid M)-\frac1{\varphi(d)}\right|\ll_K(\log Y)^{-K}. \tag{2}$$
For each $p$, the allowed $q$ lie in one reduced class modulo $4d$, combining $q\equiv3\pmod4$ and $pq\equiv\mp1\pmod d$. Apply Bombieri–Vinogradov to counts of primes in $I$ with moduli at most $4D$, using the maximum error over all reduced classes. Divide by $B$ and average over $p$. The main term is $(\operatorname{Li}(4Y/3)-\operatorname{Li}(Y))/(2B\varphi(d))$. Replacing the factor before $1/\varphi(d)$ by 1 costs an arbitrarily large logarithmic saving by the fixed-modulus prime number theorem and $\sum_{d\le D}1/\varphi(d)\ll(\log Y)^2$. [Bombieri, *On the large sieve*](https://doi.org/10.1112/S0025579300005313)

For odd prime powers $\ell^j\le D$, (2) bounds their total expected log contribution by $\log D+O(1)$. The $j=1$ terms use $\sum_{\ell\le z}\log\ell/(\ell-1)=\log z+O(1)$, and the remaining main terms sum to $O(1)$. Powers of 2 contribute $O(\log\log Y)$ in expectation: factor out at most $O(\log\log Y)$ digits unless $M$ is divisible by the first power of 2 greater than $(\log Y)^5$. The latter occurs for $O(Y^2/(\log Y)^5+1)$ pairs by uniqueness of $pq$, each contributing $O(\log Y)$.

All higher powers $\ell^j>D$, $j\ge2$, together cost $o(1)$ in expectation. For each such power $v$ there are $O(Y^2/v+1)$ pairs. On $v\in[Z,2Z]$ there are $O(\sqrt Z\log Z)$ such powers counted with their exponents. Summing the resulting dyadic bounds through $2Y^2$, including their logarithmic weights, gives $O((Y^2/\sqrt D+Y)(\log Y)^3)/W=o(1)$.

**Primes between $D$ and $D_0=Y/(\log Y)^8$.** In a dyadic interval $L<\ell\le2L$ restricted to this range, character orthogonality gives the principal term $1/(\ell-1)$ in $\Pr(\ell\mid M)$. The total absolute error over the interval is
$$\ll \frac{Y+L^2}{L\sqrt W}.$$
Indeed the nonprincipal characters modulo a prime are primitive. Apply Cauchy–Schwarz to the products of the two character sums and the multiplicative large sieve to each prime indicator. The squared norms of these indicators are $A,B$, and the character normalization is $O(1/L)$. Summing over intervals and inserting the weight $\log\ell$ still gives $o(1)$ error because $D/2\le L\le D_0$. The main term in this range is $\log D_0-\log D+O(1)$. [Montgomery and Vaughan, *The large sieve*](https://doi.org/10.1112/S0025579300004708)

**Primes between $D_0$ and $2Y$.** We give the upper bound for this interval, where the preceding large-sieve error would be too large. For each prime $\ell$ in it the number of ordered prime pairs in $I^2$ with $pq\equiv\mp1\pmod\ell$ is
$$\ll \frac{Y^2}{\ell(\log Y)^{3/2}}+Y^{1/2+o(1)}. \tag{3}$$
The implied bound is uniform. To see this, start with all integer pairs $a,b\in I$. Impose $m\mid a$, $t\mid b$ for integers $m,t$ coprime to $\ell$. The rectangle count for the modular hyperbola, applied to the intervals $I/m,I/t$ and the invertible right-hand side $\mp(mt)^{-1}$, gives
$$\frac{\ell-1}{\ell^2}\frac{|I|^2}{mt}+O(Y^{1/2+o(1)}).$$
Here $|I|=Y/3$. Split each interval into $O(1+Y/\ell)=O((\log Y)^8)$ intervals of length at most $\ell$. On each rectangle the bound is $O(\ell^{1/2+o(1)})$, uniformly in its endpoints and the invertible residue. This follows by Fourier completion from Weil's Kloosterman sum bound. Rounding interval lengths contributes $O(1+Y/\ell)$ overall when forming the main term with real lengths. [Weil, *On some exponential sums*](https://doi.org/10.1073/pnas.34.5.204) · [Shparlinski, Theorem 13](https://arxiv.org/pdf/1103.2879)

Put $z=\exp((\log Y)^{3/4})$ and let $J$ be the least even integer at least $100\log\log Y$. Prime pairs must avoid each of the events $s\mid a$, $s\mid b$ for all primes $s\le z$. Truncate inclusion–exclusion after degree $J$ for an upper bound. Every intersection is one of the preceding rectangle counts with $m,t$ products of the indicated primes. The number of terms is at most $(3z)^{J+1}=Y^{o(1)}$, and $z<D_0$. The truncated main term, relative to $(\ell-1)|I|^2/\ell^2$, is at most
$$\prod_{s\le z}(1-1/s)^2+\frac{(2\sum_{s\le z}1/s)^{J+1}}{(J+1)!}\ll(\log Y)^{-3/2}.$$
This is the Bonferroni bound applied to independent events of probabilities $1/s$; the first term uses Mertens' product estimate and the second is $O((\log Y)^{-10})$ by the choice of $J$. This proves (3). Since $\sum_{D_0<\ell\le2Y}(\log\ell)/\ell=O(\log\log Y)$, summing (3) with the log weight costs $o(W\log Y)$. This proves (1).

**Obtaining both signs.** Let $Z$ count prime factors of $M$ in $(Y^{1/40},Y^{1/5}]$. Equation (2) applies to each such prime and to every product of two distinct ones. Each product has just two orderings. Consequently
$$\mathbb EZ=\lambda+o(1),\qquad \mathbb EZ^2=\lambda^2+\lambda+o(1),\qquad \lambda=\log 8>2,$$
by the reciprocal-prime sum. Cauchy–Schwarz gives $\Pr(Z>0)\ge\lambda/(1+\lambda)-o(1)>2/3$ for large $Y$.

The unique large prime on $B_\pm$ has logarithm at most $2\log Y+O(1)$. Its expected contribution from (1) first gives $\Pr(B_\pm)\ge1/2-o(1)$. Hence at least $1/6-o(1)$ of the pairs have both $B_\pm$ and $Z>0$. On those pairs the large prime has logarithm at most $(2-1/40)\log Y+O(1)$, since its cofactor contains a distinct prime exceeding $Y^{1/40}$. Subtracting this saving from the upper bound for the expected contribution gives
$$1-o(1)\le2\Pr(B_\pm)-\tfrac1{40}(\tfrac16-o(1)).$$
It follows that $\Pr(B_\pm)\ge1/2+1/480-o(1)$ for both signs individually. Inclusion–exclusion now gives the claimed simultaneous lower bound without an independence assumption.

## The fourth prime and minimality

For all but $o(W)$ of the products $n$, the interval $[n-\lfloor X^{1/3}\rfloor,n]$ contains a prime $r$. Apply Guth–Maynard's Corollary 1.4 to the integer starts $n-\lfloor X^{1/3}\rfloor$, discarding $O(X^{1/3})$ at the lower dyadic boundary. The exceptional count $O(X\exp(-c(\log X)^{1/4}))$ is $o(W)$. Then $n+1-r<\min(p,q)$ and $2r>n+1$ for large $X$. [Guth and Maynard](https://arxiv.org/pdf/2405.20552)

A prime $b$ divides every coefficient precisely when the sum of the base-$b$ digits of the weight is at most two. Otherwise split its digit units among three nonempty parts to obtain a carry-free row. If no universal prime exists, no list of at most three primes can certify the weight: choose a row avoiding each listed prime.

At $n=pq$ an odd base cannot give digit sum two because $n$ is odd, and digit sum one would make it a prime power. Digit sum at most two in base 2 would require $n-1$ to be a power of two, excluded by $s_-$. At $n+1$, discard $O(\sqrt X(\log X)^2)$ weights as follows. Base 2 gives only $O((\log X)^2)$ sums of at most two powers. In an odd prime base $b$, an even weight of digit sum at most two has form $b^i+b^j$, $0\le i\le j$. The case $i=0<j$ would make $pq=b^j$, and $i=j$ is impossible since $4\mid n+1$. For $1\le i<j$ we have $b\le\sqrt{2X}$ and $O((\log X)^2)$ possible exponent pairs per base. This discard is $o(W)$. Thus for all sufficiently large $X$ at least $W/250$ of the chosen pairs have minimum certificate size four at both weights.

## Finite checks

Run `python3 check_four_prime_pairs.py` in the project directory. The finite test uses trial division, constructs both lists and checks minimality by digit sums at all prime divisors of the even weight's endpoint coefficient. Its direct example also computes every trinomial coefficient by exact binomial products.

| $Y$ | Prime pairs in the two classes | Large factor of $pq-1$ | Large factor of $pq+1$ | Both | Optimal adjacent pairs certified |
|---:|---:|---:|---:|---:|---:|
| 100 | 12 | 5 | 6 | 3 | 3 |
| 1,000 | 598 | 377 | 305 | 186 | 186 |
| 10,000 | 31,313 | 20,024 | 17,491 | 10,966 | 10,966 |

The two example lists at 187 and 188 pass all 5,859 exact coefficients. An independent run of the repository's C++ factorial-valuation verifier checked these lists and all 21,009,459 unordered positive triples at 11,227 and 11,228 with the following input.

```text
11227 103 109 1871 11213
11228 103 109 401 11213
```
