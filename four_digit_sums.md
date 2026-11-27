# Four digit sums at a possible exception

Write $E=N(N-1)$ and let $K_A,K_B$ denote coefficients of two positive trinomial rows of weight $N$. The following counts concern weights $3\le N\le X$ and allow the choice of rows to vary with $N$.

* There are $O(X^{0.8204})$ weights for which $\gcd(E,K_A,K_B)=1$ with each coefficient indivisible by at least one of 2 and 3. If both are odd the bound is $X^{h(1/4)+o(1)}$, with $h(1/4)=0.811278\ldots$ and $h$ the binary entropy. Consequently, outside $O(X^{0.8204})$ weights, every three trinomials whose coefficients each avoid 2 or 3 have a common prime. This last deduction uses the [all-weight endpoint theorem](all_weights_endpoint_reduction.md).
* Outside $O(X^{0.8561})$ weights, any counterexample has a nonendpoint row $B$ with **all four** inequalities
  $$
  v_2(K_B)>\tfrac1{10}\log_2 N,\quad v_3(K_B)>\tfrac1{10}\log_3 N,\quad
  v_5(K_B)>\tfrac1{20}\log_5 N,\quad v_7(K_B)>\tfrac1{40}\log_7 N .
  $$
  Its part supported on these four primes exceeds $N^{11/40}$. The other nonendpoint row must avoid every one of these primes dividing $E$.
* When each coefficient is indivisible by at least one of 2 and 3, all such pairs satisfy $\gcd(E,K_A,K_B)>N^{1/8}$ outside $O(X^{0.8805})$ weights. If every part is at least two, the bound on exceptional weights improves to $O(X^{0.8511})$. For **three** trinomials with all parts at least two, each coefficient avoiding 2 or 3, all three share a divisor of $E$ exceeding $N^{1/16}$ outside $O(X^{0.9717})$ weights.

## A short digit word identifies the assigned factors

For a collection of coefficients let $G$ be their gcd with $E$. At each $q^a\parallel N$ or $N-1$, retain $q^{a-\min(a,\min_i v_q(K_i))}$ and assign it to a row of minimum valuation. Denote by $U,S$ the products assigned to a given row from $N,N-1$. Each part $z$ satisfies $U\mid z,\ S\mid z(z-1)$, by the [partial-power residue lemma](endpoint_pair_bound.md#retaining-partial-prime-powers). The product of all $US$ across the rows is $E/G$. If $z\ge2$ it determines $U,S$ up to $X^{o(1)}$ divisor choices. They determine a single progression for $N$ modulo $US$.

First suppose there are two rows with $G=1$. No part equals one. Indeed a row containing one has coefficient divisible by $N$, and every positive row has coefficient divisible by $N/\gcd(z_1,z_2,z_3)\ge3$ by decrementing each part in turn. On each dyadic interval at least one assigned modulus $US\ge\sqrt E$ is comparable to or larger than that interval's length. A part of this row therefore determines $N$ up to $X^{o(1)}$ possibilities.

Fix any prime $p$ allowed for this selected row and put $L=\lceil\log_p(2X+1)\rceil$ on $[X,2X)$. With $s_p$ the base-$p$ digit sum, the four nonnegative integers
$$ z_1,\ z_2,\ z_3,\ p^L-1-N $$
have digit sums totaling $(p-1)(L+v_p(K))$ by Legendre's formula and digitwise complementation in the fourth integer. Thus if $v_p(K)\le\alpha\log_p N$, one of the four words has mean digit at most $\mu=(p-1)(1+\alpha)/4$. If it is the fourth, it determines $N$ exactly. The number of such words is at most $p^{L F_p(\mu;u)}$ for every $0<u<1$, where
$$
 F_p(\mu;u)=\log_p(1+u+\cdots+u^{p-1})+\mu\log_p(1/u).
$$
This follows by summing the generating weights $u^{s_p(z)}$. The same count applies to every choice of the two allowed primes. For binary carry-free rows optimize at $u=1/3$; for ternary use $u=0.4343$. The four nonzero valuation thresholds use $u=0.3793,0.48,0.5883,0.6601$ at the four primes. The resulting exponents are strictly below $0.8561$. In every potential counterexample at least one nonendpoint coefficient is odd. It therefore qualifies for a small valuation; if the other also qualifies at any of the four primes the count applies. This proves the first two assertions, summing dyadic intervals.

## Larger common divisors

Suppose a selected row with parts $z_i\ge2$ has $US\gg X^b$, $b\le1$, and avoids $p$. Choose a witness deterministically per weight, give the weights in a stratum the uniform distribution, and use natural entropy $\mathcal H$. Knowing $(U,S)$ leaves at most $O(X^{1-b})$ possible weights. The [shared-label inequality](common_factor_count.md#information-that-every-individual-part-supplies) yields
$$ \mathcal H(N)\le \tfrac13\sum_i\mathcal H(z_i)
     +\tfrac23(1-b)\log X+o(\log X). $$
Indeed each $z_i$ determines the label $(U,S)$ up to $X^{o(1)}$ choices and the three parts together determine $N$.

Split at a threshold $r$ for $s_p(N)/((p-1)L)$. At or above the threshold the complement has mean digit at most $(p-1)(1-r)$, giving exponent $F_p((p-1)(1-r);u)$. Below it the three parts have average mean digit at most $(p-1)r/3$. Entropy against the same digit-generating weights gives exponent $F_p((p-1)r/3;v)+\tfrac23(1-b)$. For a pair with $G\le N^{1/8}$, the larger label gives $b=15/16$. If this row contains a one, $N$ divides its coefficient, and in the other row $N/\gcd(N,K_{\rm other})\ge N/G$ divides every part. Use this weight divisor alone as its label with $b=7/8$. Its parts exceed one for large $N$. For three interior rows with $G\le N^{1/16}$, one label gives $b=31/48$. The verifier checks the following upper bounds with rational logarithm enclosures.

| $b$ | Threshold at 2 | Threshold at 3 | Upper exponent for both |
|---:|---:|---:|---:|
| $15/16$ | .7296 | .7287 | .8511 |
| $7/8$ | .7078 | .7058 | .8805 |
| $31/48$ | .6101 | .6016 | .9717 |

The exact-arithmetic check also verifies the four digit sums for 1,347 unordered rows through 36 and the assignments on all 45,509 unordered pairs of these rows at the same weight.

~~~~sh
python3 check_four_digit_sums.py
~~~~
