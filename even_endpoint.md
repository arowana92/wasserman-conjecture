# An odd divisor at even weights below $2^{64}$

**Finite theorem.** Let $4\le N<2^{64}$ be even. Write $N=mP$ with $m$ odd and $P$ a power of two. If three positive multiples of $m$ sum to $N$, their trinomial coefficient has a prime factor in common with $N-1$.

This excludes a route to a counterexample at all these weights. Suppose one selects the endpoint $E=(1,1,N-2)$ and a row $A$ avoiding 2. If $A$ has three powers of two as its parts, it cannot avoid any odd prime of $N$ or $N-1$. Indeed its residues at any odd prime are all nonzero, whereas a carry-free row at a divisor of $N$ or $N-1$ has residue sum 0 or 1. Thus the third row would have to avoid every prime of $m$ and $N-1$. Avoiding the primes of $m$ forces all its parts to be multiples of $m$, and the theorem excludes that possibility. The same finite theorem excludes any proposed third row avoiding *all* odd primes of $N(N-1)$, for any choice of $A$.

## Arithmetic reduction

Write the proposed row as $(m y_1,m y_2,m y_3)$ with $y_1+y_2+y_3=P$. Put $M=N-1$. If the row avoids all primes of $M$, at each full prime power $w\parallel M$ the three actual residues are $1,0,0$ in some order. Since $m$ is invertible modulo $M$, two normalized parts are divisible by $w$. Consequently
$$M^2\mid y_1y_2y_3\le(P/3)^3.$$
This rejects $P=2$ and bounds $m$ by $\sqrt{P/27}+1/P$. It leaves 1,902,253 weights below $2^{64}$.

Assign each full prime power of $M$ to the exceptional position where its residue is 1. Let $d_i$ be the product assigned to position $i$. The positive part $y_i\le P-2$ is divisible by $M/d_i$. Hence $d_i\ge M/(P-2)>m$. All three groups are nonempty and each full power is at most $M/(m+1)^2$. For each group the actual part is the unique integer in $(0,M)$ equal to 1 modulo $d_i$ and 0 at the other factors. Its divisibility by $m$ and the required sum $M+1$ are exact tests.

The program computes these residues by summing the CRT elements for individual powers. It retains every partition into three nonempty groups satisfying the tests above. It then uses Legendre’s formula to check whether each candidate avoids every prime of $M$. The row count includes each unordered row once.

| Stage | Count |
|---|---:|
| Weights passing $27M^2\le P^3$ | 1,902,253 |
| CRT rows divisible by $m$ and summing to $N$ | 18,346 |
| Such rows at weights with at least three binary ones | 6 |
| Rows avoiding all primes of $M$ | 0 |

The six rows at binary-eligible weights can also be rejected by direct evaluations of Legendre’s formula. Each pair $(k,m)$ listed below produces exactly one such unordered CRT row, with $P=2^k$.

| $k$ | $m$ | Prime divisor of $N-1$ | Valuation of the row coefficient |
|---:|---:|---:|---:|
| 26 | 7 | 3 | 13 |
| 30 | 7 | 3 | 21 |
| 39 | 607 | 5 | 17 |
| 43 | 7 | 5 | 21 |
| 45 | 25 | 17 | 10 |
| 50 | 13 | 3 | 33 |

## Reproduction

Use Python 3.9 or later with python-flint installed. The search used python-flint 0.6.0.

```sh
python3 search_even_endpoint.py 64
```

For efficiency the script strips primes through 3,000,000 simultaneously from $2^k m-1$ in each progression of odd $m$. A residual $r>1$ below the square of the sieve limit is prime. If this prime exceeds the group bound, or the product of all other factors is at most $m^2$, three groups greater than $m$ are impossible. A residual below the cube of the sieve limit has at most two prime factors with multiplicity. In that case the remaining known product must exceed $m$, and the residual cannot exceed the square of the group bound. For other factorizations the checker confirms primality of the returned bases and multiplies the full powers back to the input. [FLINT integer primality](https://flintlib.org/doc/fmpz/fmpz_is_prime.html)
