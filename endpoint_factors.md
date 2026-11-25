# Restrictions on the two rows paired with $(1,1,N-2)$

Use the [endpoint factor assignment](endpoint_reduction.md). In a counterexample containing $E=(1,1,N-2)$ the other rows are $B=(Ux_i)_{i=1}^3$ and $C=(Vy_i)_{i=1}^3$, where $UV=N$, $ST=N-1$, $\sum x_i=V$, and $\sum y_i=U$. Assign the full prime powers in $U,S$ to $B$ and those in $V,T$ to $C$. Each row avoids the primes assigned to it. All parts are positive, so $U,V\ge3$.

## Lower bound on every part

**Proposition.** Every part $a$ of either $B$ or $C$ satisfies
$$a^5>27(N-1)^2.$$

Take a part $a=Ux$ in $B$. At each prime power $w$ assigned to $S$, this part is 0 or 1 modulo $w$. Because $U$ is invertible modulo $w$, either $w\mid x$ or $w\mid a-1$. Thus $S\mid x(a-1)$. At each prime power assigned to $T$, two of the $y_i$ are multiples of that power. Hence $T^2\mid y_1y_2y_3\le (U/3)^3$. As $U\le a$ and $a\ge3$, we obtain
$$N-1=ST\le \frac{a(a-1)}{U}\left(\frac U3\right)^{3/2}
\le \frac{a(a-1)\sqrt a}{3\sqrt3}<\frac{a^{5/2}}{3\sqrt3}.$$
The same argument applies with the rows interchanged.

## When $N-1$ has three distinct prime factors

Write $M=N-1=w_1w_2w_3$, where the $w_i$ are the full powers of distinct primes. The endpoint assignment forces one of the two rows to take all three factors of $M$, with a different exceptional position for each factor. Indeed, name the rows so that $S\ge V$. Each of the three positive $x_i<V$ must receive an exceptional residue for a factor of $S$, since otherwise $S\mid x_i$.

Up to reordering, this determines the actual parts of $B$:
$$e_i=\frac{M}{w_i}\left(\left(\frac M{w_i}\right)^{-1}\bmod w_i\right),$$
where the inverse is taken in $[1,w_i-1]$. The other two parts are 0 modulo $w_i$, and $e_i$ is 1 modulo $w_i$. Each lies between 1 and $M-1$ because the parts are positive and sum to $M+1$. Thus the weight is excluded unless $\sum_i e_i=N$ and $B$ avoids every prime of $M$.

If the parts sum to $N$, their gcd $G$ is coprime to $M$, so $G$ divides each inverse multiplier in parentheses. In particular $G<\min_i w_i$. The product of the full prime powers of $N$ avoided by $B$ divides $G$ and must be at least 3, because after division by the complementary scale the other row must still have three positive parts. This proves the endpoint case at **every** such weight with $\min_i w_i\le3$, regardless of the other two factors.

**Finite theorem.** The endpoint row is excluded whenever $N-1$ has precisely three distinct prime factors and each of its three full prime powers is at most 20,000. The enumeration covers 2,099,636,978 distinct weights, a nonconsecutive set. The number surviving each successive condition is given below.

| Condition | Count |
|---|---:|
| Three distinct prime bases | 2,099,636,978 |
| CRT parts sum to $N$ | 1,049,832,196 |
| Part gcd at least 3 | 101,932,827 |
| Row avoids all three primes of $N-1$ | 63,327 |
| Avoided full powers from $N$ have product at least 3 | 0 |

The C++ program enumerates triples of powers, caches inverses, and uses base-$p$ carries. At bound 500 the separate Python program computes all three residues directly, factors the full weight at the last stage, and uses factorial valuations. The stage counts agree at 235,224, 117,559, 8,821, 2, and 0.

```sh
c++ -std=c++17 -O3 check_three_factors.cpp -o /tmp/check_three_factors
/tmp/check_three_factors 20000
/tmp/check_three_factors 500
python3 check_three_factors.py
```

## An obstruction to bounding a single row's avoided product

At $N=222{,}788{,}354$, let $B=(86{,}617{,}810,45{,}491{,}401,90{,}679{,}143)$. The factors of $N$ are $2,167,283,2357$ and those of $N-1$ are $29,127,241,251$, all with multiplicity one. The valuations of the coefficient of $B$ are as follows.

| Prime | 2 | 167 | 283 | 2357 | 29 | 127 | 241 | 251 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Valuation | 29 | 4 | 0 | 2 | 0 | 0 | 0 | 5 |

The product of the primes of $N(N-1)$ that $B$ avoids is $283\cdot29\cdot127\cdot241=251{,}191{,}649>N$. Thus a bound of $N$ on that product for an arbitrary single row is false. The Python check above verifies the factorizations and valuations.

This row cannot complete an endpoint counterexample. The other row would have to avoid $2,167,2357$, hence have every part divisible by $787{,}238$. Its three normalized parts would sum to 283. Avoidance of 251 would require two of those positive parts to be multiples of 251, an impossibility.
