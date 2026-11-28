# All full endpoint powers and arbitrarily deep tests at any fixed primes

Fix any finite set $\mathcal P$ of primes and any integer $k\ge1$. There is an arithmetic progression of weights $N$ with positive rows $A,C$ at every weight in the progression such that

1. **Every** $p\in\mathcal P$ divides $N(N-1)$ to exponent exactly one. For each of them at least one of $A,C$ passes its first $1+k$ base-$p$ digits.
2. For **every** prime power $p^e\parallel N(N-1)$, at least one of the two rows passes through digit $e$. Here the test includes primes outside $\mathcal P$, even when $e$ is large.
3. For any prescribed prime $R>52154$ larger than every member of $\mathcal P$ the progression can be chosen so that $R$ divides the endpoint coefficient and both row coefficients. The first row carries by digit one and the second by digit two at $R$.

The units digit is digit one. Thus full endpoint powers, supplemented by any fixed finite choice of deeper prime tests, miss a common prime on a set of weights of **positive natural density**. Each prime chosen for deeper testing actually occurs in the endpoint. This contrasts with the [density-one theorem using a varying set of bases](bounded_carry_depth.md).

There are also two reduced progressions of primes $s,q$ for which the **full three-row conjecture** holds at all weights $1+52154\,s q^b$, for every $b\ge1$, and the same pair of row formulas has properties 1–3. Both primes vary independently. The set of all weights in the first paragraph is an ordinary arithmetic progression. The prime-parameter subfamily has this additional all-rows certification.

## The row formulas

Put $d=26077=89\cdot293$, $N=1+dW$, $V=N/7$, and use
$$
 A=(9345W,16408W,324W+1),\qquad C=(V,3V,3V).
$$
The slopes of $A$ sum to $d$. At the following positive scalars, the **entire** row $A$ is carry-free in the indicated base. The table gives the sums in successive digit columns, starting at the units digit.

| Prime | $W$ | Column sums |
|---:|---:|---|
| 2 | 2 | The binary supports of $9345,16408,324$ are disjoint. Multiplication by 2 leaves the units place for the added 1. |
| 7 | 108 | $0,6,5,6,3,6,2,3$ |
| 89 | 9587 | $1,62,55,87,3$ |
| 293 | 1654 | $1,120,209,1$ |

Adjoin these four primes to $\mathcal P$. At each require $W$ congruent to its tabulated scalar modulo $p^{k+1}$. The table and the residues modulo $p^2$ show $7\parallel N$ and $p\parallel N-1$ for $p=2,89,293$. At each other $p\in\mathcal P$ impose
$$
                  N\equiv7p m_p\pmod{p^{k+1}},\qquad
 m_5=2,\quad m_p=1\ (p\ne5).
$$
This gives $p\parallel N$ and the tested digits of $C$ agree with those of $p(m_p,3m_p,3m_p)$. The latter row is carry-free at 3 with $m_3=1$ and at 5 with $m_5=2$. At all remaining primes here the normalized sum is $7<p$. These conditions each give one class for $W$ because $p\nmid d$.

Also impose $N\equiv R\bmod R^2$. All moduli are pairwise coprime, so CRT gives one class $W_0\bmod H$. Positive integers in that class give the progression of weights. At any prime dividing $W$ the first row has residues $(0,0,1)$ at its full power in $M=N-1$, because $(W,d)=1$. At 89 and 293 use the table. At every weight prime except 7, the entire corresponding power divides $V$. At 7 the table applies. This proves the full-power assertion. At $R$ the weight is zero modulo $R$ while the first part of $A$ is a unit. In $C/R$ all three parts are units, but their sum is one modulo $R$. A carry is necessary there. The endpoint row $(1,1,N-2)$ also carries at $R$ in its units digit.

Finally $W_0\equiv2\bmod4$ and $W_0$ is a unit at every odd prime of $H$. Thus $W_0/2$ is a unit modulo $H/2$. Take primes $s\equiv W_0/2\bmod{H/2}$ and $q\equiv1\bmod{H/2}$, available arbitrarily large by Dirichlet's theorem. Then $W=2s q^b$ lies in the same class for every $b\ge1$. At $N=1+(2d)s q^b$ the full weight power at $R$ is $R$ and $N/R\equiv1\bmod R$. Since $R>2d$ and $s q^b$ has at most two prime divisors, the [weight-prime gap theorem](predecessor_weight_gap.md) proves the assertion for all triples. The same proof works with $s^a$ replacing $s$ for every $a\equiv1\bmod{\varphi(H/2)}$.

Run the verifier with any finite cutoff and positive extra depth. It forms the CRT progression and checks all prescribed powers by Legendre floors on two consecutive representatives. The proof above covers the remaining prime divisors without factoring the resulting weights.

~~~sh
python3 check_fixed_prime_depth.py --prime-limit 101 --extra 256
~~~
