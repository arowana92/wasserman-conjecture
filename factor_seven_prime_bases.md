# Two unrestricted exponents at multiplier seven

**For every pair of odd primes $p,q\le47$ and all positive integers $a,b$, the weight $1+7p^a q^b$ satisfies Wasserman's three-row conjecture.** The exponent bounds used in the verification depend on the prime bases. A version with explicit dependence on both bases gives the following further conclusions.

* Among all unordered pairs of odd primes $p,q\le Y$, at most $O((\log Y)^{72}(1+\log\log Y)^{36})$ pairs could fail at even a single pair of exponents.
* For all odd primes $q$ outside a set of relative natural density zero, **every** odd prime $p\le\exp(q^{1/37})$ works with $q$ at **every** pair of positive exponents. Here and below the density of a set of primes is relative to all primes.

## An explicit upper bound for any exceptional weight

Write $t=p^a q^b$. Equal bases or a base equal to 7 give at most two distinct predecessor primes and are already covered. Otherwise the [even-weight reduction](predecessor_five_seven.md#full-factor-seven-at-even-weights) forces $t\equiv11\pmod{72}$, $t\equiv6\pmod7$, and requires $6t$ and $t+1$ to add without a carry in bases $2,3,7,p,q$. In particular neither prime can be 3. The same reduction gives
$$n=6(t+1)=\sum_{i=1}^6 7^{h_i},\qquad h_1\ge\cdots\ge h_6\ge1.$$
Put $X=6t$, $E=1+\lfloor\log_2 X\rfloor$, $d_i=h_1-h_i$, and $F=\sum_{i=1}^j7^{h_i}=7^{h_j}R$. Thus $R\le6\cdot7^{d_j}$. Since $F\ge n/6$, and the remaining sum is at most $6n7^{-d_{j+1}}$, both $|X/F-1|$ and $|\log_7(X/F)|$ are at most
$$128(2^{-E}+7^{-d_{j+1}}).\tag{1}$$
For the logarithm use $|\log X-\log F|\le|X-F|/\min(X,F)$ and $X<n\le2X$. Omit the gap term at a full prefix. Also $X\ne F$ since $7\mid F$ and $7\nmid X$.

Take $A=10^{16}(1+\log p)(1+\log q)$, with natural logarithms. If both primes are at most 47 one can use $A=10^{16}$. Set $C=A(1+\log(E+1))$. Apply the rational form of Matveev's bound to $p,q,7,6/R$ with exponents $a,b,-h_j,1$, omitting the last rational if it equals one. All the exponents have magnitude at most $E+1$ and the height of $6/R$ is bounded by $\log\max(6,R)<2d_j+2$. In the small-base case,
$$1.4\cdot30^7\cdot4^{4.5}(\log47)^2\log7<10^{15}.$$
In the general case the same coefficient with $(\log p)(\log q)$ in place of $(\log47)^2$ is less than $A/10$. Hence $\log|X/F-1|>-(C/10)(2d_j+2)$. [Matveev's bound, Theorem 2.1](https://archive.mpim-bonn.mpg.de/id/eprint/5025/1/mpim-preprint_2023-34.pdf)

Comparing with (1) gives $E<C(d_j+2)$ or $d_{j+1}+2<C(d_j+2)$. Indeed the smaller of $E\log2$ and $d_{j+1}\log7$ is below $\log256+(C/10)(2d_j+2)$. At the full prefix only the bound on $E$ remains. Starting at $d_1=0$, iteration gives
$$E<2C^6,\qquad E<[100A\log(2A)]^6.\tag{2}$$
For the second inequality, set $H=[100A\log(2A)]^6$. At $A\ge10^{16}$ we have $1+\log(H+1)<20\log(2A)$, which gives $2[A(1+\log(H+1))]^6<H$. The ratio $(1+\log(E+1))^6/E$ decreases thereafter. For $A=10^{16}$ the same comparison at $10^{130}$ gives $E<10^{130}$.

## Counting prime bases

At $p,q\le Y$, (2) gives $E\ll (\log Y)^{12}(1+\log\log Y)^6$. All six exponents $h_i$ are at most $E+1$. There are $O((1+E)^6)$ multisets with this bound, and each multiset fixes $t$. Unique factorization recovers its unordered pair of distinct prime bases. This proves the first counting assertion. For fixed $p$, the same reasoning bounds the number of possible exceptional partners $q\le Y$ by $O_p((\log Y)^{36}(1+\log\log Y)^{36})$.

For $Y\le q<2Y$ and $p\le\exp(q^{1/37})$, we have $A\ll Y^{1/37}\log Y$ and $E\ll Y^{6/37}(\log Y)^{12}$. Consequently at most $O(Y^{36/37}(\log Y)^{72})$ primes in that interval have even one possible exception with such a partner. Each multiset supplies at most two choices for $q$. Summing over dyadic intervals and using the prime number theorem proves the second assertion. The same proof works with any fixed exponent less than $1/36$ in place of $1/37$.

## Finite certificate for the primes through 47

The two residue conditions combine to $t\equiv83\pmod{504}$. The remaining bases have sixth power one modulo 504. Trying the exponent classes modulo 6 leaves just $(13,23),(17,43),(29,31),(37,47),(41,43)$.

For each pair, suppose $3000<E\le M$. Write $\alpha=\log_7 p$, $\beta=\log_7 q$. The verifier obtains positive denominators $D$ and integers $u,v$ with both $|D\alpha-u|$ and $|D\beta-v|$ less than $1/[128(M+3000)]$. There are two ways to bound a prefix from below.

* If $\|D\log_7(R/6)\|>1/32$, where the norm denotes distance to an integer, replacing $D\alpha,D\beta$ by $u,v$ costs less than $1/64$ after multiplication by $a,b$. Thus $|\log_7(X/F)|>1/(64D)$.
* If $R/6=p^i q^k7^w$ for integers $i,k,w$, these exponents are nonnegative because $6$ is coprime to the three bases. The verifier checks $i,k\le3000$ and uses three linearly independent rows $(u,v,D)$ satisfying the same error bounds. Their scalar products with $(a-i,b-k,-h_j-w)$ cannot all vanish. One is a nonzero integer, so $|\log_7(X/F)|>1/(64D_{\max})$ for the largest denominator in the three rows.

Let $D_*$ be the denominator in the resulting bound. At every pass $16384D_*<2^{3000}$. Thus (1) rules out a full prefix and bounds the next gap by $d_{j+1}<g=\lceil\log_7(16384D_*)\rceil$. Enumerate every next gap from $d_j$ to $g-1$, using the new coefficient $7^{d_{j+1}-d_j}R+1$. At each node a necessary modular test discards coefficients that cannot give $n\equiv72\pmod{432}$ and $n\equiv6\pmod{6pq}$. The test starts with all multiples of the target by powers of seven, then at each shorter length subtracts one and closes under multiplication by seven. This includes every possible earlier coefficient because seven is invertible modulo both moduli.

| Bound $M$ | Last prefix length $J$ | New bound |
|---:|---:|---:|
| $10^{130}$ | 2 | $10^{82}$ |
| $10^{82}$ | 3 | $10^{62}$ |
| $10^{62}$ | 4 | $10^{42}$ |
| $10^{42}$ | 5 | $10^{22}$ |
| $10^{22}$ | 6 | 3000 |

At $J<6$, if $G$ is the largest $g$ at length $J$, the untreated prefixes give $E<(G+2)C^{6-J}$ by the recurrence in (2). The verifier checks this with $C<10^{16}(2+3\log_{10}M)$. At $J=6$ all full prefixes have been excluded. Finally enumerate the prime powers at $E\le3000$ exactly. The numbers of exponent pairs satisfying $t\equiv83\pmod{504}$ in the above order are $7414,5665,5218,4314,12877$. Every one has a carry in $6t+(t+1)$ in base two.

Run the following with python-flint installed.

```sh
python3 check_seven_two_bases.py
python3 check_seven_two_bases_cross.py
```

The first program uses Arb balls. The cross-check uses LLL proposals from the first program, but independently certifies the approximations with integer logarithm intervals and determinants, and checks all prefixes using the joint modulus $432pq$. It also uses an independent sieve on exponents. The carry test at layer $q^{b+1}$ requires $z+(6z\bmod q)<q$ for $z=p^a\bmod q$, and likewise with $p$ and $q$ exchanged. These tests exclude $(13,23)$ and $(29,31)$ outright. For the finite range the cross-check adds binary digit sums. [Integer logarithm enclosures](six_sevens_near_two.md#reproduction)
