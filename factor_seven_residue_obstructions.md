# Prime divisors and quadratic residues at multiplier seven

**There is a set of odd primes $p$ of relative natural density one such that $1+7p^a q^{2b}$ satisfies the full three-row conjecture for every odd prime $q$ and all positive $a,b$.** The partner prime has no size restriction. A related uniform bound applies to an odd exponent on the varying base. For each odd prime $q$, let $B_q(Y)$ count primes $p\le Y$ at which some $1+7p^a q^b$ fails with $a$ odd and $b\ge1$. Then $\sup_q B_q(Y)=o(\pi(Y))$.

**Each of the following eleven primes $r$ works at every weight $1+7r^a q^b$ for every odd prime $q$ and all positive $a,b$.**

| Order of 7 modulo $r$ | Prime $r$ |
|---:|---:|
| 3 | 19 |
| 6 | 43 |
| 18 | 117307 |
| 54 | 1628413557556843 |
| 66 | 1357105535093947 |
| 78 | 115868130379 |
| 102 | 90920746649833 |
| 114 | 24006520939 and 275430251632699963 |
| 126 | 1506625093 |
| 162 | 3426171911854086574033958308360873411 |

## A local obstruction when a prime divides the cofactor

If $1+7t$, with $t=p^a q^b$ odd, fails, then
$$t\equiv11\pmod{72},\qquad 6(t+1)=\sum_{i=1}^6 7^{e_i},\quad e_i\ge1.$$
This follows from the [even-weight reduction](predecessor_five_seven.md#full-factor-seven-at-even-weights), including its treatment of equal or overlapping prime bases. Thus the sum is 0 modulo 72, and is 6 modulo each prime divisor $r$ of $t$. Modulo 8, three of the six exponents must be even and three odd.

For $r=19,43$, the exponents in both congruences can be reduced modulo 6. Exactly 33 multisets of six residues give sum 0 modulo 72. Their sums modulo 19 belong to $\{5,8,10,12,13,15,16,17,18\}$, and their sums modulo 43 belong to $\{0,2,3,5,8,10,12,13,14,16,17,18,21,22,25,26,27,29,30,31,33,35,38,40,41\}$. The residue 6 is missing from both sets.

For the other entries let $L=\operatorname{lcm}(6,\operatorname{ord}_r(7))$. The powers of 7 modulo $72r$ have period dividing $L$. Form the set of sums of three powers with even exponents modulo $L$. For each sum of three powers with odd exponents, its complement to the unique target modulo $72r$ misses that set. Run `python3 check_seven_residue_obstructions.py` to check every entry with the standard library. The same program verifies the exact orders and proves primality by trial division or recursive Lucas certificates.

## Quadratic characters with a shared exponent period

Take $h$ to be the product of the first $J$ odd primes. For each $d\mid h$ with $d>1$, choose a prime $r_d$ for which 7 has exact order $d$. These primes are distinct, odd, and different from 3 and 7. Here is an elementary existence argument. The cyclotomic factor $\Phi_d(7)$ has absolute value at least $6^{\varphi(d)}>d$, using $\ell-1\ge\sqrt\ell$ at every odd prime $\ell\mid d$. For $r\mid d$ dividing $\Phi_d(7)$, Fermat's congruence implies $r\mid 7^{d/r}-1$. The quotient $(7^d-1)/(7^{d/r}-1)$ has $r$-adic valuation one by binomial expansion, and it contains the factor $\Phi_d(7)$. Thus all prime factors shared with $d$ contribute at most $d$ in total. Any remaining prime $r\mid\Phi_d(7)$ has order exactly $d$, since the factorization of $X^d-1$ modulo $r$ is squarefree.

Write $k=2^J-1$ for the number of chosen primes. Folding six exponents modulo $h$ produces at most
$$K=\binom{h+5}{6}$$
possible vectors $(t\bmod r_d)_d$, since 6 is invertible at each prime. At each nonzero entry we know the Legendre symbol of $t$.

In the first theorem a square $t$ is impossible since $t\equiv3\pmod8$, so we may assume $a$ odd. Remove the finitely many $p$ among the $r_d$. At every coordinate other than possibly $r_d=q$, the sign of $p$ equals the sign of $t$ because the exponent on $q$ is even. An admissible vector can have at most one zero entry. Its location is determined by the vector, and at most two sign vectors for $p$ extend it. This accounts also for partners $q$ that are among the auxiliary primes. Altogether at most $2K$ of the $2^k$ sign vectors can occur at an exceptional $p$.

For the second assertion fix $q$ and let $b$ have either parity. Outside the coordinate $r_d=q$, if present, the vector of signs for $p$ is determined by that of $t$ multiplied coordinatewise by the signs of $q$ to the power $b$. Two parities and at most one undetermined sign give at most $4K$ eligible sign vectors for $p$, uniformly in $q$.

By the Chinese remainder theorem, each sign vector occupies exactly the same number of reduced classes modulo $\prod_d r_d$. The prime number theorem in arithmetic progressions for this fixed modulus gives relative prime density $2^{-k}$ for each. For fixed $J$ there are finitely many vectors, so its error after normalization is $o_J(1)$ also for the supremum over $q$. The upper densities in the two claims are at most $2K/2^k$ and $4K/2^k$, respectively. Bertrand's postulate gives $\log h=O(J^2)$, whereas $k=2^J-1$. Sending $J$ to infinity proves both theorems. Already $J=8$ gives $h=111546435$, $k=255$, and $2K/2^k<10^{-31}$.
