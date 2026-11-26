# Six powers at positions in an arithmetic progression

**Theorem.** Let $p$ be any odd prime and let $a,b$ be positive integers. For any six integers $0\le h_i\le20$, repetitions allowed, Wasserman's three-row conjecture holds at
$$N=1+\sum_{i=1}^6p^{a+bh_i}.$$
For instance the positions can be $a,a+b,a+4b,a+9b,a+15b,a+20b$. The bound 20 applies only to the normalized range. Both $a$ and $b$ are unrestricted.

Shift the least $h_i$ to zero and absorb their common divisor into $b$ if necessary. With at most three distinct positions use the [three-position theorem](seven_digit_units.md). A carry in the six positive units leaves $s_p(N)\le5$, which is covered by the [five-unit theorem](six_digit_gaps.md). It remains to consider $s_p(N)=7$. Put $W=p^a$, $B=p^b$, $m(X)=\sum X^{h_i}$ after normalization, and $H=\deg m$. Here $p\nmid m(B)$ because the units at the least position add freely.

By the [seven-unit reduction](seven_digit_units.md#the-inner-units), a counterexample requires an assigned row with parts $(Wx(B),Wy(B),1+Wz(B))$, where $x+y+z=m$. The row scale $U\ge7$ is odd, divides $x(B),y(B)$, and is coprime to $Bm(B)$. Its assigned odd predecessor factor $S$ divides the lcm of $\gcd(m(B),x(B))$ and $\gcd(m(B),y(B))$. The companion requires
$$27m(B)^2\le U^3S^2. \tag{1}$$
Write $x,y$ for the two unshifted sums in an order such that $x$ contains two units. Then $y$ contains two or three units. Both need at least two different exponents because $U\ge7$. Write $f,k$ for $x,y$ with their lowest monomial powers removed.

## Two unshifted pairs

Suppose the respective gaps are $d,e$ in normalized units, and $g=\gcd(d,e)$. Both $d/g,e/g$ must be odd, and $H<d+e-g$ by the [pair-gap bound](seven_pair_gaps.md). Define
$$Q(X)=\frac{X^d+1}{X^g+1}\frac{X^e+1}{X^g+1},\qquad F=B^g+1.$$
Then a further necessary inequality is
$$m(B)<B^g\gcd(m(B),Q(B)). \tag{2}$$
In fact $S\mid FQ(B)$ since the two pair binomials have gcd $F$. Thus $S/\gcd(S,F)$ divides both $m(B)$ and $Q(B)$. In the shifted triple $1+p^r+p^s$, $U$ divides $1+p^r+p^s$ and $F$. Also $\gcd(S,F)$ divides $p^{|r-s|}+1$ because both unshifted pairs vanish there and the whole row sums to $1\bmod S$. Applying the binomial–trinomial gcd inequality at $n=bg$ gives $U^3\gcd(S,F)^2<27B^{2g}$. Now (1) gives (2). The integer $Q(B)$ is odd, since both quotients have an odd number of odd summands with alternating signs.

## One unshifted pair and one unshifted triple

For each polynomial pair $(A,D)$ used below take its positive-leading gcd $h$ over $\mathbb Z[X]$. The verifier checks an integral Bézout identity for the cofactors with positive constant $c$. At $B\ge3$ this gives
$$\operatorname{odd}(\gcd(A(B),D(B)))\mid \operatorname{odd}(c)\,h(B).$$
Here $\operatorname{odd}$ denotes the largest odd divisor. Every gcd polynomial here is positive on $[3,\infty)$, since it divides a polynomial with nonnegative coefficients and has positive leading coefficient.

Use this construction for $(f,k)$, obtaining $u,C_u$ with the odd part of the constant written $C_u$. Remove from $u$ every polynomial factor it shares with $m$, repeatedly taking gcds. This preserves $U\mid C_u u(B)$ since $U$ is coprime to $m(B)$. For $(m,f)$ and $(m,k)$ obtain $h_x,C_x$ and $h_y,C_y$. Set $h=\operatorname{lcm}(h_x,h_y)$ with positive leading coefficient. Then $h\mid m$ in $\mathbb Z[X]$, and $S\mid C_x C_y h(B)$ since $h(B)$ is a multiple of each of $h_x(B),h_y(B)$. Hence (1) is impossible if
$$27(m/h)^2(B)>C_u^3(C_xC_y)^2 u(B)^3. \tag{3}$$

## Finite verification for unbounded bases

For two pairs the same Bézout construction on $(m,Q)$ bounds the gcd in (2) by $C h(B)$. The polynomial tested for positivity is $m/h-CX^g$. In the pair-and-triple case it is the difference in (3). For each allocation that reaches the polynomial test, the verifier finds a $T$ such that translating $X$ to $T+X$ gives nonnegative coefficients and a positive constant. This certifies positivity for every $B\ge T$. For each odd $3\le B<T$ it checks (2) directly or checks (1) using integer gcd and lcm upper bounds. This covers every odd prime power $B$.

Run `python3 check_seven_progressions.py` with python-flint installed. There are 48,311 normalized multisets on four to six distinct levels and 3,110,335 eligible allocations. Of these, 119,453 two-pair and 965,864 pair-and-triple allocations reach a polynomial test. The largest $T$ is 3,072. All 726,632 individual base checks below the thresholds pass. A separate enumeration of the occupied supports and their multiplicities checks the count of multisets. On smaller supports separate controls compare all three-colour assignments of the six labelled units and examine 47,571 positive carry-free rows at 405 weight/base pairs. [Verifier](check_seven_progressions.py)

## An optimal four-prime certificate at 48 digits

Take $N=1+2\cdot5^4+5^{52}+3\cdot5^{67}$ and the four primes
$$5,\qquad q=N/31,\qquad s=(N-1)/(5922\cdot5^4),\qquad r=N-38.$$
No positive row avoids two of them. At $q$, avoiding a carry requires all parts to be multiples of $q$. At $s$, two parts must be multiples of $s$. Satisfying both would require $N\ge q(2s+1)$. At $r>N/2$, the two smaller parts must sum to at most $38$, which is less than $s$ and $2q$. These inequalities exclude all three pairs among $q,s,r$. Two independent enumerations give exactly 68 unordered positive rows avoiding 5. Legendre's factorial valuation at each of $q,s,r$ is positive for every one of these rows. Thus any three rows share a prime from the list.

Every individual prime can be avoided by at least one row: its base expansion of $N$ has at least three digit units, which can be distributed into three positive parts. Indeed at an odd base the digit sum is odd and cannot equal one, since $N=31q$ has two distinct prime divisors. In binary the odd $N$ has at least three ones, since $N-1$ is divisible by 5. A list of at most three primes can therefore be defeated by choosing one row for each prime. This proves optimality.

Run `python3 check_seven_prime_list.py` to reproduce the exact valuations and recursive Lucas primality certificates. [Verifier](check_seven_prime_list.py)
