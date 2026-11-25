# A digit test for three prime powers in the predecessor

Let $F\in\mathbb Z[X]$ have nonnegative coefficients, positive constant term $b$, coefficient maximum $H$ and content 1 (the gcd of its coefficients). Let $P=p^a$, $p$ prime, $a\ge1$. **If $F(P)$ has at most two distinct prime divisors, both $P F(P)+1$ and $2P F(P)+1$ satisfy Wasserman's full three-row conjecture whenever $P>4b(2b+1)H$.** This threshold is independent of the degree. The two prime divisors may vary with $P$.

There is a second consequence for two arbitrary nonconstant polynomials $f,g$ with nonnegative integer coefficients. For every sufficiently large prime power $P$ for which $f(P)$ and $g(P)$ are prime powers, the conjecture holds at
$$N=P f(P)g(P)+1.$$
The threshold is effective and depends on the two polynomials. This includes unequal degrees, unbounded exponents and nonprimitive polynomials.

## A gap in the possible gcds of a part

Write $N=PD+1$, with $D\ge1$, and expand $D=\sum_i d_iP^i$. Put $c=d_0>0$ and $H_D=\max_i d_i$. Consider a positive trinomial row of weight $N$ whose coefficient is coprime to $p$. Exactly one part $B$ is 1 modulo $P$ and the others are 0. Set
$$s=D/\gcd(D,B).$$
If $s(c+1)H_D<P$, then all $d_i$ are multiples of $c$. Put $F(X)=\sum_i(d_i/c)X^i$, so $F(0)=1$. For an integer $0\le y\le c-2$ we have
$$B=(1+yP)F(P),\qquad s=c/\gcd(c,1+yP),\qquad F_{i+1}\le(c-y)F_i\quad(i\ge0). \tag{1}$$
Here coefficients beyond the degree are zero. If a divisor $U>1$ of $N$ divides all parts, then $U\mid 1+yP$ and $y\ge1$. Thus a counterexample requiring this row and a small $s$ needs $c\ge3$, with every digit from the units through the leading place positive. It must also have $c$ equal to the gcd of the digits.

**Proof.** Kummer's criterion grouped into blocks of $a$ digits shows that the base-$P$ digits of $B$ are bounded position by position by those of $1+PD$. Denote the digit polynomials of $D,B$ by $D_*(X),B_*(X)$. We have $B_*(0)=1$. The integer $k=sB/D$ satisfies $0<k<sP$, since $B\le N-2$. Reducing $kD=sB$ modulo $P$ gives $ck=s+jP$. Here $0\le j<sc$, since $s<P$. Consequently
$$cs B_*(P)=D_*(P)(s+jP).$$
The coefficients on the left are at most $csH_D$. On the right they are at most $(s+j)H_D<s(c+1)H_D$. Both sides give ordinary base-$P$ expansions, proving the polynomial identity. By Gauss's lemma the primitive part of $D_*$ divides $B_*$ over $\mathbb Z[X]$. Its positive constant term must divide 1, so the content of $D_*$ is $c$.

The quotient $B_*/F$ has degree at most one, constant term 1, and a nonnegative linear coefficient $y$. The strict inequality $B<PD$ implies $y<c$. Comparing coefficients in $B_*\le 1+cXF$ position by position gives the last inequality in (1). If $c-y=1$, the positive coefficients of $F$ start at 1 and weakly decrease, so $F=1+X+\cdots+X^m$, including $m=0$. The other two positive parts then sum to $1+(P-1)F(P)=P^{m+1}$. They cannot add without a carry in base $p$. This proves $y\le c-2$. The identity for $s$ follows by cancelling $F(P)$ in the gcd. Finally $F(P)\mid N-1$ and $U\mid B$, giving $U\mid1+yP$.

## Locating the part

Use the [all-weight endpoint reduction](all_weights_endpoint_reduction.md). With exactly three distinct predecessor primes, its [three-factor consequence](endpoint_factors.md) gives one row $A$ avoiding all three. Each full power is assigned its own exceptional position in $A$, and some divisor $U\ge3$ of $N$ divides all parts of $A$. At a full power $P$ the exceptional part $B$ therefore has $D\mid B$, giving $s=1$. This yields a direct digit test on any of the three full powers. If $P$ is a proper power of the same base dividing $N-1$, then $B$ still contains the entire part of $D$ coprime to $p$. In fact $s=p^{v_p(D)}$, since $B\equiv1\pmod p$.

There is also a four-factor case. Suppose $N-1=2PQR$, with $P,Q,R$ full odd powers of distinct primes. The [four-bin inequality](two_prime_powers.md) forces the three odd bases onto the same row $A$. They have different exceptional positions. Otherwise there would be a position divisible by $PQR$, and also by $U\ge3$ coprime to $PQR$, although each part is less than $2PQR$. The part exceptional at $P$ therefore contains $QR$, giving $s\le2$ for $D=2QR$.

For the first theorem take $D=dF(P)$ with $d=1,2$. Its base-$P$ digit polynomial is $dF(X)$ under the displayed threshold. In the cases with at most two predecessor primes the [two-factor result](endpoint_reduction.md) applies immediately. With three, $s=p^{v_p(D)}\le db$ since $D\equiv db\pmod P$ and $0<db<P$. With four, necessarily $d=2$, $p$ is odd, and $F(P)$ has two odd bases distinct from $p$, so $s\le2$. In every case requiring a test, $s(c+1)H_D\le (2b)(2b+1)(2H)$. But the digit lemma would require the content of $dF$ to equal $c\ge3$, contradicting its content $d\le2$.

More generally, at a three-factor predecessor with $P$ full and $(c+1)H_D<P$, a zero between two nonzero base-$P$ digits of the cofactor suffices to prove the full conjecture. For example, the conclusion holds at $1+P(vP^j+c)$ for all $j\ge2$ and positive $v,c$ with $(c+1)\max(v,c)<P$, $p\nmid c$, whenever $vP^j+c$ has at most two distinct prime divisors.

A useful case has no large-$P$ threshold. Let $Z$ have at most two distinct prime divisors and satisfy $Z\equiv1\pmod P$. Then $PZ+1$ satisfies the conjecture for every prime power $P$. If the base-$P$ digits of $Z$ are all at most $(P-1)/2$, so does $2PZ+1$. In particular both weights work at every prime power when $Z$ is a polynomial in $P$ with coefficients 0 and 1 and constant term 1, subject only to the same prime-factor condition on $Z$.

Indeed, at $PZ+1$ in the three-factor case the special part is a multiple of $Z$ less than $PZ$, and its multiplier is 1 modulo $P$. This forces the part to be $Z$, coprime to $N$. At $2PZ+1$ the preceding three- and four-factor arguments make $B$ a multiple of $Z$ as well, since $p\nmid Z$. Its only possible multipliers are 1 and $P+1$. The latter gives $B=(P+1)Z$. Under the digit bound, neither $Z+PZ$ nor $1+2PZ$ has any base-$P$ digit overflow in these formulas. Comparing digits of the part to the total forces each digit of $Z$ to be at most the preceding digit. They start at 1, giving a repunit and a complement equal to a single power of $P$, impossible for the other two positive parts. For instance, for each fixed pair of coprime prime powers $P,Q$, Dirichlet’s theorem gives infinitely many primes $r\equiv Q^{-1}\pmod P$. Every weight $PQr+1$ satisfies the conjecture by the first assertion. For the 0–1 polynomial statement at $P=2$, the predecessor of $2PZ+1$ has full power 4 and at most two odd bases, so the [small full-power theorem](three_factors_small_power.md) applies.

## Two polynomial prime-power factors

It suffices to consider three distinct bases. If $f(0)$ or $g(0)$ vanishes then its value is divisible by $p$ and cannot supply a different prime base. Otherwise write $f=h f_0$, $g=k g_0$, taking positive contents $h,k$ and primitive parts. For large $P$, $D=f(P)g(P)$ has digit polynomial $fg$, and $s=1$. The digit lemma gives $f_0(0)g_0(0)=1$ and hence both constant terms are 1. Write $c=hk$. The part exceptional at $P$ is $(1+yP)f_0(P)g_0(P)$ with $1\le y\le c-2$ and $U\mid1+yP$. Also $c\mid1+yP$ by $s=1$. If both contents are 1 the contradiction is immediate.

A content $h>1$ forces $f(P)$ to be a power of a fixed prime dividing $h$. A fixed integer polynomial with at least two distinct complex roots has greatest prime divisor tending effectively to infinity on positive integer inputs. Thus, past an effective bound, $f_0$ can only have a single distinct root. Its rational coefficients make that root rational, and primitivity, nonnegative coefficients and constant term 1 give $f_0=(1+uX)^e$ with integers $u,e\ge1$. The same argument applies to $g_0$ if $k>1$. If both contents exceed 1, the prime factors of $(fg)(P)$ lie in a fixed finite set. Two distinct roots of $fg$ again contradict the greatest-prime-divisor bound. A shared single root makes both values divisible by $1+uP>1$, contrary to distinct bases. [Shorey and Tijdeman, Theorem 1 with interval length one](https://www.numdam.org/item/CM_1976__33_2_187_0/)

It remains to consider $h>1,k=1$, interchanging $f,g$ as necessary. Here $c=h$. Since $h(1+uP)^e$ is a power of a fixed prime, $h\mid1+uP$ for large $P$. Thus $u\equiv y\pmod h$ by $h\mid1+yP$ and $\gcd(P,h)=1$. The coefficient of $X$ in $f_0g_0$ is at least $eu$, and (1) bounds it by $h-y$. Hence $1\le u,y<h$ and $u=y$. But $U$ divides $1+uP$, which divides $N-1$, a contradiction.

## Reproduction

The verifier checks 354 instances inside the gcd gap, including 223 with $s>1$ and 25 with a part sharing a divisor of at least 3 with the weight and $y>0$. It also tests 9,307 carry-free rows obtained from exact binomial products, and 1,152 CRT divisor checks. Run the controls using Python.

```sh
python3 check_cofactor_digit_gap.py
```
