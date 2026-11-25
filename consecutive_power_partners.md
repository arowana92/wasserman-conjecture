# A third prime power beside two consecutive prime powers

Let $P$, $P+1$, and $mP+h$ all be prime powers, with positive integers $m,h$. They may have repeated prime bases. Let $d$ be a power of two, including 1, or three times a power of two.

**The full conjecture holds at $N=1+dP(P+1)(mP+h)$ whenever $P>\big(1+2d(m+h)\big)^5$. For every $1\le m,h\le32$ and $d\in\{1,2,3,4,6,8,12,16\}$ it holds at every prime power $P$ satisfying the hypotheses, including those below the cutoff.** The factor $mP+h$ can have arbitrary content.

There is also a bound independent of the multiplier in many binary families. Suppose $P$ is a power of two, put $a=m/\gcd(m,h)$ and $b=h/\gcd(m,h)$, and assume $|b-a|$ has an odd prime divisor. Then **every** $d=2^j$, $j\ge0$, works as soon as $P>4(m+h+1)^3$.

## When the two rows split three prescribed bases

The following degree restriction applies to any fixed $d\ge1$ and any two nonconstant polynomials $f,g$ with nonnegative integer coefficients. Suppose $P,f(P),g(P)$ are powers of three distinct primes, and put $H(X)=1+dXf(X)g(X)$. In a hypothetical endpoint counterexample assign each prime of $H(P)(H(P)-1)$ to an avoiding row. If one of the three prescribed bases is assigned to a different row from the other two, then for sufficiently large $P$, **$H$ has a rational polynomial factor whose degree equals the degree of the singleton polynomial** among $X,f,g$.

Let the two row scales from $H(P)$ be $U$ for the row taking the pair and $V=H(P)/U$ for the row taking the singleton. Each part is a multiple of its scale. Write $W$ for the singleton power and $K,L$ for the pair of powers. A part of the pair row is divisible by $UKL$, while the singleton row has two parts divisible by $VW$. Positivity of the other parts gives
$$2W<U<dW.$$
The first inequality uses the singleton row's normalized sum $U$. The second uses $UKL\le H(P)-2<dWKL$.

The row avoiding the base of $P$ splits the coefficients of $H$ into polynomials with nonnegative integer coefficients, for $P$ larger than all those coefficients. This follows by grouping the base-prime digits into blocks of size $\log_p P$. If $W=P$, one of the singleton row's parts is $P B(P)=PV t$ with $1\le t<d$, where $B$ is chosen from a finite list. Thus $B(P)\mid tH(P)$. Since $V$ has size comparable to $f(P)g(P)$, $B$ has degree $\deg f+\deg g$ on any unbounded sequence with fixed shape and $t$. Polynomial division over the rationals and clearing denominators show $B\mid H$ over the rationals. Indeed a nonzero remainder ultimately has smaller magnitude than the divisor, even after clearing its fixed denominator. The complementary degree is one.

If $W=g(P)$, the pair includes $P$ and $f(P)$. Its part divisible by $UPf(P)$ is $P B(P)=UPf(P)t$ with $1\le t<d$. On any unbounded sequence with fixed choices, the same division argument gives $f\mid B$ as polynomials. The rational polynomial $B/f$ evaluates to $Ut$, is comparable to $g(P)$ in size, and divides $tH(P)$ at the integer inputs in question. Clear its denominator and apply the argument once more to obtain a divisor of $H$ of degree $\deg g$. Interchange $f,g$ for the remaining case. Only finitely many shapes and multipliers were involved. All assignments used here come from the [all-weight endpoint theorem](all_weights_endpoint_reduction.md).

If $H$ is irreducible, there is a simple explicit bound. Set $D=\deg H$ and $B_0=H(1)$. The scale of the row avoiding the base of $P$ divides both $H(P)$ and some nonzero $B(P)$ from a part of the form $P B(P)$. The degrees are at most $D,D-1$ and the sums of the absolute coefficients are at most $B_0$. The resultant is nonzero by irreducibility, and Hadamard's inequality on the Sylvester matrix bounds the scale by $B_0^{2D-1}$. In a split assignment this scale exceeds $P$: it is $U>2W$ if its row takes the pair, and $V>f(P)g(P)$ if its row takes the singleton. This proves the explicit exclusion for $P>B_0^{2D-1}$.

## Applying the restriction to $P,P+1,mP+h$

One of $P,P+1$ contains the base 2. If two of the three prescribed powers share a base, the predecessor has at most two bases, or has exactly three including the full power 3. The [two-factor theorem](endpoint_reduction.md) and the [small full-power theorem](three_factors_small_power.md) handle both cases.

Suppose the three bases are distinct. When $d$ is a power of two there are exactly three predecessor primes, so one row takes them all by the [CRT row theorem](endpoint_factors.md). For $d=3\cdot2^j$ the cubic $H=1+dX(X+1)(mX+h)$ is irreducible. Any rational root must be $-1/t$ with $t$ a positive integer, since the constant coefficient is 1 and all coefficients are nonnegative. It would give
$$t^3=d(t-1)(ht-m).$$
We have $t\ne1$ and $t-1\mid t^3$, giving $t=2$, but $8=d(2h-m)$ is impossible when $3\mid d$. The explicit resultant bound therefore puts all three bases on one row for $P>B_0^5$ with $B_0=1+2d(m+h)$.

On that row the three exceptional positions are distinct. Otherwise one part is divisible by the full powers of all three assigned bases and by its scale $U\ge3$. The three powers have product at least $(N-1)/3$, since the only possible additional base is a full 3. This would make the part exceed $N-2$.

Let $C(P)$ be the part exceptional at $mP+h$. Its polynomial has degree at most 3, nonnegative coefficients of sum at most $B_0$, and $C(0)=0$ because another part is exceptional at $P$. It is divisible by $P+1$, so $C(-1)=0$ for $P>B_0$. Write $aX+b$ for the primitive part of $mX+h$. The congruence $C(P)\equiv1\pmod{aP+b}$ implies that $aP+b$ divides $a^3(C(-b/a)-1)$. The latter integer has absolute value at most $(B_0+1)(m+h)^3<B_0^5$. Hence $C(-b/a)=1$.

Since $X(X+1)\mid C(X)$, this last equality forces $b\mid a^3$ and $|a-b|\mid a^3$, with $a\ne b$. Primitivity gives $b=1$, $a=2$. Thus $2P+1$ divides the prime power $mP+h$ and must itself be a prime power. But $P,P+1,2P+1$ are simultaneously prime powers only for $P=2,3,4,8$. If $P$ is odd, then $P=2^k-1$ and $2P+1=2^{k+1}-1$. An even exponent $2s\ge4$ gives coprime factors $2^s-1,2^s+1$, both exceeding 1. This leaves $k=2$. If $P=2^k$, one of $k,k+1$ is odd. For an odd $l\ge5$ we have $3^{v_3(2^l+1)}=3^{1+v_3(l)}\le3l<2^l+1$ by lifting the exponent, so $2^l+1$ cannot be a prime power. This leaves $k\le3$ and proves the asymptotic assertion.

## Binary weights with an unrestricted multiplier

Take $P=2^e$, $d=2^j$, and $Q=dP$. Only the distinct-base case needs consideration. The predecessor has full powers $Q,P+1,mP+h$ with three different exceptional positions on one binary-carry-free row. Its part exceptional at $P+1$ is $Q L(P)$, where $L$ has nonnegative integer coefficients taken digitwise from the polynomial $(X+1)(mX+h)$. Indeed, the bottom $\log_2 Q$ digits of this part are zero, and the remaining binary digits split those of $(P+1)(mP+h)$. Group them into blocks of size $e$, using $P>m+h$. In particular $\deg L\le2$ and its coefficient sum is at most $2(m+h)$, independently of $j$.

Here $aP+b\mid L(P)$ because the part is divisible by $mP+h$ and $Q$ is coprime to it. Hence $aP+b\mid a^2L(-b/a)$. Its absolute value is at most $2(m+h)^3<P$, so $aX+b\mid L$ over the integers and $b-a\mid L(-1)$. But $Q L(-1)\equiv1\pmod{P+1}$. The powers of 2 modulo $P+1$ all have representatives $\pm2^i$ for $0\le i<e$, since $2^e\equiv-1$. These representatives have absolute value at most $P/2$. Since $|L(-1)|\le2(m+h)<P/2$, it must itself be a signed power of 2. This excludes every odd prime divisor of $|b-a|$ at all multipliers simultaneously.

## Checking every remaining input in the box

For each of the 8,192 families, list every $P\le B_0^5$ for which $P$ and $P+1$ are prime powers. One is even, so it suffices to test $2^e$ and $2^e-1$. After the third prime-power test there are 19,815 specializations, including repetitions across families. Removing shared-base cases and repeated weights leaves 8,438 distinct weights. The [consecutive verification through one billion](full_range_1b.md) handles the initial interval. Beyond it the two programs independently enumerate 27,480 eligible assignments at 3,134 weights. In each assignment they examine the row for which the assigned predecessor product $S$ is at least the normalized sum $V$. Such a row always exists because the two corresponding ratios have product $(N-1)/N$ and both cannot be strictly below one for integral complementary factors (if $S<V$ and $T<U$, then $ST\le(V-1)(U-1)<N-1$). All three exceptional bins must be occupied, and the normalized parts are their unique positive CRT residues below $S$. No row passes the required carries. The largest weight tested is $5{,}070{,}602{,}398{,}182{,}799{,}483{,}249{,}597{,}874{,}177$.

The first program uses the established row search. The second generates subset sums of CRT idempotents and uses factorial floors, checking its prime factors with SymPy. Run both with python-flint and SymPy installed.

```sh
python3 check_consecutive_power_partners.py
python3 check_power_partner_assignments.py
```
