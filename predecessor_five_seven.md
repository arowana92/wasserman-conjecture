# The full factors five and seven in a three-factor predecessor

Let $N-1$ have exactly three distinct prime factors. In any counterexample the [endpoint row](all_weights_endpoint_reduction.md) occurs, and one of the other two rows $A$ avoids every prime of $N-1$. Write $UV=N$ with the primes in the unitary factor $U\ge3$ assigned to $A$, and the primes of $V$ to the last row $C$. Each part of $C$ is a positive multiple of $V$. For each full power $w$ of $N-1$ there is a part $k(N-1)/w$ of $A$, where $k$ is a positive multiple of $U$ strictly less than $w$. That part is 1 modulo $w$ and 0 modulo both other powers. [Factor assignment](endpoint_factors.md)

## Full factor five

**Every even weight $N=5p^a q^b+1$, with primes $p,q$ and $a,b\ge1$, satisfies the conjecture.** For odd weights of this form, any counterexample requires
$$N=5\cdot 2^k Q+1,\qquad k\ge3,\qquad Q=\ell^e,\quad \ell\equiv3\bmod8,\quad e\ \text{odd},\quad 3\parallel N.$$
Here $\ell$ is a prime exceeding 5. The binary expansion of $Q$ must have at least seven ones, and the ones of $3Q$ and $2Q$ must occupy disjoint positions. In particular the conjecture holds at $5\cdot2^k\ell^{2b}+1$ for all primes $\ell$ and all positive $k,b$.

If there are at most two distinct primes in $N-1$, the [two-factor theorem](endpoint_reduction.md) applies. Otherwise 5 is a full power in the three-factor predecessor. Its multiplier $k_5$ is 3 or 4 and equals $U$. The possibility $U=4$ is excluded by the [carry at a power one below the denominator](three_factors_small_power.md). Thus $U=3$ and $C=(V,V,V)$. The equal positive parts prevent $C$ from avoiding 2. If $N$ were even, the assignment would require $C$ to avoid 2 since $V=N/3$ would be even. Also $\gcd(U,V)=1$ implies $3\parallel N$.

For odd $N$ write $P=2^k$ for its predecessor's full power of two and $Q$ for the remaining odd full power; $k\ge3$ by the [small-power theorem](three_factors_small_power.md). Row $A$ has the three parts
$$3PQ,\quad zP,\quad (2Q-z)P+1,\qquad z=5l,\quad 1\le l<Q.$$
Here the second part is exceptional modulo $Q$, giving $5Pl\equiv1\bmod Q$. It is also divisible by 3. Only the first part is 1 modulo 5, so $3PQ\equiv1\bmod5$. Together with $3\mid N$ this gives $PQ\equiv7\bmod15$. Since $2\mid N-1$, $A$ must avoid 2. Hence $z$ and $2Q-z$ partition the ones of $2Q$, and $3Q$ has no one in common with $2Q$.

Because $Q$ is odd, addition of $Q+2Q$ starts with no carry into binary position 1 (positions start at 0). The 1 of $2Q$ in that position forces a 1 of $Q$ there, so that $3Q$ has a zero and the addition carries into position 2. In position 2, $2Q$ again has a 1, now with an incoming carry, forcing a zero of $Q$. Thus $Q\equiv3\bmod8$. Every odd square is 1 modulo 8, giving the restrictions on $\ell,e$.

Both $zP$ and $(2Q-z)P+1$ are positive multiples of 15. Such a multiple has at least four binary ones. A sum of two powers cannot vanish modulo both 3 and 5, since $2^i\equiv-1\bmod3$ requires odd $i$ and $2^i\equiv-1\bmod5$ requires $i\equiv2\bmod4$. For three powers to vanish modulo 3 their exponents have the same parity. Their residues modulo 5 are then all among $\{\pm1\}$ or all among $\{\pm2\}$, and their sum is nonzero. The two parts therefore use at least $4+(4-1)=7$ of the ones of $2Q$.

**Counting consequence.** The number of exceptional weights of the form $5p^a q^b+1$ below $X$ is $O(X^{3/5})$. For $Q$ alone, count odd binary strings with $3Q$ disjoint from $2Q$. Reading from low to high, keep the previous input bit and the carry in $Q+2Q$. In state order $(0,0),(1,0),(1,1),(0,1)$ the allowed successor lists are $(1,2),(3),(4),(1,3)$, indexing states from 1. This follows by forbidding output 1 whenever the previous bit is 1. For $\rho=303/200$, the positive vector $(1,\rho-1,\rho(\rho-1),\rho^2(\rho-1))$ has its successor sums bounded coordinatewise by $\rho$ times itself, since $\rho^4-\rho^3-\rho^2+\rho-1>0$. Thus there are $O(\rho^n)$ length-$n$ strings, even before the terminal carry check. Since $\rho^5<8$, the number of eligible $Q\le Y$ is $O(Y^{3/5})$. Summing over $k$ gives the claim by a convergent geometric series.

## Verification across two independently varying parameters

The full conjecture holds for $N=5\cdot2^k Q+1$ with odd prime power $Q<2^{42}$ and $1\le k\le2048$. The search allows composite odd integers $Q$ coprime to 30 as extra candidates for the necessary digit and congruence conditions. Cases with $Q$ supported on 3 or 5 are already excluded by the proof.

For each binary string passing the digit test, $2^kQ\equiv7\bmod15$ fixes at most one class of $k\bmod4$. The least positive inverse $l=(5\cdot2^k)^{-1}\bmod Q$ determines the entire row. On advancing $k$ by four, replace $l$ by the unique integer $(l+jQ)/16$ with $0\le j<16$ making the numerator divisible by 16. Retain $3\mid l$, all ones of $5l$ in $2Q$, and $9\nmid N$.

There are 10,869,329 odd strings through 42 bits satisfying the binary test, including leading zeros up to the specified length. Of these, 2,894,036 satisfy the cofactor and modulus-15 tests. The search returned 107 candidate triples $(Q,k,l)$. In every case the coefficient is divisible by **both 3 and 5**. The first overflowing prime-power layers are at most $3^8$ and $5^4$, respectively. These primes also divide the endpoint coefficient and the coefficient of $(N/3,N/3,N/3)$.

Use the following commands in the project folder. The Python verifier reconstructs the candidate rows and checks factorial valuations directly. It compares the C++ generator with a separate mask and congruence enumeration on a smaller box.

```sh
mkdir -p .work
c++ -std=c++17 -O3 check_factor_five.cpp -o .work/check_factor_five
python3 check_factor_five.py .work/check_factor_five
```

## Full factor seven at even weights

For every even $N=7p^a q^b+1$, the full conjecture holds unless $N\equiv6$ or $30\bmod36$. In either remaining residue class a counterexample with three distinct predecessor primes requires $U=6$ and a part $6(N-1)/7$.

We use a binary observation. If $x\ge y\ge0$ have disjoint ones, then $s_2(x-y)\ge s_2(x)$, where $s_2$ counts ones. Subtract the powers of two in $x-y$ individually from $x$. Each subtraction decreases just one binary digit, adding ones at lower positions when it borrows. All $s_2(x)$ originally occupied positions must become zero in $y$, so there must be at least that many subtractions. In particular for every odd $t>0$ and $u=2^h$, $h\ge1$, the positive integers $ut$ and $(u-1)t+1$ overlap in binary because their difference $t-1$ has $s_2(t)-1$ ones, fewer than $s_2(ut)=s_2(t)$.

If there are three distinct predecessor primes, the full factor 7 forces $3\le U\le k_7<7$. An odd $U$ here is 3 or 5, with $2\mid V$ forcing $C$ to avoid 2. The normalized parts $C_i/V$ of a binary-free row have distinct 2-adic valuations, since equal valuations give an overlapping bit after scaling. Their sum is therefore at least $1+2+4=7$, impossible. At $U=4$ we have $k_7=4$. With odd $t=(N-1)/7$, the part $4t$ and its complement $3t+1$ overlap by the binary observation, although $2\mid U$ requires $A$ to avoid 2. This leaves $k_7=U=6$. Since $U$ is a unitary divisor, $N/6$ is coprime to 6, which gives the two stated residues.

### A restriction on every base-seven digit

**Among even weights of the form $7p^a q^b+1$, at most $O((\log X)^6)$ weights up to $X$ could be exceptions. The full conjecture holds at every such weight through $7^{129}$.**

In the remaining case let $t=(N-1)/7$ and $z=t+1$. The distinguished part $6t$ is exceptional modulo 7, so $7\mid z$. The row must avoid 7, hence each digit of the sum $z$ of its other two parts is at most the corresponding digit of $N=7t+1$. To see the implication for $t$, write $t_i$ for its digit at position $i$, starting with the units. If $7^h\parallel z$, then $t_0=\cdots=t_{h-1}=6$ and $t_h\le5$. Above position $h$ the digits of $z$ are those of $t$. The corresponding digit of $N$ is $t_{i-1}$, so $t_i\le t_{i-1}$. This also holds through position $h$ by the description of the trailing sixes. Thus the digits of $t$ weakly decrease when read from low to high.

Every such $t<7^L$ is uniquely a sum of six repunits $(7^{h_i}-1)/6$ with $1\le h_1\le\cdots\le h_6\le L$. Indeed, the number of repunits extending beyond position $i$ equals $t_i$. This gives the count $\binom{L+5}{6}$ and the polylogarithmic bound.

Since $U=6$ is unitary, both 2 and 3 divide $z$ and $\gcd(N/6,6)=1$. In fact $z\equiv12\bmod72$. If $v_2(z)=1$, the identity $N=7z-6$ makes $N/2$ even. If $v_p(z)=j$ and $p^j\ge6$ for $p=2$ or 3, then modulo $p^{j+1}$ the least residue of $6t=6z-6$ is $p^{j+1}-6$ and that of $z$ is at least $p^j$, forcing a carry in an assigned base. Thus $v_2(z)=2$ and $v_3(z)=1$. Since $9\nmid N$, $z/3\equiv1\bmod3$, giving the stated class by CRT.

We also have $6z=\sum_{i=1}^6 7^{h_i}$. Modulo 8 the left side vanishes. The right side is 6 minus twice the number of odd exponents, so exactly three exponents are odd. For $L=128$ the checker combines the sums of three odd-length repunits with the sums of three even-length repunits in the prescribed class $t\equiv11\bmod72$. Out of 6,856,577,728 monotone possibilities, 232,664,124 pass that class. Exactly one, $t=294923$, makes $6t$ and $t+1$ disjoint in binary. Legendre's formula gives
$$v_3\binom{7t+1}{t+1}=4$$
there, contradicting avoidance of 3. This proves the finite assertion. Run `python3 check_factor_seven.py` to reproduce the enumeration and the valuation.

## Finiteness with prescribed prime bases

**For each fixed pair of primes $p,q$, all sufficiently large even weights $7p^a q^b+1$ satisfy the full conjecture, uniformly over the two positive exponents.** In fact the number of possible exceptional weights for a fixed pair has an absolute upper bound independent of $p,q$.

Here is the more general digit observation. Let $w=s^e$ be a prime power, $w\mid t+1$, and suppose the binomial coefficient $\binom{wt+1}{t+1}$ is coprime to $s$. In the base-$s$ representation of $t$, each sequence of digits whose positions have a fixed residue modulo $e$ is weakly decreasing from low to high. Indeed if $j=v_s(t+1)\ge e$, all positions below $j$ in $t$ are $s-1$. For positions $i>j$, the digit of $t+1$ is $t_i$ and that of $wt+1$ is $t_{i-e}$, so Kummer gives $t_i\le t_{i-e}$. At positions $e\le i\le j$ the earlier digit is $s-1$, giving the same inequality.

Each of the $e$ sequences starts at $s-1$, and is uniquely the sum of $s-1$ strings of consecutive ones starting at its lowest place. Consequently
$$(w-1)(t+1)=\sum_{r=0}^{e-1}\sum_{j=1}^{s-1}s^r w^{h_{rj}},\qquad h_{rj}\ge1.$$
For a fixed finite set of prime bases supporting $t$, there are only finitely many such $t$. Divide by $w-1$ and move $t$ to the left to obtain a linear equation with constant term 1. Each positive term exceeds 1. A vanishing nonempty subsum must include the negative term, leaving positive terms summing to 1, which is impossible. Thus the solution is nondegenerate. The coordinates are powers of $w$ and the single number $t$ supported on the fixed primes. They lie in a multiplicative group of rank at most $e(s-1)$ plus the number of fixed primes. The theorem of Evertse, Schlickewei and Schmidt bounds the number of nondegenerate solutions in terms of the number of terms and this rank. [Evertse, Schlickewei and Schmidt, *Annals of Mathematics* (2002)](https://annals.math.princeton.edu/2002/155-3/p04)

At $w=7$ with $t=p^a q^b$ the number of coordinates is seven and the rank is at most eight. The even-weight reduction forces the distinguished part $6t$ and avoidance of 7, so the binomial hypothesis holds. This proves the stated finiteness. It also applies to the multiplier-six case at odd weights.

For explicit bounds and all odd bases through 47, see the [two-exponent theorem](factor_seven_prime_bases.md). The [quadratic-residue theorem](factor_seven_residue_obstructions.md) gives eleven individual prime bases that work for every odd prime partner and both unrestricted exponents.
