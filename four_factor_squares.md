# Two arbitrarily large prime powers at cofactors 21 and 33

**Theorem.** For each $D\in\{21,33\}$, Wasserman's full three-row conjecture holds at all of the following weights.

- $N=1+D t^2$ for every positive integer $t$ having at most two distinct prime divisors.
- $N=1+D\,2^b r^{4c}$ for every prime $r$ and arbitrary positive $b,c$.
- $N=1+D\,2^b17^c$ for all positive $b,c$.

The second and third families have short proofs by congruences. The first uses a small enumeration of rows of sum less than 33. The weights themselves have no upper bound.

## Why the two bases must be assigned to the same row

Use the [all-weight endpoint theorem](all_weights_endpoint_reduction.md) and assign the primes of $N(N-1)$ to the two other rows that avoid them. Their weight scales are coprime integers $U,V\ge3$ with $UV=N$ and parts divisible by their respective scales. The endpoint coefficient is even, so at least one of these two rows has odd coefficient.

Write $N-1=D W$, with $W$ supported on two distinct primes. Let $Q,R$ be their **full** powers in $N-1$ and put $D_0=(N-1)/(QR)$, a divisor of $D$. If their bases were assigned to different rows, two parts of the first would be $UQ a,UQ b$ and two of the second $VR c,VR d$ for positive integers $a,b,c,d$. Positivity of the third parts gives $V\ge Q(a+b)+1$ and $U\ge R(c+d)+1$, hence
$$ (a+b)(c+d)<D_0. $$
Each full static power in $D_0$ divides one of $a,b,c,d$: two parts in its assigned row vanish modulo that power, and at least one is among the two displayed. This is the [four-factor test](two_prime_powers.md). In addition **at least one of $a,b,c,d$ is even**. If both displayed multipliers in the odd-coefficient row were odd, their parts would have the same lowest binary bit.

For a proper divisor $D_0$ of 21 or 33, the old test already excludes the split since $D_0$ is 1 or a prime. Otherwise $D_0=3s$ with $s=7$ or 11. Placing 3 and $s$ in different sums gives a product at least $4(s+1)>3s$. Placing both in one multiplier gives at least $2(3s+1)$. If they divide separate multipliers in the same sum, the least possible product with some multiplier even is $\min(2(6+s),3(3+s))=2(6+s)>3s$. This excludes every split.

Call the row taking both bases $A$, of scale $U$. Each power has two zero positions, so at least one part of $A$ is a positive multiple of $QR$. In terms of the original $W$, write that part $xW$. Then $3\le U\mid x$, $x<D$, and $\gcd(U,DW)=1$. If a static prime $p\mid D$ divides $W$, then $p\mid x$ by use of the **full** power in $QR$. Row $A$ must add its distinguished part and the rest
$$ xW,\quad (D-x)W+1 \tag{1}$$
without a carry at the two bases of $W$ and at every prime of $U$. The companion is $V(f_1,f_2,f_3)$ with $\sum f_i=U$. It cannot take a static prime $p$ if $2p+1>U$. It cannot take both 3 and $s$ since some positive $f_i$ would then be divisible by $3s=D>U$.

If $W$ has at most one prime divisor, $N-1$ has at most two prime divisors or exactly three including the full power 3. These cases are already [proved](endpoint_factors.md#when-n-1-has-three-distinct-prime-factors).

## The binary exponent is unrestricted

Set $W=2^b R$, $R=r^{4c}$, and take $r$ odd. The upper binary digits in (1) require the bits of $xR$ to fit inside $DR$. Since $R\equiv1\bmod16$, the bottom four bits of $x$ must fit inside those of $D$. Also $U$ is odd. At 21 the only possibilities $(x,U)$ are $(5,5),(20,5),(17,17)$. At 33 only $(17,17)$ remains. If $r$ divides $D$, its required divisibility of $x$ excludes every entry. If $r\mid U$, that entry is impossible because $U\mid N$.

At $D=21,U=5$, the relation $5\mid N$ gives $2^bR\equiv4\bmod5$. Since $R\equiv1\bmod5$, $b$ is even. Row $A$ must take 3, whereas $xW\equiv2\bmod3$ since $R\equiv1\bmod3$.

At $D=21,U=17$, we have $W\equiv4\bmod17$. The fourth powers of units modulo 17 form $\langle4\rangle$, the even powers of 2, and 2 has order 8. Thus $b$ is even. Now $17W\equiv2\bmod3$, excluding 3 from $A$. If $A$ took 7 it would need $W\equiv5\bmod7$, impossible because $W$ is a nonzero square there. The companion cannot take both 3 and 7.

At $D=33,U=17$, we have $W\equiv1\bmod17$, again forcing $b$ even. Since $2\cdot11+1>17$, row $A$ must take 11. This requires $W\equiv2\bmod11$, while $W$ is a nonzero square and 2 is a nonsquare modulo 11.

For $r=17$ allow instead $R=17^c$ at any $c\ge1$. The same binary test holds and $U=17$ is impossible. At the sole remaining case $D=21,U=5$, $W\equiv2^{b+c}\equiv4\bmod5$ forces $b+c$ even. Then $W\equiv1\bmod3$, giving the same contradiction at 3.

## Two arbitrary prime bases at square exponents

Now $W=t^2$. Try all $3\le U\mid x<D$ coprime to $D$ and all sorted positive triples $f$ of sum $U$. At each prime $p\mid UD$, independently try the squares $W$ modulo $p^k$, with $k=3$ increased if necessary to include one place beyond the full $p$-power in $U$. Keep these necessary conditions.

1. At $p\mid U$ the full $p$-power in $N$ equals that in $U$, and (1) adds without a carry through $k$ digits.
2. At $p\mid D$, if $p\mid W$ then $p\mid x$. Either (1) or the companion adds without carries through $k$ digits. Use $V\equiv (1+DW)U^{-1}\bmod p^k$ for the companion.

The tests allow different square inputs at each prime. Every actual counterexample would survive the relaxation. Two enumerations give exactly the following possibilities.

| $D$ | $x$ | $U$ | Normalized companion triples |
|---:|---:|---:|---|
| 21 | 15 | 5 | $(1,1,3),(1,2,2)$ |
| 21 | 19 | 19 | $(5,7,7)$ |
| 33 | 14 | 7 | $(1,3,3)$ |
| 33 | 23 | 23 | $(2,3,18),(2,6,15),(2,9,12),(6,6,11),(6,8,9)$ |

If $t$ is even, write $t=2^a v$ with $v$ odd. The upper binary addition in (1) requires $xv^2$ to fit inside $Dv^2$. Since $v^2\equiv1\bmod8$, this requires the lowest three bits of $x$ to fit in $D$. No displayed row passes.

If $t$ is odd, $N$ is even and every displayed $U$ is odd, so the companion must avoid 2. Its scale is $V\equiv(1+D)U^{-1}\bmod8$. At $D=21$ every listed triple has two equal parts with overlapping lowest bits at or below the third bit after scaling. At $D=33,U=7$ the repeated 3 gives the same obstruction. At $U=23$, $V\equiv6\bmod8$, and in every listed triple two scaled parts share a bit of value 2 or 4. This exhausts both parities.

Run `python3 check_four_factor_squares.py`. The first enumeration forms squares directly and compares digits. The second traverses companion triples first, characterizes square residues by valuation and unit part, and checks integer factorial floors. The verifier also checks 54,480 modular cases of the fourth-power proof and the separate unrestricted base-17 orbits.
