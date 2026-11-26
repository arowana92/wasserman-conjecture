# Two high powers and a freely placed outer power

Let $p$ be a prime, let $n>h\ge1$, $a\ge1$, and $c\ge1$ be integers, and put $H=p^n$, $J=p^h$. Wasserman's full trinomial conjecture holds at
$$N=1+p^a(H+J+c)$$
**whenever $H>(2c)^{64}$ and $H/J>c$.** The three exponents are free within these conditions. The factors of $H+J+c$ can be arbitrary. Thus for each fixed $p,c$, once $n$ exceeds an explicit bound, only the last $\lfloor\log_p c\rfloor$ choices of $h<n$ could allow an exception, for any $a$.

**At every size**, for $1\le c\le8$ the conjecture holds at all these weights except possibly for seven choices of $(c,p,n-h)$:
$$ (5,3,1),\ (7,2,1),\ (7,2,2),\ (7,3,1),\ (7,5,1),\ (8,3,1),\ (8,5,1). $$
In particular $c=1,2,3,4,6$ work at every prime, every $n>h\ge1$, and every outer exponent $a\ge1$.

The [close-gap calculation](close_power_gaps.md) subsequently resolves offset $8$ at every size and restricts the cases at offsets $5,7$ to six explicit exponent families.

## The row forced by the prime $p$

Write $D=H+J+c$. Take the [endpoint](all_weights_endpoint_reduction.md) in a hypothetical counterexample and assign $p$ to a row $A$ avoiding it. By Kummer its parts are $p^a x,p^a y,p^a z+1$ with $x+y+z=D$ carry-free at $p$. Write $UV=N$ for the assigned full powers from the weight, with $U$ assigned to $A$. The other row has three positive parts that are multiples of $V$, so $U\ge3$. Also $U\mid x,y$, $\gcd(U,pD)=1$ and $z\ge1$.

Write $D=ST$ with $p^aS$ assigned to $A$. Every full power in $S$ divides $x$ or $y$. For a base different from $p$, at least one unshifted part is zero modulo that full power. Any extra $p$-power in $D$ divides both $x,y$ by carry freedom. Each full power in $T$ divides two normalized parts of the other row, whose sum is $U$. Consequently
$$27D^2\le U^3 S^2,\qquad S\mid\gcd(D,xy). \tag{1}$$
This is the [endpoint factor assignment](endpoint_reduction.md).

If $J\le c$, the [single high power bound](successor_powers.md) gives $27D^2<(J+c)^5\le(2c)^5$, already impossible at the stated threshold. Suppose $J>c$. Then the digits at $H,J$ each consist of one unit, and $x,y,z$ distribute these two units and split $c$ among them. If an unshifted entry is a pure power of $p$, its divisibility by $U$ is already impossible.

## A common short relation for three gcds

For $1\le A,B\le c$ write $G(A,B)=\gcd(H+A,J+B)$. There are relatively prime integers $0\le i\le j\le8$, $j\ge1$, satisfying $|in-jh|\le n/8$. To obtain them place the nine fractional parts of $th/n$ for $0\le t\le8$ in eight equal intervals, subtract two in the same interval, and divide the resulting pair by its gcd. The possible integer approximant lies between zero and $j$. If $i=0$ then $J\le H^{1/8}$ and every $G(A,B)\le2H^{1/8}$ since $J>c$.

For $i>0$, reduce $H^i=p^{in-jh}J^j$ modulo $G(A,B)$, moving a negative power to the other side first. Unless the corresponding integer difference vanishes,
$$G(A,B)\le 2c^8 H^{1/8}.$$
Vanishing requires $i,j$ both odd and $(H/A)^i=(J/B)^j$. If $i=j$ then $i=j=1$ and $H/J=A/B\le c$, excluded by hypothesis. Otherwise $j\ge3$. By coprimality there are relatively prime positive integers $s,t$ with $H/A=(s/t)^j$ and $J/B=(s/t)^i$ (compare each prime valuation). Thus $s^j\mid H$, $t^j\mid A$, and
$$G(A,B)\le AB\gcd(s^j+t^j,s^i+t^i)=AB(s+t)\le2c^2 H^{1/3}.$$
In fact $G(A,B)$ divides $AB$ times each of the two power sums, since $H+A=(A/t^j)(s^j+t^j)$ and $J+B=(B/t^i)(s^i+t^i)$. For the equality, $s/t$ is a unit modulo the gcd of the power sums and its $i$-th and $j$-th powers are both $-1$. Its square is one by Bézout and the odd exponents force $s/t=-1$ there. Conversely $s+t$ divides both sums. Any divisor of $G(A,B)$ coprime to $p$ is also coprime to $AB$ since it divides $H+A,J+B$. Such a divisor therefore divides $s+t\le2H^{1/3}$ in the vanishing case.

## Applying the bounds together

If the high units occupy different unshifted entries, write $x=H+r$, $y=J+s$, $z=t$, with $r,s,t>0$ and $r+s+t=c$. The bounds needed in (1) are
$$U\le G(r,s),\qquad S\le G(r,s+t)\,G(r+t,s).$$
Use the *same* $i,j$ for these three gcds. At most one difference vanishes, since $r^i/(s+t)^j<r^i/s^j<(r+t)^i/s^j$ whereas $H^i/J^j$ is fixed. If none vanishes, $U^3 S^2\le2^7 c^{56}H^{7/8}$. If the first vanishes, the coprimality of $U$ gives $U^3 S^2\le 2^7 c^{32}H^{3/2}$. If either of the other two vanishes, $U^3 S^2\le2^7 c^{44}H^{31/24}$. Each is less than $27H^2$ for $H>(2c)^{64}$.

In every other allocation some unshifted entry is a positive constant at most $c$. Hence $U\le c$. If the other unshifted entry is constant or takes both high units, its gcd with $D$ is at most $c$, by subtraction in the latter case. If it takes a single high unit, it must also take a positive constant, and its gcd with $D$ has the form $G(A,B)$ with $1\le A,B\le c$. Thus $S\le c\max(c,2c^8H^{1/8},2c^2H^{1/3})$. This again contradicts (1) at the stated threshold. This proves the general theorem.

## An additional bound for finite checking

When $p>c$, write $R=p^{\gcd(n,h)}$, $H=R^j$, $J=R^i$ with $j>i\ge1$ coprime. A counterexample requires
$$R\le(2c)^5.$$
Use the exact relation $H^i=J^j$ in the preceding proof. Each nonvanishing gcd bound becomes $2c^j$. For vanishing we still have odd $j\ge3$ and now $s\le R$, $t\le c<R$, so $G\le2c^2R$, or $U\le2R$ for the first gcd. In the separate-high-unit case the three respective upper bounds for $U^3 S^2$ are $2^7c^{7j}$, $2^7c^{4j}R^3$, and $2^7c^{5j+4}R^2$. All are below $27R^{2j}$ for $R>(2c)^5$. In the last two cases use $j\ge3$. In the constant-entry case the bounds are $c^7$, $4c^{5+2j}$, or $4c^9R^2$, with $j\ge3$ for the last case. These prove the claim.

For $p>c$ both $J>c$ and $H/J>c$ hold. The remaining finite list therefore has $R\le(2c)^5$ and $R^j\le(2c)^{64}$. The first verifier lists every such prime power $R$ of a prime $p>c$ and every coprime $1\le i<j$. At each allocation let $g$ be the largest divisor of $\gcd(x,y)$ coprime to $pD$. Test $27D^2\le g^3\gcd(D,xy)^2$. For $3\le c\le8$ the only surviving allocation is
$$c=8,\quad p=R=13,\quad (j,i)=(2,1),\quad D=190,\quad (x,y,z)=(171,18,1),\quad g=9.$$
Here $U$ would be divisible by 3, but $1+190\cdot13^a\equiv2\pmod3$ at every exponent. This rules it out.

For $p\le c\le8$ and $H/J>c$, the first theorem still bounds the remaining $H$ by $(2c)^{64}$. When $J>c$, restrict the constant split to the base-$p$ digit allocations of $c$. When $J\le c$, the single-high-power bound leaves only $27D^2<(J+c)^5$ and one can test every positive triple of normalized parts directly. No assignment passes (1) with a permissible scale.

The first program enumerates 10,068,947 inner power pairs with $p>c$ across the six offsets. The second program reverses this search. For each symbolic allocation it uses the determinants $(-A)^i-(-B)^j$ or the bounds in the vanishing case to cap $R$ separately, then admits even composite integer roots in those caps. This gives 4,889,124 individual root and row checks. Both programs find the single large-base allocation above. For $p\le c$ they independently test the same 218,611 inner power pairs, one using digit addition and the other using factorial floors.

At every prime, the digit sum satisfies $s_p(N)\le 3+s_p(c)$. The [five-unit theorem](six_digit_gaps.md) suffices when $s_p(c)\le2$. For $c=1,2$ this holds at every prime. For $c=3,4,6$ it holds at each of the remaining primes $p\le c$. For $c=5,7,8$, intersecting $s_p(c)>2$ with $p^{n-h}\le c$ gives exactly the seven choices displayed above.

```sh
python3 check_two_high_powers.py
python3 check_two_high_powers_cross.py
```
