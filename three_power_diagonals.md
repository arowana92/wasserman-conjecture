# Three high powers must have fixed differences

Let $C\ge2$, let $p$ be a prime, and put
$$N=p^a+p^b+p^c+C,\qquad a\ge b\ge c\ge0,\qquad p\mid N(N-1).$$
Write $n=\lceil\log_p(2C+2)\rceil$. Suppose $a>2048(n+1)$.

**Theorem.** If this weight is a counterexample, there must be $\delta\in\{0,1\}$ such that
$$s=1+p^{b-c}+p^{a-c}\ \mid\ C-\delta.$$
Moreover, for $l=(C-\delta)/s$, the three components of $l(p^{a-c},p^{b-c},1)+e$ add without a base-$p$ carry, where $e$ is zero for $\delta=0$ and is some coordinate unit for $\delta=1$. In binary, $a>b>c$ is necessary as well.

In particular, **at every offset** the conjecture holds throughout the indicated height range whenever $p^{a-c}>C$. Even at small spacings there are at most $\tau(C)+\tau(C-1)$ possible difference pairs $(a-c,b-c)$ at all these heights combined. Indeed $s-1$ determines the two powers uniquely, including the case of equal powers. For example, at binary offset 29 the only possible pair is $(2,1)$, since the only eligible divisor of 28 or 29 is 7. All other pairs work above the given height.

## Only one of the four divisors can be large

Use the [three-power allocation](three_power_translates.md#what-three-very-high-powers-force). Write $a_1=a,a_2=b,a_3=c$. In the binary or base-$p$ carry-free row the parts are $p^{a_i}+c_i$, where the $c_i$ add freely to $C$, and $a_3>4n$. These facts follow by the individual-part bound, independently of any spacing comparison. If $p=2$ any two equal high units merge before they are distributed, contrary to each of the three parts needing a unit of its own.

The products $U,S_1,S_2,S_3$ belonging to this row satisfy
$$27(N-1)^2\le U^3(S_1S_2S_3)^2. \tag{1}$$
Associate to each its target $e$, respectively zero or one of the three coordinate units. Remove its $p$-part and call the remainder $W_e$. The total power of $p$ removed is at most $p^n$, since $C,C-1$ are nonzero and strictly less than $p^{a_3}$. Each $W_e$ divides all three $p^{a_i}-x_i$ with $x_i=e_i-c_i$.

Call $W_e$ large here if it exceeds $B=2p^{32n+a/32}$. Choose *the same* two approximations for all four targets, with $1\le u_i\le32$, $0\le v_i\le u_i$ and $|v_i a-u_i a_i|\le a/32$ for $i=2,3$. They depend only on the exponents. If $W_e>B$, none of the $x_i$ is zero. The congruence between $x_1^{v_i}$ and $p^{v_i a-u_i a_i}x_i^{u_i}$ is an exact equality, since after clearing a negative power of $p$ the difference has absolute value at most $B$. In particular
$$ |x_1|^{v_i}=p^{v_i a-u_i a_i}|x_i|^{u_i}. \tag{2}$$
Here $v_i>0$, since otherwise $|x_i|=p^{a_i}>C$. These two equations uniquely determine all three positive magnitudes as soon as *one* is given. Every pair of distinct targets has a coordinate in common. That coordinate would force all three magnitudes to agree at both targets, impossible because $|j-1|\ne|j|$ for any integer $j$. Thus **at most one** $W_e$ exceeds $B$.

## The divisor of the single common power

Suppose $W_e>B$ and let $d$ be the lcm of the reduced denominators of $v_i/u_i$. Put $d_1=d$ and $d_i=dv_i/u_i$ for $i=2,3$. These positive integers have gcd one. As in the [valuation calculation](three_power_translates.md), there is a positive integer $t$ prime to $p$ with $|x_i|=p^{j_i}t^{d_i}$. In particular $t\le C$ and $0\le j_i\le n$.

Taking the valuation at $p$ in (2) gives the exact second relation
$$ a_i-j_i=(a-j_1)d_i/d.\tag{3}$$
The reduced-denominator definition implies $d\mid a-j_1$. Thus $g=(a-j_1)/d$ is a positive integer. Modulo $W_e$, both $p$ and $t$ are units. For the residue $y=p^g/t$ in the unit group, equation (3) and the original divisibilities imply
$$y^{d_i}=\operatorname{sgn}(x_i)\quad\bmod W_e\quad (i=1,2,3).$$
Take integer Bézout coefficients for the three $d_i$. Multiplying the three equalities to those signed powers shows $y=1$ or $-1$ modulo $W_e$, including when $W_e$ is composite. If $p^g=t$, coprimality to $p$ would fail. Hence
$$ W_e\le p^g+t\le p^{a/d}+C. \tag{4}$$

If $d\ge2$, the large divisor is at most $2p^{a/2}$. Its exponent in (1) is at most three. All other divisors are bounded by $B$, and the removed $p$-parts cost at most $p^{3n}$. Consequently, even allowing all nine small-factor exponents as an overestimate,
$$U^3(S_1S_2S_3)^2\le p^{3n}(2p^{a/2})^3(2p^{32n+a/32})^9 < p^{2a}.$$
Indeed the base-$p$ logarithm is at most $(3/2+9/32)a+291n+12$, less than $2a$ in the stated range. The same upper bound holds when there is no large divisor. Both contradict (1).

It follows that there is exactly one large divisor and $d=1$. All three $d_i$ now equal one. Since $W_e>2$, the three signs in the congruences for $y$ coincide. At most one $x_i$ could be positive because $x_i=e_i-c_i$ and $\sum e_i\le1$. Therefore all are negative. From (3), $j_i-j_3=a_i-a_3$ and $l=p^{j_3}t>0$. We obtain the claimed equality $c_i=lp^{a_i-a_3}+e_i$. The low digits already add freely, completing the proof.

## Independent finite enumeration

In binary, enumerate $k>h>0$, put $s=1+2^h+2^k$, and test each positive scale $l$ at which $sl$ or $sl+1$ is in range. The independent enumeration fixes the offset, factors $C$ and $C-1$, checks their divisors and tests the low parts directly. Both enumerations count all possible difference pairs. The same pair cannot arise in both the unshifted and shifted cases since $s>1$.

~~~sh
python3 check_three_power_diagonals.py
~~~

Among $2\le C\le1,048,576$, both enumerations exclude every high exponent triple at 1,027,663 offsets. At the other 20,912 offsets they find just 26,392 offset/difference pairs. The largest number of possible difference pairs at any one offset in this range is four. These counts concern the necessary condition, so the table does not assert an exception at any of the remaining pairs.
