# Closing the small gaps between two prime powers

Write $N=1+p^a(p^n+p^h+c)$, for a prime $p$, integers $n>h\ge1$, $a,c\ge1$, and let $d=n-h$. The full trinomial conjecture holds **at every such weight for $c=1,2,3,4,6,8$**. For $c=5,7$, a counterexample requires one of the six rows of this table and both additional conditions
$$p^h>2^{256},\qquad a+v=(2j+1)(h-b)\quad\text{for some integer }j\ge0.$$

| $c$ | $p$ | $d$ | $b$ | $v$ |
|---:|---:|---:|---:|---:|
|5|3|1|0|0|
|7|2|1|0|2|
|7|2|1|1|0|
|7|2|2|0|1|
|7|3|1|0|1|
|7|5|1|0|0|

The threshold applies to the inner power alone. There is no bound on the outer exponent in the exclusion. At the listed parameters a counterexample would also have to pass the companion-row conditions below and the [square-root restrictions on prime-power factors](scaled_row_blocks.md#two-square-root-restrictions-in-the-close-gap-families).

## The seven gap cases

The [two-high-power theorem](two_high_powers.md) leaves seven triples $(c,p,d)$ at $c\le8$. In each of them $p\nmid c$. Put $k=p^d$, $Q=p^h$, $D=(k+1)Q+c$. A row $A$ avoiding $p$ must have parts $(p^a x,p^a y,p^a z+1)$ with $x,y,z>0$ adding without carries to $D$. Take $U,S$ from the [endpoint assignment](two_high_powers.md#the-row-forced-by-the-prime-p), so $U\mid x,y,N$, $U\ge3$, and
$$27D^2\le U^3S^2,\qquad S\mid\gcd(D,xy). \tag{1}$$

For $Q\le c$ one can enumerate the positive carry-free splits directly. For $Q>c$ the high units $kQ,Q$ each go to a single slot, and the three constant terms split the base-$p$ digits of $c$. Thus each slot has the form $BQ+r$. In an unshifted slot a constant must be at least three and a pure power of $p$ is impossible. For two forms use $\det((B,r),(B',r'))=Br'-B'r$ to bound their gcd whenever the determinant is nonzero.

If $x$ or $y$ is proportional to $D$, say $(k+1)x=B_xD$, then $U\mid B_x$. Across the seven cases every such split has prime-to-$p$ part of $B_x$ less than three. If neither is proportional to $D$, then $S\le|\det(D,x)\det(D,y)|$. When $x,y$ have nonzero determinant, also $U\le|\det(x,y)|$. When both are constant, $U\le\gcd(x,y)$. These bounds cut off $Q$ by (1). At each remaining $Q$, take the largest divisor of $\gcd(x,y)$ coprime to $pD$ as an upper bound on $U$ in (1).

The remaining splits have proportional, increasing $x,y$. Their high units lie in different slots and, up to permutation, they are precisely
$$x=k(Q+s),\qquad y=Q+s,\qquad z=t=c-(k+1)s>0.$$
Since $p\nmid D$, every full power in $S$ divides $y$ and hence $t$. Here $s$ is always a power of $p$ including $1$. At offset eight, $(p,s,t)$ is $(3,1,4)$ or $(5,1,2)$. In both cases $U\mid Q+1,\,t p^a+1$. Modulo $U$, reduce $a$ modulo $h$ with a sign, choosing a remainder $r$ with $|r|\le h/2$. For $r\ge0$, $U$ divides one of $tp^r\pm1$, and for $r<0$ it divides one of $p^{-r}\pm t$. None vanishes since $t$ is not a power of $p$. Therefore $U\le(t+1)Q^{1/2}$, and (1) with $S\le t$ requires
$$[27(k+1)^2]^2 Q\le t^4(t+1)^6.$$
The same finite gcd test handles this bound. Across all seven triples the finite tests above use 691 allocations and exclude every one.

At the other five triples the pairs $(s,t)$ are precisely $(p^b,p^v)$ in the table. Put $H=h-b>0$, $X=p^H$, $G=X+1$. Now $S=1$ and $U$ divides both $G$ and $p^{a+v}+1$. If this gcd exceeds two it is $p^l+1$, $l=\gcd(H,a+v)$, with both quotient exponents odd. Indeed Bézout for the exponents gives $p^l\equiv\pm1$ modulo the gcd. The sign must be negative, both quotients must be odd, and the converse divisibility follows from the factorization of an odd power plus one. Unless $a+v$ is an odd multiple of $H$, one has $l\le H/3$. Then $U^3\le(1+p^{H/3})^3\le8p^H<27D^2$. This proves the asserted restriction at every exponent.

## A companion row independent of the outer exponent

In one of those six families $G\mid N$, since $G\mid p^a t+1$ and $D=(k+1)p^b G+t$. Every prime of $D$ is assigned to the companion row $C$ since $S=1$. Its parts are $(N/U)f_i$ with $\sum f_i=U$. Put $r_i=(G/U)f_i$ and $V=N/G$. Now $C=V(r_1,r_2,r_3)$, $\sum r_i=G<D$, and the residues of these three positive integers at every full power $w\parallel D$ form a permutation of $(G,0,0)$. This follows by cancelling the unit $V$ modulo $w$, where $GV=N\equiv1$. Thus $2w+1\le G$. The CRT determines the $r_i$ uniquely for each allocation of the powers to the three exceptional slots. Moreover
$$G/\gcd(r_1,r_2,r_3)\mid U,$$
since $G/U$ divides all three $r_i$.

There are four exponent-independent ways to discard one of these rows.

- The two unshifted parts of $A$ have ratio $k:1$. Their binomial coefficient is divisible by $k+1$, at every common scale. This excludes any prime of $k+1$ forced into $U$. In general $(u+v)/\gcd(u,v)$ divides $\binom{(u+v)T}{uT}$ by the identity obtained by multiplying the binomial by $u/(u+v)$.
- Apply the same identity to pairs of $r_i$. The companion must avoid every prime of $D$. For odd $p$, $k+1$ is even and $A$ carries in base two, so the companion must avoid two as well.
- At $w=\ell^e\parallel D$ the two $r_i$ divisible by $w$ have next $e$ base-$\ell$ digits in $C$ given by $(r_i/w)G^{-1}\bmod w$. A carry between them rules out every exponent.
- For $\ell\mid G/\gcd(r_i)$, let $w=\ell^{v_\ell(G)}$. The two unshifted parts of $A$ are divisible by $w$, and their next $v_\ell(G)$ digits are $kz,z\bmod w$ with $z=-p^b(G/w)t^{-1}\bmod w$. This uses $t p^a=X^{2j+1}\equiv-1\bmod w$. A carry contradicts assignment of $\ell$ to $A$.

For $p^h\le2^{256}$ the first verifier finds 2,417 unordered CRT rows. In the order above, the four tests exclude 1,397, 232, 779, and 9. The second program reconstructs the rows by partitions of the full powers, then evaluates the fixed digit layers from actual integer parts at a representative exponent. Both verify each factorization by primality and multiplication. Neither needs to factor a trinomial coefficient or the variable weight.

With python-flint installed, run

```sh
python3 check_close_power_gaps.py
python3 check_close_power_gaps_cross.py
```
