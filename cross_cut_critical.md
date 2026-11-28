# Cuts close to the fifth-root bound

Use the notation of the [cross-cut theorem](cross_cut_ratios.md) with $D=J=1$ and $x\le y$. Then
$$hLN<108x^5/3125,\qquad m^2\equiv-L\pmod x,\qquad 0<m<x. \tag{1}$$
This gives a stronger estimate for almost all choices of $x$, as well as an explicit description when $hL=1$.

## The lowest prime-power layers

Put $X_0=x-m$, $Y_0=y-m$, $r=\gcd(m,x)$ and $s=\gcd(m,y)$. For a prime power $P=p^{v_p(N(N-1))}$ the lowest-layer condition asks that at least one cut have residue at most $N\bmod P$. This includes all smaller powers of that prime since $N\bmod P$ is 0 or 1. There may be carries at powers greater than $P$.

In the integral enumeration of the cross-cut theorem, the denominators of $a/N=m/y$ and $b/N=m/x$ in lowest terms are coprime, so every weight prime passes this condition. The two primitive numerator products available at the predecessor are $mY_0/s^2$ and $mX_0/r^2$. Thus the predecessor passes exactly when their least common multiple is divisible by $N-1=mX_0Y_0/h$. By multiplying both entries by $r^2s^2$ before taking their gcd, this is equivalent to
$$\gcd(r^2Y_0,s^2X_0)\mid h. \tag{2}$$

## When the maximum in (1) is approached

In particular $N\ge54x^5/3125$ forces $h=L=1$. In that case $r=s=1$ and
$$z=(m^2+1)/x,\quad W=m-z>0,\quad y=m^2(x-m)-1,$$
$$N=Wxy,\quad a=Wmx=y-m+1,\quad b=Wmy,\quad N-1=m(x-m)(y-m). \tag{3}$$
The term $x$ is odd since even $x$ forces odd $m$ and even $y$, contradicting $\gcd(x,y)=1$. Every prime divisor of $x$ is 1 modulo 4 by $m^2\equiv-1\bmod x$. Both $W,y$ are odd by (3). Condition (2) reduces to $\gcd(X_0,Y_0)=1$.

There is an additional obstruction if $a$ is a single part in a trinomial counterexample. For every $q^e=E\parallel x$, put $d=q^{v_q(W)}$. Then
$$E^2d<xW=m(x-m)-1. \tag{4}$$
In particular $x$ cannot be a prime power in this case. Indeed, every prime of $Y_0$ is absent from $mX_0$, the predecessor product of the other cut. Its full power in $Y_0$ must therefore be avoided by the row with part $a=Y_0+1$, so its other two parts are divisible by $Y_0$. At every prime of $x$ the other cut has reduced denominator $x$, so the entire power in $N$ divides every part of the row containing $a$. These divisibilities together make the two companions $xdY_0 j_1,xdY_0 j_2$ with positive $j_1+j_2=W/d$, using $N-a=WxY_0$.

Removing the full power $Ed$ of $q$ from the row preserves the condition of no base-$q$ carry. In the resulting row the two companions are $(x/E)Y_0 j_i$. Since $Y_0=x(z(x-m)-1)-1$, their residues modulo $E$ are those of $-(x/E)j_i$. If $(x/E)(W/d)\le E$, both $(x/E)j_i$ are strictly between zero and $E$, and the two residues sum to at least $E$. This forces a carry and proves (4). When $E=x$ the inequality fails because $W<x$.

If both $a$ and $b$ are single parts, $W$ cannot be a prime power. The argument above forces the first pair of companions to be multiples of $xY_0$ and the second to be multiples of $yX_0$ (using $b\equiv1\bmod X_0$). If a prime power $W$ is assigned to the first row its companions are multiples of $WxY_0$ and sum to $WxY_0$, which is impossible. The second row gives the same contradiction with sum $WyX_0$. The case $W=1$ is excluded by the same sums. The [four-factor restriction](opposite_part_lcm.md) also excludes a prime-power $m$ for every $h,L$ when both cuts are single parts.

At $x\le10^6$, exactly 121,202 systems in (3) pass (2) with $y\ge x$. Only two pass the binary carry test. At the following power of 3, each of their two cuts has residue greater than the residue of $N$, and $3\mid N(N-1)$.

| $x$ | $m$ | $N$ | Power of 3 |
| ---: | ---: | ---: | ---: |
| 22,405 | 17,648 | 124,380,638,800,347,396,945 | 9 |
| 24,565 | 2,928 | 11,751,878,708,987,668,945 | 27 |

## All products through 50

For $x\le10^6$ and $hL\le50$, there are exactly 5,383,328 parameter systems satisfying (2), of which 98 pass parity. The only pairs avoiding a common endpoint prime have $(N,a,b)=(46,10,23),(989,92,172)$. Neither extends to trinomial rows by the exhaustive split check in the cross-cut theorem. Thus at a counterexample, every pair with smaller reduced term $x\le10^6$ satisfies
$$159375N<108x^5. \tag{5}$$
Indeed $hL\ge51$ for such a pair, and (1) is strict.

The enumeration solves $m^2\equiv-L\bmod x$ for every modulus and $1\le L\le50$. At ordinary prime powers use both square roots modulo the prime and their unique lifts. At powers of primes dividing $2L$ test all lifts. Combine the roots by the Chinese remainder theorem. Then $k=m-(m^2+L)/x$, and each $h\le50/L$ gives at most one $y$. Check integrality, (2), binary carries and the factorial valuations at all prime factors of $N,N-1$.

An independent enumeration factors $m^2+L$, tries its divisors as $x$, reconstructs $N$ from $N-1=m(x-m)(y-m)/h$, and directly checks the least common multiple of the two primitive products. For $x\le4096$, $hL\le50$, both programs give 21,935 systems passing the lowest layers, 64 passing parity and the same two cut pairs. The divisor enumeration at $hL=1$, $x\le10^6$ reproduces both entries of the base-three table.

~~~sh
python3 check_cross_cut_products.py 1000000 50
python3 check_cross_cut_critical.py 4096 50
python3 check_cross_cut_products.py 4096 50
python3 check_cross_cut_critical.py 1000000 1
~~~

Both programs use python-flint.

## Sharpness for tests using the lowest layers

These tests alone permit $N/x^5$ arbitrarily close to $108/3125$. Fix $A=5^i$, $B=13^j$ with $A>B$ and $i,j\ge0$. Write
$$m=2t^2,\qquad x=\frac AB(2t^2+2t+1),\qquad z=\frac BA(2t^2-2t+1).$$
Choose $t$ so the first quadratic is divisible by $B$, the second by $A$, and $t\equiv0\bmod(2A^2+B^2)$. Each quadratic condition is equivalent to $(2t\pm1)^2\equiv-1$ at its odd prime power. The roots of $-1$ modulo 5 and 13 lift since their derivatives are units. The three moduli are pairwise coprime, so there are arbitrarily large such $t$.

Now $xz=m^2+1$ and, for large $t$, $z<m<x$. Both $x,z$ are odd. Define the remaining values by (3). Here $\gcd(x,y)=1$ because $y\equiv m-1\bmod x$ and a common prime would divide 2. Also
$$\gcd(X_0,Y_0)=\gcd(x-m,m+1)=\gcd(x+1,m+1).$$
Any prime in the last gcd divides $B+2At$ by $Bx=A(m+2t+1)$. Since $2t^2=m\equiv-1$ at that prime, squaring gives divisibility of $2A^2+B^2$. Our choice of $t$ then contradicts $m\equiv-1$. Thus all lowest layers pass. For example $A=5,B=1,t=102$ gives $x=105065$ and $N=63959196770257703633785$, and both cuts carry in base two.

For fixed $A,B$ let $t$ tend to infinity. With $\theta=B/A$, (3) gives
$$N/x^5\ \longrightarrow\ \theta^3(1-\theta)^2.$$
The irrationality of $\log_{13}5$ and density of irrational rotations allow $5^i/13^j$ to approach $5/3$. Taking $\theta\to3/5$ proves sharpness for the lowest-layer tests.

## Almost all reduced terms

For every $0<\epsilon<\log2$, outside a set of reduced terms $x$ of natural density zero, every endpoint-coprime pair satisfies
$$N\ll \frac{x^5}{(\log x)^{\log2-\epsilon}},\qquad
y\ll \frac{x^3}{(\log x)^{\log2-\epsilon}}. \tag{6}$$
The density here is on the smaller ratio term. There is a separate sparsity statement on weights. For fixed $C>0$, the number of weights $N\le T$ admitting such a pair with $x\le C N^{1/5}$ is $O_C(T^{1/5})$.

For a positive integer $L$ let $\chi$ be the nonprincipal quadratic character of the field $\mathbb Q(\sqrt{-L})$, of conductor at most $4L$. If $-L$ is a square modulo $x$, every prime $p\mid x$ with $p\nmid2L$ has $\chi(p)=1$. Write $x=wn$, where $w$ has all its prime factors in $2L$ and $\gcd(n,2L)=1$. If $x$ passes, $2^{\omega(n)}\le\sum_{d\mid n}\chi(d)$, where $\omega$ counts the distinct prime factors. This divisor sum is nonnegative for all $n$, by its prime-power factors. Summing first over $n=dk$ and using cancellation over a period gives, for $U\ge1$,
$$\sum_{n\le U}\sum_{d\mid n}\chi(d)
 \le\sum_{k\le U}\min(4L,U/k)\ll U(1+\log(2L)).$$
Furthermore
$$\sum_w\frac{2^{\omega(w)}}w=\prod_{p\mid2L}\left(1+\frac2{p-1}\right)
 \le\prod_{p\mid2L}(1-1/p)^{-2}\ll(1+\log(2L))^2.$$
For the last step use $\prod_{p\mid2L}(1+1/p)\le\sum_{d\le2L}1/d$ and $\prod_p(1-p^{-2})^{-1}=\zeta(2)$. The resulting bound is
$$\sum_{\substack{x\le X\\-L\ {\rm square\ modulo}\ x}}2^{\omega(x)}
 \ll X(1+\log(2L))^3.$$
Sum this over $L\le(\log X)^{\log2-\epsilon}$. Outside $o(X)$ integers $x\le X$, the Hardy–Ramanujan normal order of $\omega(x)$ gives $2^{\omega(x)}\ge(\log X)^{\log2-\epsilon/2}$. Division by this lower bound shows that only $o(X)$ integers admit such an $L$. The two bounds in (6) follow from (1) and the bound on $y$ in the cross-cut theorem. [Hardy and Ramanujan, *The normal number of prime factors of a number n*](https://ramanujan.sirinudi.org/Volumes/published/ram35.html)

More generally, for $C>0$ and $\delta\ge0$, the number of weights $N\in[T,2T]$ admitting an endpoint-coprime cut pair with $x\le C N^{1/5+\delta}$ is
$$O_{C,\delta}\!\left(T^{1/5+6\delta}\big(1+\log(2+T^{5\delta})\big)^4\right). \tag{7}$$
In particular for every fixed $\epsilon>0$, the weights admitting such pairs with $x\le N^{1/3-\epsilon}$ have density zero. The weights admitting a pair with $x\le N^{1/4}$ number $O(T^{1/2}\log^4 T)$ up to $T$.

To see this uniformly in $L$, let $\rho_L(x)$ count roots of $m^2\equiv-L\bmod x$. At odd $p\nmid2L$ there are two roots at every power for split primes and zero for inert primes. At $p\mid2L$, write $\nu=v_p(L)$. For $1\le j\le\nu$ the number of roots at $p^j$ is $p^{\lfloor j/2\rfloor}$. For $j>\nu$ it is zero if $\nu$ is odd. For $\nu=2k$ it is at most $2p^k$ at odd primes, or $4p^k$ at 2, by dividing the square by $p^{2k}$ and counting roots of a unit. Summing the geometric series gives
$$\sum_{j\ge0}\rho_L(p^j)/p^j\le 1+2/(p-1)\quad(p\text{ odd},\ p\mid L),$$
and an absolute bound at 2. Separating good and bad primes in the preceding character argument therefore gives $\sum_{x\le U}\rho_L(x)\ll U(1+\log(2L))^3$ with an absolute implied constant.

Here $x\le U= C(2T)^{1/5+\delta}$ and (1) bounds $hL$ by $H=1+(108/3125)U^5/T$. For fixed $x,m,h,L$ there is at most one weight. Summing the root estimate for $L\le H$, with at most $H/L$ choices of $h$, gives $O(UH(1+\log(2H))^4)$, as claimed. For $\delta=0$ dyadic summation recovers the stated $O_C(T^{1/5})$ bound up to $T$.
