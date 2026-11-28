# Enumerating by the least common multiple of the cuts

Take two binomial cuts $2\le a\le b<N$ with $\gcd(N(N-1),\binom Na,\binom Nb)=1$. Write $g=\gcd(a,b)$ and $a=gx,b=gy$. The positive integer $m=\operatorname{lcm}(a,b)/N$ satisfies
$$N<m^7,\qquad N<gm^5,\qquad ab(a+b)>gN(N-1). \tag{1}$$
Thus $ab>N^{6/5}g^{4/5}$. At a trinomial counterexample at least one of the two nonendpoint rows has all three parts greater than $(N(N-1)/2)^{1/3}$, by applying the last inequality to the smallest part of each row. The estimates also hold for sums of two parts used as cuts.

## A symmetric divisor equation

Use $h,k,L$ from the [ratio theorem](cross_cut_ratios.md) with $c=J=1$. In particular $k=gh/m$ and $k=m(N-a)(N-b)/(N(N-1))$. The second equality follows from $N=kxy/h$, $N-1=m(x-m)(y-m)/h$ and the two cut fractions. Write
$$u=m-k,\quad L=ux-m^2,\quad L'=uy-m^2.$$
Both $L,L'$ are positive by applying the ratio theorem in each orientation. The equation $Ly=m^2(x-m)-h$ gives
$$L L'=k m^3-u h=:F,\quad (ux,uy)=(m^2+L,m^2+L'). \tag{2}$$
Conversely the positive divisors of $F$ determine the possible cuts for fixed $m,k,h$, subject to the integral reconstruction and carry tests.

Since $1\le L\le L'$, we have $L\le\sqrt F<m^2$, $ux<2m^2$ and $L+L'\le1+F$. Hence
$$u^2xy\le(m^2+1)(m^2+F)<(m^2+1)(m^2+(m-1)m^3)<m^6.$$
For the last strict inequality, the difference from $m^6$ is $m^4(m-2)+m^2(m-1)>0$ since $m>k\ge1$. Substitution gives the stronger estimates
$$hu^2N<m^7,\qquad u^2N<gm^5. \tag{3}$$
Finally $m-k\ge1$ and $m=ab/(gN)$ yield
$$ab\left(a+b-1-\frac{ab}{N}\right)\ge gN(N-1),$$
which proves the remaining assertion in (1).

For each fixed $m$ there are $O_\epsilon(m^{1+\epsilon})$ candidate weights. There are fewer than $m$ choices of $k$, and then $h\mid km$ and $L\mid F<m^4$ give $O_\epsilon(m^\epsilon)$ choices together by the divisor bound. Thus the number of weights admitting any such pair with $m\le Y$ is $O_\epsilon(Y^{2+\epsilon})$. In particular, for almost all weights, every endpoint-coprime pair has quotient greater than $N^{1/2-\epsilon}$ for each fixed $\epsilon>0$.

## Quantitative form for every pair of cuts

For unrestricted binomial cuts retain $D,J,c$ from the ratio theorem, put $G=DJ$ and $\ell=\operatorname{lcm}(a,b)$. Then
$$G>\sqrt{\frac{N(N-1)}{\ell(a+b)}},\qquad
G>\left(\frac{N^8}{4\ell^7}\right)^{1/9},\qquad
G>\left(\frac{N^6}{4g\ell^5}\right)^{1/6}. \tag{4}$$
For example if two rows each have a part at most $N^{2/3-\epsilon}$, their coefficients and $N(N-1)$ share a divisor $\gg N^{3\epsilon/2}$ by the first bound. All estimates apply for $2\le a,b<N$.

Indeed $m=c\ell/N$, $k=cJm(N-a)(N-b)/(N(N-1))$ by the parameter identities, and $u=cJm-k$ is positive. Then $u\ge1$ implies
$$c^2J\ell\left(a+b-1-\frac{ab}N\right)\ge N(N-1).$$
For the other estimates define $L=ux-Jm^2$, $L'=uy-Jm^2$, both positive. The symmetric form of the elimination is
$$cLL'=Jkm^3-uh,\qquad LL'<(Jm^2)^2.$$
To verify it, multiply $cLy=Jm^2(cx-m)-h$ by $u$ and replace $ux$ by $L+Jm^2$. Thus
$$u^2xy\le(Jm^2+1)(Jm^2+LL')<4J^3m^6.$$
Using first $N=ckxy/h$ with $k<cJm$, and then $N=cgxy/m$, gives $N<4c^2J^4m^7$ and $N<4cgJ^3m^5$. Substitute $m=c\ell/N$ and use $c\mid D$ to obtain (4).

## Sharpness and an obstruction from a third part

The exponent seven and the leading constant one in $N<m^7$ are sharp for binomial tests using only the full powers dividing $N(N-1)$. In the [case $h=L=1$](cross_cut_critical.md), take $m$ to be a positive multiple of six and set $u=1$. Then
$$x=m^2+1,\qquad y=m^2(x-m)-1,\qquad N=(m-1)xy.$$
Here $\gcd(x,y)=1$ since $y\equiv m-1\bmod x$ and $x$ is odd. The criterion for all the lowest layers to pass is $\gcd(x+1,m+1)=1$, which here is $\gcd(3,m+1)=1$. These weights have $N/m^7\to1$.

The entire case $u=h=L=1$ is impossible in a trinomial counterexample when the larger cut $b$ is a single part. Indeed $y=(m-1)(m^3+m+1)$, so every prime of $m-1$ already divides $y$. The cut $a/N=m/y$ forces the row containing $b$ to avoid every prime of $y$. Thus its parts are divisible by the full product $(m-1)y$ from $N$. Likewise the primes of $x-m$ in $N-1$ are forced into that row, and $b\equiv1\bmod (x-m)$ by $N-b=(m-1)y(x-m)$. Its two companions would both be positive multiples of their sum.

Consequently for such a single part the bound sharpens to
$$N<\tfrac12m^7+3m^5.$$
Indeed from (2),
$$N=\frac{k}{hu^2}\left(m^4+m^2L+m^2F/L+F\right)
 <\frac{m^7}{hu^2L}+3m^5.$$
The excluded case is the only one with $hu^2L=1$.

## Reproduction by the quotient

For each $m$ try $1\le u<m$, put $k=m-u$, and take the positive divisors $h$ of $km$. For $F>0$ enumerate $L\mid F$ with $L\le\sqrt F$ and $L\equiv-m^2\bmod u$, then reconstruct $x,y$ from (2). Require $y$ integral, $m<x$, $\gcd(x,y)=1$, $\gcd(m,x)\mid h$, and $h\mid kxy$. Recover $N=kxy/h$ and $g=km/h$. Test parity, then the factorial valuations at the primes directly factored from $N$ and $N-1$. Allowing either side of each surviving binomial cut to split tests its extension to positive trinomial rows.

The comparison covers all $m\le1000$ and $x>4096$. Two enumerations gave 8,102,110 systems, with largest weight 998,002,997,000,999,999,001. Exactly 38,222 passed parity. Nineteen distinct binomial pairs survived the endpoint valuations, the largest at $N=144636344$ with cuts $1112112,1980486$. Their 20 assignments of endpoint primes all fail to extend to trinomial rows. The [earlier enumeration](cross_cut_ratios.md) handles $x\le4096$. Hence **every pair of cuts in the two nonendpoint rows of a counterexample satisfies $\operatorname{lcm}(a,b)>1000N$**.

The independent program loops over $m,u,x$, using $ux\le m^2+\lfloor\sqrt{km^3-u}\rfloor$. It tries $h\mid km$ in the residue class $h\equiv m^2(x-m)\bmod L$ and recovers $y$ via the linear equation. It uses 128-bit integers for the weight and prints every pair passing the binary test. The Python audit directly factors its weights and predecessors, checks primality and reconstructs the factorizations, and compares the counts, maximum and the complete set of surviving triples $(N,a,b)$ between the two methods.

~~~sh
mkdir -p .work
c++ -std=c++17 -O3 audit_cross_cut_quotients.cpp -o .work/audit_cross_cut_quotients
.work/audit_cross_cut_quotients 1000 4096 > .work/cross_cut_quotients_cpp.log
python3 check_cross_cut_quotients.py 1000 4096 --audit-log .work/cross_cut_quotients_cpp.log
lean CrossCutRatios.lean
~~~

The Python program uses python-flint. Lean checks the divisor identity, the extremal product estimate and the strict seventh-power bound from the integer parameters. [Lean source](CrossCutRatios.lean)
