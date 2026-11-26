# A central interval at every weight

Let $N\ge3$ and let $a,b$ be interior indices of its binomial row. If
$$\max\{(N-2a)^2,(N-2b)^2\}\le 5N,$$
then $\binom Na,\binom Nb$ share a prime dividing $N(N-1)$. This includes every choice of two cuts in this interval, with no restriction on the gap between them. Consequently the two nonendpoint rows of a trinomial counterexample cannot both have such a cut. Grouping the parts on either side of a cut shows that the corresponding binomial divides the trinomial coefficient.

For any fixed $c\ge1$, at most $O_c(\log(2X))$ weights $N\le X$ can fail the binomial assertion with $5N$ replaced by $cN$. More generally the count is $O(c^2\log(2c)\log(2X))$ with an absolute constant. The [earlier central-cut estimate](central_cuts.md) allows larger distances from the center when the gap is fixed. The [rational-proportion estimate](rational_cut_windows.md) treats other centers.

## Sum and difference of the offsets

Reflect and order the cuts so that $x=N-2a\ge y=N-2b\ge0$. Write $M=N-1$ and suppose the endpoint gcd is one. Then $x>y$, since $N/\gcd(N,a)>1$ divides both $N$ and $\binom Na$. Each full power in $N$ divides an index of a cut avoiding its base, so $N\mid ab$. Put
$$t=xy/N=N-2(a+b)+4ab/N,\qquad t\equiv N\bmod 2.$$
At each full power $w$ of $M$ an avoiding cut index is $0$ or $1\bmod w$. One of $x,y$ is thus $\pm1\bmod w$ and the other is $\pm t$ with the same sign since $xy\equiv t\bmod w$. If $t=1$, this gives $M\mid x-y$, impossible for $0<x-y<M$. In general it gives two exact identities
$$(x+y)^2=kM+(t+1)^2,\qquad (x-y)^2=(k-4t)M+(t-1)^2 \tag{1}$$
with $k$ a positive integer. In fact $x+y>t+1$ since $y<N$, $x<N-1$ and $x>y$, including the case $y=0$ with $N$ even. These local residue deductions use only the lowest prime-power layers of Kummer's criterion.

## The cutoff $5N$

Assume $x^2\le5N$, so $t<5$ by $x>y$. We check $N\le25$ directly by exact integer binomial coefficients in the verifier. From here $N>25$.

If $t=0$, then $y=0$ and $N$ is even. The exact central binomial coefficient is even. The cut with offset $x$ must be odd in coefficient, which forces $2^{v_2(N)}\mid a$ and $v_2(x)=v_2(N)$. Here $x^2=kM+1$, so $k$ is odd, $1\le k\le5$. For $k=1$ the equality $N=x^2$ contradicts the valuation. In fact **any** $k=1+2^v$ with $v\ge1$ is impossible for such a central pair. From $kN=x^2+k-1$ and $v_2(x)=v_2(N)$ we get $v=v_2(x)$. Choose $x_0=\pm 2^v$ with $x\equiv x_0\bmod 2^{v+2}$. At $x=x_0$ the weight is $N_0=2^v$, with a zero side. Otherwise put $h=v_2(x-x_0)\ge v+2$. The quadratic identity gives $N\equiv N_0\bmod 2^{h+1}$. The two actual sides $(N\pm x)/2$ are congruent modulo $2^h$ to $(N_0\pm x_0)/2+2^{h-1}$. Both representatives before the addition are in $[0,2^v]$, so both sides use bit $h-1$. This rules out $k=3,5$.

For $t\ge2$, put $\delta=k-4t$. By (1), $\delta\ge0$. Since $y=tN/x$, the increasing function $x-tN/x$ for $x\ge\sqrt{tN}$ gives
$$(x-y)^2\le (5-t)^2N/5.$$
For $t=4$ this forces $\delta=0$ and $x-y=3$, impossible since the offsets have the same parity. For $t=3$ it forces $\delta=0$ and $x-y=2$. These are adjacent cuts with $(N-2a-1)^2=3N+1$, the equality excluded at all weights by the [binary lifting argument](adjacent_cuts.md#why-equality-is-impossible).

For $t=2$, the bound gives $\delta\le1$. The offsets are even, so $\delta=0$ would require $x-y=1$. If $\delta=1$, (1) gives $(x-y)^2=N$ and $(x+y)^2=9N$. Thus $x=2y$ and $N=y^2$. At a square $N=y^2$, the cut with offset $my$ is even in coefficient whenever $m$ is a power of two and $y>m$. Indeed its two sides, doubled, are $y(y-m),y(y+m)$. Scaling both summands by two preserves binary carries. If $v_2(y)\ne v_2(m)$ their valuations agree, forcing a carry. Otherwise divide out the common power from $y,m$ and cancel the common square factor of two, leaving $z(z-1),z(z+1)$ for an odd $z>1$. Choose $\sigma=\pm1$ with $h=v_2(z-\sigma)\ge2$. Modulo $2^{h+1}$ the two products are $2^h$ and $2+2^h$, in some order. Both use bit $h$. Apply this with $m=1,2$ to finish the proof.

## Counting weights and recovering them exactly

Let $x^2\le cN$ and $H$ be the smallest power of two strictly greater than $N$. At least one cut must have odd coefficient. Its larger side must contain the leading bit $H/2$, whence
$$H-N\le x\le\sqrt{cH},\qquad H/2\le N<H.$$
The identity (1) has $0\le t\le c$ and $1\le k<4(c+1)$, since $x^2/M=x^2/N+x^2/(NM)<c+1$. For fixed $H,t,k$, the positive integer $v=x+y$ lies in an interval of length $O(1+\sqrt{ck})$. It determines $N$ by (1), and the discriminant $v^2-4tN$ must be the square of $x-y$. In particular at fixed $c$ there are $O_c(1)$ candidates per power of two.

For the uniform estimate omit the $O(c^2)$ weights $N\le (c+1)^2+1$. Then (1) gives $4t\le k$. Allowing the full range of $t$ adds $O((k+1)^{3/2}/\sqrt H)=O(1+\sqrt{ck})$ to the length of the interval for $v$. We also require $v^2\equiv(t+1)^2\bmod k$. Modulo $k$ there are at most $2k\tau(k)$ pairs of square roots of the same residue, counting ordered pairs and summing over residues. Indeed $(v,t+1)\mapsto(v-t-1,v+t+1)$ has kernel of size at most two and the number of pairs with product zero is $\sum_{j\bmod k}\gcd(j,k)\le k\tau(k)$. Consequently at each $k$ there are $O((k+\sqrt{ck})\tau(k))$ candidates in the interval. Summing over $k<4(c+1)$ using $\sum_{k\le z}\tau(k)=O(z\log(2z))$, and then over $H$, gives the claimed estimate.

## A fixed gap or an exact central cut at all sufficiently large weights

Fix $c\ge1$. Apart from effectively finitely many weights, every pair in this window with endpoint gcd one has $y=0$ or $x-y=t-1$ with $t=xy/N$. Thus when both cuts miss the exact center, their reflected gap at sufficiently large exceptional weights belongs to the finite list $(t-1)/2$ for odd $3\le t\le c$. In particular central pairs with nonzero offsets and a gap tending to infinity eventually have a common endpoint prime.

To see this, suppose $t>0$ and put $L=k-4t$. For $N>(c+1)^2+1$ we have $L\ge0$. The case $L=0$ gives the asserted gap, so take $L>0$. Set $v=x+y,w=x-y$. Equation (1) implies
$$L v^2-k w^2=C,\qquad C=4t\big(k-(t+1)^2\big).$$
If $C=0$ then $k=(t+1)^2$, $L=(t-1)^2$ and (1) gives $N=y^2,x=ty$ (recall $t\ge2$). At a square $N=y^2$ every weight prime divides the coefficient of the cut $b=y(y-1)/2$. For any prime $p\mid y$, $v_p(b)=v_p(y)$ for odd $p$ and $v_2(b)=v_2(y)-1$ for even $y$. Thus $p^{v_p(N)}\nmid b$, forcing $p\mid\binom Nb$. Every interior cut coefficient has a weight prime divisor, so this case gives a common one.

If $C\ne0$ and $D=kL$ is square, the nonzero fixed product $(Lv-\sqrt D\,w)(Lv+\sqrt D\,w)=LC$ allows only finitely many integer pairs. Otherwise take a unit $\varepsilon>1$ of norm one in $\mathbb Z[\sqrt D]$, obtained from Pell's equation. Each $Lv+\sqrt D\,w$ belongs to one of finitely many orbits $\gamma\varepsilon^j$. Indeed reduction by powers of the unit puts $\gamma$ in $[1,\varepsilon)$, and its conjugate is $LC/\gamma$, so both its integral coefficients are bounded effectively. For growing solutions $j\to+\infty$. Adding the two conjugates gives
$$ N=A\varepsilon^{2j}+O_{t,k,\gamma}(1),\qquad A=\gamma^2/(4kL^2)>0. $$
Binary carry freedom also gives $2^m=N+O_c(\sqrt N)$, with $m,j=O_{c,\gamma}(\log N)$. Hence
$$|m\log 2-2j\log\varepsilon-\log A|\ll_{c,\gamma} N^{-1/2}.$$
For a nonzero linear form in these fixed real algebraic logarithms, Matveev's theorem supplies an effective lower bound $(2+m+j)^{-C_1}$ with fixed effective $C_1$. Omit $\log A$ if $A=1$. These inequalities give an effective upper bound for $N$. If the form is zero exactly, taking norms in the real quadratic field gives $2^{2m}=\operatorname{Norm}(A)$, fixing $m$. This treats all finitely many $t,k$. [Matveev's theorem](https://www.mathnet.ru/eng/im314)

The dependence on $c$ gives a growing window. For every $\eta>0$ and all sufficiently large $N$, the same dichotomy holds whenever
$$x^2/N\le(\log N)^{1/2-\eta}.$$
Here the allowed gaps in the second alternative can grow with $t$. Equivalently all other exceptions must have maximum offset larger than $\sqrt N(\log N)^{1/4-\eta/2}$.

For details put $c=\max(2,x^2/N)$. Under the proposed bound, polynomial-size weights $N\ll c^{O(1)}$ can be dropped. We have $k,L,t=O(c)$, $D=O(c^2)$ and $LC=O(c^4)$. The square-$D$ and zero-logarithm cases above give $N\ll c^{O(1)}$. In the other cases choose a norm-one Pell unit with $\log\varepsilon\ll_\delta c^{1+\delta}$ for any $\delta>0$. Indeed the complete quotients of $\sqrt D$ in its period have the form $(u+\sqrt D)/v$ with $0<u<\sqrt D$, $0<v<2\sqrt D$, $v\mid D-u^2$. There are $O_\delta(D^{1/2+\delta})$ such states by the divisor bound. Each partial quotient is at most $2\sqrt D$. A convergent over one period gives a unit of norm $\pm1$, and squaring if necessary gives the asserted bound, after decreasing $\delta$. [Pell convergents and complete quotients](https://www.m-hikari.com/ija/ija-password-2007/ija-password5-8-2007/mollinIJA5-8-2007.pdf)

Write $\alpha=Lv+\sqrt D\,w=\gamma\varepsilon^j$ as above. Since $\gamma$ is integral, $1\le\gamma<\varepsilon$ and its conjugate is $LC/\gamma$, the absolute logarithmic heights of $\varepsilon$ and $A$ are $O(\log\varepsilon+\log(2c))$. For $A$ use the denominator $4kL^2$ and these same two conjugates. The identity $2Lv=\alpha+LC/\alpha$ now gives $N=A\varepsilon^{2j}+O(c^{O(1)})$ with an absolute polynomial error, since $\alpha\ge1$. Also $m,j\ll\log N$ since $\alpha\ll c^{O(1)}\sqrt N$ and $\varepsilon\ge2$. Thus the absolute value of the logarithmic form is $\ll\sqrt{c/N}+c^{O(1)}/N$. In degree two with at most three logarithms, the negative log of Matveev's lower bound is $O((1+\log\varepsilon+\log(2c))^2\log(2+\log N))$. The third height, that of 2, is constant. This yields $\log N\ll_\delta c^{2+\delta}\log(2+\log N)$ after adjusting $\delta$. Substitution of the proposed bound and a sufficiently small $\delta>0$ gives the contradiction. It suffices to treat $0<\eta<1/2$, which includes the larger values by containment.

## Reproduction

Run `python3 check_central_square_windows.py`. At $c=64$ the sum-of-offsets search finds 62 candidate pairs below $2^{4096}$. A second enumeration using the difference of offsets finds the same list below $2^{256}$, where all 62 occur. A direct index scan agrees through 12,000. Sixty candidates have explicit common endpoint primes confirmed by factorial valuations. The only binomial exceptions, up to independent reflection and ordering, are $(N,a,b)=(46,10,23)$ and $(32477,15886,15893)$. Exact coefficient gcds confirm both. Exhausting all positive trinomial refinements at these cuts finds a common endpoint prime in every pair. The numerical bound and the eventual fixed-gap theorem are separate conclusions.
