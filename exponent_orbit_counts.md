# Counting exponent pairs using the inverse modulo the fixed prime power

Let $P=p^a$, $Q=q^b$, $R=r^c$ for distinct primes and positive exponents. Write $D=P$. The results also hold with $D=2P$ when $p,q,r$ are all odd. We consider the full conjecture at $1+DQR$. Put
$$
 \beta_s=\log_s\!\left(\frac{1}{2\sin(\pi/(4s+2))}\right),
 \qquad \alpha=\beta_r,\quad \beta=\beta_q,\quad \tau=2-\alpha-\beta.
$$
Here $1/2<\alpha,\beta<1$. The implied constants in the two bounds below depend only on $q,r$.

**Independent progressions.** Take any arithmetic progressions of $B$ positive exponents $b$ and $H$ positive exponents $c$, with arbitrary positive steps. The number of pairs that can fail is at most
$$
 \ll BH\min\{1,(D/B)^{1-\alpha}\}\min\{1,(D/H)^{1-\beta}\}. \tag{1}
$$

**Exponents moving together.** Take $H$ positive exponent pairs $(b_0+hi,c_0+ki)$ on one arithmetic progression, with both steps positive. The number of indices that can fail is
$$
                \ll H^{\beta_*} D^{1-\beta_*},
      \qquad \beta_*=\min(\alpha,\beta). \tag{2}
$$
If just $c$ varies, (2) holds with $\beta_q$ for every fixed $b$. If just $b$ varies, use $\beta_r$.

For example, in **any two length-$H$ progressions** of binary and ternary exponents, all but $O(H^2/\log(2H))$ pairs work simultaneously at $1+D2^b3^c$ for **every prime power**
$$
                     D\le H^{\tau/(1+\tau)},\qquad
                     \tau/(1+\tau)=0.36258649\ldots .
$$
A base 2 or 3 in $D$ gives at most two predecessor primes, a case already proved. With any one fixed base $p\ge5$, the fraction of pairs that can fail at some $p^a\le H/f$ is $O_p(f^{-\tau})$ for $1\le f\le H$. On any single progression of positive indices $n$, the fraction that can fail at $1+p^a6^n$ for some $p^a\le H/f$ is $O_p(f^{-(1-\beta_2)})$.

For two intervals of **consecutive** exponents and a fixed base $p\ge5$, the conclusion extends beyond $H$. For each $\epsilon>0$ all but $O_p(H^{2-\epsilon})$ pairs work simultaneously for all $p^a\le H^{2-\beta_2-\epsilon}$ when $0<\epsilon\le1-\beta_2$. Here $2-\beta_2=1.30575808\ldots$. Starting points can be arbitrarily large in all these assertions.

## The forced numerator

At a failure with $D=P$ the [endpoint reduction](all_weights_endpoint_reduction.md) and [three-factor assignment](endpoint_factors.md#when-n-1-has-three-distinct-prime-factors) give a part $xQR$ of a row avoiding $p,q,r$, where
$$
                 0<x<D,\qquad xQR\equiv1\pmod P. \tag{3}
$$
For $D=2P$ the [four-factor assignment](two_prime_powers.md#a-four-factor-family) puts all three odd bases on one row with a positive part $xQR$. Its residue modulo $P$ is zero or one. A unitary divisor $U\ge3$ of the weight divides $x$ and is coprime to $P$. The only positive multiple of $P$ below $2P$ is $P$ itself, so the zero residue is impossible. Thus (3) holds with at most two representatives for a residue of $x$.

Cut off this part and discard the low $b$ digits at $q$ or the low $c$ digits at $r$. Kummer's criterion requires both additions
$$
     xR+(D-x)R=DR\quad\text{at }q,\qquad
     xQ+(D-x)Q=DQ\quad\text{at }r
$$
to be carry-free. For any arithmetic progression of $L$ exponents the counts for a fixed $x$ are respectively at most $O_q(W_q(x)L^\beta)$ and $O_r(W_r(x)L^\alpha)$, where
$$
              W_s(x)=s^{(1-\beta_s)v_s(x(D-x))}.
$$
The estimates are uniform in the common difference, including differences divisible by large powers of the carry prime. [Carry criterion and digit-block proof](three_factor_exponent_rectangles.md#how-many-exponent-terms-can-pass-a-cut)

For any set $X$ of $m$ distinct integers $0<x<D$ we have
$$
 \sum_{x\in X}W_s(x)\ll_s D^{1-\beta_s}m^{\beta_s},\qquad
 \sum_{x\in X}W_q(x)W_r(x)\ll_{q,r}D^\tau m^{1-\tau}. \tag{4}
$$
Indeed $s\nmid D$, so at most $2D/s^j$ numerators have $s^j\mid x(D-x)$ for $j\ge1$. The $i$th largest weight among any subset is at most $(2D/i)^{1-\beta_s}$. The sum of products is largest when both lists are sorted in the same order, by swapping inverted pairs. Sum the resulting $i^{-\tau}$ bound, using $\tau<1$.

## Tracking the numerator on the two axes

On the single progression in (2), $x\bmod P$ advances by multiplication by $(q^h r^k)^{-1}$. Write $T$ for its order. Each residue occurs in one index class modulo $T$. If $H\ge T$, each class has at most $2H/T$ terms and the possible representatives comprise at most $2T$ distinct integers. Apply the $q$ carry bound on each class and the first estimate in (4). The factor $T^\beta$ from summing the weights cancels the factor $T^{-\beta}$ from the progression length. This gives $O_q(D^{1-\beta}H^\beta)$. For $T>H$ the trivial count suffices since $T\le P\le D$. Use the $r$ cut as well when $b$ varies, taking the stronger estimate if $D\le H$ and the trivial estimate otherwise. This proves (2) and its one-axis forms.

For the rectangle first take $B,H\ge D$. Let $T,S$ be the orders of $q^h,r^k$ modulo $P$ for the two steps. The subgroup they generate has size $J$. Each residue of $x$ in the resulting coset occurs at exactly $TS/J$ pairs of index classes modulo $T,S$. This is the common fiber size of the group homomorphism $(i,j)\mapsto q^{hi}r^{kj}\bmod P$. Within each such pair of classes the two digit tests concern separate indices. Their counts multiply. Sum the second estimate in (4) over at most $2J$ numerator representatives to obtain
$$
 \ll \frac{TS}{J}(B/T)^\alpha(H/S)^\beta D^\tau J^{1-\tau}
   = B^\alpha H^\beta D^\tau
             (T/J)^{1-\alpha}(S/J)^{1-\beta}.
$$
Both ratios are at most one. If $B<D\le H$, apply the one-axis estimate to each $b$. Exchange the axes when $H<D\le B$. When both axes are shorter than $D$, use $BH$. This proves (1).

To sum over all prime powers $D\le Y\le H$ in a square of side $H$, their number is $O(Y/\log(2Y))$. The primes in $(n,2n]$ divide $\binom{2n}{n}\le4^n$, giving the prime count bound by dyadic summation. The higher powers contribute $O(\sqrt Y\log(2Y))$. Formula (1) costs $O(H^{2-\tau}Y^{1+\tau}/\log(2Y))$, proving the assertion at the displayed cutoff. For one fixed prime base, summing $D^\tau$ over its powers costs only $O_p(Y^\tau)$. The same geometric summation applies to (2).

## Longer cofactor range at fixed bases

For consecutive $b,c$ the orders of $q,r$ modulo $p^a$ are both $\gg_{p,q,r}P$. For a fixed unit base $u>1$, choose $l=p-1$ for odd $p$ and $l=2$ for binary $p$. With $v=v_p(u^l-1)$, binomial expansion gives $v_p(u^{lp^j}-1)=v+j$ for all $j\ge0$. Thus $u^l$ has order $p^{a-v}$ for $a>v$.

For fixed $x$ first use the unconstrained $r$ carry bound for $b$. For each surviving $b$, (3) confines $c$ to at most one residue modulo the order of $r$. Applying the $q$ carry bound on that progression and summing (4) with all $x$ gives
$$
           \ll_{p,q,r} D B^\alpha(1+H/P)^\beta .
$$
The analogous bound holds on reversing the two axes. For $B=H$, sum over $P\le H$ using the better orientation to obtain $O_{p,q,r}(H^{1+\beta_*})$. For $H<P\le Y$ the sum costs $O_{p,q,r}(YH^{\beta_*})$. This proves the stated extension.

## When one carry layer excludes every numerator

There is also a pointwise test working for **every** outer exponent $b$. At a layer $m=q^j$ of $xR+(D-x)R$, write $z=(R\bmod m)/m$. A carry at this layer for every $1\le x<D$ is equivalent to $\{Dz\}<\{xz\}$ for every such $x$, where braces denote fractional part. Its exact set in $0\le z<1$ consists of the disjoint intervals
$$
  \left[\frac{u}{D},\frac{u}{D}+\frac{1}{Dt}\right),
  \qquad 1\le u<D,\quad (u,D)=1,\quad
         1\le t<D,\quad ut\equiv-1\pmod D. \tag{5}
$$
The right endpoint is the next fraction after $u/D$ in the Farey sequence of order $D$. The total length is $D^{-1}\sum_{t<D,\,(t,D)=1}1/t$. These exclusions also apply to any cofactor for which the earlier assignment forces a part $xQR$, including $12,15,21,33$ when coprime to $qr$. [Assignment conditions](three_factor_exponent_rectangles.md#extension-to-four-predecessor-primes)

For a direct proof put $u=\lfloor Dz\rfloor$, $z=u/D+e$. If $g=(u,D)>1$, use $x=D/g$ to disprove the strict inequality. Otherwise write $s=xu\bmod D$. Then $D-x\le ts$ and $x\le t(D-s)$ by reduction modulo $D$. For $0\le e<1/(Dt)$ these give $De<s/D+xe<1$ for every $x<D$ with $x\ge1$. If $1/(Dt)\le e<1/D$, choose $x=t$. Its fractional part is $te-1/D<De$. This proves (5) including the endpoint conventions.

Take $l=q-1$ at odd $q$ and $l=2$ in binary. In any $lM$ consecutive terms of an arithmetic progression of positive exponents $c$, with $M=q^n$, at least
$$
           l\sum_{\substack{1\le t<D\\(t,D)=1}}
                       \left\lfloor\frac{M}{Dt}\right\rfloor \tag{6}
$$
terms work at **every** $b$. Split the index modulo $l$. In each class let $f=v_q(r^{lh}-1)$ for the original step $h$. At $m=q^{f+n}$ the $M$ residues of $R/m$ modulo one form a translate of the evenly spaced grid of mesh $1/M$, by the order calculation above. Each interval in (5) contains at least the floor of its length times $M$. For $D=5^{10}$, at least **2,195,614** terms in every progression of $2^{41}$ positive exponents $c$ thus certify all $1+5^{10}2^b r^c$, for every fixed prime $r\ne2,5$ and every $b\ge1$. The bound is uniform in the starting point, common difference, and prime $r$.

The same layer can handle a class of cofactors simultaneously. Fix $X\ge L\ge1$ and $(d,L)=1$. For each $t\mid L$ put $h/t$ equal to the reduced fraction with $1\le h\le t$ and $dh\equiv1\bmod t$ (use $1/1$ at $t=1$). Each interval $[h/t-1/(tX),h/t)$ certifies every eligible $2\le D\le X$ with $D\equiv d\bmod L$. For $t>1$ and $z=h/t-e$ in this interval, $\{Dz\}=1/t-De$. Every smaller $x$ has $\{xz\}>\{Dz\}$, since its residue numerator modulo $t$ is either positive or zero, in which case $\{xz\}=1-xe$. At $t=1$ the same comparison reads $1-De<1-xe$ if $De<1$, and also holds at $De=1$. Two distinct right endpoints with denominators $t,u$ differ by at least $1/(tu)$, so the intervals are disjoint. Their total length is $X^{-1}\sum_{t\mid L}1/t$.

The [orbit verifier](check_exponent_orbits.py) compares the residue-class enumeration to direct binomial cuts on short and long progressions, including collapsed orbits and doubled cofactors. The [billion-weight cube and its translations](periodic_exponent_cube.md) use exact carry witnesses at the bases $2,3,5$.
