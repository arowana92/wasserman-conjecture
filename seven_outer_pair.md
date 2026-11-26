# Congruences of the outer shift and prime exponent positions

Let $p$ be a prime, $p\mid N-1$, and $s_p(N)=7$. Denote the six positive exponents, with multiplicity, by $a+h_j$, where $\min h_j=0$ and $\max h_j=H$. In a hypothetical counterexample consider a row avoiding $p$ with the units-place 1 in a pair. This always covers odd $p$. Put $u=a+k$ for the other exponent in that pair, and $d$ for the gap of the other pair. Write $g=\gcd(d,u)$. Then
$$
  d/g\in\{1,3\},\qquad u/g\ \text{odd},\qquad
  g>H/4\ (p\text{ odd}),\qquad g>2H/9\ (p=2).
  \tag{1}
$$
Consequently, for a fixed allocation, a **single residue class** contains all possible exceptional outer shifts. There is at most one such shift in any interval of length $H/2$ at odd primes, or $4H/9$ at 2. This includes the other allocation at 2 as well. For an entire normalized pattern the corresponding bounds on the number of shifts are 60 at odd primes and 105 at 2.

**Prime exponents.** For every prime $p$ the conjecture holds at
$$ N=1+\sum_{j=1}^6 p^{q_j} $$
when the $q_j$ are primes, distinct if $p=2$. At odd $p$ each $q_j$ may independently be a prime or a power of two, with repetitions allowed. There is no span or exponent cutoff. We also prove the binary case in which the six exponents are distinct powers of two.

The arithmetic result behind these families is that $C=1$ in the allocation below forces $3\mid g$ at odd $p$, and forces $g$ to be neither prime nor a power of two at $p=2$. Moreover, at odd $p$, the quotient $d/g=3$ in (1) requires $9\mid d$.

## The common gap

Use the [seven-unit allocation](seven_span_and_gaps.md#which-allocations-need-checking). The parts are $p^a x,p^a y,1+p^u$ with $x=p^i(1+p^d)$ and $y=p^j(1+p^e+p^f)$; the inner singleton is $z=p^k$. The assigned weight factor $U$ divides all three parts. The predecessor powers in $A$ divide both $p^a x$ and $1+p^u$, and those in $C$ divide both $p^a y$ and $1+p^u$. The integers $U,A,C$ are odd, prime to $p$, pairwise coprime, and satisfy
$$ U\mid y,\quad A\mid y+z,\quad C\mid x+z,\quad
27p^{2H}\le U^3A^2C^2. \tag{2} $$
In particular $UA\mid \gcd(1+p^d,1+p^u)$. The odd divisor $U>1$ forces the two gap quotients $d/g,u/g$ to be odd: $p^{2g}=1\bmod U$ by Bézout and an even quotient would give $1=-1\bmod U$. The gcd of the two plus binomials is then $1+p^g$, either by the Euclidean algorithm or the same argument modulo the whole gcd. Also $\gcd(C,1+p^d)=1$ since $C\mid x+z$. In particular
$$ UA\mid1+p^g,\qquad u=g\ \Longrightarrow\ C=1. \tag{3} $$

By the [two-trinomial bound](two_trinomial_cuts.md) we have $C\le4p^{2H/3}$. At odd primes apply the [binomial–trinomial inequality](seven_pair_gaps.md) to the binomial at exponent $g$, giving $U<3p^{2g/3}$. Also $UA\le(1+p^g)/2\le2p^g/3$ by oddness. Thus
$$27p^{2H}\le U(UA)^2 C^2 < (64/3)p^{8g/3+4H/3}.$$
This proves $g>H/4$. In binary the [span-100 theorem](seven_span_and_gaps.md) allows $H>100$. Now $(1+2^g)^3\ge(UA)^3\ge(27/16)2^{2H/3}$, so $g\ge3$, and $(1+2^{-g})^3<27/16$ gives $g>2H/9$. At either base $d\le H$ excludes every odd quotient $d/g\ge5$.

Set $c=1/4$ at odd primes and $c=2/9$ at 2. Gaps $d\le cH$ are impossible. If $3\mid d$ and $d>3cH$, (1) requires $a+k\equiv d/3\bmod(2d/3)$, including the possibility $g=d$. Otherwise $a+k\equiv d\bmod 2d$. Each period is strictly greater than $2cH$. At 2, if 1 belongs to the triple, the [alignment theorem](seven_binary_alignment.md) gives one class modulo $2n$ for the outer shift, with the two unshifted pairs both having gap $n$; it also gives $27\cdot2^{2H}\le(1+2^n)^3$, hence $n>2H/3$. The shifted-pair allocations number at most $6\binom52=60$. The extra binary allocations number at most $\binom62\cdot3=45$, by choosing the pair of positive positions in the triple and matching the other four. This proves the spacing claims. Below we may take $H>100$ at all bases by the span theorem.

## What a missing companion factor requires

If $C=1$, equations (2) and (3) give $27p^{2H}\le (1+p^g)^3<27p^{3g}$, so $g>2H/3$ and $d=g$. In fact
$$ U>p^{H/6}. \tag{4} $$
For this use the circular gaps of the four exponents of $y+z$ modulo $g$. Rotating beyond a largest gap places their residues in an interval of length at most $3g/4$ starting at zero. Multiplication by a power of $p$ is invertible modulo $1+p^g$, and reducing exponents there only changes signs. The resulting signed sum is a nonzero integer: a zero would give $1+p^g\mid y+z$, contrary to $U\mid y$ and $U\nmid z$. Thus $A\le4p^{3g/4}$. Equation (2) implies $U^3\ge(27/16)p^{2H-3g/2}>p^{H/2}$, proving (4).

We use the same comparison for the index-three branch at odd primes. Set $h=\max(101,3g)$. If $d=3g$, then $H\ge h$ and the circular bound on $A$, together with $C\le4p^{2H/3}$, gives
$$ 256^6 U^{18}\ge 27^6 p^{4h-9g}. \tag{5} $$
By (1), $H<4g$ so $g\ge26$. For $C=1$ we have $g>2H/3$ and hence $g\ge68$; in this range $h=3g$. Any upper bound contradicting (5) also contradicts (4). We will exclude (5) at all odd bases where $3\nmid g$ that still have seven-unit cases, and exclude (4) in binary when $g$ is prime or a power of two.

## A lattice bound with the binomial exponent fixed

Put $J=1+p^e+p^f$. For an integer $n>0$ and integers $v,w$ define
$$ r=ne-gv,\quad s=nf-gw,\quad W=\max(0,r,s)-\min(0,r,s). $$
Suppose the two polynomials $1+T^n$ and $1+p^{r/n}T^v+p^{s/n}T^w$ have no common nonzero root. Then
$$ U\le 3^n p^W. \tag{6} $$
To see this, remove the least Laurent monomial and the least coefficient power from the second input. Its pairwise resultant with the monic binomial is an algebraic integer. At $T=p^{g/n}$ the values are divisible by $U$ in the integer ring with $p$ inverted. Thus so is the nonzero resultant. At every embedding its modulus is bounded by $(3p^{W/n})^n$ by taking the product at the $n$ unit roots of $1+T^n$. Since $U$ is prime to $p$, taking the norm gives (6).

In the following two cases this pair is coprime whenever $0<n<g$.

- At odd $p$ assume $3\nmid g$. Write $n',v',w'$ for the slopes divided by their common gcd. At a common root the [three-term phase classification](binomial_trinomial_gcd.md#locating-every-common-root) gives $n'\mid 2g$, since $e-(g/n')v'$ and $f-(g/n')w'$ are half-integers and the direction is primitive. For equal heights the primitive phase has order 6 with $n'$th power $-1$, so $3\mid n'$. The only unequal-height case has $p=3$ with a relative phase of order 12, again requiring $3\mid n'$. This contradicts $3\nmid g$.
- At $p=2$ assume $g\ge68$ is a prime or a power of two, the bit positions within $J$ are distinct, and $f<3g/2$. Again $n'\mid2g$. For prime $g$, $n'<g$ leaves only $1,2$. For a power of two, $n'<g$ is a power of two dividing $g$. Neither case allows the equal-height phase of order 6. Relative magnitude $\sqrt2$ requires $4\mid n'$ (a phase of order 8) and a noninteger half-height; the first condition fails for $n'=1,2$, and an integer scale $g/n'$ excludes the second when $g$ is a power of two. In the final binary case the two lowest terms have the same phase and the high magnitude is twice either. Primitivity forces the phase to be $-1$, hence $n'$ is odd and must be 1. The two lowest slopes differ by a nonzero even integer. Equal slopes at equal heights would repeat an actual bit. Their actual exponents then differ by at least $2g$, contrary to $f<3g/2$. These are all the [small-base phase cases](odd_seven_count.md#a-half-span-bound-on-the-first-common-divisor).

We give a uniform choice of direction. The triples $(n,r,s)=(n,ne-gv,nf-gw)$ form a lattice of determinant $g^2$. For $c>0$ the convex symmetric body $c|n|+W\le D$ has volume $2D^3/c$: the hexagon $W\le t$ has area $3t^2$. If $D^3>4cg^2$, the interior contains a nonzero lattice point. Indeed its half-sized copy has volume greater than a fundamental cell, so two of its points coincide modulo the lattice, by adding volumes in that cell. Their difference gives the assertion. With $D<g$ the point cannot have $n=0$, since then both residuals are multiples of $g$. Change its sign to take $n>0$.

For $p=2$ take $c=8/5>\log_2 3$, $D=g/6$. The volume criterion holds for $g\ge1383$. For odd $p$ take $c=1$, $D=g/6-1$. At $g=900$ we have $D^3>4g^2$ and $D^3/g^2=g(1-6/g)^3/216$ increases thereafter. In either case $n<g$. Thus (6) contradicts (4) in binary or gives $U\le p^{g/6-1}$ at odd primes. The latter contradicts (5), since $(27/256)^{1/3}>1/3$ and $h\ge3g$.

## The finite ranges

Only $p=2,3,5,7,11,13,17,19,23$ need checking by the [all-span theorem from base 29](seven_base27.md). The [modular-power verifier](check_seven_outer_pair.py) tests the following complete ranges.

- In binary take each prime or power of two $68\le g\le1382$, and $0<e<f<3g/2$. With $H_0=\max(101,g,f)$, (2) and the circular bound for $A$ require $U^3\ge(27/16)2^{2H_0-3g/2}$ when $C=1$.
- At odd bases take $26\le g\le899$, $3\nmid g$. Reduction modulo $2g$ and interchanging the last two terms of $J$ give $0\le e\le f<2g$. Use (5) with $h=\max(101,3g)$ for both branches as justified above.

Let $F$ be the odd part of $1+p^g$ and $G=\gcd(F,J)$. Since $p^f=-1-p^e\bmod G$, the integer $G$ divides $\gcd(F,(1+p^e)^g-(-1)^{g+f})$. Thus one modular power bounds all candidates of each parity of $f$. In binary 210,225 modular powers and 70,081 direct gcds give $G\le 2^{\lfloor(4H_0-3g)/6\rfloor}$ throughout. Each odd prime takes 539,566 modular powers. All but three triples meet $G\le p^{\lfloor(4h-9g)/18\rfloor-1}$. The three have $p=3,g=58$ and $(e,f)=(16,86),(30,46),(70,100)$, each with $G=12413$. Direct exponentiation gives $256^6\cdot12413^{18}<27^6\cdot3^{174}$. Thus every candidate contradicts its necessary bound.

The [independent verifier](check_seven_outer_pair_cross.py) uses ordinary integer gcds and (6). For each $n<g$ in the budget it rules out all residues $(e,f)\bmod g$ with $n\log_p3+W$ within the binary bound $g/6$ or the odd bound $(4h-9g)/18-1$, using the rational upper estimates above. All remaining pairs receive direct gcd tests. In the odd case one may translate the three residues modulo $2g$ past a largest gap, which multiplies the evaluated trinomial by a unit modulo $F$. This gives $2g-f\ge\max(e,f-e)$ and reduces duplication. The independent calculation checks 4,091,854 binary and 22,297,376 odd-base numerical gcds.

## Applying the result to the six exponents

Assume $H>100$, since every shorter span is [already certified](seven_span_and_gaps.md). If the absolute exponent $u$ is a power of two, the odd quotient $u/g$ forces $g=u$. For prime $u$, the bound $g>2H/9>1$ also forces $g=u$. Therefore $C=1$, $g=d=u>2H/3$. At odd bases the conclusion $3\mid g$ now excludes either type, since $u>3$. The binary conclusion excludes either type as well whenever the units-place digit belongs to a pair.

If all six exponents are prime, the other binary branch is also impossible. Its shifted triple has positive positions $r,r+h$ with $h=n+1$, $r\equiv0\bmod 2n$, or with $h=n-1$, $r\equiv n+1\bmod 2n$. Thus either $r$ or $r+h$ is divisible by $2n$, while $n>2H/3$ and neither prime can be an even integer greater than 2. If all six exponents are distinct powers of two, this branch would require two disjoint pairs at the same distance $n$. A difference $2^j-2^i$ with $j>i$ uniquely determines $i$ from its 2-adic valuation and then $j$. This completes the stated binary cases. At odd primes any carry in the proposed expansion reduces the digit sum to at most five and is [already covered](six_digit_gaps.md).

Another consequence needs no primality assumptions on the exponents. **At any odd base, if none of the six positive exponents is divisible by 3 and the largest is at least five times the smallest, the conjecture holds.** Indeed $3\nmid u$ rules out $d/g=3$ by (5), so $g=d>H/4$. The quotient $u/g$ is odd; it cannot equal 1 since that would give $C=1$ and $3\mid g$, and it cannot be 3 since $3\nmid u$. Thus $u\ge5g>5H/4$. Writing $a$ for the smallest exponent, the largest is $a+H$, so this forces $a>H/4$, contrary to $a+H\ge5a$.

For example the primes $2,12007,33013,55001,65537,65539$ give a weight at every prime base. The pair $2,65539$ has gap 65537, itself another exponent, so the exponent-divisibility test alone would retain that choice. The companion-factor argument rules it out.

Run the first verifier with python-flint. The second uses only the Python standard library.

~~~sh
python3 check_seven_outer_pair.py
python3 check_seven_outer_pair_cross.py
~~~
