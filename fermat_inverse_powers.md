# Two inverse lifts at a four-prime predecessor

Let $P=65537$. For every positive integer $t$ with at most two distinct prime divisors, Wasserman’s full three-row conjecture holds at both
$$
                   1+P\,t^{8192},\qquad 1+2P\,t^{8192}.
$$
There is no restriction on the exponents in $t$. The second family includes predecessors with four distinct prime divisors. A separate result with completely independent exponents at fixed bases is
$$
        N=1+2\cdot65537\cdot2359333^{\,i}\cdot23855467^{\,j},
                   \qquad i,j\ge1 .
$$

## Two possible lifts and the binary test

Write $N=1+hPW$, $h=1$ or 2, $W=t^E$ for the indicated power $E=2^k$. In a counterexample the [one- and two-numerator reduction](static_square_cofactors.md#at-most-two-possible-numerators) gives $P\nmid t$, a row $A$ with part $xW$, and a second row $C=V(f_1,f_2,f_3)$, where
$$
   1\le x<hP,\quad xW\equiv1\bmod P,\quad UV=N,\quad
   3\le U\mid x,\quad f_i>0,\quad \sum f_i=U .
$$
Each row avoids all primes of its respective weight scale. The first row avoids the primes of $t$. The endpoint is even, so at least one of $A,C$ must avoid 2.

The inverse $x\bmod P$ is an $E$th power. For each odd prime $\ell\mid U$, the class $-hP\bmod\ell$ must also be an $E$th power, since $hP t^E\equiv-1\bmod\ell$. In particular $\ell\mid x$ and $\ell\nmid hP$. The following table lists every numerator with at least one eligible odd prime, with all such primes at that numerator.

| $h$ | $E$ | $(x,\{\ell\})$ |
|---:|---:|---|
| 1 | 16384 | $(65281,\{97,673\})$ |
| 1 | 8192 | $(65281,\{97,673\})$ |
| 2 | 8192 | $(65521,\{65521\}), (65793,\{7\})$ |

The search is small because the $E$th powers modulo $P$ form a group of size four or eight. At an odd divisor $\ell$ of a lift, check $(-hP)^{(\ell-1)/\gcd(E,\ell-1)}\equiv1\bmod\ell$. An independent enumeration raises every possible nonzero input modulo $P$ to its power, inverts the lifts and tests all input residues at each proposed weight prime. These lists suffice because $4\nmid N$. Indeed $N$ is odd for $h=2$ or even $t$, and $N\equiv2\bmod4$ otherwise. Hence $U\ge3$ contains an odd prime. Each listed $x$ is odd, so $U$ is odd as well.

Let $m=4E$. Every odd integer $d$ has $d^E\equiv1\bmod m$, by $d^2\equiv1\bmod8$ and repeated squaring. If $t$ is even, write $t=2^s d$ with $d$ odd. Since $A$ avoids 2, its binomial cut $xW,(hP-x)W+1$ adds without a binary carry. Above the bottom $sE$ positions the digits of $xd^E$ must fit inside those of $hPd^E$. At the modulus $m$ this requires every bit of $x\bmod m$ to be present in $hP\bmod m$. Every listed $x$ fails.

Now take odd $t$. At $h=2,E=8192$ the two residues $x\bmod32768$ are $32753$ and $257$, whereas $N\equiv3\bmod65536$. Thus $A$ carries at 2 and $C$ must avoid 2. In the two cases $U$ divides $65521$ or $7$, respectively, using the multiplicity of 7 in $65793$. Since $U\le65536$, every $f_i$ is positive and below 65536. The scale $V$ is odd, so all three residues $Vf_i\bmod65536$ are nonzero. They cannot add without carries to the two occupied low bits of $N$.

At $h=1,E=16384$ we have $U\mid65281=97\cdot673$ and $2\parallel N$. The companion avoids 2 since its scale $V$ is even and $U$ is odd. Divide its entries by 2. The resulting scale $V/2$ is odd and the total $N/2$ is 1 modulo 32768. Hence two of the $f_i$ must be positive multiples of 32768. Their sum with the third part would exceed $65536$, whereas $U\le65281$.

At $h=1,E=8192$ the same argument with $N/2\equiv1\bmod16384$ forces $U>32768$, leaving only $x=U=65281$. The two normalized companion entries divisible by $16384$ are $16384i,16384j$, with $i+j\le3$. They cannot be equal, because two equal positive integers have a common occupied bit at any scale. This leaves the companion direction
$$
        f=(16384,32768,16129),\qquad C/2=z f,\qquad
        z=V/2\equiv65281^{-1}\equiv257\pmod{16384}.
$$

**The coefficient in this direction is odd exactly when successive binary ones of its multiplier have positions differing by at least 16.** The forced prefix $257=2^8+1$ violates this condition and finishes the proof at $8192$.

Here is a 128-case certificate for the exact spacing. Read the multiplier digits starting at the units position. If the prefix read so far is $b$ in $j$ places, its state consists of $s_i=\lfloor f_i b/2^j\rfloor$. Reading digit $d$ gives the three values $s_i+df_i$. Their units must sum to at most one. Their quotients by two form the next state. At any carry-free prefix also $\sum s_i=\lfloor U b/2^j\rfloor$, by the Legendre layer test.

Consider the 16 states $S_k=(\lfloor f_i/2^k\rfloor)_i$ for $1\le k\le16$, with $S_{16}=0$. Any valid predecessor of one of these states again belongs to the list. To verify this for target $S$, choose $d=0,1$ and one of the four residue triples $u$ with sum at most one. The only possible source is $s_i=2S_i+u_i-df_i$. Retain it only if a common real $v\in[0,1)$ can satisfy $\lfloor f_i v\rfloor=s_i$ and $\lfloor Uv\rfloor=\sum s_i$. These are four half-open intervals with rational endpoints. The 128 tests give exactly the zero-digit edges $S_k\to S_{\min(k+1,16)}$ and the one-digit edge $S_{16}\to S_1$. All are actual legal digit transitions. Any accepted integer ends at zero after enough leading zeros, so backwards induction puts its entire path in this list. Conversely, reading ones at distance at least 16 follows these same edges.

A separate forward calculation starts at $(257,514,253)$ after the prescribed 14 low digits. It reaches 3,201 states and checks that none has three pairwise disjoint binary entries when the remaining input is zero. This independently excludes every completion of the prefix. Lean verifies the 128 predecessor tests and the backward induction excluding any finite completion of this state. [Lean certificate](FermatCompanion.lean)

## A direct spacing theorem

For integers $h\ge3$ and $1\le b\le h-2$, put $A=2^{h-2}$ and $c=2^b-1$. When $h\ge2b$, the coefficient on the ray
$$
                     z(A,2A,A-c)
$$
is odd **exactly** when successive ones of $z$ have positions differing by at least $h$. The preceding direction has $h=16,b=8$. The condition $h\ge2b$ is sharp for the claimed spacing.

For the proof let $S=4A-c$, assume the addition is carry-free, and write $z=M+y$, where $M$ is its highest binary power and $0\le y<M$. All parts are below $4AM$. As $4AM$ is a power of two, their carry-free sum is also below it. Hence $Sy<cM$. Since $c^2\le S$, we get $cy<M$.

Set $a=\lfloor Ay/M\rfloor$. The three quotients $\lfloor f_i y/M\rfloor$ for $f=(A,2A,A-c)$ are $a,2a+\epsilon,a-\delta$ with $\epsilon,\delta\in\{0,1\}$. They are all below $c$ by $Sy<cM$. The high parts $\lfloor f_i z/M\rfloor$ must have disjoint bits. Since $2^b\mid A$, their low $b$ bits, interpreted as integers, are exactly
$$
                         a,\quad 2a+\epsilon,\quad 1+a-\delta .
$$
If $a>0$ and $\delta=1$, the first and last coincide. If $\delta=0$, the highest bit of $a$ is shared by $a+1$ unless $a$ consists entirely of ones up to that bit. In the latter case $a+1$ shares its occupied bit with $2a+\epsilon$. Thus $a=0$, necessarily $\delta=0$, and disjointness forces $\epsilon=0$ too. All $f_i y$ are consequently below $M$. They occupy only the low block and add freely there, so $Sy<M$. If $y>0$, its highest position lies at least $h$ places below that of $M$, since $S>2^{h-1}$. The addition on the tail $y$ remains carry-free. Induction proves the spacing condition. Conversely, the three fixed slopes have disjoint bits in the same $h$ positions, so shifted copies at the asserted distances add freely.

When $h<2b$, the support of $S$ is $\{0\}\cup\{b,\ldots,h-1\}$. Its translate by $h-b$ is disjoint from it. Consequently $z=1+2^{h-b}$ makes the coefficient odd despite a distance smaller than $h$.

## The smaller exponent 4096

The full conjecture holds at $N=1+65537t^{4096}$ for every **odd** $t$ with at most two distinct prime divisors. It also holds for even $t$ in that class if its odd part is $1$ or $7$ modulo 8, or if $t$ is a cube. In particular every $1+65537t^{12288}$ works without a parity restriction.

The odd-prime table at this exponent has exactly three entries,
$$
          (49153,\{13\}),\quad (65281,\{97,673\}),\quad(65533,\{13\}).
$$
All $x$ are odd. For odd $t$ the normalized companion must avoid 2 after removing the factor 2 in $V$. Since $N/2\equiv1\bmod8192$, at least two $f_i$ are multiples of $8192$. Thus $U>16384$, forcing $x=U=65281$. The two divisible parts are $8192i,8192j$, with $1\le i<j$ and $i+j\le7$. This gives nine directions. In each, $z=V/2\equiv U^{-1}\bmod8192$ is necessary. Forward binary searches from these prefixes visit 10,000 states in total, and none can end in a carry-free state.

For an independent reverse check let $f$ be one direction and append its sum to obtain $g=(f_1,f_2,f_3,U)$. Enumerate the states $(\lfloor g_i v\rfloor)_i$ over $0\le v<1$ by evaluating at every left breakpoint $j/g_i$. Retain states whose first three coordinates sum to the fourth. All prefixes of a carry-free integer must satisfy this equality. A digit $d$ sends a state $s$ to $(\lfloor(s_i+d g_i)/2\rfloor)_i$. Working back from zero across the retained edges gives all possible states of completable prefixes. The nine enumerations inspect 1,093,118 cells in total and exclude every prescribed prefix.

For even $t=2^s d$ with $d$ odd, the binary test on $xd^{4096}$ inside $P d^{4096}$, first taken modulo 16384, leaves only $x=49153$ and $U=13$. Thus $t^{4096}\equiv3\bmod13$ and $t^{4096}\equiv4\bmod65537$. The former is impossible for a cube $t$, since every unit raised to $3\cdot4096$ is 1 modulo 13. If $d\equiv\pm1\bmod8$, then $d^2\equiv1\bmod16$ and repeated squaring gives $d^{4096}\equiv1\bmod32768$. At this modulus $x\equiv16385$ and $P\equiv1$, again forcing a carry. These are conditions on any unresolved input at 4096.

## Two independent exponents at four fixed primes

More generally take $N=1+2Pt$ with $t$ coprime to $P$ and having at most two distinct prime divisors. Suppose $t\equiv1$ or $-1\bmod P$ and $t$ is a nonzero quadratic residue at each of $3,11,331,43691$. Then both inverse lifts are coprime to $N$.

For $t\equiv1\bmod P$ the lifts are $1,P+1$. The odd prime factors of $P+1$ are $3,11,331$, at each of which 2 is a quadratic nonresidue. Hence $N\equiv1-2t$ cannot vanish at any of them. For $t\equiv-1\bmod P$ the lifts are $P-1=65536,2P-1=3\cdot43691$. The two odd prime factors of the latter are 3 modulo 4. At either, $N\equiv1+t$ cannot vanish. In both cases $N$ is odd. The decrement identity therefore makes the binomial for each possible cut divisible by the full weight.

The primes $2359333$ and $23855467$ are $1$ and $-1$ modulo $P$, respectively, and both are squares at the four auxiliary primes. Every product of their positive powers meets the hypotheses.

## Uniform counting with four predecessor primes

The [uniform cut bound](inverse_weight_cuts.md#uniform-estimates-from-one-or-two-cuts) extends to the two lifts. In the notation of that result take $N=1+2FDw$ with the $P_i$ odd and allow each inverse multiplier in $[1,2P_i)$, still congruent to $t_i^{-1}$ modulo $P_i$. Count an input if **any** such positive row lacks full weight divisibility. For $\ell\mid N$ the map from $w\bmod m$ to $N\bmod m$ at $m=\ell^{k+1}$ is still bijective, since $\ell\nmid2FD$. A fixed residue for $a_i\bmod m$ now has at most $1+2P_i/m$ actual lifts, and each lift still specifies $w\bmod P_i$. The tail uses at most $2B/\ell$ multipliers at the smallest modulus and $\ell<2B$. These constant factors preserve (1), including both exponents in that estimate.

For $r=1$ the two-numerator reduction applies when $Dw$ has at most two distinct prime factors. For $r=2$, take $D=1$ and $P_1,P_2,w$ powers of three distinct odd primes. All three odd bases must go to one row by the [four-factor assignment](two_prime_powers.md#a-four-factor-family). They have distinct exceptional positions, since a part exceptional for none would be divisible by $P_1P_2w$ and the coprime scale $U\ge3$, exceeding $N$. Thus both lifts occur among the tested possibilities.

In particular for **any odd prime powers $P,Q$ of distinct bases** and $\kappa<2\sqrt{\log 2}$, the conjecture holds simultaneously at both $1+PQR$ and $1+2PQR$ outside a fraction $O_\kappa(\exp(-\kappa\sqrt{\log P}))$ of primes $R\in[P^2,2P^2]$, uniformly for unrestricted $Q$. This follows by adding the two one-cut bounds. If $R$ repeats a base, the doubled predecessor has at most three distinct primes including the full power 2, hence is covered by the [small-power theorem](three_factors_small_power.md).

Run `python3 check_fermat_inverse_powers.py` to reproduce the inverse lists, both companion searches, the residue classes and the prime factorizations. Run `lean +leanprover/lean4:v4.36.0 FermatCompanion.lean` for the finite-state proof at exponent 8192.
