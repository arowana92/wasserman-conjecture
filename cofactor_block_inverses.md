# Modular inverses across gaps in the digits of a cofactor

Let $p$ be prime, $T>1$ coprime to $p$, and suppose $T$ has at most two distinct prime factors. At an exception to Wasserman's conjecture at $1+p^eT$, $e\ge1$, the least positive inverse $r_e$ of $p^e$ modulo $T$ must fit digit by digit inside $T$ in base $p$. In fact $r_e$ is a carry-free sum of two positive integers, each divisible by a different prime factor of $T$. To recall the argument, if $T$ has only one prime factor the conjecture already holds. Otherwise the unique predecessor row has parts $p^e x,p^e y,p^e z+1$ with $x+y+z=T$ carry-free. The two unshifted parts have residues $(1,0)$ and $(0,1)$ at the two full powers in $T$, and $z\ge1$. Thus $r_e=x+y<T$. In particular a pure power of $p$ cannot be $r_e$. [CRT and inverse criterion](inverse_residue_cofactors.md#the-inverse-and-its-digits)

**A zero gap certifies an individual weight.** If the base-$p$ digits of $T$ have at least $e$ consecutive zeros between occupied positions, the conjecture holds at $1+p^e T$. More generally the same conclusion holds when $e$ modulo the order of $p$ modulo $T$ agrees with a positive integer no larger than such a gap. The cofactor may have any number of other occupied positions.

Here is a version applying to any $p$-avoiding row, without a restriction on the prime factors of $T$. Write $W=p^e$, and let $B=1+Wx$ be its part equal to one modulo $W$. Then $x$ fits digitwise inside $T$. Put $s=T/\gcd(B,T)$, so $k=sB/T$ is a positive integer below $sW$ and coprime to $p$. At any occupied position $n>0$, put $L=T\bmod p^n$ and $H=\lfloor T/p^n\rfloor$. Then
$$s(1+WL)\ \ge\ p^n. \tag{1}$$
Indeed if the reverse strict inequality held, split $x=Jp^n+I$ with $0\le I\le L$. In $kT=s+sWx$, both low quantities $kL$ and $s+sWI$ would be below $p^n$. Hence $kH=sWJ$, impossible since $p\nmid kH$. The three-factor row above has $s=1$. A gap of $e$ zeros before $n$ gives $L\le p^{n-e}-1$, contradicting (1). For transferring the exclusion to every congruent exponent use $x=T-r_e$ for the three-factor row: the base-$p$ digits of $x$ fit precisely when those of $r_e$ do. At a smaller congruent exponent $j$ use the same $r_e$ and $B=1+p^j(T-r_e)$. This satisfies $T\mid B$ and $B<p^jT$. The proof of (1) then applies directly to $x=T-r_e$ if the digits fit.

## Separated exponent clusters

Fix any split $T=Hp^n+L$ with $H,L>0$, $L<p^n$. If $r_e=Jp^n+I$ fits, then $J\le H$, $I\le L$, and
$$\Delta_e=HI-LJ=H r_e-JT,\qquad
 0<|\Delta_e|\le HL,\quad p^e\Delta_e\equiv H\bmod T. \tag{2}$$
The determinant cannot vanish because $\gcd(r_e,T)=1$ and $0<H<T$. It is divisible by $g=\gcd(H,L)$.

Take any integer $b\ge1$ such that $(p^b+1)HL<T$. If two exceptional exponents differ by $d\le b$, then $\Delta_e=p^d\Delta_{e+d}$: the difference is a multiple of $T$ of absolute value less than $T$. Consequently in any chain of exceptional exponents with successive gaps at most $b$ the distance from the first to the last is at most $\lfloor\log_p(HL/g)\rfloor$. In every interval of $K$ exponents the count is at most
$$\left(1+\lfloor\log_p(HL/g)\rfloor\right)
  \left(1+\left\lfloor(K-1)/(b+1)\right\rfloor\right).$$
Within an interval, take each maximal chain separately. The first exponents of distinct chains are at least $b+1$ apart.

**An initial quadratic range.** Suppose further that $L>1$, $\gcd(H,L)=1$, $p\nmid HL$ and $(2HL)^{2(e+2n)}<p^{n^2}$. Then $1+p^e T$ satisfies the conjecture. Indeed choose $q$ nearest to $e/n$, $r=e-qn$ and $r_\pm=\max(\pm r,0)$. From $Hp^n\equiv-L$ and (2),
$$T\mid(-L)^q\Delta_e p^{r_+}-H^{q+1}p^{r_-}.$$
Its absolute value is at most $2(HL)^{q+1}p^{n/2}$ by $|\Delta_e|\le HL$, $L^q\le(HL)^q$, and $H^{q+1}\le(HL)^{q+1}$. This is at most $(2HL)^{e/n+2}p^{n/2}<p^n<T$. Equality to zero would force $L\mid H^{q+1}p^{r_-}$ when $q\ge1$, which is impossible. For $q=0$ it would give $p^e\Delta_e=H$, also impossible.

## A bound depending on the lengths of the occupied blocks

Consider an arithmetic progression of exponents with positive step $d$. Separate the nonzero base-$p$ digits of $T$ at gaps of at least $d$ zeros. Suppose there are at least two resulting blocks, and let $D_i$ be their positive integer values after shifting each block to the units position. Fix any block $i$. If $k+1$ successive exponents in this progression are exceptional, then
$$ p^{dk}\le\max_j\operatorname{lcm}(D_i,D_j). \tag{3}$$
This bound is independent of the distances between the blocks.

Read the progression backwards so successive residues are obtained by multiplication by $b=p^d$ modulo $T$. For one step write $br=uT+r'$ with $0\le u<b$. Split $r,r'$ in the same positions as $T$. At any internal boundary $p^v$ where a block starts above a gap, $T\bmod p^v<p^v/b$. The two expressions $b(r\bmod p^v)$ and $u(T\bmod p^v)+(r'\bmod p^v)$ lie in $[0,p^v)$ and agree modulo $p^v$, hence coincide. Subtraction at successive boundaries (and the integer equality at the top) gives $bR_j=uD_j+R'_j$ for each block. Thus the differences of the fractions $R_j/D_j$ between any two blocks multiply by $b$ at each step. For any fixed $i$ some block $j$ has a different initial fraction, since otherwise every fraction equals $r/T$ in lowest terms and $T$ divides $D_i<T$. Its initial separation is at least $1/\operatorname{lcm}(D_i,D_j)$ and the final separation is at most one. This proves (3).

In particular take an **odd** $T$ with at most two prime factors. If its binary expansion has multiple runs of ones, let $l,L$ be the shortest and longest run lengths. Every $l+L+1$ consecutive positive exponents contain at least one weight $1+2^eT$ satisfying the full conjecture. In (3) take $d=1$ and the shortest block. The right side is below $2^{l+L}$. If all runs have the same length $L$, every $L+1$ consecutive exponents suffice because their block values all equal $2^L-1$. If there is only one run, $T=2^L-1$ and **every exponent works**: the inverse of $2^e$ modulo $T$ is a pure power of two. This supplies a gap bound at every odd two-prime cofactor without a lower bound on the order and without knowing its factorization.

Run `python3 check_cofactor_block_inverses.py` for 3,797 exact coefficient checks, 6,778 checks of the quadratic window, and direct modular orbits. The orbit test includes 921 cases attaining a run bound of at least two for the necessary digit criterion. Lean 4.36.0 checks the signed determinant congruence and the one-step multiplication of block differences. [Formal arithmetic lemmas](BlockRatioBounds.lean)
