# Short cofactor intervals above multiples of a prime power

Let $P=p^a$ and $L=p^j$ with $j\ge a\ge1$, and put $T=L/P$. Write $\omega(W)$ for the number of distinct prime factors of $W$.

**Square-root gap.** If $W=L+d$, $0<d<L$, and $\omega(W)\le2$, the full three-row conjecture holds at $1+PW$ whenever
$$
                         d^2<4W+T.
$$
A sharper necessary inequality at a failure is $d^2\ge(U+[W]_U)W+T$, where $U\ge3$ is the weight factor assigned to the row taking all three predecessor bases and $[\,\cdot\,]_U$ denotes least nonnegative residue.

**Uniform bound in aligned intervals.** Let $b=1$, or $b=2$ with $p$ odd. For any positive $h$ with $P\nmid h$ and $0<bH<L$, consider all $1\le d\le H$ with $\omega(hL+d)\le2$. The number of these $d$ for which $1+bP(hL+d)$ can fail the conjecture is
$$
 \ll_\epsilon (1+bH^2/L)\,[P(1+bH^2/L)]^\epsilon .
$$
The constant is independent of $h,p,a,j,H$. When $P\mid h$ one can absorb the entire $p$-part of $h$ into $L$ before applying the bound. For $H=C\sqrt L$ it suffices to inspect the divisors of $O(1+bC^2)$ specified nonzero integers, however large $h$ is.

## A signed divisor condition

The [static numerator theorem](static_square_cofactors.md#at-most-two-possible-numerators) gives the following necessary conditions at a failure with $W=hL+d$. We have $p\nmid W$, and a row avoiding $p$ has a part $xW$ with $3\le x<bP$, $xW\equiv1\bmod P$ and $\gcd(x,1+bPW)\ge3$. These use the [all-weight endpoint reduction](all_weights_endpoint_reduction.md). The calculation below applies for any positive $b$ whenever such a part is given.

Set $e=(xd-1)/P$, so $0\le e<bd$. The low digits of the total through modulus $PL$ form $1+Pbd$. The low digits of $xW$ form $1+Pz$, where
$$
 z\equiv e+Thx\pmod L,\qquad 0\le z<L .
$$
Write $z\preceq_p bd$ if each base-$p$ digit of $z$ is at most the corresponding digit of $bd$. This is exactly Kummer's test for the cut through these positions. [Bergman, Lemma 5](https://arxiv.org/abs/0806.0607) It implies, for some signed integer $t$,
$$
 z=e+Tt,\qquad t\equiv hx\pmod P,\qquad dt\equiv h\pmod P,\qquad
 |dt|\le bPd^2/L. \tag{1}
$$
Indeed $0\le z\le bd$ and $0\le e<bd$. Conversely, for each admissible $x$, any $t\equiv hx\bmod P$ giving $0\le e+Tt<L$ recovers the actual low digits. Thus the digitwise test after (1) is exact for the truncated binomial cut.

Take $h_0$ to be a least absolute representative of $h\bmod P$. Every surviving $d$ divides some $h_0+iP$ with $|h_0+iP|\le bPH^2/L$. These integers are nonzero by hypothesis. The elementary divisor bound gives the stated count. If $bH<T$ then $t=0$ in (1), contradicting $P\nmid h$. In particular, if $bH^2<L(1-|h_0|/P)$, only divisors of $|h_0|$ need testing. This argument also applies at interior digit positions, with unrestricted digits above $h_0$.

## The leading block and the two other parts

For the square-root gap take $b=h=1$. At a failure there are exactly three predecessor bases, so $W=QR$ with $Q,R$ full powers of distinct primes, coprime to $P$. The forced row has parts $xW,k_QPR,k_RPQ$, where $1\le k_Q<Q$, $1\le k_R<R$ and $U$ divides all three multipliers $x,k_Q,k_R$. Since $N=1+PL+Pd<2PL$, one part receives the leading power $PL$. [Three-factor CRT row](endpoint_factors.md#when-n-1-has-three-distinct-prime-factors)

If $xW$ receives it, then $xW\ge PL+1$ by its residue. The two other parts give $d\ge k_QR+k_RQ\ge2U\sqrt W$. This also implies the sharper claimed inequality since $[W]_U<U$ and $T<W$.

Otherwise exchange $Q,R$ so that $k_QPR\ge PL$ and put $\delta=Q-k_Q\ge1$. Thus $\delta R\le d$. The other two parts sum to $P\delta R+1$. Using $xd=1+Pe$ gives
$$
 d\,\delta R=eW+T+dQk_R,\qquad
 d^2\ge(e+\delta k_R)W+T.
$$
Here $dQk_R\ge\delta RQk_R$. Also $e>0$ since $x\ge3$. As $U\mid x,N$ and $N=1+PW$, the two congruences $Pe\equiv PW\equiv-1\bmod U$ give $e\equiv W\bmod U$. Since $(U,W)=1$, we have $e\ge[W]_U\ge1$. This proves both bounds.

## The next block and a square power

Let $h=1+kP$, $k\ge0$. For a cut passing through modulus $PL$ set $f=\lfloor x(L+d)/(PL)\rfloor$ and assume $bd<L$. Then $0\le f\le b$. After deleting the low digits through $PL$, the two summands and their total are
$$
 kx+f,\quad k(bP-x)+b-f,\quad bkP+b.
$$
For $b<p$, a positive $k$ divisible by $p$ is impossible. At the first nonzero position of $k$ the small offsets contribute no digit, and the two nonzero digits sum to a multiple of $p$. If $p\nmid k$, the next block therefore requires
$$
 k\equiv(i-f)d\bmod P,\qquad 0\le i\le b,\quad i\ne f. \tag{2}
$$
This leaves at most $b$ residue classes for all positive block indices. In particular every positive $k$ divisible by $p$ certifies the whole interval $1\le d<L/b$.

There is an explicit boundary at **every** square prime power $P=S^2$ with $S=p^c\ge7$. Take $L=P$, $b=1$, $1\le d\le S+\lfloor\sqrt S\rfloor$. A failure requires all of
$$
 d=S+1,\quad k>0,\quad k\equiv-S-1\bmod P,\quad
 \gcd(k-1,S^2-S+1)>1. \tag{3}
$$
In fact $|dt|\le d^2<2P-1$ in (1), with the strict inequality valid for $S\ge7$. Since $dt\equiv1\bmod P$, the value 1 would force $d=x=1$. The other possibilities give $dv=P\pm1$ for $v=|t|$. For the positive sign, $x=v$ and $e=1$. For the negative sign, $x=P-v$ and $e=d-1$. The bounds on $z$ give $1+v\le d$ in both cases. Hence $d=S+u$ with $1\le u\le\lfloor\sqrt S\rfloor$, using $p\nmid d$. Put $v=S-u+i$. Then $(S+u)i=u^2\pm1$. Its right side lies in $[0,S+u)$, so it vanishes. Only $d=S+1$, $v=S-1$ with the negative sign survives.

Here $x=P-S+1$, $x(P+d)=P^2+P+1$ and $f=1$. The truncated cut does pass at $p$, and (2) gives the class of $k$. Reduction modulo $x$ gives $1+P((1+kP)P+d)\equiv k-1\bmod x$, since $S^3\equiv-1\bmod x$. The weight-scale condition proves (3) and rules out $k=0$. Writing $k=tP-S-1$, any common scale also divides $t^2+t+7$. Indeed $(t-1)S\equiv t+2$ modulo that scale. Multiply $S^2-S+1\equiv0$ by $(t-1)^2$ and substitute.

## Two digit blocks and independently varying exponents

Put $D=p(p+1)/2$. Among unit residue classes $W\bmod P^2$, the exact number passing the necessary inverse cut through modulus $P^2$ for $b=1$ and $x\ge3$ is
$$
 (D-1)D^{a-1}-2-
 \begin{cases}
 0 & p=2,\\
 \dfrac{p+3}{2}\left(\dfrac{p+1}{2}\right)^{a-1} & p>2.
 \end{cases}
$$
Indeed write $W=d+Ph$ modulo $P^2$, $0<d<P$. For the inverse $x$ of $d$, the higher residue of $xW$ runs through the permutation $xh+(xd-1)/P\bmod P$ as $h$ varies. Precisely $\prod_i(d_i+1)$ residues fit the base-$p$ digits $d_i$ of $d$. Summing over a nonzero units digit and unrestricted higher digits gives $(D-1)D^{a-1}$. The disallowed numerator 1 contributes two. For odd $p$ the numerator 2 corresponds to $d=(P+1)/2$, giving the last term. The count is exact for this digit test.

Now fix a prime $r\ne p$ and a prime power $Q$ coprime to $p$. For $W=Qr^e$, choose $e=e_0+s i$ with $s\ge1$, and set $v=v_p(r^s-1)$. Assume $1\le v\le a$, with $v\ge2$ when $p=2$. For every integer $a\le t\le2a-v$, the fraction of possible failures at $1+bPQr^e$ among **any $p^t$ consecutive indices** with positive exponents is at most
$$
                 b\left(\frac{p+1}{2p}\right)^{t-a}. \tag{4}
$$
The estimate is uniform in $Q$, including its exponent. For instance at $1+2^a Q\,3^e$, $a\ge3$, each block of $2^{2a-3}$ exponents of one parity has possible failure fraction at most $(3/4)^{a-3}$, below $2.4\cdot10^{-8}$ at $a=64$ and $2.5\cdot10^{-16}$ at $a=128$.

To see this, raising $1+p^v u$, $p\nmid u$, to the $p$th power gains exactly one place under the stated hypotheses. This follows by binomial expansion and iteration. Thus in the block of indices $W$ runs once through all lifts of one unit residue modulo $p^v$ to modulus $p^{v+t}$. Put $\lambda=v+t-a$, so $v\le\lambda\le a$. Write each lift as $d+Ph$ with $0<d<P$ and $0\le h<p^\lambda$. There are at most $b$ choices of $x$ for each $d$. The high $\lambda$ digits of $xW$ above the $P$ shift are again a permutation in $h$, and must fit the low $\lambda$ digits of $bd$. Multiplication by $b$ permutes the lifts of the fixed units class onto those of another fixed class, since $p\nmid b$. The bottom $v$ digits of $bd$ contribute at most $p^v$ choices for the high residue. Summing over the $\lambda-v$ free target digits gives $D^{\lambda-v}$, and the remaining $a-\lambda$ digits of $d$ give $p^{a-\lambda}$. Divide the resulting bound $b p^v D^{\lambda-v}p^{a-\lambda}$ by $p^t$ to get (4). Each uncounted class certifies every eligible weight in that class.

The same index can certify many different $p$-exponents at once. Write $\eta=(p+1)/(2p)$. Fix $A\ge v$ and $t\ge2A-v$. In any $p^t$ consecutive indices, the proportion that can fail at **any** $A\le a\le t-A+v$ is at most $2b\eta^{A-v}/(1-\eta)$. For $2a-v\le t$ the interval is a union of full periods and the bound is $b\eta^{a-v}$. For the larger $a$ use (4) to get $b\eta^{t-a}$. Summing the two geometric series proves the assertion without a factor for the length of the exponent range. For instance, with $p=2,r=3,s=2$, the fraction is below $1.92\cdot10^{-7}$ simultaneously for every $64\le a\le t-61$ in every block of $2^t$ exponents of either parity, for all $t\ge125$ and any fixed odd prime power $Q$.

A single set of exponent classes can certify an unrestricted second prime base. For $a\ge3$ put $M=2^{2a-2}$. There is an explicit set of **$M-6\cdot3^{a-3}+2$ classes** $e\bmod M$ for which the conjecture holds simultaneously at
$$
                1+2^a\,3^e q^{Mj}
        \qquad(q\text{ any odd prime},\ e,j\ge1).
$$
Every odd $q$ has $q^M\equiv1\bmod 2^{2a}$, by squaring from modulus 8. The powers of 3 run once through the odd residues congruent to 1 or 3 modulo 8, since $v_2(3^2-1)=3$. Thus the classes to retain can be computed using $W=3^e\bmod 2^{2a}$ alone. The fixed three low bits of $d$ give respectively 2 or 4 allowed subwords, and the remaining $a-3$ digits sum to $3^{a-3}$ in each case. Subtract the two words with numerator 1. Exactly $6\cdot3^{a-3}-2$ classes pass the truncated cut. The complementary classes prove the assertion, with no condition on $q$ or $j$ beyond those displayed.

## Larger checked intervals

For the following parameters take all $W=(1+kP)L+d$ with $1\le d\le H$ and $\omega(W)\le2$. Every weight at $k=0$ satisfies the conjecture. At positive $k$, any failure must use one of the indicated tails and satisfy (2).

| $P$ | $b$ | $L/P$ | $H$ | Pairs passing the layer inequality | Pairs passing all low digits |
| --- | ---: | ---: | ---: | ---: | ---: |
| $2^{64}$ | 1 | 1 | $2^{39}$ | 9,008 | 3 |
| $2^{128}$ | 1 | 1 | $2^{69}$ | 2,007 | 2 |
| $3^{40}$ | 2 | 3 | 96,628,603,792 | 803 | 8 |

At $2^{64}$ the three tails are $4{,}294{,}967{,}297$, $150{,}323{,}855{,}353$ and $244{,}813{,}135{,}853$. At $2^{128}$ they are $18{,}446{,}744{,}073{,}709{,}551{,}617$ and $387{,}381{,}625{,}547{,}900{,}583{,}929$. In both cases $f=1$ at every hit. Hence at all other classes of $k\bmod P$ than the negatives of these tails, the whole window is certified.

Before the second tail at $2^{64}$, (3) continues to hold. The number $x=18{,}446{,}744{,}069{,}414{,}584{,}321$ is prime. An exception anywhere in the window $d\le150{,}323{,}855{,}352$ thus requires **one class modulo $2^{64}x$** for $k$, from $k\equiv-(2^{32}+1)\bmod 2^{64}$ and $k\equiv1\bmod x$. At $2^{128}$ the corresponding multiplier has prime factors $769$ and $442499826945303593556473164314770689$.

The verifiers enumerate divisors in (1) and, separately, positive factor pairs of $mP\pm1$ followed by exact factorial-floor tests. Both validate their factorizations. An exhaustive direct scan over the tails in smaller windows also tests prefixes other than 1. The ternary case records all eight pairs in the verifier. At $k=0$ only one ternary pair has the required common weight scale, and its $W$ has four distinct prime factors. The checks use Python with python-flint and Lean 4.36.0. [Divisor verifier](check_predecessor_cofactor_windows.py) · [Factor-pair check](check_predecessor_window_cross.py) · [Lean arithmetic lemmas](PredecessorCofactorWindows.lean)

~~~sh
python3 check_predecessor_cofactor_windows.py
python3 check_predecessor_window_cross.py
lean +leanprover/lean4:v4.36.0 PredecessorCofactorWindows.lean
~~~
