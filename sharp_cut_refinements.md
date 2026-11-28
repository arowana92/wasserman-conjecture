# The exact predecessor divisor after splitting a binomial cut

For coprime positive $a,b$, fix $D\mid b$. Take **all** positive rows
$$
  (aV,s,bV-s),\qquad V\ge1,\quad 0<s<bV,\quad \gcd(s,b)\mid D.
$$
Write $N=(a+b)V$ and let $K$ be the coefficient of the row. The least positive integer $C_D$ making $N-1\mid C_D K$ throughout this class is
$$
 C_D=\frac{C(a,b)}{\gcd(C(a,b),b/D)},                                   \tag{1}
$$
where $C(a,b)$ is the [sharp constant for a pair of parts in ratio $a:b$](part_ratio_predecessor.md). In particular $D=b$ imposes no restriction on the split, and $D=1$ removes every prime factor of $b$ from the constant. The split may change with $V$.

## Binomial cut and sharpness after refinement

Put $M=N-1$, $B=\binom{N}{aV}$ and $F=\binom{N-2}{aV-1}$. Then
$$
 \frac{B}{M}=\frac{a+b}{abV}F,\qquad K=B\binom{bV}{s}.
$$
A prime $p\mid M$ is coprime to $V(a+b)$. Thus the defect $v_p(M)-v_p(B)$ is $v_p(ab)-v_p(F)$. More generally for $p\mid ab$ the same equality holds at every unit $V$, while a nonunit $V$ gives no positive defect. This is exactly the minimization formula for $c_p=v_p(C(a,b))$. In particular $C(a,b)\mid ab$ and $M\mid C(a,b)B$ at all scales.

For $p^k\parallel b$, put $\ell=v_p(D)$. At a unit $V$, the identity $s\binom{bV}{s}=bV\binom{bV-1}{s-1}$ gives at least $k-\ell$ additional factors of $p$ when $\ell<k$. It follows that the maximum positive defect is bounded by $(c_p-k+\ell)_+$, where $x_+=\max(x,0)$. At $p\mid a$ the binomial bound $c_p$ still applies. These are the exponents on the right of (1).

To prove sharpness at any prime with a positive proposed exponent, take a unit $V$ attaining the maximum binomial defect $c_p$. Let $h$ be the multiplicative order of $p$ modulo $a+b$, and set $W=(p^h-1)/(a+b)$. Take $p^L>(a+b)V$ and replace $V$ by $V'=V+p^LW$. The base-$p$ additions for the two separated blocks are independent. At $W$ the upper binomial argument has all digits $p-1$, so the upper block adds no carries. The valuation of $N-1$ also stays fixed. If $p\mid b$, take $s=p^\ell$. Its gcd with $b$ is allowed. The inner binomial has valuation exactly $k-\ell$, since $\binom{bV'-1}{p^\ell-1}$ is a unit by its low $k$ digits. If $p\mid a$, take $s=p^L$. The digit at position $L$ in $bV'$ is nonzero since $p\nmid bW$, so $\binom{bV'}{s}$ is a unit and $\gcd(s,b)=1$. Both splits have positive remaining part for sufficiently large $L$. This gives individual witnesses for every prime exponent in (1).

Every coefficient in the class contains $a+b$ by factorial decrement in the outer binomial. If $C_D\le a+b$, then $C_D<a+b$ by coprimality. The coefficient therefore contains a divisor of $N(N-1)$ of size at least $(a+b)(N-1)/C_D\ge N$. Any two rows meeting this bound at the same weight share an endpoint prime, since $N^2>N(N-1)$. The choices of cuts, terms and splits can differ between rows. A row containing two parts in one of the previously proved [integer or geometric ratios](part_ratio_predecessor.md) also meets the bound. Indeed if its primitive pair is $(ag,bg)$ then the missing predecessor factor divides $gC(a,b)$ by applying the pair theorem at scale $V$ and cancelling $V$, which is coprime to $N-1$. This is less than its primitive sum when $C(a,b)\le a+b$.

## Two families with full predecessor divisibility

1. Suppose $P=p^k$, $k\ge1$, $\gcd(b,c)=1$ and $p\nmid bc(b+c)$. Every coefficient on the ray $(b+c,Pb,Pc)$ is divisible by $N-1$. Its primitive pair-gcd product is $P$, so only $p$ can cause loss. At every scale the row $(P(b+c),b,c)$ has exactly the same weight and $p$-valuation by the base-$p$ digit-sum formula. This second direction has pairwise coprime entries so its full predecessor divides its coefficient. For example $(2,P,P)$ works for every odd prime power, however large $P$ is.
2. Put $a=(p^m-1)/(p-1)$ and $b=p^m$, for $m\ge1$. For **every** $1\le x<b$ the ray $(a,x,b-x)$ has
   $$
      N-1\mid \gcd(a,x(b-x))K.
   $$
   The binomial cut has $c_p=1$ by the [geometric ratio formula](part_ratio_predecessor.md#shifting-prime-powers-between-the-parts). When $p\nmid V$, splitting $bV$ at $xV$ supplies at least $m-v_p(x)\ge1$ factors of $p$. At $p\mid V$ the predecessor is already a unit. For any other prime, the missing power divides the [primitive pair-gcd product](sharp_direction_predecessor.md). Its prime-to-$p$ part is $\gcd(a,x(b-x))$, since $\gcd(x,b-x)$ is a power of $p$ and $\gcd(a,b)=1$. In particular $M$ divides $K$ outright whenever $\gcd(a,x(b-x))=1$. This includes $(2^m-1,2^{m-1},2^{m-1})$ at every $m$ and $(a,p^{m-1},(p-1)p^{m-1})$ whenever $\gcd(m,p-1)=1$.

For an example where the second conclusion improves a large coarse bound, let $P=2^k$, $k\ge1$, and take
$$
 \big(P(P-1)/2,\ P(P+1)/2,\ P^2-1\big).
$$
The first two terms sum to $P^2$. The primitive pair-gcd product is $P(P^2-1)/2$, but the sharp predecessor constant divides $P^2-1$, less than half of the primitive sum $2P^2-1$.

## An arbitrary companion when the full predecessor divides

Let a primitive direction have sum $u$ and suppose its scaled coefficient is divisible by $M=N-1$ at the weight under consideration. If a second row of primitive sum $v$ has no common endpoint prime with it, then
$$
                v\ge27u^2+1,\qquad N\ge 27u^3+u .
$$
Indeed $\gcd(u,v)=1$, $uv\mid N$, and every prime of $M$ must be avoided by the second row. The retained-power lemma and AM–GM give $27M^2\le v^3$, whence $27u^3\le N^3/(N-1)^2<N+3$ for $N\ge5$. As $u\ge3$ divides both $N$ and $27u^3$, this forces $N\ge27u^3$. Next $v\ge27u^2(1-1/N)^2>27u^2-1$ because $54u^2/N\le2/u<1$. Coprimality excludes $v=27u^2$. For $N<5$, the condition $uv\mid N$ is already impossible.

Run the [checker](check_cut_refinements.py) with Python 3. It constructs sharp prime witnesses by the separated-block argument, checks (1) on exact coefficients and checks the two families against factorial valuations and the independent direction algorithm.

~~~sh
python3 check_cut_refinements.py
~~~
