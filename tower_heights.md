# Almost every prime works at every finite tower height

For an odd integer $x\ge3$ set $T_1(x)=x$ and $T_{h+1}(x)=x^{T_h(x)}$. There is a set of primes of relative natural density one such that, for each prime $q$ in it, **Wasserman's full conjecture holds simultaneously at**
$$
           1+d\,T_h(q)^{2^j},\qquad
           h\ge1,\quad j\ge0,\quad 1\le d\le q^{0.08}.
$$
If $j=0$ the cofactor range improves to $d\le q^{0.11}$, still simultaneously at all heights. Both $h$ and $j$ are unbounded at each prime.

Here is a quantitative version. Put $\alpha=\log_2\lambda=0.551463\ldots$, with $\lambda^3=\lambda^2+1$, and $K=1+\lceil\log_2 X\rceil$. Fix any odd positive $b$. Among primes $3\le q\le X$ the number of possible exceptions for $1+d T_h(q)^{b2^j}$, $1\le d\le D\le X$, is
$$
 O(D^4 X^\alpha K)\quad(j=0),\qquad
 O(D^4 X^\alpha K(K+D))\quad(j\ge0).
$$
The constants are absolute, including uniformity in $b$. Thus the second assertion also works simultaneously for all odd $b\le q^\beta$ and $d\le q^\gamma$ whenever $\gamma,\beta\ge0$ and $5\gamma+\beta<1-\alpha$. This follows by summing over $b$ and comparing the count with the prime number theorem. The first assertion has $4\gamma+\beta$ in the condition.

## How many different residues do the towers give?

Write $v=v_2$, with $v(0)=\infty$. For odd $x,y$, $v(T_h(x)-T_h(y))=v(x-y)$. This holds at height one. In comparing $x^{T_h(x)}$ with $y^{T_h(y)}$, changing the base at the fixed odd exponent gives exact valuation $v(x-y)$. Changing the exponent at the base $y$ gives valuation at least $2+v(T_h(x)-T_h(y))$ since $v(y^n-1)\ge2+v(n)$ for positive even $n$. Interchanging the exponents if necessary gives the same statement for any order. These valuation facts follow by factoring a difference of odd powers, then repeatedly squaring an odd integer. We also have
$$
                 v(T_{h+1}(x)-T_h(x))\ge2h+1 .
$$
The first difference is divisible by eight and the valuation increases by at least two at each step. Thus $S_h=T_h^b$ has both properties since $b$ is odd. In particular only $O(K)$ functions $x\mapsto S_h(x)\bmod 2^K$ occur, each distinguishing the odd inputs through $X$.

For $j\ge1$ put $R=2^{j+2}$ and $Z_{h,j}=(S_h^{2^j}-1)/R$, an integer. Within each odd class modulo four,
$$
 v(Z_{h,j}(x)-Z_{h,j}(y))=v(x-y)-2\quad(x\ne y).
$$
Indeed $S_h(x)\equiv x\pmod4$ and for odd $A,B$ in the same class modulo four, $v(A^{2^j}-B^{2^j})=v(A-B)+j$. At fixed $x$ the height increment in $Z$ has valuation at least $2h-1$. Also
$$
                  Z_{h,j+1}=Z_{h,j}+2^{j+1}Z_{h,j}^2 .
$$
It follows that modulo $2^K$ the functions cease changing in either index once that index exceeds $K+1$, uniformly in the other index. In each class modulo four every individual function still distinguishes the inputs through $X$.

## Counting rows after removing their low digits

We use the repository's [binary carry-prefix theorem](binary_affine_counts.md#counting-low-digit-strings). For three positive slopes $s_i$, include their sum $s_0$ and set $w=\lambda^{\max_i v(s_i)-\min_i v(s_i)}$. For three affine parts with these slopes and nonnegative offsets, at most $w2^{\alpha K}$ strings of $K$ input bits extend to a carry-free addition. This bound holds after any fixed low prefix. Removing that prefix just changes the offsets to the three quotients rounded down, with no incoming carry in a successful addition. Also $\sum_{s_1+s_2+s_3\le D}w\ll D^3$. Indeed $w\le\sum_{i=0}^3 2^{\alpha v(s_i)}$, each value occurs $O(D^2)$ times in each position, and $\sum_{t\le D}2^{\alpha v(t)}\ll D$.

Put $Q=T_h(q)^{b2^j}$. By the [endpoint and factor assignment](moving_power_bases.md#summing-over-growing-cofactors), a failure requires a binary carry-free row of one of the forms
$$
 A=(aQ,a'Q,cQ+1),\quad a+a'+c=d,\qquad
 C=V(f_1,f_2,f_3),\quad U=\sum f_i\le D,\quad UV=1+dQ .
$$
All slopes are positive. In using $C$ we may assume $A$ carries and $U$ is odd. The sum of the weights $w$ over all choices of $A$ is $O(D^3)$, and for all choices of $C$ and $d$ it is $O(D^4)$.

For $j=0$, each of the $O(K)$ tower functions distinguishes the $K$ low bits of $Q$ on the inputs in question. For $C$ the same bits distinguish the next $K$ bits of $V$ after its fixed $v(d)$ low bits, since $U$ is odd. This gives the first bound by the prefix theorem. Congruences with an odd denominator make sense even on inputs that fail the integrality test. An actual successful integer always gives an extendible prefix.

For $j\ge1$ write $Q=1+RZ$. Deleting $j+2$ low bits in $A$ gives slopes $a,a',c$ on input $Z$ and offsets $\lfloor a/R\rfloor,\lfloor a'/R\rfloor,\lfloor(c+1)/R\rfloor$. These are fixed at each $j$ and all vanish once $R>D+1$.

For $C$, set $v_0=[(1+d)/U]_R$ (least nonnegative residue) and $u=(1+d-Uv_0)/R$. Then $V=v_0+RY$ with $Y=(dZ+u)/U$. Deleting the same low bits gives the three parts $f_iY+\lfloor f_i v_0/R\rfloor$. Once $R>D(D+1)$, we have $-U<u\le0$ and
$$
                \left\lfloor f_i v_0/R\right\rfloor
                   =\left\lfloor -f_i u/U\right\rfloor .
$$
In fact the numerator perturbation $f_i(1+d)/R$ lies in $[0,1)$ before division by $U$. There are therefore at most $U$ choices of the whole tuple of $u$ and the offsets throughout this tail. For a fixed $u$ the bits of $Z$ distinguish the next $K$ bits of $Y$ above its common low $v(d)$ bits. The carry-prefix theorem again bounds the number of acceptable inputs per function and tuple by $w2^{\alpha K}$, up to the four residue classes of $x$.

There are $O(K)$ height choices at precision $K$. For $j$ through $\max(K+1,\lceil\log_2(D(D+1))\rceil)$ there are $O(K)$ choices with fixed offsets. Beyond this point there is a single $Z$ function per height and at most $D$ tuples for $C$, or one for $A$. Multiplying by the summed weights proves the second bound. In particular passing all low-digit tests was used only as a necessary condition for a counterexample.

Run `python3 check_tower_heights.py` for independent modular tower evaluations, distance and stabilization checks, and comparisons of the high-digit formulas with integer parts and exact coefficients. The verifier makes 172,638 distance comparisons and checks 2,594 exact multinomial coefficients.

## Every prime through 10,000

For **every** prime $q\le10000$, every $1\le d,b\le1000$ and every $a\ge0$, the conjecture holds at
$$
                           N=1+d q^{\,bq^{q^a}}.
$$
Already at $N=1+d q^{\,bq^a}$ the only triples $(d,q,b)$ still possible are
$$
        (629,4943,97),\quad(629,4943,685),\quad
        (986,2447,333),\quad(986,2447,687).
$$
In particular every $b\le96$ works at this level too. Setting $b=1$ proves the assertion at $1+dT_h(q)$ for all $h\ge1$ throughout this explicit prime range. The [square theorem](common_power_exponents.md) covers every further squaring.

We use the [13 remaining small prime–cofactor pairs](prime_base_density.md#checking-all-small-prime-bases) and their [necessary row shapes and residues](shifted_directions.md). For each shape write $A=(s_1 W,s_2 W,s_3 W+1)$, $C=V(f_1,f_2,f_3)$ with $U=\sum f_i$ and $UV=1+dW$. The initial class $W\equiv r\pmod M$ and a list of primes that $A$ must avoid come with the shape. The two possible odd shapes at 986 are both tested. The prior square theorem handles even $b$, and its cube counterpart handles $3\mid b$ when $d\ne986$.

At a tested layer $w=\ell^k$ the first row is determined modulo $w$. Set $v=\gcd(w,U)$ and take $w$ large enough that $\ell\mid w/v$ and $\ell\nmid U/v$. If $v\nmid1+dW$ the shape is impossible. Otherwise $V$ is determined modulo $w/v$ by dividing the congruence by $v$ and inverting $U/v$. If $\ell$ divides both $U$ and $V$, discard the input by coprimality. If $\ell\mid U$ or the first row must avoid $\ell$, a low carry in $A$ excludes the input. A low carry in $C$ excludes it if $\ell\mid V$. When $\ell\mid d$, both rows carrying excludes the input. Each rejection uses digits fixed by the chosen modulus.

For each initial modulus or digit test modulus $w$, let $L$ be the order of $q$ modulo $w$ and $P$ the order of $q$ modulo $L$. All required coprimalities hold in the finite table. Use order one at modulus one. At $W=q^{bq^a}$ its residue is determined by $a\bmod P$, computed by two modular powers. At $W=q^{bq^{q^a}}$, one more power suffices and the period in $a$ divides the order of $q$ modulo $P$. These periods apply starting at $a=0$.

For each $b$ the verifier discards its class if a test admits no index, or if two sets of indices are disjoint modulo the gcd of their periods. At $d=986,q=2447,b=891$ a third test is needed for the first row with slope 119. The moduli $4096,2401,24389$ have index periods $16,294,2436$ and direct enumeration of their common period $68,208$ leaves no index. The resulting four triples at the first level are all excluded by the second-level test.

Run `python3 check_small_prime_towers.py`. For each encountered residue its congruence and carry decisions are compared with an independently chosen integral lift of $W$ and factorial floors on the two reconstructed rows. The program also compares all modular orbit entries with direct exponentiation.
