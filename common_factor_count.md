# Counting weights from factors shared by the parts of a row

The number of weights $N\le X$ at which Wasserman's three-row conjecture can fail is
$$ O(X^{0.8946}). $$
The [previous global bound](binary_windows.md) was $O(X^{0.9329})$. There is also a uniform bound in intervals: if $I\subset[X,2X]$ has $D\ge X^{37/50}$ consecutive integers, at most $O(DX^{-1/100000})$ weights in $I$ can fail. A more general version below allows fixed low binary digits as well as fixed high digits.

Two stronger questions admit estimates by the same argument. Call a weight exceptional for a pair of row lengths if some two positive rows of those lengths have coefficients sharing no prime divisor of $N(N-1)$. For lengths two and three the number through $X$ is $O(X^{0.962})$. For two binomial coefficients it is $O(X^{0.978})$, improving the [exponent 0.993](binomial_endpoint_density.md). These assertions concern all choices of rows at each weight outside the exceptional set.

## Counting small assigned moduli using one predecessor factor

Put $Q=2^L$, work in $Q/2\le N<Q$, and choose rows $A,B$ with no common endpoint prime. Choose $A$ avoiding 2; assign 2 to it. Split the full prime powers of the weight and predecessor as $UV=N$, $ST=N-1$, assigning $U,S$ to $A$ and $V,T$ to $B$. These are the [endpoint assignments](endpoint_reduction.md). The [endpoint theorem](all_weights_endpoint_reduction.md) makes the case of two ternary rows equivalent to the original conjecture. No part in either row equals one: such a row has coefficient divisible by $N$, and every other positive row has coefficient sharing a prime with $N$, since $N\mid B_i K_B$ for each part.

Put $K=US$. **If $B$ has three parts, the number of weights with such an assignment and $K\le Y$ is**
$$ \ll Q^{-2/5+o(1)}Y^{8/5},\qquad 1\le Y\le Q. \tag{1} $$
Here and below $Q^{o(1)}$ absorbs divisor bounds on nonzero integers of size $Q^{O(1)}$ and powers of $\log Q$.

Divide $B$ by $V$. At each full power in $T$, two normalized parts vanish. Group the powers into three pairwise coprime bins by the index of the unique nonvanishing part. Choose a smallest bin $d$, so $T=dt$, $d\le T^{1/3}$. Its designated normalized part is $zt$ with $z>0$, and the other two sum to $dy$ with $y>0$. We have
$$ U=zt+dy,\qquad V(zt+dy)=Sdt+1,\qquad Vzt\equiv1\pmod d. \tag{2} $$
The last congruence is for the actual part of $B$. In particular $Vz$ is invertible modulo $d$ (with the condition empty at $d=1$).

Use dyadic boxes $U\asymp U_0$, $S\asymp S_0$, $d\asymp d_0$; write $V_0=Q/U_0$, $T_0=Q/S_0$, $K_0=U_0S_0$. Then $z\ll U_0d_0/T_0$ and $y\ll U_0/d_0$. A nonempty box permits these bounds on the counts of the positive parameters, even near one.

1. Fix $S,V,d,z$. Equation (2) determines $t$ modulo $d$, so $N=1+Sdt$ lies in one class modulo $Sd^2$. It is also zero modulo $V$, with $\gcd(V,Sd)=1$. There are $O(1+Q/(S_0V_0d_0^2))$ possible weights. Counting the fixed parameters gives
$$ \ll S_0V_0d_0(U_0d_0/T_0)\big(1+Q/(S_0V_0d_0^2)\big)
 \ll Q^{2/3}S_0^{4/3}+K_0. $$
2. Fix $S,d,z,y$. Eliminating $t$ from (2) gives $U\mid Sd^2 y-z$. This is positive: $(Sd-Vz)t=Vdy-1>0$ since $V\ge2$. Each divisor $U$ determines $t=(U-dy)/z$ and then $N$. Counting these parameters, using $d_0\ll T_0^{1/3}$, gives
$$ \ll Q^{o(1)}S_0d_0(U_0d_0/T_0)(U_0/d_0)
 \ll Q^{-2/3+o(1)}K_0^2 S_0^{-1/3}. $$

Combine the two bounds, taking powers $1/5,4/5$ of the terms with opposite $S_0$ exponents and retaining the additive $K_0$. Since $T^2$ divides the product of the three normalized parts, $27T^2\le U^3$ and $K\ge3(N-1)^{2/3}$. In a nonempty box the resulting $Q^{-2/5}K_0^{8/5}$ therefore absorbs $K_0$ up to a constant. Summation proves (1). For a two-part complementary row the [existing count](binomial_endpoint_density.md#count-when-the-modulus-specified-by-a-cut-is-small) is $\ll Q^{-1/3+o(1)}Y^{5/3}$; its proof only uses the assigned complementary row and the two-part factorization there.

## Information that every individual part supplies

Choose one canonical witness and assignment per exceptional weight, fixing the row lengths and the length $a\in\{2,3\}$ of the binary carry-free row. On any given stratum of weights choose $N$ uniformly; let $P=(U,S)$ and $Z_i=A_i$. Use base-two Shannon entropy, written $\mathcal H$. Each $Z_i\ge2$ satisfies $U\mid Z_i$ and $S\mid Z_i(Z_i-1)$. Thus $\mathcal H(P\mid Z_i)=o(L)$ **for every** $i$. Since $P$ is determined by $N$ and $N=\sum Z_i$, conditional subadditivity gives
$$ \mathcal H(N)\le \frac1a\sum_{i=1}^a\mathcal H(Z_i)
              +\left(1-\frac1a\right)\mathcal H(N\mid P)+o(L). \tag{3} $$
Indeed $\mathcal H(Z_i)=\mathcal H(P)+\mathcal H(Z_i\mid P)-\mathcal H(P\mid Z_i)$ and $\sum_i\mathcal H(Z_i\mid P)\ge\mathcal H(N\mid P)$. The coefficient $1-1/a$ is important when $K$ is below the length of the interval being counted.

Write $H(v)=-v\log_2v-(1-v)\log_2(1-v)$ with its continuous endpoint values. For a block of $w$ binary positions in which $N$ has $vw$ ones, the mean over $i$ of the block entropy in $Z_i$ is at most $wH(v/a)$ by subadditivity and concavity. This uses carry freedom: the parts partition exactly those ones. A block known in advance to contain $b$ prescribed ones costs at most $bH(1/a)$ in that average.

Take $K\ge2^{\lfloor kL\rfloor}$, $h=L-\lfloor kL\rfloor$. Fix $e=v_2(N(N-1))$, $\eta=N\bmod2$ and the slot of the units bit if present. All low $e$ bits of the parts are now fixed. The word $n=(N-\eta)/2^e$ has length $L-e$ with padding. Given $P$, it is in a single residue class modulo the odd number $K/2^e$. Either its first or last $h$ bits determine it: the first restrict it to an interval of length $2^{L-e-h}\le K/2^e$, and the last increase its modulus to $2^hK/2^e\ge2^{L-e}$. If $L-e<2h$ the unrestricted count is already at most $2^{2h}$.

Otherwise stratify by the densities $r_1,s,r_2$ of ones in its first $h$, middle $m=L-e-2h$, and last $h$ positions (omit an empty block). Put $r=(r_1+r_2)/2$. The smaller of the two end entropies bounds $\mathcal H(N\mid P)$ and is at most $hH(r)$. The ordinary binomial count and (3) now give, respectively,
$$
\mathcal H(N)\le 2hH(r)+mH(s),\qquad
\mathcal H(N)\le 2hH(r/a)+mH(s/a)+(1-1/a)hH(r)+o(L).
$$
Take a fraction $\lambda$ of the second bound and $1-\lambda$ of the first. Since $h=(1-k)L+O(1)$, $m\le(2k-1)L$, the coefficient of $L$ is at most
$$ (1-k)\big[(2-\lambda(1+1/a))H(r)+2\lambda H(r/a)\big]
 +(2k-1)\big[(1-\lambda)H(s)+\lambda H(s/a)\big]. \tag{4} $$

| Binary row $a$ | Other row | $k$ | $\lambda$ | Tangencies $r_0,s_0$ | Upper exponent |
|---|---|---:|---:|---|---:|
| 3 | 3 | .80910 | .70108 | .63946, .71319 | .8946 |
| 2 | 3 | .85117 | .55158 | .59869, .62045 | .962 |
| 3 | 2 | .75132 | .67204 | .63107, .69422 | .919 |
| 2 | 2 | .78673 | .46489 | .57935, .59279 | .978 |

In each row the small-modulus exponent from (1) or the binomial count is strictly below the last entry. For (4) replace each entropy by its affine upper tangent at $r_0,r_0/a,s_0,s_0/a$, as appropriate. The verifier checks the four corners of the resulting affine expression using rational logarithm enclosures. The strict margins absorb $o(L)$, the polynomial number of strata and both choices for the binary row. Summing the dyadic blocks proves the stated global bounds.

## Uniformity with prescribed bits at both ends

In any $L$-bit dyadic range, fix arbitrary high and low binary strings leaving one consecutive block of $n=\lfloor(37/50)L\rfloor$ free bits. At most $O(2^{(37/50-1/100000)L})$ weights in this set can fail the conjecture, uniformly in both strings and the location of the block. In particular the low string may prescribe an arithmetic progression modulo a power of two.

Put $\theta=37/50$, $k=.7124$. Equation (1) for the entire dyadic range bounds the small-$K$ contribution by $Q^{8k/5-2/5+o(1)}$. For the rest fix $e,\eta$ and the units slot as before. If the low prescribed string has length $b$, let $c=\min(b,e)$. The low $c$ bits of the parts are fixed by carry freedom. On the other prescribed positions, at most $L-n-c$ ones may be occupied. For fixed $P$ the weight lies in a single class modulo $\operatorname{lcm}(2^b,K)=2^bK/2^c$. Its allowed interval has length at most $2^{b+n}$. Thus $\mathcal H(N\mid P)\le n+c-kL+O(1)$, since $k<\theta$.

For any fixed density $r$ on the free bits, (3) bounds the logarithm of the count by
$$ (L-n-c)H(1/3)+nH(r/3)+(2/3)(n+c-kL)+o(L). $$
Here the coefficient of $c$ is negative. For $r\le 1/2+1/400$, use the tangent at $1/6$ for $H(r/3)$. The resulting exponent is at most
$$ (1-\theta)H(1/3)+\theta\big[H(1/6)+\log_2(5)/1200\big]+(2/3)(\theta-k)<\theta-1/100000. $$
The upper tail of the free words is $O(2^{n H(1/2+1/400)})$, whose exponent is also strictly below $\theta-1/100000$ by $H(1/2+t)\le1-2t^2/\ln2$. This proves the assertion per cylinder. For ordinary intervals in $[X,2X]$ cover the intersection with each of the at most two dyadic ranges by the aligned blocks of $2^n$ integers it meets, using an empty low string. The number of blocks is $O(1+D/X^\theta)$, which gives the interval result.

Run the integer and entropy checks from the repository root with the standard library.

~~~sh
python3 check_common_factor_count.py
~~~
