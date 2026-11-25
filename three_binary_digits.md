# Odd weights with three binary ones

**Theorem.** For every sufficiently large odd $N$ with exactly three ones in its binary expansion, any three positive trinomial coefficients of weight $N$ have a common prime divisor. The threshold is uniform in both exponents of $N=1+2^s+2^t$, $0<s<t$. Odd weights with at most two binary ones have the universal divisor 2 at every weight $N\ge3$.

## The unique odd coefficient

Suppose there is a counterexample at $N=1+2^s+2^t$. Its row avoiding 2 must be $A=(1,2^s,2^t)$, since the three binary units must go to different positive parts. The rows $B,C$ must between them avoid every prime of $N$. Assign each such prime to one of them which avoids it. The product of the full powers assigned to each divides all parts of the respective row. Both products are at least 3, since the parts in the complementary row have to sum to $N$. Their product is $N$, so one of $B,C$ has minimum at least $\sqrt N$. Neither is the endpoint $(1,1,N-2)$.

By the height 6 bound, the second part of every sorted row is at least $120^{-7/5}N^{1/5}$ for $N\ge12$. Thus, for sufficiently large $N$, both $B,C$ have second part at least 150. Applied to $A$ the same bound excludes $s\le\sqrt t$ for large $t$. We may also suppose $t-s\ge4$. The three fixed differences 1, 2 and 3 give $N=(1+2^{t-s})2^s+1$, already covered for sufficiently large $s$ by the fixed-multiplier theorem. [Height 6 bound](power_intervals.md) · [Fixed-multiplier theorem](binary_weights.md)

## A short relation among the exponents

Set $L=\lceil\sqrt t\rceil$. Place the $L+1$ fractional parts of $i s/t$, $0\le i\le L$, in $L$ equal intervals. Two in the same interval give integers $1\le b\le L$, $0\le a\le b$, and $|bs-at|\le t/L\le\sqrt t$. Since $s>\sqrt t$, we have $a\ge1$. Divide $a,b$ by their gcd and put $D=bs-at$ for the reduced pair. We retain $1\le a\le b\le L$ and $|D|\le\sqrt t$. The exact relation is
$$2^{\max(-D,0)}(2^s)^b=2^{\max(D,0)}(2^t)^a.$$

Fix $H=300$. For integers $u,v\ge0$, define
$$R(u,v)=2^{\max(-D,0)}u^b-2^{\max(D,0)}v^a.$$
Omit an offset $0\le h<H$ if it can be written $h=1+u+v$ with $R(u,v)=0$. We claim that **at most 18 offsets** are omitted.

If $a=b=1$, then $v=2^{t-s}u$ at a zero, and the values of $h$ lie in the progression $1+(1+2^{t-s})u$ of step at least 17. This gives at most 18 offsets. Otherwise $b\ge2$. There is the pair $(0,0)$. At any positive solution and every odd prime $q$, we have $b v_q(u)=a v_q(v)$, so the valuations are $(a n_q,b n_q)$ for some $n_q\ge0$. For $q=2$ the valuations solve $b i-a j=D$. The nonnegative solutions form one progression $(i_0+a n,j_0+b n)$ for $n\ge0$, since $a,b$ are coprime. Thus every positive pair has the form
$$u=2^{i_0}z^a,\qquad v=2^{j_0}z^b\qquad(z\text{ a positive integer}).$$
There are at most $\lfloor\sqrt{300}\rfloor=17$ possibilities when $v<H$, proving the claim.

## Bounding the contribution of the binary row

Let $J$ be the offsets retained. Put $K$ equal to $\operatorname{lcm}(1,\ldots,H)$ times the product of $|R(u,v)|$ over the nonnegative pairs with $1+u+v\in J$. All factors are nonzero. There are fewer than $H^2$ of them, and each is bounded by $2^{|D|+1}H^b$. In particular,
$$\log K=O_H(\sqrt t)=o(\log N)$$
uniformly over the exponent pair. Since $K$ may depend on $N$, this uniform estimate matters.

If $w=p^e>H$ divides $N-h$ and $A$ avoids $p$, its residues sum to $h$. The residue of the first part is 1, so the other residues $u,v$ satisfy $h=1+u+v$. The exact power relation shows $w\mid R(u,v)$. For $h\in J$ this implies $w\mid K$. In particular $h=0$ cannot be covered by $A$ at such a layer.

Set $L'=\sum_{h\in J}(H-h)\ge W(300)-18\cdot300=39{,}750$, with $W(x)=x(x+1)/2$. Remove $\gcd(N-h,K)$ from each $N-h$. At every remaining prime-power layer only one $h$ occurs and one of $B,C$ must avoid the prime. The weighted lemma therefore gives
$$\prod_{h\in J}\left(\frac{N-h}{\gcd(N-h,K)}\right)^{H-h}\ \mid\ P,$$
where $P$ is the product for cutoffs of total 300 on each of $B,C$. On the row of minimum at least $\sqrt N$ use $(100,100,100)$, and on the other use $(0,150,150)$ in increasing order of the parts. For sufficiently large $N$ these fit. Hence
$$P\le N^{3W(100)+2W(150)}=N^{37{,}800}.$$
The left side is at least $(N/(2K))^{L'}=N^{L'-o(1)}$ uniformly as $t\to\infty$, a contradiction. This proves the theorem. [Weighted divisibility lemma](weighted_obstruction.md)

For an arithmetic check of the short-relation step, run `python3 check_three_binary_digits.py`. For eight exponent pairs the script enumerates all 358,800 pairs of candidate residues with $1+u+v<300$. It checks that $\gcd(2^s-u,2^t-v)$ divides the proposed relation and counts the actual omitted offsets. The pair $(s,t)=(12,16)$ attains the upper bound of 18.
