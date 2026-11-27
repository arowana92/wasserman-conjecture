# Two nearly equal parts of an odd trinomial coefficient

Let a positive trinomial row have weight $N$, odd coefficient, and parts $u>v$ and $w$. Put $e=u-v$. Let $T$ be the largest power of two at most $u$. Then there are integers $j\ge0$, $c,r\ge0$, and $d\ge1$ with
$$
 u=T+c,\qquad v=T-d,\qquad w=2jT+r,\qquad c+r\le d-1 .
$$
In particular $1\le 2(j+1)T-N\le e$. The interval $(N,N+e]$ must contain a multiple of $2T$, which is larger than either of the selected parts.

By Kummer's criterion the three bit sets are disjoint. Thus $v<T$ and $w$ omits the bit $T$. Below $T$, the bits available after assigning $v=T-d$ are exactly those of $d-1$. The disjoint integers $c$ and $r=w\bmod 2T$ use only these bits. This proves the representation and the bound since $e=c+d$.

## Which limiting proportions occur

For fixed positive real numbers $a,b$, odd coefficients can approach the proportions $(a,a,b)$ with coordinate errors $o(W)$ at scale $W\to\infty$ **if and only if $b/a$ is an even integer**. In that case bounded errors already suffice.

More explicitly, suppose $u,v,w$ lie within $E\ge1$ of $aW,aW,bW$ respectively, with $u\ge v$. At an odd coefficient the two positive parts are distinct. Then $|T-aW|\le E$, and $r\le d-1<2E$. Writing $m=2j$, we obtain
$$
 |b-ma|W\le (m+3)E .
$$
If $aW\ge 2E$, then $m\le(bW+E)/(aW-E)\le2b/a+1$. When $b/a$ is not an even integer put $\delta=\operatorname{dist}(b,2a\mathbb Z)>0$. Every such row with $aW\ge2E$ and $W>(2b/a+4)E/\delta$ has even coefficient. For integer slopes the single condition $W>(2b/a+4)E$ suffices, since $\delta\ge1$ and $a\ge1$. Conversely, if $b=2ja$ with $j\ge1$, take $T=2^h$, $W$ a nearest integer to $T/a$, and the three parts $T,T-1,2jT$. Their bit sets are disjoint, and the errors are bounded independently of $h$.

This gives a single common prime for any three rows that each meet the sufficient condition, allowing different slopes and errors in each row. In particular it applies when two rows meet the condition and the third is the endpoint $(1,1,N-2)$, whose coefficient is even.

## The exact least gap at a fixed weight

Suppose $N$ has at least three binary ones, as is necessary and sufficient for an odd positive trinomial coefficient. Fix an occupied bit $T$ greater than the lowest bit $T_0$ of $N$. Among all odd coefficients and pairs of parts for which $u$ has highest bit $T$ and $v<u$, the smallest possible gap is
$$
 T-(N\bmod T)+
 \begin{cases}
 T_0,&T\text{ is the highest bit of }N,\\
 0,&\text{otherwise}.
 \end{cases}
$$
Indeed put $R=N\bmod T=c+v+r$. Then $u-v=T-R+2c+r$. If $T$ is below the highest bit, $u=T,v=R$ leaves all higher bits for $w>0$ and attains the bound. At the highest bit $w=r>0$ needs at least $T_0$, and $u=T,v=R-T_0,w=T_0$ attains the bound. There is no eligible pair at the lowest bit. For consecutive occupied bits $T<T'$, the change in $T-(N\bmod T)$ is $T'-2T\ge0$. Thus the overall minimum is exactly the difference of the two lowest occupied bits.

As a counting consequence, fix $0<\eta\le1$ and $1\le H\le X$. There are $O(H/\eta)$ weights $N\in[X,2X]$ admitting an odd row with a pair $u>v$ satisfying $u\ge\eta X$ and $u-v\le H$. In fact the next $H$ integers after $N$ contain a multiple of $2T>\eta X$, hence of the least power of two exceeding $\eta X$. At most $1+(X+H)/(\eta X)$ such multiples are relevant and each accounts for at most $H$ weights.

Run **python3 check_binary_near_equal_parts.py** for an independent enumeration using exact factorial coefficients. It checks the gap for every eligible bit at each tested weight, the interval restriction for every odd row, and examples on the unrestricted slopes. [Verifier](check_binary_near_equal_parts.py)
