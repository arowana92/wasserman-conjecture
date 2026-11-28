# Exact divisors on all integer translates

Fix a prime $p$, a positive integer $k$, and an arbitrary integer $e$. Consider every positive row on
$$
                    (p-1)X+pY-p(p^k-1)Z=e.
$$
Write $e=pv-u$ with $0\le u<p$ and $v\in\mathbb Z$. Let $s$ denote base-$p$ digit sum. The **minimum** of $v_p(K(X,Y,Z))$ over all these rows, at all weights, is independent of $k$. It is

| Condition | Minimum |
|---|---:|
| $v<0$, or $v\ge0$ and $s(v)>u$ | 0 |
| $v=0$ and $u=p-1$ | 2 |
| All other cases | 1 |

For $p=2$ the **entire integer gcd** on $X+2Y-2(2^k-1)Z=e$ is
$$
\begin{cases}
4 & e=-1,\\
2 & e=0\text{ or } e=2^j-1\ (j\ge1),\\
1 & \text{otherwise}.
\end{cases}
$$
Thus two rows of the same weight may independently choose a permutation, $k$, and any offset in the first two cases, and both coefficients then share 2 with the endpoint. For example, with $k=1$ the gcd is four on $X+2Y-2Z=-1$, two on the parallel lines with offsets $0,1,3,7,15,\ldots$, and one on every other translate. These are gcds over all weights on a line. An individual weight can force more primes.

## The lower bound

Put $X=pT+u$ and $N=X+Y+Z=T+p^kZ+v$. When $v\ge0$, the [digit-sum identity](rational_row_line_divisors.md) gives
$$
 v_p(K)=\frac{s(Y)+u-s(v)}{p-1}+v_p\frac{N!}{T!(p^kZ)!v!}.
$$
If $s(v)\le u$, the first term is a positive integer. When $v=0$ and $u=p-1$, $Y=(p^k-1)Z-(p-1)T-u$ is a positive multiple of $p-1$, so $s(Y)\ge p-1$ and the bound is two. It is attained by setting $Y=p-1$, $B=(p^k-1)/(p-1)$, taking $1\le T\le B$ congruent to $-2\bmod B$, and setting $Z=(T+2)/B$. Here $T<p^k$, so $T+p^k Z$ has no carry.

## A carry-free row when the minimum is zero

We use a digit construction. If $W\ge0$, $u<p$, and $s(W)\ge u$, let $C$ consist of the $u$ lowest digit tokens of $W$, with each digit $d_i$ counted as $d_i$ tokens of value $p^i$. Put $T=(C-u)/(p-1)$. Then
$$
             (pT+u)+(W-C)=T+W
$$
is carry-free on each side. To see this, the digit of $T$ at position $i$ counts the chosen tokens strictly above $i$, because $(p^j-1)/(p-1)$ has ones below $j$. Below the highest chosen position all tokens of $W$ were chosen and the sum of the digit and incoming tokens is at most $u$. At and above the highest chosen position $T$ has zero digits. For the left side, below the highest chosen position $W-C$ is zero, at that position the two digit contributions sum to the original digit, and above it $pT+u$ vanishes. The case $u=0$ gives $T=0$. [Lean proof of the digit-allocation lemma](TranslatedLineDigits.lean)

For $v\ge0$ with $s(v)>u$, take $W=v$, $X=pT+u$, $Y=W-C>0$ and temporarily $Z=0$. This gives the specified offset and no carry. For $v<0$, choose $J$ with $p^{J-1}>|v|$, take $Z=p^J$, and put $W=p^{J+k}+v$. Its digits at $J-1$ and $J$ are both $p-1$, so the $u$ lowest tokens are strictly below $J$. Set $X=pT+u$ and $Y=W-C-Z$. The digit at $J$ can supply $Z$ without borrowing. There is at least one further unused token if $J$ is taken sufficiently large. Thus $Y>0$ and the addition of the three parts is carry-free.

In each construction one can append a high homogeneous block with parts $(p(p^k-1),0,p-1)$ at any digit position above the entire weight. This makes $X,Z$ positive if needed and gives infinitely many weights with the same valuation. The block's first and third parts occupy disjoint digit positions.

## A row with exactly one carry

In the remaining cases let $m=u-s(v)$, so $0\le m\le p-2$. Take $Z=p^J>v$ and $H=J+k-1$. Define
$$
 (p-1)T=v-s(v)+m(p^H-1)+(p^H-p^J),
 \qquad Y=(p-1-m)p^H,\qquad X=pT+u.
$$
Then $Y>0$, and substitution gives $N=T+p^kZ+v=X+Y+Z$. To check that $T+p^k Z+v$ has no carries, the digits of $(v-s(v))/(p-1)$ count the tokens of $v$ at strictly higher positions. The other two summands in $T$ add $m$ at each position below $H$, and one at positions $J,\ldots,H-1$. Below $J$ the digits of $T+v$ are at most $m+s(v)=u$. From $J$ through $H-1$ they are at most $m+1\le p-1$. All digits from $H$ upward in $T+v$ vanish, and $p^kZ=p^{H+1}$. Since $s(Y)=p-1-m$, the identity gives valuation one. Append a high homogeneous block as above if needed.

There are exactly $\binom{a+p}{p-1}$ offsets forcing $p$ among the $p^{a+1}$ consecutive offsets $1-p\le e\le p^{a+1}-p$ for $a\ge1$. For fixed $u$, these correspond to $a$ digits of $v$ with sum at most $u$. Their number is $\binom{a+u}{u}$ since $u\le p-1$. Summing over $u$ proves the count.

## Why no odd prime divides the gcd in the binary case

Write $t=2(2^k-1)$ and fix any odd prime $q$. Take $Q=q^J>|e|$ with $J\ge1$. Let $A=t-\mathbf1_{e<0}$, take $b=1$ if $A\equiv-1\pmod q$ and $b=0$ otherwise, and set
$$
          Z=Q,\qquad Y=bQ,\qquad X=(t-2b)Q+e.
$$
If $b=0$ the digit of $X$ at $J$ is $A\bmod q\le q-2$, so it can accept the digit of $Z$. If $b=1$ that digit is $q-3$, leaving room for both $Y$ and $Z$. All their other digits vanish. The only case with $X=0$ is $(k,e,q)=(1,0,3)$, where $(4,1,3)$ works directly. In particular $X$ cannot be negative when $t=2$ and $e<0$, since then $A=1$ and $b=0$.

To supply a missing $Y$, take a multiple of the edge row $(0,2^k-1,1)$ of weight $q^f-1$, using $2^k\mid q^f-1$. Its two nonzero parts have complementary digits. Append it above the previously constructed row. This gives a positive row at the required offset with valuation zero at $q$.

Run `python3 check_translated_line_minima.py` for 1,990 comparisons with an independent shortest-path computation, 988 direct checks using exact binomial coefficients, 2,352 checks at large offsets, and 7,260 witnesses avoiding odd primes. The shortest-path computation starts its relation carry at $-e$ for the primitive normal. Its interval of possible carries is the convex hull of $-e$ and the [homogeneous bounds](rational_row_line_divisors.md). Run `python3 check_translated_line_minima.py --line 2 1 -1` to print a positive row attaining valuation two.
