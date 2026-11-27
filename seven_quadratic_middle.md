# A middle block of binary digits at the three quadratic cofactors

Let $P,R$ be prime powers, $10P+7=DR$, $D\le10000$, and put $j=\lfloor\log_2P\rfloor$. At $N=(2P+1)(5P+1)$ the full conjecture holds outside the following necessary classes for $j$.

| Cofactor $D$ | Period | Classes still possible |
|---:|---:|---|
| 4587 | 138 | $4,26,67,90$ |
| 7383 | 1166 | 683 classes, given by the eight smaller progressions below |
| 9453 | 68 | $15,23,40,64$ |

All other cofactors through 10000 are already settled. Taking the union, **1,036,672 of the 2,735,436 residue classes of binary orders certify the entire family at every size**. The complementary classes retain the [finite certification and sparse exceptional-set bound](seven_ternary_rotation.md).

## The middle scale

The [quadratic cofactor reduction](seven_quadratic_cofactors.md) supplies $(D,z,w)$ and an odd companion coefficient in every residual case. With $T=2P+1$, $X=TR$, $Y=5T$, and $A=(2z,2w,D-2z-2w)$, twice the three parts of that companion are
$$
                      A_1 X,\qquad A_2 X,\qquad A_3 X-Y .
$$
In particular every column of these three binary integers contains at most one 1. From the [upper binary words](seven_quadratic_prime_bases.md) and their [scale error](seven_ternary_rotation.md), there is a listed rational $f$ and an integer $m$ such that
$$
 |X-32f2^m|\le2080\,2^j,\qquad |P^2-c2^m|\le1100000P,\qquad c=8Df/5 .
$$
Write $m=2u+r$, $r\in\{0,1\}$, and $\alpha=\sqrt{c2^r}$. Division by $P+\alpha2^u$ gives $|P-\alpha2^u|\le1100000$, hence $|Y-10\alpha2^u|<12000000$. We may take $P>2^{65536}$ by the [earlier finite computation](seven_quadratic_suffixes.md), so $u>65000$. The twelve possible $\alpha^2$ have denominator below $2^{18}$ and $1<\alpha<150$. Their denominators in lowest terms each have odd 5-adic valuation, hence none is a rational square. In particular $j=u+\lfloor\log_2\alpha\rfloor$: the distance from $\alpha$ to either bounding power of two exceeds $1/(2^{18}\cdot512)$ by factoring the difference of squares.

Let $t=\lfloor\log_2(10\alpha)\rfloor$, and $B=2^{u+t+1}$. The highest bit of $Y$ is at position $u+t$. More precisely the bits of $Y$ down to position $u-K$ are exactly those of $s=\lfloor10\alpha2^K\rfloor$, for $K=200$ or 256. The programs check
$$
 (s+2^{-400})^2<100c\,2^{r+2K}<(s+1-2^{-400})^2 .
$$
The approximation error for $Y/2^{u-K}$ is below $2^{-400}$. Above this middle block, $Y$ has only zeros. Put $h=u+r+4-t$. Then
$$
 \left|\lfloor X/B\rfloor-f2^h\right|<210 ,
$$
since $B>Y>10P$. As $h>100$, the first 32 binary digits of $\lfloor X/B\rfloor$ agree with those of $f$ after the binary point, and its length is $h$. Indeed each listed $f$ is strictly between $1/2$ and 1 and the check measures its distance from the two bounding length-32 dyadic endpoints. The independent checker uses the first 40 digits.

## Testing the shift across both scales

At a cut above $k$ low bits of $X,Y$ the incoming offset is
$$
 s_i=\lfloor A_i x-\delta_{i3}y\rfloor,\qquad
 x=(X\bmod2^k)/2^k,\quad y=(Y\bmod2^k)/2^k .
$$
Use all such offsets for real $0\le x,y<1$. The first program enumerates the breakpoints of each $A_i x$, allowing $\lfloor A_3x\rfloor$ or $\lfloor A_3x\rfloor-1$ in the third position. For input bits $d,e$ of $X,Y$, set $v=s+dA-e(0,0,1)$. A step is allowed if at most one coordinate of $v$ is odd, with upper offset $\lfloor v/2\rfloor$ coordinatewise. At the top, the offset itself must be nonnegative and have disjoint bits.

Read the fixed high prefix of $X$ backward from the top offsets, with zero shift bits. Then iterate allowing both bits of $X$, still with zero shift. This gives a sequence of finite sets $G_i$ for the remaining $i$ free bits above $B$. It repeats from the initial set after the period in the next table; equality of the entire set proves recurrence at all lengths. Independently read the known middle prefix of $Y$ from its low end, starting with every offset and permitting both bits of $X$. Let $V$ be the resulting set at the boundary $B$. An actual row requires $G_{h-32}\cap V\ne\varnothing$. Converting $h$ to $j$ gives the table.

| $D$ | $f$ | $r$ | Period | $j$ residues |
|---:|---|---:|---:|---|
| 4587 | $256/417$ | 0 | 138 | 26 |
| 4587 | $256/417$ | 1 | 138 | 4,67,90 |
| 7383 | $19597/35328$ | 0 | 22 | 5,16 |
| 7383 | $19597/35328$ | 1 | 22 | none |
| 7383 | $486233/876544$ | 0 | 106 | 1,9,14,23,24,57,81,84 |
| 7383 | $486233/876544$ | 1 | 106 | 29,50,62,64,87,102 |
| 7383 | $836139/1507328$ | 0 | 11 | 7 |
| 7383 | $836139/1507328$ | 1 | 11 | 3,4,7 |
| 7383 | $46678369/84148224$ | 0 | 106 | 4,6,11,19,21,45,54,56,58,76,78,98,99 |
| 7383 | $46678369/84148224$ | 1 | 106 | 1,10,15,21,26,40,47,50,52,54,56,57,59,61,70,76,81,84,96,99 |
| 9453 | $256/411$ | 0 | 68 | 15,40 |
| 9453 | $256/411$ | 1 | 68 | 23,64 |

The second checker inverts each of the four allowable output columns, tests each new offset by rational interval intersection, and transports period labels from high to low through the entire prescribed shift prefix. It derives its root bits by successive squaring comparisons. The two programs agree on every residue in all twelve sets.

~~~~sh
python3 check_seven_quadratic_middle.py
python3 check_seven_quadratic_middle_cross.py
~~~~
