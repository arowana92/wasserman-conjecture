# The leading digits at cofactor five

**For every odd prime power $Q$ and every $1\le k\le45$ except $32,34,44$, the full three-row conjecture holds at $N=1+5\cdot2^kQ$.** In particular this settles the previously open family $k=31$ and all ten new binary exponents except for one class at 44. At the three remaining exponents a counterexample requires respectively
$$
 Q\equiv 5727322795\pmod {2^{33}},\quad
 Q\equiv 429493347\pmod {2^{35}},\quad
 Q\equiv 26828080362291\pmod {2^{45}} .
$$
All assertions here include arbitrary exponents of the odd prime. The [earlier finite-range and exponent-separation results](five_fixed_binary_separation.md) apply to 32 and 34.

There is also a restriction at **every** binary exponent. Put $P=2^k$, and let $c$ be the least positive inverse of $Q$ modulo $2P$. A possible failure requires the reduced dyadic fraction $r=c/(5P)$ to belong to 11,742 explicitly computable intervals in $[0,2/5]$ of combined length **less than 0.015382**. In particular
$$
             \min(r,2/5-r)<\frac{95574}{2286935}<0.0418 .
$$
The first two digits in bases 2 and 3 already require $r\le1/15$ or $r\ge1/3$. The interval restrictions follow solely from leading carries once the row and its arithmetic progression have been determined.

## A test on any three part proportions

Let the positive parts be $x_i$, with sum $N$. If their coefficient avoids $p$, all base-$p$ additions are carry-free. At the scale $M=p^j$ write $v_i=\lfloor x_i/M\rfloor$ and $n=\sum v_i=\lfloor N/M\rfloor$. Then
$$
                        \frac{v_i}{n+1}<\frac{x_i}{N}
$$
for all $i$, since $N<M(n+1)$ and the parts are positive. Starting at $(0,0,0)$, extend the leading words by $v_i\mapsto pv_i+b_i$, using nonnegative digits of sum less than $p$. Require a nonzero sum at the first step. If the proportions have upper bounds $u_i$, discard every prefix with $v_i\ge(n+1)u_i$ in some position. If no path has length $h$, the coefficient is divisible by $p$. To apply the argument at more steps than the weight has digits, multiply every part by a sufficiently large power of $p$. This preserves both the proportions and carry freedom.

At cofactor five, a possible counterexample has a row
$$
 (3P Q,\ cQ,\ (2P-c)Q+1).
$$
Its coefficient must avoid 2, 3 and 5, and $15\mid c$. The preceding conditions follow from the all-weight endpoint reduction and the [cofactor-five row theorem](factor_five_binary_exponents.md#parametrizing-the-row-by-the-low-bits-of-the-odd-power). At a surviving prefix $Q\equiv q\bmod 2P$ with $0<q<2P$, the three proportions are bounded above by $(3P,c,2P-c+1/q)/(5P)$. These caps depend only on the template.

At $k=31$ the ternary obstruction has a stronger standalone form. Put
$$
       (A,B,C)=(6442450944,4192014075,102953221).
$$
For **every positive integer $w$**, and $t=0$ or 1, the coefficient of $(Aw,Bw,Cw+t)$ is divisible by 3. For $w\ge65536$ use the bounds $(A,B+1/65536,C+1/65536)/(A+B+C)$. No ternary path has length 28. For $1\le w<65536$ the first positive factorial layer is always among the first 12, for both offsets. [Lean proof](FiveLeadingDigits.lean) · [Integer verifier](check_five_leading.py)

## Extending the check through 45

The previous [low-digit enumeration](factor_five_binary_exponents.md) leaves one template each at 31, 32 and 34 through exponent 35. The ten additional exponents give 14,206,342 templates.

| $k$ | Low-digit templates | After the low-digit graphs | After the leading tests and the short check at 44 |
|---:|---:|---:|---:|
| 36 | 117907 | 0 | 0 |
| 37 | 178495 | 0 | 0 |
| 38 | 269766 | 0 | 0 |
| 39 | 407955 | 4 | 0 |
| 40 | 617815 | 1 | 0 |
| 41 | 934867 | 0 | 0 |
| 42 | 1413575 | 4 | 0 |
| 43 | 2139305 | 6 | 0 |
| 44 | 3234656 | 8 | 1 |
| 45 | 4892001 | 9 | 0 |

The independent checker lifts the inverse $c$ together with $q$, testing the factorial layer for $3q+2q$ at every binary prefix. For each retained inverse it consumes the original certificate in order, verifies the exact template, and exhausts its claimed graph using the high quotients of the total and the first two parts. At 2 it also tracks $3Q,2Q$ directly. It compares all counts and checks for extra records at the end. The subsequent leading calculation rules out 30 of the 32 remaining shapes from these ten exponents, in at most 81 columns. For example, two binary digits already exclude $c=645587271675$ at 39, and two ternary digits exclude $c=790417059915$ at 40.

At 44 there is a further exclusion for $c=8136949755$, $q=28147172193075$. Above or at $Q=3^{58}$ use the sharper last cap with $1/3^{58}$ in place of $1/q$. No ternary prefix survives 58 steps. Below the cutoff enumerate the binary extensions of the specified 45 bits through the 92-bit ceiling. Reject overlaps of the actual parts and of $3Q,2Q$ as each bit is read. The seven resulting prefixes all fail the full binary addition when the upper multiplication digits are included. An independent leading check grows the total first and enumerates the allowed parts by their individual caps. It confirms all 31 leading exclusions and repeats the short check with factorial quotients.

## Projecting the leading cells onto the fraction $r$

Write $N=5t+1$ and $x_1=3t$, where $t=PQ$. At each scale $M=p^j\ge2$ put $a=v_1$, $b=v_2$, $z=v_3$ and $n=a+b+z$. For a row avoiding $p$ the remainders $R_i=x_i-Mv_i$ sum to $N\bmod M<M$. Hence
$$
    M(5a-3n)=-2R_1+3(R_2+R_3)-3,\qquad -2\le 5a-3n\le2 .
$$
Putting $d=5a-3n$, a further leading digit triple $(i,j,l)$ takes $(n,d,b)$ to
$$
      (pn+i+j+l,\;pd+2i-3(j+l),\;pb+j).
$$
At the last chosen column the interval supplied by $(n,d,b)$ is
$$
    \left[\frac{b}{n+1},\ \frac25-\frac{z}{n+1}\right]
    =\left[\frac{b}{n+1},\ \frac{5b+2+d}{5(n+1)}\right].                 \tag{1}
$$
Indeed $r=x_2/(N-1)>b/(n+1)$, and $r<2/5-z/(n+1)+1/(N-1)$ by applying the same prefix inequality to $x_3$. Both $r$ and $2/5$ have denominator dividing $5P$. If $Q\ge n+1$, a strictly positive excess over the upper endpoint is at least $1/(5P(n+1))\ge1/(N-1)$, a contradiction. For the computations below $n+1\le2^{20}$. All $Q<2^{20}$ satisfy the full conjecture at every $k$ by the [bounded odd-power theorem](factor_five_exponents.md). Larger $Q$ give sufficient digits at all three bases with $M\ge2$.

Intersect the unions in (1) for 20 binary places, 12 ternary places and eight quinary places. There are 12,978 components, including 1,236 singletons. Every endpoint has reduced denominator at most $5\cdot2^{20}$. In a possible exception $k\ge32$ by the result above, and $r=(c/5)/2^k$ has reduced denominator exactly $2^k$ since $c$ is odd. All endpoints are therefore excluded. The union before removing endpoints is symmetric under $r\mapsto 2/5-r$ by interchanging the second and third digit words. On the lower half of the interval the largest right endpoint of a positive-length component is $95574/2286935$. The total length in the theorem is less than 3.846 percent of the full interval of length $2/5$.

One verifier maps rational endpoints down to a grid of mesh $2^{-52}$, which preserves their exact order since distinct endpoints are separated by at least $(5\cdot2^{20})^{-2}$. An independent verifier enumerates the integer $n$ first. Its only possible $a$ is $\lfloor(3n+2)/5\rfloor$. It checks the shorter prefixes and then distributes the available digits of $n-a$ between the two other parts. It sorts and intersects intervals by exact cross multiplication and gives a combined length at most $15381646326/10^{12}$. At depth two alone, the binary intervals miss $(3/20,1/4)$, and the ternary intervals miss $(1/15,1/6)$ and $(7/30,1/3)$. Their overlap gives the stated shorter test. [Exact interval checker](check_five_intervals_cross.cpp)

The [binary direction theorem](binary_direction_measure.md) shows that refinements can make the total permitted length arbitrarily small.

## Reproduction

Run in a shell with pipeline failure propagation. The large enumeration streams its certificates between the two programs and retains only the residual shapes.

~~~sh
set -o pipefail
mkdir -p .work
c++ -O3 -std=c++17 check_factor_five_binary_exponents.cpp -o .work/check_factor_five_binary_exponents
c++ -O3 -std=c++17 check_five_leading_cross.cpp -o .work/check_five_leading_cross
.work/check_factor_five_binary_exponents 36 45 --witnesses | .work/check_five_leading_cross 36 45 > .work/five_leading_rows.txt
python3 check_five_leading.py .work/five_leading_rows.txt
python3 check_five_prefixes_cross.py .work/five_leading_rows.txt
c++ -O3 -std=c++17 check_five_intervals_cross.cpp -o .work/check_five_intervals_cross
.work/check_five_intervals_cross
lean FiveLeadingDigits.lean
~~~
