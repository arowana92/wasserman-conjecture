# Prime powers forced on rational lines of trinomial rows

Write $K(X,Y,Z)=N!/(X!Y!Z!)$ at $N=X+Y+Z$, with all parts positive. For any prime $p$ and positive exponents $h,k$, the condition
$$
             N=\lfloor X/p^h\rfloor+p^k Z                                      \tag{1}
$$
forces $p^{\gcd(h,k)}\mid K(X,Y,Z)$. This includes a family of $p^h$ parallel lines. On each line the minimum valuation in **every** weight residue class has a closed formula, and infinitely many rows attain it.

For example, $X+2Y-6Z\in\{-1,0\}$ forces parity. On the line $X+2Y=6Z$ the best guaranteed powers of two are $8$ at $N\equiv0\pmod 7$, $4$ at residues $1,2,4$, and $2$ at residues $3,5,6$. Two positive rows at the same weight can satisfy (1) with different orderings and different exponents at $p=2$. Their coefficients then share 2 with the endpoint coefficient $N(N-1)$. The same argument works at any prime $p\mid N(N-1)$. See the [endpoint reduction](all_weights_endpoint_reduction.md).

## The valuation on each parallel line

Put $d=\gcd(h,k)$, $L=h+k$, $D=p^d-1$, $M=p^L-1$, and let $s(n)$ denote the sum of the base-$p$ digits of $n$. Write $X=p^hT+u$ with $0\le u<p^h$. Condition (1) says
$$
 (p^h-1)X+p^hY-p^h(p^k-1)Z=-u,\qquad N=T+p^k Z .
$$
It gives $Y+u=(p^k-1)Z-(p^h-1)T=MZ-(p^h-1)N$. The digit-sum formula for factorial valuations therefore gives the identity
$$
 v_p(K)=\frac{s(Y)+s(u)}{p-1}+v_p\binom{T+p^k Z}{T}.                      \tag{2}
$$

Fix $u$ and a weight class $N\bmod (M/D)$. Let $R$ be the representative in $[1,M]$ of $-u-(p^h-1)N\bmod M$. The exact minimum of (2) on the class is
$$
                        \frac{s(R)+s(u)}{p-1}.                              \tag{3}
$$
To prove the lower bound, fold the base-$p$ digits of the positive integer $Y$ into $L$ bins, sending position $i$ to $i\bmod L$. Move each group of $p$ tokens to a single token in the next bin cyclically. Each move reduces the count by $p-1$ and preserves the residue modulo $M$. The final nonempty array represents an integer from 1 to $M$, necessarily $R$. Thus $s(Y)\ge s(R)$. Also $D\mid Y+u$ and $Y+u>0$. Folding $Y+u$ into $d$ bins gives $s(Y)+s(u)\ge s(Y+u)\ge d(p-1)$.

For sharpness put $Y=R$, $A=(p^h-1)/D$, $B=(p^k-1)/D$. The numbers $A,B$ are coprime, and $D\mid R+u$. Choose $1\le T\le B$ such that $A T\equiv-(R+u)/D\pmod B$, and set $Z=((R+u)/D+AT)/B$. This gives positive parts and the required weight class since $\gcd(p^h-1,M)=D$. As $T<p^k$, the binomial term in (2) has no carries. To retain the same $Y,u$ at arbitrarily large weights, add $(p^k-1)p^j$ to $T$ and $(p^h-1)p^j$ to $Z$ with $p^j>N$. The two new blocks in $T+p^k Z$ have disjoint digit positions. The weight increases by $Mp^j$ and (3) is attained again.

On the homogeneous line $u=0$, the residues $R$ are precisely $D,2D,\ldots,M$. The global minimum is $d$, attained at $R=D$. The maximum class minimum is $L$, attained only at $R=M$. Exactly
$$
                       \binom{L/d+p-2}{p-1}^{d}
$$
weight classes have minimum $d$. Indeed, equality in the $d$-bin fold means that the sum of the $L/d$ digits in each bin was exactly $p-1$ before any cyclic carrying. Stars and bars counts these independent choices. In particular, for $p=2$, $h=k$, all $2^h$ nonzero classes modulo $2^h+1$ give minimum $h$, while the zero class gives $2h$.

More generally, put $v=N-T-p^k Z$ and suppose $v\ge0$. The same calculation gives
$$
 v_p(K)=\frac{s(Y)+s(u)-s(v)}{p-1}
       +v_p\frac{N!}{T!(p^kZ)!v!}.
$$
In particular $s(v)\le s(u)$ forces $p\mid K$, since $Y>0$ and $s(Y)+s(u)-s(v)$ is a multiple of $p-1$. For each fixed $h,k$ this supplies infinitely many other parallel lines by taking $u=1$, $v=p^j$.

## The exact gcd on any homogeneous line

For arbitrary positive $r,s,t$, consider all positive rows on $rX+sY=tZ$ over all weights. The gcd of their coefficients is effectively computable, including prime powers and a minimizing row for each candidate prime. Normalize by $\gcd(r,s,t)$ and exchange $r,s$ if needed so $r\le s$.

Only primes dividing $(r+t)/\gcd(r,t)$ or $(s+t)/\gcd(s,t)$ need checking. Indeed, on the edge $Y=0$ the primitive pair has weight $m=(r+t)/\gcd(r,t)$. If $p\nmid m$, choose $e$ with $m\mid p^e-1$ and scale the pair to sum to $p^e-1$. Its two parts have complementary base-$p$ digits. The other edge admits the same construction when its primitive weight is a unit. Concatenating the two blocks at disjoint digit positions gives a positive row with no carries.

For a candidate $p$, keep the relation carry $c$, addition carry $a$, and set of occupied parts. Begin with $(c,a)=(0,0)$ and the empty set. Read digits $i,j,l$ from least to most significant and allow the transition when
$$
 c'=(c+ri+sj-tl)/p
$$
is integral. Set $a'=\lfloor(a+i+j+l)/p\rfloor$, mark the nonzero coordinates, and give the edge cost $a'$. The relation carry stays in $[-t,r+s]$, and $a$ stays in $\{0,1,2\}$. An accepting state has both carries zero and every part marked. Multiplying the relation recurrence by the successive powers of $p$ shows that each accepting path gives a positive row on the line. Every such row gives an accepting path after appending enough high zeros. The cost is its valuation by the multinomial carry criterion, counting a carry of two twice. A shortest path therefore gives the exact minimum. At least one positive row always exists, for example $(t,t,r+s)$.

## Which lines force a prime?

For a primitive normal with $p\mid b=\gcd(r+t,s+t)$, every coefficient is divisible by $p$ by the [primitive weight divisor](row_line_obstructions.md). The homogeneous lines in (1) have $b=1$ and give a different source of divisibility. The following classification has been checked in finite boxes. The height is $\max(r,s,t)$, with $1\le r\le s$, $t>0$.

| Height | Primes checked | Primitive lines | Forced lines outside $p\mid b$ |
|---:|---|---:|---|
| 256 | 2 | 7,002,070 | 29, all from (1) |
| 64 | every prime | 110,646 | all from (1) |

For the second box, primes above 128 are excluded by the two-edge argument. At odd primes the additional lines number four at 3, one at 5, one at 7, and zero at every larger prime. An unrestricted classification of homogeneous lines would require an argument beyond the census.

The [census](check_line_prime_census.cpp) restricts to digits with $i+j+l<p$. Here $c\in[-t,s]$. Once the two-edge test and $p\mid b$ have been handled, exactly one edge is blocked. Indeed, if both reduced edge weights are multiples of $p$, then $p\mid r+t,s+t$. The good edge already supplies two positive parts. It suffices to find a closed carry-free walk through zero using the missing coordinate at least once, and append the good edge. This tests for a positive row at any weight. The program checks the resulting decision against (1).

Run the following commands from the repository root.

```sh
python3 check_line_divisors.py 3 4 12
python3 check_line_divisors_cross.py
c++ -O3 -std=c++17 check_line_prime_census.cpp -o /tmp/check_line_prime_census
/tmp/check_line_prime_census 2 256
for p in 2 3 5 7 11 13 17 19 23 29 31 37 41 43 47 53 59 61 67 71 73 79 83 89 97 101 103 107 109 113 127; do
    /tmp/check_line_prime_census "$p" 64
done
```

The [independent checks](check_line_divisors_cross.py) attain 3,099 homogeneous class minima, test (2) on 25,458 rows, and attain 13,140 offset/class minima twice with different weights. A separate full three-coordinate search makes 23,832 carry-free decisions without using the two-edge reduction. For all 171 primitive lines through height seven, finite lists of exact multinomials attain the gcd computed by the shortest-path algorithm, using 250,009 coefficient evaluations.
