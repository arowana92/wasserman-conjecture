# Two predecessor powers and the digits above one of them

Write $N=1+DPR$, where $D\ge1$, $P=p^k$, $R=q^e$, the primes $p,q$ are distinct and $k,e\ge1$. The powers need not be full powers of their bases in $N-1$.

**One low digit.** If $DR\equiv1\pmod L$ for a power $L=p^j\ge D$, $j\ge1$, the full trinomial conjecture holds at $N$. This condition is independent of $k$. For example, with any fixed $p\nmid D$, choose such an $L$. Every prime $q\equiv D^{-1}\bmod L$ works for every $k\ge1$, as does every power of $q$ congruent to $q\bmod L$. There are infinitely many such primes by Dirichlet's theorem.

**Two low digit units.** Suppose $h=[DR]_L$ is the sum of exactly two base-$p$ digit units and $Dh\le L$. Write their sizes as $w,rw$ with $r$ a power of $p$, allowing $r=1$. The full conjecture holds if
$$
                         D\le(1+r)(4(1+r)+1).
$$
It also holds at arbitrary $D$ if the same two-unit hypothesis holds in the other direction, with $DP$ modulo a power of $q$. For instance, the one-sided theorem proves every weight
$$
                    1+385\cdot 2^k\cdot41^{2+1024m},
                    \qquad k\ge1,\quad m\ge0.
$$
Indeed $41^{1024}\equiv1\bmod8192$ and $385\cdot41^2\equiv17\bmod8192$. Here $385\cdot17<8192$, the ratio is 16, and $385<17\cdot69$. There are five distinct predecessor primes.

## The carry-free cuts

In a counterexample the endpoint reduction gives the row $(1,1,N-2)$. Assign each prime of $N(N-1)$ to one of the other two rows avoiding it. The products of the assigned full powers of $N$ are $U,V\ge3$ with $UV=N$. The parts in the respective rows are multiples of $U,V$. At each assigned predecessor power the three residues are $0,0,1$. Above the $P$ shift, the three quotients on division by $P$ in the row taking $p$ split the digits of $DR$ without a base-$p$ carry. These statements follow from Kummer's criterion. [Endpoint reduction](all_weights_endpoint_reduction.md) · [Two-power assignment](two_prime_powers.md)

At any layer $L=p^j$ with $1\le h=[DR]_L<L$ and $Dh\le L$, consider both assignments.

* If one row takes both primes, the zero positions modulo $P,R$ intersect. One part is $aPR$ with $U\mid a$, $0<a<D$. Above the shift, $aR,(D-a)R$ add without a carry. Their residues $v,h-v$ modulo $L$ are both positive, since $R$ is a unit and $a,D-a<L$. As $Dv\equiv ah\bmod L$ and both sides lie strictly between 0 and $L$, we have
  $$
                              Dv=ah.
  $$
  In particular $U\mid v$ because $(D,U)=1$. Also $a=(D/g)c$ with $g=(D,h)$, $0<c<g$ and $U\mid c$. Grouping is impossible whenever $g\le3$. A single digit unit cannot give two positive residues. At two digit units, $v$ is a pure $p$-power, also impossible since $U>1$ is prime to $p$.

* If the primes go to different rows, take two parts $UPx,UPy$ in the row taking $p$ and $VRz,VRw$ in the row taking $q$. The four multipliers are positive. The third normalized parts are positive, so $V>P(x+y)$ and $U>R(z+w)$. Multiplication gives $(x+y)(z+w)<D$, and the size of the two displayed parts in the second row gives $2V<DP$. Likewise $2U<DR$. The residues $i=[Ux]_L,j=[Uy]_L$ are positive because $U$ is a unit and $x,y<D/2\le L$. They use disjoint digit units of $h$. They satisfy $xj\equiv yi\bmod L$ with both products less than $L$, so
  $$
                              xj=yi.
  $$
  One digit unit is impossible. At two units we may order the multipliers as $x,rx$.

For the two-unit estimate, $D\mid xzw$. In fact each full prime power of $D$ at a base other than $p,q$ divides a displayed multiplier of its assigned row, by intersecting the zero positions. If it divides $rx$ it divides $x$. Any $p$-part of $D$ divides both multipliers in the first row by using the full predecessor power above $P$. Any $q$-part divides both multipliers in the second row. Thus
$$
            (1+r)x(z+w)<D\le xzw.
$$
Since $4zw\le(z+w)^2$, this requires $z+w\ge4(1+r)+1$, proving the claimed bound. With two sparse layers, the other two multipliers are $z,r'z$, where $r'$ is a $q$-power. Now $D\mid xz$, contradicting $(1+r)(1+r')xz<D$. The two ratios may coincide. Related results use low digits at primes of the weight. [Weight-prime argument and formal ratio lemma](two_digit_weight_cuts.md)

## A polynomial test independent of the number of digits

Use the canonical radix-$P$ expansion
$$
              N=F(P),\qquad F(X)=1+hX+f_2X^2+\cdots+f_nX^n,
              \qquad K=\max_{0\le l\le n} f_l .
$$
Suppose $h>0$, $P>D$, $Dh\le P$ and $D(h+2)K\le4P$. A split assignment requires a factorization
$$
 F(X)=(1+tX)G(X),\qquad 2\le t\le h,\qquad G\in\mathbb Z_{\ge0}[X].
 \tag{1}
$$
In particular $t\mid f_n$. **If grouping is excluded and no divisor $2\le t\le h$ of $f_n$ gives (1), the full conjecture holds.** In particular a leading digit $f_n=1$ excludes splitting. Either $(D,h)\le3$ or at most two base-$p$ digit units in $h$ suffice to exclude grouping. The bound on the digits is independent of the degree. It can therefore test an arbitrarily long expansion without factoring $N$. For example, take $D=55$, $P=1024$ and the prime $R=20{,}610{,}923{,}369{,}137$. At $N=1{,}160{,}807{,}204{,}149{,}795{,}841$ the digits from lowest to highest are $1,7,1,2,4,7,1$. The low digit test in either direction has more than two units before the required gap, but the polynomial has leading digit one and $(55,7)=1$.

To prove this, apply the split argument at $L=P$ and select the smaller residue $i\le h/2$, retaining its multiplier $x$. Since $UV=N\equiv1\bmod P$,
$$
                  iV=bP+x,\qquad 0<b<iD/2\le Dh/4.
$$
Indeed $V>2P>x$ and $2V<DP$. The radix-$P$ digit polynomial $B(X)$ of $Ux$ has $B(0)=i$, and each digit is at most the corresponding digit of $DR$ by carry freedom. In particular its coefficients are at most $K$. We have
$$
                      (bP+x)B(P)=ixF(P).
$$
Every coefficient on the left is at most $(b+x)K<P$, and on the right it is at most $ixK<P$, by $x<D/2$, $i\le h/2$ and the assumed inequality. Uniqueness of the radix expansion gives the identity in $X$. Hence $-x/b$ is a root of $F$. The reciprocal polynomial is monic, so $t=b/x$ is a positive integer. The identity becomes $(1+tX)B=iF$. Comparing coefficients from the constant term up proves $i\mid B_l$ for every coefficient, so $G=B/i$ has nonnegative integer coefficients and $F=(1+tX)G$.

In fact $Ux=iG(P)$, $Uy=jG(P)$ by the exact ratio. The actual row is
$$
                   G(P)\big(Pi,Pj,1+(t-i-j)P\big).
                   \tag{2}
$$
Its positivity gives $i+j\le t$, and $h=t+G_1\ge t$, where a missing coefficient is zero. This proves (1).

## A second range for quadratic weights

Suppose $DR=MP+h$ with $M\ge1$, $0<h<P$, $P>D$, $Dh\le P$, and
$$
                         D^2(h^2+4M)\le32P.
$$
A split assignment requires positive integers $t,u$ with $t+u=h$, $tu=M$, and disjoint positive digit subcounts $i,j$ of $h$ with $i+j\le t$. This version allows $M\ge P$ if the inequalities permit it. In particular, if the grouped assignment is excluded as above, a negative or nonsquare discriminant $h^2-4M$ certifies the weight. **Every $h\le3$ certifies the weight throughout this range**, including the square-discriminant cases.

With $b,x,i$ as before, reducing $b^2N$ by $bP\equiv-x\bmod V$ gives
$$
                    V\mid b^2-hxb+Mx^2.
$$
The expression is strictly less than $D^2(h^2+4M)/16\le2P$ since $b<Dh/4$ and $x<D/2$. Completing the square gives a lower bound $-h^2x^2/4>-D^2h^2/16\ge-2P$. As $V>2P$, it is zero. Thus $t=b/x$ is an integer root of $T^2-hT+M$ and $u=h-t>0$. Equation (2) follows from $iV=x(1+tP)$ and the exact ratio, with $G(P)=1+uP$.

At one or two digit units we have $i+j=h$ if the two residues both exist, which contradicts $i+j\le t<h$. For $h\le3$ grouping is impossible by $(D,h)\le3$. The only unresolved allocation at three is $t=2,u=i=j=1$. The row in (2) is then $(P+1)(P,P,1)$, so $U\mid P+1$. For any prime $\ell\mid U$, the first and last parts sum to $(P+1)^2$. If $\ell^c\parallel P+1$, their sum is divisible by $\ell^{2c}$ but the last part is not. They carry in base $\ell$, contradicting the assignment of $\ell$ to that row. This accounts for every positive pair of roots at $h\le3$.

Run `python3 check_predecessor_digit_splits.py` and `lean PredecessorDigitSplits.lean`. The Lean file checks the quadratic divisibility and the inequality from positivity in (2). The Python check also tests the full three-way assertion on 311,972 exact coefficients at 168 selected weights through 250, enumerates all 7,770 digit polynomials of degrees two through five with digits at most five and constant one, and checks 11,350 grouped binomial cuts. At weight 241,860, two actual rows avoid the two assigned predecessor primes and realize a permitted polynomial factor. The digit, root and assignment arguments above are ordinary mathematical proofs. [Lean lemmas](PredecessorDigitSplits.lean) · [Independent arithmetic checks](check_predecessor_digit_splits.py)
