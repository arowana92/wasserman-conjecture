# Three digit units above a predecessor power

Let $p,q$ be distinct primes, $P=p^a$, $Q=q^b$, $a,b\ge1$, and $N=1+DPQ$, where $D$ is a positive integer coprime to $pq$. For a power $L=p^k$, put $h_p=DQ\bmod L$ (least nonnegative residue). We use a layer when $Dh_p\le L$. Define $h_q$ symmetrically using $DP$ modulo a power of $q$. A base-$p$ digit of value $d$ in position $j$ contributes $d$ units of size $p^j$.

**One-prime theorem.** Suppose a layer at $p$ has exactly three units. Write their sizes in order as $1,t,tu$ and set $m=1+\min(t,u)$. If
$$
                     h_p\nmid D,\qquad D\le\min\{5h_p,\ m(4m+1)\},                 \tag{1}
$$
then every three positive trinomial coefficients of weight $N$ have a common prime divisor. This holds simultaneously for **every** $a\ge1$.

**Two-prime theorem.** Suppose there are layers at $p,q$ each having at most three units. Any counterexample requires
$$
 h_p=h_q=s,\qquad s\mid D,\qquad s-1=p^i+p^j=q^k+q^l.                    \tag{2}
$$
The powers in each sum are the two digit units remaining after removing one unit of size one from $s$. Moreover some unitary divisor $U\ge3$ of $N$ divides $s-1$. One row has a part $(s-1)(N-1)/s$ and avoids $p,q$ and every prime dividing $U$. Here a unitary divisor satisfies $\gcd(U,N/U)=1$.

For primes below 257, **1,408 of the 1,431 unordered base pairs certify the conjecture for every such weight**. The other 23 pairs require a specific common residue from the table below. There is no bound on $D$, the exponents or the high digits.

**Three-prime theorem.** Let three distinct primes divide $N-1$, with full powers $P_1,P_2,P_3$. Suppose at each prime $p_i$ a layer $L_i$ has a remainder $h_i=((N-1)/P_i)\bmod L_i$ of at most three digit units, and
$$
              L_i\ge h_i\max_{j\ne i}\frac{N-1}{P_iP_j}\qquad (i=1,2,3).
$$
Then the conjecture holds at $N$. All three primes can be arbitrarily large. The proof below uses a separate result about interior binomial cuts and does not require the finite prime table.

## Why the assignments force (1) and (2)

The [endpoint reduction](all_weights_endpoint_reduction.md) lets us assume one row is $(1,1,N-2)$. Assign each prime of $N(N-1)$ to one of the other two rows whose coefficient it does not divide. The respective assigned full weight powers have products $U,V\ge3$, with $UV=N$, and the two rows have parts divisible by $U,V$ respectively. At an assigned predecessor power the parts have residues $0,0,1$. Carry freedom above that power splits the digits of the quotient. We use the following consequences of the [low-layer identities](predecessor_digit_splits.md#the-carry-free-cuts).

If one row takes both $p,q$, one of its parts is $xPQ$ with $0<x<D$ and $U\mid x$. There is a proper, nonempty subcollection of the units of $h_p$ with sum $v$ such that $Dv=xh_p$. Since $(D,U)=1$, we have $U\mid v$. A single $p$-power is impossible for $v$. With at most three units, $v$ therefore takes two and its complement $h_p-v$ is a single $p$-power. The identity $D(h_p-v)=(D-x)h_p$ gives $h_p\mid D$ since $p\nmid h_p$. With a three-unit bound in both bases we get, for some $i',k'\ge0$,
$$
                \frac{D-x}{D}=\frac{p^{i'}}{h_p}=\frac{q^{k'}}{h_q}.
$$
Taking valuations at $p,q$ gives $i'=k'=0$ because the respective denominators are units at their own bases. Thus $h_p=h_q=s$, $v=s-1$, $x=D(s-1)/s$ and (2) follows.

If the two primes go to different rows, two parts in the first are $UPx,UPy$ and two in the second are $VQz,VQw$. All four multipliers are positive. Positivity of the third parts gives $V\ge P(x+y)+1$, $U\ge Q(z+w)+1$ and hence
$$
                             (x+y)(z+w)<D.                                  \tag{3}
$$
Each full prime power of $D$ divides a displayed multiplier of its assigned row, since the two pairs of zero positions must intersect. Write $(x,y)=g(A,B)$ in lowest ratio. The low-layer identity $xj=yi$ says that $A:B$ is the ratio of two disjoint positive subcollections $i,j$ of the units of $h_p$. Let $F$ be the part of $AB$ coprime to $p$ (remove every factor $p$). Then $D\mid gFzw$. This follows prime power by prime power, since $\operatorname{lcm}(x,y)=gAB$ and $p\nmid D$.

For (1), if $i,j$ are singletons then $F=1$ and $A+B\ge m$. From (3) and $4zw\le(z+w)^2$ we get $z+w>4(A+B)$, so $D>m(4m+1)$. Otherwise the two subcollections use all three units. Their gcd divides both their sum $h_p$ and the singleton, hence is one. Thus $A+B=h_p$, and $F<h_p$ since one reduced term is a $p$-power. Now (3) gives $z+w>4h_p/F>4$, so $D>5h_p$. This excludes either split under (1).

For (2), write $(z,w)=g'(C,E)$ in lowest ratio, and let $F'$ be the part of $CE$ coprime to $q$. At each layer at least one of the two subcollections is a singleton. Consequently $F<A+B$ and $F'<C+E$, also when one subcollection leaves some units unused. The same assignment gives $D\mid gg'FF'$, contradicting (3). The grouped row is therefore compulsory.

## An unrestricted binary exponent with arbitrarily many predecessor primes

For every $g\ge3$, set $t=2^g$, $D=t^2+1$, $h=1+t+t^2$ and $L=2^{4g+1}$. There is one even class $b_0\bmod 2^{4g-1}$ for which $D3^{b_0}\equiv h\pmod L$. Indeed $h/D\equiv1\pmod8$ as a modular quotient, and 9 generates the units congruent to one modulo eight because successive squaring gives $v_2(9^{2^j}-1)=3+j$. Since $Dh<L$, $D<h$, and (1) holds with $m=t+1$, the conjecture holds throughout
$$
             N=1+(2^{2g}+1)2^a3^b,\qquad a,b\ge1,\qquad b\equiv b_0\pmod{2^{4g-1}}.
$$
For $g=3$ this gives $N=1+65\cdot2^a3^{386+2048j}$ for all $a\ge1$, $j\ge0$. For $g=9$ the representative is $30{,}828{,}923{,}520$ modulo $34{,}359{,}738{,}368$. In these families the occupied positions $0,g,2g$ all occur before an initial one- or two-unit residue can satisfy the required size bound. In fact the respective moduli are at most $t$ and $t^2$ while $D=t^2+1$.

The number of distinct prime factors of $N-1$ is unbounded. Set $g=3^n$, $n\ge1$, and $D_n=4^{3^n}+1$. Each step factors as $D_n=D_{n-1}(A^2-A+1)$ with $A=4^{3^{n-1}}$. The two factors have gcd dividing three, and $3\nmid A+1$. The second factor exceeds one, so each step introduces a new prime. Starting at $D_0=5$ gives at least $n+3$ distinct primes after including 2 and 3.

## Classification at all 54 primes below 257

The table gives $(q,s-1)$ for each possible smaller prime $p$. An omitted pair is excluded. In a listed case, $s$ must still divide $D$, both remainders must equal $s$, and $N$ must have a unitary divisor $U\ge3$ dividing $s-1$.

| $p$ | Possible $(q,s-1)$ |
|---|---|
| 2 | $(3,10),(5,6),(5,130),(11,12),(11,132),(17,18),(19,20),(23,24),(47,48),(67,68),(71,72),(79,80),(131,132),(191,192)$ |
| 3 | $(11,12),(29,30),(83,84),(89,90),(107,108),(251,252)$ |
| 5 | $(29,30),(149,150)$ |
| 11 | $(131,132)$ |
| 13 | $(47,2210),(181,182)$ |

Here is a finite bound proving that the table covers arbitrary digit positions. Write $p^i+p^j=q^k+q^l$, with $i\le j$, $k\le l$, and $B=\max(j,l)$. For $B\ge100$, both smaller terms are at most $KB$, $K=28{,}561$. To see this, $q^k$ divides $p^{j-i}+1$. For odd $q$ and positive gap, if there is a least $d>0$ with $p^d\equiv-1\pmod q$, every eligible gap is $d$ times an odd integer. Lifting at $q$ bounds its full power by $q^{v_q(p^d+1)}B$. When no such $d$ exists, or the gap is zero, the odd prime contributes one. For $q=2$, the 2-part of $p^e+1$ is at most the larger of two and the 2-part of $p+1$, by parity of $e$. The verifier checks the maximum over ordered prime pairs. Interchanging bases gives the other bound.

Interchange sides so $j=B$. Then $|p^B-q^l|\le KB$, and $l>0$ for $B\ge100$. Apply Matveev's bound to $\Gamma=q^l p^{-B}-1\ne0$. Its two-logarithm constant is $1.4\cdot30^5\cdot2^{4.5}\log p\log q<3\cdot10^{10}$. It gives
$$
             B\log 2 < 3\cdot10^{10}(1+\log B)+\log(KB),
$$
so $B<M=10^{14}$ by monotonicity of the right-hand side divided by $B$. [Matveev's original theorem](https://www.mathnet.ru/eng/im314) · [Languasco, Luca, Moree and Togbé, Theorem 2.1](https://pure.mpg.de/rest/items/item_3607286_6/component/file_3699621/content)

For $\theta=\log p/\log q$ the verifier certifies a reduced fraction $a'/b'$ with $2M<b'<2^{70}$ and $|b'\theta-a'|<1/(2M)$ for each of the 2,862 ordered pairs. The rational log enclosures and determinant argument are detailed in the [two-logarithm check](seven_mixed_primes.md#a-bound-on-the-exponent). Since $B<b'$, the integer $Ba'-lb'$ is nonzero and $|B\theta-l|>2^{-71}$. But $\min(p^B,q^l)>p^B/2$, giving $|B\theta-l|<3KB/2^B$ by the mean value theorem for the logarithm and $2/\log2<3$. At $B=100$ this is already below $2^{-71}$ and it decreases. Thus $B<100$.

The enumeration intersects the two-power sums below this exponent bound, requires $s$ to be a unit at both primes and discards $s-1$ whose part coprime to $pq$ is below three. A separate pass extracts valuations and tests whether the remaining cofactor minus one is a prime power, including exponent zero. It returns the same 25 triples. For example, at $(2,3)$ only $s=11$ survives, with $U=5$. Write $T=(N-1)/11$. The forced part and its complement are $10T,T+1$. Carry freedom at 5 forces $T\equiv4\pmod{25}$. With $h=v_5(T+1)\ge1$, the residue of $10T$ modulo $5^{h+1}$ is $5^{h+1}-10$. The positive residue of $T+1$ must be less than ten, hence equals five.

## Three simultaneous binomial cuts

Let $P_i=p_i^{a_i}\ge2$ be distinct prime powers at distinct bases, and suppose a positive integer $W$ has the form
$$
 W=P_i+t_i P_i L_i,\qquad t_i\ge1,\quad L_i=p_i^{k_i},\quad k_i\ge1.
$$
For any pair $i,j$ satisfying
$$
                      t_i t_j(P_i+P_j+1)<W,                              \tag{3a}
$$
there is at most one unordered interior cut of $W+1$ whose binomial coefficient avoids both primes. Here interior means both parts are at least two. With $P_i>P_j$, that cut must contain $(W-P_j)/(P_i-P_j)$, and $P_i-P_j$ must divide both $t_i,t_j$. With three such powers, if (3a) holds on two edges joining all three indices, **every interior binomial coefficient in row $W+1$ is divisible by at least one of the three primes.**

Indeed a cut $a$ avoiding the specified primes has, at each of their layers,
$$
 a=k_i' P_iL_i+e_iP_i+f_i,\qquad 0\le k_i'\le t_i,\qquad e_i,f_i\in\{0,1\}.
$$
The low digits come from $P_i+1$ in $W+1$. Put $\lambda_i=k_i'/t_i$. Then
$$
                 a=\lambda_i W+(e_i-\lambda_i)P_i+f_i.                   \tag{4}
$$
This gives $W|\lambda_i-\lambda_j|\le P_i+P_j+1$. Under (3a) the two rationals are equal, since unequal rationals of those denominators differ by at least $1/(t_it_j)$. Write $\lambda$ for the common value. The case $e_i=1$, $e_j=0$ would give $(1-\lambda)P_i+\lambda P_j=f_j-f_i\le1$, whereas the left side is at least two. Hence the $e_i$ agree. If $f_i=f_j$, distinctness of the powers in (4) gives $\lambda=e_i$, and $a$ is one of $0,1,W,W+1$. Otherwise exchange the cut and its complement if necessary to arrange $e_i=e_j=0$. For $P_i>P_j$ we get $\lambda(P_i-P_j)=f_i-f_j=1$ and $f_j=0$. This gives the stated candidate and divisibility of the two denominators. The boundary $P_i-P_j=1$ is permitted. For example, $(W,P_i,P_j)=(12,4,3)$ gives the surviving cut $9+4$, and $(105,7,5)$ gives $50+56$ at difference two.

For **three** indices the connected comparisons give one $\lambda$ and one $e$. Two of the three $f_i$ agree, forcing $\lambda=e$. Equation (4) now gives an endpoint at once. For example, the three primes $2,5,17$ cover all interior cuts at $W+1=4931$, and $2,3,7$ do so at $W+1=59053$. The verifier checks every cut at these two weights by exact factorial valuations. Lean checks the common-denominator argument from (4) to an endpoint. [Lean lemmas](ThreeUnitLayers.lean)

To apply this to Wasserman's conjecture, assume a counterexample to the three-prime theorem. Apply (2) to each pair with the same fixed prime assignment. All three primes go to one row, and their remainders equal the same $s\ge3$ dividing each pair's cofactor. Since $s$ is coprime to all three bases, we can put
$$
           N-1=sW,\qquad W=E P_1P_2P_3,\qquad E\ge1.
$$
The grouped row contains $(s-1)W$, so its other two parts add to $W+1$. Both are positive multiples of $U\ge3$. Their binomial cut avoids each of the three primes, since its coefficient divides the full row coefficient. From $sW/P_i\equiv s\pmod{L_i}$ and $(s,p_i)=1$ we obtain $W/P_i=1+t_iL_i$ with $t_i\ge1$. The stipulated margin gives
$$
                L_i\ge s^2 E\max_{j\ne i}P_j,\qquad
                t_i<\min_{j\ne i}P_j/s^2.
$$
In descending order $P_1>P_2>P_3$ this implies $t_i t_j<P_2P_3/s^4$ for each pair. Since $2P_1+1<s^4 E P_1$, (3a) holds on every edge. The three-prime binomial result excludes the two remaining parts. This proves the theorem without any bound on the primes or on the number of other factors of $N-1$.

There is also a two-prime consequence for the exceptional residue in (2). Write $W=(N-1)/s=EPQ$, $E=D/s$. The two gap hypotheses give $t_p<Q/s^2$ and $t_q<P/s^2$ whenever the two quotients $W/P-1,W/Q-1$ are positive, as they are here. Thus $P+Q+1\le s^3D$ suffices for (3a). In that range a counterexample requires
$$
       |P-Q|<\min(P,Q)/s^2,\qquad U\mid\min(P,Q)+1.
$$
Indeed $|P-Q|$ divides both $t$'s. Also $U$ divides both remaining parts and their sum $W+1$. Multiplying the forced cut by $|P-Q|$ shows $U\mid W-\min(P,Q)$, giving the second assertion.

## Counts in exponent progressions at the 23 remaining pairs

Fix $D,p,q$ with $p,q<257$. Count only $(a,b)$ satisfying the hypotheses of the two-prime theorem. In a Cartesian product of two arbitrary positive-integer arithmetic progressions of length $H$, the number of possible exceptional weights is $O(H^{1.619})$. On two synchronized progressions $(a_j,b_j)$ of length $H$ the count is $O(H^{0.776})$. In any progression of $H$ values of one exponent, at most $O(H^{0.876})$ values can support even one such exception as the other exponent ranges over all positive integers. The implied constants are absolute, including uniformity in $D$, starting exponents and positive steps.

For each surviving table entry put $w=s-1$, $E=D/s$, $W=E p^a q^b$. The two cuts $wW,W+1$ add without a carry at both bases. Removing the shift at $p$ gives a carry-free addition $wm+m$ with $m=Eq^b$. Independently the $q$-condition gives the same addition for $m=Ep^a$. Thus we can bound one exponent at a time, also when the other exponent or its admissible layer varies.

At base $p$, let $c$ be the incoming carry of the multiplication by $w$, so $0\le c<w$. An input digit $d$ is allowed if $d+((wd+c)\bmod p)<p$. It leads to state $\lfloor(wd+c)/p\rfloor$. When $p\nmid w$, each state allows exactly $(p+1)/2$ digits. Both $wd+c$ and $(w+1)d+c$ permute the residues modulo $p$. Subtracting the sum of residues after addition from the sum before it, then dividing by $p$, counts exactly $(p-1)/2$ overflows. When $p\mid w$, the verifier iterates the transition for 64 digits starting with the vector of ones and uses its largest entry $A$. In either case write $n$ for the block length (one or 64) and choose the least integer $b_p$ such that $A^{1000}\le p^{n b_p}$. The number of admissible $j$-digit words from any incoming state is $O(p^{j b_p/1000})$ by splitting into blocks. The finite table gives
$$
        \max b_p=876,\qquad \max(b_p+b_q)=1619,\qquad \max\min(b_p,b_q)=776,
$$
where both orientations are included for the first maximum.

To transfer this estimate to any exponent progression, write its positive step as $d_0$ and split its indices modulo $p-1$, or modulo 2 for $p=2$. In each class the ratio of successive $m$ values is $z=q^{d_0(p-1)}$ at odd $p$, or $q^{2d_0}$ in binary. Let $f=v_p(z-1)$, so $f\ge1$ and in binary $f\ge3$. Lifting gives order exactly $p^j$ modulo $p^{f+j}$. All $m$ in this class share the low $f$ digits and the same entering multiplication carry above them. For the least $j\ge0$ with $p^j\ge H$, the next $j$ digits distinguish all terms of the class. This bounds the number of eligible exponents by $O(H^{b_p/1000})$ independently of $d_0,E$ and $f$. Multiply the two bounds on a Cartesian product or take the smaller on synchronized progressions. Sum over the finite table.

Run `python3 check_three_unit_layers.py` and `lean +leanprover/lean4:v4.36.0 ThreeUnitLayers.lean`. The Python verifier checks the rational separation bounds, computes the possible cuts in two independent ways, and certifies the exponent bounds with integer powers. It also checks every triple of coefficients via their prime supports at 123 selected weights through 250, using 217,620 exact trinomial coefficients. [Reproducible checks](check_three_unit_layers.py)
