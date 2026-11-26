# Differences between parts in the two nonendpoint rows

Suppose three positive trinomial rows of weight $N$ have coefficient gcd one. By the [endpoint theorem](all_weights_endpoint_reduction.md), one row is $(1,1,N-2)$ up to permutation. In each of the other two rows select two unequal parts, and let their positive differences be $d,e$, ordered so that $d\le e$. Write $M=N-1$.

1. The number $t=de/N$ is an integer at least three. It satisfies
   $$M\mid t\big((d-e)^2-(t-1)^2\big),\qquad M\mid t(d^2-1)(d^2-t^2). \tag{1}$$
2. The equality $e-d=t-1$ is impossible for every prime $t$ and for every $t\le2660$, at all weights.
3. For $1\le B\le X$, the number of weights $N\le X$ of counterexamples admitting any such choice with $t\le B$ is $O_\epsilon(B^2 X^{1/3+\epsilon})$ for every $\epsilon>0$. For example, if $de\le N^{5/4}$ the count is $O_\epsilon(X^{5/6+\epsilon})$.

Every nonzero difference in either nonendpoint row also satisfies $d^7>27(N-1)^2$, even when the companion row has three equal parts. The selected positions are arbitrary.

## Congruences for the two differences

Use the [prime-power assignment](endpoint_reduction.md) with $UV=N$, $ST=M$. The rows $A,C$ have scales $U,V$ and avoid the primes of $US,VT$, respectively. The six pairwise coprime bins $s_1s_2s_3=S$ and $w_1w_2w_3=T$ record the exceptional position of residue 1 at each full predecessor prime power. Set $d=A_1-A_2=Ux$, $e=C_1-C_2=Vy$. Then $xy=t$ and $d,e>t$, since each difference is smaller than $N$. The bin congruences give
$$s_1\mid d-1,\quad s_2\mid d+1,\quad s_3\mid x,\qquad
w_1\mid e-1,\quad w_2\mid e+1,\quad w_3\mid y.$$
In particular $M\mid t(d^2-1)(e^2-1)$. Use $de=tN\equiv t\bmod M$ to obtain
$$k=\frac{t(d^2+e^2-t^2-1)}{M}\in\mathbb Z_{>0},\qquad
(k-2t^2)M=t\big((d-e)^2-(t-1)^2\big). \tag{2}$$
Multiplying the initial congruence by $d^2$ proves the second divisibility in (1). Its right-hand side is nonzero.

Exactly one of $S\ge V$ and $T\ge U$ holds. If both fail then $M\le(V-1)(U-1)$, and if both hold then $M\ge UV$. On the side where the inequality holds all three bins exceed one, since a part with an empty bin would be a positive multiple of the whole predecessor modulus after removing the row scale. Temporarily name that side $A$. Then $s_3\mid x$ rules out $t=1$. Write
$$A_3/U=fS/s_3,\qquad C_3/V=gT/w_3.$$
The product of the actual third parts is smaller than $N M$, so $fg<s_3w_3$. For $t=2$ we must have $x=s_3=2$, $y=w_3=1$, and $f=g=1$. Now $N,U,V,T$ are odd because $2\mid M$ and $2$ is in $S$. The normalized $C$ row consists of two consecutive integers and the odd integer $T$. Their sum is even, contradicting its value $U$.

There is also a bound at each individual weight. Let $L$ be the longest run of ones in the binary expansion of $N$, and set $H=2^{L+1}(2^{L+1}+1)$. Then
$$M<H^3t^6. \tag{3}$$
To see this, name the row assigned 2 as $A$, in this paragraph without ordering the two differences. Its smaller selected part $a$ and its larger part $a+d$ have disjoint bits. Write $P=2^h$ for the largest bit of the larger part. Then $P-d\le a<P$. For $\ell=\lceil\log_2 d\rceil\le h$, all bits in positions $\ell,\ldots,h-1$ of $a$ are ones, also present in $N$. Thus $a<2^{L+1}d$, which also holds for $\ell>h$. At each bin of $S$ at least one selected part divided by $U$ vanishes, giving $S\mid a(a+d)/U^2$ and $S<Hx^2$. On the other side $T\le w_3(U-2)\le y(U-2)$ from $C_3/V$, hence $V<Sy$. Finally $T\mid y(e^2-1)<y^3V^2$, so $M<H^3x^6y^5\le H^3t^6$.

## Bounds for an arbitrary pair of rows

Every positive difference $d$ between parts of either nonendpoint row in a counterexample satisfies
$$d^7>27(N-1)^2. \tag{9}$$
This does not require a positive difference in the other row. More generally, for any two positive rows $A,C$ let $K_A,K_C$ denote their coefficients and set $D=\gcd(N,K_A,K_C)$, $J=\gcd(N-1,K_A,K_C)$. For a gap $d\ge2$ in $A$ one has
$$27(N-1)^2\le J^2 D^3 d^2\gcd(d,N/D)(d^2-1)^2. \tag{10}$$

Use the [retained prime powers](endpoint_pair_bound.md#retaining-partial-prime-powers) with $DUV=N$ and $JST=N-1$. At a power assigned to $S$ the actual residues are $1,0,0$ in some order, just as above. With $x=d/U$, this gives $S\mid x(d^2-1)$. The other three parts divided by $V$ sum to $DU$, and their product is divisible by $T^2$, so $27 T^2\le(DU)^3$. Multiplying proves (10), since $U\mid\gcd(d,N/D)$. Also $DU\ge3$ and $U\mid d$, so $Dd\ge3$ even for $d=1$. In a counterexample $D=J=1$, and (10) gives (9).

For any second positive gap $e$ in $C$, put $r=de/N$. Then $\tau=Dr$ is an integer, and
$$M\mid J\tau\big(D^2(d-e)^2-(\tau-D)^2\big),\qquad
M\mid J\tau(d^2-1)(D^2d^2-\tau^2). \tag{11}$$
Indeed $\tau=(d/U)(e/V)$ and $ST\mid\tau(d^2-1)(e^2-1)$ by the bin congruences. Use $Dde=\tau N\equiv\tau\bmod ST$ after multiplying the latter divisibility by $D^2$ or $D^2d^2$, respectively. In particular, when $(d-e)^2\ne(r-1)^2$ the endpoint coefficient and the two row coefficients have a common divisor of size at least
$$DJ\ \ge\ \left(\frac{N-1}{r\,|(d-e)^2-(r-1)^2|}\right)^{1/3}.$$
This follows from the first divisibility of (11) since $JD^3\le(DJ)^3$.

As a spacing consequence of (2), if $M>t(t-1)^2$ and $e-d\ne t-1$ then
$$e-d\ge\sqrt{(t-1)^2+M/t}.$$
A smaller value would make the nonzero multiple of $M$ in (2) strictly between $-M$ and $M$. This applies in particular for all prime $t$ and for all $t\le2660$ under the stated size hypothesis.

## When the right-hand side of (2) vanishes

Resume the order $d\le e$ and suppose $e=d+t-1$. Two factorizations follow:
$$tM=(d-1)(d+t)=(e-t)(e+1),\qquad M=(V-x)(U+y).$$
Both factors of the last product are positive. Substituting $d\equiv-1\bmod s_2$ and $e\equiv1\bmod w_1$ in the first products gives $s_2,w_1\mid2(t-1)$. Also $s_1\mid V-x$ and $w_2\mid U+y$, by $UV\equiv1\bmod M$. Define
$$b=(V-x)/s_1,\quad c=(U+y)/w_2,\quad bc=s_2s_3w_1w_3.$$
There are positive integers $m,h,n,j$ giving the exact normalized rows
$$(m s_1+x,\ m s_1,\ h s_1),\qquad (n w_2,\ n w_2-y,\ j w_2). \tag{4}$$
Here $b=2m+h$, $c=2n+j$, and
$$s_3\mid m,\ s_2\mid h,\ \gcd(m,s_2)=\gcd(h,s_3)=1,\qquad
w_3\mid n,\ w_1\mid j,\ \gcd(n,w_1)=\gcd(j,w_3)=1.$$
These follow by examining the exceptional position at each bin in (4). In particular $\gcd(b,s_3)=\gcd(c,w_3)=1$. If $t$ is prime, the side with all three bins greater than one must have third bin $t$, and the opposite third bin is 1. If $s_3=t$, then $b\mid s_2w_1\mid2(t-1)$ but $b\ge2t+s_2$. If $w_3=t$, use $c$ instead. This proves the prime assertion.

For fixed $t$, the quantities in (4) other than $s_1,w_2$ range over a finite list. The remaining two bins satisfy
$$xc w_2-yb s_1=t+1. \tag{5}$$
Indeed $Vy-Ux=t-1$, $V=x+b s_1$, and $U=cw_2-y$. In particular $x,y$ are coprime. Let $G=\gcd(xc,yb)$. Equation (5) gives $G\mid t+1$, so $G$ is coprime to $xy$ and $G=\gcd(b,c)$. The solutions for $(s_1,w_2)$ form an arithmetic progression of step $(xc/G,yb/G)$.

Write $L_0=s_2s_3$, $H_0=w_1w_3$. The reduced fractions $b/L_0$ and $c/H_0$ are reciprocals. Their cancellation factors $\alpha=\gcd(b,L_0)$ and $\beta=\gcd(c,H_0)$ each divide 2 by (4). Since $L_0,H_0$ are coprime, at most one is 2. Thus
$$b=\alpha H_0/\beta,\qquad c=\beta L_0/\alpha,\qquad G\le2. \tag{6}$$
For the last bound, if $\alpha=2$ then $b=2H_0$, $c=L_0/2$, with $H_0$ coprime to $c$. The other two cases follow directly.

### A binary cutoff

In the binary-avoiding row let two selected parts be $a,a+\delta$, and put $P=2^{\lfloor\log_2(a+\delta)\rfloor}$ and $Q=2^{\lceil\log_2\delta\rceil}$. Then $a<P\le a+\delta$. If $P\ge Q$, the bits occupied by $a$ from $\log_2Q$ through $\log_2P-1$ force the third part to have the form $lP+s$ with $0\le s<Q\le2\delta$.

If this is row $A$, write $a=P-u$ with $1\le u\le d$. Its third part is $ha/m$, so $(h-ml)P=hu+ms>0$. It follows that $P\le bd$, including the case $P<Q$. Hence $m s_1<b x$.

If it is row $C$, write its larger part as $P+u$, with $0\le u<e$. Its third part is $j(P+u)/n$. Thus $(nl-j)P=ju-ns$. If $P>ce$, this forces $j=nl$. An odd integer multiple of a positive integer shares its lowest occupied bit, so $l$ must be even. Consequently
$$n w_2<(c+1)y\qquad\text{or}\qquad 2n\mid j. \tag{7}$$
This also covers $P<Q$. Whenever $2n\mid j$, one has $n\mid c$.

### Exhaustion through 2660

For each $xy=t$ enumerate $s_3\mid x$, $w_3\mid y$ and coprime $s_2,w_1\mid 2(t-1)$. Retain the coprimality and nonempty-bin conditions, the factorizations $bc=s_2s_3w_1w_3$, (4), and solvability of (5). Try every possible assignment of 2. If a fixed bin of the opposite row contains 2, that choice is impossible. For row $A$ the bound $m s_1<b x$ gives a finite check. For row $C$ use the finite bound in (7) unless $2n\mid j$.

For $2\le t\le2660$ there are 54,150 parameter choices before assigning 2. The finite checks give 5,228 positive instances counting the choice of the binary row, with largest weight 48,688,385,386,086,601. Of these, 1,092 have coprime weight scales and pairwise coprime bins. Every one has a binary carry in the assigned row.

For a modulus $q=p^v$, enumerate all solutions of (5) modulo $q$ after dividing its coefficients and constant by $G$. Reject any solution assigning $p$ to multiple bins or both weight scales. If $p$ is assigned to a row, reject any carry among that row's first $v$ base-$p$ digits. These tests depend only on the residue class. Exhaustion at a single modulus excludes the entire progression.

A second digit test handles a prime divisor of a variable factor $F$ with unknown valuation. In (4) the actual row assigned to $F$ is a vector $\eta+Fv$, with $\eta=0$ for a weight factor and one entry of $\eta$ equal to 1 for a predecessor bin. Use $d=1+bys_1$ and $e=-1+cxw_2$ to obtain the following quotients.

| $F$ | $v$ |
|---|---|
| $U$ | $(m s_1+x,m s_1,h s_1)$ |
| $V$ | $(n w_2,n w_2-y,j w_2)$ |
| $s_1$ | $(by+Um,Um,Uh)$ |
| $w_2$ | $(Vn,Vn-cx,Vj)$ |

If $p\mid F$ and $p^a\parallel F$, the digits at position $a$ are the three residues $(F/p^a)v_i\bmod p$. Here $a\ge1$, so the offset digit contributes only at position zero. If every nonzero multiplier modulo $p$ gives digit sum at least $p$, that entire residue class fails at its first nonzero scale digit. Every coordinate of $v\bmod p$ is determined by the two bin residues modulo $p$. The second implementation obtains $v$ by exact division at a positive lift on (5); the table proves independence of that lift.

The remaining 114 choices give unbounded progressions. Every one is excluded by these tests at one prime. The second implementation uses the finite moduli 81, 25, 361 and 8192, together with the first-nonzero-digit test at their prime bases. Both programs raise an error on an unresolved case. The cutoff 2,660 is computational. At $t=2661$, for example, the current tests leave the shape $(x,y,b,c,m,h,n,j)=(2661,1,20,42,3,14,1,40)$ with $(s_2,s_3,w_1,w_3)=(7,3,40,1)$ unresolved. This specifies a progression on which higher carry conditions still have to be tested.

## Counting the possible weights

First suppose $k\ne2t^2$ in (2). Split into dyadic ranges $R\le e/d<2R$. The bound $de=tN$ gives $d\le\sqrt{tX/R}$. At each $d$ the second divisibility in (1) allows $O_\epsilon(X^\epsilon)$ choices for $M$, by the divisor bound for a nonzero integer of polynomial size in $X$.

Alternatively (2) gives $1\le k\ll t^2R$. Completing the square yields
$$(2t^2 e-kd)^2-(k^2-4t^4)d^2=4t^3(t^3+t-k). \tag{8}$$
If the right side vanishes, substitution of $k=t^3+t$ in the quadratic before completing the square gives $(e-td)(te-d)=0$. One gap would divide the other, making the larger divisible by both $U,V$ and hence by $N$, which is impossible.

For each remaining $k$ there are $O_\epsilon(X^\epsilon)$ possible pairs. To verify uniformity, write (8) as $u^2-Fd^2=H$ with $FH\ne0$. For positive square $F$ factor $H$ into the two linear factors. Otherwise write $F=Dv^2$ with $D$ squarefree. The algebraic integer $u+vd\sqrt D$ generates an ideal of norm $|H|$ in the full quadratic ring of integers. There are at most $\tau(|H|)$ such ideals, since a prime power $p^a$ contributes at most $a+1$ choices in the split case and at most one otherwise. Generators of a given ideal in (8) differ by a unit of norm one. In an imaginary quadratic field there are at most six units. In a real quadratic field both conjugates have absolute value at most $X^{O(1)}$ and at least $X^{-O(1)}$ by their nonzero integral product. Distinct positive norm-one units are spaced multiplicatively by at least $(3+\sqrt5)/2$, since a norm-one unit greater than one has integer trace at least 3. Including signs gives $O(\log(2X))$ choices per ideal, uniformly in the field. The elementary divisor bound absorbs this factor. [Sutherland, ideal norms and prime decomposition](https://math.mit.edu/classes/18.785/2021fa/LectureNotes6.pdf)

Taking the smaller of the two estimates per dyadic range and using $\min(a,b)\le a^{1/3}b^{2/3}$ with $a=t^2R$, $b=\sqrt{tX/R}$ gives $O_\epsilon(tX^{1/3+\epsilon})$ per $t$, absorbing the logarithmic number of ranges.

Finally suppose $k=2t^2$. For fixed $t$ there are $O_\epsilon(t^\epsilon)$ choices of $(x,y,s_2,s_3,w_1,w_3,b,c)$ by the divisor bound. In the finite binary cases count progression parameters directly, without choosing $m,n$. The bounds $s_1< bx/s_3$ and $w_2<(c+1)y/w_3$ with (5) allow at most
$$2+\frac{Gb}{s_3c}+\frac{G(c+1)}{w_3 b}\ll 1+t$$
parameters. Indeed $b\le2w_1w_3$, $c\ge2w_3+w_1$, and the interchanged inequalities hold. In the remaining case put $\lambda=c/n$, an integer divisor of $c$ by (7). Then $N+e=\lambda(P+u)$, $0\le u<e$. Since $e=d+t-1<2d$, this gives $|d^2-t\lambda P|\le5t\lambda d$, so $d$ is within $5t\lambda$ of $\sqrt{t\lambda P}$. Along (5) the step in $d$ is $tbc/G\ge t\lambda$. Thus each of the $O(\log(2X))$ powers $P$ allows only $O(1)$ values of $d$ for each $\lambda$, and $d$ determines $N$. The total in the equality case is $O_\epsilon(t^\epsilon(t+\log(2X)))$ per $t$. Summing both cases through $B$ gives the claimed bound.

## Reproduction

From this directory run the commands below. The first and third programs use python-flint. The second uses the Python standard library. It reconstructs the row sizes via (6), scans the two free bins directly in their finite ranges, and checks base-$p$ digits in the unbounded cases.

The third program checks (10) on 723,380 single gaps and (11) on 2,173,374 gap pairs from exact trinomial coefficients at weights through 45. The size bound for the common divisor exceeds one in 623,560 nonvanishing cases. It also checks the binary cutoff on 437,250 choices through weight 1,024 and (2) and (8) on 1,125 nonzero gap pairs from independently constructed low-residue assignments.

~~~sh
python3 check_resonant_row_gaps.py
python3 check_resonant_row_gaps_cross.py
python3 check_row_gap_products.py
~~~
