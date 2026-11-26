# Products of parts across the two endpoint rows

Write $K_A$ for the trinomial coefficient of a positive row $A$ of weight $N$. For **any** two rows $A,C$, choose any part $a\ge2$ of $A$ and any part $b\ge2$ of $C$. Set $D=\gcd(N,K_A,K_C)$ and $J=\gcd(N-1,K_A,K_C)$. Then
$$9D^{13}J^6(ab)^9\ \ge\ N^9(N-1)^2. \tag{1}$$
In particular the common divisor $DJ$ of the two coefficients and the endpoint coefficient grows without bound along every sequence with $N\to\infty$ and $ab=o(N^{11/9})$.

In every counterexample, the [endpoint theorem](all_weights_endpoint_reduction.md) supplies the row $(1,1,N-2)$. For **each of the nine choices** of $a,b$ from the other two rows, $t=ab/N$ is a positive integer with
$$N<t^{9/2},\qquad t>32768.$$
The first bound is analytic and gives $ab>N^{11/9}$. The second comes from the finite enumeration below, which imposes no upper cutoff on $N$. More generally, for every $\epsilon>0$ the number of possible exceptional weights with some $t\le Y$ is $O_\epsilon(Y^{2+\epsilon})$.

## The divisor identity

Use the [retained prime-power assignment](endpoint_pair_bound.md) for $A,C$. It gives positive integers $U,V,S,T$ with $DUV=N$, $JST=N-1$. The entries in $A$ are multiples of $U$, those in $C$ are multiples of $V$. At each retained prime power in $S$ the three actual residues of $A$ are $0,0,1$ in some order. The same holds for $T$ and $C$. Write $x=a/U$, $y=b/V$, so that $t=xy=Dab/N$ throughout the general proof.

Split $S=FB$ according as the chosen entry $a$ has residue 0 or 1. Likewise split $T=HL$ using $b$. We have $F\mid x$, $H\mid y$. Write the normalized rows as
$$(x,Bf_1,Bf_2),\qquad (y,Lg_1,Lg_2).$$
Each prime power in $F$ divides one of $f_1,f_2$, so $F\mid f_1f_2$. Similarly $H\mid g_1g_2$. Set $l=f_1+f_2$, $r=g_1+g_2$, $d=FH$, $G=DJ$. Then $l\ge2\sqrt F$, $r\ge2\sqrt H$, $d\mid t$, and
$$DV=x+lB,\qquad DU=y+rL,\qquad N-1=JdBL.$$
Putting $e=Gd-lr$ and multiplying the two row-sum identities gives
$$eBL=t-D+lyB+rxL,\qquad
(eB-rx)(eL-ly)=Gdt-De. \tag{2}$$
Both factors are positive. In fact $(eB-rx)L=D(b-1)$ and $(eL-ly)B=D(a-1)$. In particular $e\ge1$.

Denote these positive integers by $h_1,h_2$. Then $h_1h_2<Gt^2$ and $h_1+h_2\le h_1h_2+1$. Also
$$rx<Gd x/l\le Gt\sqrt{x}/2,\qquad ly<Gt\sqrt{y}/2.$$
For the first estimate use $F\le x$, $H\le y$ and $l\ge2\sqrt F$. Interchanging the roles gives the second. Expanding $(h_1+rx)(h_2+ly)$, using $lr t<Gt^2$, gives
$$BL\le e^2BL<2Gt^2+\tfrac12Gt^{3/2}(Gt^2+1).$$
Since $d\le t$, it follows that $N-1\le3JG^2t^{9/2}$. Substitution proves (1).

In a counterexample $D=J=1$, and the two scales $U,V$ are at least 3 because each is the sum of three positive normalized parts. Here $t=ab/N$ and $d>lr\ge4\sqrt d$, giving $t\ge d\ge17$. The sharper preceding bound gives
$$N<1+2t^3+\tfrac12t^{5/2}(t^2+1)<t^{9/2}.$$
The last inequality follows by dividing by $t^{9/2}$ and using $t\ge17$.

If $N-1$ has exactly three distinct prime divisors, the [three-factor assignment](endpoint_factors.md) puts all of $N-1$ on one row. Call it $A$, so $H=L=1$ and $N-1=FB$. Then $(F-Ul)B=Ux-1>0$. Hence $U<F/l\le\sqrt F/2$ and $N-1=FB<F^{3/2}x/2\le x^{5/2}/2$. In this case **every** cross-row pair satisfies $ab>N(2(N-1))^{2/5}$.

## Counting and testing possible weights

For fixed $t$, take $xy=t$ and $F\mid x$, $H\mid y$. These choices cost $t^{o(1)}$ by the divisor bound. The positive integers $l,r$ satisfy $lr<d\le t$, giving $O(t\log(2t))$ possibilities. Each positive factor pair of $dt-e$ in (2) fixes $B,L$, and then $U,V,N$. This costs another $t^{o(1)}$. Summing over $t\le Y$ gives the stated count of weights.

The C++ enumerator groups by $d=FH$ and orders the rows so that $F<H$. They are coprime, and equality would force $d=1$. It enumerates all unordered splits $f_1+f_2=l$ with $F\mid f_1f_2$ and similarly for $H$. Here $\gcd(l,F)=1$ since $F\mid x$, $\gcd(F,B)=1$, and $V=x+lB$ is coprime to $F\mid N-1$. Thus the allowed residues are $f_1\equiv l i\bmod F$ for all idempotents $i$ modulo $F$ (independently 0 or 1 at each full prime power). Equation (2) reconstructs all possible $B,L$ by enumerating the smaller of its two factors through $\lfloor\sqrt{dt-e}\rfloor$. The code requires $F,H,B,L$ pairwise coprime and $\gcd(U,V)=1$. It then checks Kummer's criterion for the actual second row at every prime dividing $VHL$, followed by the first row at every prime dividing $UFB$. The prime sieve suffices since $B,L\le2t^2$ by (2) and $U,V\le t+2t^3$. The weight and the parts use 128-bit integers.

The independent Python search starts from the actual parts. Since $F\mid x$, $H\mid y$, $B\mid a-1$, $L\mid b-1$, we have $N-1\mid t(a-1)(b-1)$. Reducing $ab=tN$ modulo $N-1$ shows that
$$k=\frac{t(a+b-t-1)}{N-1}$$
is an integer with $0<k<t$. For positivity use $a>t$, $b>t$, and for the upper bound use $(N-a)(N-b)>0$. In particular $ab(a+b)>N(N-1)$, and the two rows cannot both have a part at most $(N(N-1)/2)^{1/3}$. Direct expansion now gives
$$(ka-t^2)(kb-t^2)=t(t-k)(t^2-k).$$
Both factors are positive since $(ka-t^2)(N-1)=t(a-1)(a-t)$ and likewise for $b$. The verifier factors the right side for each $k$, recovers $a,b,N$, and assigns the full powers of $N-1$ by the actual residues of $a,b$. It loops directly over the possible remaining parts and uses factorial floors for its final divisibility tests.

| Maximum $t$ | Systems of scales and normalized sums | Row pairs including assigned factors | Second rows passing assigned bases | Pairs passing all assigned bases | Largest reconstructed $N$ |
|---:|---:|---:|---:|---:|---:|
| 200 | 781 | 852 | 0 | 0 | 2,355,278,661 |
| 800 | 15,755 | 20,277 | 0 | 0 | 2,832,795,555,731 |
| 10,000 | 1,507,916 | 2,601,378 | 51 | 0 | 291,872,273,441,283,151 |
| 32,768 | 10,259,193 | 20,157,972 | 396 | 0 | 59,143,769,369,569,335,781 |

The independent search matches the first two lines, including the largest weights. At $t=1542$ alone it also recovers the first second row to pass its assigned bases, among 15 systems and 16 pairs. At the full bound, an independent factorial-floor audit of all 405 pairs attached to the 396 passing second rows finds a common prime in every pair. They represent 329 distinct weights. The prime 2 suffices for every one, and the minimum of the two coefficient valuations at 2 is always at least two. A weight may occur in several systems. The largest reconstructed weight only describes this search. The Python program additionally checks (1), the integrality of $Dab/N$, and the four-part divisibility in the [two-prime-power criterion](cross_prime_power_gaps.md) using exact binomial products in 1,357,965 comparisons through weight 40.

~~~sh
mkdir -p .work
c++ -std=c++17 -O3 check_cross_parts.cpp -o .work/check_cross_parts
.work/check_cross_parts 32768 --audit > .work/cross_parts_32768_audit.txt
.work/check_cross_parts 800
# Python verifier requires python-flint.
python3 check_cross_parts.py --audit .work/cross_parts_32768_audit.txt
python3 check_cross_parts.py --limit 800 --weight 40
python3 check_cross_parts.py --first 1542 --limit 1542 --weight 3
~~~
