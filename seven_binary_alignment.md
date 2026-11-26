# Alignment of seven binary digits

Let $N$ be odd with exactly seven binary ones. In any counterexample, a nonendpoint row avoiding 2 partitions the ones into two pairs and a triple. If the units-place bit belongs to the triple, then the row, up to permutation, has the form
$$2^i(1+2^n),\qquad 2^j(1+2^n),\qquad 1+2^r(1+2^h),$$
where $n,i,j,r,h>0$. The only possibilities are
$$h=n+1,\ r\equiv0\bmod 2n,\qquad\text{or}\qquad
h=n-1,\ r\equiv n+1\bmod 2n.$$
All seven bits are distinct. The companion row must avoid every prime dividing the odd part of $N-1$.

**Consequence.** For all positive $a$ and $b\ge8$ and all six **distinct** nonnegative $h_i$, the conjecture holds at
$$N=1+\sum_{i=1}^6 2^{a+b h_i}.$$
This extends the [large-base theorem](seven_large_bases.md) to binary digit positions in a single residue class modulo $b$, including arbitrary outer shifts.

## Excluding unequal unshifted gaps

Use the [endpoint and prime-power assignment](endpoint_reduction.md), assigning 2 to the row in question. Its weight factor $U\ge3$ divides each part and is odd. Therefore no part can consist of a single binary bit. If the units-place bit lies in the triple, denote the gaps of the two unshifted pairs by $d\ge e>0$, the exponents of their lower bits by $i,j$, and put $g=\gcd(d,e)$, $F=2^g+1$. Both $d/g,e/g$ are odd and $U\mid F$. As in the odd-base proof let $S_0,S_d,S_e$ denote the factors from $m=\operatorname{odd}(N-1)$ assigned to the row with zero positions on both pairs, the first pair and the triple, or the second pair and the triple, respectively. Since $US_0\mid F$,
$$2^d\le m\le (F/3)^{3/2} S_d S_e.$$
Here $S_d\mid (2^d+1)/F$ and $S_d\mid 2^j(1+2^e)-1$, since a power assigned there has residue 1 on the other pair and is coprime to $F$. The interchanged divisibilities hold for $S_e$.

These are exactly the inequalities used in the [binary unequal-gap exclusion](equal_predecessor_gaps.md#base-two-with-large-common-gap-divisor), including its two finite verifiers. For clarity, its additional input $d<3e+g$ for unequal gaps follows here as well. If $e=1$, then $F=3$, $S_e=1$ and $m\le S_d\le(2^d+1)/3<2^d$. Otherwise $e\ge2$ and
$$S_d\ge\frac{3^{3/2}2^d}{\sqrt F(2^e+1)}>3\cdot2^{d-e-g/2}.$$
Use $\sqrt F\le\sqrt{3/2}\,2^{g/2}$ and $2^e+1\le(5/4)2^e$. Of the three residues $0,j,j+e$ on a circle of circumference $d$, the latter two delimit an arc of length $d-e$ whose longer empty piece has length at least $(d-e)/2$. Rotate the residues into an interval of length at most $(d+e)/2$. This multiplies the sum in the divisibility for $S_d$ by a unit power of 2 and introduces signs through $2^d=-1\bmod S_d$. If $d\ge3e+g$, the resulting signed sum is smaller in absolute value than $S_d$, so must vanish as an integer. Three signed binary powers can vanish only at positions $t,t,t+1$. This requires $e=0,\pm1\bmod d$, inconsistent with $2\le e\le d-2$. The unequal-gap exclusion now applies, giving $d=e=n$ and $S_d=S_e=1$.

## The gcd at equal gaps

Write $f=2^n+1$, $G=\gcd(f,1+2^r(1+2^h))$, $D=\gcd(f,1+2^h)$. Then $G,D$ are coprime divisors of $f$, $U\mid G$ and $S_0\mid D$. The last assertion follows because both pairs vanish at $S_0$, so the positive-bit terms of the triple sum to zero there. The span of the positive bits is at least $\max(n,h)$. Consequently
$$G^3D^2\ge27\cdot 4^{\max(n,h)}. \tag{1}$$
As $G^3D^2\le f^3$, this requires $h<3n/2$.

If $D=1$, rotate the residues $0,r,r+h$ modulo $n$ past their largest empty arc. The resulting span is at most $2n/3$. The signed sum has absolute value strictly less than $3\cdot2^{2n/3}$, since either its terms occupy more than one position or the span after rotation is zero. It vanishes modulo $G$ and by (1) must vanish as an integer. Hence it also vanishes modulo $f$, so $G=f$. Two of the three residues coincide. If $n\mid h$, then $h=n$ by the span restriction, contradicting $D=1$. Thus $r$ or $r+h$ is divisible by $n$. If this is an odd multiple of $n$ then the other term is a unit modulo $f$ after cancellation, impossible. If $r=0\bmod 2n$, the equation gives $2^{h-1}=-1\bmod f$. If $r+h=0\bmod 2n$, it gives $2^{r-1}=-1\bmod f$. The order of 2 modulo $f$ is $2n$, because it divides $2n$ and no positive exponent at most $n$ gives 1. This proves the stated two alternatives. Also $S_0=1$, so every prime of $m$ is assigned to the companion.

## Excluding a nontrivial $D$

If $D>1$, put $g=\gcd(n,h)$. Then $D=2^g+1$ and $n/g,h/g$ are odd. The case $n\mid h$ has $h=n$ and $G=1$, violating (1). Thus $t=n/g\ge3$. Since $2^n=-1\bmod G$, we may use the divisibility
$$G\mid 2^{r+n}(1+2^h)-1.$$
For $g\ge29$ the [binary resultant bound](equal_predecessor_gaps.md#base-two-with-large-common-gap-divisor) gives $\log_2 G\le n/3+3t$. Since $(t-2)g\ge9t$,
$$3\log_2G+2\log_2D < n+9t+2g+2 \le 2n+2,$$
contradicting (1). For $g\le28$, [inequality (6) for a single binomial gap](equal_predecessor_gaps.md#approximating-two-exponent-ratios-in-base-two) gives $\log_2 G\le160+3n/5$, using $K=10$ in both cases of that bound. With $\log_2 D<29$, (1) would imply $n<2700$.

Two finite integer verifications exclude the remaining range. For the first, if $L$ is the largest empty circular arc among the three residues, rotating as above gives $G\le3\cdot2^{n-L}$. Indeed an integer zero would imply $f\mid1+2^r(1+2^h)$, impossible since $D>1$ divides $f$ and the expression is 1 modulo $D$. Consequently $3L\le n+2(g+1)$. The verifier tests the remaining residues $r\bmod2n$ for every allowed $n,h$. It generates 270,710 residue tests, including 4,701 with $G>1$.

The second calculation omits the arc restriction. As $G$ is coprime to $D$ and $1+2^h$ is a signed power of 2 modulo $G$, it divides $\gcd(f/D,(1+2^h)^{2n}-1)$. Of 1,783,660 gap pairs, two order bounds require a full period scan, and neither yields a survivor. An additional direct enumeration tests all 623,156 residues through $n=85$ and finds 160 passing (1), all in the two stated aligned families. These survivors show that (1) by itself permits the aligned cases. [Verifier](check_seven_binary_pairs.py)

## Six positive bits in a common exponent class

If $b\ge8$ divides all differences between the six positive positions, both $n$ and $h$ in the triple-at-1 case would be multiples of $b$, contradicting $h=n\pm1$. Thus 1 belongs to a pair. The [pair-and-triple proof](seven_large_bases.md#the-remaining-allocation) now works at $B=2^b\ge256$. Here $U\ge3$ by the endpoint reduction, $A,C$ are odd because $m$ is odd, and $v_2(bd)=v_2(u)$ still follows from $U>2$. All numerical estimates and order arguments use only these properties, subject to two small points.

First, if all three terms of $y$ coincide modulo $d$ with the same sign, then $U\le3$. At any root of $X^d+1$ these terms have sum of modulus 3, so $h_A=1$, $A\le4^d$. The two trinomials have distinct actual positions because the six binary bits are distinct, giving $C\le3^H$. Thus $U^3A^2C^2\le27\cdot144^H<27B^{2H}$. Opposite signs in a collision still give $U=1$.

Second, the exceptional degree $d=3t$ needs the size bound at small $H$ as well. Here $A$ divides either $3^t+1$ or $3^t-1$ and is odd. Hence $A\le(3^t+1)/2\le3^t$. The three positions of $x+z$ are distinct, so $C\le3^H$ and $(AC)^2\le(9\sqrt3)^H<16^H$ from $4t\le H$. All other degree cases of the proof hold for every $H\ge1$. This completes the binary theorem.
