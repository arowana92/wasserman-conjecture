# Limiting values at the six-part product bound

For an endpoint counterexample let $M=N-1$ and $L=Z_A Z_C/(N^3 M^2)$, where $Z_A,Z_C$ are the products of the nonendpoint parts. The [six-part product bound](six_part_product.md) gives $M\ll L^{25}$. This note describes the weights where the two sides have comparable size.

* Along any sequence of counterexamples with $N\to\infty$, every positive accumulation point of $M/L^{25}$ is $1/j$ for an integer $j\ge13$. In particular its upper limit is at most $1/13$. The possible denominators omit a set of natural density greater than 55%. If one imposes only the assigned lowest prime-power residues, the largest possible limiting value is exactly one.
* For any fixed $C>0$ the weights with $L\le C M^{1/25}$ belong, up to a bounded initial set, to finitely many rational polynomial families of degree 25. The number of possible exceptions in this class through $X$ is $O_C(X^{(1-\gamma_C)/25})$ for some $\gamma_C>0$.

## Reconstructing the families

Follow the decreasing-fraction order in the [product proof](six_part_product.md#a-rational-product-close-to-an-integer). We have $D=J=G=1$. At each proper rectangle let $\Delta=K-XY$ denote its positive integer deficit numerator. The estimate for the next bin gains a factor $1/\Delta$, since the missing terms exceed $\Delta/K$. Write $\Delta_0,\Delta_1,\Delta_2,\Delta_3$ for the successive deficits, starting with the single cell $f_1g_1/(s_1t_1)$.

The four steps have weights $8,4,2,1$. If the two sums $a,b$ of weights on the respective sides satisfy $|a-b|<9$, the same proof gives $M\le9^{40}L^{23}$ since all multiplier exponents are then at most 23. The condition $M\ge c L^{25}$, for any fixed $c>0$, bounds $L$ and $N$ in this case. For the two remaining orders, interchange rows so that the two bins on the first side are added first. The upper bound, keeping the factors from the four deficits and $s_1/f_1<9g_1/t_1$, is
$$M<
\frac{9^{40} g_1^{25} f_1^{16} f_2^8 f_3^4 g_2^2 g_3}
 {t_1^9\Delta_0^8\Delta_1^4\Delta_2^2\Delta_3}.
$$
Comparison with $cL^{25}$ bounds each parameter here except $n=g_1$ in terms of $c$ alone. The deficits reconstruct the six bins exactly:
$$
\begin{aligned}
s_1&=(f_1n+\Delta_0)/t_1,&
s_2&=(f_2s_1n+\Delta_1)/\Delta_0,&
s_3&=(f_3s_1s_2n+\Delta_2)/\Delta_1,\\
V&=f_1s_2s_3+f_2s_1s_3+f_3s_1s_2,&
t_2&=(g_2t_1 V+\Delta_3)/\Delta_2,&
t_3&=(g_3t_1t_2 V-1)/\Delta_3 .
\end{aligned}
$$
Indeed adding a bin $v$ on the first side changes $K-XY$ to $v(K-XY)-f S'Y$ for the prior first-side product $S'$, and analogously on the second side. The final deficit numerator is $M-N=-1$.

These are rational polynomials in $n$ of degrees $1,2,4,6,6,12$ in the displayed order. Multiplying the leading coefficients in $M=\prod_i s_i t_i$ and dividing by that of $L^{25}$ gives the reciprocal of
$$
j=t_1^9\Delta_0^8\Delta_1^4\Delta_2^2\Delta_3\,
 f_1^9 f_2^{17} f_3^{21} g_2^{23}g_3^{24}.
$$
There are finitely many such families for any positive lower bound on $M/L^{25}$. A convergent sequence with positive limit has an infinite subsequence in a single family, proving the discreteness assertion. The explicit family with all these constants one attains limit one for the lowest residues.

Every entry of both actual rows and their common sum are nonconstant rational polynomials in $n$ with positive leading coefficient. Split the integer inputs into classes modulo a common denominator of their coefficients. A class supporting integral entries gives integer-coefficient polynomials after substituting its arithmetic progression, since all nonconstant coefficients after substitution are integral and the constant values are integral as well. The [low-digit polynomial estimate](polynomial_row_divisors.md#divisibility-with-a-power-saving-exceptional-set) at base two shows that each row carries outside $O(Y^{1-\gamma})$ parameters through $Y$, for a fixed $\gamma>0$ for the given polynomials. At an exception at least one of the two rows avoids 2. Here $n\le L\le C X^{1/25}$, and summing over the finitely many families and classes proves the count.

## Excluding denominators below thirteen

If $j\le12$ the five bounded numerators, $t_1,\Delta_0,\Delta_1$ must all be one. Write $u=\Delta_2$, $r=\Delta_3$. Thus $u^2r=j$, $n=L$ is even, and, in shorter notation, the bins and scales are
$$
a=n+1,\quad b=na+1,\quad c=nab+u,\quad V=ab+ac+bc,\quad
d=(V+r)/u,\quad e=(Vd-1)/r .
$$
The rows are $U(bc,ac,ab)$ and $V(nde,e,d)$, with $S=abc$ and $T=de$.

If $u=2$, then $a,b$ are odd and $c$ is even. The first two parts of the first row have equal 2-adic valuation, but its bin $c$ requires it to avoid 2. For $u=3$ only $r=1$ is possible. If $n\equiv0\bmod3$ then $V\equiv1\bmod3$, preventing integrality of $d$. In the other two classes $c$ and one of $a,b$ are both multiples of 3, violating bin coprimality.

It remains to consider $u=1$, $1\le r\le12$. The first row always carries at 2 since $a,b,c$ are odd.

- For $r=1$ the [explicit common divisor](six_part_product.md#why-the-exponent-is-attained-by-the-residue-conditions) is $V$.
- For $r=2,4,6,10,12$ we have $1\le v_2(r)\le2$. Since $V$ is odd and $e=V+(V^2-1)/r$ is an integer, $e$ is odd. Also $d$ is odd, so the second row carries at 2.
- For $r=5,8,11$ use $r\equiv-1\bmod3$. At $n\equiv0\bmod3$ we have $3\mid V$ and $e,d$ are nonzero and opposite modulo 3. As $\gcd(e,d)=1$, their sum divides the binomial coefficient of the scaled pair, forcing a carry at a prime assigned to the second row. At $n\equiv1\bmod3$, $b\equiv0$ and $a+c\equiv0$. The coprime $a,c$ similarly force 3 in the first coefficient although $3\mid b$. At $n\equiv2\bmod3$, both $a$ and $d$ are divisible by 3.
- For $r=3,7,9$, every even residue of $n$ forces a carry in the second row modulo $32,512,512$, respectively. This is a finite test of 16, 256 and 256 residues. For each, compute $a,b,c,V,d$ in the displayed formulas and $e=(Vd-1)r^{-1}$ modulo the specified power of two. The sum of the least residues of $Vnde,Ve,Vd$ exceeds their bitwise OR in every case.

This excludes all $j\le12$. Run `python3 check_product_limits.py` with SymPy installed for the modular check and exact symbolic checks on the leading-coefficient formula.

## A positive density of excluded denominators

Suppose $j$ is squarefree. Every parameter in its product formula with exponent at least two must be one. Thus $u=1$ and $r=j$. For $j$ even, the argument for $v_2(r)=1$ above excludes the family. For odd $j$ the three residue classes of $n$ modulo 3 used above exclude every $j$ congruent to 2 modulo 3. In addition the binary test can be made on the last deficit itself, since the parts modulo 512 depend only on $r$ and $n$ modulo 512 for $r$ odd. Exactly 99 of the 256 odd classes of $r$ force a binary carry for every even class of $n$. The checker tests every class directly.

A second disjoint set of exclusions consists of $j=4v$ with $v$ odd and squarefree. Now the formula allows only $u=1$ and $u=2$, with all the earlier parameters fixed at one. The first choice has $v_2(r)=2$. The second has an even bin $c$ on the first row, where the other two parts have equal 2-adic valuation. Both were excluded above.

For the density calculation the squarefree integers have density $6/\pi^2$. Among them the fraction divisible by 2 is $1/3$. Among odd squarefree integers, the fraction congruent to 2 modulo 3 is $3/8$, independent of the odd residue class modulo 512. These formulas follow by inclusion-exclusion over prime squares. At a prime $p$ the proportion of multiples of $p$ among squarefree integers is $1/(p+1)$, the unit classes are equidistributed, and conditions at distinct primes multiply. The density of $4v$ for odd squarefree $v$ is $1/\pi^2$. The union of these specified exclusions therefore has density
$$
\frac6{\pi^2}\left(\frac13+\frac23\left(\frac38+\frac58\frac{99}{256}\right)\right)+\frac1{\pi^2}
=\frac{2799}{512\pi^2}>0.5539 .
$$
