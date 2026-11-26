# Seven digit units at large predecessor primes

**Theorem.** For every prime $p\ge257$, Wasserman's three-row conjecture holds at every weight $N$ with $p\mid N-1$ and $s_p(N)=7$. There are no restrictions on the six positive exponents or their multiplicities.

More generally let $p$ be odd and $a,b$ positive. For arbitrary nonnegative integers $h_1,\ldots,h_6$ the conjecture holds at
$$N=1+\sum_{i=1}^6 p^{a+b h_i}\qquad\text{if }B=p^b\ge256.$$
Thus any six positive powers of 3 in one exponent class modulo 6, for example, give a family after adding 1. At any odd predecessor prime, including 3, an additional conclusion holds without a restriction on $B$. In a counterexample, a row avoiding $p$ must group the units-place 1 into a **pair**.

The [refined cutoff](seven_refined_cutoffs.md) subsequently lowers the bound to $B\ge64$, including every predecessor prime $p\ge67$.

## Locating the units-place digit

A counterexample contains the endpoint row. Assign $p$ to a nonendpoint row avoiding it. Write $U$ for its assigned weight factor and $S$ for the product of its assigned full predecessor powers after removing the full power of $p$. The [seven-unit reduction](seven_digit_units.md) gives $U\ge7$, both $U,S$ odd, and $27m^2\le U^3S^2$, where $m=(N-1)/p^{v_p(N-1)}$. If there is a carry among the six specified positive powers, then $s_p(N)\le5$ and the [five-unit theorem](six_digit_gaps.md) suffices. Otherwise the avoiding row partitions the seven units into two pairs and a triple.

Suppose the 1 belongs to the triple, with paired parts $p^i(1+p^d),p^j(1+p^e)$ and triple $1+p^r+p^s$. Put $g=\gcd(d,e)$ and order the pairs so $d\ge e$. Both gap quotients by $g$ are odd, and $U\mid p^g+1$. Split $S$ into $S_0,S_d,S_e$ according as the two zero residues modulo each assigned power occur on both pairs, on the first pair and the triple, or on the second pair and the triple. Thus $S_0\mid p^g+1$ and $S_0\mid p^{|r-s|}+1$. The [binomial–trinomial gcd inequality](seven_pair_gaps.md) gives $U^3 S_0^2<27p^{2g}$, so
$$p^d\le m<p^g S_d S_e. \tag{1}$$

Here $S_d$ divides both $p^d+1$ and $p^j(1+p^e)-1$, and $S_e$ has the interchanged divisibilities. We recall the [resultant bound for this pair](equal_predecessor_gaps.md#root-of-unity-products). At a root $\zeta$ of $X^d+1$, $\zeta^e$ has even order, so $|1+\zeta^e|\ne1$. The resultant with $X^j(1+X^e)-1$ is nonzero. If $0,j,j+e$ are distinct modulo $d$, orthogonality gives mean squared modulus 3 at the $d$ roots, hence $S_d\le3^{d/2}$. If $d\mid e$ then $S_d=1$. In the other collision cases the equations force $p^{\pm e}=-2\bmod S_d$ or $S_d=1$. Raising to the odd power $d/g$ gives $S_d\le 2^{d/g}-1$. In particular $S_d\le 2^d$ always, and $S_d\le3^{d/2}$ when $g\ge2$. The same bounds hold with $e$ in place of $d$.

For equal gaps both factors are 1, contradicting (1). Otherwise $d-e\ge2g$. For $p\ge5$, $S_d S_e\le2^{d+e}\le4^{d-g}<p^{d-g}$. For $p=3$ we need $g\ge2$ since $U\ge7$ divides $3^g+1$. Then $S_d S_e\le3^{(d+e)/2}\le3^{d-g}$. This contradicts (1) in every odd base.

## The remaining allocation

Translate the minimum $h_i$ to zero by increasing $a$. Set $H=\max h_i$, $W=p^a$ and $m=\sum B^{h_i}$. The [21-position theorem](seven_progressions.md) covers $H\le20$, so take $H>20$. The remaining parts have the form
$$W x(B),\quad W y(B),\quad 1+W z(B),$$
where $x=X^i(1+X^d)$, $y$ is a sum of three monomials, and $z$ is one monomial. Neither $x$ nor $y$ is supported on a single position. Write $z=X^v$, $u=a+bv$, $F=X^d+1$. The two actual pair gaps are $bd,u$. Their 2-adic valuations agree since both powers of $p$ are $-1\bmod U$ and $U>2$. Every full power of $S$ has a zero on the shifted part: otherwise it would divide $p^u$ by the residue 1 there. Write $A,C$ for the respective products of powers with their other zero on $Wx(B),Wy(B)$. Then
$$U\mid\gcd(F(B),y(B)),\quad A\mid\gcd(F(B),(y+z)(B)),\quad
C\mid\gcd(y(B),(x+z)(B)),\quad C\mid p^u+1. \tag{2}$$
We use $27B^{2H}\le27m^2\le U^3 A^2 C^2$. All three integers in (2) are odd and prime to $B$.

## Bounding the factor from two trinomials

Strip the lowest monomial individually from $y$ and $x+z$, and denote the resulting polynomials by $P,Q$. If $P=Q$, then $U$ divides both $x(B)$ and $(x+z)(B)$, hence the unit $z(B)$, impossible. Otherwise **every common root of $P,Q$ is a root of unity**.

Indeed a trinomial $f=1+X^j+X^n$ with $0<j<n$ has at most one irreducible noncyclotomic factor, with multiplicity one, by Ljunggren's trinomial theorem. If two such polynomials share this factor, their quotients by their reciprocals must coincide, since all remaining factors are cyclotomic and reciprocal (the positive sums do not vanish at 1). Write $f^*=X^n f(1/X)$. A reciprocal $f$ has $n=2j$ and is itself cyclotomic in the sense that all its roots are roots of unity. Otherwise the first nonzero coefficient of $f/f^*-1$ has degree $\min(j,n-j)$ and sign specifying which side is smaller. Reversing both polynomials if necessary, this determines $j<n/2$ in both. At $X=2$ the ratio $f/f^*=(1+2^j+2^j t)/(1+(1+2^j)t)$, with $t=2^{n-j}$, is strictly decreasing in $t$. Hence $n$ agrees as well. A repeated position gives either $X^n+2$ or its reciprocal $2X^n+1$, irreducible by Eisenstein and reversal. Their roots are not algebraic units, whereas roots of the monic trinomials with constant term 1 are units. This also handles two repeated-position inputs. [Ljunggren, Theorem 3](https://doi.org/10.7146/math.scand.a-10593)

Let $h$ be the positive-leading primitive polynomial gcd of $P,Q$. Each of its roots $\zeta$ makes the three unit terms of $x+z$ sum to zero. They form an equilateral triangle, so $\zeta^d$ is a primitive cube root. These roots are simple in $Q$: at three distinct positions the two terms of the logarithmic derivative after normalization have distinct, nonopposite directions. A repeated position has no roots on the unit circle. Thus $h\mid X^{3d}-1$. For any odd prime $\ell\mid C$, the order of $p\bmod\ell$ has 2-adic valuation $v_2(u)+1$ because $\ell\mid p^u+1$. It cannot divide $3bd$. In particular $\gcd(C,h(B))=1$.

We use an integer Bézout bound. For polynomials of degrees $m_1,m_2$ with primitive gcd $h$, the Sylvester columns are the first $m_2$ shifts of the first polynomial and the first $m_1$ shifts of the second. Their rational span contains $h$ by the Euclidean algorithm applied to the cofactors. Choose a column basis and a nonsingular minor on these columns. Cramer's rule places the determinant times $h$ in the integer column span. The absolute determinant is bounded by the product of column norms by Hadamard. The coefficient vectors here have Euclidean norm at most $\sqrt5$, and norm $\sqrt3$ for three distinct positions. Since $C$ is coprime to $h(B)$, this proves
$$ C\le5^H,\qquad C\le 15^{H/2}\ \text{if $y$ has distinct positions},\qquad
C\le3^H\ \text{if both have distinct positions}. \tag{3}$$

## Common roots at the unshifted pair

Put $h_U=\gcd(F,y)$ and $h_A=\gcd(F,y+z)$, both monic. If $F\mid y+z$, then $F(B)\mid m$ and $U=1$ since $U$ divides $N=1+Wm$. Exclude this case. For any polynomial $E$, reduction of multiplication by $E$ modulo $F$ gives a $d$-dimensional signed shift matrix. Its rank is $d-k$ for $k=\deg\gcd(F,E)$. Its rational image contains that gcd by Bézout. The same minor argument gives
$$\gcd(F(B),E(B))\le L^{d-k}\gcd(F,E)(B), \tag{4}$$
where $L$ is the norm of the reduced coefficient vector of $E$. All column norms equal $L$. If $k=d$, the gcd equals $F$ and the estimate follows directly.

Write $k=\deg h_A$. Four unit complex numbers summing to zero have an opposite pairing, since the two endpoints on the unit circle with any fixed nonzero sum are unique up to order. At a root of $F$, a pair of exponents with difference $w$ is opposite only if $v_2(w)=v_2(d)$. For a matching with two nonzero differences $w_1,w_2$, the common roots are exactly the roots of $X^j+1$ for $j=\gcd(d,|w_1|,|w_2|)$ if the three quotients by $j$ are odd, and otherwise there are none. This follows by Bézout, as for the integer pair gaps. Each nonempty matching has $j\le d/3$, since $j=d$ would give $F\mid y+z$. There are at most two admissible matchings: the allowed differences join opposite residue classes modulo $2^{v_2(d)+1}$, and the resulting bipartite graph on four labels has at most two perfect matchings. If two give distinct $j=d/n_1,d/n_2$, their union has size $d(1/n_1+1/n_2-1/\operatorname{lcm}(n_1,n_2))$. The odd integers $3\le n_1<n_2$ give at most $7d/15$, using $\operatorname{lcm}(n_1,n_2)\le n_1n_2$. Hence $k\le7d/15$.

If two terms of $y$ occupy the same residue modulo $d$ with opposite signs after reduction, then $U$ divides the remaining monomial, impossible. If all occupy the same residue with the same sign then $U\le3$, also impossible. In the remaining collision case two terms coincide modulo $2d$. They cannot match one another, and interchanging them gives the very same root conditions. Thus $k\le d/3$ in this case.

## The size inequality

First suppose $h_U=1$. For three distinct residues of $y$ modulo $d$, (4) gives $U\le3^{d/2}$ and $A\le6^{(d-k)/2}(B+1)^k$. Indeed adding $z$ to the three signed positions gives squared norm at most 6, and each root of $h_A$ has modulus 1. With a collision the corresponding bounds are $U\le5^{d/2}$ and $A\le10^{(d-k)/2}(B+1)^k$. Combining these with (3), using $B\ge256$ and $d\le H$, bounds $(U^3 A^2 C^2/B^{2H})^{1/H}$ respectively by
$$
\frac{3^{3/2}6^{8/15}15(257/256)^{14/15}}{256^{16/15}},
\qquad
\frac{5^{3/2}10^{2/3}25(257/256)^{2/3}}{256^{4/3}}.
$$
Both are less than 1 by integer powers. This contradicts the required size.

If $h_U\ne1$, the [trinomial root classification](seven_pair_gaps.md#three-distinct-positions-modulo-n) gives $h_U=X^{2t}-X^t+1$, where $d/t=3\bmod6$ and the three positions of $y$ differ from its least position by the residue classes $0,2t,4t\bmod 6t$. In particular $H\ge4t$. All three residues are distinct modulo $d$. Every difference is divisible by $2t$, whereas $v_2(d)=v_2(t)$. No two terms of $y$ can be opposite at any root of $F$. Hence $h_A=1$ and $A\le6^{d/2}$. If $d\ge9t$, (4) gives $U\le3^{d/2}B^{2t}$. Using (3), the ratio above is at most $3^{3/2}\cdot6\cdot15/256^{4/3}<1$.

It remains that $d=3t$. Then $F=(X^t+1)h_U$ and $U\le3h_U(B)<3B^{2t}$ because the three terms of $y$ agree modulo $B^t+1$ after removing their common lowest power. Also $A$ is coprime to $h_U(B)$ since $h_U\mid y$ and $z$ is a monomial. Thus $A\mid B^t+1$. Modulo $A$, $y(B)+z(B)=0$ gives $3B^j=-B^v$ for the least position $j$ of $y$. Raising to the $t$-th power yields $A\le3^t+1$.

We claim $(AC)^2\le16^H$. When $x+z$ has three distinct positions use $C\le3^H$. For $t\ge3$ the bound $(3^t+1)^{1/t}\le28^{1/3}<(16/9)^2$ and $4t\le H$ suffice. For $t=1,2$ use $100<(16/9)^{21}$ and $H>20$. If $x+z$ has a repeated position, its stripped form is $2+X^d$ or $1+2X^d$. In the variable $T=X^t$ it has degree three, while the stripped $y$ has three distinct terms and degree at most $H/t$. These polynomials are coprime by the repeated-position argument above, with the actual degrees now taken in $T$. The Sylvester bound at $T=B^t$ gives $C^2\le 27\cdot5^{H/t}$. As $(3^t+1)^2\le 2^H$, we get $(AC)^2\le27\cdot10^H<16^H$ for $H>20$. Finally $2H-6t\ge H/2$, so $16^H\le B^{2H-6t}$. With $U<3B^{2t}$ this contradicts $27B^{2H}\le U^3A^2C^2$.

Run `python3 check_seven_large_bases.py` with python-flint. It checks 246,051 distinct trinomial pairs, 31,920 four-term gcds and the cutoff inequalities using integer arithmetic. A separate digit enumeration checks 69,616 positive carry-free rows at 808 weight/base pairs, including outer shifts outside the exponent progression for the smaller primes. [Verifier](check_seven_large_bases.py)

## An optimal prime list with six distinct positive exponents

Take $p=257$ and $N=1+p+p^2+p^3+p^7+p^{17}+p^{23}$. Four primes certify every three-row selection:
$$p,\qquad q=N/3,\qquad s=(N-1)/(1598p),\qquad r=N-44.$$
A $p$-avoiding row has two positive parts with sum at most $N-p^{23}<\min(s,2q)$. It cannot avoid $q$, which requires three multiples of $q$, or $s$, which requires two multiples of $s$. Its two smallest parts also sum to at least $1+p$, larger than $N-r$, so it carries at $r$. Avoiding both $q$ and $s$ would require $N\ge q(2s+1)$. Avoiding $r>N/2$ forces two small parts of sum at most 44, incompatible with either $q$ or $s$. Hence each row avoids at most one listed prime.

The list has minimum cardinality. Since $N=3q$ with $q\ne3$ and $N$ is odd, its digit sum at an odd prime is odd and exceeds one. Its binary digit sum is at least three because $p\mid N-1$. Every prime can therefore be avoided by a positive row, by allocating three or more digit units among the parts. Three chosen primes can be avoided one per row. The verifier supplies recursive primality certificates and checks all 301 $p$-avoiding rows against the three other primes by Legendre's formula. [Four-prime verifier](check_seven_prime_list.py)
