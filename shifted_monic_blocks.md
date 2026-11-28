# A quadratic block with an independently placed leading shift

Let $b,c\ge0$ be integers, $H=b+c+3$, $G(X)=X^2+bX+c$, and $p$ a prime. Write $X=p^a$ and $Q=p^e$, where $a,e\ge1$. Each of the following conditions proves Wasserman's conjecture for **every three positive rows** of weight
$$N=1+QG(X).$$

1. If $G$ is irreducible over $\mathbb Q$ or is a square, $X>H^{15}$ suffices for **every** outer exponent $e$. For $0\le b,c\le10$ under this factorization hypothesis, both exponents and the prime are unrestricted. This gives 100 quadratics, including $X^2+10X+10$.
2. For **every** quadratic $G$ above and **every** prime $p$, the conjecture holds whenever
   $$a\ge212\log_p H,\qquad 1\le e\le \frac{a^2}{24\log_p H}.$$
   A slightly larger range is $12e\le a(a\log_H p-106)$. The two exponents may be coprime. For instance $1+p^e(p^a+6)(p^a+35)$ works for every prime when $a\ge1696$ and $1\le e\le a^2/192$, since $H=254<2^8$.
3. If $g=\gcd(a,e)$, it also suffices that $p^g>H^{15+4e/a}$. In particular every shifted monic cubic $P^3+bP^2+cP+1$ satisfies the conjecture at every prime power $P>H^{19}$, including cases in which both polynomials in the [previous resultant criterion](polynomial_prime_powers.md) are reducible.

The bound behind condition 2 can be applied directly at other exponents. For **any** integers $j\ge1$, $\ell\ge0$, it suffices that
$$p^a\ge H^{33+2\ell+6j}p^{2|je-\ell a|}. \tag{*}$$
There is an arbitrary-degree version of condition 1 below.

## Normalizing the row that avoids the base prime

By the [endpoint reduction](all_weights_endpoint_reduction.md), a counterexample supplies two rows $A,C$ between them avoiding every prime of $N(N-1)$. Assign $p$ to $A$. Write $UV=N$ and $QBT=N-1$ with $U,QB$ assigned to $A$ and $V,T$ to $C$, as products of full prime powers. In particular $BT=G(X)$. If $p\mid G(X)$, its entire contribution to $G(X)$ is in $B$. The normalized entries of $C$ after division by $V$ sum to $U$. Each full power in $T$ divides two of those entries. Their product and AM–GM give
$$U\ge3,\qquad 27G(X)^2\le U^3 B^2.\tag{1}$$

For $X>\max(b,c)$, Kummer's criterion writes $A_i=Q L_i(X)+\delta_i$, with $\delta=(0,0,1)$. The three nonnegative integral polynomials $L_i$ sum coefficientwise to $G$. The digits within each coefficient block also split without carries. All $L_i$ are nonzero because each actual part is a positive multiple of $U$. One polynomial is monic quadratic and the other two have degree at most one. The two unshifted polynomials $L_1,L_2$ both have values divisible by $U$, since $U$ is coprime to $Q$. One of them has degree at most one, giving $U\le HX$.

At a prime $q$ assigned to $B$, the residues in $A$ at the full $q$-power in $QG(X)$ are $0,0,1$ in some order. If the 1 occurs at slot $k$, every unchanged slot $i$ with $\delta_i=[i=k]$ has $L_i(X)\equiv0$ modulo the full $q$-power in $G(X)$. This cancellation also works at $q=p$, where $k=3$ and the larger modulus contains the full factor $Q$. In particular, every full power of $B$ divides at least one of $L_1(X),L_2(X)$.

## Irreducible quadratics, squares and higher degrees

For each unshifted polynomial use $L_i$ if it has degree at most one, and $G-L_i$ otherwise. This gives a nonzero line $\alpha X+\beta$ with $0\le\alpha\le b$ and $0\le\beta\le c$, with the same divisibility at factors of $G(X)$. For $\alpha=0$ use the bound $\beta$. For $\alpha>0$, any common divisor of $G(X)$ and the line divides
$$R=\beta^2-b\alpha\beta+c\alpha^2.$$
If $R\ne0$, its absolute value is at most $H^3$. If $G$ is irreducible, $R$ cannot vanish. If $G=(X+h)^2$ and $R=0$, the line is $\alpha(X+h)$. For $q^v\parallel X+h$, the condition $q^{2v}\mid\alpha(X+h)$ implies $q^v\mid\alpha$. Thus the product of the full powers assigned to that line divides $\alpha^2$. These arguments bound $B$ by $H^6$ on taking the product of the two line bounds. Now (1) gives $27X^4\le H^{15}X^3$, proving the first threshold.

More generally, take a monic nonnegative polynomial $G=X^n+J$ with $J\ne0$ of degree $r<n$. Factor $G$ over $\mathbb Q$ into full powers $f^k$ of distinct monic irreducibles. Let $D$ be the sum of the degrees of the full powers with $\deg(f^k)>r$. **If $2D>3r$, then $1+p^eG(p^a)$ satisfies the conjecture for every $e$ at every sufficiently large $p^a$, with an effective threshold depending only on $G$.** For instance $X^5+2X^3+6X+10$ is Eisenstein at 2, and $2D=10>9$.

To prove this, split the coefficient blocks with $X$ larger than all coefficients and replace the top-degree part by its complement as above for predecessor divisibility. Each of the resulting two unshifted tests is a nonzero polynomial $L$ of degree at most $r$. For a counted $f^k$, write $L=f^t R$ with $t<k$. If $q^v\parallel f(X)$ and $q$ is assigned to this test, then $q^{kv}\mid L(X)$ and $q^v\mid R(X)$. The assigned portion of $f(X)$ divides $|\operatorname{Res}(f,R)|$, a fixed nonzero integer. This holds even at primes shared by different polynomial factors after specialization. Using both tests, the portion of the entire degree-$D$ product at primes of $B$ is bounded. Its remaining portion divides $T$, so $T\gg_G X^D$. At least one actual unshifted polynomial itself has degree at most $r$, giving $U\ll_G X^r$. Equation (1) excludes $2D>3r$. All constants come from finitely many coefficient splits, polynomial evaluations and integer resultants.

## Comparing the exponents at a linear factor

Write $G=(X+h)(X+k)$, $0\le h<k$. Choose integers $j\ge1$, $\ell\ge0$, and put $v=je-\ell a$, $v_+=\max(v,0)$, $v_-=\max(-v,0)$. Thus $p^{v_-}Q^j=p^{v_+}X^\ell$. Set $W=p^{|v|}H^{\ell+3j}$ and assume $X>H$. For a low-degree unshifted part $\alpha X+\beta$ with $\alpha>0$, the congruence $QG(X)\equiv-1\pmod U$ gives
$$
 U\mid p^{v_+}(-\beta)^\ell (\beta^2-b\alpha\beta+c\alpha^2)^j
           -(-1)^j p^{v_-}\alpha^{\ell+2j}.                         \tag{2}
$$
Multiply by $\alpha^{\ell+2j}$ and use $\alpha X\equiv-\beta$. If nonzero, the integer on the right has absolute value at most $2W$. A constant part already bounds $U$ by $H$.

For either $d\in\{h,k\}$, let $B_d$ be the portion of $X+d$ supported on primes of $B$. Thus $B=B_hB_k$, including at primes dividing both evaluated factors. Split $B_d$ according to the exceptional slot at each prime. The group with unchanged exceptional slot divides every $L_i(-d)$. For a switch to slot $i=1$ or 2, let $i'$ denote the other unshifted slot. The group divides both $L_{i'}(-d)$ and
$$p^{v_+}(-d)^\ell L_i(-d)^j-p^{v_-},                                 \tag{3}$$
since $QL_i(X)\equiv1$ at its assigned powers. We use the convention $x^0=1$. The evaluations of the $L_i$ have absolute value at most $H^3$, and (3) has absolute value at most $2W$.

When the $L_i(-d)$ do not all vanish, at least one of the two unshifted evaluations is nonzero, since their sum with the third is $G(-d)=0$. Thus at most one switched slot needs (3) to bound its group. If that integer is nonzero whenever the opposite evaluation vanishes, multiplying the three group bounds gives $B_d\le2H^6W$. When **both** unshifted evaluations are nonzero, the stronger $B_d\le H^9$ follows directly. If (3) vanishes, then $\ell>0$ (at $\ell=0$ one has $v=je>0$), and $d$ must be a power of $p$, including 1.

**At least one of $B_h,B_k$ is at most $2H^6W$.** If neither root annihilates all three parts, failure of the preceding test at both roots would make both parameters powers of $p$. Then $G(X)$ has at most four base-$p$ digit units. Splitting them among three nonempty parts gives at least two singleton parts. One is unshifted, forcing $U$ to divide a power of $p$.

Now suppose $X+d$ divides all $L_i$ as polynomials. The other linear factor $X+d'$ cannot do so as well. For $d=0$, if $d'$ were a power of $p$, $G(X)=X(X+d')$ would have only two digit units. For $d>0$ the three quotients are two positive integral constants and one monic line of nonnegative constant term, by evaluation at zero. Their sum is $X+d'$, so $d'>0$. If $d'=p^t$, the full $p$-part of the constant block $dd'$ divides each part of that block by carry freedom. Thus $d'$ divides both positive constant quotients, whose sum is at most $d'$, an impossibility. In every remaining case $d'>0$ is not a power of $p$, so (3) is nonzero there. This proves the assertion. In the case $d=0$, $d'>0$ already follows from distinctness.

## Why the two integer bounds cannot fail together

Assume $X\ge H^{33}W^2$ as in (*). If $U\le H^3$, use the bound just proved at one factor and $X+d\le2X$ at the other. Equation (1) gives $27X^2\le16H^{21}W^2$, a contradiction.

For $U>H^3$, an unshifted part of degree at most one must be a line $L=\alpha X+\beta$ with $\alpha,\beta>0$. Indeed $U$ is coprime to $X$. Put $r=\beta/\alpha$. The other unshifted polynomial $M$ satisfies $M(-r)=0$, since otherwise $U$ divides the nonzero integer $\alpha^{\deg M}M(-r)$ of absolute value at most $H^3$. Also $r$ differs from $h,k$. Equality would give $U\mid\alpha$ by coprimality with $G(X)$. Hence at least one of $-h,-k$ is a zero of **neither** unshifted part. Only $M$ could have another root besides $-r$, and it has degree at most two. At this parameter we have $B_d\le H^9$.

If (2) is nonzero, then $U\le2W$ and $B\le2H^9 X$. Equation (1) gives $27X^2\le32H^{18}W^3$, contradicting (*). If (2) is zero, we claim that the other parameter $d'$ satisfies $B_{d'}\le2H^6W$. The parts cannot all vanish at $-d'$, since $L(-d')\ne0$. The only remaining obstruction would require both (3) at a switched slot and the opposite unshifted evaluation to vanish.

In this event $d'$ is a power of $p$, $\ell>0$, and (2) gives the rational identity
$$p^v(-r)^\ell G(-r)^j=(-1)^j.$$
At a denominator prime of $r$, monicity gives $v_q(G(-r))=2v_q(r)<0$. At a numerator prime, $G(-r)$ is integral locally. Thus all numerator and denominator primes of $r$ equal $p$, and $|G(-r)|$ is a power of $p$ with integer exponent. The unshifted polynomial vanishing at $-d'$ must also vanish at $-r$. It is necessarily the monic quadratic $(X+r)(X+d')$, so $r$ is a positive integer. Let $d''$ be the other root parameter. Subtracting this quadratic from $G$ leaves
$$(d''-r)(X+d').$$
This polynomial is the sum of the other two nonzero parts, so $d''>r$. Its positive multiplier $d''-r$ divides $|G(-r)|$, hence is a power of $p$. Both its linear and constant coefficients therefore consist of a single base-$p$ digit unit. The remaining unshifted line $\alpha(X+r)$ is positive in **both** coefficient blocks. By carry freedom it must take both units in full, leaving the shifted polynomial zero. This contradiction proves the claim. [Lean verifies](ShiftedQuadraticRoots.lean) the recovery of the monic quadratic from its two roots and the final two-unit contradiction directly at the moduli $(d''-r)d'$ and $(d''-r)X$. The first divides $X$ because it is a power of $p$ smaller than $X$.

We now have $B\le2H^{15}W$ and $U\le HX$, giving $27X\le4H^{33}W^2$, again impossible. This finishes (*). It applies to irreducible quadratics and squares by condition 1 as well.

**Quadratic range of outer exponents.** Among $0,e,2e,3e$ on the circle modulo $a$, either two residues coincide or a circular gap is at most $a/4$. Taking the corresponding difference gives $1\le j\le3$ and $0\le\ell\le3e/a+1$ with $|je-\ell a|\le a/4$. Consequently (*) follows from $(53+6e/a)\log_p H\le a/2$. This gives both forms of condition 2.

**Exact exponent ratio.** Take $j=a/g$ and $\ell=e/g$. If $c=0$, the low unshifted line is $\alpha X$ and $U\mid\alpha$. With $c>0$, a zero in (2) would require a rational root $-r$ of the monic polynomial $Y^\ell G(Y)^j-(-1)^j$ with constant term $\pm1$. The only possible nonpositive rational root is $-1$, requiring $|G(-1)|=1$. For distinct positive $h,k$, the product $(h-1)(k-1)$ is zero or at least two. Thus $U\le2W$ in either case. Using $B\le4H^6WX$ gives $27X^2\le128H^{12}W^5$. Under condition 3, $X>H^{4\ell+15j}$ and $W=H^{\ell+3j}$, contradicting that inequality. This condition implies the first threshold in the irreducible and square cases.

## Verification at every exponent in the coefficient box

When $0\le b,c\le10$, exactly 100 choices give a square or irreducible quadratic. The two checkers cover every inner power below the symbolic cutoff while leaving the outer exponent free.

The first checker bounds $U$ by the numerical gcd of the two unshifted polynomials using the [polynomial specialization routines](polynomial_prime_powers.md). It bounds $B$ by their full-power resultant bounds, enumerates its possible divisors, and uses polynomial translation to obtain a cutoff in $X$ for each divisor. It finds roots of $G(X)$ modulo $B$ by square roots modulo $4B$. The second checker uses the linear remainders in the proof, elementary linear congruences and CRT instead. For its cutoff it uses $27X^4\le C^3B^2$ when the two unshifted parts have gcd bounded by a constant $C$, and $27X\le(\alpha+1)^3B^2$ whenever the low line is $\alpha X+\beta$ with $X>\beta$. For overlapping blocks both checkers enumerate the carry-free splits of $G(X)$ directly.

For each numerical split $(x,y,z)$, $U\mid\gcd(x,y)$ and $U$ is coprime to $pG(X)$. The [shifted denominator theorem](shifted_directions.md) excludes $G(X)/\gcd(x,y,z)<364$ and any denominator with a full power at least one third its size. After (1), each full power $w$ of $T$ requires two multiples of $w$ in the normalized companion of sum $U$ and the third entry congruent to $U\pmod w$. One implementation enumerates exceptional slots by CRT. The other enumerates the individual admissible entries and tests their sums. For each remaining companion, a single prime $q\ne p$ rules out every residue of $Q=p^e$ in its finite power orbit modulo $Uq^m$ at the tested depth. The row assigned $q$ must have no carry at every layer. Also the sum of the entries of any subset divided by their gcd divides its multinomial coefficient, by the integer identities obtained by multiplying that coefficient by each entry divided by the sum. This gives fixed divisors for the two unshifted parts and all cuts of the normalized companion.

Both programs check 90,805 polynomial shapes. The first checks 19,287 overlapping-block rows and excludes 269 companion occurrences by local tests. The independent checker directly tests 184 companion occurrences with sorted entries. These counts can differ through permutations and preliminary tests. It also checks 70,803 split-case determinant systems, including 12,687 nontrivial switched-slot divisibilities. Direct binomial products verify 2,313 carry-free coefficients and every endpoint pair at 12 small weights. A further 3,835 assignments with two prescribed roots are all rejected by direct factorial floors, independently of the last two-unit argument.

~~~sh
python3 check_shifted_monic_blocks.py
python3 check_shifted_monic_blocks_cross.py
lean ShiftedQuadraticRoots.lean
~~~

The Python programs use python-flint and SymPy. The Lean file uses Lean 4.36 and its standard library.
