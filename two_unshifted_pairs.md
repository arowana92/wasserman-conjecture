# Two binomial parts in a candidate row, at any digit sum

Let $p\mid N-1$ be prime. In a counterexample suppose a nonendpoint row avoids $p$ and has two parts
$$p^i(1+p^d),\qquad p^j(1+p^e),\qquad i,j\ge1.$$
Then **$d=e>0$**. This conclusion places no upper bound on the digit sum of $N$ or on the number of base-$p$ units in the third part. If $F=p^d+1$, $a=v_p(N-1)$, $m=(N-1)/p^a$, then the total product of the full powers of $N(N-1)$ at all primes this row avoids is at most $p^a F<(N-1)/2$.

In particular at most one of the two nonendpoint rows can have two such parts, even if one uses a different prime for the other row. The two rows together must cover all the full powers of $N(N-1)$ by the endpoint reduction.

## The unequal-gap proof with an arbitrary third part

Assign all avoidable endpoint prime powers to this row, and write $U$ for the factor from $N$. Carry avoidance at $p$ gives $i,j\ge a$ and identifies the third part as 1 plus further positive-position units. The third part exceeds 1 since $U\ge3$ divides it. Each gap is positive because otherwise $U$ would divide $2p^i$ or $2p^j$ while coprime to $p$.

Order $d\ge e$, put $g=\gcd(d,e)$ and $G=p^g+1$. The two quotients by $g$ are odd and $U\mid G$. Let $S_0,S_d,S_e$ denote the assigned factors of $m$ with both pairs zero, only the first pair zero, and only the second pair zero, respectively. At an assigned full power precisely two actual parts are zero and the third is 1. This gives
$$US_0\mid G,\quad S_d\mid (p^d+1)/G,\quad S_d\mid p^j(1+p^e)-1,$$
and the two interchanged divisibilities for $S_e$. In particular $S_d,S_e$ are coprime to $G$. At odd $p$ the paired parts have the same 2-adic valuation, so this row carries at 2 and $US_0\mid G/2$. The [companion product bound](endpoint_reduction.md) yields
$$p^d\le m\le(U/3)^{3/2}S_0 S_d S_e
\le (G/(3c))^{3/2} S_d S_e,\qquad
c=\begin{cases}2,&p\text{ odd},\\1,&p=2.\end{cases}$$

These inequalities and the bin divisibilities are exactly the inputs of the unequal-gap proof at six units. At $p\ge5$ its two resultant estimates depend only on $d,e$ and the powers outside the binomials. At $p=3$ the cases $g=1$, $2\le g\le4$ and $g\ge5$ use just the same inputs, including the factor of two in $G/2$. At $p=2$ the only extra spacing input for the finite calculation is $d<3e+g$ when $d>e$. The rotation argument in the binary seven-unit proof derives it from the displayed inequalities alone. Thus the unequal case is excluded for every base. [Odd-base estimates and finite binary exclusion](equal_predecessor_gaps.md#root-of-unity-products) · [Binary spacing deduction](seven_binary_alignment.md#excluding-unequal-unshifted-gaps)

For equal gaps, $S_d=S_e=1$ and $US_0\mid F$. Both unshifted parts contribute to $m$, and the third contributes a positive amount beyond its units digit. It follows that
$$m>F(p^{i-a}+p^{j-a})\ge2F.$$
Multiplying by $p^a$ proves the claimed bound on all full powers avoided by the row. If $D$ is the span of all positive base-$p$ digit positions of $N$, then $p^D\le m\le(F/3)^{3/2}<p^{3d/2}$. Thus each equal gap uses more than two thirds of the entire span, at every digit sum.

## The companion needs at least three further prime powers

Put $T=m/S_0$. The three positive companion parts normalized by their assigned weight scale sum to $U$. Each full prime power of $T$ has one exceptional position, at which the normalized part is $U$ modulo that power, and the other two parts are zero modulo it. Group the powers by exceptional position, with products $t_1,t_2,t_3$. The part in position $k$ is at least $T/t_k$. Thus
$$\frac1{t_1}+\frac1{t_2}+\frac1{t_3}\le\frac{U}{T}
\le\frac{F}{m}<\frac1{p^{i-a}+p^{j-a}}.$$
In particular $N-1$ needs at least four distinct prime factors, including $p$. If it has exactly four, the three other full prime powers are the three bins and must satisfy the displayed strict inequality. At odd $p$, this rules out $2\parallel N-1$ in that case. At $p=2$ the two lower bit positions differ, so $2^{i-a}+2^{j-a}\ge3$; in particular a full power equal to 3 rules out the case.

For example, the full conjecture holds at every odd six-bit weight with $N-1=2^a\cdot3\cdot q^v r^w$ for distinct primes $q,r>3$ and positive $a,v,w$. Every possible nonendpoint row avoiding 2 at six units has two unshifted pairs, so the exclusion applies to the weight itself.

## The endpoint divisor of a single coefficient

For arbitrary positive parts $x,y,z$ of sum $N$, put $K=N!/(x!y!z!)$ and $H_K=N(N-1)/\gcd(N(N-1),K)$. An unconditional test using only the two chosen parts is
$$H_K\mid L(x,y),\qquad L(x,y)=\gcd(x,y)\gcd(x,y-1)\gcd(x-1,y).$$
Indeed $Kxy/[N(N-1)]$ is the trinomial coefficient after deleting one unit from each chosen part. Deleting two units from the same part gives the same integrality for $x(x-1)$ and $y(y-1)$, including zero if that part is one. The gcd of these three numerators before multiplication by $K$ equals $L(x,y)$. Prime by prime, a divisor of $xy$ divides at least one of $x,y$. If it divides both, the gcd takes their minimum valuation. A prime dividing $x$ alone contributes the minimum valuation of $x$ and $y-1$, and the other case is symmetric.

The divisibility $H_K\mid L(x,y)$ is also [formalized in Lean](PairDivisor.lean). The proof starts from Pascal’s recurrence, proves the factorial formula for the resulting trinomial coefficient, and then deletes units and proves the gcd divisibility within Lean. It checks with Lean 4.12.0 using `lean PairDivisor.lean`.

In particular, if $x,y$ have the same prime support, then $H_K\mid\gcd(x,y)$. For the equal-gap parts above we obtain $H_K\mid p^aF$, even though their common factor is $p^{\min(i,j)}F$: $p^a$ is the full power of $p$ in $N(N-1)$. This strengthens the bound on avoided prime powers to one that retains partial valuations of $K$ as well.

For any other row of the same weight with chosen parts $v,w$ and coefficient $J$,
$$\frac{N(N-1)}{\gcd(N(N-1),K,J)}
=\operatorname{lcm}(H_K,H_J)\mid\operatorname{lcm}(L(x,y),L(v,w)).$$
For example, if the first pair has identical prime support and one part of the second pair divides the other, the two rows share an endpoint divisor of size at least three. Indeed $L(x,y)\le(N-1)/2$, and $L(v,w)$ divides the larger part when the smaller divides the larger. Each could also have a pair of different, nondividing parts with the same prime support, in which case the shared divisor is at least five.

## A sparse prime dividing the weight

Suppose $q\mid N$ is prime and $s_q(N)\le6$. The row with the two unshifted binomials must itself avoid $q$: otherwise the companion would need to avoid it. The [single-coefficient bound](six_units_weight_prime.md) at $q$ gives $H_J<(N-1)/2$ for such a companion. Since $H_K<(N-1)/2$, their least common multiple is too small to cover $N(N-1)$.

In particular for odd $p$, $q=2$ is impossible since the two binomials have the same 2-adic valuation. For $q$ odd with $s_q(N)\le5$, neither binomial can be a pure $q$-power as both are divisible by $p$. They must therefore both have exactly two $q$-units, with the third part a pure $q$-power. We claim that **$i=j$ is necessary**. For odd $p$ this requires multiplicity at least two at both $i$ and $i+d$. For $p=2$ it is impossible. Thus at every odd weight with a prime $q\mid N$ of digit sum at most six, a binary carry-free row in a counterexample cannot have two even parts with two bits apiece. Here $s_q(N)$ is odd, so the five-unit argument applies.

To prove the claim, suppose $i<j$ and write $p^iF=q^b(1+q^f)$, $p^jF=q^b(1+q^l)$, with $l>f\ge0$. Both valuations equal $b=v_q(F)$ since $q$ is odd. If $f>0$, divisibility by $1+q^f$ gives $l=fT$ for some odd $T\ge3$: the order of $q$ modulo $1+q^f$ is $2f$, since it divides $2f$ and exceeds $f$. The quotient $(1+q^{fT})/(1+q^f)=p^{j-i}$ is odd, excluding $p=2$. At odd $p$, lifting at $q^f=-1\bmod p$ makes its valuation $v_p(T)$. But the quotient is greater than $T$ for $q^f\ge3$ and $T\ge3$, contradicting its being a pure $p$-power.

If $f=0$, $p^iF=2q^b$ requires $p=2$, $i=1$ and $F=q^b$. Now $1+q^l=2^j$. Even $l$ gives $2\bmod8$, and an odd $l>1$ introduces a nontrivial odd factor by the geometric sum. Hence $l=1$ and $q=2^j-1$. Similarly $q^b-1=2^d$ forces $b=1$ or $q^{b/2}-1,q^{b/2}+1$ to be powers of two differing by two. In the former case $2^j-2^d=2$, giving $d=1,j=2$. In the latter case $q=3,b=2,d=3$. Both give $F\le9$, whereas the companion inequality and $m>2F$ require $F>108$.

## Verification


Run `python3 check_two_unshifted_pairs.py`. It checks the two-part divisor on 284,025 pairs of positions using exact multinomial coefficients. For 585,852 carry-free choices of the binomials and third part, 550,399 have digit sum above seven. Of these choices over all digit sums, 23,019 have unequal gaps and a possible scale of at least three, and every one fails the companion size bound. Four equal-gap choices pass that bound. The sparse weight-prime comparison independently recovers the two small binary exceptions in the proof. [Verifier](check_two_unshifted_pairs.py)
