# Adding three powers to a fixed offset

Let $C\ge2$ and let $p$ be prime. Add at most three powers of $p$ to $C$, with arbitrary nonnegative exponents and repetitions allowed. Write $N$ for the result. Throughout the prime-base statements assume $p\mid N(N-1)$, a condition automatic at $p=2$.

**Theorem.** If the conjecture fails and the largest added exponent satisfies $a>90\lceil\log_p(2C+2)\rceil+9$, there are exactly three added powers, with exponents $a\ge b\ge c$. For some positive integer $l$,
$$ C=l(1+p^{b-c}+p^{a-c}),\qquad l,\ lp^{b-c},\ lp^{a-c}\ \text{add without a base-}p\text{ carry}. \tag{1}$$
One of the two rows paired with the endpoint must be $(p^c+l)(p^{a-c},p^{b-c},1)$ up to permutation. In particular $s_p(C)=3s_p(l)$, where $s_p$ is the base-$p$ digit sum.

**All-exponent corollaries.**

- The conjecture holds at **every weight obtained by adding up to three powers of two to 785**. For every prime $p$, it holds for the corresponding translates of 785 satisfying $p\mid N(N-1)$.
- At every offset $2\le C\le63$, every prime $p$, and every exponent choice, a counterexample would require (1) and the indicated row. In binary, this proves the conjecture at all translates of 41 of these 62 offsets. The remaining 21 offsets admit 22 offset/difference pairs altogether. They are tabulated in the [shifted-row verifier](check_shifted_power_rows.py).

## Comparing the other three common divisors

Put $n=\lceil\log_p(2C+2)\rceil$. By the [endpoint theorem](all_weights_endpoint_reduction.md), in a counterexample one can take a nonendpoint row $A$ avoiding $p$. The [three-power allocation and approximation](three_power_translates.md#what-three-very-high-powers-force) in this height range gives three parts $p^{a_i}+c_i$ with $(a_1,a_2,a_3)=(a,b,c)$, $a_3>4n$, and $c_i$ adding freely to $C$. Let $U$ be the product assigned to $A$ from $N$, and split its assigned product from $N-1$ as $S=S_1S_2S_3$ according to the coordinate congruent to one. Write $W_e$ for the prime-to-$p$ part of the divisor at target $e$, with $e=0$ for $U$ and a coordinate unit for each $S_i$. Then
$$27(N-1)^2\le U^3S^2,\qquad W_e\mid p^{a_i}+c_i-e_i\quad (1\le i\le3). \tag{2}$$
The total $p$-part is at most $p^n$, since $N\equiv C\bmod p^{a_3}$. Some $W_e$ therefore satisfies $W_e\ge p^{2a/9-n}>2p^{9n+a/9}$. Fix this target and put $x_i=e_i-c_i$.

For clarity, approximate $a_i/a$ for $i=2,3$ by $v_i/u_i$ with $1\le u_i\le9$, $0\le v_i\le u_i$ and $|v_i a-u_i a_i|\le a/9$. Compare $x_1^{v_i}$ with $p^{v_i a-u_i a_i}x_i^{u_i}$. After clearing the power of $p$ the difference is divisible by $W_e$ and smaller in absolute value than $W_e$. None of the $x_i$ is zero, and $v_i=0$ would give $|x_i|=p^{a_i}>C$. Let $d$ be the least common multiple of the reduced denominators, $d_1=d$ and $d_i=dv_i/u_i$ for $i=2,3$. They have gcd one. As in the [single common power argument](three_power_diagonals.md#the-divisor-of-the-single-common-power), the valuation identities and Bézout give
$$a_i=j_i+d_i g,\quad x_i=p^{j_i}(\sigma t)^{d_i},\quad W_e\mid p^g-\sigma t,\quad W_e\le2p^{a/d}. \tag{3}$$
Here $j_i=v_p(|x_i|)$, $t\ge1$ is prime to $p$, $g=(a-j_1)/d\ge1$ and $\sigma\in\{-1,1\}$. In detail, the signed congruences $(p^g/t)^{d_i}=\operatorname{sgn}(x_i)\bmod W_e$ show by Bézout that $p^g/t=\sigma\bmod W_e$. Since $W_e>2$, the signs then agree as integers. For the bound, $p^g\ne t$, $t^d\le C<p^a$ and $g\le a/d$.

Take any other target $e'$. It agrees with $e$ at some coordinate $j$ and differs by one at some coordinate $i$. Choose these coordinates to include coordinate 1. Set
$$ w=d_j/\gcd(d_i,d_j),\qquad v=d_i/\gcd(d_i,d_j).$$
Then $w\le \min(d,9)$, since the reduced ratio involving $d_1$ has numerator and denominator at most 9. Put $x'_i=e'_i-c_i$. Modulo $W_{e'}$, the unchanged coordinate gives $p^{g d_j}\equiv(\sigma t)^{d_j}$, since $p$ is invertible. At the changed coordinate, raise $p^{a_i}\equiv x'_i$ to power $w$ and substitute $d_i w=d_jv$. Thus
$$W_{e'}\mid (x'_i)^w-x_i^w,\qquad W_{e'}\le w C^{w-1}. \tag{4}$$
The difference is nonzero because adjacent integers have unequal absolute values. The bound follows by factoring the difference of powers, using $|x_i|,|x'_i|\le C$.

If $d\ge3$, (2)–(4) give
$$\log_p(U^3S^2)\le 3n+3(1+a/d)+7(8n+\log_p9)\le a+59n+26<2a.$$
Here the large factor has exponent at most three and the remaining three have total exponent at most seven. If $d=2$, the same bound with $w\le2$ gives at most $3a/2+10n+10<2a$. Both contradict (2). Thus $d=d_i=1$ for all $i$, and the common sign must be negative since at most one $x_i$ is positive. Equation (3) gives $c_i=lp^{a_i-c}+e_i$, with $l>0$. A unit target now gives a scaled triple of powers plus one unit, excluded by the [shifted-row theorem](shifted_power_rows.md). Target zero gives (1).

## Finite heights

Use the [verified initial interval](full_range_1b.md) for $N\le10^9$. Above it each part $y$ of $A$ must satisfy $y^5>27(N-1)^2$ by the [individual-part bound](endpoint_factors.md#lower-bound-on-every-part). For the offsets checked here $(2C)^5<27\cdot10^{18}$. Each of three parts must therefore receive an added unit, and if the smallest added unit is $p^c$ then $(p^c+C)^5>27\cdot10^{18}$ and $p^c>C$. This allocation remains valid when units merge during addition by assigning every constituent to the part receiving the merged unit. Merging two of the three high units is impossible since all three parts need one. Thus in binary $a>b>c$, and at $p=3$ the three cannot all be equal. The digits of $C$ split freely among the three parts. Also $p\mid C(C-1)$, so all possible primes can be enumerated.

The main verifier enumerates the two smaller high units first. Write their actual parts as $y,z$, put $P=p^{v_p(C(C-1))}$, and take $g,L,R$ to be the prime-to-$p$ parts of $\gcd(y,z),\gcd(y-1,z),\gcd(y,z-1)$. Outside $p$, the disjoint factors $U,S_1$ both divide $g$, and $S_2,S_3$ divide $L,R$. The $p$-part occurs at only one target. Consequently
$$27(N-1)^2\le P^3g^3 L^2 R^2.$$
This bounds the largest exponent before forming a full row. For each remaining row, let $G_e=\gcd_i(A_i-e_i)$. Check whether $U^3S^2$ can reach (2) using the largest divisors made of full prime powers of $N,N-1$ allowed by $G_0,\prod_{e\ne0}G_e$. These three nonzero-target caps are pairwise coprime. A row with all three prime-to-$p$ parts equal gives exactly the displayed unshifted form: since $p^c>C$, the leading base-$p$ positions of its entries are $a,b,c$. Subtracting the high terms recovers (1). A row with this shape after subtracting a coordinate unit is covered by the shifted-row theorem.

A second program enumerates the largest exponent first. It uses a different necessary condition: at least one of the four $G_e$ must have prime-to-$p$ part at least $(27(N-1)^2)^{1/9}/P$. It intersects two lists of exponents giving sufficiently large pairwise gcds, generates low rows independently by digit-sum equality and checks the resulting rows with integer arithmetic.

~~~sh
python3 check_three_power_offsets.py --limit 63
python3 check_three_power_offsets.py --one 785
python3 check_three_power_offsets_cross.py 29 83 785
~~~

At offsets 2 through 63 the main verifier checks 17,580,825 binary rows and 5,746,106 odd-base rows after its pair filter. Only unshifted rows survive the full row tests. At 785 the possible bases are $2,5,7,157$, and its respective digit sums are $4,5,5,5$. Therefore (1) is impossible at every one. In the finite branch, 964,179 binary rows and 477,347 odd-base rows pass the pair filter. None survives the full row tests. The independent method tests 10,001, 19, 428 and 70 candidate rows respectively in these four bases and leaves zero possible exceptions.
