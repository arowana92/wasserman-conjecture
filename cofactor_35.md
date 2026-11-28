# Two unrestricted prime powers at cofactor 35

Consider $N=1+35\cdot2^a q^m$, where $q$ is prime and $a,m\ge1$. The full three-row conjecture holds throughout the following families.

- $4\mid m$.
- $m$ is even and either $a$ is odd or $q\equiv\pm1\bmod8$.
- $q\equiv-1\bmod64$, for every $a,m$. This includes $q=127$.

For every fixed odd prime $q$, a set of exponents $m$ of **lower natural density greater than 99.15%** certifies every binary exponent $a$ simultaneously. When $q\not\equiv5\bmod8$ the lower density is one.

There is also a result with two arbitrary prime bases. The conjecture holds at **every $N=1+35t^8$ with $t$ having at most two distinct prime factors**. For a fourth power $N=1+35t^4$ with the same restriction, any failure requires $t$ odd, divisible by 3, and congruent to $\pm3$ modulo both 7 and 8. Every prime-power $t$ works at the fourth power.

## The only possible split

Apply the [endpoint reduction](all_weights_endpoint_reduction.md). Write the two variable powers as $P,Q$ and their product as $W$. If their bases are assigned to different nonendpoint rows, call the row of odd coefficient $A$ and let $P$ denote its assigned power. Both cannot be even in coefficient since the endpoint coefficient is even. Let $U,V$ be the weight factors on the two rows. As in the [general leading-digit argument](leading_predecessor_assignment.md), there are pairs of actual parts $UPi,UPj$ and $VQk,VQl$ with
$$
                  (i+j)(k+l)<35,\qquad 35\mid ijkl .
$$
The odd coefficient requires different binary valuations for $i,j$. In particular at least one is even. The two factors 5 and 7 cannot divide the same multiplier or lie in opposite pairs: these placements give product bounds of at least 72 and 48 respectively. Putting them both in the second pair gives $(i+j)(k+l)\ge3\cdot12>35$. Thus one divides each multiplier of the first pair. The only possibility with an even multiplier and $i+j\le17$ is $\{7,10\}$. The other pair is $\{1,1\}$. Its equal parts make its coefficient even, so the odd row is unique.

For a binary variable power $P=2^a$, the assigned row must be the odd one. Its normalized slopes include $7U,10U$ and partition $35Q$ without carries. They force
$$
       2Q<U<35Q/17,\qquad Q\equiv5\bmod8.                 \tag{1}
$$
The inequalities follow from the positive third parts. For the congruence, $U$ is odd and $7U,10U$ have distinct occupied bits in the lowest three places, all contained in $35Q$. Trying $U=1,3,5,7\bmod8$ leaves only $U=7$, giving residues $1,6$ and $35Q=7\bmod8$. In particular splitting is impossible for square $Q$.

For odd $W=t^4$ the weight satisfies $N\equiv4\bmod16$. Its odd row must avoid 2, so $4\parallel U$ because the equal-parts companion is even in coefficient. After dividing the odd row by four, the two parts $7(U/4)P$ and $10(U/4)P$ need both of the lowest two binary positions, while $N/4\equiv1\bmod4$ supplies only one. This excludes a split. A single variable prime needs no split argument.

## When one row takes both bases

Such a row has a part $xW$ with $3\le U\mid x<35$ and $\gcd(U,35W)=1$. Take a position zero at the full predecessor powers of both variable bases. If a static prime $\ell=5,7$ divides $W$, this choice ensures $\ell\mid x$. With only one variable base use either zero position at its full power. Its binomial cut $xW,(35-x)W+1$ must add freely at every prime assigned to this row. The companion is $V(f_1,f_2,f_3)$ with $V=N/U$ and the positive $f_i$ summing to $U$. These deductions follow from Kummer's criterion. [Bergman, Lemma 5](https://arxiv.org/abs/0806.0607)

Here is the full finite local elimination. At each $\ell\mid35U$, work modulo $\ell^3$, except at 2 use $2^8$. Enumerate the possible scale residues $w$ independently at each prime. With $W=2^aQ$, $Q$ an odd square, these are $2^{a\bmod2}z^2$ at odd primes. With odd $W=t^e$, use $z^e$, $e=4,8$. Use unit $z$ except that a nonunit is allowed at $\ell=5,7$ if $\ell\mid x$. At $\ell\mid U$ require that the valuations of $1+35w$ and $U$ agree up to the working precision, and that the binomial cut is carry-free in that precision. At $\ell=5,7$ test either the cut or the companion, computing $V=(1+35w)/U$ modulo the prime power. At 2 for odd $t$ the row assigned the weight prime 2 must pass the test. For $W=2^aQ$, $q$ odd, instead require $U$ odd and $xQ,(35-x)Q$ carry-free in their lowest four bits, testing all odd square residues for $Q$. This condition follows by dropping the first $a$ positions of the cut. A real assignment must pass every test.

| Input | Surviving $(x,U)$ | Remaining restrictions |
|---|---|---|
| $2^a Q$, $Q$ an odd square | $(9,3),(9,9),(11,11)$ | $a$ even in all 18 assignments |
| Odd $t^4$ | $(16,4)$ | Unique companion proportions $(1,1,2)/4$ |
| Odd $t^8$ | None | All excluded |

At $x=9,11$ the binary cut fails for $Q\equiv1\bmod16$. Every fourth power of an odd integer is one modulo 16, as is the square of an integer congruent to $\pm1\bmod8$. This proves the even-exponent claims. For $q\equiv-1\bmod64$ and odd $m$, the cut residues modulo 64 are $64-x,64-(35-x)$; their sum exceeds 64. Equation (1) also excludes splitting. At $q=2$, the [fixed binary cofactor theorem](binary_cofactor_rays.md) covers factor 35 at all exponents.

For odd $t^4$ in the remaining case $U=4$, the companion cannot take 7 since its normalized sum is four. Hence $16W\equiv1\bmod7$ (if $7\mid W$ then already $7\mid x$ is impossible). This gives $t\equiv\pm3\bmod7$. The first part $16W$ requires the bit of value four after dividing the row by four, forcing $W\equiv17\bmod32$ since $W\equiv1\bmod16$. Thus $t\equiv\pm3\bmod8$. If $3\nmid t$ then $3\mid N$ and $3\mid V$; the companion must avoid 3. After removing its common power of three, all three parts are units, so their lowest digits carry. Hence $3\mid t$ is also necessary. An even $t$ with at most two bases has only one odd prime base, or none, so the binary-variable theorem covers all even fourth and eighth powers. The input $t=1$ is in the [checked initial range](full_range_1b.md).

In the surviving fourth-power case write $t=3^c s$, $3\nmid s$, $c\ge1$. The grouped row also avoids 3. Dividing its binomial cut by $3^{4c}$ and dropping the single units digit in the complement gives the necessary carry-free addition $16s^4+19s^4$ in base 3. Through the first three digits this forces
$$
                  s^4\bmod27\in\{10,13,19,22\}.
$$
For $s=1$ this excludes the last possible prime-power input. More generally if $s=r^e$ for a fixed prime $r\ne3$, the two slopes $16,19$ and their sum are units modulo 3. The [equal-valuation exponent bound](three_factor_exponent_rectangles.md#how-many-exponent-terms-can-pass-a-cut) then leaves only $O(H^{\log_3 2})$ candidate indices $e$ on any $H$-term arithmetic progression, uniformly in $r$, at every $c$. In particular every $9\mid e$ works, since $s^4=1\bmod27$ by Euler's theorem.

## Combining the leading and trailing binary digits

In the grouped case the two multiples $xQ,(35-x)Q$ must add freely in binary, independently of $a$. For each $1\le x<35$ the [exponent digit estimate](three_factor_exponent_rectangles.md#how-many-exponent-terms-can-pass-a-cut) allows only $O(H^{\log_2((1+\sqrt5)/2)})$ exponents $m$ in any interval of length $H$, uniformly over the odd prime $q$. Summing over $x$ gives density zero for these exponents. It remains to bound the split using (1). If $q\not\equiv5\bmod8$, its powers never equal 5 modulo 8, giving the density-one assertion.

Take $s=16$, $L=2^s$, and $q\equiv5\bmod8$. Its powers run once through all $L/4$ classes equal to one modulo four in each period $L/4$: repeated squaring gives $v_2(q^{2^k}-1)=k+2$. For each odd $u\bmod L$ keep only disjoint bit sets of $7u,10u$ reduced modulo $L$. A necessary condition on $Q$ is that $35Q\bmod L$ contain their union for at least one $u$. Exactly 2,743 of the 16,384 scale classes pass.

For $M=35Q$ large, let $F$ be the power of two for which $L\le M/F<2L$. Set $h=\lfloor M/F\rfloor$ and $w=\lfloor70U/F\rfloor$. From (1), $4M/F<70U/F<(70/17)M/F$. Thus
$$
 4L\le w<140L/17,\qquad 17w/70\le M/F<(w+1)/4 .
$$
The integers $\lfloor w/10\rfloor$ and $\lfloor w/7\rfloor$ are the high digits of $7U,10U$ after discarding the places below $F$. They must be disjoint and their union must be contained in $h$. Enumerating these bounds gives 4,243 possible $h$.

On each of the exponent residue classes from the low-digit test, the fractional parts of $\log_2 35+m\log_2 q$ are uniformly distributed: the progression step times $\log_2 q$ is irrational. This follows from the geometric sum for each nonconstant Fourier mode and the equidistribution criterion. The frequency of a given $h$ on that class is $\log_2(1+1/h)$. Therefore the upper density of split exponents is bounded by the product of the low-class proportion and the sum of these frequencies. With $J=2^{30}$, the sum of $\lceil J/h\rceil$ over the 4,243 prefixes equals 37,294,869. Using $\log(1+1/h)<1/h$ and $\log2>69/100$ gives
$$
 \frac{2743}{16384}\frac{100\cdot37294869}{69\cdot2^{30}}<0.0085 .
$$
The logarithm lower bound follows from the first two terms of $2\operatorname{artanh}(1/3)$. Together with the density-zero grouped cases, this proves the stated lower density, with the quantifier over all $a$ inside the conclusion for each $m$.

Run **python3 check_cofactor_35.py**. It repeats the local elimination using factorial quotients and individual digit sums, reverses the leading-prefix enumeration through nine bits, and checks the trailing masks against direct three-summand additions through the same depth.
