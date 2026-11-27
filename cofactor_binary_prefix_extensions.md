# Simultaneous certification above every prescribed binary prefix

Let $D_0,B,J$ be positive integers and put $S=BD_0$. Fix any residue modulo $2^t$, with $t\ge1$ and $2^t$ at least the largest power of two at most $S$. Put
$$
                         n=8+20\lceil\log_2 S\rceil .
$$
At least a fraction $2^{-J}$ of its lifts modulo $2^{t+Jn}$ have the following property **for every positive integer $Q$ in the lift and every $1\le D\le D_0$ simultaneously**. For any two positive trinomial rows $A,C$ at $N=1+DQ$, writing $K_A,K_C$ for their coefficients,
$$
 Q\mid A_i\text{ for some part of }A,\quad \gcd(N,K_A,K_C)\le B
       \quad\Longrightarrow\quad 2^J\mid\gcd(K_A,K_C). \tag{1}
$$
At $B=J=1$ at least half the lifts certify the **full three-row conjecture** in both of the following cases.

- $Q$ is any prime power, for every $1\le D\le D_0$.
- $Q$ is a product of powers of two distinct primes, and $D\le D_0$ is a prime power, twice an odd prime power, 12, 15 or 1. This also holds at 21 and 33 by the parity refinement.

In the first case the two nonendpoint rows of a counterexample must include a row avoiding the base prime of $Q$. Its parts at modulus $Q$ have residues $0,0,1$, so it has a positive part divisible by $Q$. In the second case one nonendpoint row must avoid both bases. Use the full powers of the two primes in $N-1$. The residual cofactor is a divisor of $D$, and every divisor of a listed cofactor is again on the list. The two sets of zero positions intersect, giving a part divisible by $Q$ even if $D,Q$ share prime factors. A counterexample has weight gcd one between its two nonendpoint rows, but (1) makes both coefficients even, and the endpoint coefficient is even. These conclusions use the all-weight endpoint theorem. Shorter prefixes may be extended first to the stated length. [Endpoint reduction](all_weights_endpoint_reduction.md) · [Two-power assignment](two_prime_powers.md) · [Cofactors 21 and 33](four_factor_squares.md)

The same argument works at any other cofactor where the two bases must be avoided by the same row. For example, write $Q=PR$ as powers of distinct primes coprime to $D$. A split assignment forces positive divisors $U,V$ with $UV=N$, $V\ge2P+1$ and $U\ge2R+1$. Thus it is enough that $N$ have no divisor in $[2P+1,N/(2R+1)]$.

## Which pairs of parts need a binary carry

Group the other two parts of $A$ to obtain the cut $aQ,(D-a)Q+1$ with $1\le a<D$. A carry in this cut contributes at least as much to the full row, since splitting a summand cannot decrease a Legendre layer. Its part gcd $u$ divides $N$ and is coprime to $Q$, hence divides $a$ and is at most $D$. Put $m=N/\gcd_i C_i$. The factorial decrement identity shows $N/u\mid K_A$ and $m\mid K_C$. As $m\mid N$, its quotient by $\gcd(m,N/u)$ divides $u$. The gcd in this last denominator divides $\gcd(N,K_A,K_C)$. Therefore $3\le m\le BD\le S$ and $C=(N/m)(f_1,f_2,f_3)$ for positive $f_i$ summing to $m$. These bounds require no assumption on which primes divide $D$. [Normalized denominator argument](digit_gap_transfer.md)

It suffices to force a shared occupied bit for every relevant pair of the forms
$$
            aQ,(D-a)Q+1;\qquad f_1 V,f_2 V,\quad V=(DQ+1)/m .
$$
An overlap forces a positive Legendre contribution in the full row at the corresponding layer. We use at most $D_0^2$ pairs of the first kind and $D_0 S^3$ choices $(D,m,f_1,f_2)$ of the second kind. Throughout, the sum of the two slopes is at most $S$.

## A bound for one new block of input bits

For positive integer slopes $a,b$ with $a+b\le S$, prescribe the first $t'$ bits of an integer scalar $X$. Allow fixed nonnegative offsets $c_1,c_2$ on the two parts. Set $d=\min(v_2(a),v_2(b))$. Among the $2^l$ choices of the next $l$ scalar bits, at most $S\varphi^l$ choices allow the bits of $aX+c_1,bX+c_2$ to be disjoint at positions $t'+d,\ldots,t'+d+l-1$, where $\varphi=(1+\sqrt5)/2$. The units position is 0.

To prove the bound, divide the parts by $2^{t'+d}$ and take floors. The high parts have slopes $a/2^d,b/2^d$ in the new scalar and fixed nonnegative offsets. Take their own sum as the third affine expression. Use the ordinary digits of the first two expressions and the complemented digits of the third. Disjointness of the two $l$-bit words requires the three digits to sum to one in each column. Two of the three slopes are odd and the third has delay $h\ge1$ equal to its 2-adic valuation. Its digit at column $i+h$ is determined when input bit $i$ is read and toggles between the two choices. The two undelayed digits in column $i$ both toggle.

Keep the already determined ones of the delayed output in the next $h$ columns, including the fixed initial digits, with weight $\varphi^{-b'}$ for $b'$ ones. If the current column contains such a one, at most one input can make both current undelayed digits zero. Removing the one multiplies weight by at most $\varphi$. If there is none, the combined relative weight of the possible successors is at most $1+1/\varphi=\varphi$, since their new delayed digits are opposite. Thus the count after $l$ inputs is at most $\varphi^{l+h}$. Now $2^h\le S$ and $\varphi<2$. This is the two-summand case of the [affine digit estimate](binary_affine_counts.md), applied to a finite window regardless of earlier carries.

For the cut in row $A$ use $t'=t$. For $V$, write $m=2^e m'$ with $m'$ odd and $g=v_2(D)$. Since $t\ge e$, the fixed low bits determine whether $2^e\mid DQ+1$. When it divides, $V$ has a fixed prefix of length $t'=t+g-e$. As the next $l$ bits of $Q$ vary, the next $l$ bits of $V$ range bijectively over all words. Indeed, after the fixed prefix, a change by $2^t z$ in $Q$ changes $V$ by $2^{t'}(D/2^g)(m')^{-1}z$ in 2-adic arithmetic, with an odd multiplier. The argument also covers words for which $m'$ does not divide the actual $N$, by using its modular inverse. This only enlarges the test.

The fraction of next words at which even one pair has no overlap in its window is bounded by
$$
          2S^5(\varphi/2)^l . \tag{2}
$$
For $l=n$ this is at most $1/2$, since $\varphi^4<8$. Repeat the construction $J$ times, keeping all universally successful words at each step. Each fixed pair uses disjoint consecutive windows as the starting prefix length increases by $n$. Thus there are at least $J$ different valuation layers in both rows, proving (1).

## Consequences and an explicit example

* For any odd prime $q$, any progression $a+bi$ of positive exponents contains an entire congruence class of indices of bounded step on which (1) holds for all $D\le D_0$. Restrict first to even $i$. Put $v=v_2(q^{2b}-1)\ge3$ and $t=\max(v,\lfloor\log_2 S\rfloor)$. On fixing $i=0\bmod 2^{1+t-v}$ the next $Jn$ digits of $q^{a+bi}$ run through every lift of the resulting prefix, with period $2^{Jn}$ in the new index. This follows by squaring $1+2^v u$ with $u$ odd, increasing the valuation by one at each step. One may therefore take a step at most $2^{1+\lfloor\log_2 S\rfloor+Jn}$ in the original indices, uniformly in $a,b,q$. This gives a residue class whose every member works, complementing the [counts in finite exponent intervals](simultaneous_exponent_intervals.md).
* In $[X,2X)$, outside $O(X^{24/25})$ integer values of $Q$, every positive integer $Q$ in that interval has this property simultaneously for $D\le X^{1/50}$: if either of two rows has a part divisible by $Q$ and at least one coefficient is odd, their coefficients have a common divisor of $N=1+DQ$ greater than $X^{1/50}$. Take $B=D_0=\lfloor X^{1/50}\rfloor$, $t=\max(1,\lfloor\log_2 S\rfloor)$, $l=\lceil\log_2(2X)\rceil-t$, and count the failed words over all initial prefixes. By (2) their number is $O(2^t S^5\varphi^l)=O(X^{3/4} S^{21/4})=O(X^{24/25})$, using $S\le X^{1/25}$. A good word forces both rows even under the small weight-gcd hypothesis.
* For the single class $Q\equiv1058742659284095799\pmod{2^{64}}$, formula (1) holds with $D_0=32$, $B=2$, $J=3$. The verifier checks the overlapping bits directly throughout this class and includes denominators with a factor of two. It also checks the finite-word bound with exact integer counts from an independent enumeration. [Verifier](check_cofactor_binary_prefix_extensions.py)

## A direct lift at any prescribed prime

For a single cofactor $D$ there is a direct construction valid at every prime $p$. It gives actual residues by CRT and forces two contributions from each tested layer of the companion row. Start from any unit $r\bmod p^t$, $t\ge1$. Set $S=BD$, $g=\lfloor\log_p S\rfloor$, and let $L$ be the lcm of all positive integers at most $S$ coprime to $pD$. Choose the unique $v$ in $[2,1+p^tL]$ with
$$
                 v\equiv-r\pmod{p^t},\qquad Dv\equiv1\pmod L .
$$
Let $j$ be the least integer with $p^j>Dv$ and put $H=\max(t,j+J-1+g)$. On the whole class $Q\equiv-v\pmod{p^H}$, the hypotheses on the left of (1) imply
$$
                       p^J\mid K_A,\qquad p^{2J}\mid K_C . \tag{3}
$$
In particular $H-t$ is bounded independently of $r,t$ by $1+\lceil\log_p(D(L+1))\rceil+J+g$.

In fact at each layer $p^i$ with $j\le i<j+J$, the residues of the two parts of the cut are $p^i-av,p^i-(D-a)v+1$, whose sum exceeds $p^i$. For $m$ as above, its prime-to-$p$ part divides $L$ by $m\mid N$ and $(N,D)=1$. Its $p$-part divides $Dv-1$ since $N\equiv1-Dv\pmod{p^H}$. Thus $m\mid Dv-1$. Dividing the weight congruence by $m$ gives $V\equiv-(Dv-1)/m\pmod{p^{H-g}}$. Every $f_i(Dv-1)/m$ is positive, and their sum $Dv-1$ is less than $p^j$. Hence the three residues of $C$ at each layer sum to $3p^i-(Dv-1)$, giving a contribution of two.

For the prime-power and grouped two-power families the lift proves the full conjecture whenever $p$ divides $N(N-1)$, by the same endpoint argument. It extends the infinite-subprogression conclusion to the even prime base $q=2$. For each $D$ with an odd prime divisor $p$, use that $p$. Along $Q=2^{a+bi}$ take first $i$ divisible by $p-1$. If $t=v_p(2^{b(p-1)}-1)$, the powers exhaust all lifts of the initial $p^t$-prefix modulo $p^H$ with period $p^{H-t}$ in the new index. This follows from $v_p((1+p^t u)^{p^n}-1)=t+n$ for $p\nmid u$, by the binomial theorem. This restricts the index to a class of step at most $(p-1)p^{H-t}$. Cofactors with no odd prime divisor already give a predecessor with just one prime factor. Repeating for each $D\le D_0$ yields an infinite subprogression certifying the full conjecture simultaneously at all these cofactors, with step bounded solely in terms of $D_0$. For odd bases the binary word count gives the sharper bound above.
