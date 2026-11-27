# Powers of one prime at cofactor five

Write $N(k,Q)=1+5\cdot2^kQ$, with $k\ge1$ and $Q$ a positive prime power. Call $Q$ exceptional if the full three-row conjecture fails at **some** such $k$.

**For every prime power $R$, at least 998 of the first 1,000 odd powers $R,R^3,\ldots,R^{1999}$ satisfy the conjecture at every binary exponent $k$.** The possible two exceptional powers may depend on $R$. Within the three-element set $R,R^{121},R^{251}$ there is at most one. The same holds for the pair $R,R^{19}$. These statements allow a different $k$ at every power. In particular, three distinct exceptional powers of the same prime, of exponents $e_1,e_2,e_3$, must satisfy $\max(e_i)>1999\gcd(e_1,e_2,e_3)$ by dividing out the common exponent.

**For any $m$ one can choose $m$ distinct positive odd exponents of size $O(m^{10})$ so that, for every prime power $R$, at most two of the resulting powers of $R$ are exceptional.** The exponents can be determined by finite searches independently of the prime base.

When $1\le k\le35$ there is a stronger [separation theorem for powers of one prime](five_fixed_binary_separation.md). Among the first $10^{18}$ powers of any prime power $R$, at most one pair consisting of a power and a binary exponent can fail. The sequence of odd exponents in that interval extends to one infinite sequence retaining this property.

## The test in the odd prime power alone

The [cofactor-five reduction](predecessor_five_seven.md#full-factor-five) shows that an exceptional $Q$ must be an odd power of a prime congruent to 3 modulo 8. A mandatory row in a counterexample has an odd coefficient and a part $3\cdot2^kQ$. Its other parts sum to $2\cdot2^kQ+1$. A carry in adding these two quantities would make the coefficient even by Kummer's criterion, so the ones of $3Q$ and $2Q$ must be disjoint. For a prime power $R$ with $R\not\equiv3\pmod8$, no odd power can be exceptional. We may restrict to $R\equiv3\pmod8$.

For $n\ge3$ let $S_n$ consist of the odd exponent classes $a\bmod 2^{n-2}$ for which $3q$ and $2q$ have disjoint ones in the lowest $n$ places, with $q=3^a\bmod 2^n$. The powers of 3 generate the group of classes 1 or 3 modulo eight in the units modulo $2^n$. Indeed $v_2(3^{2^i}-1)=i+2$ for $i\ge1$, by successively factoring differences of squares. The odd powers are therefore in bijection with the residues equal to 3 modulo eight. If $R=3^u\bmod 2^n$ then $u$ is odd, and an exceptional $R^a$ requires $ua\in S_n$.

At nine bits $|S_9|=16$. The ratios $xy^{-1}\bmod128$, with $x,y\in S_9$, cover every odd class except 19 and 27. In particular, if positive odd $a,b$ satisfy $b\equiv19a$ or $27a\bmod128$, at least one of $R^a,R^b$ satisfies the conjecture at every $k$. For the three exponents $1,121,251$, the three pairs first become incompatible at respectively 15, 14 and 9 bits, in the order $(1,121),(1,251),(121,251)$. [Lean verification of the nine- and fifteen-bit implications](FivePowerSeparation.lean)

## Checking a thousand entire families

Start with $r=3\bmod8$ and all 1,000 odd exponents from 1 through 1999. At $n$ bits keep the exponents for which $3r^a$ and $2r^a$ have disjoint low bits. A prefix with fewer than three exponents needs no descendants. Otherwise test the lifts $r$ and $r+2^n$ at $n+1$ bits. For odd $r,a$, the two $a$th powers differ by $2^n$ modulo $2^{n+1}$ by the binomial theorem. Thus the new digit of $2r^a$ stays fixed and the new digit of $3r^a$ toggles. The search visits 73,328,239 surviving prefixes from 3 through 120 bits. Its largest level has 4,918,691 prefixes at 30 bits. Only 18 remain at 120 bits and none at 139.

A second traversal indexes residues by $r=3^u$ and lifts $u$ by $2^{n-2}$ instead. It propagates the powers by multiplication and tests the factorial-layer inequality $[3r^a]_{2^{n+1}}\le[5r^a]_{2^{n+1}}$, writing $[x]_m$ for the least nonnegative residue. This is equivalent to no overflow when the residues of $3r^a,2r^a$ are added at that layer. The two independent traversals agree on every level count and on the exact 18 terminal prefixes and their exponent lists at 120 bits. Python rechecks each of those lists against all 1,000 powers from scratch and exhausts the remaining lifts with arbitrary-precision arithmetic. The use of unsigned 128-bit multiplication in the two short traversals preserves all requested low digits.

## Counting exponents and constructing larger sets

Let $\rho$ be the Perron root of the binary transition matrix with successor lists $(1,2),(3),(4),(1,3)$. The states record the previous digit of $q$ and the incoming carry in $q+2q$. The digit of this sum must be zero whenever the previous digit is one. This gives the four transitions and $|S_n|\ll\rho^n$; the matrix is primitive. The [explicit matrix estimate](predecessor_five_seven.md#full-factor-five) gives $\rho<303/200<2^{3/5}$. Set $\sigma=\log_2\rho$.

For each prime $\ell$, at most $O(H^\sigma)$ exponents $e$ in **any interval of $H$ consecutive exponents** can make $\ell^e$ exceptional. The constant is absolute. If $\ell\equiv3\bmod8$, its powers have period $2^{n-2}$ modulo $2^n$. Choose the period between $H$ and $2H$. All exponents in the interval then have different residues, and at most $|S_n|$ pass. The other prime classes give no exceptions. For any fixed positive exponent, the same bound holds on the number of possible prime bases in any interval of $H$ consecutive integers, uniformly in the exponent. For odd exponents the power map permutes the units modulo $2^n$, and for even exponents all cases are already excluded. Thus for every fixed prime, a density-one set of exponents gives the full conjecture at all $k$ simultaneously.

The count also holds for rapidly growing inputs. For $d\ge1$, take $X_0(j)=j$ and $X_i(j)=b_i^{X_{i-1}(j)}$, where $b_i>1$ are fixed integers equal to 3 or 5 modulo eight and $b_d$ is prime. In any $H\ge2$ consecutive positive $j$, at most $O(\rho^{2d}H^\sigma)$ values give an exceptional $X_d(j)$, with an absolute implied constant. In fact $v_2(b^v-b^w)\le v_2(v-w)+2$ for unequal positive exponents, by the same difference-of-squares factorization and the odd-exponent geometric series. Across $d$ iterations the loss is at most $2d$. Distinct indices in the interval therefore give distinct residues modulo $2^{\lceil\log_2 H\rceil+2d}$. There are $O(\rho^{\lceil\log_2 H\rceil+2d})$ passing residues. In particular if $d\le(\log_2 H)/4$, the number of possible indices is $O(H^{9/10})$ even as the height grows.

Finally put $M=2^{n-3}$, the size of the group of odd exponent classes modulo $2^{n-2}$. Across all bases there are at most $M$ translates $u^{-1}S_n$ and each has size $O(M^{3/5})$. Choose $2m$ distinct elements of the group uniformly. The expected sum, over all translates, of the number of triples lying in a translate is
$$
 M\binom{2m}{3}\frac{\binom{|S_n|}{3}}{\binom M3}
 \ll m^3 M^{-1/5}.
$$
Choosing the power of two $M$ between $C m^{10}$ and $2C m^{10}$ for a sufficiently large absolute constant $C$ makes this expectation less than $m$. Some choice therefore has fewer than $m$ bad triples, counting repetitions across translates. Delete one member of each surviving bad triple. At least $m$ members remain and no translate contains three. Their least positive representatives are bounded by $2M$. Exhaustive enumeration of the finite group and its subsets will produce such exponents.

Run the two enumerations and the Lean finite proof from the project folder.

~~~sh
mkdir -p .work
c++ -std=c++17 -O3 check_five_odd_power_blocks.cpp -o .work/check_five_odd_power_blocks
python3 check_five_odd_power_blocks.py .work/check_five_odd_power_blocks
lean FivePowerSeparation.lean
~~~
