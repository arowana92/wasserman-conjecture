# Every power index admits a cofactor with arbitrarily delayed common carries

Let $\mathcal H_D$ be the [finite divisibility list's union](cofactor_index_principle.md) for the bounded-digit endpoint test on $N=1+dq^e$, $1\le d\le D$, with $q>d+1$ prime. **No positive exponent belongs to all $\mathcal H_D$.** For each $e\ge1$ one may take $d=247s\ll e^5$ with $s>247$ prime and a single pair of candidate row shapes at that cofactor. For every prime cutoff $B$ and digit depth $k$, there are arbitrarily large prime bases $q$ such that each prime through $B$ clears its first $k$ digits in at least one of the three rows. The conjecture holds for every triple at the constructed weights. Consequently every exponent in $\mathcal H_D$ is $\gg D^{1/5}$, with an absolute implied constant.

There is also a result **simultaneous in every exponent** when the precision is fixed first. For any finite set of positive multipliers $c$ and any $B,k$, there are two independent reduced arithmetic progressions of primes $s,q$ such that the full conjecture holds at **all** $1+c s q^b$ for every $b\ge1$. At $c=247$ the selected three rows hide all common carries from the indicated low digits for every $b$, although both 2 and 3 divide their three coefficients. A prescribed prime larger than $\max(B,247)$ can certify this weight and divide all three coefficients within its first two digits. The precision, cutoff and finite set of multipliers can all be arbitrarily large.

## Choosing one prime cofactor for an unbounded digit depth

Write $W=sP$, with $P=q^e$. Use the rows
$$
 E=(1,1,N-2),\qquad A=(114W,117W,16W+1),\qquad C=(V,V,V),\qquad 3V=N=1+247W .
$$
The assigned factor of $N$ on $A$ is $U=3$. At $r=2,3,13,19$ respectively the following rational targets for $W$ are available:
$$
 t_2=-1/247,\qquad t_3=-16/121,\qquad t_{13}=43,\qquad t_{19}=127 .
$$
At 2 the companion vanishes locally. At 3 the first row has carry-free eventually periodic digits and $v_3(N)=1$. At 13 and 19 the first row adds without a carry at the respective finite integer targets. [Digit expansions and checks](cofactor_first_local_failure.md)

At these four primes put $h_r=1+v_r(e)+{\bf1}_{r=2}$ and require $s\equiv t_r\bmod r^{h_r}$. Then $t_r/s$ is an $e$th power in the local units. Indeed all units congruent to one at this precision are $e$th powers by the cyclic unit groups at odd primes and the sign and principal-unit factors at 2. The same holds for every subsequent extension of that prefix. For the exact first binary carry prescribe
$$
                    1+247sP\equiv 2^j\pmod{2^{j+1}},\qquad j\ge\max(k,h_2).
$$
The scalar target is unchanged modulo $2^{h_2}$. The endpoint and $A$ carry in the units column and $C$ first carries at digit $j+1$, numbering the units digit as 1.

At each remaining prime divisor of $e$ require $s\equiv1\bmod r$, except use $-1$ at 31 and 251. Impose $s\equiv-1\bmod251$ in every case. These choices keep $r$ out of $1+247s$, since $248=8\cdot31$ and neither exceptional odd prime divides $246$. The CRT class consists of units and its modulus is at most $251\cdot4\cdot3\cdot13\cdot19\,e$. Its least prime exceeds 247 by the condition at 251 and is $O(e^5)$ by Xylouris's Theorem 2.1. [Xylouris, 2011 dissertation, Theorem 2.1](https://d-nb.info/1017437416/34)

Here are the remaining local choices for $P$. At $r=s$ take $P=1$. After removing the $s$ shift the slopes $114,117,16$ sum to $247<s$, so $A$ avoids $s$. At every other prime outside the four specified above take $P=1$ whenever $r\nmid1+247s$. If $r\mid1+247s$, then $r\nmid e$. Use $P=-1/(247s)$ in the $r$-adic units. Its residue is one modulo $r$, so it has local $e$th roots and gives $V=0$ to every precision. At the actual varying prime $q>247s+1$ the three slopes $114s,117s,16s$ sum below $q$, so $A$ avoids $q$.

These data pass the endpoint tests at every depth. Choose local root residues at any desired precisions and use CRT and Dirichlet's theorem to obtain arbitrarily large actual prime bases $q$. To certify the full conjecture on that progression as well, refine the binary scalar prefix for $P$ using the [binary extension theorem](cofactor_binary_prefix_extensions.md) for the cofactor bound $247s$, or for any larger fixed bound. The original binary prefix is long enough that every refinement is locally an $e$th power. Its congruence for the first carry is preserved. This proves the assertion for each fixed $e$. The [simultaneous cofactor result](simultaneous_cofactor_obstructions.md) gives one polynomial-size cofactor for every $e$ in a prescribed initial range.

## Two mandatory common primes with independently delayed carries

For every positive integer $W$ the first-row coefficient $K_A$ is divisible by 3. At an integral companion the weight is divisible by 3, so the endpoint coefficient is divisible by 3. The three equal positive parts in $C$ always carry in ternary at their first nonzero digit. Thus 3 is a genuine common divisor even when the first row clears arbitrarily many ternary digits.

Here is a finite proof of the first assertion. When the incoming multiplication carries are $(x,y,z)$, an input digit $i$ forms $(x+114i,y+117i,z+16i)$. If its three units digits sum to at least 3, the row carries. Otherwise the successor consists of the three integer quotients by 3. Start at $(0,0,1)$; input zeros lead to $(0,0,0)$. A first nonzero digit from either state goes into a list of 25 triples. Each permitted transition from the list stays in it. With scalar input exhausted, every listed triple itself has a ternary carry within five digits. These two finite assertions and the induction on the entire positive integer $W$ are checked in Lean. Its conclusion exhibits a finite Legendre layer at which the three residues sum to at least the modulus. [Lean proof of the universal ternary carry](Ternary247.lean)

One can place the first ternary carry in $A$ exactly. Fix $j_3\ge\max(k,h_3,4)$ with $5\nmid j_3$, and use $W\equiv t_3\bmod3^{j_3}$. At ternary position $j_3$ (units position zero) the first two output digits have positive sum by the five-column period in the [reference calculation](cofactor_first_local_failure.md). Both slopes are divisible by 3, so the next scalar digit does not change those two digits. The slope 16 is a ternary unit, so that digit can set the third output to 2 and force a carry. The first $j_3$ digits were carry-free. The weight still has valuation exactly 1 at 3 and both $C$ and the endpoint carry at their units. Consequently the gcd has ternary valuation exactly 1, first witnessed in all rows by digit $j_3+1$.

Independently, use $1+247W\equiv-3\cdot2^{j_2}\bmod 2^{2j_2}$ with $j_2\ge\max(k,h_2)$. The first common binary carry is still witnessed at digit $j_2+1$. Now $V/2^{j_2}\equiv-1\bmod2^{j_2}$, so each of the next $j_2$ layers contributes a carry to $K_C$. Each of the first $j_2$ layers contributes to $K_A$ because the total is zero there and two of its parts are odd. Since $N(N-1)$ has binary valuation exactly $j_2$, the gcd has binary valuation exactly $j_2$ too. All the new scalar residues preserve the local power images after division by $s$.

## Changing the order of the quantifiers

Fix $B,k\ge1$ and a finite set of positive multipliers. Include 247 among them. Prescribe a residue $w_0$ for $W$ as follows. At 2 use $1+247W\equiv-3\cdot2^k\bmod 2^{2k}$. At 3 use the preceding exact construction with any $j_3\ge\max(k,4)$ not divisible by five; no power restriction is needed on $W$ here. At 13 and 19 use their targets through $k$ digits. At all other primes through $B$ take $W\equiv1$, except take $-1$ at 31. At those other primes $N(N-1)$ is a unit. These conditions alone establish the claimed low-digit behavior of the selected rows whenever $W$ is a sufficiently large positive integer in the resulting class.

For each multiplier $c$ choose a distinct prime $p_c>\max(B,247,c)$, and add the congruence
$$
                   W\equiv(p_c-1)c^{-1}\pmod {p_c^2}.
$$
Thus $p_c\parallel 1+cW$ and $(1+cW)/p_c\equiv1\bmod p_c$. The [weight-prime criterion](predecessor_weight_gap.md) proves the full conjecture at $1+cW$ for **every** such $W$ with at most two distinct prime factors. At $c=247$ the prime $p_c$ divides the endpoint coefficient. The first row carries at the units position since it has nonzero residues summing to zero modulo $p_c$. For $C$, the three identical entries divided by $p_c$ have units digit $3^{-1}\bmod p_c$, whose triple is greater than $p_c$. This prime therefore witnesses a common carry within two digits while all primes through $B$ still miss one. Its valuation in the gcd is exactly one since it divides the weight exactly once.

Let $M$ be the product of the required prime-power moduli. The combined class $w_0\bmod M$ is reduced. Take **any** primes $s\equiv w_0$ and $q\equiv1\bmod M$. Then $W=s q^b$ satisfies the same congruences at every positive exponent $b$. One can also take $W=s^a q^b$ for every positive $a\equiv1\bmod\varphi(M)$. Both prime classes have positive relative density, independently of one another.

These progressions also permit all three parameters to grow. For an integer $n\to\infty$ take $B=k=n$ and use all $c\le n/(2\log n)$ together with 247. Choose the distinct $p_c$ in $(n,2n]$, which contains enough primes by the prime number theorem, and use $j_3=n$ or $n+1$. The high-precision powers in $M$ contribute at most $n\log(2^2\cdot3\cdot13\cdot19)+O(1)$ to $\log M$. The remaining primes through $n$ contribute at most $(1+o(1))n$, and the auxiliary squares at most $(1+o(1))n$. Thus $\log M<11n$ for large $n$. Both least primes in the two reduced classes exceed $n$ and are $O(M^5)$. With $L=\log\max(s,q)$ this gives $L<56n$. For infinitely many such pairs the full conjecture holds simultaneously for all $b\ge1$ and $c\le L/(120\log L)$, while at multiplier 247 all primes through $\lfloor L/56\rfloor$ miss a common carry through that many digits. The stated range of $c$ follows from monotonicity of $n/\log n$ for large $n$. [Least primes in progressions, Theorem 2.1](https://d-nb.info/1017437416/34)

For a fixed cofactor $d$ the exact first binary carry of three equal parts cannot be delayed without limit simultaneously at exponents $b$ and $2b$ of the same odd prime. In fact
$$
 \gcd(1+dq^b,1+dq^{2b})=\gcd(1+d,q^b-1).
$$
Subtract the second entry from $q^b$ times the first, then reduce the first modulo $q^b-1$. Hence the smaller of the two binary weight valuations is at most $v_2(d+1)$ for odd $d$. This explains why the cofactor can be fixed before the precision in the first theorem, and changes with the precision in the simultaneous construction.

## A certified pair that works at every exponent

The example uses $B=47$, $k=8$, and $c=1,2,3,4,5,6,247$. At 247 it prescribes $p_{247}=251$. For $c=1$ through 6 it uses $257,263,269,271,277,281$ in order. The prime bases are
~~~text
s = 1306232213153877917593933367823536683066159885737114363997573716868895663450681
q = 322880541893073397389196821002188355954113056549133075089295162936641374385930241
~~~
For every $b\ge1$ the full conjecture holds at all seven weights. At $1+247s q^b$, the first common binary carry in the sense of having seen a carry in all three rows occurs at digit nine. The first common ternary carry is also witnessed at digit nine. The gcd has valuations 8 at 2 and 1 at 3. The prime 251 works by digit two and divides the gcd exactly once. Lean checks the uniform weight congruences at every exponent. The independent verifier recursively certifies both primes and directly checks the Legendre floors for six exponents, including a 66,327-bit weight.

~~~sh
python3 check_all_exponent_cofactors.py
python3 check_all_exponent_cofactors_cross.py
lean UniformCofactorPowers.lean
lean Ternary247.lean
~~~

The Python programs use python-flint. The Lean files use Lean 4.36.0.
