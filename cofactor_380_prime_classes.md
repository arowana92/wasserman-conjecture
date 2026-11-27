# Prime bases that work at every exponent for cofactor 380

An explicit set of prime bases of relative natural density greater than **99.85366%** satisfies the full conjecture at $1+380q^e$ for every $e\ge1$. For each prime $q$, **at most one residue class of exponents modulo 58,500** can contain any failures. At the three previously [unresolved small bases](prime_base_density.md), the necessary exponent classes are 51,541 for 1,733, 52,433 for 1,973 and 42,971 for 5,717.

## The exponent and four carry bases

By the [classification of primitive denominators](shifted_directions.md), the only possibility has odd scale $W=q^e\equiv823\pmod{2470}$ and rows
$$
 A=(130W,247W,3W+1),\qquad C=(V,4V,8V),\qquad 13V=380W+1.
$$
The endpoint is $(1,1,380W-1)$. The first row must avoid 5, 13 and 19, and the second must avoid every prime of $V$. At the binary layer four the residues of $A$ sum to five for either odd residue of $W$. Hence $C$ must avoid 2. The shifted bit strings of $V,4V,8V$ are disjoint exactly when any two ones of $V$ are at least four positions apart. The coefficient of $C$ is also divisible by 3, since it contains $\binom{12V}{4V}=3\binom{12V-1}{4V-1}$. Therefore $3\nmid V$ and $W\not\equiv1\pmod3$.

Starting at $W\equiv3\pmod5$ in the addition for $A$, the only residue modulo $5^4$ that extends to a finite carry-free input is 458. Starting at $W\equiv4\pmod{13}$, the only such residue modulo $13^2$ is 56. Starting at $W\equiv6\pmod{19}$, the residue $234=6^{19}\bmod361$ cannot extend. Both checkers use the [carry-graph criterion](carry_intersections.md). One retains every state from which some word terminates without a carry, then reads the low digits. The other tests each candidate prefix by factorial floors, then searches its continuation states for zero.

Put $M=625\cdot169\cdot19=2{,}006{,}875$. The three constraints fix $W\equiv R=1{,}471{,}708\pmod M$. Its three coordinate orders are $500,78,9$. Their least common multiple $L=58{,}500$ equals the exponent of the whole unit group (take the lcm of $500,156,18$). Thus $q^e=R$ in the group requires $q$ of order $L$ and $e$ coprime to $L$. Write $q=R^j$ modulo $M$ with $\gcd(j,L)=1$. Necessarily $e\equiv j^{-1}\pmod L$. Conversely every such base class attains $R$. Exactly 14,400 of the $500\cdot156\cdot18$ unit classes pass. This proves the single-exponent-class assertion and implies $q\equiv2\pmod3$ unless $q=3$. At $q=3$ the residue modulo 13 has order three and cannot produce 4 of order six. The bases dividing $2M$ are excluded by the odd unit scale.

If $q^9\equiv1\pmod{361}$, then $W$ must be the unique ninth root of unity modulo 361 above 6, namely 234. Uniqueness follows by expanding $(6+19t)^9$ modulo 361 and using $19\nmid9\cdot6^8$. Among the nineteen lifts of each surviving base modulo 19, exactly one is lost this way.

Binary spacing forces $V\equiv1\pmod{16}$ and $W\equiv1\pmod4$ since $V\equiv1\pmod4$ to begin with. If $W\equiv1\pmod{16}$ then $V\equiv49\pmod{64}$, which already violates the spacing. Since $e$ is odd, $h=v_2(q-1)=v_2(W-1)$ is 2 or 3. When $e$ is specified modulo four, $W$ is fixed modulo $2^{h+2}$, because $2^{h+2}\mid q^4-1$. The relation for $V$ fixes another two bits. Testing these bits gives

| $e\bmod4$ | Allowed $q\bmod32$ |
|---:|---|
| 1 | 5, 9, 13, 21, 29 |
| 3 | 5, 13, 21, 25, 29 |

For each of the 14,400 classes the three additional proportions are $18/19$ at 19, $1/2$ at 3, and $5/16$ at 2. Multiplication by the Chinese remainder theorem leaves a fraction $3/1976$ of the unit classes. The prime number theorem in fixed arithmetic progressions gives the same relative density among primes.

## Common primes from the companion row

There are further explicit common divisors at primes outside the forced carry bases. For example, if $W\equiv1\pmod{127^2}$, then $127\mid 380W+1$ and $V/127\equiv3/13\equiv10\pmod{127}$. The digits $10,40,80$ sum to 130, forcing a carry in $C$. The units digits of $A$ are $3,120,4$, summing to 127. Thus 127 divides all three coefficients in the only candidate configuration.

More generally, choose a prime $r\nmid380L$ with $(-380)^L\equiv1\pmod r$. Put $z=-380^{-1}\bmod r$. At a layer $m=r^k$ let $T=z^{r^{k-1}}\bmod m$. This reduces to $z$ modulo $r$ and satisfies $T^L\equiv1\pmod m$ by expanding an $r$th power successively. If $W\equiv T\pmod m$, then $r\mid N$ and $A$ carries at $r$: its two unshifted slopes have gcd 13, and $r\ne13$ cannot divide both. Compute $v=(380T+1)/13$ modulo $m$, using the modular inverse. If
$$ (v\bmod m)+(4v\bmod m)+(8v\bmod m)\ \ge m,$$
then $r$ is a common prime for the candidate configuration. For any surviving base class modulo $M$ its exponent is fixed at a unit $e_0$ modulo $L$. The single unit class $q\equiv T^{e_0^{-1}\bmod L}\pmod m$ forces $W\equiv T$ at **every** such exponent. It removes a fraction $1/\varphi(m)$ regardless of $e_0$.

There are 53 distinct primes passing this test through 100,000 at some $2\le k\le5$. Take the first passing $k_r$ for each. Their moduli are relatively prime to one another and to the earlier coordinates. The surviving fraction for the listed tests is exactly
$$ \frac{3}{1976}\prod_r\left(1-\frac1{(r-1)r^{k_r-1}}\right)=0.001463365295974897\ldots .$$
This leaves the stated certified density. The finite set of primes dividing the combined modulus does not affect the density. Passing these necessary tests does not exhibit a counterexample.

Run both checks with standard-library Python 3.9 or later. The second computes the lifted roots by their derivatives and checks actual factorial floors at positive integral lifts.

~~~~sh
python3 check_380_prime_classes.py
python3 check_380_prime_classes_cross.py
~~~~
