# Prime bases that work at every cubic exponent

There is an explicit union of prime residue classes of relative natural density
$$
 \frac{188285}{188384}=0.999474477662646\ldots
$$
on which **every $1+986q^n$ satisfies the full three-row conjecture, for every exponent $n\ge1$ simultaneously**. On these same prime bases every $1+dq^{3e}$ works for every $d\le1000$ and $e\ge1$, by the [earlier cube reduction](denominator_residue_hulls.md). In particular an exceptional base must satisfy $q\equiv-1\pmod{204}$.

The surviving fraction $99/188384$ is exact for the simultaneous local tests below. For each prime class in that fraction the tests can be met to **any specified finite number of digits** at all four carry bases, with a common exponent divisible by 3. These classes still require the actual coefficient divisibility tests at the full integer.

## The local conditions

The two [rows at 986](cofactor_986_powers.md) have common companion $(V,2V,4V)$, $7V=986W+1$, $W=q^n$ odd. Denote the rows by their first slopes 119 and 203. The binary carry graph for each distinguished row has exactly two productive states, with a single outgoing edge labeled zero in each. Its only accepted nonnegative input is zero. Thus the companion must avoid 2. Equivalently the binary ones of $V$ are separated by at least three positions.

The companion coefficient is divisible by 3 and 5 at every positive $V$. It contains the factors $\binom{3V}{V}$ and $\binom{5V}{V}$ individually, respectively multiples of 3 and 5 by $\binom{hV}{V}=h\binom{hV-1}{V-1}$. All primes dividing $V$ are assigned to the companion in this row classification. Hence $3,5\nmid986q^n+1$, requiring $q\equiv2\bmod3$ and $q\not\equiv-1\bmod5$, since $n$ is odd. The primes dividing the finite modulus in this note are already [settled at all exponents](prime_base_density.md).

Use the following tests to identify the remaining classes. At 29 write $q\equiv10^g$ with $g$ a unit modulo 28. For the first row let $t=1$, and for the second let $t=19$. Put $n_0\equiv t/g\bmod28$. The positive exponent must have this residue. For the finite tests take its unique representative divisible by 3 modulo 84.

| Coordinate of $q$ | Necessary test |
|---|---|
| $7^2$ | $q^{n_0}$ modulo 49 lies in $\{22,29\}$ for the first row or $\{22\}$ for the second. |
| $17^2$ | $q\equiv-1\bmod17$ and $q\not\equiv-1\bmod289$. |
| $29^3$ | In addition to primitivity modulo 29, $v_{29}(q^{28}-1)\le1$ for the first row or $\le2$ for the second. |
| $2^7$ | For $n_0\equiv1\bmod4$ take $q\bmod128$ in $\{7,11,15,27,39,43,47,59,71,75,79,91,95,103,107,111,123\}$. For $n_0\equiv3\bmod4$ use $\{3,15,19,23,31,35,47,51,55,67,79,83,87,99,111,115,119\}$. |

For the 7-test, $q$ has order 1 or 3 modulo 7. If the order is 3 every possible exponent is divisible by 3. If the order is 1 its residue modulo 7 suffices to determine $q^n\bmod49$, so the chosen representative works in both cases. At 17 the prescribed $W\equiv84\bmod289$ is attainable exactly for the indicated $q$: write $q=-1+17v$ and impose $vn\equiv5\bmod17$. The low-prefix and primitivity assertions are proved in the [exponent comparison](cofactor_986_powers.md#a-single-row-and-exponent-class-for-each-prime).

At 29 let $h=v_{29}(q^{28}-1)$. For fixed $n_0\bmod28$, $q^{n_0}\bmod29^h$ is the unique lift of $w=10$ or 21 with $w^{28}=1$ at that modulus. Uniqueness follows by lifting one digit, since $28w^{27}$ is a unit. The lift at precision $i$ is $w^{29^{i-1}}\bmod29^i$. For the first row the lift 416 modulo $29^2$ forces a carry. For the second, 137 modulo $29^2$ extends to a carry-free input, whereas 9388 modulo $29^3$ forces a carry. [Factorial-layer witnesses](common_power_exponents.md)

At 2 put $h=v_2(q^4-1)$. The powers with the chosen exponent class occupy $W\equiv q^{n_0}\bmod2^h$ at the first $h$ bits. Accordingly $V$ has a fixed prefix of $h+1$ bits since $v_2(986)=1$ and 7 is odd. Test its binary spacing. A passing prefix extends by zero bits. For $h\ge10$ the odd prime $q$ is $\pm1\bmod256$, since $h=2+\max(v_2(q-1),v_2(q+1))$. In either case the necessary nine-bit prefix of $V$ fails the spacing test. In the other cases $h\le9$. Exhaustion modulo 1024 gives exactly the indicated classes modulo 128. In particular all of them are 3 modulo 4, which combines with the conditions at 3 and 17 to give the simple modulus 204.

## Simultaneous attainability and density

To see that the surviving classes really pass the stated local tests to arbitrary depth, vary the exponent as $n=n_0+84z$. For an odd $p$, if $v_p(s-1)=h\ge1$, the powers of $s$ modulo $p^k$ fill the subgroup of units congruent to 1 modulo $p^h$. The same holds at 2 for $h\ge2$. The binomial expansion gives $v_p(s^j-1)=h+v_p(j)$, so the order equals the size of that subgroup. Use $s=q^{84}$.

At 7 this subgroup consists of all lifts from modulus 49, since $v_7(q^3-1)=1$ and $84=3\cdot28$. At 17 all lifts from the first digit can be reached. At 29 all lifts from modulus $29^h$ can be reached, and at 2 all lifts from modulus $2^{v_2(q^4-1)}$ can be reached. Use one accepted continuation of each distinguished-row prefix at the odd primes. For the binary companion choose the spaced integer that extends its prefix by zeros. The map $W\mapsto(986W+1)/7$ sends each binary residue class modulo $2^h$ onto the corresponding class modulo $2^{h+1}$ at finer 2-adic precision. Each chosen finite continuation therefore restricts $z$ to a class modulo a power of its own prime. The Chinese remainder theorem realizes all the continuations at once. At 3 and 5 the exclusion of prime divisors of $V$ holds automatically on the whole exponent class. This proves the local attainability assertion.

For each of the twelve primitive residues at 29, the fraction of its lifts with $h\le b$ is $1-29^{-b}$, by unique lifting of the torsion root. At 7 there are three possible residues modulo 7 with $q^3=1$. For each one and each allowed target at 49 there is one lift, since $7\nmid n_0$. The two rows' fractions at 7 are therefore $6/42$ and $3/42$. The fractions at 17, 2 and the combined modulus 15 are $1/17,17/64,3/8$. Multiply across the coprime moduli by CRT and use the prime number theorem in fixed arithmetic progressions. The resulting masses are $9/25984$ and $135/753536$. They are disjoint by the modulo-7 and modulo-29 ratio test for a common base, so their sum is $99/188384$.

The independent calculation enumerates the binary orbit, determines the order-29 lifts using a discrete-log table at $29^3$, and enumerates the possibilities at 7, 17 and 15. Both programs use only standard-library Python.

~~~~sh
python3 check_986_prime_classes.py
python3 check_986_prime_classes_cross.py
~~~~
