# All prime powers at cofactor 1651

**The full three-row conjecture holds at $N=1+1651q^b$ for every prime $q$ and every integer $b\ge1$.** This also settles three previously open rays in the [successor-exponent theorem](successor_exponents.md), giving the full conjecture at
$$N=1+p^a(p^j+d),\qquad p\text{ prime},\quad a,j\ge1,\quad 1\le d\le248.$$
Through $d=256$ the only remaining ray is $d=249$, $p=4877$, $j=1$, already checked through $a=10^{24}$.

At cofactor 2222 the full conjecture holds whenever $Q=q^b$ is $1\bmod3$, $2\bmod5$, a power of 3, or a power of 5. In particular **every square prime power $Q$ works at $1+2222Q$**. These results give the unrestricted cases of the [second-interval classification](inverse_pair_digit_gap.md).

The same enumeration proves the full conjecture at **every $N=1+2563q^b$**. None of the 393 primitive candidates at $2563=11\cdot233$ survives the local tests, and there are no candidates at either proper prime divisor.

## Enumerating the rows at these cofactors

Write $N=1+dQ$ with $d=1651$ or 2222 and $Q=q^b$. In a counterexample containing the [necessary endpoint](all_weights_endpoint_reduction.md), let $A$ be a row avoiding $q$. By Kummer its parts are $aQ,b'Q,cQ+1$ with $a+b'+c=d$. All three multipliers are positive, since every part is divisible by the row scale $U\ge3$. After removing $g=\gcd(a,b',c)$ use $W=gQ$ and the primitive denominator $m=d/g$. The [shifted-direction test](shifted_directions.md) applies even when $W$ is composite or shares primes with $m$. The companion parts are $V(f,g',h)$, with $UV=N$ and $f+g'+h=U$.

Two enumerations, one starting with the slopes and one reconstructing the companion by CRT, agree on every primitive candidate at all the divisors $m$ of 1651 and 2222. The only denominators with locally surviving rows and their necessary shapes are below. The three numbers in the first-row column mean the actual row is $(aW,b'W,cW+1)$.

| $m$ | $U$ | First-row slopes | Companion divided by $V$ | Necessary $W$ |
|---:|---:|:---|:---|:---|
|1111|7|$(308,707,96)$|$(1,2,4)$|$7858$ or $81\bmod15554$|
|1651|4|$(104,1524,23)$|$(1,1,2)$|$2297\bmod6604$|
|2222|7|$(77,2121,24)$|$(1,2,4)$|$8101\bmod15554$|

At 1651 only the middle row applies and $W=Q$. It must avoid 2, 13 and 127. There are 297 primitive candidates before these local carry tests. At 2222 the other rows have $W=2Q$ and $W=Q$ respectively, so only the even residue in the first line is relevant. This computation uses no digit bound on $Q$.

## The two cycles in base 13

For the middle shape let $(a_1,a_2,a_3)=(104,1524,23)$. Read the base-13 digits of $Q$ from right to left. Start with the three carries $(0,0,1)$. At input digit $d$ and carries $c_i$, form $y_i=a_i d+c_i$. Allow the transition precisely if $\sum_i(y_i\bmod13)<13$, then replace $c_i$ by $\lfloor y_i/13\rfloor$. A state can end a word if its three remaining integers add freely in base 13. Induction on the number of digits read proves the test from ordinary multiplication. There are 558 reachable states. After the required initial digit $2297\bmod13=9$, exactly 268 states both can be reached and can lead to acceptance.

In this latter graph all directed cycles are self-loops. At $(0,0,0)$ the only allowed continuation leading to acceptance is a zero. At $(69,1016,15)$ the loop digit is 8. Leaving this state permits only zeros afterward. Removing the self-loops leaves an acyclic graph with longest path length 72 from the state after the initial digit. These statements are checked by topological deletion and a reachability test. A canonical input missing the nonzero loop has fewer than 74 digits. Two independent [simultaneous carry searches](inverse_pair_digit_gap.md#the-classification-through-128) already exclude all $2\le Q\le10^{1000}$, far above $13^{74}$.

There are exactly eight prefixes which first reach $(69,1016,15)$. If their low-digit value is $v$ and their length is $t$, appending digits 8 up to length $n\ge t$ gives
$$3Q=2\cdot13^n+C,\qquad C=3v-2\cdot13^t.$$
For each of the eight prefixes the code compares all 63 possible exponent classes with $3Q\equiv3\cdot2297\bmod(3\cdot6604)$. The period divides 63 after the first digit, verified by modular exponentiation. Exactly one prefix survives, with $t=3$, $C=1$, $n\equiv54\bmod63$.

For $3Q=2\cdot13^n+1$ a final four-residue check suffices. At $n=0,1,2\bmod4$ the residues of the three parts sum to at least the modulus at the following layers:

| $n\bmod4$ | Modulus | Part residues |
|---:|---:|:---|
|0|16|$(8,4,8)$|
|1|32|$(8,20,16)$|
|2|16|$(8,4,8)$|

This follows since $2(13^4-1)$ is divisible by 32. In the remaining class $n\equiv3\bmod4$, $5\mid Q$ since $13^3\equiv2\bmod5$. Thus $Q$ must be a power of 5. But then at the first column above the $Q$ shift the three base-five digits are $104,1524,23\bmod5$, namely $4,4,3$. They carry, whereas $A$ was chosen to avoid the prime below $Q$. This exhausts the possibilities.

## The second cofactor

At $N=1+2222Q$, both possible companions are $(V,2V,4V)$ with $U=7$. For $s=3$ or 5, at the first nonzero base-$s$ digit of $V$ the three residues sum to at least $s$, for every nonzero input digit. If $s\mid N$, then $s\mid V$, and the companion must avoid $s$ by the assignment. This is impossible. In particular $Q\equiv1\bmod3$ and $Q\equiv2\bmod5$ are excluded from any counterexample. For $Q$ a power of 3 or 5, the first-row slopes at scale $Q$ are either $(77,2121,24)$ or $(616,1414,192)$. At each of the two primes both slope triples have a carry. There is then no possible $A$ avoiding the prime below $Q$.

## Reproduction

With python-flint installed, run `python3 check_cofactor_1651.py` to compare the enumerations over every divisor and reproduce the table. Run `python3 check_1651_carry_graph.py` to reconstruct the base-13 graph, check its cycles, enumerate the eight prefixes, and verify the modular exclusions. Run `python3 check_inverse_pair_scales.py` for the independent digit-graph and factorial-floor checks of the finite portion.
