# Fixed binary exponents at cofactor five, with unbounded odd prime powers

For every odd prime power $Q$ the full three-row conjecture holds at $N=1+5\cdot2^k Q$ for every $1\le k\le30$ and for $k=33,35$. For all $k\le35$ the only residue classes still possible are the following. The table gives necessary conditions for an exception.

| $k$ | $Q\bmod 2^{k+1}$ | $c$ |
|---:|---:|---:|
| 31 | 741526579 | 4192014075 |
| 32 | 5727322795 | 6285315 |
| 34 | 429493347 | 34312377675 |

The exponent and prime base of $Q$ are unrestricted. The [earlier calculation in the other direction](factor_five_exponents.md) fixes the odd power and allows the binary exponent to vary.

## Parametrizing the row by the low bits of the odd power

Write $P=2^k$ and $B=2P$. The [cofactor-five reduction](predecessor_five_seven.md#full-factor-five) proves the result for $k<3$, and otherwise forces one row of a counterexample to have an odd coefficient also coprime to 3 and 5. Its parts, in some order, are $3PQ,zP,(2Q-z)P+1$, where $0<z<2Q$, $zP\equiv1\pmod Q$, $15\mid z$, and the bits of $z$ lie among those of $2Q$. Hence $z$ is even. The last part is a positive multiple $cQ$ with $0<c<B$ (the alternative $zP=1$ is impossible), and is divisible by 15. The reduction gives $(Q,15)=1$, $PQ\equiv7\pmod {15}$ and $3\parallel N$. Therefore
$$
 c=[Q^{-1}]_B,\qquad 15\mid c,\qquad
 A=(3PQ,cQ,(2P-c)Q+1).
$$
Here brackets denote the least nonnegative residue. Further, $3Q$ and $2Q$ have disjoint bits. Enumerate the odd $q_0<B$ for which those two products are disjoint in their first $k+1$ places and $15\mid[q_0^{-1}]_B$. There are 29,258 pairs $(k,q_0)$ for $3\le k\le30$, and 230,427 for $3\le k\le35$.

## Searching all higher digits

For the binary test put $Y=(cQ-1)/B$. The bits of $Y$ must be a subset of those of $Q$, since $Y=(2Q-z)/2$ and $z/2$ split the bits of $Q$. Start with the prefix $q_0$ of $Q$. After reading $k+1+j$ bits, let $v$ be the value of this prefix and retain
$$
 u=\lfloor v/2^j\rfloor,\qquad f=\lfloor(c v-1)/(B2^j)\rfloor .
$$
The next bit $b$ of $Q$ gives bit $(f+bc)\bmod2$ of $Y$ at position $j$. Discard the step if it exceeds $u\bmod 2$. The new pair is $(\lfloor u/2\rfloor+bP,\lfloor(f+bc)/2\rfloor)$. In parallel multiply $Q$ by three as $Q+2Q$ and check disjointness of $3Q,2Q$ at the newly read place. This needs only the previous bit and the incoming carry of that multiplication. We can terminate exactly when padding with zero bits passes both checks. For the subset check this means $f$ is a submask of $u$.

In bases $p=3,5$ the input scalar is again $Q$, with slopes $(3P,c,2P-c)$ and offsets $(0,0,1)$. At 3 initialize with every residue modulo 9 giving $3\parallel 5PQ+1$. At 5 initialize with $PQ\equiv2\bmod5$. Retain a residue only if the low parts add without a carry. From the resulting state of upper parts, a new scalar digit $b$ adds $b$ times the three slopes to that state. Reject if the new units digits sum to at least $p$ and otherwise take the three quotients by $p$. A state admits termination exactly when its three entries add without a carry. Each quotient stays bounded by its slope or initial offset. Thus exhaustion of the reachable states with no termination excludes every finite $Q$ in that base.

The C++ search tries the three bases with increasing resource limits. Only exhausted graphs give exclusions. For five cases reaching the final stage use the complete binary graph and the required residues of $Q$ modulo 45. Starting at each terminating state with suffix $y=0$, traverse a reverse edge labeled $b$ by replacing $y$ with $2y+b\bmod45$. At an initial state test $Q\equiv q_0+B y\bmod45$. This discards two further classes, $c=1402106595$ at 31 and $c=1757311275$ at 32. The other three appear in the table.

The Python verifier independently generates the prefixes and checks every claimed exclusion. It tracks the total and two parts, checking that the two digits fit in the total. At 2 it processes the original row directly, with a separate multiplication test on $3Q,2Q$. It also exhausts the reverse residue search when needed. The comparison checks that every allowed prefix has been accounted for.

~~~sh
mkdir -p .work
c++ -std=c++17 -O3 check_factor_five_binary_exponents.cpp -o .work/check_factor_five_binary_exponents
python3 check_factor_five_binary_exponents.py .work/check_factor_five_binary_exponents 35
~~~
