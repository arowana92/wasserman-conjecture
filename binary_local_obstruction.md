# Simultaneous carries can escape every prescribed digit depth on a fixed binary ray

Let $q=2^a$, $N=12031q+1$ and take three rows
$$E=(1,1,N-2),\qquad A=(681q,9222q,2128q+1),\qquad C=(N/3,N/3,N/3).$$

**Theorem.** Given any $L\ge1$, there is an infinite arithmetic progression of positive exponents $a$ with $3\parallel N$ on which the following hold **simultaneously for every prime $r<6.7\times10^{15}$**.

- At $r=2$ the coefficient of $A$ is odd.
- At odd $r\nmid N(N-1)$ the coefficient of $E$ is a unit.
- At odd $r\mid N(N-1)$ we have $v_r(N(N-1))=1$. At least one of $A,C$ adds without a carry through the layer $r^{L+1}$, including all lower layers.

More generally the same assertion holds for any finite set of primes in which every odd $r$ satisfies $2^{r-1}\not\equiv1\pmod{r^2}$ or belongs to $\{1093,3511\}$. The numerical range follows from the published exhaustive search of Dorais and Klyve. The verifier in this repository independently checks that condition for all primes through one million. [Dorais and Klyve, *A Wieferich Prime Search up to* $6.7\times10^{15}$](https://cs.uwaterloo.ca/journals/JIS/VOL14/Klyve/klyve3.html)

Thus no test that requires a simultaneous Kummer carry within a fixed number of layers beyond the full endpoint powers at those primes can exclude this entire ray. An actual common divisor can arise from higher digits or another prime. All exponents through $10^{24}$ on this ray satisfy Wasserman's conjecture. [Bounded-cofactor theorem](binary_cofactor_rays.md)

## Three finite witnesses

Set $d=12031=53\cdot227$. The slopes $681,9222,2128$ sum to $d$ and have pairwise disjoint binary expansions. The extra 1 occupies a fresh position for positive $a$, proving the binary assertion.

Start with $a=b\pmod D$, where
$$b=106675,\qquad D=22211280=2^4 3^3\cdot5\cdot7\cdot13\cdot113.$$
This progression lies inside $a\equiv907\pmod{5876}$, the necessary progression in the cofactor theorem. For $p=3,53,227$ the following positive scalar $t_p$ makes $(681t_p,9222t_p,2128t_p+1)$ carry-free in base $p$. The three zero valuations follow by summing the usual integer factorial floors.

| $p$ | $t_p$ | $m=p^{1+v_p(D)}$ | $2^b\bmod m=t_p\bmod m$ | $2^D\bmod pm$ |
|---:|---:|---:|---:|---:|
| 3 | 17140160990152208 | 81 | 20 | 163 |
| 53 | 33 | 53 | 33 | 2333 |
| 227 | 8 | 227 | 8 | 13394 |

In each case $2^D-1$ has exact $p$-adic valuation $v_p(m)$. Hence $2^{b+Dn}$ can agree with $t_p$ to any prescribed precision. Indeed if $2^D=1+c p^j$ with $p\nmid c$, multiplication by $(2^D)^k$ permits every next digit as $k$ runs from 0 to $p-1$. After fixing that digit the exponent step is multiplied by $p$, increasing the valuation by exactly one by binomial expansion. The primes here are odd. The three refinements are compatible because each only multiplies the step by powers of its own prime. Match through at least $p^{L+1}$ at each prime. Since $dt_3+1\equiv6\pmod9$, every resulting $N$ has $3\parallel N$. The two odd factors of $N-1$ are always exactly 53 and 227, each to the first power.

For $r=5,7,13,113,1093,3511$, respectively, $r-1\mid D$ and the residues of $d2^b+1$ are $4,4,2,98,851,3458$. Hence none divides $N$ anywhere in the original progression. These six primes account for every odd prime factor of $D$ outside the three witnesses and both exceptional Wieferich primes in the numerical range.

## Adjoining the other primes in increasing order

Suppose the current progression is $a=b'\pmod M$. For the next odd prime $r$ still to treat, $r\nmid M$. Initially all prime factors of the step belong to $D$ or to the three witness primes. At each subsequent step the only possible new prime factors are $r$ and divisors of $r-1$, which justifies the increasing order.

Replace $M$ by $\operatorname{lcm}(M,r-1)$ while retaining $b'$. This fixes $2^a\bmod r$. If $r\nmid d2^{b'}+1$, the endpoint already avoids $r$ and all future refinements preserve this. Otherwise choose the target
$$2^a\equiv (3r-1)d^{-1}\pmod {r^{L+1}}.$$
It is the existing residue modulo $r$. Since $r$ is ordinary at base 2 and $r\nmid M$, the binomial expansion gives $v_r(2^M-1)=1$ from $(r-1)\mid M$. The same one-digit lifting argument realizes the target on a subprogression. Now
$$N\equiv 3r\pmod {r^{L+1}},\qquad N/3\equiv r\pmod {r^{L+1}}.$$
The prime $r$ is at least 5. Therefore $r\parallel N$ and the three equal parts of $C$ add without a carry through that layer. Processing a finite list leaves an infinite progression. This proves the theorem even when the depth is measured after stripping the forced factor $r$ out of each part of $C$.

Run `python3 check_binary_local_obstruction.py` to verify the seeds, the Wieferich scan through one million, and the example at cofactor 1903. The program constructs full progressions for all primes through 200 at depth 32 and for all primes through 4096 at depth 12. It checks actual row residues at two positive exponents in each resulting progression. At the latter cutoff the extra tested weight primes are 19 and 229, each dividing $N$ exactly once.
