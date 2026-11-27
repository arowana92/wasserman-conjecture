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

## The leading digits at one selected prime

Fix any one of $p=3,53,227$ and any $J\ge1$, in addition to the preceding finite prime list and lower depth. Infinitely many exponents from the constructed progression also have **no carry in $A$ at any of the $J$ highest layers in base $p$**. The choice of the single leading base is made before choosing the exponents. Thus inspecting both ends in that base, together with the prescribed low digits in all the other bases, still leaves infinitely many exponents.

Let $s=(681,9222,2128)$ and $F(x)=\lfloor dx\rfloor-\sum_i\lfloor s_ix\rfloor$. It is periodic modulo one. Each rational $x$ below satisfies $1<dx<p$ and $F(p^j x)=0$ for **all** $j\ge0$. The calculation is finite because the forward orbit of $x$ under multiplication by $p$ modulo one is eventually periodic.

| $p$ | $x$ | distinct residues | eventual fixed residue |
|---:|---:|---:|---:|
| 3 | $423213851608691/2701703435345984178$ | 39 | $1/2$ |
| 53 | $21/146068$ | 3 | $21/52$ |
| 227 | $3/226$ | 1 | $3/226$ |

For a fixed $J$, right continuity of the finitely many floor functions gives a small open interval immediately to the right of $x$ on which the same floors hold for $0\le j<J$. Take a closed subinterval of positive length in its interior, separated from every discontinuity and still satisfying $1<dy<p$. Along any arithmetic progression of exponents, $a\log_p 2$ is dense modulo one, since $2$ and $p$ are distinct primes. Hence infinitely often $q/p^n$ lies in this subinterval for some integer $n$ tending to infinity. Then $p^n\le N<p^{n+1}$. For division at the layer $p^{n-j}$ the shifted part and the total contribute additional real offsets $p^{j-n}$ in their floors. For large $n$ these offsets preserve every floor by the separation margin. The vanishing of $F$ proves the asserted upper-layer tests.

Run `python3 check_binary_local_obstruction.py` to verify the seeds, the Wieferich scan through one million, and the example at cofactor 1903. The program constructs full progressions for all primes through 200 at depth 32 and for all primes through 4096 at depth 12. It checks actual row residues at two positive exponents in each resulting progression. At the latter cutoff the extra tested weight primes are 19 and 229, each dividing $N$ exactly once.

## Prime exponents with individually prescribed first carries

The common prime divisors can themselves be located on this ray. Fix $L,J,m\ge1$ and integers $H_3,H_{53},H_{227}\ge\max(45,L+3)$. Choose any one of the three bases for a leading window of $J$ layers. There are **arbitrarily long arithmetic progressions of prime exponents** with all of the following properties at every term.

- The full conjecture holds at $N=12031\cdot2^a+1$. For the displayed triple $E,A,C$, all of $3,53,227$ are common divisors, each of exact valuation one in the common gcd.
- For each $p=3,53,227$ the **first** carry in $A$ occurs precisely at $p^{H_p}$. At each of $p^{H_p+8},\ldots,p^{H_p+7+m}$ the factorial-floor contribution is two. In the chosen leading base, all $J$ highest layers of $A$ have zero contribution.
- The lower-window assertions above still hold simultaneously at every tested prime below $6.7\times10^{15}$, or for any admissible finite set as in the theorem.

Consequently the three first carries can be placed at the same arbitrarily prescribed index $H\ge\max(45,L+3)$ in their respective bases. The positive conclusion for **all** triples at those weights follows from the bounded-cofactor classification. Any counterexample at 12031 would have the specific row $A$ and require it to avoid 3. Since 12031 is the only remaining odd cofactor through 65535 at prime exponents, each term of these progressions also proves the conjecture simultaneously at every weight $1+d2^a$ for odd $d\le65535$. [Necessary row at 12031](binary_cofactor_rays.md)

**Prescribing the middle digits.** Use $b=265327$ in the initial progression with the same $D=22211280$. Now $\gcd(b,D)=1$. In the order $5,7,13,113,1093,3511$ the six protected residues of $12031\cdot2^b+1$ are $4,4,2,84,1085,3254$. The three local witnesses still match their initial moduli. After matching through at least $p^6$, the exponent residues modulo $3,53,227$ are respectively $1,31,3$, all units.

For the ternary prime put $t=H_3-2$, and at the other two primes put $t=H_p-1$. Impose
$$2^a\equiv t_p-p^t\pmod {p^{H_p+8+m}}.$$
This retains the witness's low digits and its nonzero exponent residue modulo $p$. Each $dt_p+1<p^t$ (at 3 already $dt_3+1<3^{43}$). The lower block therefore adds exactly without carries. At $p^{t+i}$, throughout the imposed range $i\ge1$, the contribution is exactly
$$\left\lfloor\frac{\sum_{s\in\{681,9222,2128\}}((-s)\bmod p^i)}{p^i}\right\rfloor.$$
Indeed a nonnegative representative for the scalar is $t_p+p^t(p^i-1)$, and the lower block has total below $p^t$. At 3 the contributions for $i=1,2$ are 0 and 1. At 53 and 227 the contribution at $i=1$ is 1. Whenever $p^i>d$ the three residues sum to $3p^i-d$, giving contribution two. This proves all the claimed middle layers. The endpoint has valuation one at each prime. The equal parts of $C$ give a ternary carry at the first occupied digit. At 53 and 227 they are units whose three equal residues sum to more than one and are 1 modulo $p$ in total.

**Keeping a primitive progression.** When adjoining a further prime $r$, extend the current exponent class to a unit class modulo $\operatorname{lcm}(M,r-1)$. This is always possible by CRT, preserving the previous class and choosing a unit at new prime factors. If $r$ divides $N$ on this progression, use either $N/3\equiv r$ or $N/3\equiv 2r\pmod{r^{L+1}}$. Both targets give no carry through that layer, since $r\ge11$. At most one gives an exponent divisible by $r$: at the lift from $r$ to $r^2$ the two targets give distinct digits, and the exponent step is a unit modulo $r$. Choose the unit lift. Further lifts preserve it. This constructs a primitive progression of exponents with every prescribed finite congruence.

Primes in this progression have positive relative density among all primes. For the chosen leading base $p$, the values $a\log_p2$ are uniformly distributed modulo one over the primes in the progression, by Vinogradov's theorem. Selecting the interior interval from the leading-window proof leaves a set of primes of positive relative density. Green and Tao's Theorem 1.2 gives arbitrarily long arithmetic progressions inside that set. [Vinogradov, distribution on primes in an arithmetic progression](https://www.mathnet.ru/php/archive.phtml?jrnid=im&option_lang=eng&paperid=3453&wshow=paper) · [Green and Tao, Theorem 1.2](https://annals.math.princeton.edu/wp-content/uploads/annals-v167-n2-p03.pdf)

The same verifier constructs primitive progressions both with first layers $(45,45,45)$ and with $(60,53,71)$. In the second case it checks the low layers at every prime through 4096, the first carry at each of the three prescribed layers, and eight consecutive double contributions in each base, directly from the scalar residues.

### An explicit prime exponent

The following 1,694-bit prime $a$ is written in hexadecimal, most significant digits first. Join the lines to recover the integer.

```text
32e06c8f77e5efab07dc1f49fa23b1377d0fe98110ffa8ef122bd06a73fabad8d154710f
ae6f359ac0f456938941f05d18a85a32c0db0e571b768bcda0e94290f8e089bb580a9e79
284e20c68e9a4786aa80dce8bb4e2adcc06538425428b7dc6ad03b2be91ac6284d7cff52
faa31cf040c6c78784f86fc291aa8f76ddcf1c1545edf9d12e8d6088de903bb41be13ec1
2fb8ff1b6489beb648c0d42e9c5cd93a16d66622afc387cc474715a99d9973359bdf05c0
f2d8944fbc4a24beed8240e114f2a8df43ab1845bd35019dd55b7a6b72fd6e8f
```

Here the first carry of $A$ occurs at layer 45 in each of the three bases. Layer 53 has contribution two in each base. The four highest ternary layers have contribution zero. The low-window construction with $L=2$ is verified through prime 19 for this particular exponent.

For a primality certificate put $F=11^{246}$. Direct modular exponentiation gives
$$F\mid a-1,\qquad F^2>a,\qquad 2^{a-1}\equiv1\pmod a,\qquad
\gcd(2^{(a-1)/11}-1,a)=1.$$
For every prime divisor $r$ of $a$, the multiplicative order of 2 modulo $r$ divides $a-1$ and does not divide $(a-1)/11$. Hence $F$ divides that order and $r-1$. Every prime divisor of $a$ exceeds its square root, proving primality.

The verifier reconstructs the same integer from the exponent progression. It checks the low and middle layers by modular powers. For the four highest layers it encloses $2^a$ between integer mantissas times a common power of 3. Each binary exponent digit squares these bounds and multiplies by 2 if the bit is one. Rescaling rounds both endpoints outward. The upper and lower quotient bounds agree for all three parts and the total at each tested layer, including the unit offsets in the third part and the total. Thus computing the top layers needs no floating-point logarithm and never constructs $2^a$.
