# Exponent counts and two new families at the 4096th power

Let $P=65537$, let $q$ be any prime, and put $N=1+P\,2^b q^{4096j}$ for positive integers $b,j$. The full three-row conjecture holds throughout both families below.

- $169\mid b$ and $169\mid j$, with both exponents unbounded.
- $b\le10^6$ and $13^{20}\mid j$, with no bound on $q$ or $j$.

For arbitrary $j$ and fixed $q,j$, every interval of $H$ consecutive binary exponents $b$ contains at most $O(H^{0.759})$ possible exceptions. The implied constant is absolute, independent of $q,j$ and of the starting exponent. There is also an explicit simultaneous class count below. These bounds concern a necessary condition for failure.

## The cut at two primes

At a failure the [half-exponent reduction](fermat_half_exponent.md) gives $4\mid b$, a row avoiding $13$ and $P$, and one of its parts $xW$, where
$$
 W=2^b q^{4096j},\qquad x=49153,\qquad y=P-x=16384.
$$
It gives $xW\equiv1\bmod P$ and $13\mid N$. Consequently $W\equiv4\bmod P$ and $W\equiv3\bmod13$. The binomial coefficient with parts $xW,yW+1$ divides the trinomial coefficient of that row, so the two parts must add without carries in both bases. If $q=2,P$, the two-prime predecessor result already applies. If $q=13$, then $13\mid N$ is impossible. [Endpoint theorem](all_weights_endpoint_reduction.md)

Let $C_L$ count scale prefixes passing the first $L+1$ base-$P$ places with units digit 4, and $D_K$ count those passing the first $K+1$ base-13 places with units digit 3. $C_0=D_0=1$. The [nested intervals](predecessor_digit_matrix.md) give $C_L=(A^L\mathbf1)_5$. In base 13 start with the two multiplication carries $(a,c)=(11343,3781)$ obtained from input digit 3 and the offset $+1$ in the second part. At each next input digit $d$ keep it exactly when
$$
 ((xd+a)\bmod13)+((yd+c)\bmod13)<13 .
$$
The child state is $(\lfloor(xd+a)/13\rfloor,\lfloor(yd+c)/13\rfloor)$. This is ordinary long multiplication, so it counts precisely the clear prefixes, including ones ending in zero digits. The graph reached from the starting pair has 65,536 states and 458,752 digit transitions. Its adjacency matrix $B$ counts transitions with multiplicity. Then $D_K=(B^K\mathbf1)_{\rm start}$.

The first predecessor counts are $C_1=5$, $C_2=196596$ and $C_3=6442549034$. At 13 one obtains $D_1=6$, $D_2=40$, $D_4=1956$ and $D_{10}=231470218$. For an all-depth bound the checker verifies the coordinatewise integer inequalities
$$
 8A^9\mathbf1\le262145 A^8\mathbf1,\qquad
 10^6 B^{20}\mathbf1\le7000001 B^{19}\mathbf1.
$$
The vectors on the right before multiplication by their scalar have positive entries: every base-$P$ state admits a zero digit, and at 13 the digit of $xW$ is predetermined since $13\mid x$, while the other digit runs through all 13 possibilities. Applying the nonnegative matrices repeatedly therefore gives $C_L=O((262145/8)^L)$ and $D_K=O((7000001/10^6)^K)$. The constants refer only to these fixed matrices.

## Exact number of classes passing the two prefixes

Write $T=q^{4096j}$. For each fixed unit $T$ coprime to $13P$, the conditions $2^bT=4\bmod P$ and $2^bT=3\bmod13$ either disagree or specify one class modulo 96. Here the orders of 2 modulo $P,13$ are respectively 32 and 12. If a starting congruence has no solution, the count is zero.

In the compatible case, exactly $C_LD_K$ classes modulo $96P^L13^K$ pass the two prefix tests. Indeed $v_P(2^{32}-1)=1$ and $2^{12}\equiv40\bmod169$. Thus $v_P(2^{96}-1)=v_{13}(2^{96}-1)=1$. For either odd prime $r=P,13$, binomial expansion shows that raising $1+ru$, with $r\nmid u$, to $r^t$ gives valuation $t+1$ after subtraction of one. Powers of $2^{96}$ therefore permute all lifts of the chosen low digit in each base. The two periods $P^L$ and $13^K$ are coprime, so every pair of lifts occurs once. At $(L,K)=(2,10)$ the count is
$$
 \frac{45\,506\,118\,977\,928}{56\,843\,115\,215\,370\,210\,650\,976}
       <8.01\cdot10^{-10}
$$
of all exponent classes. The joint count describes just these two tests; additional conditions on a counterexample can reduce it further.

Taking $L=0$ and choosing the least $K\ge0$ with $96\cdot13^K\ge H$ bounds the number of possible exponents in the interval by $D_K$. Integer comparison gives $(7000001/10^6)^{1000}<13^{759}$. This proves the uniform interval bound. At $P$ alone the same argument uses period $32P^L$ and gives exponent less than $15/16$ since $(262145/8)^{16}<P^{15}$.

## A synchronized power at 13

For $13^s\mid j$ and $q\ne13$, $T^3\equiv1\bmod 13^{s+1}$ by Euler's theorem, since $12\cdot13^s\mid3\cdot4096j$. For each unit $a$ put $\tau(a)=a^{13^s}\bmod13^{s+1}$. Its order divides 12 and it is congruent to $a\bmod13$. Reduction modulo 13 is injective on these roots of unity, since the kernel has order $13^s$. In particular $T=\tau(T)$ at this modulus.

Consequently the necessary low digit $W\equiv3\bmod13$ forces the unique lift
$$
 W=\tau(3)\big(2/\tau(2)\big)^b\bmod13^{s+1}.
$$
As $b$ runs modulo $13^s$, this permutes all lifts of 3, by the order calculation at 13. Exactly $D_s$ classes pass the 13 test, simultaneously for all $q$ and $13^s\mid j$, in the sense that every failure must lie in this one set independent of $q,j$.

For $169\mid b,j$, take $s=2$. The lift is $W=1160\bmod2197$. Then $xW\bmod2197=936$ and $(yW+1)\bmod2197=1391$; their sum exceeds 2197. This proves the first family.

For $13^{20}\mid j$, test all $4\mid b$ from 72 through one million. Modulo $13^{21}$ the forced lift starts at $\tau(3)$ for $b=0$ and is multiplied by $16/\tau(16)$ at each step of four. Only $b=50088,248052$ have no carry through the first twenty places, and both carry at place twenty-one. The [previous finite theorem](fermat_half_exponent.md) supplies $b<72$ and all $4\nmid b$. This proves the second family.

Run **python3 check_fermat_exponent_digits.py** from the repository root. Besides both coordinatewise bounds, the program checks the counts by direct Legendre quotients in small bases and through six places at 13. It independently verifies the million-exponent range by enumerating $b/4\bmod13^3$, computing modular powers directly for its surviving classes, and testing individual digits.
