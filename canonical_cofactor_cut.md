# One binomial cut for two static prime powers

Let $a=p^i$, $b=q^j$ be powers of distinct primes, $i,j\ge1$. Set $D=ab$, $n=D+1$. Let $c$ be the least positive inverse of $b$ modulo $a$, and $d$ the least positive inverse of $a$ modulo $b$. Then
$$
 bc+ad=n,\qquad 1\le c<a,\quad 1\le d<b.
$$
Indeed the sum is one modulo $ab$ and lies strictly between 1 and $2ab+1$.

**Theorem.** Put $g=\gcd(c,d)$. If the conjecture fails even once in the family
$$
 N=1+Dt^e,\qquad t\ge1,\quad\omega(t)\le2,\quad e\ge1,\quad L_{n^3}\mid e,
$$
then all three conditions below hold:
$$
 g\ge3,\qquad p\nmid\binom b d,\qquad q\nmid\binom a c . \tag{1}
$$
Here $L_v$ is the exponent of the unit group modulo $\operatorname{lcm}(1,\ldots,v)$. The same conclusion holds with independently chosen positive multiples of $L_{n^3}$ as the two varying exponents, including shared bases with $D$. Equivalently, the single binomial $\binom n{bc}=\binom n{ad}$ must avoid both static primes and its two parts must have gcd at least three. Indeed $g=\gcd(bc,ad)$ since their sum $ab+1$ is coprime to $ab$. At $p$ the cut at $ad$ removes $i$ trailing digits of the variable portion of the weight, leaving precisely the cut at $d$ of $b$. The corresponding cut at $q$ is at $c$ of $a$.

The period can also be replaced by an explicit integer smaller than $D^4(D+1)^2$, simultaneously for all varying prime bases. A product congruence gives families with both varying exponents equal to one. [Polynomial period and ordinary prime factors](static_product_congruence.md)

In particular **$b\equiv1\bmod a$ or $2b\equiv1\bmod a$ suffices to certify the entire family**. A distribution consequence is that among $D\le X$ having exactly two distinct prime divisors, at most $O(X/\log X)$ can admit any failure in the family. Its proportion among these cofactors tends to zero. A second consequence for $b\equiv-1\bmod a$ appears [below](#one-residue-class-has-polynomially-many-digit-patterns).

## Why both primes use the same cut

Use the [reference-weight reduction](two_power_seed.md) and its [two-static-power specialization](static_power_transfer.md). In a grouped assignment let $B_0=(x,n-x)$ be the cut, $3\le U\le x<ab$ the unitary divisor of $n$ dividing $x$, and $C_0=(n/U)(f_1,f_2,f_3)$ the positive triple with $\sum f_i=U$. At each static prime at least one of $B_0,C_0$ avoids that prime.

The companion $C_0$ cannot avoid both: two of its entries must be zero modulo $a$ and two zero modulo $b$. The intersection contains an entry divisible by $ab$, impossible for a positive triple of total $ab+1$. Thus $B_0$ avoids at least one prime, say $q$. The endpoints $x,n-x$ are at least two since $U\ge3$ divides both. Write $x=mb+\delta$ with $\delta=0,1$ and $1\le m<a$, using its residue at $b$. Then $U$ divides $ax-mn=a\delta-m$, which has positive absolute value less than $a$. This forces $U<a$. Were $C_0$ to avoid $p$, two of its entries would be multiples of $a$, hence so would their normalized parts $f_i$ since $(n/U,a)=1$. Their sum would already exceed $U$. Thus $B_0$ avoids both static primes. Among its parts the multiple of $b$ must be one modulo $a$: an entry divisible by both would be at least $ab$. It is therefore $bc$, and the other is $ad$. In particular $U\mid g$, giving $g\ge3$.

In a split assignment the static primes are taken by the same row and divide distinct multipliers there. Shared prime bases rule out a split by absorption. The [cross-prime equality](static_power_transfer.md#two-different-varying-bases-in-different-rows) gives primitive multipliers $X,Y$ with $a\mid X$, $b\mid Y$ and reference parts $kX,kY$ summing to $n$, carry-free at both $p$ and $q$. Necessarily $kX=ad$, $kY=bc$. Primitivity gives $k=g$, while the strict product-of-sums bound gives $k\ge3$. This proves (1). It also specifies the primitive pair $(ad/g,bc/g)$ before considering any weight primes.

## Counting all cofactors with two prime divisors

First take $a=p<b=q$ prime. Fix the units digit $b_0$ of $b$ at $p$, so $c$ is fixed. Write $b=b_0+pt$ and
$$
 d=(p-c)t+((p-c)b_0+1)/p .
$$
For the cut to pass, every digit of $d$ must be at most the corresponding digit of $b$. Given the preceding digits of $t$, the current digit of $t$ permutes the current digit of $d$, since $p-c$ is a unit modulo $p$. Consider two consecutive digits of $t$. If the first is $u$, the second digit of $d$ at this pair must be at most $u$ (the corresponding digit of $b$). Even allowing all first digits, at most $\sum_{u=0}^{p-1}(u+1)=p(p+1)/2$ ordered pairs pass. This bound holds after every prefix. Among the unit classes of $b\bmod p^L$ the fraction that can pass is therefore at most
$$
 ((p+1)/(2p))^{\lfloor(L-1)/2\rfloor}\quad(L\ge1). \tag{2}
$$

For $p\le\sqrt X$ put $Y=X/p$, and take $L=\lfloor\log Y/(2\log p)\rfloor$. For $L\ge1$, Brun–Titchmarsh at modulus $p^L\le\sqrt Y$ bounds the number of primes up to $Y$ in each unit class by $4Y/(\varphi(p^L)\log Y)$. Multiply by the number of allowed classes from (2). Since $(p+1)/(2p)\le3/4$ and $\log Y\ge(\log X)/2$, the resulting bound is
$$
 \ll \frac{X}{p\log X}\exp\big(-\eta\log X/\log p\big)
$$
for some absolute $\eta>0$. For $L=0$ the same bound follows from $\pi(Y)\ll Y/\log Y$ since $\log X/\log p<3$. [Montgomery and Vaughan, *The large sieve*, equation (1.12)](https://personal.science.psu.edu/rcv4/personal/Publications/large_sieve.pdf)

The sum over $p$ is $O(X/\log X)$. Indeed on each range $2^s\le\log X/\log p<2^{s+1}$ the prime reciprocal sum is $O(1)$ by partial summation from $\pi(z)\ll z/\log z$, and $\sum_s\exp(-\eta2^s)$ converges. The number of cofactors with some proper power $u=\ell^f$, $f\ge2$, also has size $O(X/\log X)$. For $u\le\sqrt X$ use the bound $O(X/(u\log X))$ on the number of prime powers up to $X/u$ and sum the convergent reciprocals of proper prime powers. For $u>\sqrt X$, sum $O(\sqrt{X/b}\log X)$ possible proper powers over all integers $b\le\sqrt X$, for $O(X^{3/4}\log X)$. Finally, the total number of two-prime cofactors is $\gg X\log\log X/\log X$: for each $p\le X^{1/4}$ take primes $q$ between $X/(2p)$ and $X/p$ and sum the prime number theorem lower bound. This proves the relative-density assertion.

## One residue class has polynomially many digit patterns

Suppose $b\equiv-1\bmod a$. Then $c=a-1$, $d=(b+1)/a$, and $g=\gcd(a-1,b+1)$. Write $b+1=p^h u$ with $h\ge i$ and $p\nmid u$. Let $\preceq_p$ denote digitwise comparison in base $p$. Kummer at $p$ requires $d\preceq_p b$, equivalently
$$
 u\preceq_p p^i u-1.
$$
In fact $d=p^{h-i}u$ and $\lfloor b/p^{h-i}\rfloor=p^i u-1$. The low $i$ digits of the right side are all $p-1$ and its upper word is $u-1$. Thus $u_i<u_0$ and $u_j\le u_{j-i}$ for $j>i$, indexing digits from zero. Each of the $i$ chains of digits a distance $i$ apart is nonincreasing. In the chain starting at zero the digits after the first are at most $p-2$. For $b<p^L$, each nonincreasing string over digits $0,\ldots,s$ is determined by at most $s$ positive-digit counts, with $O(L+1)$ choices each. Counting $h$ as one more parameter gives only $O_a((1+L)^{i(p-1)})$ possible $b$. At $a=4$ the [exact formula](four_times_sparse_powers.md) also imposes the gcd condition.

**For fixed $p,i,q$ there is an effective threshold beyond which every $b=q^j\equiv-1\bmod a$ certifies the family.** In particular when $q\equiv-1\bmod a$ this covers every sufficiently large exponent $j$: even exponents give the already certified residue $+1$.

Here are details of effectivity. Regard $u$ as a word in base $a$ with digits $d_0,\ldots,d_s>0$ up to the last nonzero digit, in order from the units place. Coordinatewise monotonicity makes this sequence nonincreasing. Subtract $u$ from $au$, and put $d_{s+1}=0$. The nonzero differences $d_l-d_{l+1}$ give a representation
$$
 (a-1)q^j=\sum_{\nu=1}^T c_\nu p^{l_\nu}-d_0p^h-(a-1),\qquad
 l_1>\cdots>l_T\ge h+i,\quad \sum c_\nu=d_0<a .
$$
In particular $T\le a-1$ and $q^j+1<p^{l_1}\le a(q^j+1)$ from the top digit of $u$. Also $h=O_{p,q}(1+\log j)$ because $p^h\mid q^{2j}-1$, using the [prime-power order bound](prime_power_ray_density.md#counting-exponents).

For a prefix through $\nu$, write its sum as $p^{l_\nu}B_\nu$; then $B_\nu\le a p^{l_1-l_\nu}$. Divide by $(a-1)q^j$ and subtract one. The difference is nonzero since the prefix is divisible by $p$ and the denominator is a unit at $p$. Its absolute value is $O_a(p^{l_{\nu+1}-l_1})$, setting $l_{T+1}=h$, by the omitted terms and the correction of size $O_a(p^h)$. Apply Matveev to the positive rationals $p,q,B_\nu/(a-1)$ with exponents $l_\nu,-j,1$, omitting the last term when it is one. Their heights are bounded by constants depending on $p,q,a$ for the first two, and $O_a(1+l_1-l_\nu)$ for the third. Since $l_1=O_{p,q,a}(1+j)$, this gives effectively
$$
 l_1-l_{\nu+1}\ll_{p,q,a}(1+l_1-l_\nu)\log(2+j).
$$
Iteration yields $l_1-h\ll_{p,q,a}(\log(2+j))^{a-1}$. Together with the bound on $h$ and $q^j<p^{l_1}$, this bounds $j$ effectively in terms of the three fixed parameters. [Matveev's theorem](https://www.mathnet.ru/eng/im314)

## Verification

Run **python3 check_canonical_cofactor_cut.py**. The verifier checks the two smaller binomials against the exact coefficient at the reference weight and compares the exclusions with the general grouped and split reference algorithms. Through $D=1{,}000{,}000$, 13,863 of 288,726 unordered pairs survive, so 274,863 cofactors certify all the indicated weights. Through $D=500$, the general reference algorithms independently confirm all 250 canonical exclusions. The verifier also checks the descending-digit criterion on 130,182 integers. [Verifier](check_canonical_cofactor_cut.py)
