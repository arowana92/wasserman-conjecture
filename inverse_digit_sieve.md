# A digit bound when a part contains most of the predecessor

Call a weight exceptional if some three positive trinomial rows have coefficients with gcd one. Fix a prime $p$ and define
$$
 \Lambda_p=\frac{1}{2\sin(\pi/(4p+2))},\qquad \beta_p=\log_p\Lambda_p.
$$

**Theorem.** In any arithmetic progression of $H$ positive weights with common difference coprime to $p$, at most $O_p(H^{\beta_p})$ are exceptional, satisfy $p\mid N-1$, and have at most three distinct prime divisors of $N-1$. The bound is uniform in the starting point and the common difference. In ordinary intervals of $H$ consecutive integers the bound improves to $o_p(H^{\beta_p})$, uniformly in the starting point.

For $p=2$ this counts all odd exceptional weights with at most three predecessor primes, with exponent $\log_2((1+\sqrt5)/2)=0.6942419\ldots$. For each fixed odd $p$, both bounds also apply to exceptional weights with $2\parallel N-1$, exactly three distinct odd predecessor primes, and $p$ one of those primes.

For an arbitrary positive common difference $d$ and first term $A$, the progression bound becomes $O_p(H^{\beta_p}p^{(1-\beta_p)a})$, where $a=v_p(\gcd(d,A-1))$. Thus even an arbitrarily large power of $p$ in the step costs only a constant when the full $p$-power of $A-1$ stays fixed.

## Encoding the digits of the cofactor

Write $M=N-1=WT$, with $W=p^e$, $e\ge1$, and $p\nmid T$. If a row avoids $p$ (its coefficient is coprime to $p$), Kummer's criterion says that its parts add without carries in base $p$. Exactly one part $B$ is $1$ modulo $W$. The other two are multiples of $W$. Put
$$s=\frac{T}{\gcd(T,B)},\quad t=\frac Ts,\quad B=kt,\quad x=\frac{B-1}{W}=\frac{kT-s}{sW}.$$
Here $1\le k<sW$ since $B<M$, and $kT\equiv s\bmod W$. The base-$p$ digits $x_j$ are at most the corresponding digits $T_j$ because the digits of $B$ fit within those of $N=WT+1$.

Fix $s,e$ and one of the $s$ possible lifts of $k$ from its least positive residue modulo $W$. For $L\ge e$, the first $e$ digits of $T$ determine that residue. Given the first $j+e$ digits, varying digit $j+e$ permutes digit $j$ of $x$. In fact a change of $a$ at place $p^{j+e}$ changes $x$ by $(k/s)ap^j$ modulo $p^{j+1}$, where $s$ and $k$ are units. This congruence defines the output digit even for prefixes considered before testing $s\mid T$. Exactly $T_j+1$ next digits are allowed.

Keep the first $e$ digits unchanged. Encode each later true digit by mapping the sorted allowed digits in order onto the largest $T_j+1$ members of $\{0,\ldots,p-1\}$. The encoded digit is at least the true digit and at least $p-1-T_j$. Thus encoded digits at distance $e$ have sum at least $p-1$. The map is injective because the initial digits determine $k$ for the chosen lift, and each successive rank recovers the next true digit from the known prefix.

Let $A_{ij}=1$ if $i+j\ge p-1$ and zero otherwise, for $0\le i,j<p$. Its largest eigenvalue is $\Lambda_p$. To see this put $\theta=\pi/(2p+1)$ and $v_j=\sin((j+1)\theta)$. With $v_{-1}=0$, consecutive row sums of $Av$ differ by $v_{p-1-j}=\cos((j+1/2)\theta)=\Lambda_p(v_j-v_{j-1})$, including the first row compared with zero. The vector $v$ is positive, and $A$ is symmetric and nonnegative, so its operator norm is $\Lambda_p$. Each chain of $m\ge1$ digit positions spaced $e$ apart admits at most $\mathbf1^\mathsf{T} A^{m-1}\mathbf1\le p\Lambda_p^{m-1}$ words. Multiplying over the chains gives, per lift, at most
$$ \Lambda_p^L(p/\Lambda_p)^e \tag{1}$$
possible length-$L$ prefixes. This also covers $L<e$ since then the trivial count $p^L$ suffices and $\Lambda_p\le p$. The all-ones Rayleigh quotient gives $\Lambda_p\ge(p+1)/2>\sqrt p$. Also $\Lambda_p<p$ by the explicit formula or strict row-sum comparison. In particular $\beta_p<1$ and $p/\Lambda_p^2<1$.

## Progressions and the factor $s$

Consider $H$ successive terms of $N=A+di$, with $d\ge1$ and $p\nmid d$. For a fixed $s\le H$ coprime to $p$, let $g=\gcd(d,s)$ and $n=\lceil\log_p(Hg/s)\rceil$. The congruence $sp^e\mid N-1$ is either empty or cuts out one class of indices modulo $sp^e/g$. For $e<n$ the corresponding values of $T$ before imposing the exact valuation have step $sd/g$, a $p$-adic unit. There are at most $\lceil Hg/(sp^e)\rceil\le p^{n-e}$ such indices. Their $T$ values therefore have distinct residues modulo $p^{n-e}$. Applying (1) with $L=n-e$ and summing over the lifts bounds the contribution at $e$ by
$$s\Lambda_p^n(p/\Lambda_p^2)^e. \tag{2}$$
All $e\ge n$ together contribute at most one index, by the single congruence $sp^n\mid N-1$. Summation of (2) gives $O_p(s(Hg/s)^{\beta_p})$ for the existence of a $p$-avoiding row with this $s$. In particular the counts for any $s\le Y\le H$ are $O_p(H^{\beta_p}Y^{2-\beta_p})$ on ordinary intervals and $O_p(H^{\beta_p}Y^2)$ on progressions with $p\nmid d$.

With at most two predecessor primes the conjecture already holds. With exactly three, the endpoint reduction forces one row to avoid all three bases, with three different positions for the parts congruent to one. Its $B$ is divisible by $T$, giving $s=1$. When $2\parallel M$ and $M$ has three odd bases, all three odd bases must again go to one row. Its parts are $U$ times positive integers of sum $N/U$ with $U\ge3$ coprime to $M$. Each normalized part is less than $M/2$, so no position can be zero modulo all three odd powers. Their three exceptional positions are distinct. At an odd $p$ in this row $s$ divides 2. This proves the progression bounds. [Endpoint reduction](all_weights_endpoint_reduction.md) · [Three predecessor primes](endpoint_factors.md) · [Four predecessor primes](two_prime_powers.md)

For the extension to $p\mid d$, restrict to $s\mid N-1$ first. On the remaining indices $(N-1)/s=B+Di$ has step $D=d/\gcd(d,s)$ and length $H'\le\lceil H\gcd(d,s)/s\rceil$. Put $f=v_p(d)$ and $n=\lceil\log_p H'\rceil$, omitting empty progressions. If $a=f$, divide $B+Di$ by $p^f$ to obtain a progression of unit step modulo $p$. At valuation $e=f+b$, $b\ge0$, the candidates for $T$ occupy distinct prefixes of length $n-b$ when $b<n$, by the same congruence argument. Formula (1) contributes at most $s\Lambda_p^n(p/\Lambda_p)^f(p/\Lambda_p^2)^b$. The entire tail $b\ge n$ contributes at most one.

If $a<f$, the valuation $e=a$ is fixed. The case $a=0$ contributes nothing. Otherwise all the values of $T=(N-1)/p^a$ have the same first $f-a$ digits. The next $n$ digits distinguish the candidates, since after deleting the common low digits their step is a $p$-adic unit. Apply the encoding inside this new block. Keep its first $e$ digits free, at which point the earlier fixed digits and these digits determine $k$ for each lift. Varying a later digit permutes the corresponding digit of $kt$, which must fit beneath the true $T$ digit $e$ places earlier. This gives at most $s\Lambda_p^n(p/\Lambda_p)^a$ candidates, also for $n<e$ by the trivial count. In either case the fixed-$s$ bound has acquired at most the factor $(p/\Lambda_p)^a=p^{(1-\beta_p)a}$. This proves the extension after summing over the bounded list of $s$ in the theorem.

## A second prime on ordinary intervals

In either application the row containing $B$ avoids every prime of its assigned scale $U\ge3$, a divisor of $N$ dividing each of its parts. Since $\gcd(t,N)=1$, we have $U\mid k$. Hence there is a prime $r\mid k$ dividing $U$. It is coprime to $sW$. For fixed $s,e,k,r$ both $p$ and $r$ must allow the carry-free cut $kt,(sW-k)t+1$. We bound the possible blocks of $h$ consecutive base digits of $t$ at these two primes, from any digit position.

At $p$, fix the incoming multiplication carries for $kt$ and $st$ at the start of the block. There are finitely many possibilities depending only on $k,s$. For these carries the $h$ input digits correspond bijectively to the $h$ output digits of $st$, since $s$ is a unit. At positions $j\ge e$ within the block the digit of $kt$ must be at most the digit of $st$ at $j-e$, by $N=Ws t+1$. Once the lower digits are fixed, varying the current $st$ digit permutes the current $kt$ digit, with multiplier $k/s\bmod p$. Thus after leaving the first $e$ digits free, the same rank encoding gives $O_{p,e,s,k}(\Lambda_p^h)$ possible blocks.

At $r$, fix incoming multiplication carries for $kt$ and $(sW-k)t+1$ and write $b=v_r(k)\ge1$. In the next $b+1$ digits at least one input word is forbidden. The first input digit can force digit $b$ of the first expression to equal $r-1$, independently of later inputs. Fix any intermediate digits. Since $sW-k$ is a unit, the last input digit can then force the other expression's digit $b$ also to equal $r-1$. The collision violates the cut test. Iterating on disjoint blocks and allowing all finitely many carry states gives $O_{s,e,k,r}(r^{\sigma h})$ possible length-$h$ blocks, with $\sigma=\log_r(r^{b+1}-1)/(b+1)<1$.

Use these blocks as digit alphabets in bases $p^h,r^h$, adjoining zero if needed. Their compact digit sets have dimensions at most $\beta_p+O(1/h)$ and $\sigma+O(1/h)$. The bases are multiplicatively independent. The uniform affine intersection bound in Corollary 1.17 of Corso and Shmerkin gives $O(H^{\max(0,\beta_p+\sigma-1)+\varepsilon})$ possible inputs in any interval of length $H$, after choosing $h$ large in terms of $\varepsilon$ and the fixed parameters. Indeed in each base truncate at the least block-base power at least $H$. The higher quotient takes at most two values, the scaled low digits belong to the compact digit set, and the resulting dilation is between 1 and that block base. The intersection estimate at scale $1/H$ then counts the separated inputs. An interval of $H$ weights also gives an interval of at most $H$ inputs $t$. [Corso and Shmerkin, Corollary 1.17](https://arxiv.org/html/2409.04608v1#S1.SS6) · [Details of the interval rescaling](five_seven_pairs.md#from-the-digit-sets-to-moving-intervals)

The last exponent can be chosen strictly below $\beta_p$ for every fixed $s,e,k,r$. Fix $E$. There are only finitely many such choices for $s\in\{1,2\}$ and $e\le E$. The remaining $e>E$ contribute $O_p(1+H^{\beta_p}\sum_{e>E}(p/\Lambda_p^2)^e)$ by (2). Divide by $H^{\beta_p}$, first let $H$ and then $E$ tend to infinity, and obtain the claimed uniform little-$o$ estimate.

## Finite residue certificates

For $s=1$, take a modulus $p^m$ and write $L=m-e$. When $L\ge e$ the known low digits of $T$ determine $k=T^{-1}\bmod W$, taken in $[1,W-1]$. Discard $k<3$ since $U\mid k$. Discard $W\le4$ by the proved small-power case. Test the known digits of $kT$ against those of $N$ through position $L-1$. When $L<e$ retain every unit prefix after the small-power test. Retain the class $N\equiv1\bmod p^m$ to allow larger valuations. [Small predecessor powers](three_factors_small_power.md)

| Modulus | Classes with $p\mid N-1$ | Classes surviving | Classes excluded |
|---|---:|---:|---:|
| $2^{24}$ | 8,388,608 | 37,939 | 8,350,669 |
| $3^{13}$ | 531,441 | 8,725 | 522,716 |
| $5^9$ | 390,625 | 10,733 | 379,892 |

Every excluded class satisfies the full conjecture at all its weights with at most three distinct predecessor primes. The surviving counts are upper bounds. Run `python3 check_inverse_digit_encoding.py` to reproduce the table by direct modular enumeration and independently by growing allowed digit prefixes. The checker also compares 500,338 finite prefixes with direct modular digit comparisons, recovers every accepted prefix, and tests the encoding on 2,478 actual rows using exact binomial products.

At larger moduli the path count gives a certificate without enumerating every residue. Write $P_\ell(i)$ for the number of length-$\ell$ chains with first digit $i$, so $P_1(i)=1$ and $P_{\ell+1}(i)=\sum_{j\ge p-1-i}P_\ell(j)$. For $L\ge e$ the chain lengths are $\ell_i=1+\lfloor(L-1-i)/e\rfloor$ for $0\le i<e$. Multiply the chain counts, summing over all first digits with digit zero of $T$ nonzero. Subtract the products for the specific initial $e$ digits giving $k=1$ or $2$, when these are units. These initial digits are preserved by the injection. Sum over $e$ with the same boundary cases as in the table. This gives the following rigorous upper bounds on classes that can contain exceptions under the three-prime hypothesis.

| Modulus | Upper bound on surviving classes | Upper fraction among classes $N\equiv1\bmod p$ |
|---|---:|---:|
| $2^{128}$ | 450,309,764,530,657,106,546,034,206 | $2.647\times10^{-12}$ |
| $3^{81}$ | 13,162,362,795,641,191,516,829,074,211 | $8.905\times10^{-11}$ |
| $5^{55}$ | 362,545,042,571,185,944,680,949,795,777 | $6.532\times10^{-9}$ |

The checker uses integer arithmetic and verifies the recurrence by directly listing encoded words at three smaller moduli. A separate binary calculation reproduces the bound at $2^{128}$ using Fibonacci numbers. With $F_1=F_2=1$, a chain of length $\ell$ starting at 0 or 1 admits $F_\ell$ or $F_{\ell+1}$ words respectively. Thus fewer than three in a trillion odd residue classes modulo $2^{128}$ can contain even one counterexample with at most three predecessor primes, independently of how large representatives of the class are.

## Transferring the sieve to powers

For each $m\ge1$ the candidate classes modulo $p^m$ form a set of size $O_p(p^{m\beta_p})$. This follows by summing (1) with $L=m-e$, $s=1$, over $1\le e<m$, and adding the single class $N\equiv1\bmod p^m$. Thus the savings transfer through maps with bounded fibers on units. For example, fix integers $c,b$ with $p\nmid c(1-b)$, a positive integer $k$, and take $N=c u^k+b\ge3$ on any $H$-term progression of integer inputs $u$ of step coprime to $p$. The number of exceptional inputs with $p\mid N-1$ and at most three distinct predecessor primes is $O_p(k H^{\beta_p})$, uniformly in the coefficients, initial input and step. Take the least $p^m\ge H$ for $H\ge2$. The inputs have distinct residues. Every candidate weight gives a unit value for $u^k\bmod p^m$, with at most $2k$ possible $u$ by the cyclic unit group for odd $p$ and its cyclic subgroup of index at most two for $p=2$.

When $\gcd(k,p(p-1))=1$, exponentiation permutes all units and the factor $k$ drops out. In particular the binary table applies with exactly the same residue counts to $c u^k+b$ for each odd exponent $k$, every odd $c$ and every even $b$, with no upper bound on $k$. The computable surviving classes of $u$ are the preimages of the original sieve.
