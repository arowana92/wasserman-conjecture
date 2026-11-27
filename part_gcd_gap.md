# The gap between one part and the gcd of the other two

**Theorem.** Suppose three positive trinomial rows of weight $N$ have coefficient gcd one. One row is a permutation of the [endpoint](all_weights_endpoint_reduction.md). In either other row $(a,b,c)$, put $M=\gcd(b,c)$. Then
$$a=M \qquad\text{or}\qquad |a-M|>1024.$$
Every choice of the distinguished part is allowed. At most one of the two nonendpoint rows can realize the equality at even one part. The bound 1024 uses a finite computation below. Analytically, every nonzero such gap is $\gg N^{1/5}$.

More generally let $d,e$ be positive integers, and put $h=da-eM$. At any fixed rational slope $e/d$, a nonendpoint row in a counterexample must satisfy
$$h=0,\quad h=d,\qquad\text{or}\qquad |h|\gg_{d,e} N^{1/5}.$$
When $d,e\le8$ are coprime, every $|h|\le256$ other than $0,d$ is excluded at all weights. For $d=1$ the case $h=d$ is also impossible, since then $\gcd(a,M)=1$.

## A bound for arbitrary companions

Take arbitrary positive rows $A=(a,b,c)$ and $A'$ of common weight $N$ and write $G=\gcd(N(N-1),K_A,K_{A'})$. Put $x=b/M$, $y=c/M$, $r=x+y\ge2$, $H=|h|$, $E=|d-h|$, and assume $h\ne0,d$. If $e+hr\ne0$, then
$$ N-1\ \le\ \frac{G^3 e E^2 H^2\max(d,H)}{27d}
              +\frac{G^{3/2} E(E+H)}{d}\sqrt{\frac H{27}}. \tag{1}$$
In the exceptional case $e+hr=0$ there is a separate, smaller bound
$$ N-1\le \frac{G^{3/2}E e^2}{4\sqrt{27H}}. \tag{2}$$
Thus for fixed $d,e,h$ outside the two exact lines, the endpoint gcd grows at least as a constant times $N^{1/3}$, uniformly in the companion row.

As in the [retained-power argument](endpoint_pair_bound.md#retaining-partial-prime-powers), write $G=DJ$, let $U,S$ be the retained products from $N,N-1$ assigned to $A$ and let $T$ be the complementary predecessor product. Then $N-1=JST$, $DU\ge3$, and $T^2\le(DU/3)^3$. In particular $K=JT\le G^{3/2}(U/3)^{3/2}$. The weight product $U$ divides $a,b,c$, hence divides $M$ and $h$. Also $g=\gcd(M,a-1)=\gcd(M,N-1)$ divides $E$.

At the powers of $S$ where $a$ is one modulo the power, both $b,c$ vanish, so their product of powers divides $g$. At every other assigned power $a$ vanishes. Their product divides $a/U$ since $S$ is coprime to $U$. It also divides $e+hr$ by the identity
$$ e+hr=(dr+e)a-e(N-1).$$
Finally each power of $S$ divides at least one of $b,c$, so $S\mid\operatorname{lcm}(b,c)=Mxy$. Since $S\mid N-1$, this implies $S\mid gxy$. In sum,
$$ S\mid g(a/U),\qquad S\le g|e+hr|\ \ (e+hr\ne0),\qquad S\le gr^2/4. \tag{3}$$

If $e+hr=0$ then $r=e/H$, so (2) follows directly from (3). Otherwise, $N-1=KS\le gK(e+Hr)$ and $d(N-1)=M(e+dr)+h-d$. Consequently
$$ M\le gK\max(d,H)+\frac{E}{e+dr}.$$
Using $N-1\le gK a/U$ and $da=eM+h$ now gives (1): apply $K^2/U\le G^3 H^2/27$, $K/U\le G^{3/2}\sqrt{H/27}$, and $g\le E$. At $G=1$ these two inequalities prove the fifth-root assertion. The [Lean file](PartGcdGap.lean) verifies the resultant divisibility and the exact identities used in the finite search.

If $a=M$, then $A=M(1,x,y)$ with $\gcd(x,y)=1$. Every pair of actual parts has gcd $M$, which is coprime to $N-1$. Thus $S=1$. Both nonendpoint rows cannot have such a part at $G=1$ because their predecessor products multiply to $N-1>1$.

## Comparing the two nonendpoint rows

At a **single reduced positive rational slope** $e/d$, at most one of the two nonendpoint rows can have any part within $o(N^{1/5})$ of that slope times the gcd of its companions along a sequence of counterexamples. More precisely, choose any positive constant smaller than the asymptotic lower gap implicit in (1), (2) at that slope. For all sufficiently large $N$, at least one of the rows has all three errors $|da-eM|$ above that constant times $N^{1/5}$.

To see this, a small error at either row must put it on the exact line $da=eM$ or the shifted line $d(a-1)=eM$. On the exact line write $M=dL$, $a=eL$. The primitive row is $(e,dx,dy)$, with $\gcd(x,y)=1$. Its coordinate gcd before normalization is $L$ and $N=L(e+dr)$. Hence $N/U$ is a multiple of $e+dr\ge e+2d$. Its assigned predecessor product divides $de$: the primes of $S$ are coprime to $L$, the group with residue one at $a$ divides $d$, and the other two groups together divide $e$. On the shifted line $U\mid h=d$.

For slopes $e/d$ and $e'/d'$ at the two rows, both exact would give $N-1\le de\,d'e'$. Both shifted would give $N\le dd'$. Exact at the first and shifted at the second requires $d'\ge e+2d$ because $N/U$ is the second weight product. The reverse orientation requires $d\ge e'+2d'$. Neither mixed orientation is possible at a common slope. The same exclusion works for every pair of slopes whose larger reduced denominator is at most twice the smaller one.

In particular, take any subset of the reduced slopes $e/d$ with $1\le e,d\le8$ whose largest denominator is at most twice the smallest. In any counterexample **one nonendpoint row simultaneously has $|da-eM|>256$ at all three positions and all slopes in the subset**. The verified initial interval disposes of the two small exact-case bounds. The finite range calculation below supplies the other cases. The row can depend on the subset.

There is also a single-row consequence uniform over *all* positive integer numerators $t$, with no bound on $t$. Put $h=a-tM$ and $H=|h|>0$. At $G=1$ the case $h=1$ is impossible. The scale bound in the proof of (1) gives, whenever $t+hr\ne0$,
$$ M\le (H+1) H (H/3)^{3/2}+(H+1)/3. $$
Thus a nonzero remainder at any such multiple has size $\gg M^{2/7}$, independently of the multiple. The omitted case has $h=-H$, $t=Hr$ and the exact value $a=H(N-1)/(H+1)$.

## Small remainders at many slopes

Use (1), (2) at $G=1$ and the [verified initial interval](full_range_1b.md). Each row in this table excludes all $0<|da-eM|\le H$ with $da-eM\ne d$ for every $1\le d,e\le L$ simultaneously.

| $L$ | 1 | 8 | 16 | 32 | 64 | 100 | 128 | 256 | 512 | 1000 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| $H$ | 121 | 80 | 69 | 60 | 52 | 48 | 45 | 35 | 13 | 5 |

The script verifies each bound by exact integer squaring after substituting $E\le d+H$, $e\le L$. It checks (2) separately, using the actual $|h|\ge1$ in its denominator.

## Exhausting the remaining scales

At $G=1$ choose $U\mid H$, $U\ge3$, $g\mid E$. Since $U\mid M$ and $\gcd(U,g)=1$, write $M=gUm$. Set $\beta=(h-d)/g$. The normalized parts of the companion row sum to $U$. Every full power in $T$ divides two normalized parts. Thus $T$ divides the lcm of their three pairwise gcds and is coprime to $U$. Enumerating the positive partitions of $U$ gives a short finite list of possible $T$. By (3),
$$ k=\frac{g T a}{U(N-1)}\ \in\mathbb Z_{\ge1}.$$
The case $e+hr=0$ lies below $10^9$ throughout the two searches by (2), as the program also checks. In the other cases the bound on $M$ above with $K=T$ gives a finite upper bound on $k$ for $N>10^9$. Clearing denominators gives
$$ m\big(kU(e+dr)-egT\big)=T(h/U)-k\beta. \tag{4}$$
The right side is nonzero. If it vanished, $m>0$ would give $kU(e+dr)=egT$. Multiply by $h$ and substitute $Th=kU\beta$ and $g\beta=h-d$ to obtain $e+hr=0$. Thus one enumerates $m$ by positive divisors of the absolute right side of (4), then reconstructs $r,M,a,N$. Check integrality, positivity, the two unitary splits at $N,N-1$ and the binomial cut $(a,N-a)$. This cut must avoid every prime assigned from $U$ because the trinomial coefficient is a multiple of $\binom{N}{a}$.

| Range | Assignments after unitary tests | After binomial cut | Distinct weights | Largest weight |
|---|---:|---:|---:|---:|
| $d=e=1$, $|h|\le512$ | 89 | 5 | 5 | 191,271,031,381 |
| $d=e=1$, $|h|\le1024$ | 819 | 111 | 101 | 249,249,002,003 |
| Coprime $d,e\le8$, $|h|\le256$ | 248 | 5 | 5 | 10,233,442,571 |

The [full endpoint solver](check_binary_gap_finite.py) excludes all these weights, testing 12,377 size-feasible assignments in the second range and 273 in the third. The independent slope-one check enumerates $M$ first, tries divisors $S$ of $ga/U$, then obtains $T$ by solving $ST\equiv h-1\bmod M$. It uses a factorial valuation for the binomial cut. At the four shifts giving the five weights in the first range ($-467,-421,-337,454$), this looser test leaves respectively 17, 12, 5 and 2 weights. All are excluded by the full endpoint solver as well. Direct exact trinomial computations through weight 70 check (1) in 251,719 cases and (2) in 165 cases, using the minimum actual endpoint gcd over all companions of each designated row.

~~~sh
# Requires python-flint for the endpoint checks above 10^9.
python3 check_part_gcd_gap.py --limit 1024
python3 check_part_gcd_gap.py --limit 256 --denominators 8 --numerators 8
python3 check_part_gcd_gap_cross.py -467 -421 -337 454 1024
~~~
