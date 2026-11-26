# Central cuts in the two nonendpoint rows

A cut is the sum of a nonempty proper selection of the parts of a row. Its binomial coefficient divides the row coefficient, since the parts on each side can be merged. Reflect each of two cuts independently by replacing it with its complement if necessary, and order the resulting indices as $1\le a\le b\le N/2$. Write
$$M=N-1,\qquad s=b-a,\qquad d=N-a-b,\qquad x=\frac{(N-2a)(N-2b)}{N}=\frac{d^2-s^2}{N}.$$

**Theorem.** Fix $s\ge1$ and $C>0$. Among weights $N\le X$ with $d\le C N^{3/4}$, at most $O_{s,C}(\log(2X))$ admit cuts of gap $s$ for which
$$G=\gcd\left(NM,\binom Na,\binom Nb\right)=1.$$
Except at finitely many weights, these cuts must satisfy
$$4b(N-a)=NM,\qquad d^2=(2s+1)N+s^2. \tag{1}$$
In the explicit range $|(x-1)^2-4s^2|<M$, the same identities hold at **every** weight with $G=1$. Equal cuts always give $G>1$. The [endpoint theorem](all_weights_endpoint_reduction.md) applies these restrictions to any cuts of the two nonendpoint rows of a trinomial counterexample.

For $1\le s\le512$ and $3\le N<2^{4096}$, no pair of **positive trinomial** rows with cuts satisfying (1) has endpoint gcd one. The binomial version has exactly two exceptions in this box, and the refinements of both can be checked separately.

## A congruence for arbitrary cuts

Put $D=\gcd(N,\binom Na,\binom Nb)$ and $J=\gcd(M,\binom Na,\binom Nb)$. Then $H=Dx$ is an integer and
$$M/J\ \mid\ (H-D)^2-4D^2s^2. \tag{2}$$

Use the [retained prime powers](endpoint_pair_bound.md), which work also for two-part rows by the same factorial valuation proof. The prime powers of $N/D$ can be assigned to cuts so that each assigned power divides its cut. Thus $N\mid Dab$, proving integrality by expansion. At each assigned power $w$ of $M/J$, at least one cut is $0$ or $1$ modulo $w$. Since $N\equiv1\bmod w$, we have $H\equiv D(1-2a)(1-2b)\bmod w$. Substituting $a=0,1$ or $b=0,1$ in this product gives $D(1+2s)$ or $D(1-2s)$ modulo $w$. The powers $w$ are pairwise coprime, which proves (2).

Now take $G=1$, so $D=J=1$. If $s=0$ the weight powers all divide the same cut, which is strictly between zero and $N$. Hence $s\ge1$. For $|(x-1)^2-4s^2|<M$, (2) gives $x=1\pm2s$. Only the positive sign is possible since $x\ge0$. The identity for $d^2$ follows at once. Expanding $N(x-1-2s)$ gives $NM-4b(N-a)$, proving (1). For example, the explicit range includes $2s<\sqrt M$ and $|x-1|<\sqrt M$.

## Integral points outside the equality case

Here is a stronger form that also controls bounded common divisors. For fixed positive $C,B$ and fixed integer $s\ge1$, consider all pairs with $d\le C N^{3/4}$ and $G\le B$. With finitely many exceptions, they satisfy (1), or they have $s=1$ and
$$9(N-1)=4z^2(z^2-4) \quad\text{for an integer }z. \tag{3}$$
The second possibility requires $G>1$.

Indeed,
$$\lambda=\frac{(x-1)^2-4s^2}{M}$$
has denominator dividing $JD^2$ by (2), and $0\le x\le C^2\sqrt N$. Hence $\lambda$ takes only finitely many values. Negative values give bounded $M$, since $\lambda M\ge-4s^2$. Zero gives (1). For each fixed positive $\lambda$, eliminating $N$ yields
$$\lambda d^2=f(x),\qquad f(x)=x^3-2x^2+(1-4s^2+\lambda)x+\lambda s^2. \tag{4}$$
When $f$ has distinct roots this is a nonsingular affine curve of genus one. Since $Dx$ and $d$ are integers and $D\le B$, Siegel's integral-point theorem gives finitely many solutions. For instance, for each fixed $D$ choose a positive integer $L$ clearing the coefficients of $\lambda f(X/D)$. Set $X=Dx$ and $Y=L\lambda d$, also clearing the denominator of $\lambda$ in $L$. This gives $Y^2=L^2\lambda f(X/D)$ with a squarefree cubic over the rationals and integral $X,Y$. Each $x$ determines $N$ through $\lambda M=(x-1)^2-4s^2$. [Siegel's original paper and English translation](https://doi.org/10.1007/978-88-7642-520-2)

A repeated root $r$ of the rational cubic is rational, since a nonrational repeated root and its distinct conjugate would use at least four roots, counting multiplicities. Equations $f(r)=f'(r)=0$ give
$$\lambda=4s^2-3r^2+4r-1,\qquad \lambda s^2=2r^2(r-1).$$
Thus $r>1$ and $2r$ is integral, by the rational root theorem applied to
$$2r^3+(3s^2-2)r^2-4s^2r+s^2(1-4s^2)=0.$$
Viewed as a quadratic in $s^2$, the eliminated equation has discriminant
$$(3r^2-4r+1)^2+32r^2(r-1)=(r+1)^2(r-1)(9r-1).$$
For the integer $k=2r\ge3$, this requires $9k^2-20k+4$ to be an integer square. For $k>3$ it lies strictly between $(3k-4)^2$ and $(3k-3)^2$. At $k=3$ the positive root for $s^2$ is $1$, giving $\lambda=9/4$ and $f(x)=(x-3/2)^2(x+1)$. Here $d>0$, so $x+1=z^2$ for a rational $z$ by (4). Substitution gives (3). A reduced denominator of $z$ has fourth power dividing $4$, so $z$ is integral. If $G=1$, (2) makes $\lambda$ integral and excludes $9/4$. This proves the finiteness assertion outside (1). The argument supplies no numerical cutoff for the finite set.

## One candidate per residue in a binary block

In (1) with $G=1$, $N\mid ab$ implies $N\mid b(N-a)$, so $4\mid M$. Let $P$ be the unique power of two with $P/2\le N<P$. Since $2\mid NM$, one binomial is odd. Its two parts add without a binary carry, so one of those parts contains the leading bit $P/2$. The largest of the four parts is $N-a$, giving $N+d+s=2(N-a)\ge P$. For $t=2s+1$, the identities in (1) now give
$$(d+s)(d+s+1)\ge tP,\qquad d^2<tP+s^2,\qquad d^2\equiv s^2\pmod t.$$
Consequently
$$\lfloor\sqrt{tP}\rfloor-s\le d\le \lfloor\sqrt{tP+s^2-1}\rfloor.$$
The difference of the endpoints is at most $2s<t$. Since $t$ is odd and coprime to $s$, the last congruence has exactly $2^{\omega(t)}$ roots by factoring $(d-s)(d+s)$ prime power by prime power. Each root gives at most one $d$ in the block, hence at most one weight. This gives the stated logarithmic count after the finitely many other weights are added. In the explicit range, summing over $1\le s\le S$ gives $O(S\log(2S)\log(2X))$ possible pairs, since $2^{\omega(t)}$ is at most the divisor count and $\sum_{t\le T}\sum_{h\mid t}1\le T(1+\log T)$.

## Reproduction

Run `python3 check_central_cuts.py`. The first check uses exact binomial coefficients for all 338,349 oriented cut pairs through weight 200. The search then enumerates the residue representatives in every binary block for gaps 1 through 512 and weights below $2^{4096}$. Exactly 244 triples $(N,a,b)$ pass the positivity, parity and binary tests. Of these, 169 pass the full low-power tests at $N$ and $N-1$. The tests require the lcm of $\gcd(N,a)$ and $\gcd(N,b)$ to be $N$, and the lcm of $\gcd(M,a(N-a))$ and $\gcd(M,b(N-b))$ to be $M$, as follows from the assigned powers. For 167 triples the checker finds an explicit common prime and verifies it using factorial valuations.

The checker also reconstructs the 169 possibilities from the weight factors. Write $N=uw$ with $u\mid b$ and $w\mid a$ by assigning each full weight power to an avoiding cut. In (1) the positive integers $i=(N-2b)/u$ and $j=(N-2a)/w$ satisfy $ij=t$ and $wj-ui=2s$, hence $\gcd(i,j)=1$. For each ordered factorization $ij=t$ with $\gcd(i,j)=1$, $w$ lies in the single class $wj\equiv-1\bmod i$. Passing to the next pair $(u+j,w+i)$ increases $N$ by $2wj+1$. The binary test requires $P\le N+wj$, so only the largest $w$ in this progression with $uw<P$ can work. This independent enumeration, with parity checked by binary digit sums, returns the same set after the low-power tests.

The remaining binomial pairs have exact gcd one with $NM$.

| $N$ | $a$ | $b$ | distinct prime-divisibility masks for refinements of $a$, $b$ |
|---:|---:|---:|---:|
| 477 | 106 | 153 | 3, 6 |
| 32,477 | 15,886 | 15,893 | 3, 6 |

For each cut the program splits either side into every unordered pair of positive parts. It computes the resulting trinomial valuations at all primes dividing $NM$. At each of the two weights every mask from the first cut intersects every mask from the second. This verifies the trinomial conclusion for the whole stated range. Outside the central-cut hypotheses the general conjecture remains open.
