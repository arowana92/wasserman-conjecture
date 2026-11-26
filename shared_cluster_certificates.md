# The same four primes at consecutive weights

For every fixed positive integer $D$, there are $\gg_D X^{2/3}/(\log X)^3$ distinct products $M=PQR\le X$ of three primes with the following property. **For each $1\le d\le D$ there is a single list of four primes certifying every three positive trinomial rows at both $dM$ and $dM+1$. Four is the minimum size at each weight.** The list consists of $P,Q,R$ and a prime that may depend on $d$.

There is a version with a growing number of consecutive multipliers. Fix $\epsilon>0$. Uniformly for $1\le D\le X^{31/1269-\epsilon}$ there are $\gg_\epsilon X^{2/3}/(D\log^3 X)$ products $M\le X$ with shared four-prime certificates for every $1\le d\le D$. For these same products, both minimum sizes equal four at every $d\equiv1\bmod4$ in the range. The certificates prove the full conjecture directly at the indicated weights.

## Adding a fourth prime

Use the [three incompatible base primes](clustered_predecessors.md) at $n=dM$ and $n+1$. If a prime $\ell\le n$ satisfies $2\ell>n+1$ and $n+1-\ell\le P$, then $\{P,Q,R,\ell\}$ certifies **both** weights (with $P,Q,R$ prime). At either weight a row avoiding $\ell$ has two positive parts whose sum is at most $n+1-\ell$. At $n$ a row avoiding a core prime $W\in\{P,Q,R\}$ has all parts divisible by $W$. At $n+1$ at least one of the two small parts is divisible by $W$. Either possibility violates the sum bound. Thus all four primes are pairwise incompatible at each weight; three rows can avoid at most three of the four.

More generally a fourth prime base can come from a power $L$ with $n=mL+a$, $m\ge1$, $0\le a<L-1$, $\gcd(m,a)=1$, and $m(a+1)<P$. Require its base to be distinct from the three core primes. At $n+i$, $i=0,1$, a row avoiding the base of $L$ has parts $b_j L+c_j$ with $\sum b_j=m$, $\sum c_j=a+i$. If it also avoids $W\mid n$ from the core, its residues $\delta_j$ modulo $W$ sum to $i$. The integer $m(c_j-\delta_j)-a b_j$ is a multiple of $W$ of absolute value at most $m(a+1)$, so it vanishes. Hence $m(A_j-\delta_j)=n b_j$, where $A_j=b_jL+c_j$. Since $\gcd(n,m)=1$, $m\mid b_j$. At most one $b_j$ is positive. Of the other two positions, at least one has $\delta_j=0$ and therefore $A_j=0$. This proves the incompatibility. With $g=\gcd(m,a)>1$ the same calculation forces every $A_j-\delta_j$ to be a multiple of $n/g$.

## Counting close products

Take $Y=X^{1/3}/2$ and $H=\lfloor\sqrt{Y/D}/20\rfloor$. As in the [prime-bin count](clustered_predecessors.md), take $Y\le P<Q<R\le2Y$ with $R-P\le H$ and $R-Q>Q-P$. We may restrict all primes to the same fixed reduced residue class. The prime number theorem in arithmetic progressions still gives $\gg Y/\log Y$ primes with a constant depending on the modulus, and the upper bound on primes in short intervals is unchanged. There are $\gg YH^2/\log^3Y$ distinct products.

The number of ordered pairs of **integer** triples of this width whose products differ by at most $Y$ is
$$\ll_\eta Y^{1+\eta}H^2 \qquad(\eta>0). \tag{1}$$
To see this, write $S=P+Q+R$, $s=Q-P$, $t=R-Q$. Then
$$27M=S^3-3SB+C,\quad B=s^2+st+t^2,\quad C=-(2s+t)(s-t)(s+2t).$$
Also $0\le S-3M^{1/3}\le6H^2/Y$. Indeed, writing each factor as $S/3$ plus a deviation of size at most $H$ gives $0\le(S/3)^3-M\le3YH^2+H^3$. Divide by the sum of the three quadratic terms in a difference of cubes. Products within $Y$ have cube roots differing by at most $1/(3Y)$. Our choice of $H$ therefore forces their integer sums $S$ to agree for large $Y$.

Each $B\le3H^2$ has $Y^{o(1)}$ representations as $s^2+st+t^2$. For coprime $s,t$, $t$ is a unit modulo $B$ and $s/t$ is a root of $x^2+x+1$ there. Each such root determines the pair since the determinant of two representations has absolute value less than $B$. There are at most two roots at each prime power: roots are simple away from 3, and the only root modulo 3 does not lift modulo 9. For nonprimitive pairs sum over square divisors of $B$ and use the divisor bound.

If $B=B'$, the gap choices cost $O(H^2Y^{o(1)})$, with $O(Y)$ choices for $S$. Otherwise for each choice of gap pairs the product inequality is $|3S(B-B')-(C-C')|\le27Y$. This allows $O(1+Y/|B-B'|)$ integers $S$. Sum over $B,B'\le3H^2$, using $H^2\le Y$ and the harmonic sum. This proves (1).

## Finding the fourth prime at all multipliers

For each $d$, partition the range of $n=dM$ into bins of length $\lfloor Y/4\rfloor$. If $[n-Y/2,n]$ contains no prime, its bin meets a prime gap of length at least $Y/2$. A gap of length $g$ meets $O(1+g/Y)$ bins. The number of bins that can contain a failure is therefore
$$\ll_\eta Y^{-2}(dY^3)^{123/100+\eta},$$
by the mean-square prime gap theorem. A gap crossing the edge of the range is included by summing over all prime starts through $8dY^3$. [Stadlmann, *On the mean square gap between primes*](https://arxiv.org/abs/2212.10867v2)

The sum of squared product counts over all bins is at most (1), since two products $dM,dM'$ in a common bin have $|M-M'|\le Y$. Cauchy–Schwarz with the count of bad bins, followed by summing over $d\le D$, bounds the number of products failing at any multiplier by
$$\ll_\eta Y^{369/200+O(\eta)}D^{223/200+O(\eta)}.$$
Here $H^2\ll Y/D$. The bound is $o(Y^2/(D\log^3Y))$ for $D\le X^{31/1269-\epsilon}$ by taking $\eta$ sufficiently small; the necessary strict inequality in terms of $Y$ is $D^{423}<Y^{31}$ with a power saving. Each remaining prime within $Y/2$ below $n$ gives the desired single list at both weights.

## Minimum size

A prime $b$ divides every coefficient at weight $N$ exactly when the base-$b$ digit sum of $N$ is at most two. Otherwise its digit units can be distributed to three nonempty parts without carries. If every prime is individually avoidable, no list of at most three primes suffices, by choosing a row for each prime.

For growing $D$ choose all of $P,Q,R$ congruent to 3 modulo 4. At $d\equiv1\bmod4$, $n=dM$ is $3\bmod4$, composite, and has at least three distinct prime factors. Its odd-base digit sum is odd and cannot be one. Its binary digit sum cannot be one or two since $n-1>2$ has exactly one factor of 2. At $n+1$ we have $4\mid n+1$ and $n$ is composite. Thus $n+1$ is neither prime, twice a prime, nor one more than a prime. For any $N$ outside these three forms, a universal base $b$ satisfies $b\le\sqrt N$ (the assertion holds for large $N$). Up to a bound $T$ this leaves $O(\sqrt T\log^2T)$ possible weights by enumerating digit positions. For all $d\le D$ together, the resulting discard is $O(Y^{3/2}D^{3/2}\log^2(DY))$, which is smaller than the main term throughout the stated range.

For fixed $D$, impose a fixed reduced residue class on all three primes. For each $1\le d\le D$ choose two new odd primes $u_d,v_d$ congruent to 2 modulo 3 and coprime to $d$. Choose the common class $c$ with $dc^3\equiv1\bmod u_d$ and $dc^3\equiv-1\bmod v_d$ for every $d$. Cubing permutes the units at these primes, so CRT supplies the class. Now both $dM-1$ and $dM+1$ have bounded odd prime divisors. For large $Y$, neither $dM$ nor $dM+1$ can be prime, twice a prime, or one more than a prime. For $dM$, use its three distinct large prime factors and the divisor of $dM-1$; for $dM+1$, use its bounded odd divisor and the compositeness of $dM$. Discarding the $O_D(Y^{3/2}\log^2Y)$ remaining candidate weights having a universal base proves the fixed-$D$ statement.

## Finite results and reproduction

The first pair of consecutive weights with a shared list of minimum size four at both weights is **51 and 52**. The shared lists there are $\{2,3,5,17\}$ and $\{2,3,7,17\}$. The checker computes every coefficient through weight 52, tests for a universal prime, and enumerates all four-cliques among primes pairwise incompatible at both consecutive weights. When no prime is universal, four primes certify a weight exactly when they are pairwise incompatible: a row avoiding two would leave only two primes for the other two rows to avoid.

Take $(P,Q,R)=(1000000007,1000000009,1000000021)$. The three-prime theorem applies through $d=\lfloor P/(2\cdot14^2)\rfloor=2{,}551{,}020$. For every multiplier in this entire range the largest prime at most $dM$ gives a shared fourth prime. The largest distance to this prime is 1,115, at $d=1{,}958{,}064$.

| Minimum at $dM$ | Minimum at $dM+1$ | Multipliers |
|---:|---:|---:|
| 4 | 4 | 2,467,252 |
| 4 | 1 | 50,145 |
| 1 | 4 | 33,045 |
| 1 | 1 | 578 |

At $d=1$ the two weights start at $1{,}000{,}000{,}037{,}000{,}000{,}399{,}000{,}001{,}323$, and the fourth prime is 52 smaller. The four-prime certificate itself is valid also when one weight has a universal prime.

The enumerator uses FLINT's proved primality predicate. To test for a universal odd base of an even $N$, it checks $N-1$ and all possibilities $N=2b^j$ or $N=b^j+b^k$ with $0\le k<j$, using integer roots. For $j\ge2$ in the last case $b=\lfloor N^{1/j}\rfloor$ because $b^j+b^{j-1}<(b+1)^j$. At odd $N$ the only odd candidates are bases of pure prime powers. The binary check uses the digit sum directly.

The independent checker factors both endpoint coefficients at 502 sampled multipliers and tests all their prime bases. It also proves the sample fourth primes at $d=1$ and $d=2{,}551{,}020$ by recursive modular order certificates. Small exact trinomial products give 186,292 successor comparisons including 42 successful transfers outside the incompatible cases. Factorial floors independently check 1,558,582 normalized binomial rows across 178 increasing-gap triples of prime powers up to 1,200. The close-product identity and the count by common sum are checked on 5,137 close pairs of integer triples, including 144 with different quadratic norms. There are also 105 explicit jointly avoiding trinomial rows in the opposite gap ordering.

~~~sh
# Both commands use python-flint. The parallel job count can be changed.
python3 check_clustered_predecessors.py
python3 check_cluster_certificates.py --jobs 10
~~~
