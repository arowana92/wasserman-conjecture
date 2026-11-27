# Common divisors at rationally spaced exponents and along power towers

Write $\alpha=\log_2\lambda=0.551463\ldots$, where $\lambda^3=\lambda^2+1$ and $\lambda>1$. Let $I$ be any interval of $K\ge2$ consecutive integer indices. Take a nonconstant rational function $R\in\mathbb Q(X)$ and consider only $j\in I$ at which $R(j)$ is a positive integer. Let $s$ be the number of distinct finite zeros and poles of $R'$ over an algebraic closure, counting a point once regardless of multiplicity. In particular $s\le1$ for every nonconstant shifted power $u+v(w+xj)^m$, uniformly in the positive degree $m$.

**Odd prime powers.** Start with $R(j)$ and apply any fixed finite list of operations $n\mapsto B^n$ with odd integer bases $B\ge3$, obtaining $Q_j$. The list may be empty. The number of indices where $Q_j$ is an odd prime power and some weight $1+dQ_j$ with $1\le d\le Y$ fails Wasserman's conjecture is
$$O\big((1+s)Y^4 K^\alpha\big). \tag{1}$$
The constant is absolute. It is uniform in the rational coefficients, the interval, and the bases and height of the tower. The prime underlying $Q_j$ need not be fixed when the tower is empty.

For example, for any odd prime $q$ and any positive $a,b,m$, take
$$Q_j=q^{5^{7^{(a+bj)^m}}},\qquad 1\le j\le K.$$
At all but $O(K^{\alpha+4/9})$ indices, **all** cofactors through $K^{1/9}$ satisfy the full three-row conjecture. Here $\alpha+4/9<1$. The bound still holds if the inner 5 and 7 are replaced by an arbitrarily long fixed list of odd bases. With no tower, (1) also applies directly to prime-power values of a rational function, including in intervals arbitrarily far from the origin.

**Binary prime powers.** Set $Q_j=2^{R(j)}$. For $Y\ge3$, the number of indices with a failure at any $1+dQ_j$, $1\le d\le Y$, is at most
$$O\!\left((1+s)Y^9 K^{1-\kappa/\log Y}\right),\qquad \kappa=\log(9/7). \tag{2}$$
For example, all cofactors through $\exp(\sqrt{\log K}/8)$ satisfy the conjecture simultaneously outside $O((1+s)K\exp(-0.8\sqrt{\log K}))$ indices. This follows from $8\log(9/7)-9/8>0.8$. Logarithms without a subscript are natural.

For the 32,768 odd cofactors through 65535, the sharper count is $O((1+s)K^{0.36})$. At $d=12031$ alone it is $O((1+s)K^{0.14})$. Both hold even after inserting a power tower between $R(j)$ and the outer base 2, with no dependence on its height. For the joint statement all inserted bases must be coprime to $3\cdot5\cdot7\cdot11$. At 12031 it suffices that none be divisible by 3. The earlier interval counts covered consecutive exponents. In particular these new bounds apply to every arithmetic progression of exponents, with a constant independent of its step. [Binary cofactor and carry-state classification](binary_cofactor_rays.md)

## How a rational function acts on residue classes

Normalize $v=v_p$ by $v(p)=1$ and put $v(0)=+\infty$. Use its extension to a finite splitting field of the derivative numerator and denominator. On a ball $x+p^n\mathbb Z_p$ suppose every finite zero or pole $z$ of $R'$ satisfies $v(x-z)<n-1$. Then throughout the ball $v(R')$ is constant, and for distinct rational points $x_1,x_2$ in it
$$v(R(x_2)-R(x_1))=v(R'(x_1))+v(x_2-x_1). \tag{3}$$
All finite poles of $R$ are poles of $R'$. The valuation extension used here is the standard one for finite extensions of $\mathbb Q_p$. [Brian Conrad, notes on completion of an algebraic closure](https://math.stanford.edu/~conrad/248APage/handouts/algclosurecomp.pdf)

For a proof put $M=\max_z v(x_1-z)$. If there are no such points then $R'$ is a nonzero constant and (3) is immediate. Otherwise the derivative ratio has the expansion
$$R'(x_1+T)/R'(x_1)=\prod_z(1+T/(x_1-z))^{m_z}.$$
The multiplicities $m_z$ are signed integers, so all binomial coefficients in these factors are integers. The coefficient of $T^{k-1}$ has valuation at least $-(k-1)M$. For $T=x_2-x_1$, the Taylor term of degree $k\ge2$ for $R$, divided by $R'(x_1)T$, therefore has valuation at least
$$(k-1)(v(T)-M)-v_p(k)>0.$$
Indeed $v_p(k)\le k-1$ and $v(T)-M>1$. These bounds tend to infinity, so the expansion converges. It is the Taylor expansion of the rational function on this ball, since every pole is farther away. Its linear term has uniquely least valuation. This proves (3) without a factor depending on the multiplicities.

At depth $n\ge1$, a potentially bad residue ball has a point with $v(x-z)\ge n-1$ for some $z$. For each $z$, all such integer points agree modulo $p^{n-1}$, so there are at most $ps$ potentially bad balls of depth $n$. Start at depth zero and split only potentially bad balls. At most $p^2(s+1)$ new good balls occur at each depth. Stop at $\lceil\log_p K\rceil$, when every remaining ball contains at most one index of $I$. The estimate is independent of where the zeros and poles lie and whether they are rational.

## Common low digits in a tower

Fix a prime $p$ and assume all tower bases are units at $p$. Within each good ball divide the eligible indices into a bounded number of classes so that for distinct $x,y$ in a class the final outputs satisfy
$$v_p(Q_x-Q_y)=c+v_p(x-y) \tag{4}$$
with one constant $c$ on the class. The following construction makes the number of classes independent of the height.

Use $L=\operatorname{lcm}(1,\ldots,p-1)$ for odd $p$, and $L=2$ for $p=2$. Take $R(j)$ modulo $L$, separating individually its positive values below $1+\lfloor\log_2 L\rfloor$. By (3) each individual value occurs at most once in the ball. For $n,m$ at or above the cutoff with $n\equiv m\pmod L$, every $B\ge2$ satisfies $B^n\equiv B^m\pmod L$. At a prime power $r^e\parallel L$, both residues vanish for $r\mid B$ since $n,m\ge e$. Otherwise $\varphi(r^e)\mid L$: its relatively prime factors $r^{e-1}$ and $r-1$ both divide $L$ (also directly for $L=2$). The exponents stay above the cutoff throughout the tower.

For odd $p$, differences in each class are multiples of $p-1$, giving
$$v_p(B^n-B^m)=v_p(B^{p-1}-1)+v_p(n-m).$$
At $p=2$, differences are even and the added constant is $v_2(B^2-1)-1$. These identities follow by expanding powers of a principal unit. Iterating and applying (3) proves (4). With only the outermost base, one can simply color by the exponent modulo $p-1$ (modulo 2 at 2). With no tower, (3) suffices.

Suppose the scalar tested for carry avoidance is $T_j=(A Q_j+B)/U$ for fixed integers $A,U>0$ and $B$, considering only nonnegative integral values. On the good ball write $j=j_0+p^n i$. The index $i$ lies in an interval of length $k\le\lceil K/p^n\rceil$. Equation (4), with the affine rescaling included, gives
$$v_p(T_i-T_{i'})=t+v_p(i-i')$$
for some integer $t$ on a class with at least two indices. For $t\ge0$ the scalars share their first $t$ digits. Their next $\ell=\lceil\log_p k\rceil$ digits distinguish them. For $t<0$ the first $\ell$ digits already distinguish them. Thus any uniform upper bound $C p^{\delta\ell}$ on the number of readable blocks after a common prefix bounds the number of eligible indices here by $Cp^\delta k^\delta$. A singleton is covered as well when $C\ge1$.

Sum over the good balls, classes and remaining singleton balls. For $0<\delta<1$ the geometric sum over depths costs $O(K^\delta/(1-p^{-\delta}))$. Constants for a fixed $p$ and $\delta$ are independent of coefficients, starting position and tower height. For a single outer base, the explicit cost is $O((1+s)p^4 C K^\delta/(1-p^{-\delta}))$, including the at most $p$ phase classes. This is where a very long common digit prefix, arising from a large step or tower height, drops out of the count.

## Applying the carry bounds

A counterexample contains the endpoint row. For a prime power $Q$ dividing $N-1=dQ$, a row avoiding its base has the form $(xQ,yQ,zQ+1)$ with $x+y+z=d$. Assign full powers of $N$ to the two nonendpoint rows. The part assigned to this row has product $U\ge3$, divides $x,y$, and also divides $zQ+1$. Thus $z>0$. The other row is $V(f,g,h)$ with $UV=1+dQ$, $f+g+h=U\le d$. One of the two rows avoids 2. These necessary conditions work even when the base of $Q$ divides $d$. [Endpoint reduction](all_weights_endpoint_reduction.md) · [Affine carry-prefix bounds](binary_affine_counts.md)

For odd $Q$ use base 2 and the row avoiding it. With slopes $s_1,s_2,s_3$ and total $s_0$, set $D=\max_i v_2(s_i)-\min_i v_2(s_i)$. The binary block bound is $C2^{\alpha\ell}$ with $C=2^{\alpha D}$. It applies after deleting any common prefix because the new offsets are the corresponding nonnegative upper carries. Over all positive triples of sum at most $Y$, the sum of $C$ is $O(Y^3)$. Indeed $C\le\sum_{i=0}^3 2^{\alpha v_2(s_i)}$, each value at one position has multiplicity $O(Y^2)$, and $\sum_{m\le Y}2^{\alpha v_2(m)}=O(Y)$. The first row uses scalar $Q$, and its slopes determine $d$. The other uses scalar $(dQ+1)/U$ and costs at most a further factor $Y$ for the choice of $d$. This proves (1). [Sum over slopes](tower_heights.md#counting-rows-after-removing-their-low-digits)

For $Q=2^{R(j)}$ take an odd prime $p\mid U$ (the weight is odd). The first row must avoid it. The odd-prime block bound is $C((2p+1)/3)^\ell$, and the weighted-prefix proof gives $C\le e^D$ for the analogous spread of valuations at $p$. Indeed the weights cost at most $\exp(D(p-1)/p)$ for $p\ge5$ and $(5/4)^{2D}$ at 3. Here $D\le\log_p Y$ implies $C\le Y$. The exponent $\delta=\log_p((2p+1)/3)$ satisfies $1-\delta\ge\kappa/\log Y$. There are $O(Y^3)$ slope triples and at most $Y$ prime choices for each. Using $p^4\le Y^4$ and $1-p^{-\delta}\ge4/7$ proves (2). For the eleven specific rows the finite transition estimates bound blocks after any reachable live state following the first allowed digit. An empty common prefix costs at most one additional factor $p$. Applying (4) with their tabulated exponents proves the sharper finite-cofactor claims.

Run `python3 check_rational_exponent_counts.py`. The program checks (3) on 6,843 exact rational comparisons, including derivatives with nonrational zeros, and checks 305 valuation identities for two-layer towers. The binary exponent in the general block estimate is sharp for a prefix test alone. Put $Q=9^j$ and use the three parts $(Q,2Q,4Q)$. Across $2^k$ consecutive indices the residues modulo $2^{k+3}$ run through every value equal to 1 modulo 8. A passing prefix has ones separated by at least three positions. The free $k$ bits give $F_k=F_{k-1}+F_{k-3}$ with $F_0=1,F_1=2,F_2=3$, hence $F_k\asymp 2^{\alpha k}$. The program checks the first fifteen lengths directly. This lower bound concerns the binary test on those parts.
