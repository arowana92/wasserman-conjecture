# Counting exceptional weights with several prescribed endpoint primes

Let $E=N(N-1)$. Count a weight $N$ once if some three positive trinomial rows of sum $N$ have coefficients with gcd one. The following upper bounds hold through $X$. The density column describes all weights satisfying the indicated congruence.

| Condition | Natural density | Number of possible exceptions |
|---|---:|---:|
| $30\mid E$ | $4/15$ | $O(X^{0.7593})$ |
| $210\mid E$ | $8/105$ | $O(X^{0.7265})$ |
| $2310\mid E$ | $16/1155$ | $O(X^{0.7106})$ |

For every $\epsilon>0$ there is a fixed primorial $Q$ such that the number of possible exceptions with $Q\mid E$ is $O_\epsilon(X^{2/3+\epsilon})$. Explicit exponents are $0.70$, $0.68$ and $0.67$ for the primorials ending at $79$, $283$ and $2113$, respectively. Each is a positive-density set of weights. The implied constants may depend on the primes.

## A weighted digit count in any fixed collection of bases

For distinct primes $p$ in a nonempty finite set $I$, fix $0<w_p\le1$. Write
$$
 D_p(w)=\log_p\sum_{j=0}^{p-1}w^j,\qquad A_p(w)=(p-1)\log_p(1/w),\qquad
 \Delta_I(w)=\max\big(0,1+\textstyle\sum_{p\in I}(D_p(w_p)-1)\big).
$$
Let $M_p=p^{L_p}\asymp X$ and choose arbitrary integer shifts $e_p$. Uniformly in these shifts,
$$
 \sum_{z:\ 0\le z+e_p<M_p\ (p\in I)}
       \prod_{p\in I}w_p^{s_p(z+e_p)}\ll_\eta X^{\Delta_I(w)+\eta}
 \qquad(\eta>0).                                                   \tag{1}
$$
Here $s_p$ is the digit sum in base $p$. A single base gives (1) directly by multiplication over digit positions. For an empty set of bases use $\Delta_\varnothing=1$ for variables supported on an interval of length $O(X)$.

For $d=|I|\ge2$, take the product of the Bernoulli digit measures on $[0,1]$ with digit probabilities $w_p^j/\sum_{l=0}^{p-1}w_p^l$. The rectangles specified by the words for $z+e_p$ have disjoint interiors and null boundaries. Each has mass equal to its summand in (1) divided by $\prod_p M_p^{D_p(w_p)}$. The rectangles lie within $O(1/X)$ of a line with direction $(1/M_p)_p$. The normal $(d-1)$-plane contains no coordinate axis, and its directions range in a compact set with this property.

For $q>1$ the $L^q$ dimension of each factor is at least $D_p(w_p)$, since $\sum_j\Pr(j)^q\le(\max_j\Pr(j))^{q-1}$. The dimension formula for the projected product gives the minimum of $d-1$ and the sum of the factor dimensions. Proposition 3.2 gives a locally uniform finite-scale lower estimate for the projected spectrum, also uniform in the model parameter. A finite cover gives uniformity over the normals above. A ball of radius $h$ in the projection meets boundedly many mesh cells of side comparable to $h$. Bound their masses by the $q$-norm of all cell masses, then choose $q$ large enough to absorb the factor $(q-1)/q$ into an arbitrarily small loss. The mass of the tube is $O_\eta(X^{-\min(d-1,\sum D_p)+\eta})$. Multiplication by the normalizing factors proves (1). The sources are Lemma 5.1, Proposition 3.2, Theorem 5.4 and the integer-base application in Section 6.2 of [Corso and Shmerkin](https://arxiv.org/pdf/2409.04608). Details about the scale hypothesis appear below.

For any random integer $Z$ in the summation range, relative entropy against the normalized summands in (1) gives
$$
 \frac{\mathcal H(Z)}{\log X}\le \Delta_I(w)+
 \sum_{p\in I}\frac{\mathbb E[s_p(Z+e_p)]\log(1/w_p)}{\log X}+o(1).       \tag{2}
$$
Here and below an upper bound with $o(1)$ allows an arbitrarily small fixed exponent loss. All prime sets and parameters are fixed before $X$ tends to infinity.

## Applying the digit count to both rows

The [endpoint reduction](all_weights_endpoint_reduction.md) supplies two rows $A,B$ with $\gcd(E,K_A,K_B)=1$, where $K_A,K_B$ are their coefficients. For each full prime power dividing $N$ or $N-1$, choose a row avoiding its prime, resolving ties deterministically. In each row the product of its assignments is $P_i=U_i S_i$ from $N$ and $N-1$, respectively. Every part $z$ satisfies $U_i\mid z$ and $S_i\mid z(z-1)$. All parts exceed one. Indeed a part equal to one forces its coefficient to be divisible by $N$, and every positive trinomial coefficient has a nontrivial factor of $N$ by the decrement identity. We work on $[X,2X)$, label the rows so $P_A\le P_B$, and take a canonical witness for each weight.

Weights with $P_A\le X^k$ number $O(X^{k+o(1)})$ for $0<k<1$, by [reconstruction from the complementary row](two_predecessor_congruences.md). Among the other weights $P_A>X^k$ and $P_B\ge\sqrt E\gg X$. Stratify by the assignment of the prescribed primes to the two rows. Give the weights in each stratum the uniform law. Write $H=\mathcal H(N)/\log X$. Each part determines its row's label $(U_i,S_i)$ up to $X^{o(1)}$ divisor choices. That label determines $N$ modulo $P_i$. Conditional subadditivity for the three parts therefore gives
$$
 H\le \frac1{3\log X}\sum_{z\in i}\mathcal H(z)+
      \begin{cases}2(1-k)/3 & i=A,\\0&i=B\end{cases}
      +o(1).
$$
The [entropy argument](common_factor_count.md#information-that-every-individual-part-supplies) counts the parts by position.

Let $P$ now denote the fixed prime set dividing $E$, partitioned as $I$ on $A$ and $J$ on $B$. Put $L_p=\lceil\log_p(2X+1)\rceil$ and $r_p=\mathbb E[s_p(N)]/((p-1)L_p)$. The three parts of the row assigned $p$ have digit sums totaling $s_p(N)$. Their complements in the argument are the words $M_p-1-N$, with digit sums $(p-1)L_p-s_p(N)$. These words are all shifts of the same variable $-N$, exactly as required by (2).

Use tilts $x_p$ for the individual parts and $u_p$ for the complements. Applying (2) to the three positions of each row and to $-N$ gives three upper bounds for $H$, up to $o(1)$:
$$
 \Delta_I(x)+\tfrac13\sum_{p\in I}A_p(x_p)r_p+\tfrac23(1-k),\quad
 \Delta_J(x)+\tfrac13\sum_{p\in J}A_p(x_p)r_p,\quad
 \Delta_P(u)+\sum_{p\in P}A_p(u_p)(1-r_p).
$$
Give these bounds positive weights $a,b,c$ of sum one. Put $v_p=a$ on $I$ and $v_p=b$ on $J$. Maximizing over $0\le r_p\le1$ gives $H\le T+2a(1-k)/3+o(1)$, where
$$
 T=a\Delta_I(x)+b\Delta_J(x)+c\Delta_P(u)
       +\sum_{p\in P}\max\{v_p A_p(x_p)/3,cA_p(u_p)\}.
$$
At $k=(T+2a/3)/(1+2a/3)$ the two counting bounds coincide. For the three moduli the checker certifies across all $8,16,32$ assignments that these exponents are below $.75921,.72644,.71054$, respectively. It uses rational tilts, outward bounds for logarithms and the displayed coordinatewise maximum. The strict margins give the table after summation over dyadic intervals. Its density entries follow by choosing $0$ or $1$ independently at each prime.

## An explicit criterion for arbitrarily many bases

Fix $2/3<t<1$. For any fixed finite pool of primes $P_0$, restrict to weights for which
$$
 \sum_{p\in P_0,\ p\mid E}\frac{1}{\log p}>
 \frac{2t^2}{2t-1}+\frac{(3t-1)^2}{2(3t-2)}.                          \tag{3}
$$
There are $O_{P_0,t}(X^t)$ possible exceptions in this set.

Indeed take any nonempty subset of assigned primes $I$ on a single row, write $M=\sum_{p\in I}1/\log p$, and put $w_p=\exp(-\lambda/(p-1))$. For a uniform base-$p$ digit $D$, $Y=D/(p-1)$ has mean $1/2$ and lies in $[0,1]$. Convexity gives $\mathbb E e^{-\lambda Y}\le(1+e^{-\lambda})/2\le e^{-\lambda/2+\lambda^2/8}$, by $\cosh x\le e^{x^2/2}$. Hence $\sum(1-D_p(w_p))\ge M(\lambda/2-\lambda^2/8)$ and $\sum A_p(w_p)=M\lambda$.

For $M>2$ choose $\lambda=2(1-\sqrt{1-2/M})$. Then $\Delta_I(w)=0$. Use the same tilts on the three parts and the complementary words. Their four digit sums total $(p-1)L_p$ at each $p$. Weight the row bound by $3/4$ and the complement bound by $1/4$. This gives
$$
 H\le F(M)+\begin{cases}(1-k)/2&\text{on the smaller product},\\0&\text{on the larger product},\end{cases}
 \quad F(M)=\frac{M\lambda}{4}=\frac1{1+\sqrt{1-2/M}},
$$
up to $o(1)$. For $1/2<s<1$, $F(M)<s$ whenever $M>2s^2/(2s-1)$.

If the larger product receives reciprocal-log mass exceeding the first term in (3), its bound is strictly below $t$. Otherwise the smaller product receives mass exceeding the second term, giving $F(M)<(3t-1)/2$. In the latter case balance the small-product count with the one-row entropy bound at $k=(1+2F(M))/3<t$. In the former case the larger product bound holds without any removal of small products. Finitely many assignments and intersections $P_0\cap\{p:p\mid E\}$ suffice. This proves (3).

The sum of $1/\log p$ over primes diverges, for example by comparison with the divergent sum of reciprocal primes. Thus taking all primes through a sufficiently large fixed prime in (3) proves the primorial assertion. Its classes have density $2^{\omega(Q)}/Q$. The numerical examples in the introduction use certified lower bounds on this sum.

## Checking the logarithmic scale hypothesis in the projection formula

The preprint's display (P1) omits logarithms in the ratio of contractions. Its Lemma 5.1 uses the phase increments $\log p_d/\log p_i$, and Section 6.2 applies the projection formula to integer bases. Here is a direct verification of the projected exponential separation used in the proof of its Theorem 5.4 for the digit grids in (1). It applies to each coordinate subcollection as required by the induction there.

Order the $d\ge2$ bases as in Lemma 5.1. In the phase torus of the first $d-1$ coordinates take the closed subgroup generated by $(\log p_d/\log p_i)_{i<d}$. Every coordinate projection maps onto the circle, by irrationality. The identity component is a rational subtorus of finite index. Its rational tangent space thus contains a vector $v$ with all entries nonzero, by avoidance of finitely many proper hyperplanes. On bounded local lifts, differentiating along $v$ gives the distinct rates $v_i\log p_i$ for $i<d$, and rate zero for the final coordinate. Distinctness follows from the irrationality of every pairwise log ratio.

Fix a line projection whose direction has all entries nonzero. For two different level-$n$ atoms, on a local lift their projected difference is of the form $f(t)=\sum_{i=1}^d c_i e^{\beta_i t}$ along parallel segments in direction $v$. The rates are those above. At level $n$ each coordinate has $O(n)$ integer-base digits. Distinct digit words on the same grid differ by at least $\exp(-C n)$, and the phase scalings on the lifts are bounded and bounded away from zero. Consequently $\exp(-Cn)\le\max_i|c_i|\ll1$ whenever the atoms differ. One can include all grid lengths occurring at that level, and all pairs from each grid, in a list of at most $\exp(Cn)$ candidate differences. The coordinates of a candidate are fixed grid values before multiplication by phase scalings. The estimates hold on any part of the chart where that grid is used.

After division by $\max_i|c_i|$, the derivatives of orders $0,\ldots,d-1$ cannot all have absolute value below a fixed positive constant at any point of a bounded segment. This follows by inverting the Vandermonde matrix of the distinct $\beta_i$. Derivatives through order $d$ are bounded above uniformly. Divide the segment into boundedly many sufficiently short intervals. On each, some derivative of order $j\le d-1$ stays bounded away from zero. The set where the normalized function has absolute value below $\delta$ has length $O(\delta^{1/(d-1)})$. For $j=0$ it is empty for small $\delta$. For $j\ge1$, a set of length $l$ contains $j+1$ points separated by $\gg_j l$, and the divided-difference mean-value formula at these points gives $l^j\ll_j\delta$.

Haar measure is Lebesgue measure in finitely many such lifted coordinate charts, up to constant factors. Cutting charts at coordinate discontinuities gives boundedly many intervals along each line. Integrate the sublevel estimate in the transverse parameters. The union of all candidate events $|f|<\exp(-Rn)$ has summable measure in $n$ when $R$ is large enough. Borel–Cantelli proves exponential separation for almost every phase. It also gives injectivity on the one-step support almost everywhere, or this follows directly by applying the same sublevel estimate at a single level with $\delta\downarrow0$. A one-coordinate subcollection already has the separated digit grid. This verifies the input of Lemma 5.2 in the integer-base case, with arbitrary positive digit probabilities. The induction for Theorem 5.4 and the uniform estimate of Proposition 3.2 then apply. [Corso and Shmerkin, Sections 3.1, 5 and 6.2](https://arxiv.org/pdf/2409.04608)

Reproduce the rational certificates and the independent checks of the four digit sums against exact binomial products with the standard library.

~~~~sh
python3 check_joint_prime_counts.py
~~~~
