# Exceptional orders with two unrestricted prime-power factors

Let $t=p^a q^b$ for arbitrary primes $p,q$ and positive integers $a,b$. Call a binary order $j$ exceptional if some $2^j\le t<2^{j+1}$ gives a failure of Wasserman's full trinomial conjecture at $5t+1$, at $7t+1$, or at $9t+1$ with $t$ even.

**Every interval of $H$ consecutive nonnegative binary orders contains at most $O(H^{1-1/10000000})$ exceptional orders, uniformly in the starting order.** The same statement holds if we also allow any failing weight $1+dt$ with $d\le1000$ and $t$ a prime power, by the [two-base theorem](logarithmic_intervals.md).

The following bounds for the two-factor family are effective at the stated lengths. In particular every window of $10^{1000}$ consecutive binary orders contains at least $10^{28}$ consecutive orders throughout which every weight of these forms satisfies the conjecture. The same run can be chosen to work for all the bounded cofactors at prime powers as well.

| Window length $H$ | Upper bound for exceptional orders |
|---:|---:|
| $10^{100}$ | $2.32\cdot10^{99}$ |
| $10^{300}$ | $1.60\cdot10^{293}$ |
| $10^{1000}$ | $5.02\cdot10^{971}$ |

## Leading digits in three bases

The [factor-five and factor-seven reductions](five_seven_pairs.md) and the [factor-nine reduction](factor_nine_pairs.md) leave the cuts below. Each requires $kt$ and $(d-k)t+1$ to add without carries in the three specified bases. The numbers $s_i$ are upper bounds on the dimensions of the corresponding sets of allowed digit blocks, checked by two path-count algorithms at block length 256.

| $d$ | $k$ | Bases | $(s_1,s_2,s_3)$ |
|---:|---:|:---|:---|
| 5 | 3 | 2,3,5 | $(.61,.69,.51)$ |
| 7 | 3 | 2,3,7 | $(.57,.66,.62)$ |
| 7 | 5 | 2,5,7 | $(.65,.74,.59)$ |
| 9 | 5 | 2,3,5 | $(.59,.56,.71)$ |
| 9 | 7 | 2,3,7 | $(.58,.50,.77)$ |

The remaining multiplier $k=6$ at $d=7$ has weakly monotone base-seven digits for $t$, and the binary cut dimension is at most $569/1024$. The [two-base proof](logarithmic_intervals.md) applies at either parity whenever this multiplier is forced. In particular its saving $(1-.56)/332000$ exceeds the one needed here, since the polynomial count of base-seven prefixes is absorbed in the margin between $.56$ and $569/1024$.

Write the three bases in a row as $2,p,r$, and put $B_1=2^{256}, B_2=p^{256}, B_3=r^{256}$. In coordinate $i$ take the compact digit set $K_i\subset[0,1]$ formed by allowing all the $B_i$-digits readable in that carry graph starting and ending at any reachable states. The zero block is available. These sets are invariant under multiplication by $B_i$ modulo one and have upper box dimensions at most $s_i$ by the elementary count of cylinders. Their bases are pairwise multiplicatively independent.

For a surviving $t$ put $n_i=1+\lfloor\log_{B_i}t\rfloor$ and $x_i=t/B_i^{n_i}$. Each $x_i\in K_i\cap[B_i^{-1},1]$: pad the highest block with zeros, read all blocks along the actual accepting path, then append fractional zero digits. On the compact positive rectangle $D=\prod_i[B_i^{-1},1]$ define
$$
 \Phi(x)=\left(\frac{\log x_2-\log x_1}{\log B_2},
                \frac{\log x_3-\log x_1}{\log B_3}\right).
$$
The phase $n_1\theta$ modulo $\mathbb Z^2$, where $\theta=(\log_p2,\log_r2)$, belongs to the image of $F=\Phi(D\cap\prod K_i)$ modulo $\mathbb Z^2$. Indeed its coordinates before reducing are $n_2+\Phi_1(x),n_3+\Phi_2(x)$. Each possible $n_1$ covers at most 256 ordinary binary orders, and a block of $H$ binary orders meets at most $H+1$ consecutive block orders.

## A uniform bound for straight lines in the phase plane

Let $\tau=\max(s_1+s_2+s_3-1,s_1,s_2,s_3)$. For every $\sigma>\tau$, there is a constant $C_\sigma$ such that, uniformly over every affine line $\ell$ and $0<\delta<1$, the set $F\cap\ell$ can be covered by at most $C_\sigma\delta^{-\sigma}$ squares of side $\delta$.

Here is the proof, including directions with a vanishing component. Write a line as $v\cdot y=c$ with unit normal $v$. The hypersurface in the digit coordinates has equation $g_v(x)=v\cdot\Phi(x)=c$. Set $h=v_1/\log B_2$, $k=v_2/\log B_3$. Its gradient is
$$
 (-(h+k)/x_1,\ h/x_2,\ k/x_3).
$$
At least two components are nonzero. By compactness the second largest absolute component is bounded below uniformly on $D$ and the unit circle. All second derivatives are uniformly bounded on a positive neighborhood of $D$.

For an affine plane whose normal has three nonzero components, Theorem 1.15 of Corso and Shmerkin with codimension one bounds the upper box dimension of its intersection with $\prod K_i$ by $\max(s_1+s_2+s_3-1,0)$. If the normal has exactly two nonzero components, apply their two-coordinate line theorem on those coordinates and take the product with the free coordinate. This gives $s_k+\max(s_i+s_j-1,0)\le\tau$. Both conclusions also hold for translated diagonal scalings of the product with positive bounded scale factors: pull the plane back under the diagonal map. [Corso and Shmerkin, *Dynamical self-similarity, $L^q$-dimensions and Furstenberg slicing in $\mathbb R^d$*, Theorem 1.15](https://arxiv.org/pdf/2409.04608)

We give the compactness argument needed to use these affine statements uniformly for the curved surfaces. In a cube $Q$ of side $R$, cover each $K_i$ by its digit cylinders of length between $R$ and $B_i R$. The projection of $Q$ meets at most three such cylinders per coordinate. Normalize $Q$ to the unit cube. Each product cylinder then lies in a translated diagonal copy of $\prod K_i$ with factors in $[1,\max B_i]$ and bounded translations. There are at most 27 copies.

Take any sequence of nonempty normalized pieces of $D\cap\prod K_i\cap\{g_v=c\}$ with $R\to0$. Pass to subsequences on which the pieces converge in Hausdorff distance, the copies converge in their parameters, and the tangent planes at a chosen point of each piece converge. The Taylor error in these coordinates is $O(R)$ after dividing by $R$, so the limit lies on a single plane with at least two nonzero normal components. It lies in the finite union of the limiting product copies. Its upper box dimension is therefore at most $\tau$.

Fix $\sigma>\tau$. There are $R_0>0$ and an integer $m_0$ such that every nonempty piece in a dyadic cube of side $R\le R_0$ meets at most $2^{m\sigma}$ dyadic children of relative side $2^{-m}$ for **some** $1\le m\le m_0$. Otherwise one could take $R\to0$ with this assertion failing for every $m$ up to a growing bound. The preceding compact limit has upper box dimension at most $\tau$. For some fixed sufficiently large $m$ it is covered by $O(2^{m(\tau+\sigma)/2})$ cubes of side $2^{-m}$. Enlarging each by a fixed number of neighboring cubes covers all sufficiently close pieces and still uses at most $2^{m\sigma}$ children, a contradiction.

Now subdivide each piece at its own such scale. Each subdivision increases neither the sum of the $\sigma$-powers of the side lengths nor the largest side. Stop individual branches when the side first drops to at most $\delta$. Every terminal side is at least $2^{-m_0}\delta$ for small $\delta$. The sum bound gives $O(\delta^{-\sigma})$ terminal cubes from a fixed initial cover at scale $R_0$. Finally $\Phi$ is Lipschitz on $D$, so the required squares cover its projected slice as well. All constants are independent of the line. This argument handles limits of normals in which a component tends to zero.

## Counting visits without joint logarithmic independence

Every row of the table has $\tau\le.98$. Fix $\sigma=.985$. For any integer $Q\ge2$, pigeonhole $Q^2+1$ successive multiples of $\theta$ into $Q^2$ squares in the unit torus. There is a step $1\le m\le Q^2$ for which
$$
 0<\rho=\max(\|m\theta_1\|,\|m\theta_2\|)\le1/Q .
$$
On the other hand $\rho\gg m^{-60000}$. The first coordinate has $p=3$ or 5. Apply Gouillon's two-logarithm bound to $m\log2-a\log p$, where $a$ is the nearest integer to $m\log_p2$. The coefficient of $\log m$ in the negative log bound is at most $7200(3.409+1.705)\cdot1\cdot\log5<60000$; the other terms are bounded independently of $m$. Small $m$ with $a=0$ are absorbed by the constant. [Details of the rational-height substitution](logarithmic_intervals.md#two-leading-prefixes-and-a-rotation) · [Gouillon, Corollary 2.3](https://www.numdam.org/article/JTNB_2006__18_1_125_0.pdf)

Split $H$ consecutive indices into progressions modulo $m$. Along a progression the orbit advances by a fixed vector of sup norm $\rho$ on the torus. Split it into segments with at most $\lfloor1/\rho\rfloor$ points each. Each segment lifts to a straight segment of bounded length, with successive points spaced by $\rho$ in sup norm. Only boundedly many integer translates of the bounded set $F$ can be relevant to this segment. The line estimate at scale $\rho/2$ therefore allows $O(\rho^{-\sigma})$ visits per segment. There are $O(H\rho+m)$ segments, independently of the initial phase. The count is
$$
 O(H\rho^{1-\sigma}+m\rho^{-\sigma})
  = O(HQ^{-(1-\sigma)}+Q^{120002}).
$$
Take $Q=\lfloor H^{1/120004}\rfloor$ for sufficiently large $H$. Both terms are $O(H^{1-1/10000000})$. This only uses a lower bound for **one** logarithm ratio and works at any starting phase. Summing over the finitely many necessary cases proves the uniform theorem.

## Computing much stronger bounds at specified window lengths

For the table work with ordinary base-$2,p,r$ digits and the index $i=\lfloor\log_2 t\rfloor$. First certify that any $H$ consecutive orbit points are pairwise separated in torus sup norm by more than $1/D$ for an explicit $D$ near $\sqrt H$. This is a finite lattice calculation in dimension three. Set $T=H^4$, choose integers $X,Y$ within 2 of $T\theta_1,T\theta_2$ using rigorous log intervals and put $w=\lceil T/(HD)\rceil$. Generate the lattice with rows $(w,X,Y),(0,T,0),(0,0,T)$. A displacement $0<n<H$ violating separation would have a nonzero lattice vector of squared norm at most
$$
 (Hw)^2+2(\lceil T/D\rceil+2H)^2 .
$$
LLL supplies a candidate reduced basis. The checker verifies the integer change of basis and its determinant, then verifies that **each** squared Gram–Schmidt length exceeds this threshold. The last nonzero coefficient of any lattice vector in this basis proves the asserted lower bound. Starting with $D=\sqrt H$, double $D$ until the certificate succeeds.

A triple of leading words of lengths $n_1,n_2,n_3$ covers a rectangle in phase space. Its sides have lengths at most $(P_1^{-1}+P_i^{-1})/\log b_i$ for $i=2,3$, where $P_i=b_i^{n_i-1}$ and $(b_1,b_2,b_3)=(2,p,r)$. This follows from the [two-prefix inequality](logarithmic_intervals.md#two-leading-prefixes-and-a-rotation). A side of length $l$ is covered by $1+\lfloor Dl\rfloor$ subarcs of length strictly less than $1/D$. Each product subrectangle holds at most one orbit point. Multiply their number by the number of prefix triples and add the initial binary orders where the words might be too short. Optimize over nearby word lengths. The sum over all six cuts gives the table. For $H=10^{1000}$ the pigeonhole bound on the longest gap is $(H-M)/(M+1)>10^{28}$, where $M$ is the combined count. The two-base checker adds fewer than $1.45\cdot10^{780}$ occupied orders for all the bounded cofactors at prime powers, preserving the same guaranteed run length. The finite checker also verifies six small rotation-separation bounds independently by comparing ordinary integer powers.

Reproduce the dimension bounds, then the finite orbit and prefix counts, with Python and python-flint.

~~~~sh
python3 check_paired_power_orders.py
python3 check_three_base_orders.py
~~~~
