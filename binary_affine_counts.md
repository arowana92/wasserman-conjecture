# How often can affine parts add without a carry?

Let $k\ge2$ and let $a_i t+b_i$ ($1\le i\le k$) have positive integral slopes and nonnegative integral offsets. Write $a_0=\sum_i a_i$, $b_0=\sum_i b_i$. Let $D$ be the difference between the largest and smallest 2-adic valuations of $a_0,a_1,\ldots,a_k$. Let $\lambda_k>1$ be the real root of $\lambda_k=1+\lambda_k^{-(k-1)}$.

**At most $\lambda_k^{L+D}$ residues of $t$ modulo $2^L$ extend to any nonnegative input at which the $k$ parts add without a binary carry.** This is sharp in its exponential rate. For triples the count exponent is $\alpha=\log_2\lambda_3=0.551463\ldots$. For two summands it is $\log_2((1+\sqrt5)/2)=0.694241\ldots$.

## Counting low digit strings

Use the binary digits of the $k$ parts and the *complement* of each digit of their total. A complete addition is carry-free exactly when these $k+1$ digits sum to one in every column. Put $v=\min_i v_2(a_i)$, including $i=0$. Divide all $k+1$ expressions by $2^v$, rounding the offsets down separately. The omitted low columns are independent of $t$ and the remaining columns still obey the same digit test. The new slopes have delays $s_i=v_2(a_i)-v$ between 0 and $D$. Reading input bit $j$ determines the digit of expression $i$ at position $j+s_i$. Given the previous input, the two choices of this bit give opposite values of the output digit, also when it is complemented. Lower digits at positions below $s_i$ are fixed initially.

The number $m$ of delays equal to zero is even and at least two, by the sum identity for the slopes reduced modulo two. If $m\ge4$, at most one next input bit is possible at every node, since the two choices together put $m$ ones in the current column. Suppose $m=2$. Keep a record of future columns already containing a determined one. Reject any collision of two ones in such a column. Immediately before input bit $j$, the recorded columns lie between $j$ and $j+D-1$. If there are $b$ of them, give the prefix weight $\lambda_k^{-b}$.

If the current column is recorded, both undelayed outputs must be zero, which permits at most one choice. Removing the current column multiplies the weight by $\lambda_k$ and all new ones can only decrease it. If the current column is free, the two choices create $z$ and $k-1-z$ future ones before collision checks. Their combined relative weight is at most
$$\lambda_k^{-z}+\lambda_k^{-(k-1-z)}\le1+\lambda_k^{-(k-1)}=\lambda_k,$$
by convexity in $z$. The initial weight is at most one and each surviving prefix has weight at least $\lambda_k^{-D}$. This proves the bound, even when the terminal digits have yet to be tested. For slopes $1,2,\ldots,2^{k-1}$ (successive powers of two) and zero offsets, the admissible inputs have ones spaced at least $k$ places apart. Their exact counts satisfy $F_L=F_{L-1}+F_{L-k}$ for $L\ge k$, giving the claimed sharpness.

## An odd prime base

For three affine parts as above, any odd prime $r$, and $D=\max_i v_r(a_i)-\min_i v_r(a_i)$ including the total, the number of extendible residues modulo $r^L$ is
$$O_{r,D}\big(((2r+1)/3)^L\big).$$
The constant is independent of the offsets. Use the digits of the three parts together with $r-1$ minus the total's digit. Their sum must equal $r-1$ at every column. After removing the common delay there are $m\ge2$ undelayed expressions. At each of these, the next input digit permutes all $r$ output values. If the current column already has assigned digit sum $S$, and $n$ choices pass, then $m n(n-1)/2\le n(r-1-S)$. This follows by summing the current digits over the $n$ choices. For $m\ge3$ the branching is at most $(2r+1)/3$.

For $m=2$, give each prefix weight $z^{-b}$ where $z>1$ and $b$ is the sum of pending digit values in future columns, including the current column. Prune sums above $r-1$, so $b\le D(r-1)$. There are $n\le r-S$ admissible choices. Each of the two delayed expressions permutes its new digit values. Cauchy–Schwarz therefore bounds the sum of relative weights after one step by
$$z^S\sum_{j=0}^{r-1-S}z^{-2j}.$$
For $r=3$ take $z=5/4$. The three bounds are at most $41/20<7/3$. For $r\ge5$ take $z=\exp(1/r)$. The geometric sum is at most
$$\frac{2e}{3\sqrt3(1-e^{-2/r})}\le\frac{e(r+2)}{3\sqrt3}\le\frac{2r+1}{3}.$$
The first bound maximizes $x(1-x^2)$ at $x=1/\sqrt3$ with $x=z^{S-r}$. The middle inequality uses $e^{2/r}\ge1+2/r$. The last follows from $e/\sqrt3<11/7$ and $(2r+1)/(r+2)\ge11/7$. The weighted prefix count now gives the assertion.

## Consequences for Wasserman's conjecture

Write $N=dQ+1$ with $Q=q^e$ a prime power. In a possible counterexample the [endpoint](all_weights_endpoint_reduction.md) occurs and one other row $A$ avoids $q$. Its parts are $aQ,bQ,cQ+1$ with $a+b+c=d$. Assign the prime powers of $N$ to the two rows avoiding their bases, giving $UV=N$ with $U$ assigned to $A$. The other row is $C=V(f,g,h)$ with $f+g+h=U$. The [endpoint assignment](endpoint_reduction.md) gives $U\ge3$ and $U\mid a,b$. Also $c\ge1$, since otherwise $A$ has a part equal to one. Thus all three slopes on each row are positive and $U\le d$. One of $A,C$ is binary carry-free, since the endpoint coefficient is even.

- **Uniform count with growing cofactors.** For $X\ge3$ and $Y\ge1$, at most $O(Y^3 X^\alpha)$ weights $N\le X$ with $N-1=dQ$, $d\le Y$, can be exceptions. The implied constant is absolute. For binary $A$ fix $d$ and one of $O(d^2)$ triples of slopes. Use $L=\lceil\log_2(X/d)\rceil$ for $Q<X/d$ and $D\le\log_2 d$. For binary $C$ do the same with $U\le Y$, the normalized triple, and $V\le X/U$, allowing one more digit at equality. In both cases the count per shape is $O(X^\alpha)$. Each parameter and shape determines its weight. This also bounds exceptions with $N-1$ having *any* prime-power divisor at least $(N-1)/Y$.
- **Along each fixed odd prime base, uniformly in the exponent interval.** For fixed $d$ and odd prime $q$, every interval of $K$ consecutive positive exponents contains at most $O_{d,q}(K^\alpha)$ possible exceptions at $d q^e+1$. The order of $q$ modulo $2^L$ is $\gg_q 2^L$, by [lifting the exponent](prime_power_ray_density.md). Choose $L=\log_2 K+O_q(1)$ with order exceeding $K$. For a fixed shape, the residues $Q\bmod2^L$ distinguish the exponents in this interval. In the binary $C$ case the residues of $V$ modulo $2^{L+v_2(d)}$ distinguish them as well by $UV=dQ+1$. Apply the prefix theorem and sum over the fixed finite lists of shapes.
- **Along the binary base.** At $d2^e+1$ every interval of $K$ positive exponents has $O_d(K^{1-\log(9/7)/\log d})$ possible exceptions for $d\ge3$. Now $U$ is odd, and $A$ must avoid some odd $r\mid U$ with $r\le d$. Apply the odd-prime bound at $Q=2^e$ and use $\operatorname{ord}_{r^L}(2)\gg_r r^L$. The exponent is at most $\log_r((2r+1)/3)\le1-\log(9/7)/\log d$. For $d\le2$ there is no assignment with $U\ge3$.
- **Three factors in the predecessor.** If $N-1$ consists of three full prime powers of distinct bases and one is a fixed $w$, the number of possible exceptions below $X$ is $O_w(X^{\log_2((1+\sqrt5)/2)})$. The [unique CRT row](endpoint_factors.md) $A$ has a part $kt$, where $t=(N-1)/w$ and $3\le U\mid k<w$. If $A$ avoids 2, use the two-summand theorem on $kt,(w-k)t+1$. Otherwise $C$ avoids 2 and its normalized sum $U<w$ is bounded, giving the stronger count $O_w(X^\alpha)$ for that case. Both counts allow the other two prime powers to vary without any upper bound.

## Two independent bases at every fixed cofactor

For each fixed $d\ge3$ the count of possible exceptional weights $dQ+1\le X$ improves to
$$O_{d,\epsilon}\!\left(X^{\max(0,\alpha-\log(9/7)/\log d)+\epsilon}\right)$$
for every $\epsilon>0$. Besides the binary test, one of the two affine rows must avoid an odd prime $r\le d$ from a fixed finite list. Use $r\mid U$ if $U$ has an odd factor. Otherwise $N$ is even and $d$ is odd, with $d\ge 2U+1$ by the sizes of $a,b,c$, so take any $r\mid d$. The endpoint coefficient is divisible by $r$ in either case.

Here are the details needed to combine the two counts. In a fixed row and base $p$, delete $j$ low digits of an accepted input $t$ and write $t_0=t\bmod p^j$. The remaining input $\lfloor t/p^j\rfloor$ gives the same three slopes with new offsets $\lfloor(a_i t_0+b_i)/p^j\rfloor$ bounded by $a_i+b_i$. Its row is carry-free by the absence of a carry from the deleted digits. Thus for each block length $h$, every aligned $h$-digit block of an accepted input belongs to a single alphabet of size at most $C2^{h\alpha}$ at 2, or $C((2r+1)/3)^h$ at $r$, after summing over the finitely many offset states. Leading zero padding is allowed. The integer sets with digits in these alphabets, in bases $2^h$ and $r^h$, have dimensions bounded by the corresponding logarithms of their alphabet sizes divided by $h\log p$.

Theorem B of Glasscock, Moreira and Richter bounds the upper counting dimension of the intersection of affine images of these integer missing-digit sets by the sum of their dimensions minus one, or zero if that is larger. The bases are multiplicatively independent. For different rows the parameters are related by $V=(dQ+1)/U$, an affine map with positive slope. For the same row use the identity map. Taking $h$ sufficiently large absorbs the fixed alphabet constants in $\epsilon$, and the sets of shapes and primes are finite. [Glasscock, Moreira and Richter, Theorem B](https://doi.org/10.1112/jlms.12902)

Run `python3 check_binary_affine_counts.py` for finite-prefix checks using factorial valuations and the sharp examples.
