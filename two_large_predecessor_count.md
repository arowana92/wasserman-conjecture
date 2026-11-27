# Counting exceptional weights with two large prime-power factors

**Theorem.** Put $\beta=\log_2((1+\sqrt5)/2)$. For $X\ge3$ and $1\le Y\le X$, at most
$$O\!\left(X^\beta Y^{3-\beta}\right)$$
weights $N\le X$ can fail Wasserman's trinomial conjecture and have a representation $N-1=dQR$ with $d\le Y$ and $Q,R$ prime powers of positive exponent. The constant is absolute. The primes and exponents may all vary, and the two bases may coincide.

In particular for each fixed positive $d$ the exponent is $0.6942419\ldots$ independently of the two prime powers. Even allowing every cofactor $d\le X^{1/16}$ at once gives $O(X^{0.839})$ exceptional weights. With a single prime-power factor, the bound is $O(X^\alpha Y^{3-\alpha})$, where $\alpha=\log_2\lambda_3=0.551463\ldots$ and $\lambda_3^3=\lambda_3^2+1$.

## Summing over the coefficients of affine parts

For $j=2,3$, let $\gamma_j$ be $\beta,\alpha$ respectively. For a list of positive integers $a_1,\dots,a_j$ of sum $L$, put $\Delta=\max v_2(a_i)-\min v_2(a_i)$ with $a_0=L$ included in both extrema. For fixed nonnegative offsets $b_i$, the number of positive integers $t\le X/L$ for which the $a_i t+b_i$ add without binary carries is
$$O\big((X/L)^{\gamma_j}2^{\gamma_j\Delta}\big)\qquad (L\le X).$$
Indeed use a binary modulus of size strictly greater than $X/L$ and at most $2X/L$ in the affine prefix bound. [Binary affine count](binary_affine_counts.md)

For $Z\le L<2Z$, $Z\ge1$, the sum of $2^{\gamma_j\Delta}$ over all these ordered slope lists is $O(Z^j)$. If $\Delta\ge h$, some $a_i$ or $L$ is a multiple of $2^h$. For $2^h\le2Z$ there are $O(Z^j/2^h)$ lists satisfying this condition. For an individual slope fix the other $j-1$ slopes, then count its possible multiples. For the sum, fix the other $j-1$ slopes and count the last in one residue class. Larger $2^h$ cannot occur. Weight these tail bounds by $2^{h\gamma_j}$ and sum, using $\gamma_j<1$. Summing also over dyadic $Z$ gives
$$\sum_{L\le Y}\ \sum_{a_1+\cdots+a_j=L}
     (X/L)^{\gamma_j}2^{\gamma_j\Delta}
   \ll X^{\gamma_j}Y^{j-\gamma_j}. \tag{1}$$

At an exception use the endpoint row $E=(1,1,N-2)$. Assign each prime of $N(N-1)$ to one of the other two rows avoiding it. Write $UV=N$ for the resulting coprime scales, so the two rows have parts divisible by $U,V$ respectively, and both scales are at least 3. The endpoint coefficient $N(N-1)$ is even, so at least one of the two rows adds without binary carries. [All-weight endpoint reduction](all_weights_endpoint_reduction.md) · [Factor assignments](endpoint_reduction.md)

## One power and two powers assigned to one row

First consider $M=N-1=dQ$ with $d\le Y$ and $Q$ a prime power. A row $A$ avoiding its base has parts $aQ,bQ,cQ+1$, with positive $a,b,c$ and $a+b+c=d$. Positivity of $c$ follows since a part equal to one is incompatible with scale $U\ge3$. Also $U\mid a,b$ by coprimality, so $U\le d$. If $A$ avoids 2, use (1) with three slopes and $t=Q$. Otherwise the other row $C$ is $V(f,g,h)$ with $f+g+h=U\le Y$. Use (1) for those slopes and $t=V$, which determine $N=UV$. This proves the single-power bound, allowing $Q$ to divide $d$ as well.

For a representation involving two powers of distinct bases, absorb any remaining factors of these bases from $d$ into $Q,R$, so they are full powers and the cofactor only decreases. Coincident bases are already covered by the single-power bound. If one row $A$ takes both bases, the zero positions modulo the two powers overlap. A positive part at their intersection is $kQR$ with $3\le U\mid k<d$. If $A$ avoids 2, its cut $kt,(d-k)t+1$ is binary carry-free at $t=QR$. Use (1) for two slopes of sum $d\le Y$. If $C$ avoids 2 use the same three-slope argument as above at $U<d$. These bounds are both $O(X^\beta Y^{3-\beta})$ since $\alpha<\beta$ and $Y\le X$. [Two-power assignment](two_prime_powers.md)

## When the two powers go to different rows

Name a binary carry-free row $A$ and the full power assigned to it $Q$. The other row $C$ takes $R$. Two positive parts of $A$ have the form $mx,my$ with $m=UQ$. In $C$, two normalized parts are positive multiples of $R$ and the three normalized parts sum to $U$. Thus $U\ge2R+1$ and $m>2QR$. Since $m(x+y)\le N-1=dQR$, we have $x+y<d/2$. The binary digits of $mx,my$ are disjoint, and $m\le X/(x+y)$. Count possible $(x,y,m)$ by (1) to get $O(X^\beta Y^{2-\beta})$.

It remains to count how many weights can come from one such triple. We may take $Y\le\sqrt X$, since for larger $Y$ the claimed bound exceeds $X$. First count all weights with $Q\le(X/Y)^{1/20}$ by the single-power result applied to $R$, with cofactor $dQ\le Y(X/Y)^{1/20}$. Their count is at most
$$O\!\left(X^\alpha Y^{3-\alpha}(X/Y)^{(3-\alpha)/20}\right)
   \ll X^\beta Y^{3-\beta}.$$
The exponent comparison is $\alpha+(3-\alpha)/20<\beta$. For example $\alpha<0.552$ and $\beta>0.694$ follow by bounding their defining roots.

Now $Q>(X/Y)^{1/20}\ge X^{1/40}$. Since $\gcd(U,Q)=1$, $Q$ is the full power of its base in $m$. The eligible full powers of $m$ are pairwise coprime, so their product is at most $m\le X$ and there are at most 40 choices for $Q$. Each determines $U=m/Q$. Then $N\equiv0\bmod U$ and $N\equiv1\bmod Q$ determine one class modulo $m$. The inequality $N=1+dQR<1+Ym/2$ gives $O(Y)$ weights in this class. Multiplication by the preceding count of $(x,y,m)$ completes the proof.

## A consequence for any odd trinomial coefficient

The split argument gives a restriction independent of the number of prime factors of the predecessor. For $1\le Y\le X$ and every $\varepsilon>0$, at most $O_\varepsilon(X^{\beta+\varepsilon}Y^{3-\beta})$ weights $N\le X$ admit an odd trinomial coefficient with two parts $a,b$ satisfying
$$\gcd(a,b,N(N-1))\ \ge\ (N-1)/Y.$$
Indeed take that gcd as $m$, put $a=mx$, $b=my$, and use $x+y\le Y$ and the disjoint bits of $mx,my$ in (1). For fixed $m$, each full prime power dividing $m$ divides either $N$ or $N-1$ since $m\mid N(N-1)$. This gives at most $2^{\omega(m)}=O_\varepsilon(X^\varepsilon)$ classes for $N$ modulo $m$, with $O(Y)$ representatives per class up to $1+Ym$. The elementary estimate follows by charging a factor $2$ to $p^\varepsilon$ for every $p\ge 2^{1/\varepsilon}$ and absorbing the smaller primes into the constant. For example, for every $\theta>2/(3-\beta)$ at almost all weights every odd trinomial satisfies $\gcd(a,b,N(N-1))<N^\theta$ for each pair of parts. Apply the count on dyadic ranges with $Y=O(X^{1-\theta})$ when $\theta<1$. The assertion for $\theta\ge1$ follows directly from positivity.

Run `python3 check_two_large_count.py` to check the delay tails for 41,664 ordered lists of slopes and the split inequalities and CRT reconstruction for 324,116 oriented pairs of assignments. Direct binomial products identify 74 pairs that avoid both assigned primes, including 12 with an odd coefficient on the first row.
