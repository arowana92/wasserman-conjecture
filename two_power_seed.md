# Checking two unbounded prime powers at one reference weight

Let $D\ge1$, $n=D+1$, and let $L_j$ be the exponent of the unit group modulo $\operatorname{lcm}(1,\ldots,j)$. The following finite test uses only $D$. If it excludes both assignments described below, then **every** weight
$$
 N=1+D t^e,\qquad \omega(t)\le2,\qquad L_{n^3}\mid e,\qquad e\ge1
 \tag{1}
$$
satisfies Wasserman's three-row conjecture. Here $t\ge1$ is arbitrary. The same exponent works simultaneously for any two prime bases and for independently chosen positive multiples of that exponent on each base.

Through $D=1500$, the test certifies **1,341 cofactors** for all bases. They include every $1\le D\le43$. If $\gcd(t,D+1)=1$, it certifies **1,452 cofactors** in this range. The first unresolved cases for unrestricted bases are $44,70,78,90,92,126$. These are survivors of necessary conditions. [Unbounded cofactor families and a shorter test at products of two prime powers](static_power_transfer.md)

Write $W=t^e$. Every power of a base in $W$, taken with its exponent in $W$, is one modulo each prime power at most $n^3$ of every *other* base. This follows from the definition of $L_{n^3}$ and the Chinese remainder theorem. We use a layer of size strictly above $Dn$ at each prime $p\mid Dn$ when needed. The first such $p$-power is at most $pDn<n^3$. The period is a sufficient uniform choice.

By the endpoint reduction, a failure has two nonendpoint rows which together avoid all primes of $N(N-1)$. Assign the primes to avoiding rows. Their weight scales $U,V$ satisfy $UV=N$, $\gcd(U,V)=1$, $U,V\ge3$, and divide all parts of their respective rows. At every assigned full predecessor power the parts are $0,0,1$ in some order. [Endpoint theorem](all_weights_endpoint_reduction.md) · [Prime assignments](two_prime_powers.md)

## Both bases in the same row

Call the row taking the bases of $W$ row $A$, of scale $U$. Intersect its zero positions at the full predecessor powers of these bases. One part is $xW$ with $U\mid x$, $3\le U\le x<D$, and $x$ is divisible by the full power in $D$ of each base shared with $W$. This includes one base, and when $W=1$ one may just select a part. Write the companion as $V(f_1,f_2,f_3)$ with positive $f_i$ of sum $U$.

Since $U<D$ is coprime to $W$, $W\equiv1\bmod U$. Thus $U\mid n$. In fact $U$ is a unitary divisor of $n$. For $p\mid U$, the weights $N,n$ agree modulo a $p$-power greater than $n$ and $v_p(U)=v_p(N)$. Set $v=n/U$. The reference cut and triple are
$$
 B_0=(x,n-x),\qquad C_0=v(f_1,f_2,f_3).
$$

At $p\nmid W$ dividing $Dn$, carry freedom in the actual cut transfers to $B_0$ if $A$ takes $p$, by reducing modulo a layer greater than $n$. If $C$ takes $p$, the triple transfers to $C_0$ since $p\nmid U$ and $V\equiv v$ at that layer. In either case the nonnegative reference entries sum to $n$, so all their digits are visible at the chosen layer.

At a base $p\mid W$, remove the full $p$-power in $W$ from $W$ to give $R$. Above this shift the actual cut at $xW$ gives a carry-free addition of $xR,(D-x)R$. Since $R\equiv1$ at the needed layers, the cut $(x,D-x)$ avoids $p$. If $p\mid D$, both cut entries have units digit zero, so adding one to the second entry preserves carry freedom. Consequently the reference conditions are

- $B_0$ avoids every prime of $U$.
- For each prime of $D$, at least one of $B_0,C_0$ avoids it.
- The primes of $v$ that $C_0$ fails to avoid must be bases of $W$, so there are at most two. The cut $(x,D-x)$ avoids each of them. Under $\gcd(t,n)=1$ there are none.

## Different bases in different rows

Write $W=QR$ for the two prime powers as they occur in $W$. Two parts of the first row are $UQx,UQy$ and two of the other are $VRz,VRt'$, with positive multipliers. As in the two-power assignment, $(x+y)(z+t')<D$. We keep $D$ as given, including when it shares primes with $W$.

For a pair $x,y$, let $g=\gcd(x,y)$ and let $\alpha=x/g,\ \beta=y/g$, sorted in increasing order. At a prime $p$ assigned to this row the following is a necessary reference test at a total $T$:
$$
 \exists h\ge1:\quad v_p(h)=v_p(g),\quad
 (h\alpha,h\beta,T-h(\alpha+\beta))
 \text{ has nonnegative parts and adds without a base-}p\text{ carry}. \tag{2}
$$
The first two parts are positive. The last may vanish.

For a regular prime $p\mid Dn$ with $p\nmid W$, use $T=n/p^{v_p(n)}$. Indeed $N\equiv n\bmod B$ with $B$ a $p$-power exceeding $Dn$. The scale of the assigned row has $p$-part $P=p^{v_p(n)}$. Its two displayed parts $mx,my$ have positive residues $i,j\bmod B$, since $p^{v_p(mx)},p^{v_p(my)}<B$ for $m=UQ$ or $VR$. These residues fit without carries into $n$. The congruence $xj\equiv yi\bmod B$ is exact equality since both products are smaller than $Dn$. Hence $(i,j)$ is an integral multiple of $(\alpha,\beta)$. Dividing all reference entries by $P$ gives (2), with the stated valuation of $h$ because the scalar $m/P$ is a $p$-adic unit. [Exact ratio lemma](TwoDigitWeightCuts.lean)

At a base $p$ of $Q$ or $R$ which divides $Dn$, its assigned row is the one with that displayed power. Above the shift, the two terms on the $Q$ side are $Ux,Uy$ in a total $DR$, and analogously on the other side. Reduction gives exactly the same ratio argument but with $T=D$ and a unit scalar. If $p\mid D$, adding one to the third reference entry transfers this triple to $T=n$, as all three units digits were zero. Thus all primes of $D$ still satisfy the regular test. For a prime of $n$ that is also a base of $W$, its assigned row may use $T=D$ instead of $n/p^{v_p(n)}$. Each row has at most one such substitution. There are none when $\gcd(t,n)=1$.

## Computing the least possible sums

Assign the primes of $Dn$ to two rows. For a fixed assignment and choice of substitutions, work separately on each row. For each of its primes and each coprime positive pair $(\alpha,\beta)$ tabulate the *least* possible $p^{v_p(h)}$ in (2), discarding the pair when no $h$ exists. Multiply the resulting costs over those primes. Their product is the smallest $g$ that could pass all local tests for that primitive pair, since they constrain valuations at distinct primes. Multiply by $\alpha+\beta$ and minimize over pairs. The minimum on an empty assignment is two. If the product of the two row minima is at least $D$ for every assignment, a split is impossible.

Run `python3 check_two_power_seed.py`. The program generates the regular tables by allocating base-$p$ digits. Through $D=150$ it checks the split test independently by enumerating all four multipliers and evaluating exact binomial products. It independently checks grouped assignments by exact coefficients through 150 and at the positive control 475. At $D=44$ the pairs $(1,1)$ and $(4,11)$ pass the local tests with product of sums 30. Through 1500 the only cofactor admitting a grouped reference pair is 475, where $x=76,\ U=4,\ (f_1,f_2,f_3)=(1,1,2)$ permits the two missing weight primes 7 and 17. The program prints the full list of surviving cofactors. [Standard-library verifier](check_two_power_seed.py)
