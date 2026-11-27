# When one large weight prime certifies several families

Consider finitely many weights
$$
 N_i=1+c_i\prod_{j=1}^m x_j^{a_{ij}},
$$
where $c_i\ge1$, the $x_j$ are prime bases, and each nonzero vector $a_i$ of nonnegative integers uses at most two bases. We classify when one arbitrarily large prime $p$ can certify the full trinomial conjecture simultaneously at all these weights by its first nonzero weight digit. For $p>\max_i c_i$ the exact condition used is
$$
 N_i\equiv p^{h_i}\pmod{p^{h_i+1}}\qquad(h_i\ge1).
 \tag{1}
$$
The [weight-prime theorem](predecessor_weight_gap.md) then implies that **every choice of three rows at each weight has a common prime divisor**. That divisor may depend on the rows. The prime $p$ provides the valuation condition proving the full assertion.

Let $A$ be the matrix with rows $a_i$. The following are necessary and sufficient for (1) to be achievable at infinitely many $p$.

* Every integer relation $\sum_i\lambda_i a_i=0$ satisfies $\prod_i(-1/c_i)^{\lambda_i}=1$ in $\mathbb Q$.
* The rows admit an ordered partition $I_1,\ldots,I_t$ into nonempty blocks. At step $j$, some rational linear functional on the exponent vectors is 1 on $I_j$ and 0 on all later blocks. Its values on earlier blocks are unrestricted.

In a partition that works, the blocks can be assigned **any strictly increasing positive valuation levels**, including arbitrarily large gaps. There are at most $\operatorname{rank} A$ blocks. For sufficiently large $p$, recognition requires only the first condition, the desired partition test, and solvability of the original monomial equations with targets $-1/c_i$ in $\mathbb F_p^\times$. Below we give the exact prime-base density at each such $p$. The set of successful weight primes then has positive relative natural density.

The base primes may already be restricted to any specified reduced progressions. Choose $p$ outside their moduli, then combine the new conditions by CRT. Each resulting tuple of reduced progressions gives the full conjecture for every choice of prime bases in those classes. Every exponent entry can also be increased by a multiple of $\varphi(p^{1+\max h_i})$ while keeping the support of its row of size at most two.

## Why the two conditions are necessary

Suppose arbitrarily large primes work. There are finitely many ordered partitions by the sizes of the $h_i$, so one occurs infinitely often. For any relation $\lambda$ one has
$$
 \prod_i c_i^{\lambda_i}=\prod_i (N_i-1)^{\lambda_i}.
$$
All terms are units modulo $p$. Reduction at infinitely many primes gives the first rational identity, including its sign.

At a step of the partition use a relation supported on the remaining rows and let $h$ be their smallest valuation for the particular $p$. Divide out the exact signs and constants in the identity just proved. Expand $\prod_i(1-N_i)^{\lambda_i}=1$ modulo $p^{h+1}$. Products of two weight terms vanish at this precision, also for negative $\lambda_i$ by unit inversion. This gives $p\mid\sum_{i\in I_j}\lambda_i$. Infinitely many such $p$ force the sum to be zero as an integer. Thus the block indicator on the remaining rows annihilates their rational relation space. It is a column-space vector, exactly the claimed linear functional. The rank drops when the selected rows are removed, since the functional vanishes on the later span.

This argument excludes all sufficiently large primes if a fixed partition fails: clear the denominators of a violated relation, or take a nonzero integer block sum. The possible exceptional primes divide a fixed nonzero integer. There are finitely many partitions.

## Lifting a successful partition

Put the exponent matrix in Smith form $S=UAV$ using unimodular integer matrices, with nonzero invariant factors $d_1,\ldots,d_r$. Row operations multiply and invert the monomial equations; column operations give bijective changes of variables in a multiplicative group of units. Apply the row operations to the targets, obtaining $T_j=\prod_i(-1/c_i)^{U_{ji}}$. Below the diagonal rank these numbers are exactly 1 by the first condition. The remaining equations are $y_j^{d_j}=T_j$.

There are infinitely many primes at which all these equations have unit roots. For an elementary proof, adjoin one algebraic root of each equation to $\mathbb Q$, choose an integral primitive element $\theta$ of this finite extension, and express the roots and their inverses in $\mathbb Z[\theta,1/M]$ for some integer $M$. If the extension is trivial, take $\theta=0$. Its monic minimal polynomial $f$ has a root modulo infinitely many primes. In fact, if a finite set of primes is prescribed, choose $b$ with $f(b)\ne0$ and a multiple $D$ of all those primes. Then $f(b+tDf(b))/f(b)$ is an integer congruent to 1 modulo $D$ and unbounded in absolute value. Some new prime divisor gives a root. Away from $M$ one can reduce $\theta$ to that root and obtain all required unit roots at once. [Primitive element theorem](https://pi.math.cornell.edu/~kbrown/6310/primitive.pdf)

At any such $p$ outside the primes dividing the invariant factors or the targets' numerators and denominators, the roots lift to any requested precision. Each derivative $d_j y_j^{d_j-1}$ is a unit, so the next digit is found by a linear congruence. Inverting the Smith operations gives initial unit inputs with $c_i x^{a_i}=-1$ modulo $p^{1+\max h_i}$ for every row.

At block $I_j$ write its functional as $v/L$ with $v\in\mathbb Z^m$ and $L>0$. Avoid the finitely many primes dividing these denominators. For its prescribed height $h$, multiply each current $x_\ell$ by $(1-p^h/L)^{v_\ell}$ in the units modulo $p^{1+\max h_i}$. At a later row the change is identically one because its dot product with $v$ vanishes. At a selected row the multiplier is $(1-p^h/L)^L\equiv1-p^h\bmod p^{h+1}$, giving (1). Earlier rows preserve their already fixed digits. This proves sufficiency and the freedom of the heights.

There is an exact density for the possible auxiliary primes. In the joint splitting field of the binomials $X^{d_j}-T_j$, take the fraction of Galois group elements which fix at least one root of each binomial. Outside finitely many primes, Frobenius has this property precisely when all the binomials have roots in the prime field: its permutation on distinct reduced roots is the $p$-power map. The identity belongs to this conjugacy-stable set. Chebotarev's theorem gives the stated fraction as the relative natural density. [Chebotarev's theorem, Theorem 28.9 and Remark 28.12](https://ocw.mit.edu/courses/18-785-number-theory-i-fall-2021/mit18_785f21_lec28.pdf)

## Counting base tuples and minimizing prime lists

For fixed attainable heights and sufficiently large $p$ as above, the relative density among independent tuples of prime bases is
$$
 \frac{\prod_{j=1}^r\gcd(d_j,p-1)}{(p-1)^r}\,
 p^{-\sum_{k=1}^{\max h_i} r_k},
 \qquad r_k=\operatorname{rank}_{\mathbb Q}\{a_i:h_i\ge k\}.
 \tag{2}
$$
The Smith diagonal and cyclicity of $\mathbb F_p^\times$ give the first factor. At higher precision the solutions form a coset of the homogeneous kernel. The unit group at an odd prime is the direct product of its order-$(p-1)$ subgroup and its principal units. On writing the latter as powers of $1+p$, the homogeneous conditions are $a_i\cdot y=0\bmod p^{h_i}$. Select independent rows successively in decreasing order of height. Every dependent row is a combination of earlier selected rows with denominators prime to $p$ outside finitely many minors, and so adds no constraint. The selected linear forms are jointly surjective at full precision since a maximal minor is a unit. Each selected row costs a factor $p^{-h_i}$, and their heights sum to $\sum_k r_k$. This proves (2) by counting reduced classes and using the prime number theorem in progressions.

For any prescribed partition, the column-space tests use rational arithmetic. To find one when the heights are free, select independent rows of the current matrix. Try the zero-one vectors on these basis rows and interpolate the rational functional on all remaining rows. Remove any nonempty group of ones from a successful zero-one vector. This greedy choice is valid: if a successful sequence existed before removal, its restrictions still work afterward, skipping empty groups.

For any sufficiently large lower cutoff on the available auxiliary primes, the minimum number needed by criterion (1) is exactly the minimum number of compatible subcollections partitioning the rows. A separate prime for each compatible subcollection suffices by CRT. In the other direction, assign each row to one prime that certifies it; beyond the finitely many exceptions, every resulting subcollection satisfies the two conditions. Single rows always qualify. For example, different multipliers $c_i$ on the same exponent vector cannot share arbitrarily large certificate primes, since their ratio would have to be one. Thus the one-prime-per-multiplier construction for [simultaneous progressions in the same product](all_exponent_cofactors.md) uses the fewest such primes once the cutoff is sufficiently high.

## Three weights requiring the same valuation

Take $X=sq$, $Y=s^2q^3$ and the three weights
$$
                 1+4X^2,\qquad 1+6XY,\qquad 1+9Y^2 .
$$
The relation on exponent rows is $(1,-2,1)$. Its dot product with a nonempty zero-one block indicator is zero only when all three entries are one. For $p>9$ the same reasoning already holds modulo $p$, so all three valuations must coincide. Such a $p$ works at any prescribed common height $h$ exactly when $p\equiv1\bmod4$. Necessity follows from the first square modulo $p$. For sufficiency lift a square root $u$ of $p^h-1$ modulo $p^{h+1}$ and set $2X=3Y=u$ there. The change from $(s,q)$ to $(X,Y)$ is invertible in the units, with $s=X^3/Y$ and $q=Y/X^2$. The fraction of eligible weight primes is $1/2$. At each one the base density for height $h$ is $2/((p-1)^2p^{2h})$.

At $p=13$ take $s=89^a,q=193^b$. For **every $a\equiv1\bmod12$, $b\equiv1\bmod78$**, all three weights satisfy the full conjecture. At $a=b=1$ they are $1\,180\,197\,317$, $5\,868\,816\,171\,461\,815$, and $29\,184\,105\,731\,699\,278\,143\,082$. They are all congruent to 13 modulo 169. The verifier exhausts the exponent residues modulo 12 and 156 and recovers precisely the two pairs $(1,1),(1,79)$ for this height-one certificate. The inputs here are powers of the two indicated prime bases, so the same weight-prime theorem applies.

## Families requiring arbitrarily many different valuation levels at one prime

Fix $m\ge2$ and positive integers $c_1,\ldots,c_m$. Write $X_i=c_i x_i$ and use the weights
$$
                 1+X_i\quad(1\le i\le m),\qquad
                 1+X_i X_{i+1}^2\quad(1\le i<m).
$$
Choose any prime $p>3$ exceeding all the multipliers $c_i,c_ic_{i+1}^2$. One such prime certifies **all $2m-1$ weights simultaneously**, but it must use strictly increasing valuations $h_1<\cdots<h_m$ on the first $m$ weights. Every such choice works by taking $c_i x_i\equiv p^{h_i}-1\bmod p^{h_i+1}$. This gives independent prime progressions for all $m$ bases and allows all their exponents to vary through $1\bmod\varphi(p^{1+h_m})$.

Indeed put $u=1+X_i$, $v=1+X_{i+1}$, and $h=\min(v_p(u),v_p(v))$. The next weight satisfies
$$
             1+X_i X_{i+1}^2\equiv u+2v\pmod {p^{h+1}} .
$$
The two original weights have first nonzero digit 1, so the resulting digit is respectively 1, 2 or 3 according as the first valuation is smaller, larger or equal. At $p>3$ none is zero. This proves both necessity and sufficiency of the strict ordering. The prime-base density at the specified heights is $(p-1)^{-m}p^{-\sum h_i}$.

For $m=2$, take $c_1=5,c_2=7$ and $p=251$. The primes $s=567059$ and $q=648352291$ give all three certified families
$$
      1+5s^a,\qquad 1+7q^b,\qquad 1+245s^a q^{2b},
      \qquad a,b\equiv1\pmod{15750250}.
$$
Their valuations at 251 are $1,2,1$. For fixed heights $1,2,1$ the exact density (2) is $1/(250^2\,251^3)$. Requiring valuation one at all three weights fails at every $p>3$.

The verifier uses Smith decomposition over the integers so that parity conditions on relations retain the full integer kernel. It compares (2) and the ordered partition test with exhaustive unit-tuple counts for 64 combinations of a system, a prime and a height pattern, checks the forced ordering through seven levels, and recursively certifies the example primes. Run with Python, SymPy 1.14 and python-flint.

~~~sh
python3 check_shared_weight_primes.py
~~~
