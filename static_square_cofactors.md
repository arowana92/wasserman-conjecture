# Square exponents at fifteen cofactors, and cubes at ten

**The three-row conjecture holds at $N=1+Dt^2$ for every positive $t$ having at most two distinct prime factors and every**
$$
 D\in\{2,3,4,5,6,8,9,10,13,16,18,22,25,37,58\}.
$$
**It also holds at $N=1+10t^3$ for all such $t$.** The two prime powers in a square may have independently chosen positive even exponents, and in a cube independently chosen positive multiples of three. The [existing theorem at six](two_prime_powers.md) permits arbitrary exponents there.

## At most two possible numerators

Let $b=q^j$ be a prime power, $j\ge1$. Put $D=hb$ with $h=1$ or $2$, requiring $b$ odd for $h=2$. Let $N=1+DW$ with $W\ge1$, $\omega(W)\le2$. At a failure necessarily $q\nmid W$. If $r$ denotes the least positive inverse of $W$ modulo $b$, then some $x$ satisfies
$$
 x\equiv r\bmod b,\qquad 3\le x<hb,\qquad \gcd(x,N)\ge3. \tag{1}
$$
The sole possible numerator for $h=1$ is $r$. For $h=2$ the two possibilities are $r,r+b$. In fact $xW$ is a part of a row whose unitary weight scale $U\ge3$ divides both $N$ and $x$.

To prove this, use the [endpoint reduction](all_weights_endpoint_reduction.md). The row assignment cannot split the two varying bases, because after their full powers are removed the remaining cofactor is twice an odd prime power, a prime power, or one. Each forbids a split by the [four-factor criterion](two_prime_powers.md). With at most one varying base grouping is immediate. Intersect the zero positions of the assigned full powers to get a part $xW$ in row $A$, with $3\le U\le x<hb$. When $W=1$ just select a part. The companion has normalized positive entries of sum $U$. It cannot take $q$: two normalized entries would be positive multiples of $b$, impossible at sum $U<2b$. If $q\mid W$, the full-power intersection gives $b\mid x$. The only possible positive multiple below $hb$ is $b$, and $U\mid b$ contradicts $\gcd(U,b)=1$. Otherwise $xW$ is zero or one modulo $b$, with zero excluded for the same reason. This proves (1).

For $h=1$, $W\equiv1\bmod b$ immediately excludes a failure. Thus for every prime power $b$, taking $W=t^e$ with $e$ any positive multiple of the unit-group exponent $\lambda(b)$ certifies the entire family, including when $q\mid t$.

For $h=2$ and $W\equiv1\bmod b$ the sole remaining numerator is $b+1$, so it suffices that $\gcd(2W-1,b+1)=1$. For every odd prime power $b$ the family $W=t^e$ is certified whenever $e$ is a positive multiple of
$$
 \operatorname{lcm}\big(\varphi(b),\ \ell-1:\ \ell\text{ an odd prime dividing }b+1\big).
$$
Indeed if $q\nmid W$, we have $W\equiv1\bmod b$. At each odd $\ell\mid b+1$, $W$ is zero or one modulo $\ell$, so $2W-1$ is a unit. The period is at most $b(b+1)$ by the product bound.

## Checking squares

When $W=t^2$ and $q\nmid W$, (1) requires $x$ to be a unit square modulo $b$. Any odd prime $\ell\mid\gcd(x,N)$ must have $-hb$ a nonzero square modulo $\ell$ because $hb t^2\equiv-1\bmod\ell$. The only other way to have $\gcd(x,N)\ge3$ is for $4$ to divide it, which requires $4\mid x$ and $hb\equiv3\bmod4$. These requirements involve only $x<hb$ and its prime divisors.

The finite checks have a second formulation avoiding factorization or quadratic reciprocity. For each integer $3\le x<hb$ for which $xz^2\equiv1\bmod b$ has a solution, check $\gcd(x,hb y^2+1)\le2$ for every $0\le y<x$. These include all possible residues of $t$ modulo $x$. At $h=2,b=5$, for example, $x$ can only be $4,6,9$, and $10t^2+1$ is coprime to both 2 and 3. The kernel-checked tables use this formulation. Among the admissible prime powers $b\le1000$, the local test excludes exactly the following cases.

| $h$ | Certified $b$ |
| --- | --- |
| 1 | $2,3,4,5,8,9,13,16,25,37$ |
| 2 | $3,5,9,11,29$ |

For example the test leaves $x=9,\ell=3$ at $h=2,b=7$. [Lean verification](StaticSquareCofactors.lean) · [Independent Python verifier](check_static_square_cofactors.py)



## The cube at ten

At $b=5$, put $N=1+10W$. The odd scale $U<10$ is coprime to 5, so it must be 3, 7, or 9.

- If $U=3$, unitarity gives $3\parallel N$, so $W\equiv2$ or $5\bmod9$.
- If $U=7$, then $x=7$. The divisibilities $7\mid N$ and $5\mid7W-1$ give $W\equiv2\bmod7$ and $W\equiv3\bmod5$.
- If $U=9$, then $x=9$ and $9\mid W+1$. Write $P=3^{v_3(W+1)}\ge9$. Modulo $9P$, the residue of $9W$ is $9P-9$, and the residue of $W+1$ is a positive multiple of $P$. They sum to at least $9P$. This forces a base-3 carry in the cut $(9W,W+1)$ even though $A$ must avoid 3. [Prime-power carry lemma](three_factors_small_power.md)

Cubes modulo 9 are $0,1,8$ and cubes modulo 7 are $0,1,6$. This excludes the first two cases for every cube $W$. Run **lean StaticSquareCofactors.lean** with Lean 4.36 and **python3 check_static_square_cofactors.py** to reproduce the scale implications, residue exclusions, and square tables.
