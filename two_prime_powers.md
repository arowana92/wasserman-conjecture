# Two prime powers dividing the predecessor

Let $Q,R$ be the full powers of two distinct primes in $N-1$, and set $D=(N-1)/(QR)$. The [all-weight reduction](all_weights_endpoint_reduction.md) allows us to take one row to be $E=(1,1,N-2)$. In a counterexample let the two other rows be $A,C$. Assign each prime of $N(N-1)$ to a row avoiding it and let $U,V$ be the products of the full powers of $N$ assigned to $A,C$, respectively. Thus $UV=N$, $\gcd(U,V)=1$, $U,V\ge3$, and the parts of $A,C$ are multiples of $U,V$, respectively.

**Theorem.** If the bases of $Q$ and $R$ are assigned to different rows, $D$ must admit a factorization into four pairwise coprime positive integers with
$$D=abcd,\qquad (a+b)(c+d)<D. \tag{1}$$
This is impossible exactly when $D$ is 1, a prime power, twice an odd prime power, 12, or 15. The classification concerns (1), a necessary condition for the assignment. If the two bases go to the same row $A$, then some part of $A$ is $kQR$, with $U\mid k$ and $3\le U\le k<D$. In particular, in every case where (1) fails, a counterexample requires a unitary divisor of $N$ in $[3,D)$.

**Proof.** At a full power of $N-1$ the residues of an avoiding row are $0,0,1$ in some order. If $Q$ and $R$ are both assigned to $A$, at least one position is zero modulo both. Its positive value is at most $N-2<DQR$, proving the second statement.

For a split assignment suppose $A$ takes $Q$ and $C$ takes $R$. Two parts of $A$ have the form $UQx,UQy$ and two of $C$ have the form $VRz,VRt$, with all four multipliers positive. The remaining parts are positive multiples of their row scales, so $V\ge Q(x+y)+1$ and $U\ge R(z+t)+1$. Multiplying gives $(x+y)(z+t)<D$. Every full prime power $w$ of $D$ is assigned to a row in which at most one part is nonzero modulo $w$. At least one of the two displayed parts of that row is divisible by $w$. The factors $UQ$ and $VR$ are each coprime to $D$, so $w$ divides one of $x,y,z,t$. Allocate each $w$ to such a multiplier. The resulting four products $a,b,c,d$ satisfy (1).

Equivalently (1) requires $(1/a+1/b)(1/c+1/d)<1$. With at most one prime base in $D$ the product is greater than one. With two full powers $s\le t$ of distinct primes the only potentially successful pattern places them in separate bins of the same pair. Its test is $2(1/s+1/t)<1$. This fails precisely for $s=2$ and for $(s,t)=(3,4),(3,5)$. If $D$ has at least three distinct prime factors, put its largest full power and the product of all remaining powers in the first pair, and put $1,1$ in the second. The two entries in the first pair are at least 5 and 6, respectively, which suffices.

## A four-factor family

**For every pair of primes $p,q$ and positive integers $i,j$, the full conjecture holds at $N=6p^i q^j+1$.**

First let $Q=p^i,R=q^j$ have distinct bases exceeding 3. Here $D=6$, so the two bases go to the same row. With $3\le U\le k<6$ and $\gcd(U,6)=1$, we must have $U=k=5$. Write $t=QR$. The distinguished part is $5t$ and its complement is $t+1$, a multiple of 5. If $5^h\parallel t+1$, then $5t\equiv-5\bmod 5^{h+1}$ and $t+1$ has least residue at least $5^h\ge5$. Thus the addition of this part to its complement carries in base 5, which the row was required to avoid. This uses the [prime-power carry argument](three_factors_small_power.md). If the bases coincide or one divides 6, $N-1$ has at most two distinct primes, or it has exactly three with a full power 2 or 3. The [two-factor](endpoint_reduction.md) and [three-factor](endpoint_factors.md) results cover these cases.

More generally, for $N-1=2wQR$ with four full powers of distinct bases, **all three odd bases must go to the same row**. Apply (1) to each pair of odd bases, leaving twice the third odd power as $D$. With $w$ the smallest odd power, the scale of this row is a unitary divisor $U$ of $N$ satisfying $3\le U<2w$. If $U=2w-1$ is a prime power, $k=U$ is forced in the same-row statement, and the same carry argument excludes it.

Run `python3 check_two_prime_powers.py` to exhaust the four-bin criterion through $D=5000$. The checker tests 175,282 pairs of low-residue assignments, including 14,594 covering a further power. Direct binomial products identify 37 pairs that avoid both assigned prime bases in all digit positions.
