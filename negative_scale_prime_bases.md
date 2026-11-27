# Prime bases with a growing cofactor range at every odd exponent

For all sufficiently large $X$, there are **at least $X^{0.09}$ primes $q\le X$** for which Wasserman's full conjecture holds at every weight
$$
 1+dq^e,\qquad 1\le d\le (\log X)/5,\qquad e\ge1\text{ odd}.
$$
For these same primes all positive exponents work outside a single set of $O((\log X)^{0.8689})$ cofactors in that interval. Here $\log$ is natural. Neither conclusion sets an upper bound on the exponent.

## A negative residue for the scale

For an integer $D\ge5$ set $H=\lceil(D+1)/2\rceil$ and
$$
 M_D=\prod_{p\le H\text{ prime}}p^{\lfloor\log_p(D+1)\rfloor}.
$$
Choose a prime $q>D+1$ with $q\equiv-1\pmod{M_D}$. The following hold simultaneously for $1\le d\le D$.

- Every odd exponent works.
- Every even exponent works if the conjecture holds at $d+1$ (the cofactor $d=1$ works unconditionally).

In particular every exponent works throughout $d\le D$ when $D<10^9$, by the [finite theorem](full_range_1b.md).

To prove the odd case, suppose a counterexample has scale $Q=q^e$. The [endpoint assignment](small_weight_transfer.md#specializing-the-two-rows) supplies an avoiding row $(aQ,bQ,cQ+1)$ with positive slopes summing to $d$. A unitary factor $U\ge3$ of $N$ divides both $a,b$, and the row avoids every prime $p\mid U$. Choose one. Since $U\le(d-1)/2$, $p\le H$. Write $w$ for the largest power of $p$ at most $D+1$. Then $Q\equiv-1\pmod w$ and $pw>D+1$. Since $p\mid a,b$, the two part residues modulo $pw$ are exactly $pw-a,pw-b$. Their sum exceeds $pw$ because $a+b<d$. This carry contradicts the assignment. The argument also applies directly to arbitrary positive scales on a row satisfying these hypotheses and the congruence.

For an even exponent $Q\equiv1\pmod{M_D}$. Write $n=d+1$. As in the [positive-scale transfer](small_weight_transfer.md#specializing-the-two-rows), $U$ is a unitary divisor of $n$, and the two small rows have parts divisible by $U$ and $n/U$, respectively. The proof of this fact uses only primes $p\mid U$, all at most $H$. If $n$ is prime this is impossible since $3\le U<n$. If $d$ is prime then $N-1$ has at most two distinct prime divisors and the [endpoint theorem](all_weights_endpoint_reduction.md#consequences) already settles it. In all other cases every prime of $n(n-1)$ is at most $H$, since a prime exceeding $H$ cannot have a second positive multiple through $D+1$. Thus all prime-power layers required to transfer carry freedom to weight $n$ are present in $M_D$. Also $q>n$. A counterexample at the large weight would give one at $n$. The [global exceptional-weight bound](two_predecessor_congruences.md#the-unrestricted-three-row-count) now supplies the exceptional cofactor count.

## Counting distinct prime bases

By the prime number theorem, $\log M_D=(1/2+o(1))D$. Indeed the first powers contribute $\vartheta(H)\sim H$, while all higher powers together contribute $O(\sqrt D\log D)$ by summing over primes at most $\sqrt{D+1}$.

Use the bound $C m^{5.2}$ for the least prime in each reduced residue class modulo $m$, with an absolute constant $C\ge1$. For $0<c<2/5.2$ take $D=\lfloor c\log X\rfloor$ and let $R$ be the largest power of two with $M_D R\le(X/C)^{1/5.2}$. For large $X$ it exists and
$$
 R=X^{1/5.2-c/2+o(1)}.
$$
The class $-1\bmod M_D$ has exactly $R$ reduced lifts modulo $M_D R$, since $M_D$ is even. The least primes in these lifts are distinct, exceed $D+1$ since they are at least $M_D-1$, and are at most $X$. For $c=1/5$ the count exponent is $6/65>0.09$. This proves the stated result. [Xylouris, *On Linnik's constant*](https://arxiv.org/abs/0906.2749)

## A certified 447-digit example

At $D=2026$, the number $q=12 M_{2026}-1$ is prime. Thus all **2,026 cofactors and every positive exponent** work at this explicit base.

Here is a short integer primality certificate. Set $m=q+1$ and $t=45{,}001{,}160$. Write $T_j=z^j+z^{-j}$ formally for $z+z^{-1}=t$. The recurrences $T_{2j}=T_j^2-2$, $T_{2j+1}=T_j T_{j+1}-t$, with $T_0=2,T_1=t$, compute the traces modulo $q$. The checker verifies $T_m\equiv2\pmod q$ and $\gcd(T_{m/p}-2,q)=1$ at each of the 170 distinct primes $p\mid m$. The factorization is complete by the definition of $M_{2026}$ and $12=2^2\cdot3$.

For any prime divisor $r$ of $q$, take a root $z$ of $z^2-tz+1$ in $\mathbb F_{r^2}$. The trace equals 2 exactly when the corresponding power of $z$ is one. The checks force its order to be $m$. Frobenius either fixes $z$ or exchanges it with $z^{-1}$, so that order divides $r-1$ or $r+1$. Hence $q+1=m\le r+1$ and $r=q$. This proves primality.

Run `python3 check_negative_scale.py` for the certificate and the independent small exact-binomial tests of the negative carry.
