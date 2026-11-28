# A gap above the lowest base-$p$ digit at a three-prime predecessor

Write $N-1=PW$ with $P$ a prime power and $W\ge1$ having at most two distinct prime factors. Take any prime $p\mid N(N-1)$, let $M$ be the member of $\{N-1,N\}$ divisible by $p$, and put $R=p^{v_p(M)}$.

**The full three-row conjecture holds at $N$ if $M/R\equiv1\bmod A$ for some power $A=p^m\ge P$ with $m\ge1$.** This means either $M$ is a power of $p$, or its lowest occupied digit is one and the next occupied position is at distance at least $\lceil\log_p P\rceil$.

In particular, in a counterexample with exactly three distinct prime factors of $N-1$, any gap $M/R\equiv1\bmod A$ would require $A$ to be strictly smaller than **each** full prime power in $N-1$. Hence $A^3<N-1$. A full prime power is the largest power of its prime dividing $N-1$. The prime $p$ may divide the weight itself.

## The cut and its companion

Use the [one-numerator reduction](static_square_cofactors.md#at-most-two-possible-numerators). At a failure $\gcd(P,W)=1$ and the least positive inverse $x$ of $W$ modulo $P$ occurs as a part $xW$ of a row $B$. Here $1\le x<P$ and $B$ avoids all primes of $W$. There are coprime positive scales $UV=N$ with $3\le U\mid x$ and a companion row $V(g_1,g_2,g_3)$ where the $g_i$ are positive and sum to $U$. Each row avoids all the primes of its own scale ($U$ for $B$, $V$ for the companion). Avoidance means that the prime does not divide the coefficient. By Kummer, each part of an avoiding row uses only digits of the total, with at most the available multiplicity at each position. Write $e\preceq_p d$ for this digitwise inequality. [Bergman, Lemma 5](https://arxiv.org/abs/0806.0607)

The following version allows a larger low residue. Put $d=(M/R)\bmod A$. It lies in $[1,A)$ because $M/R$ is a unit at $p$. Suppose
$$
                              dP\le A. \tag{1}
$$

**When $p\mid N-1$ and $p$ is not the base of $P$**, failure forces $P\mid d$. Indeed $W=RT$ and $PT\equiv d\bmod A$. Above the bottom $v_p(R)$ digits, the digits of the part $xRT$ fit inside those of $PT$. This follows by taking the integer quotients of $N=1+RPT$ and $xRT$ by $R$. Thus $xT\bmod A=e\preceq_p d$. We have $xd\equiv Pe\bmod A$. The left side is positive and below $A$. The right side is at most $A$, and equality to $A$ would contradict the congruence. Hence $xd=Pe$, and $\gcd(x,P)=1$ gives $P\mid d$.

**When $p\mid N$**, one of the following conditions is necessary for failure.

1. $P\mid Rd-1$ and there is an integer $0<e<d$ with $e\preceq_p d$ and $(Rd-1)x=RPe$.
2. The base-$p$ digit sum of $d$ is at least three, $\gcd(d,U)\ge3$, and the actual companion parts are $N e_i/d$ for three positive integers $e_i$ adding to $d$ without a carry.

If $p\mid U$, row $B$ avoids $p$, so $xW\bmod RA=Re$ for some $e\preceq_p d$, since $N\bmod RA=Rd$. Using $PW=N-1$ gives $(Rd-1)x\equiv RPe\bmod RA$. By (1) both sides are at most $RA$ and the left side is strictly between zero and $RA$. They are equal. If $e=d$ the right side would exceed the left, since $x<P$. This proves the first condition.

If $p\mid V$, then $R\mid V$ and $z=V/R$ is a unit at $p$. Removing the common power $R$ gives a carry-free companion with parts $z g_i$ and total $zU=N/R$. Its residues $e_i=z g_i\bmod A$ therefore add to $d$ without a carry. For each $i$,
$$
                    d g_i\equiv U e_i\bmod A .
$$
Both sides lie in $[0,A)$ since $g_i<U\le x<P$. Thus $d g_i=U e_i$ and $e_i>0$, giving the asserted parts and digit sum. Each $g_i$ is a multiple of $U/\gcd(d,U)$. Their sum is $U$, giving $\gcd(d,U)\ge3$.

These implications give the gap theorem with $d=1$. In the predecessor case with $p$ the base of $P$, coprimality implies $P=R$ and $W=M/R\equiv1\bmod A$. Then $P\mid A$ gives $x=1$, contradicting $U\ge3$. The other cases exclude $d=1$ directly.

**Additional consequences.** In the weight case (1) suffices when $P\nmid Rd-1$ and either $\gcd(d,N)\le2$ or $s_p(d)\le2$. In the predecessor case with distinct prime bases for $P$ and $R$, (1) suffices whenever $P\nmid d$, with no restriction on the digit sum of $d$. In particular every $d<P$ works there if $dP\le A$.

For example, take $P=11$, $p=5$, $R=5^i$, and any odd prime $q\equiv23\bmod125$. Let $j\equiv1\bmod100$ with $i,j\ge1$. Then $11q^j\equiv3\bmod125$ by Euler's theorem. For $N=1+11\cdot5^i q^j$ we have $d=3$, $3\cdot11\le125$, and $11\nmid3$. Thus all these weights satisfy the conjecture. Both exponents can grow independently.

The [binary parity theorem](binary_run_packing.md) gives complementary information when $P$ exceeds the gap modulus. Run **python3 check_low_digit_inverse_cuts.py** to check the modular implications and **lean +leanprover/lean4:v4.36.0 LowDigitInverseCuts.lean** for the congruence and proportion lemmas.
