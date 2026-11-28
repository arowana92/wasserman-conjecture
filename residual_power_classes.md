# Unrestricted binary and ternary exponents at cofactors 7, 11 and 19

For every $D\in\{5,7,11,19\}$ and **all** positive $a,b$, the three-row conjecture holds at
$$
                              N=1+D\,2^a3^b.                               \tag{1}
$$
Only three binary digits and one ternary digit are needed in the cut test. For the [two-prime three-unit theorem](three_unit_layers.md), a prime-power common remainder now gives only finitely many possible exceptional weights at each fixed base pair, even as the cofactor varies. If the cofactor $D$ itself is a prime power, all 23 outstanding small-base pairs are settled.

## Proof of (1)

In a counterexample the [endpoint reduction](all_weights_endpoint_reduction.md) and the [assignment for three full predecessor powers](endpoint_factors.md#when-n-1-has-three-distinct-prime-factors) give a row avoiding all three predecessor primes. Write $W=2^a3^b$. The row has a part $xW$, where $1\le x<D$ and $xW\equiv1\pmod D$. Every part is divisible by a unitary divisor $U\ge3$ of $N$. Since $(W,N)=1$ and $N$ is prime to six, an odd prime $\ell\ge5$ divides $U,x,N$. Thus $W\equiv-D^{-1}\pmod\ell$.

At 2, above the shift by $2^a$, the numbers $x3^b,(D-x)3^b$ add without a carry. At 3 the corresponding addition uses $x2^a,(D-x)2^a$. This follows by grouping the two remaining parts as $(D-x)W+1$. Let $A,B\in\{0,1\}$ be the parities of $a,b$. Modulo eight the binary summands use $3^B$, and modulo three the ternary summands use $2^A$. Also $W=2^A3^B y^2$. Consequently $x2^A3^B$ must be a nonzero square modulo $D$, and $-D2^A3^B$ a nonzero square modulo $\ell$.

At $D=5$ there is no eligible $x$. At $D=7$ only $x=5$ is eligible. The binary cut forces $b$ even, so $W$ is a square modulo seven since 2 is a square there. The required inverse of 5 is a nonsquare.

At $D=11$, the eligible multipliers are $5,7,10$. For $x=5,7$ the two summands overlap in their lowest three bits for both $3^b\equiv1,3\pmod8$. At $x=10$, the binary cut forces $b$ even. In ternary the two units digits are both $2^a\bmod3$, forcing $a$ even. Now $W$ is a square and $10W\equiv1\pmod{11}$ would require $W\equiv-1\pmod{11}$.

At $D=19$ the binary test discards $x=5,7,13,14,15$. It leaves $10,11,17$.

* At $x=10$, $\ell=5$ and $W\equiv1\pmod5$. Both 2 and 3 are nonsquares at 5 and at 19. Thus $a+b$ is even, although $W\equiv2\pmod{19}$ requires it odd.
* At $x=11$ the ternary test gives $a$ odd. But $W\equiv4\pmod{11}$, and 2 is a nonsquare and 3 a square at 11.
* At $x=17$ the ternary test gives $a$ odd and the binary test gives $b$ even. This makes $W$ a nonsquare at 19, whereas $17$ is a square there.

The Lean kernel checks all four finite parity certificates. The Python verifier independently walks the full power residues modulo $3D\ell$ and $8D\ell$ without using quadratic characters. [Lean certificate](ResidualPowerClasses.lean) · [Python verifier](check_residual_power_classes.py)

## An arbitrary moving prime

For every prime $r$, every $a,b\ge1$, and each cofactor $D$ in a row of the table, the conjecture holds at $1+D\,2^a r^{Jb}$. Each row lists the additional cofactors certified at that exponent, retaining those from earlier rows.

| $J$ | Additional $D$ |
|---:|---|
| 2 | 3, 5, 7, 9 |
| 4 | 17, 49, 81 |
| 8 | 25, 73, 89, 97, 233 |
| 16 | 113, 193 |
| 32 | 137 |

Indeed let $D=d^f$ be an odd prime power. A repeated prime base leaves at most two predecessor primes and is already settled. Otherwise the same assignment gives $x$, $U$ and an odd prime $\ell\mid (U,x)$ with $\ell\nmid Dr$. Every odd $J$th power is 1 modulo $4J$ for the displayed powers of two, by squaring starting from odd squares modulo eight. Thus the lowest $\log_2(4J)$ bits of $x$ and $D-x$ must be disjoint. With $A=a\bmod J$, both $x2^A$ modulo $d$ and $-D2^A$ modulo $\ell$ must be $J$th powers of units, by $xW\equiv1\pmod D$ and $DW\equiv-1\pmod\ell$. Trying $1\le x<D$, the eligible prime divisors of $x$ and $0\le A<J$ gives the table among odd prime powers through 250. These tests can be applied to any odd prime power $D$.

## Two sparse layers with a power cofactor

Take $N=1+Dp^a q^b$, $p<q$ distinct primes, $(D,pq)=1$. Assume at each of $p,q$ there is a layer as follows. At $p$ let $h_p=Dq^b\bmod L_p$ where $L_p$ is a positive power of $p$, $Dh_p\le L_p$, and $h_p$ has at most three base-$p$ digit units. Use $Dp^a$ at $q$ with the same conditions in base $q$. The [two-prime theorem](three_unit_layers.md) already settles all cases except $h_p=h_q=s\mid D$. In that case write $w=s-1$, $E=D/s$, $W=Ep^a q^b$. A failure requires a row with part $wW$, complementary sum $W+1$, and a unitary divisor $3\le U\mid w$ of $N=sW+1$ dividing all its parts. Its companion is $(N/U)(i,j,k)$ with $i+j+k=U$ and positive entries.

**For $p,q<257$, every eighth power $E$ suffices.** For a square $E$ the only remaining pairs are $(2,17),(2,71),(2,191)$, and for a fourth power only $(2,71)$ remains. In each case $s=q+2$, $a$ must be odd and $b$ even. The restrictions here are necessary conditions for failure.

Suppose $E=z^J$, $J=2,4,8$, and write $A=a\bmod J$, $B=b\bmod J$, $\sigma=p^A q^B$. Then $W$ is $\sigma$ times a $J$th power. We use the following congruences.

1. Since $p\nmid s$, $E q^b\equiv1\pmod {L_p}$, so $q^B$ must be a unit $J$th power modulo $p$. In binary use $q^B\equiv1\pmod{4J}$ when $L_p\ge4J$. The latter follows from $L_p\ge Ds\ge s^2$ in every binary table case. At $q$ use the corresponding condition on $p^A$.
2. If $r$ is an odd prime divisor of $U$, then $W\equiv-1\pmod r$ because $U\mid W+1$. Hence $-\sigma$ must be a $J$th power modulo $r$. The same holds for every odd $r^e\parallel s$ with $2r^e\ge U$. In fact the companion cannot take $r$, since two positive normalized entries would have to be multiples of $r^e$ and their sum is less than $U$. If $r\mid W$, the valuation of $wW$ is strictly between zero and $v_r(N-1)$, contradicting the residues $0,0,1$ at the assigned full power. Otherwise $wW$ is the unit position and $W\equiv-1\pmod{r^e}$.
3. If $u$ is the full 2-part of $U$ and $u>1$, then $p,q,z$ are odd. Unitarity gives $sW+1\equiv u\pmod{2u}$ and $W\equiv-1\pmod u$. Test these against $W\equiv\sigma y^J\pmod{2u}$ for an odd $y$.

Starting with the 25 possible $(p,q,w)$ triples in the earlier theorem, enumerate $U\mid w$, $U\ge3$, $(U,pqs)=1$ and $A,B\bmod J$. The following are the surviving $(p,q;w;U)$, collecting the possibilities for $U$ in braces.

| $J$ | After the congruences |
|---:|---|
| 2 | $(2,17;18;\{3,9\}),(2,71;72;\{3,9\}),(2,191;192;3),(5,29;30;3),(5,149;150;3),(13,181;182;7)$ |
| 4 | $(2,71;72;\{3,9\}),(13,181;182;7)$ |
| 8 | $(13,181;182;7)$ |

The verifier compares a cyclic-group power test with direct enumeration of modular powers. It also checks the parities. They are $(0,1)$ at $p=5$ and $(1,0)$ in the other surviving rows.

At $(5,29)$ and $(5,149)$ the companion with $U=3$ has three equal parts, so its coefficient is even. The grouped row must avoid 2 since the endpoint coefficient is even. If $W$ is odd, $N$ is even and the companion's scale $N/3$ is even, forcing the companion itself to avoid 2. If $W$ is even, strip its full binary power to get the odd integer $m$. The required addition at 2 includes $wm+m$ above the shift. The odd part of the square $E$ is 1 modulo eight, both moving primes are 5 modulo eight and $a+b$ is odd. Thus $m\equiv5\pmod8$ and $wm\equiv6\pmod8$, which share a bit.

At $(13,181)$ we have $s=183$, $w=182$, $U=7$. If $3\mid W$, the grouped row cannot take the full predecessor power of 3 by the valuation argument above. That power is at least nine and is too large for the companion of normalized sum seven. Otherwise $W\equiv1\pmod3$ because $E$ is a square and both bases are one modulo three. The distinguished part is 2 modulo three. Thus the companion takes 3, and its normalized entries must be $3,3,1$. But $V=(183W+1)/7\equiv7\pmod9$, so the residues of $3V,3V,V$ modulo nine are $3,3,7$, giving a carry. This removes the last pair at an eighth power.

### The cofactor classes at bases 2 and 3

At $(p,q)=(2,3)$ the only possible common remainder is $s=11$, and $U=5$. Of the 440 residue classes modulo 1320 prime to six, **at most 40 can contain $E=D/11$ in a failure**. The smallest positive representative is 43.

Indeed $W=E2^a3^b$ must satisfy $W\equiv-1\pmod{11}$ by condition 2 and $W\equiv-1\pmod5$ since $U\mid W+1$. In particular $(E,55)=1$. The layers give $E2^a\equiv1\pmod3$ and $E3^b\equiv1\pmod8$. Thus $(-1)^a=(E/3)$ and $E\bmod8$ is 1 or 3 according as $b$ is even or odd. Put $\delta=1$ in the first case and $-1$ in the second. Taking Legendre symbols of $W$ at 11 and 5 gives
$$
                      (E/3)(E/11)=-1,\qquad (E/3)(E/5)\delta=1.
$$
There are two choices for $E\bmod8$, two for $E\bmod3$, two for $E\bmod5$ and five for $E\bmod11$ under these restrictions. The verifier also enumerates exponent residues directly using the four congruences for $W,E2^a,E3^b$.

## When the third factor is a prime power

Under the two-layer hypotheses, **every weight with $D$ a prime power satisfies the conjecture if either $p,q<257$ or $D\le1032$**. In the second alternative the bases are unrestricted.

Indeed (2) of the [layer theorem](three_unit_layers.md) requires $s\mid D$, so $D$ and $s$ are powers of the same prime $r$. Since $U\le s-1$, the companion cannot take $r$. If $s<D$, then $r\mid E$ and the grouped row also cannot take $r$ by the valuation of $wW$. Thus $D=s$ and $E=1$, an eighth power. This proves the first alternative.

For the second, enumerate the prime powers $s\le1032$. For each choose all prime bases in which $w=s-1$ has precisely two digit units, requiring $s$ prime to each base. Each eligible base is at most $w$ since $w=2$ cannot have $U\ge3$. Take the pairs for which the part of $w$ coprime to $pq$ is at least three. There are 17 triples $(p,q,w)$, of which only $(3,269,270)$ and $(3,809,810)$ are outside the previous table. None of the 17 pass the eighth-power congruences. The calculation enumerates only $s$ and the two-unit representations and places no bound on either exponent of the weight.

At arbitrary $D$ there is also a size criterion in this prime-power case. At least one of $p,q$ divides $D-1=w$. Otherwise both sums of two units must have a units digit and $w-1$ would be a positive power of both bases. Repeated equal units are impossible because $w=2p^i$ would leave at most two after stripping $p$. Thus a base dividing $w$ contributes two distinct positive-exponent powers. In particular $D\ge 1+t+t^2$ with $t=\min(p,q)$ is necessary. Any moving base not dividing $D-1$ has to be the base of a prime power $D-2$.

## An effective exponent bound at a prime-power remainder

Keep the two-layer hypotheses and allow arbitrary $E=D/s$. If $s=r^d$ is a prime power, every possible exception satisfies
$$
 a+b\ \ll_{p,q,s}\ \log(2E)\big(1+\log\log(3E)\big)^{d(r-1)}.             \tag{2}
$$
The constant is effective. In particular, for **any fixed** moving primes $p,q$, only finitely many possible exceptions with a prime-power $D$ and the two layers remain, allowing the base and exponent of $D$ to vary. There is an effective finite search for that set.

As before, $r$ must go to the grouped row, $r\nmid W$ and $s\mid W+1$. From carry freedom of $wW,W+1$ in base $r$, the digits of $W$ decrease weakly along each residue class of positions modulo $d$. Indeed, up to position $v_r(W+1)-1$ all digits of $W$ are maximal. Above position $v_r(W+1)$ the digits of $W+1$ equal those of $W$ and must fit in the digits of $sW+1$, which are shifted by $d$. At positions from $d$ through the transition, the comparison digit $d$ places earlier is maximal. Splitting each column into strings of ones therefore gives
$$
           w(W+1)=\sum_{i=1}^m c_i s^{k_i},\qquad
           m=d(r-1),\quad k_1\ge\cdots\ge k_m\ge1,
           \quad c_i\in\{1,r,\ldots,r^{d-1}\},\quad \sum c_i=w.
$$
This is the [monotone-digit representation](predecessor_five_seven.md#finiteness-with-prescribed-prime-bases).

Put $k_{m+1}=0$, $F_j=\sum_{i\le j}c_i s^{k_i-k_j}$ and $B=\max(3,a,b,k_1)$. For $1\le j\le m$ define
$$
          \Gamma_j=\frac{wE}{F_j}\,p^a q^b s^{-k_j}-1.
$$
Its numerator over $F_j s^{k_j}$ is the omitted tail minus $w$. Each nonempty tail contains a term at least $s>w$, and at the last step the numerator is $-w$. Hence $0<|\Gamma_j|\le w s^{k_{j+1}-k_1}$. Apply Matveev's theorem for at most four positive rational factors, omitting a factor equal to one. The changing rational has logarithmic height at most $\log(wE)+(k_1-k_j)\log s$, since $F_j\le w s^{k_1-k_j}$. The other heights are fixed. Comparing the bounds gives, with an effective $C_{p,q,s}$,
$$
   (k_1-k_{j+1})\log s\le C_{p,q,s}\log(2B)
       \big(\log(2E)+(k_1-k_j)\log s\big).
$$
Starting at $j=1$ and iterating $m$ times bounds $k_1$ by $O_{p,q,s}(\log(2E)(\log(2B))^m)$. Since $W+1\le s^{k_1}$, we have $a+b\ll_s k_1$ and $B\ll_s k_1$. To invert the bound, first use $(\log(2B))^m\le\sqrt B$ above an effective threshold depending on $m$. This gives $B\ll_{p,q,s}(\log(2E))^2$. Substitution into the logarithm now gives (2). [Matveev's theorem in rational form, Languasco, Luca, Moree and Togbé, Theorem 2.1](https://doi.org/10.1007/s12188-025-00293-9)

More generally fix a finite set of primes $T$ and let $F$ be the part of $E$ left after removing all its prime factors in $T$ with multiplicity. Then
$$
        \log(2W)\ll_{p,q,s,T}\log(2F)\big(1+\log\log(3F)\big)^{d(r-1)}
$$
with an effective constant. In the same linear form, replace $wE/F_j$ by $wF/F_j$ and one fixed rational factor for each prime of $T$ occurring in $E$. Increase $B$ to include their exponents. They are all at most $\log_2 W\ll_s k_1$, so the iteration and inversion apply with $F$ in place of $E$. In particular there is an effective finite search even when the exponents of arbitrarily many prescribed prime bases in the cofactor vary simultaneously.

For the finiteness claim, first bound $s$ for fixed $p,q$. In $w=p^i+p^j=q^k+q^l$, with $i\le j$, $k\le l$, put $H=\max(3,j,l)$. The [lifting argument for two sums](three_unit_layers.md#classification-at-all-54-primes-below-257) gives $p^i,q^k\le C_{p,q}H$, using the multiplicative orders at these two fixed primes. Thus $|p^j-q^l|\le C_{p,q}H$. At large $H$ both exponents are positive and the larger-exponent power is at least $2^H$. The same two-logarithm bound gives $H\log 2\ll_{p,q}1+\log H$, an effective upper bound for $H$. Enumerate these two-power sums, retain the prime-power values of $s=w+1$, and use $D=s$, $E=1$ in (2). This gives a finite exponent search even for moving bases beyond the classified range.

### An effective bound from the smaller moving power alone

Let $M_0=\min(p^a,q^b)$. Even with unrestricted $E$, at a prime-power remainder $s=r^d$ there is the effective bound
$$
          \log(2W)\ll_{p,q,s}\log(2M_0)
                    \big(1+\log\log(3M_0)\big)^{d(r-1)+2}.
$$
Thus a prescribed upper bound on just one of the two moving powers gives an effective finite search over all cofactors and over the entire other exponent.

Interchange the bases so $Q=q^b=M_0$ and $P=p^a$. The layer gives $EQ=1+t p^k$ with $1\le t<Q/s^2$. In the monotone expansion above this yields an equality between two sums in coprime bases,
$$
                \sum_{i=1}^m c_i s^{k_i}=w+w p^a+w t p^{a+k}.
$$
All coefficients are positive and bounded by $T=2wQ$. Partition the equality into minimal equal subsums, repeatedly removing a minimal one from both sides. In each component all proper nonempty signed subsums are nonzero. Let $x_1\ge\cdots\ge x_n$ be its pure powers of $s$ or $p$, before multiplication by their coefficients. Here $n\le m+3$. Write $D_j=\log(x_1/x_j)$ and let $B$ be the maximum of three and all exponents in the original equality.

For $j<n$, the signed prefix has absolute value at most $nT x_{j+1}$ by its omitted tail. If all its signs agree, its absolute value is at least $x_1$. Otherwise group the two bases to write it as $A s^u-C p^v$, where $u,v$ are the smallest respective exponents in the prefix. Factoring those powers gives $1\le A,C\le nT\exp(D_j)$. Matveev applied to $(A/C)s^u p^{-v}-1\ne0$ gives a lower bound for its absolute value of $\exp(-K\log(2B)(\log(2T)+D_j))$, with effective $K$ depending on the fixed bases and number of terms. Omit any trivial factors. If only a rational constant remains, the same bound follows directly from its denominator. Multiplying by $C p^v\ge x_j$ and comparing to the tail gives in both sign cases
$$
              D_{j+1}\ll_{p,q,s}\log(2B)\big(\log(2T)+D_j\big).
$$
Starting at $D_1=0$ bounds $D_n$ by $O_{p,q,s}(\log(2T)(\log(2B))^{n-1})$. At the full equality $A s^u=C p^v$, so $s^u\mid C$ because $s,p$ are coprime. As $C\le nT\exp(D_n)$ and $x_n\le s^u$, we have $\log x_1\le\log(nT)+2D_n$. This applies to every component, including the one containing the largest power in the full sum. The exponents satisfy $B\ll_s\log(2W)$, since each term is at most $w(W+1)$. Conversely $W+1\le s^{k_1}$, so bounding the largest power bounds $W$. Invert the resulting inequality as in (2), now with $n-1\le m+2$ and $T=2wQ$.

## Finiteness when the cofactor and both exponents vary

**Fix distinct primes $p,q$ and impose the two sparse layers. Only finitely many exceptional weights can have a common remainder $s$ that is a prime power.** Here $E$ may vary over arbitrary positive integers coprime to $pq$. In particular this finiteness holds for the entire paired case at bases 2 and 3. Among the 1,431 pairs below 257, only ten pairs can still admit an unbounded sequence of failures under the layer conditions. Their composite remainders are given below.

| $p$ | $(q,s)$ |
|---|---|
| 2 | $(11,133),(19,21),(67,69),(131,133)$ |
| 3 | $(83,85),(89,91),(251,253)$ |
| 11 | $(131,133)$ |
| 13 | $(47,2211),(181,183)$ |

For fixed $p,q$ there are only finitely many possible $s$, by the two-logarithm argument above. Fix one such $s$. Let $P=p^a$, $Q=q^b$ and choose the larger of $P,Q$, calling it $R$ with prime base $v$. Because $W=EPQ$, the layer at $v$ gives
$$
            R L_v\mid W-R,\qquad R L_v\ge R s^2 E\ge s^2\sqrt W.           \tag{3}
$$
In particular $W>R$, as the other prime power is at least two.

We use Evertse's Theorem 2 over the rationals, equation (14), with $T=\{v\}$. For a fixed finite set $S$ containing infinity and $v$, let $x_1,\ldots,x_n$ be nonzero integers whose prime factors belong to $S$. If every nonempty subsum is nonzero, then for each $\epsilon>0$
$$
       \left|\sum x_i\right|_v
          \ge C_{S,n,\epsilon}\max_i |x_i|_v\,(\max_i |x_i|)^{-\epsilon},
          \qquad |v|_v=v^{-1}.
$$
The factors in the product formula for each $S$-unit cancel. [Evertse, *On sums of S-units and linear recurrences*, p. 228](https://www.numdam.org/item/CM_1984__53_2_225_0.pdf#page=5)

For $s=r^d$ use the monotone expansion above to write
$$
                    w(W-R)=\sum_i c_i s^{k_i}-w-wR .
$$
These are integer $S$-units with $S$ consisting of infinity and the primes dividing $pqrw$. Each has size at most $w(W+1)$. Delete any zero subsums until none remain. The total is positive, so some positive term $c_i s^{k_i}$ remains. It is a $v$-adic unit because $v\ne r$. The number of remaining terms has a fixed upper bound, so the minimum of the finitely many constants in Evertse's inequality is positive. At $\epsilon=1/4$ the local norm of the total is $\gg_{p,q,s} W^{-1/4}$. By (3) it is at most $s^{-2}W^{-1/2}$. This bounds $W$ and proves finiteness. The theorem used here gives no numerical cutoff. The effective bound (2) remains available when $E$ is specified.

**Digit consequence for every common remainder, including composite $s$.** Fix also an integer base $t\ge2$ prime to $pq$, integers $u\ge1$ and $z$, and a positive integer $K$. Among the weights satisfying the two-layer hypotheses, only finitely many can fail if the expansion of $Y=uN+z>0$ in base $t$ has at most $K$ maximal runs of equal nonzero digits. Thus along any increasing sequence of failures at a fixed pair, the number of these runs must tend to infinity in every such base, even after a fixed affine change of the weight.

To see this, write $(t-1)Y$ as a sum of at most $2K$ signed terms. A run of digit $d$ in positions $i$ through $j$ contributes $d(t^{j+1}-t^i)$. Subtract $(t-1)(u+z)$ and $(t-1)usR$ to obtain $(t-1)us(W-R)$. Omit zero terms. All prime factors of the terms belong to a fixed finite set, and each term is $O_{t,u,z,s}(W)$. Remove zero subsums again. The only term whose $v$-adic norm could tend to zero is the negative term involving $R$. The total is positive, so at least one of the other terms remains, with $v$-adic norm bounded below in terms of the fixed data. The same estimates at $1/4$ and $1/2$ apply. This argument allows an unbounded total number of nonzero digits within the runs.

Run `python3 check_residual_power_classes.py` and `lean +leanprover/lean4:v4.36.0 ResidualPowerClasses.lean` to reproduce the finite checks.
