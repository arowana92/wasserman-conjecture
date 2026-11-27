# Almost all cofactors work simultaneously at superpolynomial prime powers

There are absolute constants $c,\delta>0$ such that, outside $O(H^{1-\delta})$ integers $n\in[H,2H]$, Wasserman's full conjecture holds at **every**
$$
 1+n p^e,\qquad p\le1249\text{ prime},\quad e\ge1,\quad e\log p\le c\log H\log\log H.
$$
The exceptional cofactors account for all these primes and exponents together. In particular, any fixed power of $H$ is eventually below the upper bound for the scale $p^e$. The cofactor may itself be divisible by $p$.

## Counting with the large scale removed

Fix $B=p^e$. In an exception of weight $N=Bn+1$, assign $p$ to a row $A$ avoiding it. The [endpoint assignment](small_weight_transfer.md#specializing-the-two-rows) gives $A=(Ba_1,Ba_2,Ba_3+1)$ with positive slopes summing to $n$, and a second row $C=V(y_1,y_2,y_3)$ with $y_1+y_2+y_3=U$, $UV=N$, $U,V\ge3$. The slopes $a_i$ add without a base-$p$ carry, by reading above the bottom $e$ places of $A$.

The predecessor product assigned to $A$ is $Bs$ for some $s\mid n$, and the remaining product is $T=n/s$. Put $K=Us$. The following bound for $1\le Y\le H$ is uniform in $B$:
$$
 \#\{n\in[H,2H]\text{ admitting an assignment with }K\le Y\}
 \ \ll\ H^{-2/3}Y^2 (\log(2H))^{O(1)}
 \max_{1\le j\le B(2H)^{10}}\tau(j). \tag{1}
$$
Only the normalized complementary row and the scale congruence are needed here. To see this, group the prime powers of $T$ according to which part of $C$ has residue 1. The resulting three pairwise coprime products include a smallest one $d\le T^{1/3}$; empty products are allowed. Write $T=dt$. The designated normalized part is $zt$, and its two companions sum to $dy$, with $z,y>0$. Thus
$$ U=zt+dy,\qquad UV=1+Bsdt,\qquad Bs d^2 y-z=(Bs d-Vz)U>0. $$
Positivity follows from $(Bs d-Vz)t=Vdy-1>0$. The [Lean proof](ScaledCofactorBins.lean) verifies the identity and positivity in natural numbers. Once $s,d,z,y$ are given, the positive divisor $U$ determines $t$ and $n$. This is the [one-bin divisor argument](common_factor_count.md#counting-small-assigned-moduli-using-one-predecessor-factor) with the known factor $B$ removed from the count of choices for $s$.

Put $U,s,d$ in dyadic ranges and use their size symbols on the right up to absolute constants. Since $T\asymp H/s$, $z\le Ud/T$ and $y\le U/d$, the number of parameter choices in a nonempty box is at most
$$
 O\big(s d (Ud/T)(U/d)\big)\ll H^{-2/3}U^2 s^{5/3}\le H^{-2/3}K^2
$$
up to the fixed dyadic factors. The bounds count positive integers even when the allowed interval is short. All candidate dividends in an actual assignment lie within the range in (1). Summing the boxes proves (1).

For sufficiently small fixed $c>0$ and $\log B\le c\log H\log\log H$, the divisor-function factor in (1) is $H^{O(c)+o(1)}$ with an absolute implied constant. This follows from $\log\tau(m)\ll\log(2m)/\log\log(3m)$. One elementary proof splits primes at $\sqrt{\log m}$: the small primes contribute $O(\sqrt{\log m}\log\log m)$ to $\log\tau(m)$ and, above the split, the sum of exponents is at most $2\log m/\log\log m$.

## Reading digits at the size of the cofactor

Choose $k<5/6$ close to $5/6$. By (1), assignments with $K<H^k$ contribute at most $H^{2k-2/3+O(c)+o(1)}$. At greater $K$ a specified pair $P=(U,s)$ leaves one residue for $n$ modulo $Us$, since $s\mid n$, $Bn\equiv-1\pmod U$, and $\gcd(U,Bs)=1$. Each individual slope $a_i$ permits at most $H^{O(c)+o(1)}$ pairs $P$: its actual part $z=Ba_i+\mathbf1_{i=3}$ satisfies $U\mid z$ and $Bs\mid z(z-1)$, and $z\ge3$. Apply the same divisor bound to $z$ and $z(z-1)$.

Let $l=\lceil\log_p(2H+1)\rceil$. Outside $O_p(H^{1-\gamma_p})$ weights, the sum of the $l$ base-$p$ digits of $n$ is at most $((p-1)/2+\epsilon)l$ for any fixed $\epsilon>0$ and some $\gamma_p>0$. This follows from the bounded independent-digit tail on the complete range of $l$-digit strings, which has size $O_p(H)$. Choose a canonical witness for each remaining exceptional cofactor and give the cofactors the uniform distribution. Use natural-log entropy $\mathcal H$. As the slopes determine $n$ and each slope almost determines $P$, the [conditional entropy inequality](common_factor_count.md#information-that-every-individual-part-supplies) yields
$$
 \mathcal H(n)\le \tfrac13\sum_i\mathcal H(a_i)
 +\tfrac23 (1-k)\log H+(O(c)+o(1))\log H .
$$

For $0<x<1$, a digit distribution of mean $m$ has entropy at most $\log(1+x+\cdots+x^{p-1})-m\log x$, by comparison with the distribution proportional to $x^j$. The total of all digit sums of the three slopes equals the digit sum of $n$. Thus with $x=(p-1)/(p+5)$ the first term is at most $(A_p+O_p(\epsilon)+o(1))\log H$, where
$$ A_p=\log_p\!\left((1+x+\cdots+x^{p-1})x^{-(p-1)/6}\right). $$
For every prime through 1249 the verifier checks $A_p<8/9$ by raising the rational inequality to the eighteenth power. At $k=5/6$ the extra term $2(1-k)/3$ is $1/9$. The strict digit margins allow a common $k<5/6$, then sufficiently small positive $\epsilon,c$, so both the small-$K$ and large-$K$ exponents are strictly below one. The sets of primes and hence of relevant margins are finite. The $O(\log H\log\log H)$ scale exponents per prime cost only $H^{o(1)}$ in the final union. This proves the theorem.

Run `python3 check_growing_steps.py` to reproduce the integer inequalities at all 204 primes, 31,920 independent factorial-valuation comparisons, 568 complementary row and scale tests from exact coefficients, and 250 constructed reconstructions with three nonempty bins. The exponent bound is an asymptotic theorem. The finite tests check the identities and strict inequalities used in its proof.
