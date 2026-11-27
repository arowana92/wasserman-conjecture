# Repunits of twice-prime length

Let $R_n(X)=1+X+\cdots+X^{n-1}$ and $Q=p^a$ for a prime $p$ and $a\ge1$. At **each of the composite lengths $n=10,14,22,26$**, Wasserman's conjecture holds for every three positive trinomial rows at $N=R_n(Q)$, for every $Q$.

There is also an effective large-base theorem. Let $h$ be an odd prime and $m=2h-1$. If
$$ 2\varphi(m)>3h\qquad\text{or}\qquad m=3r\text{ with }r\text{ prime}, $$
then the conjecture holds at $R_{2h}(Q)$ for all sufficiently large $Q$, with the threshold computable from $h$. This includes $m$ prime whenever $h\ge5$.

## Which digit partitions can have a large common factor?

Use the [endpoint reduction and digit allocation](repunit_weights.md#what-a-single-digit-subset-must-satisfy). In a counterexample a row avoiding $p$ partitions the $2h$ powers of $Q$ into three parts $F_i(Q)$. Some $U\ge3$ divides $N$ and every part. Write $ST=N-1$, where primes in $S$ are assigned to this row and primes in $T$ to the other nonendpoint row. We have $27T^2\le U^3$. For every full power assigned to $S$ the three residues of this row are $0,0,1$.

The weight polynomial is $\Phi_h\Phi_{2h}\Phi_2$. A part with a smallest number $k$ of terms satisfies $2\le k\le 2h/3<h$ since $U$ excludes a monomial. Its polynomial is coprime to $\Phi_h$ because $\Phi_h(1)=h$ does not divide $k$. Thus $\gcd(U,\Phi_h(Q))$ is bounded by a nonzero integer resultant depending only on the partition.

Call a partition paired if both positions in each $\{i,i+h\}$ belong to the same part. If all three polynomials are divisible by $\Phi_{2h}$, the partition must be paired. In fact the difference between the indicator coefficients at $i$ and $i+h$, for $0\le i<h$, gives a polynomial of degree at most $h-1$ vanishing at a primitive $2h$th root. It is zero or $\pm\Phi_{2h}$, since its coefficients are $0,\pm1$. In the latter case this part occupies exactly one position of every pair. A second nonempty part with the same divisibility must use the other position in every pair, leaving no third part. This proves the claim.

In an unpaired partition choose a part with nonzero resultant against $\Phi_{2h}$. This bounds $\gcd(U,\Phi_{2h}(Q))$. Hence $U\ll_h Q+1$. Here one may bound a divisor of a product by the product of its gcds with the factors, even when the evaluated factors share primes. For the smallest part there is at least one bad nontrivial $m$th root, meaning $F(\zeta)\notin\{0,1\}$, by the Parseval argument in the [prime-length proof](repunit_weights.md#the-two-resultant-bounds), which uses only $2\le k\le n/3$. Since $m$ is odd, its conjugate is distinct and bad too. The resultant bound gives $T\gg_h (Q-1)^2$. This contradicts $27T^2\le U^3$ at sufficiently large $Q$, at **every** twice-prime length.

## Paired partitions

Here $F_i=(1+X^h)G_i$, where the $G_i$ partition positions $0,\ldots,h-1$. We have $U\ll_h Q^h$ using the bound at $\Phi_h$ above. Also $1+Q^h$ is coprime to $N-1$, so two of the $G_i(Q)$ must vanish at each full power in $S$.

At a divisor $d>1$ of $m$, suppose no two $G_i$ vanish at a primitive $d$th root. For each pair at least one of its two resultants with $\Phi_d$ is nonzero. Consequently $\gcd(S,\Phi_d(Q))$ is bounded uniformly in $Q$. More precisely, write $a_i=|\operatorname{Res}(\Phi_d,G_i)|$, allowing zero. An upper bound divisible by that gcd is
$$ c_d=\operatorname{lcm}\big(\gcd(a_1,a_2),\gcd(a_2,a_3),\gcd(a_3,a_1)\big). $$
Each prime power in the gcd divides two of the norms by the integral resultant identity. The condition on the pair is exactly $c_d\ne0$. Multiply $\Phi_d$ over all such $d$ to form $H$, and set $e=\deg H$. Since $H(Q)\mid N-1$, prime-by-prime comparison gives $H(Q)\le T\prod c_d$, the product taken over those same orders. Therefore $T\gg_h Q^e$.

The order $m$ always contributes. At $\zeta=\exp(2\pi i/m)$ all the powers from $1$ through $\zeta^{h-1}$ lie in an open semicircle since $h-1<m/2$. No nonempty sum with nonnegative coefficients can be zero there. Thus $e\ge\varphi(m)$ and $2\varphi(m)>3h$ suffices.

If $m=3r$ with $r$ prime, the order $r$ contributes as well. Two sums vanishing at a primitive $r$th root would each need a positive multiple of $r$ terms, yet $h=(3r+1)/2<2r$. For $r\ne3$ this gives $e\ge3(r-1)=2h-4$, enough for $h>8$. The only prime $h\le8$ in this case is $h=5$, $r=3$. There $\varphi(9)+\varphi(3)=8$ suffices. This proves the large-base theorem. All choices of supports and their integer resultants are finite and effective.

## Checking every base at the four stated lengths

In both cases the verifier proves $U\le C(Q^z+1)$ or a constant bound, and $H(Q)\le T K$ for a product of cyclotomic polynomials of degree $e$ with $2e>3z$. It certifies $Q^e\le vH(Q)$ at $Q\ge2$ by checking the coefficients after the substitution $Q=Y+2$. These estimates give an integer upper bound for $Q^{2e-3z}$.

For unpaired partitions the program enumerates only the smallest subset. When this subset vanishes at order $2h$ it uses the uniform bound $(h^2/(h-1))^{(h-1)/2}$ for the absolute nonzero norm of another part there. Indeed the signed differences from opposite positions have squared coefficients summing to at most $h$. Substitution $X=-Y$ and Parseval at the $h$th roots bound the sum of the squared absolute values at the primitive $2h$th roots by $h^2$. The other bounds use the actual subset. The program recovers the exact integer norms in a splitting prime field whose characteristic exceeds twice the elementary absolute-value bound.

For paired partitions it enumerates the unordered triples of nonempty supports on $h$ positions. The constant multiplying $Q^h+1$ in the bound on $U$ is the gcd of the three norms at order $h$. For the predecessor it uses the $c_d$ above. At large cutoffs it refines these norm bounds at each candidate prime power. For the weight all three $G_i$ must vanish at the same root. For a predecessor factor two must vanish. If $s\nmid d$, the roots modulo $s^j$ are the primitive $d$th roots obtained by lifting from $s$. If $s\mid d$ in this refinement, $s$ is odd, the roots modulo $s$ have order $d/s^{v_s(d)}$, and $v_s(\Phi_d(Q))\le1$. The last assertion follows by applying the lifting-the-exponent formula to successive $s$-power multiples of the order of $Q$ modulo $s$. For each local contribution the program includes its possible roots and the option of exponent zero, then combines the congruences by CRT. This bounds the base separately in each class.

For each remaining evaluation the gcd of the parts bounds $U$. In an unpaired test the gcd of the weight and the single part suffices. In a paired test the largest possible $S$ divides the gcd of $N-1$ with the lcm of the three pairwise gcds of the $G_i(Q)$. Further necessary conditions come from prime factors up to 500. If a prime produces a carry, the chosen row cannot receive its power on either side. For an unpaired test, a carry when adding $F(Q)$ and $N-F(Q)$ already suffices. Every candidate prime-power base is excluded.

| Length | Unpaired subsets | Largest cutoff | Base and subset tests | Paired partitions | Largest cutoff | Base and partition tests across classes |
|---:|---:|---:|---:|---:|---:|---:|
| 10 | 165 | 4 | 78 | 25 | 0 | 0 |
| 14 | 1,456 | 7 | 1,504 | 301 | 4 | 66 |
| 22 | 280,577 | 57 | 1,226,244 | 28,501 | 344 | 20,217 |
| 26 | 2,533,960 | 2,297 | 16,021,932 | 261,625 | 1,406,713,432 | 227,688 |

A separate standard-library calculation checks a partition attaining the largest cutoff. Its half-supports are $\{0,1,3,4,7\},\{2,10,11\},\{5,6,8,9,12\}$. Exact integer determinants give the three norms $53,53,53$ at order 13 and $63125,4651,63125$ at order 25. Since $63125=5^4\cdot101$ and $25\nmid\Phi_{25}(Q)$, the latter common bound reduces to 505. Directly visiting all residue classes modulo $53\cdot505$ leaves 61,586 integer bases to check. Of these, 10,614 need a carry test beyond the numerical gcd bound. None survives.

~~~sh
python3 check_double_repunit_weights.py
python3 check_double_repunit_extreme.py
~~~
