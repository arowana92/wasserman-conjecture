# A gap in the weight and two prime powers in its predecessor

**Theorem.** Write $N=1+DW$, with $D,W$ positive and $W$ having at most two distinct prime divisors. Let $P=p^{v_p(N)}$ for a prime $p\mid N$. If for some $m\ge1$
$$
N/P\equiv 1\pmod {p^m},\qquad p^m\ge D,
$$
then any three positive trinomial coefficients of weight $N$ share a prime. $D$ and $W$ may share prime factors.

For instance, fix **any odd $D$**, any odd prime $q$, any positive exponent $a$, and any prime $r\equiv5\bmod8$. Set $m=\max(1,\lceil\log_2 D\rceil)$. At least one in every $2^m$ consecutive positive exponents $b$ satisfies the full conjecture at $1+Dq^a r^b$. The binary congruence alone certifies a set of exponents of natural density $2^{1-m}$. Both statements are uniform in $q,a,r$.

## Proof of the theorem

A coefficient avoids a prime if the prime does not divide it. Suppose a counterexample exists. By the all-weight endpoint reduction it contains the row $(1,1,N-2)$. Its two remaining rows $A,C$ must between them avoid all primes dividing $N(N-1)$. Assign each such prime to one row avoiding it. Let $U,V$ be the respective products of the assigned full prime powers of $N$. Then $UV=N$, $\gcd(U,V)=1$, and $A,C$ are multiples of $U,V$ respectively, by the zero low digits of $N$. Their normalized triples sum to $V,U$, respectively. [Endpoint reduction](all_weights_endpoint_reduction.md) · [Two-power assignment](two_prime_powers.md)

Write $B=p^m$. For the row assigned $p$, division of its parts by $P$ leaves a carry-free addition in base $p$, by Kummer's criterion. Its three least nonnegative residues modulo $B$ sum to exactly one. In particular, two parts of that row divided by $P$ are multiples of $B$.

**Both predecessor bases in one row.** Name that row $A$. At each prime power dividing $W$, taken to its full exponent in $W$, the residues of the three parts are $0,0,1$ in some order. This follows from $N\equiv1$ at that power and carry-free addition. The two sets of zero positions intersect, so some part is $kW$, with $U\mid k$ and $0<k<D$. Here a part is at most $N-2<DW$ and $\gcd(U,W)=1$. The same reasoning works for one base or for $W=1$.

If $C$ takes $p$, its parts are $V f_i$ with $f_1+f_2+f_3=U\le k<D$. Since $V/P$ is a unit modulo $B$, two of these positive $f_i$ must be divisible by $B$, which is impossible.

If $A$ takes $p$, write $k=Pz$. The residue of $zW$ cannot be zero modulo $B$ since $0<z<D\le B$ and $p\nmid W$. It must be one. Now $DW=N-1\equiv P-1\bmod B$, so multiplication by $z$ gives
$$
z(P-1)\equiv D\pmod B .
$$
But $0<z(P-1)<Pz=k<D\le B$. The positive difference is strictly less than $B$, a contradiction. [Lean verification of this congruence obstruction](PredecessorWeightGap.lean)

**The predecessor bases in different rows.** Factor $W=QR$ into the powers of its two distinct primes. Suppose $A$ takes the base of $Q$, and $C$ the base of $R$. Two parts of $A$ are $UQx,UQy$, and two of $C$ are $VRs,VRt$, for positive multipliers $x,y,s,t$. The positive remaining parts give
$$
V\ge Q(x+y)+1,\qquad U\ge R(s+t)+1,\qquad (x+y)(s+t)<D.
$$
The last inequality follows by multiplication and $UV=1+DQR$. All four multipliers are less than $D$. If $A$ takes $p$, at least one of $x,y$ must be zero modulo $B$ because $(U/P)Q$ is a unit. If $C$ takes $p$, the same holds for $s,t$. Both possibilities contradict positivity and $D\le B$.

## Combining different weight primes

Here take $Q,R$ to be the full powers of two distinct primes in $N-1$, so $D=(N-1)/(QR)$ is coprime to $W=QR$. Prime-power weights already have a universal prime. At other weights, for every $p\mid N$ put
$$
P_p=p^{v_p(N)},\quad B_p=p^{v_p(N/P_p-1)},\quad G=\prod_{p\mid N}B_p .
$$
When $B_p>1$ let $j_p$ be the least positive inverse of $W$ modulo $B_p$. Put $j_p=1$ when $B_p=1$. Define $F=\{p\mid N:P_pj_p\ge D\}$.

1. In a split assignment, $DG\mid xyst$. Each full power of $D$ divides at least one displayed multiplier in its assigned row, since at least two of its three parts are zero there. Each $B_p$ also divides a displayed multiplier by the single-low-digit argument. All these factors are pairwise coprime. The arithmetic-geometric mean inequality and $(x+y)(s+t)<D$ now require $16G<D$.
2. In a same-row assignment, every prime in $F$ must be assigned to $C$. If assigned to $A$, it would give $k=P_p z$ with $zW\equiv0$ or $1\bmod B_p$. Either residue requires $z\ge j_p$, contradicting $k<D$. The case $B_p=1$ uses only $z\ge1$.
3. For $p\in F$, two of the three $f_i$ are multiples of $B_p$. Thus $2B_p+1\le U<D$. For distinct $p,q\in F$, one common position is a multiple of $B_pB_q$, so $B_p B_q<D$. Also
$$
27\left(\prod_{p\in F} B_p\right)^2\le 27 f_1 f_2 f_3\le U^3\le(D-1)^3.
$$

If $16G\ge D$ and any one of the requirements in step 3 fails, the conjecture follows. For example at $N=1{,}001{,}234{,}650{,}368$ the predecessor is $7\cdot107\cdot1861\cdot718303$. Choose the last two powers, so $D=749$. The gaps at 2 and 17 are 4 and 289, and $P_pj_p$ at those primes is 768 and 1717. Both are forced onto $C$ in a same-row assignment, but $4\cdot289>749$. The same product excludes a split.

## Families of exponents and primes

For the binary family in the introduction put $c=Dq^a$. At every $h\ge2$ the powers of $r\equiv5\bmod8$ run through all residues equal to one modulo four, with period $2^{h-2}$ modulo $2^h$. Indeed $v_2(r^{2^i}-1)=i+2$ for $i\ge0$, by factoring the difference of squares at each step.

If $c\equiv1\bmod4$, choose $v=1$. If $c\equiv3\bmod4$, choose $v=2$. Then $cr^b\equiv2^v-1\bmod 2^{v+m}$ has exactly one solution class modulo $2^{v+m-2}$. At such an exponent $2^v\parallel N$ and $N/2^v\equiv1\bmod 2^m$, proving the claimed bound on gaps. When $c\equiv1\bmod4$, every weight has $v=1$ and the resulting density is $2^{1-m}$. When $c\equiv3\bmod4$, every $v\ge2$ gives one such class, and their densities sum to $\sum_{v\ge2}2^{2-v-m}=2^{1-m}$. The omitted tail $v>H$ has upper density at most $2^{1-H}$ for $H\ge2$, by the same order calculation modulo $2^{H+1}$. This proves the natural density assertion.

For an arbitrary positive $D$, take any primes $p,q$ with $p\nmid Dq$, any $a\ge1$, and $p^m\ge D$ with $m\ge1$. Every prime $r$ in the class
$$
r\equiv(p-1)(Dq^a)^{-1}\pmod{p^{m+1}}
$$
works at $N=1+Dq^a r$, by taking $P=p$. The class contains infinitely many primes by Dirichlet's theorem. [Proof of Dirichlet's theorem, MIT 18.785](https://ocw.mit.edu/courses/18-785-number-theory-i-fall-2021/mit18_785f21_lec18.pdf)

## Reproduction

Run `python3 check_predecessor_weight_gap.py` and `lean PredecessorWeightGap.lean` (Lean 4.36). The exact small check computes binomial products, takes their gcds with the endpoint coefficient, then checks every pair of distinct or repeated gcd values. Through 300 it independently checks 431,408 rows and 4,646 pairs at the 147 weights selected by the criterion. A separate modular check verifies 896 exponent progressions.

On the 100,000 consecutive weights starting at $1{,}001{,}234{,}567{,}890$, the program chooses the two largest full predecessor powers. It certifies 22,997 weights by the one-gap theorem and another 8,256 by combining the criteria above. Of the latter, 8,232 use a single forced prime, 18 use two, and 6 use the product inequality after the preceding tests. Another 7,250 weights have a prime-power weight or fewer than two predecessor bases and are omitted from these counts. The tests give no conclusion for the other 61,497 weights. [Certificate and independent coefficient check](check_predecessor_weight_gap.py)

### All 100,000 weights in the interval

The full conjecture holds at every weight in the interval above. Run `python3 check_weight_gap_assignments.py`. For each assignment of primes to the two rows, impose one modulus for each assigned predecessor power and each assigned $B_p>1$. In a normalized row of sum $H$, each modulus has residues $H,0,0$ in some order. For a predecessor power this follows from the row scale times $H$ being $N$. For a weight gap it follows by dividing the single-low-digit condition by the unit part of the scale.

For each allocation of exceptional positions, the Chinese remainder theorem gives one positive minimum for each part modulo the product $L$ of the row's moduli. The three minima fit precisely when their sum is at most $H$, since their sum is congruent to $H$ modulo $L$. This independently rechecks all 31,253 inequality certificates. Exactly 10,207 of the 100,000 weights admit some pair at this stage.

At least one row has $L\ge H$. Indeed the two products multiply to $(N-1)G\ge N-1$, while the two sums multiply to $N=UV$. If both inequalities went the other way, the products would multiply to at most $(U-1)(V-1)<N-1$. In a row with $L\ge H$ the positive representatives are therefore the actual normalized parts. The verifier multiplies by the row scale and evaluates factorial floors at every assigned prime. No weight in the interval survives.

As a second full-interval check, the earlier subset-sum verifier factors the inputs using FLINT, works with the predecessor moduli alone and enumerates individual normalized parts by subsets rather than assigning CRT positions. It tested 1,590,756 assignments after its size bounds, with zero surviving rows. The following command reproduces that count with python-flint and SymPy installed. [Subset-sum verifier](check_power_partner_assignments.py) · [Weight-gap assignment verifier](check_weight_gap_assignments.py)

~~~sh
python3 - <<'PYTHON'
from check_power_partner_assignments import check_weight

totals = [0, 0]
for n in range(1001234567890, 1001234667890):
    assignments, rows = check_weight(n)
    totals[0] += assignments
    totals[1] += rows
print(totals)
PYTHON
~~~
