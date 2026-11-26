# Five weights with the same two prime-power factors

Let $t=p^i q^j$ with $p,q$ prime and $i,j\ge1$, allowing $p=q$. Call a weight exceptional if three of its positive trinomial coefficients have gcd one.

* At least three of $5t+1,7t+1,8t+1,9t+1$ satisfy the full conjecture. At least one of $8t+1,16t+1$ satisfies it.
* For $t\le10^{1000}$ at least **four** of the five weights satisfy it. This range was exhausted by two independent integer algorithms.
* For unrestricted $t$ the only possible simultaneous exceptions are the pairs of multipliers $(7,16)$ and $(9,16)$. The first pair has only finitely many possible inputs, with an effective upper bound. The second requires $t\equiv1267\pmod{1584}$ and identical parts in its two distinguished rows. In every interval of $H\ge2$ consecutive inputs, only $O_\epsilon(H^\epsilon)$ inputs can have two exceptions, for every $\epsilon>0$.
* The number of possible exceptions at $8t+1$ alone is $O(H^{0.198})$ in every interval of $H\ge2$ consecutive inputs.

## Which pairs can occur

The [three-weight theorem](factor_nine_pairs.md) excludes a pair among the multipliers 5, 7 and 9. For even $t$, both $8t$ and $16t$ have at most two prime factors and the [two-factor theorem](all_weights_endpoint_reduction.md#consequences) applies. For odd $t$, $5t+1$ satisfies the conjecture. If multiplying $t$ by any multiplier under consideration gives fewer than three distinct prime factors, that weight is covered as well.

For $w=8,16$ a possible exceptional assignment gives a row with part $at$, with $3\le U\mid a<w$, $at\equiv1\pmod w$, $U\mid wt+1$. Here $U$ is the product of the weight powers assigned to that row. The row avoids the predecessor primes and all primes dividing $U$. In particular its cut $at,(w-a)t+1$ is binary carry-free. At $w=8$, only $a=U=3,5$ remain, since $a=7$ is [excluded](three_factors_small_power.md). At $w=16$ we allow every odd $a$ and every odd divisor $U\ge3$ of $a$. [Three-factor assignment](endpoint_factors.md#when-n-1-has-three-distinct-prime-factors)

The odd-input restrictions on a simultaneous multiplier $v$ are as follows. In each case $bt\equiv1\pmod v$.

| $v$ | $b$ | Necessary conditions |
|---:|---:|---|
| 7 | 6 | $t\equiv11\pmod{72}$, and $6t,t+1$ add without carries at 2, 3 and 7 |
| 9 | 4 | $t\equiv115\pmod{144}$, and $4t,5t+1$ add without carries at 2 and 3 |
| 9 | 7 | $V=(9t+1)/7$ is integral and $V,2V,4V$ add without binary carries |
| 9 | 8 | $4\mid9t+1$ and $8t,t+1$ add without binary carries |

At $b=7$ the companion avoids 2. Its normalized parts have distinct 2-adic valuations and sum to 7, forcing $1,2,4$. The other cases follow from the [factor-seven restriction](predecessor_five_seven.md#a-restriction-on-every-base-seven-digit) and the [factor-nine assignment](factor_nine_pairs.md#the-allowed-rows-at-factor-nine).

The finite binary test below leaves exactly two assignments, written $(w,a,U;v,b)$:
$$
 (16,3,3;7,6),\qquad (16,11,11;9,4).
$$
Here are details for reproducing the all-length exclusion. For each binary cut read the input digits from low to high. At a digit $d$ and carries $c_i$ for affine parts of slopes $s_i$, retain the edge exactly when the bits of $c_i+s_i d$ sum to at most one. The next carries are the integer quotients by two. The initial carries are the offsets, and stopping is allowed precisely when the remaining carries have disjoint bits. Run both necessary cuts on the same input. In the fractional companion case use $t=7z+3$ and $V=9z+4$. The low $\log_2 w$ digits of the input are determined by $at\equiv1\pmod w$. Reject incompatible low classes in the table. For the remaining odd congruences use $M=\operatorname{lcm}(U,v,9)$ at $v=7$, or $M=\operatorname{lcm}(U,v)$ at $v=9$. Starting from each stopping state with suffix residue zero, close backwards along edges, sending a residue $h$ to $2h+d\bmod M$. This computes exactly the attainable residues of finite suffixes after the prescribed prefix. Test the congruences on the resulting $t$.

For multipliers 8 and 16 together, the multiplier of the distinguished part at 16 is $a$ or $a+8$, with $a=3,5$ at 8. The first case would give $a\mid 8t+1,16t+1$, whose gcd is one. In each of the lifted cases the two binary cuts and $a\mid8t+1$ alone have no accepted suffix. A second verifier traverses the states **forward**, also tracking the input residue and its current binary place modulo $M$. It works directly with rational forms $(st+b)/d$ of odd denominator, including $d=7$. For carry $c$, input bit $x$, and output bit $e=(sx+c)\bmod2$, its next carry is $(sx+c-de)/2$. The terminal high parts must be nonnegative integers with disjoint bits within each row. These carries are bounded because each step halves the preceding carry up to a bounded additive term. Both verifiers recover the same two assignments.

## Multipliers 7 and 16

The binary cuts are $6t,t+1$ and $3t,13t+1$, with $t\equiv11\pmod{16}$. After the four low bits, the complete accepted language consists of $m$ repetitions of $(1,0,1)$, one zero, $k$ repetitions of $(1,0)$ and zero padding, read least significant bit first. Finite words can be extended by zeros without changing their value. The primary verifier gives and checks the full nine-state transition table reached after the initial bits. Consequently
$$
 t=11+80(8^m-1)/7+32\cdot8^m(4^k-1)/3,\qquad m,k\ge0 .
$$
The identity $42(t+1)=448\cdot8^m4^k+32\cdot8^m+24$ gives at most six binary ones. The factor-seven cut forces the digits of $t$ to decrease weakly from low to high in base 7, starting at 6. Thus $6(t+1)$ is a sum of six powers of 7, as is $42(t+1)$. Its digit sum in base 7 is at most six. The effective two-base digit-sum theorem proves finiteness. More explicitly, the [six-unit logarithm bound](seven_binary_mixed_range.md#an-absolute-upper-cutoff-from-the-two-bases) applied to $n=42(t+1)$ at bases 2 and 7 gives $n<2^{10^{180}+1}$. That argument only uses the two expansions and the bound on the bases. [Stewart, *On the representation of an integer in two different bases*](https://uwaterloo.ca/pure-mathematics/sites/default/files/uploads/documents/j-reine-ange-math-1980.pdf)

There is also a restriction **$2k>3m+2$** coming from the actual trinomial rows. Write $t=PQ$ for its two distinct full prime powers. Pair the two nondistinguished parts at 7 with those at 16 by their residues modulo $P,Q$. Each has residue one modulo one power and zero modulo the other, so paired parts are congruent modulo $t$. Write $x,y$ for the parts at 7, with $x+y=t+1$. Both are strictly less than $t$ because both are positive multiples of 42. The other two at 16 are then $x+it,y+jt$ with $i,j\ge0$, $i+j=12$. These parts are multiples of 48, whereas $42\mid x,y$. Since $t$ is coprime to 6, both shifts are divisible by 6. At $i=j=6$, both sums would contain every bit of $6t$, since $6t,x,y$ are binary-disjoint. They would overlap. Interchanging the parts if necessary, the two rows therefore share $x$, with remaining parts $y$ and $y+12t$. In particular $336\mid x$.

Put $A=8^m$, $B=4^k$, $L=32A$. The binary ones of $x$ must occur in both $t+1$ and $13t+1$. Their residues modulo $L$ are $(80A+4)/7$ and $(144A-32)/7$. These are multiples of 4 summing to $L-4$. On division by 4 their sum is an all-ones binary integer, so the two residues have disjoint bits. Hence $x=LD$ where $1\le D\le(B-1)/3$, the latter bound following from $\lfloor(t+1)/L\rfloor=(B-1)/3$. Relabel $P,Q$ so $P\mid x$ and $Q\mid x-1$. Then $P\mid D$, and $42t=L(14B+1)-18$ gives $Q\mid18D-(14B+1)$. Therefore
$$
 t\le D(14B+1-18D)\le (B-1)(8B+7)/3 .
$$
The quadratic is increasing up to $(B-1)/3$. Substitution into the identity for $42t$ gives $L(14B+1)\le112B^2-14B-80<8B(14B+1)$, hence $32\cdot8^m<8\cdot4^k$ as claimed.

## Multipliers 9 and 16

The congruences $t\equiv115\pmod{144}$ and $11\mid16t+1$ give $t\equiv1267\pmod{1584}$. The two distinguished rows must be
$$
 (4t,x,y),\qquad (11t,x,y),\qquad x+y=5t+1,\quad 1584\mid x,y.
$$
Indeed, pair the other parts by the same CRT residues modulo $P,Q$. They sum to $5t+1$ in both rows. Paired parts differ by $it$ with $|i|<5$, because each part is at least 4. The two scales make the parts multiples of 36 and 176, respectively, so $4\mid i$. If a difference of $4t$ occurs, its smaller part is binary-disjoint from $4t$: every part under comparison is a submask of $5t+1$, which is disjoint from $4t$ by the cut at 9. Adding $4t$ to the smaller part therefore introduces its bits into the larger part, which is also a submask of $5t+1$. This contradiction forces $i=0$, giving the divisibility by $\operatorname{lcm}(36,176)=1584$.

Use the binary graph for the two cuts together, the ternary graph for $4t,5t+1$, and the base-11 graph for $11t,5t+1$. The total path counts from all states for blocks of 256 digits have bit lengths 78, 225 and 679. The exact integer inequalities give respective dimension bounds $0.31,0.56,0.77$ for the block digit sets. Their sum is strictly below 2. The [moving-interval argument](five_seven_pairs.md#from-the-digit-sets-to-moving-intervals) applies to these three multiplicatively independent bases, proving the claimed $O_\epsilon(H^\epsilon)$ bound. For factor 8 alone, use the cut $at,(8-a)t+1$ in bases 2 and $a$. At block length 4096 its two dimension bounds are strictly below $(0.488,0.644)$ for $a=3$ and $(0.488,0.710)$ for $a=5$. The strict margins permit the exponent $0.198$ in the two-base covering bound. [Corso and Shmerkin, Corollary 1.17](https://arxiv.org/pdf/2409.04608)

## Independent finite checks

The first search intersects accepted high-digit prefixes in the required bases. It visits 103,498 nodes for the 7/16 pair and 100,600 for 9/16 through $t=10^{1000}$. An independent search bisects ordinary integer intervals and uses factorial floors. For a cut $at,(w-a)t+1$ at a prime-power layer $q$, it excludes all $t\in[l,h]$ when
$$
 \left\lfloor(wl+1)/q\right\rfloor>
 \left\lfloor ah/q\right\rfloor+\left\lfloor((w-a)h+1)/q\right\rfloor .
$$
Each individual floor difference is nonnegative at an actual input, so this lower bound certifies a carry. At singleton intervals the test uses every layer. The search restricts to residues 587 modulo 1008 and 1267 modulo 1584, respectively. It visits 72,263 and 50,379 nonempty intervals and returns no survivors. Small positive controls compare the floor search with exact binomial coefficients. Neither finite algorithm needs to factor $t$.

~~~~sh
python3 check_five_related_weights.py
python3 check_five_related_weights_cross.py
~~~~
