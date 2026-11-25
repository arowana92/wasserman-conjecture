# A large prime power dividing the predecessor of the weight

Let $Q=q^e$ be a full prime power in $N-1=dQ$, so $e\ge1$ and $q\nmid d$. Write $E=(1,1,N-2)$, and, at odd weights, $M=(1,(N-1)/2,(N-1)/2)$. Rows may be permuted throughout.

**Theorem 1.** For all sufficiently large $N$ with $d\le N^{3/50}$, a counterexample must contain $E$, or must have even $d$ and contain $M$. The threshold here is absolute. In particular the full conjecture holds for all sufficiently large such weights when $d$ is an odd prime power or 1. Indeed $N-1$ then has at most two prime divisors, which [excludes the endpoint](endpoint_reduction.md).

**Theorem 2 (finite automata, unbounded weights).** For every $d\le1000$, the endpoint cannot occur in a counterexample at $N=dQ+1$ when $Q$ is a power of two. For odd $Q$ the same holds for every $d\le1000$ except possibly
$$364,380,589,629,646,703,710,780,850,874,880,946,950,962,986.$$
In particular, for **each odd $d<589$**, the full conjecture holds at every sufficiently large $dQ+1$ with $Q$ a prime power coprime to $d$. It also holds for every sufficiently large $d2^e+1$ with odd $d\le1000$. For odd $d\le1000$ the full conjecture also holds for all sufficiently large such weights whenever $Q$ is a square or cube, or whenever its prime base is below 2,897. The endpoint assertion itself has no size threshold. The listed cofactors pass necessary conditions. For every odd $d\le1000$, including 589, 629 and 703, the endpoint is also excluded throughout $N<7\cdot2^{64}$.

## Proof of Theorem 1

Suppose three coefficients have gcd one. Some row $A$ avoids $q$. Kummer gives actual parts $A_i=Qx_i+\delta_i$ where the nonnegative $x_i$ sum to $d$ and the $\delta_i$ are $0,0,1$ in some order. At least two $x_i$ are positive. Put $d_0=\gcd(x_1,x_2,x_3)$, $D=d/d_0$, and $z_i=x_i/d_0$. Then $D\ge2$, $\gcd(z_1,z_2,z_3)=1$ and
$$D A_i=z_i(N-1)+D\delta_i.$$
If $D=2$, the normalized parts are $0,1,1$, the zero position must receive the extra 1, and this is exactly $M$ with $2\mid d$. Henceforth suppose $D\ge3$ and the selection omits $E$.

Use height 12 with offsets $J=\{0\le h<12:D\nmid h-1\}$. For each $h\in J$ fix $i$ with $D\nmid z_i(h-1)$, which exists by coprimality. Define
$$K_h=\operatorname{lcm}\big(1,\ldots,12,\ |z_i(h-1)+D(\delta_i-j)|\ (0\le j\le h)\big).$$
These integers inside the absolute values are nonzero. If $w=p^k>12$ divides $N-h$ and $A$ avoids $p$, the residues of $A$ sum to $h$. Taking $j=A_i\bmod w$ in the defining identity shows that $w\mid K_h$. Thus the layers above $v_p(K_h)$ must be covered by one of the two other rows.

For either remaining row with parts $y_i$, choose cutoffs $0\le t_i\le y_i$ summing to 12. Let $P$ be the product over both rows of $B_{t_i}(y_i)$, where $B_t(y)=\prod_{j=0}^{t-1}(y-j)^{t-j}$. The [weighted divisibility lemma](weighted_obstruction.md) gives
$$\prod_{h\in J}\left(\frac{N-h}{\gcd(N-h,K_h)}\right)^{12-h}\mid P.$$
Indeed at a remaining layer the residues in a covering row sum to $h$, supplying at least $12-h$ copies at that layer. With $L=\sum_{h\in J}(12-h)$, we get $N^L\ll P\prod_{h\in J}K_h^{12-h}$.

Write $g=N-r$ for the distance to the largest prime at most $N$. The Baker–Harman–Pintz prime-gap theorem gives $g=O(N^{0.525})$. For large $N$ the row avoiding $r$ has its two smallest parts at most $g$. It cannot be $A$, which has at least two parts at least $Q=(N-1)/d$. Call this the designated row and the other remaining row ordinary. Bergman’s Proposition 10 gives a sum of at least 11 for the two smallest parts of each row since $E$ is absent. [Baker, Harman and Pintz](https://doi.org/10.1112/plms/83.3.532) · [Bergman](https://arxiv.org/abs/0806.0607)

At most one row can contain a 1 or 2. To see this, assign the full prime powers of $N$ to avoiding rows. The product assigned to each row divides all its parts. The products assigned to two rows each having a 1 or 2 would be coprime divisors of 2, of combined product at most 2. The remaining row, divided by its assigned product, would have three positive parts of sum at most 2.

Use the [height-12 cutoff table](binary_weights.md). The ordinary cost is $O(N^{30})$ for minimum at least 3 and $O(N^{36})$ otherwise. The designated cost is at most $O(N^{18.75})$ for minimum at least 3, using $g^{30}N^3$, $g^{15}N^6$ or $g^{21}N^6$. For minimum 1 use $(1,7,4)$, costing $O(g^{28}N^{10})$, and for 2 use $(2,7,3)$, costing $O(g^{28}N^6)$. In either case the ordinary minimum is at least 3. Thus always $P\ll N^{54.75}$.

- For $4\le D\le\operatorname{lcm}(1,\ldots,11)$, $L\ge78-(11+7+3)=57$. All $K_h$ are absolutely bounded, giving a contradiction for large $N$.
- For $D=3$, $L=52$ and $A$ is primitive. Indeed $(z_i)$ is either $(0,1,2)$, giving a part of 1, or $(1,1,1)$, giving consecutive parts. The other rows are imprimitive, since two primitive rows would force the third to have all parts divisible by $N$. Both therefore have minimum at least 2. With ordinary minimum 2 use $(2,5,5)$ at cost $O(N^{30})$. The maximum combined exponent is $30+6+28(0.525)=50.7<52$, and the $K_h$ are bounded.
- For larger $D$, $J$ omits only 1, giving $L=67$. The following refinement bounds the loss from layers covered by $A$ by $O(D^{202})$. It gives $N^{67}\ll N^{54.75}D^{202}\le N^{66.87}$, a contradiction.

For the refinement, choose a single $i$ for which the denominator of $z_i/D$ in lowest terms exceeds 11. Otherwise $D$ divides $\operatorname{lcm}(1,\ldots,11)$ since the $z_i$ have gcd one. For $h\ne1$, $0\le j\le h$, put $m/n=(j-\delta_i)/(h-1)$ in lowest terms with $n>0$. Then $n\le11$ and $n z_i-mD$ is nonzero. The original linear form equals $(h-1)(n z_i-mD)/n$, whose integer multiplier has absolute value at most 11.

For each distinct fraction $m/n$ use the single factor $|n z_i-mD|$ with exponent equal to the maximum $12-h$ over its occurrences. This product covers all large layers requiring $A$, up to an absolute factor. In detail, put $M_0=\operatorname{lcm}(1,\ldots,12)$. For a prime $p$ with $p^k>12$ dividing some $N-h$, all layers above $v_p(M_0)$ occur at just one offset. If $A$ avoids $p$, use its actual residue $j=A_i\bmod p^k$ at the largest such layer. Then $p^{k-v_p(M_0)}$ divides the corresponding $n z_i-mD$. Its exponent in the product is at least $12-h$. For the layers where $A$ does not avoid $p$, use $P$ as before. Layers at most 12 cost an absolute factor. This shows $N^{67}\ll P$ times the indicated product.

For each $\delta_i=0,1$, there are 43 distinct fractions and their maximum exponents sum to 202. Each factor is $O(D)$. These integer counts are reproduced by the verifier.

## Proof of Theorem 2

In a counterexample containing $E$, assign every prime of $N(N-1)$ to an avoiding row other than $E$. Let $A$ receive $q$. Write $UV=N$ for the products of full powers of $N$ assigned to $A$ and the remaining row $C$ respectively. Then $\gcd(U,V)=1$, and $U,V\ge3$ because the two rows, after division by their assigned scales, have three positive parts. Also $\gcd(U,dQ)=1$. The actual rows have the form
$$A=(aQ,bQ,cQ+1),\qquad C=(Vf,Vg,Vh),\qquad a+b+c=d,\quad f+g+h=U.$$
Here $a,b,c\ge1$ because even the shifted part is a multiple of $U\ge3$. We may assume $a\le b$ and $f\le g\le h$. Both $a,b$ are multiples of $U$.

For each full prime power $w\mid d$, a row avoiding its base has residues $0,0,1$ in some order modulo $w$. For $A$ this requires both $a,b$ divisible by $w$, or one divisible by $w$ and the other a unit. All other such $w$ must be handled by $C$, so each divides two of $f,g,h$. If $T$ is their product then $27T^2\le U^3$ and $2w+1\le U$ for each. These are the filters in the finite enumeration of $(U,a,b,c,f,g,h)$.

For a retained shape, $Q=Q_0+Ut$, with $Q_0$ the least positive solution to $dQ_0\equiv-1\pmod U$ and $t\ge0$. Both rows are triples of affine functions of $t$. Test the primes dividing $2dU$. A prime of $U$ requires $A$ to avoid it. A prime of $d$ may be avoided by either row. For odd $Q$ and odd $d,U$, the factor 2 of $N$ is assigned to $C$. For $Q=2^e$, avoidance of 2 by $A$ requires the binary digits of $a,b,c$ to be disjoint. The first digit of $t$ in base $p$ enforces $p\nmid Q$ when required and $p\nmid V$ when $p\mid U$.

For affine parts $k_i t+b_i$, start the carry state at $(b_i)$. Read digits $x$ of $t$ from least significant to most significant. The three output digits are $(k_i x+s_i)\bmod p$. Discard a transition when their sum reaches $p$ and otherwise pass to state $(\lfloor(k_i x+s_i)/p\rfloor)_i$. An input may end when the remaining state adds without a carry. Each coordinate stays within $[0,\max(k_i,b_i)]$, so breadth-first search decides existence for inputs of every length. Different primes may have different witnesses. Failure at one required base excludes the shape.

One further constraint couples two bases. If $C$ lacks two positions divisible by an odd $w\mid d$ and exactly one of $a,b$ is divisible by $w$, we must have $(a+b)Q\equiv1\pmod w$. Intersect this progression with $Q\equiv Q_0\pmod U$ before testing binary carry freedom. It is sufficient to impose one such $w$ at a time.

The enumeration tests 41,029 normalized row pairs. The standard-library verifier prints the surviving cofactors for either parity of $Q$. Run it from the repository root.

```sh
python3 check_prime_power_cofactors.py 1000
```

## The three remaining odd cofactors through $7\cdot2^{64}$

For each surviving odd $d$, the first verifier finds a single shape, always with $U=7$ and $C=(V,2V,4V)$, with the following multipliers for $A$.

| $d$ | $(a,b,c)$ | Two prime factors of $d$ |
|---|---|---|
| 589 | $(133,434,22)$ | 19, 31 |
| 629 | $(259,357,13)$ | 17, 37 |
| 703 | $(133,518,52)$ | 19, 37 |

Here $Q$ is odd. The row $C$ must avoid 2, so the ones of $V$ are at least three binary positions apart. Both primes of $d$ must be avoided by $A$ since none divides two of $1,2,4$. For each prime one of $a,b$ is a multiple and the other is a unit. Therefore $(a+b)Q\equiv1\pmod d$. With $7V=dQ+1$ this fixes $V$ modulo $d^2$.

Split the 64 binary digits of $V$ into halves and enumerate half-strings with ones separated by at least three positions. Match the halves by residue modulo $d^2$, check separation across the split, and require $V$ coprime to 7 and $Q>1$ odd and coprime to $d$. This leaves 93,187, 81,511 and 65,702 candidates respectively. For every one, Legendre’s factorial formula gives a positive valuation of the coefficient of $A$ at 7 or one of the two primes of $d$. The shape requires $A$ to avoid all three. This proves the claimed finite exclusion. The check even allows composite $Q$ satisfying the indicated conditions.

```sh
python3 check_cofactor_exceptions.py
```

## Unbounded exponents at the remaining three cofactors

Write $Q=q^e$. The congruences $dQ\equiv-1\pmod7$ and $(a+b)Q\equiv1\pmod d$ determine a unit $r$ modulo $7d$. For 589, $r\equiv6\pmod7$ is a nonsquare and $r\equiv6\pmod{19}$ is a noncube. For 629 and 703 the residues modulo 37 are 17 and 32 respectively, both of multiplicative order 36. Hence $e$ is coprime to 6 in all cases. This already excludes the endpoint for all square or cube $Q$.

The verifier also enumerates all units $v\bmod7d$ for which some positive power equals $r$. There are 24 of 3,240 unit classes at $d=589$, 12 of 3,456 at 629 and 12 of 3,888 at 703. It suffices to try exponents up to $\operatorname{lcm}(6,p_1-1,p_2-1)$, where $p_1,p_2$ are the primes dividing $d$. At the base $q$ itself, the rows $(aQ,bQ,cQ+1)$ and $(a,b,c)$ have exactly the same carry condition, by shifting the digits $e$ positions. This removes $(d,q)=(589,503)$, since $133+434+22$ carries at 503. Among prime bases $q<2897$, the only surviving pairs $(d,q)$ are $(589,937),(589,1973),(629,2129)$. In all three, $e$ must be odd. Checking $e=1,3\pmod4$ shows $dq^e+1\pmod{16}$ differs from every $7V\pmod{16}$ with the required binary separation. Thus the endpoint cannot occur for these bases at any exponent.

By the prime number theorem for arithmetic progressions, the upper relative density among primes of bases permitting a remaining endpoint at any odd $d\le1000$ is at most $1/135+1/288+1/324<0.014$. For every prime base outside those residue classes, Theorem 1 therefore gives the full conjecture at every sufficiently large $dq^e+1$ with odd $d\le1000$ coprime to $q$.
