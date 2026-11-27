# Two digit blocks that cannot support the same ratio

Let $p$ be prime, $B,c,h,a$ positive integers, $P=p^h>c$, $D=BP+c$, and $E=Bc$. Say that $r$ occurs in a block $b$ if there are positive integers $i,j$ with $r=i/j$ such that $i,j,b-i-j$ add without a carry in base $p$. A zero third summand is allowed. All digit tests below are in base $p$.

**Theorem.** Suppose no positive rational occurs in both $B$ and $c$. If $27D^2>E^7$, the full trinomial conjecture holds at
$$N=1+p^a(Bp^h+c)$$
for **every** $a\ge1$. The factors of $D$ are unrestricted. If also $p\nmid c$, the preceding weight $N_0=p^aD$ satisfies the conjecture whenever $27(N_0-1)^2>p^{3a} E^6$.

A convenient sufficient test uses the distances between positions of two digit units of each block. Include distance zero if a digit is at least two. If the two sets of distances are disjoint, no ratio can occur in both blocks. Indeed, in a pair realizing $i/j$ the distance between the lowest occupied positions of $i,j$ is $|v_p(i/j)|$. In particular $c$ can be **any** positive integer below $p^k$ if $B$ has digits zero or one with its occupied positions separated by at least $k$.

For example take $B=1+p^k+\cdots+p^{(s-1)k}$, $k,s\ge1$, and any $1\le c<p^k$. Every $a$ works whenever $2h>(5s+2)k$. Here $B<2p^{(s-1)k}$ implies $B^5c^7<32p^{(5s+2)k}<27p^{2h}$. Both blocks can contain arbitrarily many digit units. The finite theorem below removes the gap condition for 1,674 binary families.

## Three determinants

In a hypothetical counterexample reduce to the endpoint and assign each prime of $N(N-1)$ to a nonendpoint row avoiding it. Take the row $A$ assigned $p$. Its three parts are $p^a x,p^a y,p^a z+1$, where $x,y,z$ are positive and add without carries to $D$. Their positivity follows because the assigned factor $U\ge3$ of $N$ divides each actual part. Write
$$x=iP+j,\quad y=kP+l,\quad z=(B-i-k)P+(c-j-l).$$
The three upper coefficients and three lower coefficients distribute their respective blocks without carries. Put $S$ for the factor of $D$ assigned to $A$, counting any extra $p$-power in $D$ in $S$. The other nonendpoint row, after division by its scale $N/U$, has three positive parts of sum $U$. Each full power in $D/S$ divides two of these normalized parts. Each full power in $S$ divides $x$ or $y$, since at most one actual part is nonzero modulo an assigned predecessor power. Extra $p$-powers in $D$ divide both by Kummer's criterion. Consequently
$$U\mid\gcd(x,y),\qquad S\mid\gcd(D,xy),\qquad 27D^2\le U^3S^2. \tag{1}$$
The last inequality uses the arithmetic–geometric mean bound on the product of the three normalized companion parts. Also $U$ is coprime to $pD$. [All-weight reduction](all_weights_endpoint_reduction.md) · [Assignment inequalities](two_block_spacing.md)

The two determinants $Bj-ci$ and $Bl-ck$ are nonzero. For example, if $Bj=ci$, the positive proper part $x$ is proportional to $D$ in both blocks. All four numbers $i,j,B-i,c-j$ are then positive. The ratio $i/(B-i)=j/(c-j)$ occurs in both blocks, a contradiction. It follows that $\gcd(D,x)$ and $\gcd(D,y)$ are at most $E$, since each divides the absolute value of the corresponding determinant. Thus $S\le E^2$.

The remaining determinant $il-kj$ is divisible by $U$. Its absolute value is at most $E$. If it vanishes with $i,j,k,l>0$, the ratio $i/k=j/l$ again occurs in both blocks. In the other vanishing cases, both $x,y$ are constants or both are pure multiples of $P$ with upper coefficients at most $B$. As $p\nmid U$, either case gives $U\le\max(B,c)\le E$. Now (1) gives $27D^2\le E^7$.

At $N_0=p^aD$ with $p\nmid c$, all three parts of its $p$-avoiding row are $p^a x_i$ with the $x_i$ splitting the same $D$. For each pair the argument for $il-kj$ bounds the part of $\gcd(x_i,x_j)$ coprime to $p$ by $E$. Write $p^a U$ for the assigned weight scale and $S$ for the product of the assigned full powers of $N_0-1$. Then $p\nmid U S$. Allocate each such predecessor power to one of the pairs of parts it divides. With $S_1S_2S_3=S$ the resulting pair products, $US_i$ divides the corresponding gcd because $\gcd(U,S)=1$. Thus $U^3S\le E^3$ and $S\le E^3$. The normalized companion now has sum $p^a U$, giving $27(N_0-1)^2\le(p^aU)^3 S^2\le p^{3a}E^6$. This proves the second claim. Lean 4.36.0 checks the affine determinant divisibility and the last size comparison for the successor. Run `python3 check_block_ratio_lemmas.py` for 15,426 exact carry-free splits, single-distance controls, and direct checks on feasible and impossible companion rows. [Formal arithmetic lemmas](BlockRatioBounds.lean)

The same argument bounds the common divisor of a specified pair quantitatively. For any positive row $A$ avoiding $p$ and any positive row $C$, let $G=\gcd(N(N-1),K_A,K_C)$ where $K$ denotes the row coefficient. Under the disjoint-ratio hypothesis,
$$G^3E^7\ge27D^2.$$
In fact write $G=G_N G_M$ for its parts from $N,N-1$. Remove those common prime-power layers and assign the retained powers to a row of lesser valuation. A retained power from $N$ divides every part of its assigned row, and residues at a retained power from $N-1$ sum to one. The full $p$-power is retained at $A$ since $A$ avoids $p$. The factorizations now have the form $N=G_NUV$, $D=G_MST$. The two unshifted $x,y$ remain positive proper parts of $D$, including when the shifted normalized entry $z$ is zero. The same determinants give $U\le E$, $S\le E^2$. The normalized parts of $C$ sum to $G_NU$, so $27T^2\le(G_NU)^3$, proving the claim. With $p\nmid c$ the unshifted weight likewise satisfies $G^3p^{3a}E^6\ge27(N_0-1)^2$, using $G=\gcd(N_0(N_0-1),K_A,K_C)$ there. [Partial-power assignment](endpoint_pair_bound.md#retaining-partial-prime-powers)

## One shared distance can also work

Suppose $B=1+p^u$, $u\ge1$, and $c$ has digits zero or one. Assume its support has exactly one pair at distance $u$, at positions $v,v+u$, and at least two other positions. Put $R=c-(1+p^u)p^v$. All outer exponents $a$ work whenever
$$27D^2>E^7,\qquad 27D^2>P^{3/2}(R+1)^3R^2. \tag{2}$$
In particular this holds for every sufficiently large $h$ with no condition on $a$. More generally, allow larger digits in $c$ provided its unique occupied pair at distance $u$ has smaller digit one, and $R>0$ is not a power of $p$.

A proper part cannot be proportional to $D$ in both blocks. Otherwise the two upper units split one to each side and the lower parts must have ratio $p^u$. The smaller lower part can occupy only positions from which a second occupied position lies $u$ places above. This forces those lower parts to be $p^v,p^{v+u}$, which leave $R>0$. Thus $S\le E^2$ as before. By the first inequality in (2), an exception would require $U>E$. Its two unshifted forms must then be proportional and use both blocks. Up to order they are
$$x=P+p^v,\qquad y=p^u x.$$
Since $P>c\ge p^v$, put $n=h-v>0$. We have $U\mid p^n+1$ and $p^a R\equiv-1\bmod U$, using $D=Bx+R$. Every full power in $S$ divides $x$, including the power of $p$ in $D$ if present, so $S\mid R$.

Choose a nearest integer $q$ to $a/n$ and write $a=qn+r$ with $|r|\le n/2$. Reduction modulo $U$ gives $U\mid p^r R\pm1$ if $r\ge0$, or $U\mid R\pm p^{-r}$ if $r<0$. Each integer is nonzero because $R$ is positive and is not a power of $p$. Hence $U\le P^{1/2}(R+1)$. Substituting in (1) contradicts the second inequality in (2).

For an application with both upper positions free, let $c$ have at least four occupied base-$p$ positions, all with digit one, and assume their positive pairwise distances are distinct. For all $n>h\ge1$ and $a\ge1$ the conjecture holds at
$$1+p^a(p^n+p^h+c)\qquad\text{if}\qquad p^n>(2c)^{64}.$$
Put $u=n-h$. When $p^u>c$ this is the earlier [two-high-power theorem](two_high_powers.md). Otherwise $B=1+p^u\le 2c$, $E\le2c^2$, and $P=p^h>(2c)^{64}/c>(2c)^{12}$. There is either no pair at distance $u$ in the low block, or exactly one with at least two units left over. In the former case $27D^2>E^7$. In the latter (2) also holds since $R\le c$ and $P^{1/2}>(2c)^6>(R+1)^3R^2$. This removes the earlier separation condition $p^{n-h}>c$ for the stated offsets.

## Every binary gap in a finite set of block pairs

Take $2\le k,s\le6$. Write $M_{k,s}=\sum_{j=0}^{s-1}2^{jk}$. Use the following pairs $(B,c)$.

* If $k+s\le10$, take any odd $B$ whose bits are a subset of those of $M_{k,s}$ and any $1\le c<2^k$.
* For every $k,s$ in the stated range, also take $B=M_{k,s}$ and $c=2^k-1$.

**For all 1,674 distinct pairs and every $a,h\ge1$, the full conjecture holds at $1+2^a(B2^h+c)$.** For example this includes $B=1+2^6+2^{12}+2^{18}+2^{24}+2^{30}$ and $c=63$, with twelve freely allocated binary units in the cofactor when $h\ge6$.

The distance test and the large-gap inequality leave 71,129 distinct $D$ to check, including overlaps $2^h\le c$. At each such $D$, enumerate disjoint positive binary submasks $x,y$ with $x+y<D$. Let $G$ be the largest divisor of $\gcd(x,y)$ coprime to $2D$. Using $U\mid G$ and $S\mid\gcd(D,xy)$ in (1) leaves 1,351 pairs from 92,658,606. For a full power $w=q^j$ of $D$ at an odd prime, if $A$ avoids $q$ then each unshifted normalized entry is zero modulo $w$ or is a unit modulo $q$, and at least one of $x,y$ is zero modulo $w$. The power at $2$ divides both unshifted entries. Put all powers failing these necessary conditions into a product $T$. Each must be taken by the companion. Hence $27T^2\le G^3$ and $2w+1\le G$ for every $w\mid T$ on this list. Exactly 835 pairs survive.

Multiply the normalized companion by $G/U$. Its three positive coordinates sum to $G$ and at each required $w$ their residues in some order are $G,0,0$. This follows since the original three actual parts are $1,0,0$ in some order and their scale is $N/U\equiv U^{-1}\bmod w$. Enumerate the slot taking $G$ at each $w$ and use the least positive CRT representative modulo $T$ in each coordinate. The three representatives always sum to more than $G$, excluding all 835 cases.

Run `python3 check_block_ratio_families.py` for this computation. Run `python3 check_block_ratio_families_cross.py` for an independent check. The second program first absorbs all powers of two in $D$ into the free outer exponent. It builds each unshifted pair from its sum, then constructs possible single companion coordinates from subsets of prime powers. Across 39,304 distinct odd cofactors it reaches 558 final CRT tests and excludes every case.
