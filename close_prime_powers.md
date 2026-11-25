# Weights just above a multiple of a prime power

Write $P=p^a$ with $p$ prime and $a\ge1$. Let $Q=H+d$, $H=p^j$ and $N=PQ$, where $d>0$ and $p\nmid d$. Thus $P$ is the full $p$-power in $N$.

1. **Unbounded family.** The conjecture holds for every $j\ge a$ and every $1\le d\le1000$, with no upper limit on $p,a,j$. In particular it holds for $P(P+d)$ throughout this range with an unrestricted cofactor $P+d$ coprime to $p$.
2. **Size bound.** For any $j\ge1$ with $d<H$, a counterexample requires $432(N-1)^2\le P^3 d^4$. If $H=P$, it requires the stronger inequality
   $$27(N-1)^2\le P^3(4d^2/25+1)^2.$$
3. **Counting near a square.** For $Y\le P\le2Y$ and $1\le d\le D<Y$ with $H=P$, the number of possible exceptional weights is at most $Y^{1/2+o(1)}D^{5/2}$. When $Q=P+d$ is also a prime power the bound improves to $Y^{1/2+o(1)}D^2$. The $o(1)$ is uniform in $D$ as $Y$ tends to infinity. In particular the conjecture holds for a proportion tending to one in the two respective classes for $D\le Y^{1/3-\epsilon}$ and $D\le Y^{1/2-\epsilon}$, for every fixed $\epsilon>0$ whenever the classes are nonempty.

## A row at the prime dividing the weight

By the [all-weight reduction](all_weights_endpoint_reduction.md), a counterexample contains the endpoint $(1,1,N-2)$, whose coefficient is $N(N-1)$. Assign the bases of $N(N-1)$ to the two other rows, each prime to a row avoiding it. The row $R$ taking $p$ has all parts multiples of $P$. After division its parts add without carries in base $p$ to $Q=H+d$. When $d<H$, the base-$p$ expansion of $Q$ has digit 1 in position $j$. Exactly one part receives that unit $H$. Write the normalized row as $(x,y,Q-t)$ where $t=x+y\le d$, $x,y>0$.

Assign each distinct prime with its full power in $N$ or $M=N-1$. Let $Pv$ be the product assigned to $R$ from $N$. Then $v$ divides $Q$, $\gcd(v,Q/v)=1$, and $v\mid g=\gcd(x,y)$. Such a divisor $v$ is called unitary. Write $S$ for the product assigned to $R$ from $M$, and $T=M/S$. The other row has scale $Q/v$ and its three normalized parts sum to $Pv$. Each full power in $T$ divides two of them, since the no-carry condition forces two original residues to vanish and the scale $Q/v$ is coprime to $T$. Their product is therefore divisible by $T^2$, and the arithmetic–geometric mean inequality gives $27T^2\le(Pv)^3$. Thus
$$27M^2\le (Pv)^3S^2. \tag{1}$$
Every full power in $S$ divides $x$ or $y$. Since $S$ is coprime to $v$, it divides $\operatorname{lcm}(x/v,y/v)=xy/(vg)$. In particular
$$v^{3/2}S\le xy\sqrt v/g\le d^2/4. \tag{2}$$
This proves the general size bound. These arguments retain the assigned scale in the [two-part gcd bound](endpoint_pair_bound.md).

## Two congruences near a square

Now let $H=P$. Split $S=AB$ according to the positions divisible by each full prime power in $S$. Put a power in $A$ if it divides both $x$ and $y$, and in $B$ otherwise. In the latter case it divides $Q-t$. Hence
$$A\mid g,\quad B\mid P+d-t,\quad B\mid1+t(d-t),\quad \gcd(vA,B)=\gcd(v,A)=1. \tag{3}$$
For the third divisibility use $PQ\equiv1\pmod B$, $Q\equiv t\pmod B$ and $P\equiv-(d-t)\pmod B$. Also $B\mid\operatorname{lcm}(x,y)$, and $B$ is coprime to $t$.

We claim $v^{3/2} S\le 4d^2/25+1$ for every $Q$. From (2) and $v\mid g$ we have $v^{3/2}S\le t^2/(4\sqrt v)$. For $v\ge3$ this already gives the claim. If $g\ge2v$, (2) gives $v^{3/2}S\le d^2/(8\sqrt v)$. This leaves $g=v\in\{1,2\}$, so $A=1$ and $S=B\le1+t(d-t)$.

- For $v=1$ split at $t=4d/5$. Below or at this point $t^2/4\le4d^2/25$. Above it $t(d-t)\le4d^2/25$.
- For $v=2$ split at $t=19d/20$. Below or at this point $t^2/(4\sqrt2)\le361d^2/(1600\sqrt2)<4d^2/25$. Above it, if $t=d$ then $S=1$ and $d\ge4$ suffices. Otherwise $1\le d-t<d/20$, so $d>20$ and $2^{3/2}S\le 3(1+d^2/20)\le1+4d^2/25$.

Inserting the claim into (1) proves the stronger bound in item 2. When $Q$ is a power of a prime distinct from $p$, assigning it to the first row would make all three parts multiples of $N$. Thus $v=1$ in this subfamily.

The constant $4/25$ is sharp for these two numerical bounds with $g=1$. For $k\ge1$ put $z=k^2+k+1$, $x=2z-k$, $y=2z+k+1$, $t=x+y$ and $d=t+z$. Then $(x-z)(y-z)=z^2+1$, so $xy=1+tz$, $\gcd(x,y)=1$, and $xy/d^2\to4/25$. Taking $P\equiv-z\pmod{xy}$ satisfies the low congruences for $B=xy$. This example concerns the congruences and size estimates only.

## Counting all possible weights in a short range

Fix $d,t$ with $2\le t\le d$. The choices $v\mid t$, $A\mid t/v$, $B\mid 1+t(d-t)$ give at most $Y^{o(1)}$ combinations by the divisor bound. Here $A\mid t/v$ because $A\mid t$ and $A,v$ are coprime. For each combination the necessary residue conditions on $P$ are
$$P(P+d)\equiv1\pmod A,\qquad P\equiv-(d-t)\pmod B,\qquad P\equiv-d\pmod v.$$
The three moduli are pairwise coprime. The first congruence has at most $O(2^{\omega(A)}\sqrt{\gcd(A,d^2+4)})$ roots. For an odd prime power complete the square to $(2P+d)^2\equiv d^2+4$. If the right side vanishes modulo the power $r^e$ there are $r^{\lfloor e/2\rfloor}$ square roots. Otherwise a solution requires an even valuation $2b$, and lifting the two possible unit roots gives at most $2r^b$ solutions. For powers of two there are at most four unit roots and at most two preimages of $2P+d$ per residue, which changes only the absolute constant.

From (1) we have $v^{3/2}S\gg\sqrt Y$. Thus each of the allowed classes modulo $vS$ has $O(\sqrt{Yv})$ representatives in $[Y,2Y]$. Put $F=d^2+4$. For a fixed $v\mid t$ all the possibilities together therefore contribute at most $Y^{1/2+o(1)}\sqrt{v\gcd(t/v,F)}$ representatives. Writing $t=vh$ and summing uses
$$\sum_{h\le Z}\sqrt{\gcd(h,F)}\le\sum_{b\mid F}\sqrt b\lfloor Z/b\rfloor\le Z\tau(F).$$
For fixed $d$, summing over $v\le d$ now gives $Y^{1/2+o(1)}d^{3/2}$. If $v=1$ is forced, it gives $Y^{1/2+o(1)}d$. Sum over $d\le D$ to obtain item 3.

For the relative estimates, (2) excludes all offsets $d\le Y^{1/4}$ in both classes for large $Y$. At larger $D$, primes $P\in[Y,2Y]$ give $\gg YD/\log Y$ weights in the unrestricted class. Each has at most one full prime-power divisor in $(\sqrt{N/2},\sqrt N)$, since the product of two distinct such powers would divide $N$ strictly between $N/2$ and $N$. In the second class, divide $[Y,2Y]$ into $O(Y/D)$ intervals of length at most $D$. If their prime counts are $a_i$, Cauchy–Schwarz and the prime number theorem give $\sum_i\binom{a_i}2\gg YD/(\log Y)^2$ distinct prime products for $D>Y^{1/4}$. Dividing the two upper bounds by these respective lower bounds proves the stated proportions.

## A finite computation for every exponent

Use the [verified initial interval through one billion](full_range_1b.md). If $d\ge H$ and $j\ge a$, then $N\le2d^2$, inside this interval for $d\le1000$.

For $j=a$ enumerate $2\le t\le1000$, $x\le y$ with $x+y=t$, $v,A\mid g$ coprime and $B$ dividing the part of $\operatorname{lcm}(x,y)$ coprime to $t$. By (3) the only possible offsets satisfy $d\equiv t-t^{-1}\pmod B$, including all offsets when $B=1$. Enumerate $P$ in the three CRT conditions above. Equation (1) and $M>P^2$ give the finite cutoff $27P<v^3S^2$. Retain full prime powers $P$ and unitary factors $v$ and $S$, test the actual row in the bases assigned to it, and search the complementary normalized rows by their CRT residues and all positive lifts.

For $j>a$, (2) and (1) require $432p^{2j-a}<d^4$ since $M>PH$. As $2j-a\ge a+2\ge3$, this bounds all of $p,a,j$. For each remaining weight the program factors $Q$ and $M$, tries each unitary $v\mid Q$ with $2v\le d$, and all assignments to $S$ of full powers $w$ of $M$ with $vw\le d$. This last condition follows because $w$ divides $x/v$ or $y/v$. It checks (1), $4v^2S\le d^2$, and $2w+1\le Pv$ for each power left to the complementary row. The last condition uses two positive multiples of $w$ in its normalized parts. Every possible small part $x$ is a multiple of $v$ with $S\mid Px(Px-1)$. Pair these values with $x+y\le d$ and apply the full carry tests.

Run `python3 check_close_powers.py` with python-flint installed. The optional integer argument changes the maximum offset. All factorizations used by the verifier are checked by multiplication and primality proofs.

At maximum offset 1,000, the square-case progressions produced 546,364,102 points after the initial interval and size tests. Of these, 87,509,176 gave full prime-power assignments, counted with their row and assignment choices. The first row passed all of its carry tests in 257,927 cases. After testing the remaining full powers individually, 157,543 cases reached the CRT search. **None admitted a complementary row even at the required residues.** For $j>a$ the size bound left 68,184 weights above one billion and 12,324 factor assignments after the individual size tests. None admitted a first row passing its carries.

The verifier also checks the inferred divisibilities against 4,845 carry-free first rows using exact trinomial coefficients at small weights. A separate direct residue check includes 22 solutions with both $A$ and $B$ nontrivial.

The independent weight-first check enumerates the possible $P$ before the rows. Combining (1) with $4v^2S\le d^2$ gives $432Pv<d^4$. For each eligible $P,v,d$ it takes candidate divisors $S$ from $\gcd(N-1,\operatorname{lcm}(1,\ldots,\lfloor d/v\rfloor))$. It then pairs parts satisfying $S\mid x(x-Q)$ and tests divisibility by $S$ of the least common multiple of the three pairwise gcds of $x,y,Q-x-y$. Through offset 160 this gives exactly the same 1,111 preliminary prime-power assignments as the progression enumeration. Through offset 300, 87 assignments reach the complementary CRT search. Enumerating each allowed single part and testing sums gives the same empty set as assigning zero positions to the three parts. The two CRT algorithms also agree with direct enumeration at normalized sums 101 and 641 for the powers 4, 5 and 9, with 1 and 30 solutions respectively, including positive lifts beyond the modulus. Run `python3 check_close_powers_cross.py` for these comparisons.
