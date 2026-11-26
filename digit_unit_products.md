# Prime certificates from the product of three parts

Let $K_A$ be the trinomial coefficient of a positive triple $A$ of sum $N$, and let $Z_A$ be the product of its parts. For a prime $p$, list the base-$p$ digit units of $N$ in increasing order as powers of $p$ with multiplicity. If there are fewer than three units, $p$ divides every coefficient. Otherwise let $b_1,b_2,b_3$ be the three lowest units. Define
$$D_p=b_1b_2b_3,\qquad L_p=b_1b_2(N-b_1-b_2).$$
To define $H_p$, assign the units from largest to smallest, putting the next one into a currently smallest part. Take the product of the resulting three parts.

**Exact one-prime bounds.** Among all rows avoiding $p$, the greatest power of $p$ dividing every $Z_A$ is $D_p$. The smallest and largest values of $Z_A$ are exactly $L_p$ and $H_p$.

**A family with three separated digit blocks.** The full conjecture holds for every
$$N=p^a(1+p^b+t p^c),\qquad p\text{ prime},\quad 1\le a\le b<c,\quad 1\le t\le800p^{b-a}.$$
The last block can have arbitrarily many digit units. There is no upper bound on any exponent. Also, for every $b\ge1$, the conjecture holds throughout the simultaneous residue class
$$N\equiv2\pmod{2^{b+1}},\quad N\equiv3\pmod{3^{b+1}},\qquad 3\le N\le864\cdot 6^{2b}-3.$$

**An optimal four-prime certificate.** Every three rows of weight $1{,}000{,}000{,}000{,}232{,}000{,}000{,}007{,}527$ have a common prime among the four given below. Here four is necessary even if the proposed list can use any primes.

## Proof of the extrema

By Kummer's criterion a row avoiding $p$ partitions the digit units among three nonempty parts. The $p$-valuation of each part is the exponent of its lowest assigned unit. These are three distinct occurrences, so their product is divisible by $D_p$. Placing the first two units alone and giving the third part everything else attains the claimed valuation. It also minimizes $Z_A$. For any sorted positive partition the smallest part is at least $b_1$ and the second at least $b_2$, since at most one unit can be strictly smaller than $b_2$. Transferring excess from these smaller parts to the largest only lowers the product.

For the maximum, process one digit position at a time. After all units of size $w$ are assigned, the total of all smaller units is less than $w$. If a position holding a newly assigned $w$ has partial load at least $2w$ above another position, its final load exceeds the other by more than $w$. Moving this unit increases the product in any positive optimal partition, a contradiction. With the previous larger-unit loads fixed, this condition determines the sorted partial loads uniquely. In units of $w$, let $m$ be the final minimum at this stage. Each position receiving units ends at $m$ or $m+1$ and any higher position keeps its initial load. Thus the final multiset is obtained by raising the smallest initial loads together until the current digit is exhausted. Induction proves the greedy maximum. For $p=2$, if $Q,R$ are the two highest units and $r=N-Q-R>0$, the rest sum to less than $R$. In this case $H_2=QRr<N^3/32$, since the largest part exceeds $N/2$.

## A test using all primes of the endpoint

Put $M=N-1$. In any counterexample, the [endpoint theorem](all_weights_endpoint_reduction.md) supplies the row $(1,1,N-2)$, whose coefficient is $NM$. Its other two rows $A,C$ must between them avoid every prime of $NM$. Assume none of these primes is universal and let
$$\Lambda(N)=\frac{\prod_{p\mid NM}D_p}{N^3 M^2}.$$
This is a positive integer. For $p^e\parallel N$ all three units used in $D_p$ are at least $p^e$. For $p^e\parallel M$ the first is one and the next two are at least $p^e$. The remaining factor is exactly the local contribution to $\Lambda$. Equivalently it is the product of the three lowest units of $N/p^e$ in the first case, or the two lowest units of $M/p^e$ in the second.

By the one-prime bounds, $N^3 M^2\Lambda(N)\mid Z_A Z_C$. One row avoids 2 and the other has product at most $N^3/27$. Thus a counterexample requires
$$27M^2\Lambda(N)\le H_2=QRr<N^3/32. \tag{1}$$
The contributions to $\Lambda$ from *any subset* of endpoint primes already give valid lower bounds on its value. In particular if their product is at least $W$, then $N+3\le864W$ suffices to prove the conjecture. For $N\ge5$, $N^3/(N-1)^2<N+3$. The two smaller weights have universal primes.

For $N=p^a(1+p^b+t p^c)$ the first three units after removing $p^a$ are $1,p^b$ and a unit at least $p^c$. Consequently $\Lambda\ge p^{b+c}$. Under the displayed conditions $N+3\le864p^{b+c}$ since $p^a t\le800p^b$ and $p^a(1+p^b)+3\le64p^{b+c}$. For the simultaneous residue class, the quotients $N/2,N/3$ are both one modulo the corresponding $b$-th prime powers. Thus $\Lambda\ge2^{2b}3^{2b}$ unless a prime is universal. The period of the combined class is $6^{b+1}$.

In the [fixed-part notation](million_part_bound.md), if the companion $C$ is required to avoid $p$, its normalized product $m_0=y_1y_2y_3/T^2$ contains the corresponding local factor of $\Lambda$. Indeed $Z_C=V^3T^2m_0$ and the full baseline powers are in $V^3$ or $T^2$ for primes assigned to $C$. In particular, if $C$ avoids 2 and takes that prime in the assignment, then $8\mid m_0$ at even $N$ and $2\mid m_0$ at odd $N$. In binary the three smallest distinct positions begin no lower than $0,1,2$ in $N/2^{v_2(N)}$; the two lowest in the odd quotient of $M$ at an odd weight begin no lower than $0,1$. Since $m_0 T^2\le U^3/32$ by the binary upper product and $S\mid(a/U)(a-1)$ for every part $a$ of the opposite row, this gives $a^5>256M^2$ at even weights and $a^5>64M^2$ at odd weights for those parts.

## Lists of candidate common primes

For a set of bases assigned to one avoiding row its product must belong to
$$\big[\max_p L_p,\ \min_p H_p\big]\ \cap\ \big(\prod_p D_p\big)\mathbb Z.$$
An empty intersection excludes the assignment. When looking for a certificate on the original three rows, every prime in the proposed list would have to be covered by one of the three avoiding rows. Thus if every partition of the list into at most three groups has an impossible group, it certifies a common prime for the original problem directly.

At the displayed weight the four prime bases are
$$P=1{,}000{,}000{,}000{,}039,\quad Q=1{,}000{,}000{,}000{,}193,\quad R=(N-1)/2,\quad J=N-400.$$
Here $N=PQ$. The bases $P,Q$ force three multiples of the corresponding prime and $R$ forces two multiples of $R$. Thus $D_P\ge P^3$, $D_Q\ge Q^3$, $D_R\ge R^2$. The product of any two of these lower bounds exceeds $(N/3)^3$. A row avoiding $J$ has product at most $N\cdot200^2$ since its other two parts sum to at most 400. This is smaller than each of the three lower bounds. No row can avoid a pair of the four primes.

Every prime can individually be avoided at this weight. In base $P$ or $Q$ the digit sum of the other odd prime is odd and cannot be one, since the primes are distinct. A base prime not dividing $N$ with at most two units would require $N-1$ to be a prime power, whereas $N-1=2R$. Thus any proposed set of at most three primes can be avoided by choosing one row for each. The verifier recursively certifies the four primalities by modular exponentiation and gcds.

## Reproduction and scope of the test

Through weight 120 the verifier compares all three one-prime extrema with actual trinomial coefficients, including 1,163,432 prime and row comparisons. It checks 2,880 three-block specializations, 80 points in the simultaneous residue classes and 2,922 binary maxima.

At the 10,000 consecutive weights immediately after $10^{12}$, 839 have a universal prime. Among the others, inequality (1) settles 39 and the two-group interval test on all endpoint primes settles another 318. For the remaining 8,804 the product restrictions alone admit an assignment. A separate [CRT and factorial-floor verifier](check_power_partner_assignments.py) checks the entire interval using the [two-row reduction](all_weights_endpoint_reduction.md). At 7,008 weights, 153,837 assignments pass its initial size bounds. None has the required row avoiding all its assigned primes. This proves the full conjecture throughout this interval.

Run `python3 check_digit_unit_products.py --crt` with python-flint and SymPy installed.
