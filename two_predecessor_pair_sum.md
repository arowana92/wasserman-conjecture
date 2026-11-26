# Two predecessor prime powers control every pair of parts

Let $P,Q>1$ be powers of **different** primes dividing $M=N-1$, with $PQ\mid M$. In a counterexample the sum $h$ of any two parts of either row other than the endpoint satisfies
$$h>PQ.$$
For $h\le2PQ$ the only possibility is **$h=PQ+1$**, and $\gcd(N,PQ+1)$ must contain a unitary divisor $U$ of $N$ with $3\le U\le(PQ+1)/(P+Q)$. In particular this possibility is excluded if the gcd is at most two or $PQ+1$ is a prime power. There is no restriction on the remaining factors of $N-1$.

Consequently take any prime power $L$ whose base divides $N(N-1)$. If $0\le N-L\le PQ$, the full three-row conjecture holds at $N$. The offset may be as large as $2PQ$ under either of the two exclusions just given, provided $L>N/2$. For the first assertion, if $L\le N/2$ then $PQ\ge N/2$ and $M=PQ$, a weight [already settled](all_weights_endpoint_reduction.md#consequences). For $L>N/2$ a row avoiding its base puts $L$ in one part, and the other two have sum at most $N-L$.

In particular, for fixed $P,Q$, **the first integer at or above each power $2^k$, $k\ge1$, that is $1\bmod PQ$ satisfies the conjecture.** If $PQ+1$ is a prime power, the next integer in that progression works as well when the power of two exceeds $2PQ$. The same progression result holds for a leading power of either prime below $P,Q$.

The [modular-inverse analysis](inverse_pair_digit_gap.md) also settles the entire second interval above a binary power whenever $\min(P,Q)\le128$, with the larger power unrestricted. At every eligible marker the cofactor $(N-1)/(PQ)$ being one above a prime also suffices.

## The two assignments

Let $A,C$ be arbitrary positive rows and put $D=\gcd(N,K_A,K_C)$, $J=\gcd(M,K_A,K_C)$ and $G=DJ$. For every pair sum $h$ of $A$,
$$Gh>PQ. \tag{1}$$
Remove from each of $P,Q$ its gcd with $J$, giving $P',Q'$ and $P'Q'\ge PQ/J$. Under the [retained-power assignment](endpoint_pair_bound.md#retaining-partial-prime-powers) each surviving power is assigned to one of the two rows. At that row exactly one actual residue is 1 and two are 0. A reduced power equal to one can be assigned to either row and given any exceptional position. Retaining from $N$ gives scales $U,V$ with $DUV=N$. A row's entries are multiples of its scale; each scale is coprime to $PQ$.

* When both powers go to $A$, one selected part is divisible by $UP'Q'$ unless the two exceptional positions are precisely the two selected positions, in opposite assignments. In that remaining case $h-1$ is a positive multiple of $P'Q'$. Thus $h>P'Q'$.
* When both go to $C$, one of its three normalized entries is divisible by $P'Q'$. Their sum $DU$ strictly exceeds $P'Q'$, and $h\ge2U$.
* When $P'$ goes to $A$ and $Q'$ to $C$, two companion entries after normalization are multiples of $Q'$, giving $DU\ge2Q'+1$. At least one selected part of $A$ is divisible by $UP'$, so $h\ge U(P'+1)$. Thus $Dh>2P'Q'$, also with the roles of the powers exchanged.

This proves (1). In a counterexample apply the [endpoint reduction](all_weights_endpoint_reduction.md) to get $D=J=1$, $U,V\ge3$. If $h\le2PQ$, only the opposite assignments in the first case remain, forcing $h=PQ+1$. Here both selected parts are divisible by $U$, and their quotients by $U$ are positive multiples of $P,Q$ respectively. This proves the bound on $U$ and the gcd exclusion. In any nonendpoint counterexample row a pair sum cannot be a prime power: if $h=s^t$, then $U\mid h$ requires $U$ to be a positive power of $s$, whose prime the row avoids. But $K_A$ is a multiple of $\binom{h}{x}$ for either selected part $0<x<h$, and $s$ divides every such binomial coefficient by the carry into the top digit.

More precisely, for every prime $s\mid U$, the full $s$-power of $h$ divides both selected parts. Indeed their addition to $h$ has no carry in base $s$, since its binomial coefficient divides $K_A$. Write $H$ for the product of these full powers. In the opposite assignments one needs $h/H\ge P+Q$. If $A$ avoids the base of a leading power $L>N/2$, the sum of the two parts missing $L$ must also use only digits from the remainder $N-L$. These tests can exclude the second offset even with a larger gcd.

For $L>N/2$ either offset criterion gives an explicit prime list for every **pair** of arbitrary rows. Take the primes dividing $N$ together with the bases of $P,Q,L$. If the pair had no common prime in the list, assign each prime in the list to a row avoiding it. The weight primes again give $UV=N$ with $U,V\ge3$, and one row avoids the base of $L$. The proof above uses just the two indicated powers of the predecessor and therefore rules out that pair. In particular (1) itself holds with the sharper $G=\gcd(NPQ,K_A,K_C)$.

There is a direct way to use this list for **three** rows. Suppose $R>N/2$ is a prime at most $N$. For each listed prime $s\mid N$ require $N-R<2s^{v_s(N)}$; for each listed $s\mid N-1$ require $N-R\le s^{v_s(N-1)}$. An $R$-avoiding row has two parts with sum at most $N-R$. Avoiding a listed weight prime would require both to be multiples of its full power, and avoiding a listed predecessor prime would require at least one such multiple. Thus all primes on the list divide the $R$-avoiding row, and adjoining $R$ certifies all triples.

## Entire rays of exponents

Fix distinct primes $p,q$ and an integer $r\ge2$. Consider $N=1+r p^a q^b$ with $a,b\ge1$. Whenever $q^b(r-1)<p^k\le q^b r$ for an integer $k$, **every $a$ works**. Indeed $L=p^{a+k}$ lies at offset at most $PQ$, with $P=p^a,Q=q^b$. The same conclusion holds on interchanging $p,q$.

The limiting fraction of these certified values of $b$ from the displayed condition is $\delta_p=\min(1,\log_p(r/(r-1)))$. This follows from equidistribution of the fractional parts of $b\log_p q+\log_p r$, an irrational rotation: the exponential sums for each nonzero integer frequency are bounded geometric sums, giving equidistribution by Weyl's criterion. The criterion in $a$ similarly gives $\delta_q$. In rectangular ranges with both exponent bounds tending to infinity, the lower limiting proportion proved by the two conditions together is at least $1-(1-\delta_p)(1-\delta_q)$.

If $r\ge4$ and all prime factors of $r-1$ divide $pq$, **replace $r-1$ by $r-2$ in both conditions and both proportions**. In fact $\gcd(N,PQ+1)=\gcd(r-1,PQ+1)=1$, so the second-offset criterion applies. Its leading-power hypothesis follows from the strict inequality $p^k>(r-2)q^b$, which gives $L>N/2$ even for $r=4$. For example $r=385$, $p=2$, $q=3$ gives unrestricted rows of exponents at weights $1+385\cdot2^a3^b$ with five distinct predecessor primes. All $a$ work on the stated fraction $\log_2(385/383)$ of $b$, and all $b$ on the fraction $\log_3(385/383)$ of $a$.

For this example **every block of 253 consecutive positive exponents $b$ contains one that works for all $a$**. Every block of 401 consecutive positive $a$ contains one that works for all $b$. Both block lengths are sharp for their respective rotation tests. More generally the covered exponents have bounded gaps for every fixed $p,q,r$: the rotations by $\log_p q$ are dense, so finitely many translates of any open acceptance interval cover the circle by compactness.

For an exact block-length test, sort the numbers $q^j/p^{\lfloor j\log_p q\rfloor}$, $0\le j<B$, in $[1,p)$. Take the largest quotient of successive entries, including the wraparound to $p$ times the first one. If it is at most $385/383$, every translated block hits the half-open acceptance interval on the logarithmic circle. If it is strictly larger, some open interval of starting phases misses, and irrational rotation visits that interval. For $(p,q)=(2,3)$ the largest quotient at $B=252$ is $2^{233}/3^{147}$, and at 253 it is $2^{317}/3^{200}$. For $(p,q)=(3,2)$ the same two quotients occur at $B=400,401$. Exact rational comparison puts $385/383$ strictly between them.

## A 181-bit example in the second interval

Take $P=2305843009213693951$, $Q=147573952589676412931$, $r=4503599627370498$. Then
$$N=1+rPQ=1532495540865889538289620652285986858360875510602399739.$$
Both $P,Q$ are prime, and $r=2\cdot3^2\cdot307\cdot2857\cdot6529\cdot43691$. Here $N-2^{180}$ lies strictly between $PQ$ and $2PQ$, has 116 binary ones and exceeds $N^{2/3}$. The gcd of $N$ and $PQ+1$ is one. Thus the second-offset theorem proves the full conjecture at this weight. The verifier uses exact integers and recursive primality certificates.

## Two minimum four-prime certificates at the same weight

Take $N=1208925867235455872460352=64S$ and $R=N-35$, with $P=134217689$, $Q=268947953$. The four numbers $S,R,P,Q$ are prime and $N-1=3^5\cdot283\cdot487\cdot P Q$. Here $PQ<N-2^{80}<2PQ$ and $\gcd(N,PQ+1)=2$. The general construction gives five primes. In fact **both $\{2,S,R,P\}$ and $\{2,S,R,Q\}$ certify every three rows**, and four is the minimum possible size.

Fix $w=P$ or $Q$. No row avoids two primes of $\{2,S,R,w\}$. For a row avoiding $R$ the two small parts sum to at most 35, whereas avoiding 2 or $S$ forces multiples of 64 or $S$, and avoiding $w$ requires a small multiple of $w$. A row avoiding $S$ has three positive entries $Sx_i$ with $\sum x_i=64$. If it also avoided $w$ then two $x_i$ would be multiples of $w>64$. Avoiding 2 and $S$ together would force every part to be a multiple of $64S$.

For the pair $2,w$, write $c=N-2^{80}$. The two parts missing the top binary digit are positive disjoint submasks of $c$. If their row avoided $w$, one submask would be 0 modulo $w$ and the other 0 or 1. There are four positive submasks in residue 0 for either prime. There are zero in residue 1 for $P$ and one for $Q$. Every required pair overlaps. The Python verifier finds these submasks by splitting the 29 occupied positions into groups of 14 and 15. An independent C++ verifier tests all $2^{29}-1$ positive submasks directly.

Four primes therefore suffice. For necessity, every prime is individually avoidable. The digit sums of $N$ at 2 and $S$ are 30 and 64. At any other prime a digit sum at most two would require $N=1+s^j$, since the units digit is nonzero and $N>2$, but $N-1$ has several distinct prime factors. Distribute at least three digit units into three nonempty parts to avoid any prescribed prime. A list of at most three can therefore be avoided by three rows between them.

With python-flint installed, run `python3 check_two_predecessor_pair_sum.py`. It compares (1) against exact coefficients in 299,683 row-pair and prime-power choices through weight 48, and checks all 275 resulting weight-and-prime-list choices through 240 using exact coefficients and prime masks. These include 105 lists omitting a predecessor prime. It also verifies both numerical examples. To repeat the second submask check, compile `check_pair_sum_certificate.cpp` with C++17 and run the executable.
