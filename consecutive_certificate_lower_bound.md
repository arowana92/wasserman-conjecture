# How many primes must one list contain to certify consecutive weights?

There are arbitrarily long blocks of **consecutive weights**, all satisfying Wasserman's three-row conjecture, such that any one prime list certifying the whole block must contain at least half as many primes as there are weights. All the primes forced by the explicit witnesses exceed the length of the block. For instance, the construction below gives **128 consecutive weights above $2^{4096}$ requiring at least 64 distinct primes greater than 1,000** in any common certificate.

For a quantitative version fix $0<b<a$ and $\theta>a+4b$. For every sufficiently large $X$, the interval $[X,X+X^\theta]$ contains $H=\lfloor b\log X\rfloor$ consecutive weights with the following property. At each weight there is an explicit selection of three rows whose gcd exceeds one and has no prime divisor at most $y=\lfloor a\log X\rfloor$. Choosing a common prime for each of these $H$ selections uses at least $\lceil H/2\rceil$ distinct primes. In particular $\theta=1$ works when $a+4b<1$. When $X$ is a power of two and $a+4b<4/7$, one can place the whole block in an interval where the **full conjecture** is proved for every weight. Logarithms are natural.

## The exact gcd at each position

Take any integers $H\ge1$ and $y\ge H+1$. We construct a single congruence class for $B$. At each of the $H$ weights $N=B+h$, $2\le h\le H+1$, the three rows are
$$
 (1,1,N-2),\quad(h,A_1,B-A_1),\quad(h,A_2,B-A_2).
$$
The $A_i$ depend on $h$. For every sufficiently large positive $B$ in the class their coefficient gcd is exactly
$$
               G_h=\frac{B+h}{h}\,\frac{B+h-1}{h-1}.                     \tag{1}
$$
Each displayed factor is an integer greater than one with all prime factors greater than $y$. They are coprime. These triples therefore exhibit at least two distinct common primes at each weight.

Put $T=H+1$. For $p\le T$ let $Q_p$ be the least power of $p$ strictly greater than $T$. At offset $h$ only the primes dividing $h(h-1)$ need to be assigned. Put those with $p^3\le T^2$ in the first group and the rest in the second. Set $A_i=\prod_{p\text{ in group }i}Q_p$, using 1 for an empty product. For each $p\le T$ let $J_p$ be the least power of $p$ strictly exceeding all $A_i$ of groups containing $p$ over the whole block. Every $p\le T$ occurs at $h=p$. Impose
$$
 B\equiv-Q_p\pmod{J_p}\quad(p\le T),\qquad
 B\equiv0\pmod p\quad(T<p\le y).
$$
All moduli are coprime. Take $B$ greater than every $A_i$ to ensure positivity. Since $Q_p>T$, the powers of every $p\le y$ in $B+h$ and $B+h-1$ are exactly its powers in $h$ and $h-1$.

At an assigned prime $p$, $Q_p\mid A_i$ and $A_i<J_p$. The digits of $B$ below $Q_p$ are zero and all its digits from $Q_p$ up to $J_p$ are $p-1$ (positions with place value less than $J_p$). Subtracting $A_i$ from $B$ needs no borrow. Adding $h<Q_p$ uses only the empty low positions. Thus $(h,A_i,B-A_i)$ is carry-free there. Every prime through $y$ is either absent from the endpoint or avoided by one of the other rows.

For a coefficient $K$ on a row containing $h$, twice decrementing that part gives $N(N-1)\mid h(h-1)K$. This proves that the integer in (1) divides both coefficients. The endpoint coefficient is $N(N-1)$, and the factor $h(h-1)$ is coprime to the joint gcd of the other two by the carry test. This proves equality. [Lean verification of the factorial and gcd deduction](PairDivisor.lean)

## Why the list has to grow

Every prime in (1) exceeds $y\ge H+1$. Such a prime divides a weight or its predecessor at at most two positions in the block, and these must be adjacent. At adjacent offsets $h,h+1$, the exact gcd of the two witness gcds is $(B+h)/h$. This factor is common, and any possible further factor would divide both $B+h-1$ and $B+h+1$ and hence 2. Both gcds are odd. At nonadjacent offsets their gcd is one. Pairing consecutive offsets therefore shows that the minimum number of primes hitting **these specific witnesses** is exactly $\lceil H/2\rceil$. Any list certifying all triples at every weight in the block must be at least this large.

## Length and placement of the block

Let $L$ be the product of the CRT moduli. For $p\le T^{2/3}$, each group at a given offset has $O(\log T)$ members because their distinct primes divide $h(h-1)\le T^2$. Each $\log Q_p\le2\log T$. This gives $\log J_p=O((\log T)^2)$ and a combined cost $o(T)$ from these bases. For $p>T^{2/3}$ a group contains at most one prime of $h$ and at most one of $h-1$. Here $Q_p=p^2$, and $A_i\le p^2T^2<p^5$, so $J_p\le p^5$. With $\vartheta(x)=\sum_{p\le x}\log p$ we obtain
$$
       \log L\le\vartheta(y)+4\vartheta(T)+o(T).
$$
The prime number theorem gives $\vartheta(x)\sim x$. With the stated $y,H$, $L\le X^{a+4b+o(1)}$, and every $A_i=\exp(O((\log T)^2))=o(X)$. Choosing the first permitted $B\ge X$ puts $B+2,\ldots,B+H+1$ in $[X,X+X^\theta]$. For $X=2^m$ and $a+4b<4/7$, choose $\theta$ strictly between these numbers. Then $(N-X)^7\le6^6(N-1)^4$ throughout the block for large $m$, so the [proved prime-power interval criterion](all_weights_endpoint_reduction.md#consequences) applies. [Selberg, *An elementary proof of the prime-number theorem*](https://www.math.lsu.edu/~mahlburg/teaching/handouts/2014-7230/Selberg-ElemPNT1949.pdf)

## Finite checks

The [checker](check_consecutive_certificate_lower_bound.py) constructs the congruence class, tests the separate factorial valuations, checks (1) by multiplying exact coefficients in 234 small cases, and checks every pairwise gcd of the 128 witnesses in the large example. For $H=128$, $y=1000$ it finds a 2,117-bit modulus. Taking the first permitted $B\ge2^{4096}$ gives a 2,116-bit offset $B-2^{4096}$. The checker verifies the explicit interval criterion for every one of the 128 weights. The Lean theorem checks the coefficient identity underlying (1) and its gcd consequence given the factorization and coprimality hypotheses. The digit construction and interval estimate are proved above.

~~~sh
python3 check_consecutive_certificate_lower_bound.py
lean +leanprover/lean4:v4.36.0 PairDivisor.lean
~~~
