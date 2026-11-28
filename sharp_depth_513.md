# Four primes attain sharp depth 513 at a 490-digit weight

At the weight defined below, four specified primes certify every triple of positive trinomial rows. Four is the minimum possible size of a certificate. Each witness can be detected within 513 digits, and some triple requires digit 513 even when every prime is available.

## A criterion using powers of two and three

Let $h\ge1$, $U=2^{h+2}$, $D=3^{h+1}$ and $N=Uv=1+Dq$, with $v>U$ and $q>D$ primes. Suppose a prime $r>N/2$ has $0\le N-r<\min(D,q,2v)$ and is distinct from $3,q,v$. Then $3,q,v,r$ certify every triple at depth at most $h+1$.

For each of these primes consider the rows passing the first $h+1$ digits at 3, or the units digit at each of $q,v,r$. The four sets are pairwise disjoint by the following residue tests.

* At $D$ and $q$ the residues of a passing row must be $0,0,1$. A row passing both would have a positive part divisible by $Dq=N-1$, too large.
* A row passing $v$ has all three parts divisible by $v$. If it also passes at $P=D$ or $q$, its sum is at least $v(2P+1)>N$. Here $2D+1>U$ and $q>D$.
* A row passing $r$ has two positive parts with sum at most $N-r$. Passing $v$ would require both to be multiples of $v$. Passing $D$ or $q$ would require at least one to be a multiple of the corresponding modulus. The gap excludes each case.

Thus three rows can pass at most three of these four primes between them. The depth bound is attained by
$$
 E=(1,1,N-2),\quad A=(3^h q,2\cdot3^h q,1),\quad C=(2^h v,2^h v,2^{h+1}v).
$$
The endpoint limits common primes to $2,3,q,v$. Row $A$ is carry-free at $q$ and $C$ is carry-free at $v$. The nonzero first digits of $q,2q$ at 3 sum to three, so $A$ first carries at 3 in position $h+1$. Row $C$ first carries at 2 in that same position. At 2 the other two rows carry at the units digit, as do the other two at 3. The first common depth of this triple is exactly $h+1$.

A single prime can divide every coefficient only if the base-$p$ digit sum of $N$ is at most two. Otherwise one can distribute the digit units of $N$ among three nonempty parts and obtain a carry-free row. If every endpoint prime has digit sum at least three, each possible list of at most three primes can be defeated by one avoiding row per prime. Under this additional condition four is the minimum.

## Exact parameters

Take $h=512$. Set $q_0$ to the unique integer in $[D,D+U)$ satisfying $Dq_0\equiv-1\pmod U$, and put
$$
 q=q_0+233992U,\qquad N=1+Dq,\qquad v=N/U,\qquad r=N-273 .
$$
The three numbers $q,v,r$ are prime. These choices satisfy the gap inequalities above and all four endpoint digit sums exceed two. This proves both minima at the stated weight. The displayed triple has coefficient gcd **$2^{514}3^{167}$**.

| Prime | First digit carrying in $E,A,C$ | Factorial valuations in $E,A,C$ |
|---|---|---|
| $2$ | $1,1,513$ | $514,1071,1116$ |
| $3$ | $1,513,1$ | $513,167,1370$ |
| $q$ | $1$, none, $1$ | $1,0,2$ |
| $v$ | $1,1$, none | $1,1,0$ |

PARI/GP generated and validated elliptic-curve certificates of 28, 42 and 60 steps for $q,v,r$. FLINT independently proved all three numbers prime. The standard-library checker below checks those certificates using the portable Hasse-bound and Lucas verifier, compares direct digit addition with factorial floors, and checks the inequalities for all possible rows. The prime search sieved the two progressions $q_0+iU$ and $(1+Dq_0)/U+iD$ by primes up to 50,000 and tested the survivors. [Certificate data](depth_513_certificates.json.gz) · [Portable primality verifier](check_all_prime_depth_witness.py)

~~~sh
python3 check_sharp_depth_513.py
~~~

The [unconditional infinitude theorem](unbounded_absolute_carry_depth.md) gives triples requiring arbitrarily many digits. The four-prime criterion here additionally determines the exact depth for all triples at an individual weight.
