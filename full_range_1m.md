# The conjecture through one million

**Finite theorem.** For every $3\le N\le1{,}000{,}000$, any three positive trinomial coefficients of weight $N$ have a common prime divisor. Repetition of rows is allowed.

The [prime tables](RESEARCH.md) already settle $N\le5{,}000$. The [endpoint enumeration](endpoint_reduction.md) settles all selections through one million which contain $E=(1,1,N-2)$. Here we exclude every remaining selection at weights $5{,}001$ through one million.

## A single row near the largest prime

Suppose three rows have coefficient gcd one and none is $E$. Let $r$ be the largest prime at most $N$, and put $g=N-r$. Since $r>N/2$, a row avoiding $r$ has a part at least $r$. Write it as $R=(a,b,N-a-b)$ with $1\le a\le b$, $a+b\le g$. Bergman’s Proposition 10 excludes $3\le a+b\le10$, and $a+b=2$ would give $E$. This leaves 36,615,311 choices of $(N,a,b)$ across the entire additional range. [Bergman, Proposition 10](https://arxiv.org/abs/0806.0607)

For $0\le h\le5$ let $K_h$ be the lcm of 60 and all integers $w$ from 7 through $b$ satisfying
$$(a\bmod w)+(b\bmod w)\le h.$$
If a prime power $w>6$ divides $N-h$ and $R$ avoids its prime, its three residues sum to $h$. In particular the two displayed residues sum to at most $h$. For $w>b$ this would require $a+b\le h$, which is impossible. Thus $w\mid K_h$. All layers up to 6 are covered by 60.

Consequently the integer
$$F=\prod_{h=0}^5\left(\frac{N-h}{\gcd(N-h,K_h)}\right)^{6-h}$$
divides the product of the weighted factorials from the other two rows for any valid cutoffs of total 6 on each. Indeed each layer remaining above $K_h$ occurs for a unique offset and must be covered by one of those rows. [Weighted divisibility lemma](weighted_obstruction.md)

If an ordinary row has every part at least 2, use cutoffs $(2,2,2)$. For row $(x,y,z)$ the cost is at most $(xyz)^3\le N^9/3^9$. If the smallest part is 1, use $(1,2,3)$ in increasing order. The second part is at least 2 since $E$ is absent. Maximizing $y^3z^6$ with $y+z\le N$ gives cost at most $64N^9/3^9$. Two distinct rows cannot both contain a part at most 2. For the primes dividing $N$, the product of the full powers either such row could avoid divides the lcm of those two parts and hence is at most 2. The third row would have three positive multiples of at least $N/2$.

It follows that a counterexample requires
$$3^{18}F\le c N^{18},\qquad c=1\text{ if }a\le2,\quad c=64\text{ otherwise}.$$
Exactly **986** of the 36,615,311 candidates satisfy this necessary inequality.

## Certificates for the remaining rows

Each line of the [certificate file](gap_certificates_1m.txt) gives $N,a,b$ followed by either one or three primes. Eight rows have a universal prime, for which the base-$p$ digits of $N$ sum to at most two. Each of the other 978 rows has three primes dividing its coefficient, with the following property: no positive triple of weight $N$ avoids any two of them. Consequently two further rows cannot between them avoid all three.

The verifier exhausts all $p$-free triples by taking $x$ digitwise below or equal to $N$ in base $p$, then $y$ digitwise below or equal to $N-x$, with $z=N-x-y$. Restricting to $1\le x\le y\le z$ accounts for permutations. Legendre’s factorial formula tests whether a second prime divides each such coefficient. A sieve checks all certificate primes.

The enumeration that produced the file instead computed $K_h$ from gcds of $a-i,b-j$ over $i,j\ge0$, $i+j\le h$, checked pairwise incompatibility by allocating digits among three parts, and found the same 986 exceptional rows. The two descriptions of $K_h$ agree because an integer divides such a gcd exactly when the two least residues sum to at most $h$.

## Reproduction

Use the commands below, together with the [independent coefficient check through 5,000](RESEARCH.md).

```sh
python3 search_endpoint.py 3 1000000
python3 check_endpoint_search.py
python3 check_gap_range.py
```

The last command regenerates the entire integer inequality sweep, verifies each prime certificate, and checks that none is missing or unused.
