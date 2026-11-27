# Separating powers of one prime through binary exponent 35

Let $N(k,Q)=1+5\cdot2^k Q$, where $Q$ is a positive prime power. Write a pair $(k,e)$ as exceptional for a prime $\ell$ when the full three-row conjecture fails at $N(k,\ell^e)$.

**For any prime power $R$, among all $1\le k\le35$ and $1\le a\le10^{18}$, at most one pair can fail at $N(k,R^a)$.** More precisely, if $(k,e)$ and $(j,f)$ are distinct exceptional pairs for the same prime, with $k,j\le35$, then
$$
                      \max(e,f)>10^{18}\gcd(e,f).
$$
There is an infinite increasing sequence of odd multipliers $a_i=O(i^{100/79})$, starting with all the odd numbers through $10^{18}$, such that for **every** prime power $R$ at most one pair $(k,i)$ with $1\le k\le35$ can fail at $N(k,R^{a_i})$.

Two finite ranges have **no exceptions**. With $1\le k\le35$ one may take any prime power $Q<2^{384}$, or $Q=\ell^e$ with $\ell\le1{,}000{,}000$ prime and $1\le e\le2^{128}$. At $k=31$ the first limit improves to $2^{1152}$.

## The surviving binary inputs

The [previous exhaustive reduction](factor_five_binary_exponents.md) proves the full assertion for every other binary exponent through 35. For each of the three remaining exponents it gives a unique low residue $q$ of $Q$ modulo $2^{k+1}$ and a constant $c$. In a counterexample $Q\equiv3\bmod8$, and the coefficient of the row
$$
     (3P Q,\ cQ,\ (2P-c)Q+1),\qquad P=2^k,
$$
must avoid 2, 3, 5 and the prime base of $Q$. In addition $PQ\equiv7\bmod15$ and $3\parallel1+5PQ$.

The binary input graph starts after the prescribed $k+1$ digits of $Q$. Its transitions read one further digit and retain the multiplication quotients for the total and parts, rejecting a carry. A state is productive when some sequence of digits followed by high zeros finishes a valid binary input. Discard all other states by reverse traversal from the stopping states. A stopping state can always read another zero. A prefix of any possible $Q$, including one padded with high zeros, thus belongs to a path in this graph. An independent graph keeps the pending bits of $Q$ and the upper quotient of $(cQ-1)/2^{k+1}$ and tests a binary subset condition. [Derivations of both graphs](factor_five_binary_exponents.md#searching-all-higher-digits)

| $k$ | $q$ | $c$ | Distinct prefixes at 124 bits | At 176 bits | Exponent $\gamma$ |
|---:|---:|---:|---:|---:|---:|
| 31 | 741526579 | 4192014075 | 3 | 6 | .019 |
| 32 | 5727322795 | 6285315 | 324 | 4772 | .074 |
| 34 | 429493347 | 34312377675 | 8 | 33 | .105 |

Starting with one path at every productive state, sum the path counts at its productive successors for 2,048 steps. For each of the three graphs the maximum $M$ obeys $M^{1000}<2^{2048(1000\gamma)}$, checked as an integer inequality. Each shorter path can be extended since every productive state has an outgoing productive edge. Splitting a long path into blocks bounds the number of distinct $n$-digit prefixes by $O(2^{\gamma n})$. The constants depend only on these three fixed graphs.

Consequently in any interval of $H$ consecutive exponents $e$, uniformly in $\ell$, the number of exceptional pairs at the corresponding $k$ is $O(H^\gamma)$. An exceptional $\ell$ is 3 modulo eight and has order $2^{n-2}$ modulo $2^n$, so the exponents have distinct residues when $H\le2^{n-2}<2H$. For a fixed odd exponent $e$ the same bound applies to prime bases in any $H$ consecutive integers, since odd powers permute the units modulo a power of two. Even exponents are always excluded.

## The exponent ratio test

Suppose $Q_1=R^a$, $Q_2=R^b$ are exceptional. On dividing the original exponents by their gcd we may take $a,b$ positive, odd and coprime, with $R$ a prime power congruent to 3 modulo eight. For $s$ bits write $x,y$ for the discrete logarithms of $Q_1,Q_2$ to base 3 modulo $2^s$. They are odd classes modulo $2^{s-2}$ and satisfy $bx=ay\bmod 2^{s-2}$. The logarithm exists by the order calculation in the [unrestricted binary-exponent theorem](five_power_separation.md#the-test-in-the-odd-prime-power-alone).

Put $B=10^{18}$ and $D=2^{122}>2B^2$. The 335 prefixes at 124 bits have distinct logarithms, including across the three values of $k$. If $x=y$ at this precision and $a,b\le B$, invertibility gives $a=b$. The prefixes also force the same $k$, which gives just one pair. For distinct $x,y$, put $h=y/x\bmod D$ and compute a shortest nonzero vector in the integer lattice $b\equiv ha\bmod D$. Its determinant is $D$. If a vector with $1\le a,b\le B$ exists, that vector and a shortest vector both have norm at most $\sqrt2B$. Two independent vectors in the lattice have parallelogram area at least $D>2B^2$, so these two must be parallel. Since $a,b$ are coprime, a shortest lattice vector is the candidate itself up to sign.

For completeness the verifier starts with $(1,h),(0,D)$, repeatedly putting the shorter vector first and subtracting the nearest integer multiple of it from the longer. Each reduction outside the stopping condition strictly reduces a squared length. At termination the basis $u,v$ satisfies $|u|\le|v|$ and $|u\cdot v|\le |u|^2/2$. Every $nu+mv$ with $|m|\ge2$ has perpendicular component at least $\sqrt3|u|$. For $|m|=1$ its length is at least $|v|$, by the nearest-integer condition. This proves the shortest-vector claim.

Of the $\binom{335}{2}$ pairs of coarse prefixes, 3,947 admit a positive primitive odd ratio within the bounds. At 176 bits **none** of their descendant prefixes satisfy $bx=ay\bmod 2^{174}$. The independent checker takes coarse logarithms to base 11 and recovers candidates from continued fractions. Indeed if $0<b=ha-Dt\le B$, then $|h/D-t/a|<1/(2a^2)$ and $\gcd(t,a)=1$, so $t/a$ is a convergent. At the final precision it compares the actual powers $Q_1^b,Q_2^a\bmod 2^{176}$ directly. It recovers the same counts.

## One infinite sequence, independent of the prime base

At $n$ bits let $S$ be the union of the possible odd logarithms across the three graphs and set $D=2^{n-2}$. Then $|S|=O(D^{.105})$, and the set of possible ratios $S/S$ has size $O(D^{.210})$. This estimate is uniform in $n$. If two distinct positive odd exponent multipliers $a,b$ have $b/a\bmod D$ outside $S/S$ then they cannot both work for an exceptional pair over the same prime power $R$. For any one multiplier the three disjoint low prefixes already exclude two different exceptional binary exponents.

Start with the stated finite initial sequence. Process the remaining positive odd integers in increasing dyadic blocks $[X,2X)$ with $X$ a power of two. Use $D=2X$ and insert $b$ when $b/a\bmod D$ lies outside $S/S$ for every earlier inserted $a$. This is a finite calculation at each step. In a full block past the initial seed, each earlier or newly inserted $a$ can account for at most $|S/S|$ rejected candidates because all candidates have distinct residues modulo $D$. If $T$ is the total sequence length up to the end of that block, all its $X/2$ odd integers have been accounted for by at most $T(1+|S/S|)$ choices. Thus $T\gg X^{.790}$. Applying this at successive powers of two gives $a_i=O(i^{100/79})$. Any pair of selected exponents has a finite exclusion at the stage it was formed, or belongs to the already checked initial sequence. This proves the infinite assertion.

## Two verified finite ranges

For a fixed prime $\ell\equiv3\bmod8$, $\ell>5$, each possible 130-bit prefix of $Q$ determines a single odd $e$ in $[1,2^{128}]$ by dividing its logarithm by that of $\ell$ modulo $2^{128}$. First require the slopes $3P,c,2P-c$ to add without a carry at $\ell$. This is necessary since multiplication by $Q=\ell^e$ shifts the slopes and the extra 1 lies below the shift. Next test $PQ\equiv7\bmod15$, $9\nmid1+5PQ$, then the base-3 and base-5 factorial layers. Every candidate is excluded by 32 places at one of these two bases.

| $k$ | Prime bases through one million passing the slope test within the required prime class | Exponents checked | After the modulo-45 conditions |
|---:|---:|---:|---:|
| 31 | 2859 | 8577 | 426 |
| 32 | 2785 | 1275530 | 122211 |
| 34 | 2468 | 27148 | 1532 |

For the bound on $Q$ itself one can allow all odd integers. Count every accepted padded binary word below the specified power of two, then apply the odd-base and residue tests. The counts at $k=31,32,34$ are respectively 24,845 words below $2^{1152}$, 5,321,850 below $2^{384}$ and 430,247 below $2^{384}$. After the modulo-45 conditions and the first seven ternary and three quinary places the remaining counts are 0, 227 and 5. Full factorial floors at 3 and 5 eliminate all 232. The first verifier joins low and high binary words at a middle state. The second traverses forced segments of the independent graph and checks base digits directly, starting with a different low-layer sieve.

~~~sh
python3 check_five_power_separation.py
python3 check_five_power_separation_cross.py
~~~
