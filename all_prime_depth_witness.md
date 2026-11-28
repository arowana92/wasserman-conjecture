# A weight where the sharp universal carry depth is 130

There is a **728-digit weight** at which four primes certify Wasserman's conjecture for every triple, and four is the minimum. Measuring the needed number of digits over **all** primes gives a different invariant. Its exact value at this weight is **130**. Some triple has no common witness within the first 129 digits at any prime. Every triple has a common witness within 130 digits. The latter upper bound uses five specified primes.

This weight is an explicit exception to extending the [128-digit density-one theorem](bounded_carry_depth.md) to every weight. The conjecture itself holds here.

## Exact weight and three coefficients

Let $d=26077$, $H=365078^{129}$, and let $W_0$ be the unique integer in $[0,H)$ with the indicated four residues.

| $p$ | 2 | 7 | 89 | 293 |
|---|---:|---:|---:|---:|
| $W_0\bmod p^{129}$ | 2 | 108 | 9587 | 1654 |

Set $W=W_0+704588H$. Then
$$
 q=W/2,\qquad N=1+dW,\qquad v=N/7,\qquad r=(N-24)/151
$$
are exact integers. The three numbers $q,v,r$ are primes. Thus $N=7v$ and $N-1=2\cdot89\cdot293\cdot q$, with the six factors all distinct. The four local residues come from the [fixed-prime construction](fixed_prime_depth_obstruction.md). The index was found by sieving the two linear forms $q,v$ for small factors and then testing primality. The formulas above reproduce the weight without the search.

Consider the endpoint row $E=(1,1,N-2)$ and
$$
 A=(9345W,16408W,324W+1),\qquad C=(v,3v,3v).
$$
Their coefficient gcd is exactly **365078**. The table uses digit one for the units digit. A dash means the entire addition is carry-free.

| Prime | First carry in $A$ | First carry in $C$ | $v_p(K_A)$ | $v_p(K_C)$ |
|---:|---:|---:|---:|---:|
| 2 | 138 | 1 | 2267 | 2584 |
| 7 | 131 | 1 | 734 | 739 |
| 89 | 130 | 1 | 237 | 430 |
| 293 | 130 | 1 | 172 | 339 |
| $q$ | – | 1 | 0 | 2 |
| $v$ | 1 | – | 1 | 0 |

The six primes exhaust the divisors of $K_E=N(N-1)$, and $E$ carries at the units digit at each. The table therefore gives the lower bound 130 for the universal depth, with the endpoint included in the example.

## A four-prime certificate for arbitrary triples

Write $F_p$ for the rows avoiding $p$. In fact the four avoidance sets for **$v,q,293,r$** are pairwise disjoint. A positive row in $F_v$ is $v$ times a positive triple of sum seven. It cannot avoid $q$ or 293, since two of its parts would have to be zero modulo that prime and none of its normalized parts can be. A row in $F_q$ has two positive unshifted parts $kq,lq$ and a third part congruent to one modulo $q$. In particular $1\le k<52154$, since $N=52154q+1$ and $l>0$.

For any such unshifted part, avoidance of 293 requires $kq\bmod 293^j\le N\bmod293^j$ at every layer. Enumerating these 52,153 integers leaves 120 through layer 129, then $60,16,14,8,8,2,2,0$ through successive layers ending at 137. This proves disjointness of $F_q,F_{293}$.

It remains to test $r$. A row avoiding $r$ has parts $k_i r+b_i$ with nonnegative upper counts of sum 151 and low residues of sum 24, since $N=151r+24<r^2$. For $P=v,q,293^{129}$ put $n_0=N\bmod P$. If the row also passes the corresponding layers in that base, its part residues $z_i\bmod P$ sum to $n_0$. Each $0\le z_i\le n_0$ satisfies
$$
 151z_i\equiv k_i(n_0-24)+151b_i\pmod P .
$$
This is an integer equality because $2\cdot151(n_0+24)<P$. In the three cases $n_0$ is $0,1,43131359$ respectively, with $\gcd(151,n_0-24)=1$. Thus $151\mid k_i$ for all $i$, and two positive parts are at most 24. They exceed the permitted units residue sum of zero or one in the other base. The same argument works at $P=89^{129}$ with $n_0=1+d\cdot9587$.

Disjointness means that three rows can avoid at most three of the four primes between them. For minimality, a universal prime would have to divide $K_E$. Each of its six prime divisors has base-$p$ digit sum of $N$ at least three. Distributing the digit units into three nonempty boxes gives an avoiding row. Any list of at most three primes can therefore be defeated by choosing one avoiding row per prime.

The four displayed primes alone require depth 137. The row $(9669q,42485q,1)$ avoids $q$ and its first carry at 293 is at digit 137. Together with any row avoiding $v$ and any row avoiding $r$ it defeats that list through digit 136. The pairwise exclusions above use at most 137 digits throughout.

## Five primes attain the sharp depth

Use $v,q,r,89,293$ through 130 digits. The sets of rows passing $v,q,r$ are pairwise disjoint, so a triple defeating the list would need one row for each. The row passing $v$ cannot pass 89 or 293 even at their units digits. The preceding residue argument excludes both for the row passing $r$. Therefore the row passing $q$ would have to pass **both** 89 and 293. Enumerating its possible unshifted part $kq$ as above gives **no** $k$ passing those two bases through 130 digits. This proves the matching upper bound for all triples, independently of the general endpoint reduction.

There is a second check for pairs with the endpoint. In a failure the row passing $v$ cannot take $7,q,89,293$. Its companion must have all parts divisible by 7. The normalized parts sum to $v$ and must distribute the exceptional positions for $q,89,293$ among all three parts, because $89\cdot293\cdot q>v$. CRT gives just the row $A$ up to permutation. Its carries at 89 and 293 suffice.

## Primality and arithmetic verification

Run the following standard-library command. It verifies the three stored elliptic-curve primality certificates, all factorial valuations, the CRT row enumeration, and the four- and five-prime checks.

~~~sh
python3 check_all_prime_depth_witness.py
~~~

PARI/GP generated and validated certificates of 100, 93 and 94 elliptic steps for $q,v,r$. FLINT independently returned a rigorous primality result for each. The portable checker reads the [certificate data](all_prime_depth_certificates.json.gz). At a step on $n$ it checks nonsingularity, computes an affine point $Q$ over the residue ring, and verifies that $Q$ is nonidentity and $q'Q$ is the identity for the next certified prime $q'>(n^{1/4}+1)^2$. All denominators used in addition must be invertible. Thus every prime divisor $p$ of $n$ has a reduction of $Q$ of order $q'$. Hasse's bound excludes $p\le\sqrt n$. The leaves use Lucas primality witnesses with completely factored $n-1$, ending in trial division. Each factorial-floor contribution is also compared with a direct addition of the digits, including its incoming carry. [PARI certificate format](https://pari.math.u-bordeaux.fr/dochtml/html/Arithmetic_functions.html)

## The simultaneous-prime question

The [unconditional infinitude theorem](unbounded_absolute_carry_depth.md) gives unbounded absolute depth using growing endpoint powers. Here the endpoint remains squarefree, so the delay also counts digits beyond the full endpoint powers.

For any chosen $h\ge8$ one can instead impose the four template residues through $p^h$ and refine the residue at 89 to $W\equiv9587+89^h\bmod 89^{h+1}$. If $q=W/2$ and $v=(1+dW)/7$ are both primes greater than 52154, their bases are genuinely avoided by $A,C$ respectively. All four other endpoint primes clear $A$ through depth $h$. At 89 the next column of $A$ carries, because the seed parts at 9587 lie below $89^h$ and the three slope residues include two nonzero complementary residues. The $v$-row argument forces $A$ as above, so the endpoint pair depth at the weight is exactly $h+1$. The all-weight endpoint reduction then gives the ordinary conjecture at that weight.

The two linear forms for $q,v$ on this progression are admissible. At the four modulus primes both are units, including $v\equiv6\bmod7$. At each remaining prime, which is odd, at most two residue classes are forbidden. Infinitely many simultaneous prime values for each $h$ would prove unbounded depth over all weights. That last infinitude is a prime-tuples hypothesis. The 728-digit example and its sharp depth use the explicit prime certificates above.

One can preserve the **exact depth for all triples and the minimum list of four** in this conditional statement. For $h\ge16$ keep the binary and base-seven residues through $p^h$, and use
$$
 \begin{aligned}
 W&\equiv9587+75\cdot89^h&&\pmod{89^{h+1}},\\
 W&\equiv1654+2726850956040050170\cdot293^h&&\pmod{293^{h+8}},\\
 1+dW&\equiv24+151&&\pmod{151^2}.
 \end{aligned}
$$
These are the upper digits of the numerical witness at the two tested bases, now placed after an arbitrary prefix. The last condition makes $r=(N-24)/151$ integral and a unit at 151. Whenever $q,v,r$ are all sufficiently large primes on this progression, the least universal depth is $h+1$, and the four primes $v,q,293,r$ work by depth $h+8$.

Here is why the finite cut computations transfer. At 293 the low part of $q$ is 827. Every $827k$ in the enumeration and the corresponding low weight fit in $h$ digits. Thus the higher positions of $kq$ are exactly the digits from multiplying $k$ by the fixed upper block of $q$ modulo the required power. At 89 an odd $k$ already fails the low cut test: its residue is $(89^h+9587k)/2$, which lies strictly between $1+9587d$ and $89^h$. For even $k$ use $kq=(k/2)W$, whose low product again fits in $h$ digits. The tests through depth $h$ on the small products simply pad fixed integers by zeros. Consequently the surviving cut sets at $h+j$ are the same as at $129+j$ in the example, at all used $j$. The low residues for the comparison with $r$ are unchanged. The proof with five primes therefore applies at $h+1$ and the four-prime proof at $h+8$. The rows $E,A,C$ pass all primes between them through $h$. The digit sums excluding a universal prime are already at least three in the fixed low prefixes at the four small endpoint primes.

The three linear forms on the progression are admissible. At $2,7,89,293$ the added form $r$ is also a unit. At 151 the weight is 24, $q$ and $v$ are units, and $r$ is one modulo 151. At 3 the forbidden residues for $v$ and $r$ coincide since 24 is zero there. At every other prime at most three classes are forbidden. The forms have positive slopes and distinct rational roots. Thus the prime-tuples conjecture for three linear forms gives infinitely many such weights for each $h$; no infinitude assertion is needed for the computed weight.
