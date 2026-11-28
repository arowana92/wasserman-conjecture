# Logarithmic blocks requiring different large common primes

Fix $0<a<\theta$. For every sufficiently large $X$, put $y=\lfloor a\log X\rfloor$. The interval $[X,X+X^\theta]$ contains **$y-1$ consecutive weights** at which the explicit triples below require at least $\lceil(y-1)/2\rceil$ different common primes, all greater than $y$. For $X=2^m$ and $a<4/7$, the block can be chosen entirely within weights where the full three-row conjecture holds. The number $4/7$ comes from the project's [endpoint interval theorem](all_weights_endpoint_reduction.md#consequences).

An explicit example gives **1,000 consecutive weights above $2^{16384}$** at which the conjecture holds and any one list certifying the whole block must include **at least 500 distinct primes greater than 1,009**. This concerns a single list used across weights. A different short list can work at each individual weight.

The [elementary block construction](consecutive_certificate_lower_bound.md) gives the exact witness gcd at $N=B+h$, $2\le h\le T$,
$$
                 \gcd(K_1,K_2,K_3)=\frac{N(N-1)}{h(h-1)}.                 \tag{1}
$$
We retain (1) and reduce the log of the required modulus to $\vartheta(y)+o(T)$ for $T\to\infty$ and $y\ge T$. In the asymptotic claim take $T=y$. The method also permits smaller block lengths.

## Large primes on opposite sides of two consecutive integers

Let $D=T/(\log T)^4$. Call a prime $p>D$ paired if some $2\le h\le T$ has two distinct primes exceeding $D$ dividing $h(h-1)$, one of them $p$. The number of paired primes is $o(T/\log T)$.

Here are the sieve details. For large $T$, $D^2>T$. Such a pair gives $mp-nq=\pm1$ with $m,n\le(\log T)^4$, $\gcd(m,n)=1$, and $mp,nq\le T$. For fixed $m,n$ and sign, parameterize the integer solutions by $p=p_0+nt$, $q=q_0+mt$ on an integer interval of length at most $T/(mn)+1$. For every prime $\ell\nmid2mn$ the two forms have different roots modulo $\ell$, by the equation with determinant $\pm1$. For squarefree $d$ supported on these primes, there are $2^{\omega(d)}$ simultaneous bad residue classes modulo $d$ for divisibility of the product by every prime of $d$.

Put $z=\exp(\log T/(1000\log\log T))$ and let $k=\lceil20\log\log T\rceil$. The target primes exceed $z$. Apply the upper Bonferroni bound through order $2k$ to the primes $\ell\le z$, $\ell\nmid2mn$. Each squarefree modulus contributes its predicted interval density with error at most $2^{\omega(d)}$. Summing the errors costs at most $(2z+2)^{2k+1}=O(T^{1/10})$. The truncated density exceeds the full product by at most the elementary symmetric sum at order $2k+1$, hence by at most
$$
 \frac{(\sum_{\ell\le z}2/\ell)^{2k+1}}{(2k+1)!}=O((\log T)^{-10}).
$$
This follows from Mertens' estimate and $j!\ge(j/e)^j$. The full product is
$$
 \prod_{\substack{2<\ell\le z\\\ell\nmid mn}}(1-2/\ell)
           \ll\frac{(\log\log T)^4}{(\log z)^2}.
$$
Indeed $1-2/\ell\le(1-1/\ell)^2$, and restoring the omitted factors costs at most $\prod_{2<\ell\mid mn}(1-1/\ell)^{-4}\ll(\log\log T)^4$, since $mn\le(\log T)^8$. These uses of Mertens' product formula are uniform in $m,n$. Summing the interval lengths costs $O(T(\log\log T)^2+(\log T)^8)$. Thus the number of pairs is $O(T(\log\log T)^8/(\log T)^2)+O(T^{1/10}(\log T)^8)$, as claimed. [Mertens' product formula, Theorem 1.1](https://math.colgate.edu/~integers/graham17/graham17.pdf)

## Choosing the blocks of digits

Put $F=T/\sqrt{\log T}$. At each offset $h$ assign active primes $p\mid h(h-1)$ to two groups according as $p\le D$ or $p>D$. Put $Q_p$ equal to the least power of $p$ exceeding $T$.

* For $p\le D$ choose the first part $A_1=\prod Q_p$ over the small group at $h$. At each such prime prescribe $B\equiv-Q_p\bmod J_p$, where $J_p$ is a power of $p$ greater than $Q_p$ and every $A_1$ of a group containing $p$. Take the least such power. As in the original construction $\log J_p=O((\log T)^2)$, so the cost for all these primes is $O(D(\log T)^2)=o(T)$.
* For $D<p\le F$, and for all paired primes, prescribe $B\equiv-p^2\bmod p^5$. The second group has at most two members. Whenever it has two, both are paired. Its part $A_2=\prod p^2$ is strictly below $p^5$ at each of its primes since $p>D>T^{2/3}$ eventually. For a singleton with this prescription take $A_2=p^2$. The combined log cost is $O(\vartheta(F))+o(T/\log T)\,O(\log T)=o(T)$.
* At every remaining prime through $y$, prescribe only $B\equiv0\bmod p$. For an unpaired $F<p\le T$, also require that the second base-$p$ digit of $B$ be at most $p-1-\lfloor T/p\rfloor$. The next paragraph shows how to satisfy these inequalities without putting $p^2$ into the modulus.

The product $L$ of the CRT moduli has $\log L=\vartheta(y)+o(T)$. Starting from the first CRT solution at least $X$, try $\lceil T^3\rceil$ consecutive translates by $L$. For an unpaired $F<p\le T$ with just the first digit prescribed, $B/p\bmod p$ runs once through all residues per $p$ translates. At most $\lfloor T/p\rfloor$ residues fail. The sum of the excluded proportions is at most $T\sum_{p>F}p^{-2}\ll T/(F\log F)=o(1)$ by partial summation and the prime number theorem. The total end errors are at most $\sum_{F<p\le T}\lfloor T/p\rfloor\le T^2$. Some translate satisfies every inequality.

If such an unpaired prime is active at $h$, it is the sole member of the second group. Set $A_2=p^{v_p(B)}$. This is smaller than $B$ because $B$ also has other prime divisors. Splitting $B$ into $A_2$ and $B-A_2$ is carry-free. Adding $h$ is also carry-free since $h<p^2$, the units digit of $B$ is zero, and its second digit leaves room for $\lfloor h/p\rfloor$. If $p$ divides $h$ or $h-1$, the corresponding tens digit of $B+h$ or $B+h-1$ is between 1 and $p-1$, since the active multiple of $p$ is positive. This keeps the predecessor and weight valuations exactly as at $h-1$ and $h$.

Use 1 for an empty group. For all the prescribed higher blocks the [borrow-free subtraction argument](consecutive_certificate_lower_bound.md#the-exact-gcd-at-each-position) applies directly. The two rows $(h,A_i,B-A_i)$ and the endpoint now have gcd (1), with all its primes greater than $y$. Nonadjacent offsets have coprime witness gcds, and adjacent offsets share $(B+h)/h$. Therefore exactly $\lceil(T-1)/2\rceil$ primes are needed to hit these witnesses. For $T=y=\lfloor a\log X\rfloor$, the search displacement is $O(T^3 L)=X^{a+o(1)}$. This proves the interval and power-of-two assertions.

## Reproducible example

The [checker](check_sparse_certificate_modulus.py) computes the paired primes directly and chooses the exact power $J_p$ sufficient for the finite input, also on prescribed large-prime blocks. It checks independently that both cuts $\binom{N}{h}$ and $\binom{B}{A_i}$ are units at their assigned bases, that the integer in (1) has no prime factor through the cutoff, and that every weight meets the full-conjecture interval inequality.

For length 1,000 and cutoff 1,009 the modulus has 6,214 bits and the first translate above $2^{16384}$ works. Its offset has 6,213 bits. There are 18 unpaired primes using the free second digit, and 38 offsets use their singleton splits. For these finite sizes the checker takes the small threshold at least $\lfloor\sqrt T\rfloor+1$ to keep the two-member bound valid. This has no effect on the asymptotics. [Lean verification of the common-divisor identity](PairDivisor.lean)

~~~sh
python3 check_sparse_certificate_modulus.py
~~~
