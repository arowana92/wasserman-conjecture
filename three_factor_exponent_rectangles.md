# Two unrestricted exponents in a three-factor predecessor

Let $D>1$ be a prime power, and let $q,r$ be distinct primes coprime to $D$. Consider
$$N=1+Dq^b r^c,\qquad b,c\ge1.$$
Let $\Lambda_p=(2\sin(\pi/(4p+2)))^{-1}$ and $\beta_p=\log_p\Lambda_p<1$.

**Theorem.** Take any arithmetic progression of $H$ positive exponents $c$, with any positive common difference. Outside $O_q(DH^{\beta_q})$ of those exponents, **every** $b\ge1$ gives a weight satisfying the full trinomial conjecture. The bound is independent of the prime $r$, even when $r$ varies with the progression.

If $b$ also ranges over an arbitrary arithmetic progression of length $L$, the number of exceptional pairs is
$$ O_{q,r}\big(D L^{\beta_r} H^{\beta_q}\big). \tag{1}$$
Both common differences and starting points are arbitrary. If $q>D$, one can replace $\beta_q$ in these conclusions by $\log_q((q+1)/2)$, and similarly for $r>D$.

For example, $(\beta_2,\beta_3)=(0.69424\ldots,0.73691\ldots)$ to the indicated approximate precision. In two progressions of length $H$ each, all but $o(H^2)$ exponent pairs $(b,c)$ work simultaneously for every prime power $D\le H^{0.28}$ coprime to 6 at $N=1+D\,2^b3^c$. Indeed, summing (1) over $D\le Y$ costs $O(Y^2)$, and $0.56+\beta_2+\beta_3<2$. In one progression of length $H$, all but $o(H)$ exponents $c$ likewise work for **all** outer binary exponents and all eligible $D\le H^{0.15}$, at any one odd prime $r$.

## The same cut in both bases

By the [endpoint theorem](all_weights_endpoint_reduction.md) and the [three-factor assignment](endpoint_factors.md#when-n-1-has-three-distinct-prime-factors), a counterexample has a positive row adding without a carry at every base of $D,q^b,r^c$. Each of these three full powers has its exceptional residue one in a different part. Otherwise a positive part would be divisible by $N-1$ and at most $N-2$. The part with residue one at $D$ is $x q^b r^c$, with $1\le x<D$ and $\gcd(x,D)=1$. There is also a unitary weight divisor of size at least three dividing $x$, though the counting proof does not need this.

Taking the binomial cut at this part preserves carry freedom. At base $q$ the complementary part is $q^b(D-x)r^c+1$. Above the bottom $b$ places the cut is exactly
$$ x r^c+(D-x)r^c \quad\text{in base }q. \tag{2}$$
At base $r$ the same argument gives $x q^b+(D-x)q^b$ without a carry. In particular, failure of (2) for every $x$ certifies the full conjecture for all $b$. This argument also certifies the weights with coincident prime bases by the earlier result for predecessors with at most two prime divisors.

## Extension to four predecessor primes

Both bounds, their sharper exponents when the base exceeds $D$, and the explicit progression construction below also hold when $D$ is twice an odd prime power, 12, 15, 21 or 33. Keep $q,r$ distinct and coprime to $D$. For the first three cases use the [two-power assignment criterion](two_prime_powers.md), and for 21 and 33 its [parity refinement](four_factor_squares.md#why-the-two-bases-must-be-assigned-to-the-same-row). In every case a counterexample assigns both $q,r$ to one row. Their zero positions overlap, giving a positive part $xq^br^c$ with $1\le x<D$. This gives exactly the two cuts in (2). The counting and construction use no coprimality between $x$ and $D$.

## How many exponent terms can pass a cut

Fix positive $s,t$, a prime $p\nmid s+t$, and an integer $A>1$ coprime to $p$. Write $j=v_p(st)$. The number of exponents $e$ in any progression of $H$ positive terms such that $sA^e,tA^e$ add without a carry in base $p$ is
$$ \ll_p (p/\Lambda_p)^j H^{\beta_p}. \tag{3}$$
For $j=0$ the stronger exponent $\log_p((p+1)/2)$ applies.

First count possible blocks of $n$ consecutive digits of a multiplier $z$, with all lower digits prescribed. If $s,t,s+t$ are units at $p$, varying the next input digit makes each output digit of $sz$ and $tz$ range over all $p$ digits. Their sum modulo $p$ also ranges over all digits, regardless of the two incoming multiplication carries. The difference between the sum of their two totals and the total of their reduced sum is $p(p-1)/2$. Each rejected input accounts for $p$, giving exactly $(p+1)/2$ choices at a carry-free step. Here $p$ is odd.

Otherwise, interchange $s,t$ to write $s=p^j u$ with $u,t$ units and $j\ge1$. Use the digits of $uz$ as coordinates for the input block. They determine it bijectively once the lower input digits are known. Leave the first $\min(j,n)$ digits free. At each later position the digit of $tz$ cannot exceed $p-1-d$, where $d$ was the true digit of $uz$ at distance $j$ below. The digit of $tz$ varies bijectively with the current digit of $uz$. Encode the latter by its zero-based rank among the $p-d$ possible values, in increasing order. The code is injective, and the encoded digit is at most the true digit. Hence encoded digits at distance $j$ have sum at most $p-1$.

The symmetric matrix with entries $\mathbf1_{a+b\le p-1}$ for $0\le a,b<p$ has operator norm $\Lambda_p$ by the [sine-vector calculation](inverse_digit_sieve.md#encoding-the-digits-of-the-cofactor), after reversing both index orders. Its all-ones Rayleigh quotient is $(p+1)/2$. Each chain of $v\ge1$ positions therefore has at most $p\Lambda_p^{v-1}$ words. Multiplying over residue classes modulo $j$ gives at most $(p/\Lambda_p)^j\Lambda_p^n$ blocks, including $n<j$ by the trivial count $p^n$.

For exponents $e=e_0+h i$, split the index $i$ into classes modulo $p-1$ at odd $p$, or modulo 2 at $p=2$. Within a class the scale advances by $B=A^{h(p-1)}$ or $A^{2h}$. Let $f=v_p(B-1)$, so $f\ge1$ and $f\ge2$ in binary. The order of $B$ modulo $p^{f+n}$ is $p^n$. This follows by expanding the $p$-th power of $1+v p^f$ with $p\nmid v$ and iterating the one-step valuation increase. In a class of length $H'$ choose $n=\lceil\log_p H'\rceil$. The bottom $f$ digits stay fixed and the next $n$ digits distinguish all $H'$ terms. The block bounds prove (3), including its uniformity in the common difference.

For $1\le x<D$ and $p\nmid D$, set $W_p(x)=(p/\Lambda_p)^{v_p(x(D-x))}$. For $j\ge1$ at most $2D/p^j$ choices have valuation at least $j$, since $p$ divides at most one of $x,D-x$. The inequality $\Lambda_p\ge(p+1)/2>\sqrt p$ makes the resulting geometric series converge, so $\sum_x W_p(x)^2\ll_p D$. For one free axis sum (3) by Cauchy–Schwarz. For two axes (2) confines possible pairs to a union of Cartesian products indexed by the same $x$. Apply (3) to the two factors and use $\sum_x W_q(x)W_r(x)\ll_{q,r}D$ by Cauchy–Schwarz again. This proves (1).

## Finding a whole family in a finite block

There is an explicit bound when $D$ is allowed to grow as fast as the length of the search. Fix $q,r$ and an arithmetic progression of positive exponents of $r$. Put $l=q-1$ if $q$ is odd and $l=2$ otherwise. For every $n\ge1$, **any $l(q^n-1)+1$ consecutive terms contain an exponent $c$ that works at all $b\ge1$ and all eligible prime powers $D\le q^n$ simultaneously**.

Starting at the first term, use the $q^n$ indices separated by $l$ steps. With $f$ as in the order proof, $r^c$ attains the largest lift of its fixed low residue modulo $W=q^{f+n}$. This lift is $W-u$ with $0<u<q^f$ because the low residue is a unit. For every $1\le x<D\le q^n$ the two residues in (2) are $W-ux$ and $W-u(D-x)$. Their sum is greater than $W$ since $uD<W$. This forces a carry at that layer for every cut. The same exponent recurs every $lq^n$ progression steps.

For example, $3^{1291}\equiv-5\pmod{8192}$ and $3^{2048}\equiv1\pmod{8192}$. Thus **every** $1+D\,2^b3^{1291+2048t}$ satisfies the conjecture for $b\ge1$, $t\ge0$, and every eligible prime power $D\le1638$. The extension from 1024 to 1638 uses $5D<8192$ directly. The verification below also tests the digit encoding at nonzero lower prefixes, compares both cuts with exact binomial coefficients, and checks the progression construction for irregular exponent steps.

~~~sh
python3 check_three_factor_exponents.py
~~~
