# Transferring a certificate between digit gaps

Fix a prime $p$ and positive integers $d,c$. For $Q=p^a>c$ write $N=dQ+c$. A trinomial row avoiding $p$ splits into $(Qx_i+b_i)_{i=1}^3$. The nonnegative $x_i$ sum to $d$ without a base-$p$ carry and the $b_i$ sum to $c$ without a carry. Each pair $(x_i,b_i)$ has a positive entry. These conditions describe the same shapes at every separated exponent.

Suppose no such shape has $d b_i=c x_i$ for all $i$. This holds in particular if $\gcd(d,c)\le2$. **If $p\mid c(c-1)$, every separated exponent $a_0$ at which Wasserman's conjecture holds belongs to an explicit arithmetic progression of exponents at which it holds**, restricting the progression to $p^a>c$. This includes every good positive exponent of $1+dp^a$, for every cofactor and every prime.

More precisely, put $Q_0=p^{a_0}>c$ and $n=dQ_0+c$. For an integer $B\ge1$ let
$$
 L_1=\text{the prime-to-}p\text{ part of }\operatorname{lcm}(1,\ldots,Bdc),
 \qquad L_2=\prod_{\substack{r\mid n(n-1)\\r\ {\rm prime},\ r\ne p}}
 r^{\lfloor\log_r n\rfloor}.
$$
Choose $a$ with $Q\equiv Q_0\bmod L_1L_2$, for instance any exponent in the class of $a_0$ modulo the multiplicative order of $p$ at this modulus (use period one at modulus one). Let $A,C$ be positive trinomial rows at $N$ and let $K_A,K_C$ denote their coefficients. Suppose $p\nmid K_A$, its shape is nonproportional as above, and
$$
                  \gcd(N,K_A,K_C)\le B.
$$
Then there are positive rows $A_0,C_0$ at $n$ such that
$$
 \gcd(n(n-1),K_{A_0},K_{C_0})\ \mid\ \gcd(N(N-1),K_A,K_C).
$$

## Bounding the second row before specializing it

Write $u=\gcd_i A_i$. For each $i$ we have
$$
 u\mid d A_i-x_i N=d b_i-cx_i.
$$
At least one of these integers is nonzero and has absolute value at most $dc$, so $u\le dc$. If all vanish, write the positive upper vector as $gz$ with $z$ primitive. Positivity follows since a vanishing upper coordinate would also have vanishing lower coordinate. The lower vector is an integer multiple of $z$ by primitivity, hence $\sum z_i\ge3$ divides both $d$ and $c$. This proves the simple sufficient hypothesis.

Write $m=N/\gcd_i C_i$. The factorial decrement identity $A_i K_A=N K_{A-e_i}$ gives $N/u\mid K_A$, and likewise $m\mid K_C$. Since $m\mid N$, we get
$$
             m\mid u\,\gcd(N,K_A,K_C).
$$
Indeed $m/\gcd(m,N/u)$ divides $u$ and $\gcd(m,N/u)$ divides the three-way gcd. In particular $3\le m\le Bdc$.

The prime-to-$p$ part of $m$ divides $L_1$. Its $p$-part divides $Q,Q_0$, since $m\mid N$ and $v_p(N)=v_p(c)<\min(a,a_0)$. It follows that $mL_2\mid Q-Q_0$, so $m\mid n$ and $N/m\equiv n/m\bmod L_2$. Define
$$
 A_{0,i}=Q_0 x_i+b_i,\qquad C_{0,i}=(n/m)\,\frac{C_i}{N/m}.
$$
Both rows are positive and sum to $n$. The shape of $A_0$ guarantees $p\nmid K_{A_0}$. At every $r\ne p$ dividing $n(n-1)$ all three coordinates in each row agree with their original counterparts modulo every $r^j\le n$. The contribution to Legendre's valuation at this layer is the floor of the sum of the three least residues divided by $r^j$. All layers are nonnegative and the layers larger than $n$ vanish at the reference weight. The valuations in both actual rows therefore dominate those in the two reference rows. The endpoint has the required multiplicities too, since the full $r$-power of $n(n-1)$ divides $L_2$. This proves the divisibility.

For the conjecture take $B=1$. A counterexample would contain the [endpoint row](all_weights_endpoint_reduction.md), whose coefficient is divisible by $p$ because $p\mid c(c-1)$. Name the other row avoiding $p$ as $A$. The division just proved transfers its pair with $C$ to a counterexample at $n$.

## Examples

* For $N=5\cdot2^a+7$ use $a_0=3$, giving $n=47$. All three-row choices work on the resulting progression. Quantitatively, take $B=46$. If a pair with one odd coefficient has weight gcd at most 46, the transfer supplies 47 in the common endpoint gcd, since every positive trinomial at weight 47 is divisible by 47. The modulus also forces $47\mid N$. Hence **every pair with one odd coefficient has weight gcd at least 47** on this progression.
* The same argument works whenever the reference weight is a prime $n$, with $B=n-1$. There is also a shorter sufficient modulus for this conclusion. If $n>Bdc$ and is prime, only $L_1$ is needed: $3\le m\le Bdc<n$ and $m\mid n$ are already incompatible. In this shorter-modulus version the conclusion is $\gcd(N,K_A,K_C)>B$.
* For every prime $p$ and cofactor $d$, each positive exponent $a_0$ with $1+dp^{a_0}\le10^9$ now anchors an infinite progression, by the [contiguous computation](full_range_1b.md). The earlier [small-weight transfer](small_weight_transfer.md) uses a different specialization at exponent zero and a relaxed carry test there.

Run the standard-library checker below. It verifies the scale divisibility on 159,519 pairs by exact integer coefficients. It also compares 2,700 pairs of valuation layers at reference exponents and lifted exponents using modular powering, including 2,070 where the second-row scale denominator is divisible by the test prime.

~~~sh
python3 check_digit_gap_transfer.py
~~~
