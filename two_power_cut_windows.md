# Two prime powers and the smallest part

Let $N\ge3$ and $M=N-1$. Two positive trinomial rows $A,B$ with no common prime factor of $NM$ would give a counterexample with $(1,1,N-2)$. Conversely the [all-weight endpoint reduction](all_weights_endpoint_reduction.md) shows that ruling out such a pair proves the full three-row conjecture at $N$.

## A single predecessor power bounds all six parts

Let $Q=q^j>1$ divide $M$. In such a pair every part exceeds $Q$. At most one of the six parts is at most $2Q$. If this part occurs it equals $Q+1$ and $\gcd(N,Q+1)$ contains a unitary divisor of $N$ of size at least three. A unitary divisor is coprime to its complementary divisor.

To prove this, assign every full prime power of $N$ to a row avoiding its base, with products $U,V$ respectively. Then $UV=N$, the parts of the first row are multiples of $U$ and those of the second are multiples of $V$. Both $U,V\ge3$, since the normalized opposite row has three positive parts summing to the scale. Name the rows so the first avoids $q$. Two of its parts are multiples of $UQ$ and the other is $1\bmod Q$, by carry freedom at the layer $Q$. The latter is a multiple of $U$ and hence at least $Q+1$. Dividing by $U$ gives $V\ge2Q+1$, so each part of the opposite row is at least $2Q+1$. The sole possible value in $[Q+1,2Q]$ is $Q+1$ in the first row, forcing $U\mid\gcd(N,Q+1)$.

There is a quantitative version for any two positive rows, even when there is a common endpoint prime. Write $D=\gcd(N,K_A,K_B)$ and $J=\gcd(M,K_A,K_B)$ for their trinomial coefficients $K_A,K_B$. Every part $z\ge2$ satisfies $DJz\ge Q$. Retain the power $Q'=Q/\gcd(Q,J)$ and use the [retained-power assignment](endpoint_pair_bound.md). For $Q'>1$ the assigned row has parts $0,0,1\bmod Q'$, so any of its parts above one is at least $Q'$. The other row's parts are multiples of its scale $V$, with $DV\ge2Q'+1$ by dividing the assigned row through by $U$. Either alternative gives the inequality, and it also holds when $Q'=1$.

**Weight test.** Let $P=p^i>1$, where $p\mid NM$. Write $N=dP+h$ with $0\le h<P$. Suppose the base-$p$ digit sum of $d$ is at most two. The full conjecture holds if $h\le Q$. It also holds for $h\le2Q$ if $\gcd(N,Q+1)$ has no unitary divisor of $N$ above two. In a row avoiding $p$, the upper blocks add without a carry to $d$. At most two are nonzero, so at least one positive part is at most $h$. This contradicts the preceding bounds. The two primes may coincide. The units in $d$ can be arbitrarily far apart.

In particular, write $N=eQ+1$. For the extended interval the gcd to test is just $\gcd(e-1,Q+1)$. For fixed $e\ge2$, any odd prime $p\nmid e$, and any fixed positive multiple $L$ of $p$ coprime to $e$, the first interval proves the full conjecture at $\gg_{e,L}X/\log X$ weights $N\le X$ **with $L\mid N$**. Choose $P=p^i$ in $[X/(3p),X/3]$. Take primes $Q\equiv-e^{-1}\bmod L$ in the interval where
$$P/(2e)\le eQ+1-2P\le P/e.$$
Then $h<P$, $h<Q$, and the leading quotient is two. The interval for $Q$ has length $P/(2e^2)$ with endpoints proportional to $P$, so the stated count follows from the prime number theorem in a fixed progression. Both $N$ and $N-1$ may thus be required to have arbitrarily many distinct prime factors, by first choosing $e$ and $L$. [Kedlaya, Theorem 4.12](https://kskedlaya.org/ant/chap-primes-in-ap.html)

For example take $p=3$, $i=80$, $e=11\cdot13\cdot19$ and
$$Q=108817983047835986093107686275454247.$$
The resulting certified weight is $295658459940970374214973583610409189100$. Its offset above $2\cdot3^{80}$ is $40801112278527582807163197642593898$. Here $3\cdot5\cdot7\cdot17\mid N$ and $Q$ is prime. The offset lies strictly between $P/(2e)$ and $P/e$ and its square exceeds $N$.

## Leading quotients at two arbitrary endpoint bases

Let $P=p^i,Q=q^j$ have distinct prime bases dividing $NM$. They need not themselves divide $N$ or $M$. Write
$$N=dP+h=eQ+k,\qquad 0\le h<P,\quad 0\le k<Q,\quad d,e\ge1.$$
Suppose
$$\gcd(d,e)=1,\qquad de(h+k)<N,\qquad dehk<N.$$
Then any putative counterexample forces one of the two nonendpoint rows to be an **exact** weight proportion. Its three parts are $N/E$ times three positive integers of sum $E$, with $E\mid\gcd(d,N)$ or $E\mid\gcd(e,N)$. In particular the conjecture holds if both of these gcds are at most two. When $h,k\ge2$, the product inequality already implies the sum inequality.

A row avoiding $p$ has entries $a_i=j_iP+r_i$, where the nonnegative $j_i$ sum to $d$ and the $r_i$ sum to $h$. Thus $x_i=da_i-j_iN=dr_i-j_i h$ has absolute value at most $dh$. Similarly a row avoiding $q$ has leading blocks $l_i$, remainders of sum $k$, and $y_i=eb_i-l_iN$ of absolute value at most $ek$.

First suppose the same row avoids both bases. For every part we have
$$|(ej_i-dl_i)N|\le e|x_i|+d|y_i|\le de(h+k)<N.$$
Thus $ej_i=dl_i$, forcing all the leading blocks in each base into the same single position by coprimality. The other two positive parts have sum at most $\min(h,k)$, whose square is at most $hk<N$. This contradicts the [pair-sum bound](endpoint_pair_bound.md), since a pair sum $s$ in these rows has $s^7>6^6 M^4>N^{7/2}$ for $N\ge3$.

The two avoiding rows must therefore be different. From their scales $UV=N$ we have $N\mid a_i b_m$ for every $i,m$. Consequently $N\mid x_i y_m$. The product bound forces all these products to vanish. Either all $x_i=0$ or all $y_m=0$. In the first case normalize $(j_1,j_2,j_3)$ by its gcd. Its normalized sum $E$ divides $d$, and $E\mid N$ because the normalized entries are coprime and the actual entries are integral. All three are positive, so $E\ge3$. The second case gives the other alternative.

## A large predecessor power at an arbitrary leading quotient

Keep $Q\mid M$. For $N=dP+h$, $1\le h<P=p^i$ and $p\mid NM$, suppose $dh<Q$. In any counterexample, a nonendpoint row avoiding $p$ must have one of the forms
$$ A=(N/E)z,\qquad A=(M/E)z+e_j, $$
where $z$ is a primitive positive triple of sum $E$. In the first case $E\mid\gcd(d,h)$, and in the second $E\mid\gcd(d,h-1)$. Thus if both gcds are at most two, the conjecture holds. The two alternatives and their digit restrictions were previously established when $(dh)^5<27M^2$. The condition here can include offsets proportional to $N$. [Determinant criterion](projected_rows.md#two-blocks-of-digits-in-any-prime-base)

Indeed write $A_i=j_iP+r_i$ with the upper and lower blocks summing to $d,h$, and put $x_i=dA_i-j_iN$. If the opposite row avoids $q$, the scale dividing $N$ and every $A_i$ is at least $2Q+1$ by the single-power argument. It divides all $x_i$, whose magnitudes are at most $dh$, so all vanish. Otherwise $A$ avoids $q$, and $A_i\equiv\delta_i\bmod Q$ for a unit vector $\delta$. Now
$$d(A_i-\delta_i)-j_iM=d(r_i-\delta_i)-j_i(h-1)$$
is divisible by $Q$ and has magnitude at most $dh$. Hence it vanishes. In either case $j_i>0$ because every part exceeds $Q$ (and in the shifted case exceeds one). Divide the leading triple by its gcd to obtain $z$. Integrality gives the two divisibilities. This also preserves all further restrictions proved for those exact forms in the determinant criterion.

With **two** predecessor powers $Q,R$ of distinct bases, the same conclusions hold if
$$\min(Q,R)>d,\qquad dh<2\min(Q,R),\qquad |Q-R|\ne d.$$
If the opposite row avoids either base, its scale argument still forces all $x_i$ to vanish. Otherwise write the two residue vectors as $\delta,\epsilon$. For each position,
$$dA_i-j_i M=d\delta_i+m_iQ=d\epsilon_i+n_iR,\qquad m_i,n_i\in\{-1,0,1\}.$$
The bounds follow from the preceding display. If one multiplier is zero both must be zero since $Q,R>d$. If both are nonzero, opposite signs give a difference greater than $d$ and equal signs require $|Q-R|=d$ or $Q=R$. This proves the claim.

Here is an infinite family with an arbitrary binary prefix $d\ge3$. Fix $e\ge2$ coprime to $d$ and $L\ge1$ coprime to $de$. For $P=2^i$, take $N=eQ+1$ with $Q$ prime in the interval
$$P/(2e)\le h=N-dP\le P/e.$$
Require $eQ\equiv1\bmod d$ and $eQ\equiv-1\bmod L$. These conditions specify a reduced residue class modulo $dL$. For large $P$ we have $dh<Q$ and the two gcds are at most two, since $N\equiv2\bmod d$. Thus all these weights satisfy the conjecture, including when the binary prefix has arbitrarily many ones. For every sufficiently large $X$, choose a dyadic $P$ between $X/(2(d+1))$ and $X/(d+1)$. The same prime number theorem gives $\gg_{d,e,L}X/\log X$ such weights through $X$. Their offsets $h$ grow proportionally to $N$.

For $d=127$, $P=2^{256}$, $e=2717$, $L=1785$ one such weight has the prime predecessor factor
$$Q=5412450237894181046016723476830419859727020033384827332387109404041963682257.$$
It has a seven-one binary prefix. Both $N$ and $N-1$ have at least four distinct prime divisors. Its offset has $(dh)^5>27M^2$.

At the consecutive weights 5,635 and 5,636 the two-quotient theorem applies with $P=512$, $Q=625$, $d=11$, $e=9$. Their remainder pairs are $(3,10)$ and $(4,11)$. The verifier also checks both weights directly by the independent endpoint search.

## Reproduction

Use python-flint for the integer verifier and Lean 4.36.0 for the formal lemmas. The verifier checks 1,112,188 exact gcd inequalities, including 77,863 at partially retained powers, and 24,102 two-power residue intersections. It certifies primality recursively in both large examples. Lean checks the vanishing of an entire row of discrepancies from the cross-products, and the residue argument in the two-power extension.

```sh
python3 check_two_power_cut_windows.py
lean TwoPowerCutWindows.lean
```
