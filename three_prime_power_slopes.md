# Three prime-power factors in the predecessor

Let $p,q,r$ be distinct primes and let $P=p^a$, $Q=q^b$, $R=r^c$, with positive exponents. We study the full three-row conjecture at $N=1+PQR$.

**Unrestricted third exponent.** For any three of the primes $2,3,5,7$, the conjecture holds whenever at least two of their full powers are at most $10^7$. The third exponent is arbitrary. For the prime sets $\{2,3,11\}$ and $\{2,3,13\}$ the same assertion holds with bound $10^6$. These conclusions cover 2,902 choices of the varying base and the two fixed powers. The families can overlap.

**Three exponents varying.** Fix any prime $r\ge5$. In every interval of $H$ consecutive positive exponents $c$, all but $O_r(H^{0.8492})$ of them satisfy the conjecture **simultaneously** at all $1+2^a3^b r^c$ with
$$
             2^a3^b\le H^{0.7755}
       \quad\text{or}\quad \min(2^a,3^b)\le H^{0.1122}.
$$
In the second range the other power has no size bound. A separate full-period estimate controls all $2^a,q^b$ above a given threshold for any fixed odd prime $q$.

## Necessary row and two local tests

In a failure, the [all-weight endpoint reduction](all_weights_endpoint_reduction.md) and [three-factor assignment](endpoint_factors.md#when-n-1-has-three-distinct-prime-factors) force a row avoiding $p,q,r$. The parts exceptional at $P$ and $Q$ determine it as
$$
              (xQR,\ yPR,\ 1+(PQ-xQ-yP)R),
 \qquad x=(QR)^{-1}\bmod P,\quad y=(PR)^{-1}\bmod Q. \tag{1}
$$
Use least positive inverses, so $x<P$ and $y<Q$. The slope $PQ-xQ-yP$ is positive. A negative value is impossible and zero would give $P\mid xQ$. The assignment also gives a unitary divisor $U\ge3$ of $N$ dividing $x,y$ whose prime divisors the row avoids. Indeed $U$ divides each part and is coprime to $PQR$.

At base $r$, removal of the low $c$ places of (1) leaves the three slopes adding to $PQ$ without carries. At base $p$, grouping the first and third parts and removing the low $a$ places gives the necessary carry-free sum $yR+(Q-y)R=QR$. For the first $a$ digits $R\equiv(xQ)^{-1}\bmod P$. The corresponding sum at $q$ is $xR+(P-x)R=PR$, with $R\equiv(yP)^{-1}\bmod Q$. Thus for fixed $P,Q$ we can test all three bases without knowing $c$. These arguments use the multinomial carry criterion. [Bergman, Lemma 5](https://arxiv.org/abs/0806.0607)

## Exhaustion at the six sets of bases

Enumerate the pairs $x,y$ for which the three slopes add freely at $r$, the inverse numerators are units at their respective primes, and the factor of $\gcd(x,y)$ coprime to $r$ is at least three. The next two columns successively require the low digits of the grouped sums at $p$ and $q$, with $p<q$. The last column uses the original row and an assigned prime of $N$.

| Prime set | Fixed-power choices and varying base | Slope pairs | After $p$ | After $q$ | Remaining |
| --- | ---: | ---: | ---: | ---: | ---: |
| 2, 3, 5 | 692 | 184 | 8 | 0 | 0 |
| 2, 3, 7 | 618 | 1,524 | 54 | 0 | 0 |
| 2, 5, 7 | 494 | 3,647 | 198 | 12 | 0 |
| 3, 5, 7 | 332 | 59,830 | 511 | 16 | 0 |
| 2, 3, 11 | 383 | 3,945 | 113 | 6 | 0 |
| 2, 3, 13 | 383 | 11,824 | 303 | 10 | 0 |

For the last column use exponent period $L=12$ in the first four cases and $L=60$ in the last two. Some odd prime $\ell$ of $\gcd(x,y)$ coprime to $r$ must be assigned to the row, or its assigned divisor is a power of two divisible by four. In the latter case take $\ell=2$ and require $4\mid N,x,y$. For each candidate and each $e\bmod L$, $R\equiv r^e\bmod m$ at every power $m$ of $p,q,\ell$ dividing $r^L-1$. Use these layers to test the known residues $xQR\equiv1\bmod P$, $yPR\equiv1\bmod Q$, $N\equiv0\bmod\ell$ (modulo four for the binary alternative), and carry freedom of the three original parts. Every residue of $c$ is included.

The two searches use different enumerations. The C++ program splits the base-$r$ digits of $PQ$ into two sets and matches residues to generate its digit subwords divisible by $P$ or $Q$. The Python program generates the two multipliers from low digits modulo powers of $r$, groups them by exponent residue first, and checks additions by factorial floors. Both run entirely with integer arithmetic. [Digit-subset search](check_three_prime_power_slopes.cpp) · [Independent residue and floor search](check_three_prime_power_slopes_cross.py)

## Counting both ungrouped additions

Fix an exponent progression $c=c_0+s i$, with positive exponents, and write $v=v_p(r^s-1)$, $w=v_q(r^s-1)$. Assume $v,w\ge1$ and require valuation at least two at the binary prime. Choose $v\le j\le a$, $w\le k\le b$. In any interval of
$$
                       T=p^{a+j-v}q^{b+k-w}
$$
consecutive indices, $R$ ranges once over the product of its fixed unit classes modulo $p^v,q^w$, lifted to $Pp^j,Qq^k$. Indeed taking the $p$th power of $1+p^v u$, with $p\nmid u$, increases its valuation after subtracting one by exactly one. Binomial expansion and iteration give order $p^{a+j-v}$ at the first modulus. The second order is $q^{b+k-w}$ and the orders are coprime.

For $h>1$ define
$$
 M_p(h)=p^{-2}\sum_{d=0}^{p-1}\sum_{m=0}^d
             \left(\frac{d-m+1}{p}\right)^h,
 \qquad C(n,m)=\frac{m\lceil n/m\rceil}{n}.
$$
For conjugate exponents $h,h'>1$, the fraction of possible failures in any such block is at most
$$
 \left[C(Q/q^w,p^j)M_p(h)^{j-v}\right]^{1/h}
 \left[C(P/p^v,q^k)M_q(h')^{k-w}\right]^{1/h'}. \tag{2}
$$

To prove this, first fix $R\bmod P,Q$, or equivalently $x,y$. Before conditioning, $x$ and $y$ are independent and each is uniform on a single unit class modulo $p^v$ and $q^w$ in its range. Immediately above the $P$ shift the total has digits of $QR$, the second part has digits of $yR$, and the first part has digits of $z=(xQR-1)/P$. Varying the next $j$ digits of $R$ makes $z$ uniform modulo $p^j$, since its multiplying coefficient $xQ$ is a unit. Write $d_i,m_i$ for the first $j$ digits of $QR,yR$. The conditional fraction passing the ungrouped addition at $p$ is exactly
$$
                 F_p(x,y)=\prod_{i=0}^{j-1}\frac{(d_i-m_i+1)_+}{p}.
$$
The bottom $a$ digits already fit. At $q$ define $F_q$ by interchanging the two powers and multipliers. The upper lifts at the two primes are independent, so the exact fraction passing both digit tests is $\mathbb E(F_pF_q)$. Positivity of the slopes and the weight divisor are additional requirements for a failure.

For fixed $x$, the residue $R\bmod p^j$ is fixed and invertible. The $Q/q^w$ possible values of $y$ are consecutive terms of a progression of step $q^w$. Thus the probability of each value of $yR\bmod p^j$ is at most $C(Q/q^w,p^j)/p^j$. Replace that variable by a uniform digit word at this cost. The total $QR\bmod p^j$ equals $x^{-1}\bmod p^j$ and ranges uniformly over the lifts of its fixed class modulo $p^v$. Bound the bottom $v$ digit factors by one. Averaging the other $j-v$ independent digit pairs yields $\mathbb E F_p^h\le C(Q/q^w,p^j)M_p(h)^{j-v}$. Applying Hölder to $F_p,F_q$ gives (2), including their dependence through $x,y$.

In particular $M_p(2)=(p+1)^2(p+2)/(12p^3)$ by summing squares. Take $p=2,q=3,r=5,s=2$, giving $v=3,w=1$, and inspect $j=a,k=b$. At $(a,b)=(64,40)$ the possible failure fraction in every full period of either parity is less than $7.13\cdot10^{-25}$. At $(128,81)$ it is less than $5.88\cdot10^{-51}$. At the smaller $(5,3)$ the exact counts passing just the two ungrouped digit tests are 12 of 31,104 odd and 120 of 31,104 even exponent classes. The verifier compares exact lift counts with direct floor tests, including those with unequal inspected lengths. [Exact moment and period checks](check_two_prime_digit_moments.py)

## Simultaneous estimates over an interval of exponents

Fix $p,q,r$, and split an interval of $H$ consecutive positive exponents into classes with step $s=\operatorname{lcm}(2,p-1,q-1)$. Constants below may depend on the three primes and the chosen Hölder exponents. Put
$$
  \alpha_p=-\log_p M_p(h),\quad \alpha_q=-\log_q M_q(h'),\qquad
  \kappa=\alpha_p/h+\alpha_q/h'.
$$
Whenever $T\le H$, covering the interval in each progression by full blocks uses at most $O_{p,q,r}(H)$ indices in total.

First suppose $\alpha_p,\alpha_q>1$. In any interval of $H$ positive exponents, the fraction permitting a failure at **any** $P,Q\ge Y$ with $PQ\le\sqrt H$ is $O(Y^{-\kappa})$. Take $j=a,k=b$. If $P\ge Q$, (2) is
$$
              \ll Q^{-\kappa}(Q/P)^{(\alpha_p-1)/h},
$$
since $C(n,m)\le1+m/n$. Summing first over the geometrically spaced $P\ge Q$ and then over $Q\ge Y$ gives the bound. Interchange the primes when $Q\ge P$. Bounded $Y$ can be absorbed in the constant.

For any distinct $p,q\le7$ one may take $h=h'=2$, since $M_p(2)<1/p$ at each of these four primes. In particular $\kappa=(\log_2(8/3)+\log_3(81/20))/2=1.344102\ldots$ at $2,3$. For $p=2$ and **any odd prime $q$** one can obtain
$$
                \kappa>1+\log_2(4/3)(q^{-1/2}-q^{-1}).
$$
Take $h'=q/(\sqrt q-1)$ and $h=h'/(h'-1)$. Writing $i=q-(d-m+1)$ in the moment gives
$$
 M_q(h')\le q^{-2}\sum_{i\ge0}(i+1)e^{-ih'/q}
           =\frac{1}{q^2(1-e^{-h'/q})^2}
           <(1/q+1/h')^2=1/q.
$$
The strict inequality follows from $e^u>1+u$. Here $M_2(h)=(1+2^{1-h})/4$. On $0\le t\le1$ the concave function $-\log_2((1+2^{-t})/2)$ is at least $t\log_2(4/3)$ by its endpoint values. Apply this with $t=h-1$ to bound $(\alpha_2-1)/h$ from below by $\log_2(4/3)/h'$. This establishes convergence in both orders of summation.

Finally put $\beta_u=\log_u(1/(2\sin(\pi/(4u+2))))$ and $\beta=\max(\beta_p,\beta_q)$. For any $1\le Y\le\sqrt H$, except at
$$
            O\!\left(YH^\beta+H Y^{-\kappa}\log^2(2H)\right) \tag{3}
$$
exponents, the conjecture holds simultaneously when $\min(P,Q)\le Y$ **or** $PQ\le H/Y^2$. For the first condition, fix $P\le Y$. The earlier [one-free-exponent estimate](three_factor_exponent_rectangles.md) bounds possible failures for any $Q$ by $O_q(PH^{\beta_q})$. Summing over powers $P\le Y$ costs $O_p(Y)$, and the same reasoning covers $Q\le Y$. For both powers larger than $Y$, take $j=\lfloor\log_pY\rfloor$, $k=\lfloor\log_qY\rfloor$. We can assume $j\ge v,k\ge w$, since bounded $Y$ is trivial. Now $T\le PQY^2\le H$ and the two concentration factors are bounded by $1+q^w$, $1+p^v$. Formula (2) costs $O(Y^{-\kappa})$ per pair and there are $O(\log^2(2H))$ pairs.

For $2,3$ use the second moments and $Y=H^u$, $u=(1-\beta)/(1+\kappa)=0.11223155\ldots$. The two powers of $H$ in (3) coincide at $\beta+u=0.84914932\ldots$ and $1-2u=0.77553689\ldots$. Absorbing the logarithms proves the stated bound.

~~~sh
mkdir -p .work
c++ -std=c++17 -O3 check_three_prime_power_slopes.cpp -o .work/check_three_prime_power_slopes
.work/check_three_prime_power_slopes
python3 check_three_prime_power_slopes_cross.py
python3 check_two_prime_digit_moments.py
~~~
