# Intervals around smooth numbers and recurrence terms

Suppose $P,R$ are prime powers, $10P+7=DR$, and $1\le D\le10000$. Put $N=(2P+1)(5P+1)$. The following criteria establish Wasserman's full trinomial conjecture at $N$. **Each criterion can be used with any one of the four integers**
$$Q=P,\quad R,\quad 2P+1,\quad 5P+1.$$

* If $|Q-Z|\le Q^{0.9999}$ for a 17-smooth positive integer $Z$, the conjecture holds. At 19-smooth integers the radius is $Q^{0.99}$, and at 23-smooth integers it is $Q^{1/5}$. An integer is $y$-smooth if every prime divisor is at most $y$. The same radii work for a positive rational $Z$ whose reduced numerator and denominator are both $y$-smooth and at most $Q^2$.
* The radius $Q^{0.99}$ also works for products $Z=\prod_{i=1}^k g_i^{e_i}$, with arbitrary positive exponents, integer generators $g_i\ge2$, and the following bounds. The generators may vary with the input.

  | Number of generators | Largest allowed generator |
  |---:|---:|
  | 1 | $2^{10^{13}}$ |
  | 2 | $2^{500000}$ |
  | 3 | $2^{2000}$ |
  | 4 | $2^{100}$ |

* The radius $Q^{0.99}$ works around **every Fibonacci, Lucas and Pell number**. More generally fix the recurrence $Z_{n+2}=aZ_{n+1}-eZ_n$ with $e=-1,\ 1\le a\le 2^{5000000}$ or $e=1,\ 3\le a\le 2^{5000000}$. Take either initial pair $(0,1)$ or $(2,a)$. Every term of every one of these sequences, at every index $n\ge0$, gives such an interval.

Thus, for example, a possible failure requires a prime factor greater than 23 in **each** of $2P+1$ and $5P+1$, and in every positive integer within the stated fifth-root radius of either. The two factors themselves can have arbitrarily many prime divisors. The bounds hold from the first input by using the [independently checked finite range](seven_ternary_rotation.md).

## Binary approximation and an explicit gap

In a failure the [binary ray theorem](seven_quadratic_prime_bases.md) gives $|P^2-c2^m|\le1100000P$ with $c=8Df/5$ for one of six rational $f$. For the selected $Q=sP+b$ use $\gamma=s^2c$, where $(s,b)$ is $(1,0),(10/D,7/D),(2,1)$ or $(5,1)$. Expansion of the square gives
$$
 |Q^2-\gamma2^m|<10^7 Q,\qquad 10^{-3}<\gamma<10^6,\qquad h(\gamma)<26 .
$$
Here $h$ of a reduced positive rational is the logarithm of the larger of its numerator and denominator. All twenty-four $\gamma$ have odd 5-adic valuation. Thus a square of any positive rational is different from $\gamma2^m$. The checks of height and valuation use the actual ray fractions.

By the finite theorem we may assume $P>2^{10^{30}}$. Since $Q>P/1000$, write $L=\log Q>L_0=6.8\cdot10^{29}$; all logarithms here are natural. The scale satisfies $0<m<4L$. If $|Z-Q|\le Q^{1-\delta}$ for any $1/10000\le\delta\le4/5$, the square approximation gives
$$
 0<\left|Z^2/(\gamma2^m)-1\right|\le 8Q^{-\delta}
 \tag{1}
$$
for rational $Z>0$. Indeed $|Z^2-Q^2|\le3Q^{2-\delta}$ and $10^7Q\le Q^{2-\delta}$ in this range, while $\gamma2^m>Q^2/2$.

For $t$ logarithms in a real number field of degree $d$, put
$$
 C=\tfrac75\,30^{t+3}t^{9/2}d^2(1+\log d)\prod_i A_i ,
 \qquad A_i\ge\max(dh(\alpha_i),|\log\alpha_i|,0.16).
$$
Matveev's multiplicative bound for a nonzero $\Lambda=\prod\alpha_i^{b_i}-1$ gives $\log|\Lambda|>-C(1+\log B)$ whenever $B\ge\max_i|b_i|$ and the exponents are integers. All numbers here are positive reals. Zero exponents or bases equal to one can be omitted; the displayed upper estimates for $C$ still dominate the shortened products since $30A_i>1$. [Matveev's paper](https://www.mathnet.ru/eng/im314) · [Explicit multiplicative formulation, Theorem 2.1](https://doi.org/10.1016/j.jnt.2020.05.005)

For smooth rationals in (1) take $d=1$ with the bases $\gamma$ and all primes at most $y$, combining $2^{-m}$ with the 2-exponent of $Z^2$. Use heights $26$ and $\log p$. Since the numerator and denominator are at most $Q^2$, all exponents are at most $10L$ in magnitude. For an integer product of $k$ varying generators, use $\gamma,2,g_1,\ldots,g_k$ and heights $26,\log2,\log g_1,\ldots,\log g_k$. Here $Z\le2Q$ gives the same coefficient bound. In all rows the checker verifies
$$
                 10+75 C<\delta L_0.
$$
Since $1+\log(10L_0)<75$ and $(10+C(1+\log(10L)))/L$ is decreasing for $L\ge L_0$, the upper estimate (1) contradicts the lower estimate. The certificate uses rational upper bounds for each logarithm and for $\sqrt{t}$.

## Recurrence terms

Let $\Delta=a^2-4e$, $\lambda=(a+\sqrt\Delta)/2$, $\mu=e/\lambda$. The two allowed ranges give a nonsquare $\Delta$, $\lambda>3/2$ and $|\mu|<1$. The two sequences are respectively $(\lambda^n-\mu^n)/\sqrt\Delta$ and $\lambda^n+\mu^n$. Put $g=\sqrt\Delta$ in the first case and $g=1$ in the second. In both cases $|Z_n-\lambda^n/g|\le1$.

For $|Z_n-Q|\le Q^{0.99}$ and the large range above, $n>0$ and
$$
 0<\left|\frac{\lambda^{2n}}{g^2\gamma2^m}-1\right|\le16Q^{-0.01}.
$$
The bound follows from $|\lambda^n/g-Q|\le Q^{0.99}+1$ as in (1). Equality with one before subtracting is impossible: conjugation would give $\lambda^{2n}=\mu^{2n}$ since $g^2\gamma2^m$ is rational. Also $\lambda^n<2gQ$ gives $2n<10L$.

Use the same logarithm bound in the real quadratic field, with bases $2,g^2\gamma,\lambda$. For $v=5000000$ and $H=(v+2)\log2$, permissible height bounds are $2\log2,\ 2(26+2H),\ H$. Indeed $\sqrt\Delta,\lambda<a+2\le\exp(H)$, and the quadratic unit has $2h(\lambda)=\log\lambda$. The integer check again gives $10+75C<L_0/100$. This proves the stated recurrence criterion.

## A growing lower bound for prime factors

For every fixed $0<\delta<1$, if $Q$ comes from a counterexample and $Z$ is a positive integer with $|Z-Q|\le Q^{1-\delta}$, then, as $Q$ tends to infinity through such inputs, its largest prime factor satisfies
$$
 P^+(Z)\ \ge\ (1-o(1))\,
 \frac{(\log\log Q)(\log\log\log Q)}{\log\log\log\log Q}.
$$
This holds uniformly in $Z$ for fixed $\delta$ and the four choices of $Q$. At sufficiently large $Q$ the derivation of (1) works for every fixed $\delta<1$. Put $y=\max(2,P^+(Z))$ and use $\gamma$ together with all primes up to $y$ as before. Bounded $y$ contradicts $\delta L\le O(1)+C(1+\log(10L))$ for large $L$. For $y\to\infty$, the prime number theorem gives
$\log C\le(1+o(1))(y/\log y)\log\log y$, by bounding the product of the prime heights by $(\log y)^{\pi(y)}$. Taking logarithms in the same inequality yields $\log L\le(1+o(1))y\log\log y/\log y$; the discarded term is $O(\log\log L)=o(\log L)$. The function on the right is eventually increasing, and substitution of $x\log x/\log\log x$ into it gives $(1+o(1))x$. Invert with $x=\log L$ to obtain the formula.

~~~~sh
python3 check_quadratic_smooth_gaps.py
~~~~
