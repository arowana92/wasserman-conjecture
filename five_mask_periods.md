# Counting possible binary exponent classes by the inverse residue

For an odd prime power $Q=\ell^e$, let $A(Q)$ count the classes modulo $L_Q=\operatorname{ord}_Q(2)$ containing at least one positive $k$ where the trinomial conjecture fails at $1+5\cdot2^kQ$. The following bounds apply also to the explicitly computable set of classes left by the inverse and binary tests below.
$$
 A(Q)\le \sqrt{6Q-2},\qquad
 \sum_{\substack{Q\le X\\ Q\text{ odd prime power}}} A(Q)
       \ll \frac{X}{\sqrt{\log X}} .
$$
The second bound has an absolute constant and sums over all exponents $e\ge1$. The [binary direction theorem](binary_direction_measure.md#whole-progressions-of-weights-at-cofactor-five) gives a qualitative limit for each fixed prime. Here the fraction of classes left is **$O_\ell(Q^{-1/2})$** for each fixed $\ell$. The same rate holds inside any fixed progression of $k$, with its step allowed in the constant.

Write $d=\operatorname{ord}_\ell(2)$ and $v=v_\ell(2^d-1)$. Whenever $e\ge3v$ the fraction is at most $\sqrt{6}\,Q^{-1/6}/d$, uniformly in the prime. This applies in particular to every $e\ge3$ at every odd prime with $2^{\ell-1}\not\equiv1\pmod{\ell^2}$.

## A pointwise binary bound

A possible failure requires that $3Q$ and $2Q$ have disjoint binary ones. It also requires an integer $x$ with $0<x<Q$, $15\mid x$, $2^{k+1}x\equiv1\bmod Q$, and every one of $x$ in a one position of $Q$. This is $x=z/2$ in the [inverse-residue row](factor_five_exponents.md#determining-the-row-from-a-binary-subset). For fixed $Q$, each class of $k$ has a different inverse $x$. Thus $A(Q)\le2^{s_2(Q)}$ whenever the first binary test passes, and is zero otherwise. Here $s_2$ is the binary digit sum.

Read $Q$ from low to high, remembering the preceding bit and the carry of $Q+2Q$. The initial state is 0. A digit is allowed if the output bit and the preceding bit are never both 1. The final two columns can be checked by appending two zeros.

| State | Previous bit, carry | Next on 0 | Next on 1 | Potential $H$ | Weight $w$ |
|---|---|---|---|---:|---:|
| 0 | $(0,0)$ | 0 | 1 | 0 | 2 |
| 1 | $(1,0)$ | — | 2 | 1 | 1 |
| 2 | $(1,1)$ | 3 | — | 2 | 1 |
| 3 | $(0,1)$ | 0 | 2 | 1 | 2 |

On every transition from $i$ to $j$ with digit $b$, $2b-1\le H(j)-H(i)$. Thus the position of the $r$-th one, numbered from zero, is at least $2r-3$ for $r\ge2$. If there are $m$ ones, then
$$
  Q\ge 1+\sum_{r=2}^m 2^{2r-3}=(2\cdot4^{m-1}+1)/3 .
$$
This gives the square-root bound. Equality in the digit inequality is attained by $(2^{2m-1}+1)/3$ for every $m\ge2$, whose bits occupy exactly those positions. No assertion of equality for the number of exceptional classes follows.

For $e\ge v$, lifting gives $L_Q=d\ell^{e-v}$. Indeed raising $1+u\ell^j$, $\ell\nmid u$, to the $\ell$-th power increases the valuation of its difference from 1 by one, while an exponent coprime to $\ell$ preserves the valuation. The fraction of classes is therefore at most $\sqrt6\,\ell^v/(d\sqrt Q)$. If $e\ge3v$ then $\ell^v\le Q^{1/3}$. On a progression of step $s$ the period is $L_Q/\gcd(s,L_Q)$, so the relative bound increases by at most $s$.

## Distribution in odd residue classes

The table satisfies $\sum_{i\xrightarrow{b}j}2^b w(j)=2w(i)$. Give each outgoing edge probability $2^b w(j)/(2w(i))$. Start at state 0 and let $Y_n$ be the integer represented by the first $n$ digits. The probability of a particular path ending at $j$ is $2^{s_2(Y_n)}w(j)/(2^nw(0))$. In particular the sum of $2^{s_2(Q)}$ over passing $Q<2^n$ in any set is at most $2^{n+1}$ times the probability of that set.

There are absolute constants $C,c>0$ such that for every odd $D\ge3$ and every residue $b$,
$$
  \left|\Pr(Y_n\equiv b\bmod D)-\frac1D\right|
       \le C\exp\left(-\frac{cn}{\log(2D)}\right).                 \tag{1}
$$
To prove this, twist the transition probabilities at position $h$ by the phase $\exp(2\pi i a2^h b'/D)$ for digit $b'$, where $D\nmid a$. The maximum row sum of absolute values is at most 1. From every initial state a prescribed three-digit prefix can reach 0. After it, the six-digit blocks representing $0,3,11$, padded at their high ends, all return to 0. Choose $u\in\{3,11\}$ with $D\nmid au$, which exists since these integers are coprime.

In any block starting at position $h$, there is $0\le j\le\lceil\log_2 D\rceil$ such that $2^{h+j+3}au/D$ has distance at least $1/4$ from an integer. Its initial distance is at least $1/D$ since $D$ is odd, and doubling doubles the distance while it is below $1/4$. Use the three-digit prefix at position $h+j$, followed by either the word for 0 or the word for $u$. From each starting state the two length-nine paths reach the same state, each with probability at least $2^{-9}$ by the telescoping formula. Their phases are separated by at least a right angle. Removing equal mass $2^{-9}$ from each contribution, bounding its sum by $\sqrt2\,2^{-9}$, and applying the triangle inequality to the rest gives maximum row sum at most $1-(2-\sqrt2)2^{-9}$ for these nine positions. All other factors are contractions. Every consecutive block of length $\lceil\log_2 D\rceil+9$ therefore gives this bound on the matrix product. Iterating bounds every nonconstant Fourier coefficient, and averaging characters proves (1).

## Sieving the prime powers

Use $n=\lceil\log_2(X+1)\rceil$ and $z=\exp(\gamma\sqrt n)$ for a sufficiently small absolute $\gamma>0$. We give the upper-bound sieve calculation to specify the required distribution estimate. Take the odd squarefree integers $r\le z$, including 1, and put
$$
   G=\sum_r 1/\varphi(r),\qquad
   y_r=\frac{\mu(r)}{G\varphi(r)},\qquad
   \lambda_d=d\sum_{d\mid r}\mu(r/d)y_r
$$
on the same set of $d$. Divisibility inversion gives $y_r=\sum_{r\mid d}\lambda_d/d$, $\lambda_1=1$ and $|\lambda_d|\le d$. The divisor sums range only over this finite set. For an integer having no odd prime divisor at most $z$, the square of $\sum_{d\mid Y_n}\lambda_d$ is one. Its expectation bounds the probability of this event from above. The main term obtained using (1) is
$$
   \sum_{d,e}\frac{\lambda_d\lambda_e}{[d,e]}
     =\sum_r\varphi(r)y_r^2=1/G ,
$$
by the identity $\sum_{r\mid\gcd(d,e)}\varphi(r)=\gcd(d,e)$. The total error is $O(z^4\exp(-cn/\log(2z^2)))$. For small enough $\gamma$ this is $o(n^{-1/2})$. Also $G\gg\log z$. Indeed the odd numbers through $y$ number at least $y/2-1$, and at most $y\sum_{j\ge3}j^{-2}\le4y/9$ of them have a square factor greater than 1. Partial summation of the resulting lower bound for odd squarefree counts gives the claim, since $1/\varphi(r)\ge1/r$. The probability is consequently $O(n^{-1/2})$.

All odd prime powers with base greater than $z$ pass this sieve. For the smaller prime bases there are $O(zn)$ powers through $X$, and each contributes at most $\sqrt{6X}$. This gives $O(zn\sqrt X)+O(2^n/\sqrt n)=O(X/\sqrt{\log X})$ as asserted.

Run **python3 check_five_mask_orbits.py** for direct multiplication checks of the entire binary graph through 17 bits, the extremal strings, two enumerations of the inverse classes at odd moduli below 40,000, and checks of the short paths and the exact sieve identities. [Integer verifier](check_five_mask_orbits.py)
