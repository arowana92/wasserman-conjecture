# Six powers of seven near a power of two

**The largest integer $E$ for which a sum $n$ of one to six positive powers of seven satisfies $|n-2^E|\le4\sqrt{2^E}$ is 17.** Repetitions are allowed. The last example is $n=5\cdot7^4+7^6=129654$, at distance 1418 from $2^{17}$.

**Consequently, for every product $t=p^a q^b$ of two prime powers, at least one of $7t+1$ and $16t+1$ satisfies Wasserman's full three-row conjecture.** In the [five-weight theorem](five_related_weights.md), the only pair of multipliers that can still fail together is 9 and 16. This removes the earlier finite but unexhausted range for 7 and 16.

## A bound at one prefix

Assume $E>200$, put $X=2^E$, and write $n=\sum_{i=1}^s7^{b_i}$ with $s\le6$, $b_1\ge\cdots\ge b_s\ge1$. Put $d_i=b_1-b_i$. The prefix of length $j$ is $F=7^{b_j}R$, with $R$ a positive integer at most $6\cdot7^{d_j}$. Both $X/F-1$ and $\log_7(X/F)$ have absolute value at most
$$100\big(2^{-E/2}+7^{-d_{j+1}}\big), \tag{1}$$
omitting the second term at a full prefix. Indeed $X/2<n<2X$, and the sharper $F\ge n/6\ge (1-4/\sqrt X)X/6>X/7$ holds. The unused tail is at most $6n7^{-d_{j+1}}$. Use the mean value inequality for the logarithm and $\log7>1$. Also $F\ne X$ since $7\mid F$.

Matveev's lower bound for a nonzero product of rational powers minus one gives
$$ \log|X/F-1|>-10^{13}(1+\log(E+1))(2d_j+3). $$
Use the positive rationals $2,7,R$ with exponents $E,-b_j,-1$, dropping $R$ when it equals 1. Here $b_j\le E+1$, $\log R\le2d_j+2$, and $1.4\cdot30^6\cdot3^{4.5}\log2\log7<10^{13}$. [Matveev's bound in rational form, Theorem 2.1](https://archive.mpim-bonn.mpg.de/id/eprint/5025/1/mpim-preprint_2023-34.pdf)

With $C=10^{14}(1+\log(E+1))$, comparison with (1) gives either $E<C(d_j+2)$ or $d_{j+1}+2<C(d_j+2)$. At a full prefix only the first alternative is available. In fact the smaller of $(E/2)\log2$ and $d_{j+1}\log7$ is less than $\log200+10^{13}(1+\log(E+1))(2d_j+3)$. Since $d_1=0$, iteration yields $E<2C^6$. This implies $E<10^{180}$ by evaluation at that point and monotonic decrease of $(1+\log(E+1))^6/E$ thereafter.

## Separating prefixes by continued fractions

Suppose $E\le M$. Set $\alpha=\log_7 2$. Find coprime $h,Q$ with $Q>1000(M+3000)$ and $64(M+3000)|Q\alpha-h|<1$.

* If $R=2^v7^w$ with $v<E$, then $0<E-v<Q$ and $(E-v)h-(b_j+w)Q$ is a nonzero integer. Consequently $|\log_7(X/F)|>1/(64Q)$.
* Otherwise we check $\|Q\log_7 R\|>1/32$, where the norm is the distance to the nearest integer. Subtracting $E(Q\alpha-h)$ costs less than $1/64$. This gives the same strict lower bound.

The programs select a denominator separately at every prefix and check its inequalities. Put $g=\lceil\log_7(12800Q)\rceil$. Whenever $2^{E/2}>12800Q$, (1) shows that this prefix cannot be the full sum and that its next gap must satisfy $d_{j+1}<g$. Enumerate all integers $d_j\le d_{j+1}<g$, updating the coefficient to $7^{d_{j+1}-d_j}R+1$. The non-strict ordering includes repeated powers. The initial coefficient is 1.

| Current bound $M$ | Last tested length $J$ | Prefixes checked | Maximum $g$ at length $J$ | New upper bound on $E$ |
|---:|---:|---:|---:|---:|
| $10^{180}$ | 2 | 223 | 225 | $10^{74}$ |
| $10^{74}$ | 3 | 4,858 | 101 | $10^{55}$ |
| $10^{55}$ | 4 | 73,390 | 79 | $10^{38}$ |
| $10^{38}$ | 5 | 457,416 | 63 | $10^{21}$ |
| $10^{21}$ | 6 | 581,312 | 38 | 3000 |
| 3000 | 6 | 12,140 | 18 | 200 |

In the first five passes assume $E>3000$, and in the last assume $E>200$. Each test verifies $v$ is at most the assumed lower threshold in the homogeneous case and $12800Q$ is below two to the power of half that threshold. At $J<6$, the next gap is strictly below the tabulated bound $G$. Applying the preceding Matveev recurrence at the remaining at most $6-J$ prefixes gives $E<(G+2)C^{6-J}$. For $M=10^b$ the verifier substitutes $C<10^{14}(2+3b)$ and checks the new upper bound by integers. At $J=6$ no next term is available, which proves the asserted reduction directly.

Through $E=200$, the first checker enumerates the powers in the exact integer window of radius $\lfloor\sqrt{16X}\rfloor$. The second independently finds the first positive multiple of seven above or equal to its lower endpoint whose base-seven digit sum is at most six. After division by seven, either the lower bound already qualifies or the smallest qualifying number agrees up to a first raised digit and then has zeros below it. Rounding up at each successive power of seven therefore finds this minimum.

## The two weights

If both weights failed, the [binary row calculation and the shared-part bound](five_related_weights.md#multipliers-7-and-16) would give $t\equiv587\pmod{1008}$ and integers $m,k\ge0$ with $2k>3m+2$ and
$$42t=32\cdot8^m(14\cdot4^k+1)-18.$$
The digits of $t$ decrease weakly from low to high in base 7, starting at 6. Thus $6(t+1)$ is a sum of six positive powers of seven. Put $E=6+3m+2k$. Then
$$6(t+1)=2^E+(2^{5+3m}+24)/7.$$
Here $E\ge9+6m$, so the positive error is below $4\sqrt{2^E}$. The separation theorem forces $E\le17$. The resulting possible $t$ are $171,683,2731,10923,5467,21851$. None is 587 modulo 1008. This proves the assertion.

## Reproduction

Run `python3 check_six_sevens.py` with python-flint and `python3 check_six_sevens_cross.py` with the standard library. The first uses Arb balls. The second verifies all 1,129,339 prefixes with integer logarithm enclosures and checks the six terminal candidates. For $1\le n/d\le2$, write $t=(n-d)/(n+d)$. At fixed integer scale $S$, approximate $St^{2j+1}$ from below by multiplying the previous approximation by $\lfloor St^2\rfloor/S$ and rounding down. The error remains at most three since $t^2\le1/9$ and each update adds at most two. The series $2\sum t^{2j+1}/(2j+1)$ then costs at most eight units of scaled error per step and eight for its tail when the lower term first reaches zero. Range reduction by powers of two and a table of rationals with denominator 256 keeps the series short.
