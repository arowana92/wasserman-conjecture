# Every counterexample must contain the endpoint row

**Theorem.** At every integer weight $N\ge3$, any three positive trinomial rows of weight $N$ whose coefficients have gcd one must include a permutation of $E=(1,1,N-2)$. Consequently Wasserman's full three-row conjecture is equivalent to the following two-row assertion:
$$\gcd\big(N(N-1),K_A,K_B\big)>1\quad\hbox{for every pair of positive rows }A,B\hbox{ of weight }N.$$

The earlier [reduction through $10^{20}$](critical_small_parts.md) covers the initial range. We prove the reduction above that weight using products of binomial coefficients with cutoffs depending on the parts. The two-row assertion remains open in general.

## A normalized divisibility lemma

Put $F_t(x)=\prod_{j=1}^t\binom{x}{j}$ for integers $0\le t\le x$, so $F_0(x)=1$. Let $1\le H\le N$ and $M_H=\operatorname{lcm}(1,\ldots,H)$. For each row $X$ independently choose cutoffs $t_{X,i}\le X_i$ with sum $H$. If the coefficients of the three rows have gcd one, then
$$F_H(N)\ \mid\ M_H^H\prod_{X,i} F_{t_{X,i}}(X_i). \tag{1}$$
In fact it suffices that every prime dividing $F_H(N)$ is avoided by at least one row.

For a prime-power layer $w=p^e$, Legendre's formula gives its contribution to the valuation of $\binom{x}{j}$ as
$$\lfloor x/w\rfloor-\lfloor j/w\rfloor-\lfloor(x-j)/w\rfloor\ \in\{0,1\}.$$
For $w\le H$, the total contribution to $F_H(N)$ is at most $H$, supplied by $M_H^H$ at that layer. All factors on the right contribute nonnegatively. For $w>H$, write $h=N\bmod w$. The left contribution is $(H-h)_+$. Choose a row avoiding $p$, whose three residues $r_i$ modulo $w$ sum to $h$ by Kummer's criterion. Its contribution is $\sum_i(t_{X,i}-r_i)_+\ge (H-h)_+$. Summing over layers proves (1).

## Choosing a height from the parts

Suppose a counterexample omits $E$ at $N>10^{20}$. The [explicit small-part estimates](critical_small_parts.md) give $s>2000$ for the sum of the two smallest parts of every row, and exclude all individual parts from 4 through 1000. At most one of the nine parts is at most 3:

- Two in the same row would violate $s>2000$.
- If two different rows have parts at most 2, the powers of $N$ assigned to those rows together divide 2. This leaves the third row with three positive normalized parts of sum at most 2.
- If both selected parts are at least 2, the [product bound for parts in different rows](parts_in_different_rows.md) requires their product to exceed $27^{1/5}(N-1)^{2/5}$.
- A 1 and a 3 in different rows are excluded by [the three-consecutive-integers argument](near_row.md#which-rows-can-have-small-parts).

Let $m$ be the least of the parts exceeding 3, and set $H=6\lfloor m/6\rfloor$. Thus $H\ge996$ and $H\le m\le H+5$. If $r$ is the largest prime at most $N$, a row avoiding $r$ has two smaller parts summing to at most $N-r$, since $r>N/2$. At least one exceeds 3, so $m\le N-r$. Kadiri and Lumley's Theorem 1.1, Table 2 at $\log x_0=46$, gives $N-r<N/148306019$ for $N>10^{20}>e^{46}$. Put $L=\log(N/H)$, so $H/N<10^{-8}$ and $L>18$. [Kadiri and Lumley](https://www.cs.uleth.ca/~lumley/shortintervals.pdf)

Give the part $m$ cutoff $H$ and give its two companions zero. Of the remaining two rows one has all three parts at least $H$, and the other has at least two such parts. Use three equal cutoffs in the former and two equal cutoffs in the latter. Write $W(t)=t(t+1)/2$. All logarithms below are natural.

## Explicit estimates in (1)

**The selected part.** The binary entropy function $b(u)=-u\log u-(1-u)\log(1-u)$, extended continuously to the endpoints, satisfies $\log\binom mj\le m b(j/m)$ by the binomial theorem. Concavity gives $\sum_{j=1}^m b(j/m)\le m\int_0^1 b(u)\,du=m/2$ by the trapezoid rule. Thus
$$\log F_H(m)\le m^2/2<0.51H^2.$$

**The two ordinary rows.** For a row with $q=2$ or 3 selected parts $x_i$ put $t=H/q$. Their product is at most $(N/q)^q$. Using $\binom{x_i}{j}\le (e x_i/j)^j$ gives
$$\log\prod_i F_t(x_i)\le qW(t)(L+1)+q t^2/4.$$
Here $\sum_{j=1}^t j\log(t/j)\le\int_0^t x\log(t/x)\,dx=t^2/4$ by concavity and zero endpoint values. Thus the sum for both rows, divided by $H^2$, is at most
$$(5/12+1/H)L+5/8+1/H\le(5/12+1/H)L+0.63.$$

**The small layers.** Theorem 12 of Rosser and Schoenfeld gives $\log M_H=\psi(H)<1.04H$. Their bound applies at every positive $H$. [Rosser and Schoenfeld](https://doi.org/10.1215/ijm/1255631807)

**The left side.** For $j\le H$, use $\binom Nj\ge (N-H)^j/j!$ and $\log(j!)\le j\log j-j+1+\log j$, by integration of $\log x$. The sum of $j\log(H/j)$, $1\le j\le H$, is at least $H^2/4-H$. Indeed its error from the integral is at most the total variation $2H/e<H$. Hence
$$\log F_H(N)\ge W(H)(L+1+\log(1-H/N))+H^2/4-H(2+\log H).$$
Since $H\ge996$, $(2+\log H)/H<9/996$ by monotonicity, and $\log(1-H/N)>-2\cdot10^{-8}$. The right side exceeds $H^2(L/2+0.73)$.

Adding the upper bounds in (1) and comparing with the lower bound gives
$$(1/12-1/H)L < 1.04+0.51+0.63-0.73=1.45.$$
But $(1/12-1/996)18>1.48$. This contradiction completes the proof.

## Consequences

The [endpoint pair-sum bound](endpoint_pair_bound.md) now applies to the two other rows at **every** possible exception. In particular, if $Q=p^j\le N<2Q$, $p\mid N(N-1)$, and
$$(N-Q)^7\le6^6(N-1)^4,$$
the full conjecture holds at $N$. A row avoiding $p$ must differ from $E$ and must have two parts summing to at most $N-Q$. This contradicts the endpoint bound. With $p=2$ this proves the full conjecture throughout intervals of order $2^{4k/7}$ above every power $2^k$, with the stated explicit constant. Similarly the [two-prime-factor predecessor theorem](endpoint_reduction.md) now proves the full conjecture at **every** weight $N$ for which $N-1$ has at most two distinct prime divisors.

Run `python3 check_normalized_layers.py` to reproduce the prime-power checks and rational margins. The argument uses published explicit bounds and the previously documented finite reduction.
