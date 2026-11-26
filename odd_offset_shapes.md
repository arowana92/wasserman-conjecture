# Odd-prime offsets uniform in both upper exponents

Let $p\ge5$ be prime, $k>l\ge1$, $Q=p^k$, $P=p^l$, $D\ge1$, $D^{25}\le Q$, and $p\mid D(D-1)$. For the weight $N=3Q+P+D$, **the full conjecture holds at every such pair of exponents whenever $D\equiv1\bmod 2p$.**

More generally, among all eligible offsets $1\le D\le H$ with $p\mid D(D-1)$, the number which could admit an exception for even one pair of exponents is
$$O_p\!\left(H^{\log_p((p+1)/2)}\right). \tag{1}$$
In particular this counts offsets once across both unbounded exponents. The explicit necessary conditions below also test offsets individually. The [binary result](binary_offset_shapes.md) gives a separate bound for $p=2$. At $p=3$ the leading term $3Q$ is itself a prime power. The [two-unit remainder bound](million_part_bound.md) already applies since the two designated units $3Q,P$ leave at most $D$ and $D^5\le27(N-1)^2$.

## The two possible shapes

If $P\le D$, in a counterexample the row avoiding $p$ must have **three equal parts** and all the base-$p$ digits of $P+D$ must be divisible by three. In particular $p\mid D$ is necessary. Indeed $P+D\le2D<Q$. If a part received none of the three leading units $Q$, it would be too small by the [part bound](endpoint_factors.md#lower-bound-on-every-part). Otherwise their spread is at most $2D$ and $(2D)^5<27(N-1)^2$. The [spread bound](row_spreads.md#when-exactly-two-parts-agree) forces equality. Three equal parts adding without a carry have total digits divisible by three, and the units digit cannot be 1.

If $P>D$, put $R=Q/P$. The [two-power theorem](dense_power_offsets.md#independently-varying-exponents) excludes $D<R$. Otherwise put $c=1$ or $2$, $e=2-c$. A counterexample requires
$$D=(3R+1)u,\qquad u\ge1,\qquad
(Ru,\ cRu,\ (eR+1)u)\ \text{adds without a carry in base }p. \tag{2}$$
The actual row is $(RF,cRF,(eR+1)F)$ for $F=P+u$, up to permutation.

In detail, the [nonzero-difference bound](dense_power_offsets.md#the-coefficient-three) and the pure-power exclusion force the row to be $(Q+b,c(Q+b),eQ+P+f)$ with $b>0$. In the simultaneous gcd argument all bounds still apply except when the approximation has both indices 1 and its determinant vanishes. This gives either $b=R(f-eb)$ or $b=R(f-eb-1)$. Write $b=Ru$ and $F=P+u$. The parts are $(RF,cRF,(eR+1)F+\delta)$ with $\delta=0$ or 1. At $\delta=1$ the assigned scale $U\ge3$ divides $RF$ and is coprime to $F$, so $U\mid R$. Now $p\mid U$ implies $F\equiv-1\bmod p$ from the third part. The units digits of $F,cF$ have sum $2p-1-c\ge p$. Multiplying both by the same power $R$ shifts this carry to a higher digit, impossible. Thus $\delta=0$. The low portions have sum $D<P$ and carry freedom gives (2). Since $3R+1$ is even and the equal case requires $p\mid D$, every $D\equiv1\bmod 2p$ is excluded from both cases.

## Counting offsets with a possible row

At $P\le D$, the number of integers $P+D\le2H$ with all digits multiples of three is $O_p(H^{\log_p(1+\lfloor(p-1)/3\rfloor)})$. There are $O_p(\log(2H))$ choices of $P$. The digit alphabet here is strictly smaller than $(p+1)/2$, so this contribution satisfies (1).

In (2) the sum $u+cu$ must be carry-free. For $c=1$ there are $(p+1)/2$ choices per digit of $u$. For $c=2$, given the incoming multiplication carry $v=0$ or 1, a digit $x$ is allowed when
$$x+((2x+v)\bmod p)<p.$$
With no outgoing carry this allows $0\le x\le\lfloor(p-1-v)/3\rfloor$. With outgoing carry one it allows $\lceil(p-v)/2\rceil\le x\le\lfloor(2p-1-v)/3\rfloor$. At $p=6t+1$ the counts in the two intervals are $(2t+1,t)$ for $v=0$ and $(2t,t+1)$ for $v=1$. At $p=6t-1$ they are $(2t,t)$ for both incoming carries. Thus also for $c=2$ there are precisely $(p+1)/2$ choices per digit. For $u\le H/(3R+1)$ this gives $O_p((H/R)^{\log_p((p+1)/2)})$ choices, summed only over $R=p^h$ with $3R+1\le H$. The sum is geometric. This proves (1).

## Reproduction

Run `python3 check_odd_offset_shapes.py` to count the possible offsets and compare the two shapes with direct digit allocations and exact row gcds.

| Prime | Eligible offsets through one million | Equal shape | Scaled shape | Distinct offsets passing either test |
|---:|---:|---:|---:|---:|
| 5 | 400,000 | 833 | 570 | 1,395 |
| 7 | 285,715 | 2,916 | 610 | 3,511 |
| 11 | 181,819 | 2,753 | 546 | 3,278 |

The verifier constructs the two offset lists by the displayed digit conditions and recovers them independently through 5,000 by scanning offsets, divisors and residue sums. At 505 sampled weights it enumerates all 41,919 unordered positive carry-free rows. All 29 passing the size tests have one of the two stated shapes.
