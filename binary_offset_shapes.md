# Lower offsets that work at every secondary binary exponent

Let $D\ge1$, $k>l\ge1$ and $D^{25}\le2^k$. If Wasserman's conjecture fails at
$$N=3\cdot2^k+2^l+D,$$
put $h=k-l$, $R=2^h$. There must exist a positive integer $u$ such that
$$D=(3R+1)u,\qquad u,\ Ru,\ 2Ru\ \text{have pairwise disjoint binary supports}. \tag{1}$$
In that case the row avoiding 2 has to be a permutation of $(F,RF,2RF)$ with $F=2^l+u$. In particular **every offset whose binary digit sum is not divisible by three works for every $l$**. For other offsets the divisibility and disjointness in (1) are a direct test. The exceptional offsets allowed by the test are counted only once even if several exponent gaps are possible.

For all but $O(H^{\log_2\rho})$ offsets $1\le D\le H$, every weight above satisfies the conjecture for every pair of exponents in the stated range. Here $\rho$ is the positive root of $\rho^3=\rho^2+1$, so $\log_2\rho<0.552$. If one restricts to exponent gaps $h\ge2$, the bound improves to $O(H^{0.53})$. The exponent $\log_2\rho$ is sharp for the set passing (1). The counting argument gives no lower bound on the number of actual counterexamples.

**At every $1\le D\le256$ the same shape theorem holds for all $k>l\ge1$, without the size hypothesis.** In particular 225 of these offsets give the full conjecture for every pair of exponents.

## Proof of the binary shape

Write $Q=2^k$ and $P=2^l$. If $D<R$ the [independent-exponent theorem](dense_power_offsets.md#independently-varying-exponents) already excludes a counterexample. Otherwise $R\le D$ gives $P=Q/R\ge D^{24}>D$, since $R\ge2$. In the proof of that theorem the two upper bits $Q,2Q$ must go separately into the two parts missing $P$: if either part lacked a high unit it would be at most $D$. Thus in its notation $a=1$, $c=2$, $e=0$. The nonzero-$L$ bound and the $b=0$ exclusion do not use $D<R$. Therefore the row has to be $(Q+b,2Q+2b,P+f)$ with $b>0$, $3b+f=D$.

In (5) of the [gcd argument](dense_power_offsets.md#a-gcd-bound-for-two-shifted-powers), the only failure of its upper bounds without $D<R$ is an exact equality for the approximant with both indices 1. That equality is $b=Rs$, where $s=f$ or $f-1$, so $s$ is positive. Taking $u=s$ and $F=P+u$, the row is $(RF,2RF,F+\delta)$ for $\delta=0$ or 1. If $\delta=1$, its weight scale $U\ge3$ divides $F+1$ and $RF$, hence divides $R$. In particular $4\mid U$ and $F\equiv3\bmod4$. The two lowest set bits of $F$ are consecutive, so $RF$ and $2RF$ have a common bit, contradicting carry freedom. Thus $\delta=0$. The lower portions $u,Ru,2Ru$ sum to $D$ without carries and give (1). Their binary digit sums are equal.

## Counting the offsets

For fixed $h$, condition (1) says that the set bits of $u$ have no two positions differing by $1,h,h+1$. When $h=1$ these are differences 1 and 2, and $D=7u$. For strings of length $n$ the count $a_n$ including zero satisfies $a_0=1$, $a_1=2$, $a_2=3$, $a_n=a_{n-1}+a_{n-3}$. Split by the last digit, which if 1 requires the two previous digits to be zero. Induction using the equation for $\rho$ gives $a_n=\Theta(\rho^n)$. Thus the $h=1$ candidate list alone has $\Theta(H^{\log_2\rho})$ members.

For all $h\ge2$ put $n=\lceil\log_2(H+1)\rceil$. Since $D=(3\cdot2^h+1)u$, the $n$-bit word for $u$ has at least $h+1$ top zeros. Its set bits still avoid differences $1,h,h+1$ when positions are taken cyclically modulo $n$. For $2\le h\le11$, partition the ordinary word into blocks of 128 bits and a final shorter block. The largest number of allowed 128-bit blocks over these ten gaps is $13{,}289{,}854{,}365{,}121{,}163{,}790$, less than $2^{128/2}$. The recurrence retains the last $h+1$ bits and appends either 0 or a 1 allowed by the three distances. This gives $O(2^{n/2})$ possibilities at these fixed gaps.

For $h\ge12$ and $n>4h+11$, use all cyclic translates of the 60 positions
$$\{ah+b:0\le a<5,\ 0\le b<12\}.$$
They are distinct modulo $n$. Within each of the five rows the occupied bits cannot be adjacent. From one row to the next, a bit forbids the bit above it and the one diagonally above it at distance $h+1$. This permits at most $A=3{,}396{,}136{,}349$ patterns. To compute $A$, start with count one for every 12-bit mask with no adjacent ones. Four times, replace each mask's count by the sum of its predecessors' counts, allowing the upper mask $v$ when $v$ is disjoint from both the lower mask $w$ and $2w$. Sum the final counts.

Each cyclic position appears in exactly 60 translates. The number of possible words is at most $A^{n/60}$. For completeness, take a uniformly random valid word and expand its base-two entropy into successive conditional entropies. The entropy of the restriction to any translate is at least the sum of the corresponding conditional entropies, since conditioning on fewer preceding bits can only increase entropy. Summing over translates gives 60 times the entropy of the entire word as a lower bound. Each translate has entropy at most $\log_2 A$.

In the remaining case $n\le4h+11$ the potentially nonzero part of the word has length at most $n-h\le3n/4+O(1)$. The number of words with no adjacent ones is $O(\varphi^{3n/4})$, where $\varphi=(1+\sqrt5)/2$. Each of the three bounds has base-two exponent strictly below $0.53$: for the grid $A^{100}<2^{3180}$, and $\log_2\varphi<0.7$. Summing over at most $n$ choices of $h$ gives $O(H^{0.53})$. Since $\log_2\rho>0.53$, adding the first gap proves the assertions.

## Reproduction

Run `python3 check_binary_offset_shapes.py` for the integer counts, a direct comparison of the offset lists through one million, and a weight-by-weight audit at several exponents.

Through one million only 2,988 offsets meet (1) for any gap, including 1,840 at gaps of at least two. Thus the theorem covers the other 997,012 offsets for all secondary exponents. The two independent lists agree entry by entry. A separate sample of 3,460 weights tests all 2,659,304 positive unordered binary carry-free rows. All 46 passing the two size tests have the stated unshifted form.

## All exponents at the first 256 offsets

For $D\le256$ it remains to check $2^k<D^{25}$ and $N>10^9$, since [all smaller weights are covered](full_range_1b.md). Here $(2D+2)^7<27(N-1)^2$ and $2D<Q$. If $P\le D$, one binary part is at most $2D$. Otherwise the three upper units must occupy separate parts by the individual-part bound. A nonzero $L=2b-d$ is impossible by the bound $27(N-1)^2\le (|L|+2)^7$ in the [coefficient-three proof](dense_power_offsets.md#the-coefficient-three). The same pure-power argument eliminates $b=0$. Thus enumerate the lower allocations $(b,2b,f)$ with $3b+f=D$, $b>0$ and pairwise disjoint supports. For each actual row outside (1) check its gcd and the three predecessor bins directly.

Run `python3 check_small_binary_offsets.py`. There were 3,476,636 weight and exponent choices, 3,254,974 with three separate upper bits and 10,775,442 lower allocations with zero difference. Of these 4,140 had the permitted form and every other allocation failed the size test. A direct enumeration of all 153,373 positive binary-free rows at 712 control weights recovered exactly the 20 predicted rows passing the test. The offsets through 256 admitting (1) are

~~~text
7 13 14 25 26 28 49 50 52 56 63 97 98 100 104 112 119 125 126
193 194 196 200 208 221 224 231 238 245 250 252
~~~
