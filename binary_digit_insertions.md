# Stability under inserting a growing number of binary gaps

Write an integer below $2^n$ as a word of exactly $n$ binary digits, allowing leading zeros. Put
$$r=\lfloor\sqrt{\log_2 n}\rfloor,\qquad h=\lfloor n(\log_2 n)^{1/4}\rfloor.$$

**Theorem.** With at most $O(2^{\,n-n/(20\sqrt{\log_2 n})})$ exceptional input words, the following holds simultaneously for every allowed insertion. Cut the word into at most $r$ nonempty contiguous blocks, keeping its least significant one in the lowest block. Insert zero strings of arbitrary lengths between blocks, with final length at most $h$. Every resulting weight satisfies Wasserman's three-row conjecture.

One may also permute, reverse, or complement the blocks above the lowest one. In the gaps between blocks or above the highest block one may insert a total of at most $r$ runs of ones, of combined length at most $n/100$, in addition to zeros. The same exceptional-set bound covers all these choices together. If only odd outputs are considered, even the lowest block may be moved and changed, and insertions may start at the units digit. These unrestricted rearrangements also work for every even output with $v_2(N)\le n/6$. All statements about weights concern $N\ge3$.

For example, with probability tending to one, an $n$-bit input remains good for *every* placement of $\lfloor\sqrt{\log_2 n}\rfloor-1$ zero gaps above its lowest one, even when the inserted zeros together are longer than the original word. A version for any fixed $r\ge1$ and any fixed $A\ge1$ permits output length $h\le An$ and has at most $O_{r,A}(2^{(1-1/(12r+33))n})$ exceptional words.

## Counting possible outputs at fixed positions

Put $V=2^n$, $\epsilon=1/(12r+32)$. Write $h$ for the allowed length. All estimates below apply both to the growing parameters and to the fixed-parameter statement. We use a factor
$$E=2^{O(h/\log(2h)+r\log(2rh))}$$
with absolute implied constants. Products of a bounded number of such factors have the same form. A positive integer of at most $3h$ bits has at most $2^{O(h/\log(2h))}$ divisors. To see this, factors at primes up to $\sqrt h$ contribute at most $(3h+1)^{\sqrt h}$ choices, while the sum of the exponents at larger primes is $O(h/\log(2h))$.

For the assertion with a fixed lowest block, discard $m=0,1$ and inputs with $2^{v_2(m)}>V^\epsilon$. This costs $O(1+V^{1-\epsilon})$ words. An even output in that assertion retains the location of the lowest one and thus $2^{v_2(N)}\le V^\epsilon$. Odd outputs need no restriction on the lowest bits.

Fix the original cut, the chosen permutations and orientations, the positions of the moved blocks and the inserted ones. There are at most $h^{O(r)}(4^r r!)$ descriptions, absorbed in $E$. At each description the map from input words to output weights is injective. It gives disjoint variable binary fields of widths $w_i\ge1$ with $\sum w_i=n$, starting bit positions $a_i$ (counted from zero), and a fixed background $B$ of at most $n/100$ ones outside the fields. In increasing position order, $a_{i+1}\ge a_i+w_i$. Put $L_i=2^{w_i}$. Before enforcing the restrictions on the lowest input block, all field values vary freely, with $N=B+\sum 2^{a_i}x_i$, $0\le x_i<L_i$.

At a counterexample take the binary carry-free row $A$ paired with the endpoint and its complementary row. Let $U\mid N,S\mid N-1$ be the full prime-power products assigned to $A$, and $K=US$. Every part $y$ of $A$ is at least three, $U\mid y$ and $S\mid y(y-1)$. Thus $y$ determines the pair $(U,S)$ up to $E$ choices, and the pair determines $N\bmod K$. These consequences of Kummer's criterion are proved in the [endpoint reduction](all_weights_endpoint_reduction.md) and the [assigned-product bounds](projected_rows.md).

Write $K_o$ for the odd part of $K$ and put $Y=V^{4/5}$. We first count the outputs with $K_o\le Y$. At even outputs with $v_2(N)\le n/6$ we have $K\le YV^{1/6}\le V^{1-\epsilon}$. Applied to the complementary row, the [two-congruence reconstruction](two_predecessor_congruences.md#recovering-the-row-from-a-product) counts $O(V^{1-\epsilon}E)$ such weights up to $2^h$. In detail, sum its dyadic estimates (4) and (5), taking the smaller in each box. The losses consist of a fixed number of divisor counts on integers of $O(h)$ bits and polynomially many dyadic boxes. For intervals entirely below $V^{1-\epsilon}$ simply count their integers.

For an odd output fix $P=2^{v_2(N-1)}$, one of at most $h$ choices, and put $t=(N-1)/P$. The binary carry-free row receives the full power $P$, so $K_o=U(S/P)$. For $t\le V^{1-\epsilon}$ count $t$ directly. On each dyadic interval $T\le t<2T$ above that threshold the [one-bin reconstruction with the scale removed](growing_prime_power_steps.md#counting-with-the-large-scale-removed) gives
$$ O(T^{-2/3} Y^2 E) \ll V^{8/5-(2/3)(1-\epsilon)} E \le V^{1-\epsilon}E. $$
Here $Y<T$ for sufficiently large $V$ (including a constant factor at the first dyadic interval), all dividends in the reconstruction have $O(h)$ bits, and the last bound uses $\epsilon\le1/44<1/25$. The sum over intervals and $P$ is absorbed in $E$. Crucially, this estimate allows $N-1$ to contain a power of two of any size in the output range.

For greater $K_o$ use the lattice
$$\Lambda=\{z\in\mathbb Z^d:\ \sum 2^{a_i}z_i\equiv0\bmod K_o\},\qquad
\|z\|=\max_i |z_i|/L_i,$$
where $d\le r$ is the actual number of fields. Write $\lambda$ for its minimum nonzero length. A label $(U,S)$ restricts the weights to one translated lattice in the field box.

Outside at most $V^{H(3/5)}$ field vectors, the fields contain at most $3n/5$ ones, where $H(t)=-t\log_2 t-(1-t)\log_2(1-t)$. Some part of the carry-free row receives at most $n/5$ of these variable ones. It therefore has at most $2^{n/100}V^{H(1/5)}$ possible values, including its choice of fixed ones. Both estimates follow from the binomial tail bound. Retain only these typical field vectors for the remainder of the argument. Their total number of possible odd labels is at most $V^{H(1/5)+1/100}E$. If $\lambda>V^{-3\epsilon}$ there are at most $3^r V^{3r\epsilon}$ points per translate, by subdividing the unit box at scale $\lambda/2$ (for $\lambda>1$ there is at most one point). Since
$$H(3/5)<39/40,\qquad H(1/5)<29/40,\qquad
 (3r+1)\epsilon<1/4,$$
this counts at most $O(V^{1-\epsilon}E)$ weights in this case as well.

## Very short lattice vectors

Take dyadic scales $\rho\le\lambda<2\rho$ with $\lambda\le V^{-3\epsilon}$. There are $O(n+1)$ scales since $\lambda\ge2^{-n}$. Let $I=\{i:2\rho L_i\ge1\}$, $k=|I|\ge1$ and $V_I=\prod_{i\in I}L_i$. Every shortest vector is supported on $I$, with at most $6^k\rho^k V_I$ choices. For each such nonzero vector $z$ the integer $W=\sum 2^{a_i}z_i$ is nonzero. Indeed the first nonzero coordinate has absolute value below its field width $L_i$ for sufficiently large $n$, since $2\rho<1$. Its term cannot be cancelled by higher fields, by divisibility by $2^{a_i+w_i}$. Now $K_o\mid W$, and $|W|\le r2^h$. This gives at most $E$ moduli and CRT labels per vector. For a fixed odd modulus distribute its prime powers to $U$ or $S$, giving at most its divisor count in labels modulo $K_o$.

Fix the coordinates outside $I$. In the remaining box each translate contains at most
$$\exp(O(k\log(2k)))\big(\rho^{-(k-1)}+V_I/K_o\big) \tag{1}$$
points. If the affine dimension of these points is $t<k$, center disjoint Euclidean balls of radius $\rho/2$ in their affine hull after scaling the box to unit sides. The centers have Euclidean separation at least $\rho$ and diameter at most $\sqrt k$. The ratio of volumes gives $(1+2\sqrt k/\rho)^t$, which implies (1). In full dimension select one point and $k$ differences from it of maximal absolute determinant. Every other difference from that point has all coordinates of absolute value at most one in this basis, by Cramer's rule. The kernel lattice on $I$ has index $K_o$, since every coefficient $2^{a_i}$ is invertible modulo $K_o$. The integral differences therefore span a sublattice of relative index at most $k!V_I/K_o$. Each coset has at most $3^k$ points in the basis box. This proves (1) uniformly in the dimension.

For the first term of (1) count labels using the shortest vectors. Multiplying by the vector count and $V/V_I$ outside-coordinate choices gives at most $VE\rho\le V^{1-3\epsilon}E$ at each scale. For the volume term use the count of typical labels from the sparse part instead. For each label the total volume term over the outside coordinates is at most $EV/K_o$, so summing it once per label gives
$$ V^{H(1/5)+1/100}E\, V/Y\le V^{1-\epsilon}E. $$
One may assign each modulus its least dyadic scale and one shortest vector to avoid any duplication in this second count. This separate estimate of the volume term is what allows the threshold $Y=V^{4/5}$. The factor $E$ also absorbs all scales and descriptions. In the growing case $\log_2 E= o(n/\sqrt{\log_2 n})$, while $\epsilon\sqrt{\log_2 n}\to1/12$. In the fixed case $\log_2 E=o(n)$. This proves both asserted exceptional-set bounds.

The lowest-one condition in the first assertion supplies the required even-output valuation bound simultaneously at all its even outputs.

Run the standard-library checks of the binary operations, valuation, counting constants and the small-vector argument.

~~~sh
python3 check_binary_insertions.py
~~~
