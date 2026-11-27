# A ternary restriction on the exponents of the remaining binary powers

Let $P,R$ be prime powers with $10P+7=DR$ and $1\le D\le10{,}000$. The full three-row conjecture at $(2P+1)(5P+1)$ holds for every $P\le 2^{10^{30}}$. It also holds for every exponent whenever either underlying prime base is at most $2^{3\cdot10^{15}}$.

For any $H\ge2$ consecutive binary orders, meaning intervals $[2^j,2^{j+1})$, **at most $O(H^{5/38})$ individual inputs $P$ in their union can fail**, uniformly in the starting order. In particular there are $O((\log X)^{5/38})$ possible exceptions through $X$, improving the [logarithmic estimate](seven_quadratic_suffixes.md). The exponent $5/38$ is less than 0.132. The assertions allow arbitrary prime bases except where one base is explicitly bounded.

## Binary scale and ternary words

Use the [three residual shapes](seven_quadratic_cofactors.md) and put $T=2P+1$, $V=TR/2$, $A=(2z,2w,D-2z-2w)$. The candidate companion has parts $A_1V,A_2V,A_3V-5T/2$ and must avoid both 2 and 3 since $3\mid D$. If $j=\lfloor\log_2 P\rfloor$, the [binary word theorem](seven_quadratic_prime_bases.md) gives an integer $m=j+h$ and one of six explicit rational numbers $f$ with
$$
             |V-16 f2^m|\le 1040\,2^j .
$$
Here $h$ is the word length, and $j-d-1\le h\le j-d+2$ with $d=\lfloor\log_2 D\rfloor$. Thus $m=2j+O(1)$ and, at each fixed $m$, at most a bounded number of $P$ can occur: only boundedly many $j,h$ and binary words are possible, and each word admits at most one input in its block. For each fixed $f$, the circular distance between $\log_3 V$ and $m\log_3 2+\log_3(16f)$ modulo one is $O(2^{-m/2})$.

Choose the smallest ternary power $L\ge5T/2$. For $v=\lfloor V/L\rfloor$, the upper parts after division by $L$ are $vA+s$, where
$$
 s=(\lfloor A_1 x\rfloor,\lfloor A_2 x\rfloor,\lfloor A_3 x\rfloor-e),
 \qquad 0\le x<1,\quad e\in\{0,1\}.
$$
Accept every such offset. When reading the next trit $b$ of $v$ from the low end, keep the step if the sum of the trits of $s+bA$ is less than 3, then pass to $\lfloor(s+bA)/3\rfloor$ coordinatewise. This again belongs to the offset set, using $(x+b)/3$ and a loss of zero or one. To end an accepted word the remaining offset must be nonnegative and add without carries. The leading digits of $V$ agree with those of $v$ through the entire length of $v$, since $L$ is a power of 3. This gives at least $(m/2)\log_3 2-O(1)$ constrained leading trits.

Keep every offset from which one can reach an acceptable end. Write $M$ for the matrix counting edges between these offsets, including their multiplicities. The vector $v_t=M^t\mathbf1$ counts length-$t$ paths at each starting offset. Componentwise $10^6v_{513}\le b v_{512}$ for the following integer certificates.

| $D$ | Offsets reaching an acceptable end | Nonempty states reading from the high end | $b$ |
|---:|---:|---:|---:|
| 4,587 | 2,498 | 354 | 1,115,844 |
| 7,383 | 3,304 | 654 | 1,107,472 |
| 9,453 | 3,550 | 759 | 1,155,459 |

Iterating the vector inequality bounds the number of *any* admissible internal words of length $t$ by $C(b/10^6)^t$. A word in the middle of an accepted path starts at one of the retained offsets and stays in the graph. The constant absorbs small $t$ and all possible starting offsets. In particular each fixed prefix has at most $C3^{\delta t}$ extensions by $t$ trits for $\delta=5/38$, since $b^{38}<3^5(10^6)^{38}$ in each case. [Integer path certificates](check_seven_ternary_words.py)

## Counting indices even between consecutive convergent denominators

A ternary prefix of length $k$ specifies an interval of normalized mantissas of length $3^{1-k}$ in $[1,3)$. Applying $\log_3$ gives an interval on the circle of length comparable to $3^{-k}$, with absolute upper and lower constants. An arc of length comparable to $3^{-l}$ therefore meets $O(1)$ of the full length-$l$ grid intervals. Inside it the eligible set at depth $k\ge l$ is covered by $O(3^{\delta(k-l)})$ length-$k$ intervals, by the extension bound. The global count is $O(3^{\delta k})$.

Write $\alpha=\log_3 2$. Consider $J$ consecutive indices $m$ starting at $m_0\ge4J$, for a fixed limit $f$. Take a convergent $a/q$ to $\alpha$ with $q\le J$ and the next denominator $q'>J$. The standard convergent inequality gives $\sigma=|q\alpha-a|<1/q'$. Also $\sigma\gg 2^{-q}$, since the unequal integers $2^q,3^a$ differ by at least one and $3^a\asymp2^q$. Thus $O(2^{-m_0/2})$ is $O(\sigma)$, and every depth up to $\lceil\log_3(1/\sigma)\rceil+O(1)$ is available at these indices. Finitely many small $J,m_0$ are absorbed into the constant.

Group the indices by $m-m_0$ modulo $q$. The phases in a group lie within $1/q$ of a translated point of the $q$-grid because $J|\alpha-a/q|<1/q$. If the group contains a possible failure this grid point is within $O(1/q)$ of an admissible leading interval at depth $\lfloor\log_3 q\rfloor$. The global count and the grid spacing leave $O(q^\delta)$ groups.

Within a fixed group the phases progress monotonically in steps of $\sigma$ before reduction modulo one. They are contained in an arc of length at most $(1+J/q)\sigma$, splitting at the circle boundary if needed. Enlarge it by the $O(2^{-m_0/2})$ error. Use prefix scales comparable to $(1+J/q)\sigma$ and $\sigma$ in the local extension bound, using the whole circle if the larger scale exceeds one. At most $O((1+J/q)^\delta)$ fine arcs are needed, each meeting $O(1)$ group elements even after enlargement by the phase error. Multiplication gives $O(J^\delta)$ possible indices. Both the grid of groups and the spacing within each group are used here.

For an arbitrary block $[m_0,m_0+J)$ with $m_0<4J$, cover the range up to $5J$ by dyadic index intervals and split each into a bounded number of blocks satisfying the starting condition. The resulting geometric sum is $O(J^\delta)$. Finally $H$ consecutive input orders give a block of $2H+O(1)$ possible indices $m$. Sum over the six $f$ and use bounded input multiplicity at fixed $m$.

## A finite search on binary exponents

For $2^{65536}<P\le2^{10^{30}}$ we have $65536\le m\le2\cdot10^{30}$. The earlier theorem handles the smaller inputs. The ternary word length exceeds 40,000: $L<15T/2$ gives $V/L>R/15>P/15000$. The relative error in using $16f2^m$ for $V$ is less than $2^{-65515}$, from $V\ge10P^2/D$ and the binary scale inequality. This is smaller than one unit at every fixed precision used below, including the few units needed to convert to a logarithm or a normalized mantissa.

The first search certifies $|\alpha-a/q|\le1/q^2$ by integer logarithm enclosures, with
$$
 a=2180796053156940756896192173706,\qquad
 q=3456479965974452268626125476129>2\cdot10^{30}.
$$
Every eligible phase lies within $2/q$ of $m a/q+\log_3(16f)$ modulo one. Enumerate the accepted 66-trit leading words. For each associated logarithmic interval, enumerate the grid integers it can meet with that margin and invert $a$ modulo $q$. This finds every eligible exponent in the desired range. For each candidate exponent the program computes the normalized mantissa of $16f2^m$ by repeated squaring, with outward integer rounding at precision $3^{220}$. It includes the projection error and tests all possible prefixes of 125 trits. No exponent survives.

| $D$ | First-cut words | Grid candidates for each $f$, in the order shown in the binary word theorem | Final survivors |
|---:|---:|---|---:|
| 4,587 | 4,756 | 10,619 | 0 |
| 7,383 | 2,371 | 5,913, 5,934, 5,935, 5,938 | 0 |
| 9,453 | 72,799 | 156,767 | 0 |

An independent program generates high-end transitions by inverting all possible output columns. It brackets $\alpha$ and the offsets with directed binary squaring, brackets the logarithmic interval endpoints by a rational Taylor series, and counts the integer points at each prefix using Euclidean floor sums. With slope $A/S$ and integer interval endpoints $l,u$ broadened to cover every rounding error, each index contributes
$$
        \left\lfloor(Ai+b-l)/S\right\rfloor-
        \left\lfloor(Ai+b-u)/S\right\rfloor .
$$
Zero total excludes all descendants. This search uses the full interval of indices directly. It exhausts the six cases at prefix lengths 72, 71, 73, 73, 74 and 79 in table order.

For the base bound, an exponent at least $3\cdot10^{14}$ is already settled by the [uniform exponent theorem](seven_quadratic_prime_bases.md). For smaller exponents a base at most $2^{3\cdot10^{15}}$ gives the corresponding power at most $2^{9\cdot10^{29}}$. If that power is $R$, then $P<1000R$. In either case the input lies in the finite range just proved.

```sh
python3 check_seven_ternary_words.py
python3 check_seven_ternary_rotation.py 30
python3 check_seven_ternary_floor_sums.py 30
```
