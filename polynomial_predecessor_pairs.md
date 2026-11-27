# Counting possible exceptions along a polynomial relation between predecessor powers

Fix a positive integer $D$ and an integer polynomial $H$ with $H(0)=h>0$ and positive leading coefficient. Consider positive prime powers $P=p^k$ for which $R=H(P)/D$ is a positive power of a different prime. Put $N=1+P H(P)$ and $\beta=\log_2((1+\sqrt5)/2)$.

1. In **every arithmetic progression of $K$ positive integer inputs $P$**, at most $O_{D,H}(K^\beta)$ indices give a failure of the full conjecture with the indicated two prime powers. The constant is independent of the start and step.
2. For **every fixed odd prime $p$**, the same bound holds in every arithmetic progression of $K$ positive **exponents** $k$, with a constant independent of $p$, the start and the step.
3. For $p=2$, if $r$ is an odd prime dividing $D$, the exponent bound is $O_{D,H}(K^{\log_r(5r/6)})$. If $D$ is a power of two, including 1, there are no failures at any exponent.

These count possible failures among all indices. No estimate for how often $H(P)/D$ itself is a prime power is required. For instance the second assertion applies to $P=p^k$, $R=(10P+7)/D$ even in the two surviving residue classes of the [quadratic root test](quadratic_predecessor_roots.md).

There is also a pointwise result for signed polynomials. **If $H$ has a negative coefficient, the two predecessor bases cannot be split between the two rows of an endpoint counterexample at any sufficiently large $P$.** In this case the full conjecture holds for all sufficiently large eligible $P$ whenever $\gcd(D,h)\in\{1,2,3,4,6\}$. All thresholds below are effective. For example, $H(X)=2-X+6X^2$ and $D=119$ give the full conjecture at every eligible $P>821{,}516$. Here the finite maximum in the proof is $1{,}643{,}032$, and the only positive reciprocal root is $t=3$, with quotient $1-X+2X^2$. One actual pair is $P=823{,}663$, $R=34{,}206{,}080{,}687$, both prime, giving $N=3{,}352{,}739{,}681{,}390{,}681{,}240$.

## Recovering a polynomial from a split

Put $F=1+XH$, of degree $n$. Take $P>D(h+1)$. The two positive parts $UPx,UPy$ in the row taking $p$ and $VRz,VRw$ in the row taking the other base satisfy $(x+y)(z+w)<D$. At the layer $P$ above the shift, the low digits in the first row have total $h$. The [predecessor digit argument](predecessor_digit_splits.md) gives positive integers $i,j$ with $i+j\le h$, $i\le j$, $xj=yi$, and $iV=bP+x$ with $0<b<Dh/4$. In particular $b,x$ range over a fixed finite set.

For $F=\sum_{l=0}^n f_l X^l$ form the integer
$$
             J_{b,x}=\sum_{l=0}^n f_l(-x)^l b^{n-l}.
$$
Since $V\mid F(P)$ and $bP\equiv-x\bmod V$, we have $V\mid J_{b,x}$. For $2P$ larger than the maximum absolute value of this finite list, $V>2P$ implies $J_{b,x}=0$. The reciprocal polynomial of $F$ is monic, so $t=b/x$ is a positive integer and $F=(1+tX)G$ for some $G\in\mathbb Z[X]$. This comparison uses signed coefficients throughout. As in the digit argument, $Ux=iG(P)$, $Uy=jG(P)$, and the other part is $(1+(t-i-j)P)G(P)$, giving $t\ge i+j$.

In fact $G$ has nonnegative coefficients once $P$ is sufficiently large. For each of the finitely many possible roots take $P>2h\max_l|G_l|$. If $m$ is the first negative coefficient, all lower coefficients of $iG,jG$ are already ordinary radix-$P$ digits. At position $m$ the digits would be $P-i|G_m|,P-j|G_m|$. Their sum exceeds $P$ since $i+j\le h$. This contradicts the addition without a base-$p$ carry above the shift. Thus $F$ has nonnegative coefficients too, proving the pointwise split obstruction. In the grouped case the relation $Dv=ah$ and the small-denominator exclusions give the stated condition on $\gcd(D,h)$. [Grouped obstruction](quadratic_predecessor_roots.md#grouping-and-the-small-remainders)

## A finite list of pairs of parts

At every eligible exception beyond the threshold, one of the two nonendpoint rows must avoid 2. If $r\mid D$, one of them must avoid $r$ as well. At either prime its avoidance requires a carry-free pair of positive parts (or a part and its complement) from the following list.

| Assignment and row | Scalar $S(P)$ | Pair whose addition must be carry-free |
|---|---|---|
| Grouped, row taking both bases | $P H(P)/D$ | $a S,(D-a)S+1$, $1\le a<D$ |
| Grouped, companion | $F(P)/U$ | $f S,g S$, $3\le U<D$, $f,g>0$, $f+g<U$ |
| Split, row taking $p$ | $P G(P)$ | $i S,j S$ |
| Split, companion | $(x/i)(1+tP)H(P)/D$ | $z S,w S$ |

For the grouped alternatives the part at the intersection of zero positions is $aPR$, and the companion's positive normalized parts sum to $U$. In the split alternatives the previous section determines the actual first pair and gives $V=(x/i)(1+tP)$. All listed scalars are nonconstant rational polynomials in $P$. They are positive integers whenever the corresponding data occur. All slopes, denominators and choices of polynomial are bounded in terms of $D,H$. It remains to count carry-free additions from a fixed line of the table.

## Two summands at an odd prime

For any fixed positive slopes $a,b$ and nonnegative integer offsets, at most $O_{r,a,b}((5r/6)^l)$ strings of $l$ scalar digits after a common prefix can extend to a carry-free pair at an odd prime $r$. The constant is uniform in the offsets and the prefix. For $r=2$ the [two-summand binary estimate](binary_affine_counts.md) gives $O_{a,b}(\varphi^l)$ with $\varphi=(1+\sqrt5)/2$.

For the odd-prime estimate use the two digit streams and the complemented digits of their sum. Carry freedom forces their digits to sum to $r-1$ in every column. After removing the common $r$-adic valuation of $a,b,a+b$, at least two slopes are units. Dividing out the common power in the digit streams only changes the fixed offsets. The digit at delay $v_r(a_i)$ for slope $a_i$ is a permutation of the next input digit once earlier inputs are fixed. At three unit slopes, summing over $m$ admissible next inputs gives $3m(m-1)/2\le m(r-1)$ and hence $m\le(2r+1)/3$.

At two unit slopes let the remaining delay be $d>0$. At each stage, the delayed digits in the next $d$ columns have been determined. Give the prefix weight $z^{-T}$, where $T$ is their sum and $z>1$. If the current delayed digit is $s$, at most $r-s$ choices can pass, since either unit digit must be at most $r-1-s$. The new delayed digits are all distinct. The total relative weight of the children is at most
$$
                z^s\sum_{j=0}^{r-1-s} z^{-j}.
$$
For $r=3$, take $z=4/3$. The three bounds are $37/16,7/3,16/9$. For $r\ge5$, take $z=\exp(1/r)$. With $y=z^s$ the bound is $(y-y^2/e)/(1-e^{-1/r})\le e(r+1)/4<5r/6$. This uses $(1-e^{-1/r})^{-1}\le r+1$. The three-unit bound is also below $5r/6$. Prefix weights lie between $z^{-(r-1)d}$ and one, giving the result by induction. A prescribed low scalar prefix only changes the fixed offsets in the remaining streams.

## From scalar digits to parameter indices

Here are details for the uniformity in the arithmetic progressions. Use $\ell=2$ for integer inputs and odd prime bases, and $\ell=r$ for a binary prime-power input. Write $\lambda=\varphi$ or $5r/6$ accordingly. The preceding argument bounds readable blocks of length $l$ after a common prefix by $C\lambda^l$.

For integer inputs in a progression, the scalar is a nonconstant rational polynomial of the index. The [local polynomial estimate and ball count](rational_exponent_counts.md#how-a-rational-function-acts-on-residue-classes) transfer a uniform block bound to $O(K^{\log_\ell\lambda})$ indices. For the exponential inputs the polynomial acts *after* exponentiation. To include both cases, write $P_j$ for the input progression or powers in a progression of exponents. In the latter case split indices into their two parity classes for $\ell=2$, or their classes modulo $\ell-1$ for odd $\ell$. Reindex each class by consecutive integers. In each class, for some integer $c\ge0$ independent of the two indices,
$$
          v_\ell(P_j-P_{j'})=c+v_\ell(j-j')\qquad(j\ne j').
$$
For an arithmetic progression this follows by subtraction. For the powers it follows by expanding powers of a principal unit, using a square at 2. The phase split ensures that this applies for every positive step of the exponent progression.

Let $S$ be a candidate scalar polynomial and take the finite roots of $S'$ in a splitting field. For an index ball modulo $\ell^d$, all inputs have a common prefix of length $c+d$. If a present input lies within valuation at least $c+d-1$ of a root, subdivide the ball. At depth $d\ge1$ each root can affect at most $\ell$ balls since any two such indices agree modulo $\ell^{d-1}$. On each other ball, the polynomial estimate gives
$$
 v_\ell(S(P_j)-S(P_{j'}))=v_\ell(S'(P_j))+v_\ell(P_j-P_{j'}),
$$
with constant derivative valuation throughout. Thus, in the reindexed coordinates of the ball, output differences have valuation equal to index-difference valuation plus one fixed integer. If the ball contains at most $k=\lceil K/\ell^d\rceil$ consecutive indices in those coordinates, the outputs of interest share a prefix and the next $\lceil\log_\ell k\rceil$ digits distinguish them. For a negative shift, the first $\lceil\log_\ell k\rceil$ digits already distinguish all integral outputs. There are $O(C k^{\log_\ell\lambda})$ allowed indices in the ball.

There are at most $O_\ell(1+\deg S)$ newly accepted balls at each depth. Stop at depth $\lceil\log_\ell K\rceil$ and count the remaining singletons individually. The geometric series over depths gives the asserted bound, with no dependence on $c$ or the starting index. Summing over the finite list of rows and adding the bounded initial prime powers proves all three claims. In the binary-input case with no odd divisor of $D$, at most two distinct primes divide the predecessor, so the [two-prime theorem](all_weights_endpoint_reduction.md) handles every input directly.

Run `python3 check_polynomial_predecessor_pairs.py` for signed and nonnegative polynomial controls, odd-base block counts with offsets, and local distance checks for polynomials evaluated after exponentiation. [Verifier](check_polynomial_predecessor_pairs.py)
