# Exceptional weights globally and in moving intervals

Let $E(I)$ count weights $N\ge3$ in $I$ for which three triples of positive parts summing to $N$ have trinomial coefficients with gcd one.

1. $E([3,X])=O(X^{0.9329})$, in particular $O(X^{14/15})$. This improves the [earlier exponent of $15/16$](endpoint_count.md).
2. There is an absolute constant $C$ such that every interval $I$ of $D\ge X^{19/25}$ consecutive integers contained in $[X,2X]$, for $X\ge 3$, satisfies
   $$
   E(I)\le C D X^{-10^{-8}}.
   $$
   The bound is uniform in the location of the interval, including the binary digits of its starting point.

## The congruence determined by a part

Use a dyadic block $2^{L-1}\le N<2^L=Q$. By the [endpoint reduction](all_weights_endpoint_reduction.md), a counterexample includes the row $(1,1,N-2)$ with coefficient $N(N-1)$. Choose one of the other rows that avoids 2 and call it $A$. Assign 2 to $A$ and assign every remaining prime of $N(N-1)$ to an avoiding row among the two. Choose witnesses deterministically if several work. Write $UV=N$ and $ST=N-1$, where $U,S$ are the products of the respective full prime powers assigned to $A$. Let $K=US$.

Each part $z$ of $A$ satisfies $U\mid z$ and $S\mid z(z-1)$ by the carry criterion at the assigned powers. The other row consists of three positive multiples of $V$, so $z\ge U\ge3$. Therefore any one part specifies at most $Q^{o(1)}$ divisor pairs $(U,S)$, uniformly in $z$. Each pair specifies a single class for $N$ modulo $K$ by the Chinese remainder theorem. We also use the [small-modulus count](endpoint_count.md#counting-assignments-with-small-k), which shows that at most $Q^{-1+o(1)}Y^{5/2}$ weights in the block admit an assignment with $K\le Y$. The divisor estimates in both statements are uniform.

## Recovering the weight from either end

Put $k=773/1000$, $l=\lfloor kL\rfloor$ and $h=L-l$. The contribution with $K<2^l$ is at most $Q^{0.9325+o(1)}$. Otherwise put $e=v_2(N(N-1))=v_2(K)$. Write $N=2^e n+\eta$, where $\eta\in\{0,1\}$, and represent $n$ by $L-e$ bits. The low $e$ bits of the three parts of $A$ are zero, except for one units bit if $\eta=1$, since their ones partition the ones of $N$. Fix $e,\eta$ and the index of that part when present. These choices cost $O(L)$. If $L-e\le2h$, even the unrestricted number of choices for $n$ is at most $2^{2h}$.

Suppose $L-e>2h$. Given a part and one of its candidate divisor pairs, $n$ lies in a single class modulo the odd integer $M=K/2^e$. Either end of the word recovers $n$.

- Its top $h$ bits put it in an interval of $2^{L-e-h}\le M$ consecutive integers.
- Its bottom $h$ bits specify a second class modulo $2^h$. The product $2^h M\ge2^{L-e}$ is at least the size of the entire range.

Thus **any one part**, together with **either end** of the word, specifies at most $Q^{o(1)}$ possible weights. The oddness of $M$ gives the coprimality in the second argument.

## Counting with both ends

Split $n$ into blocks of lengths $h,m,h$, where $m=L-e-2h$. Fix their numbers of ones, with respective densities $r_1,s,r_2$. The number of such strata is polynomial in $L$. Write $H(t)=-t\log_2 t-(1-t)\log_2(1-t)$, with its usual endpoint values, and $B=H(1/3)$. For the number $W>0$ of exceptional weights in a stratum, the binomial count gives
$$
\log_2 W\le hH(r_1)+mH(s)+hH(r_2).
$$

Choose a uniform weight in the stratum. Let $Y$ be its first block, and let $Z_i$ be the $i$-th part of its chosen row. For each $i$, the conditional entropy of the weight given $Y,Z_i$ is $o(L)$ by recovery. The entropy chain rule bounds $\log_2 W$ by the entropy of $Y$ plus the average of the three conditional entropies of $Z_i$ given $Y$, plus $o(L)$.

The entropy of $Y$ is at most $hH(r_1)$. On its $h r_1$ occupied positions, conditional on a specified $Y$, the three probabilities that the respective parts occupy a given position sum to one. By concavity their average entropy is at most $B$. All other positions of $Y$ cost zero. On the middle block the occupation probability averaged over parts and positions is $s/3$, even conditional on $Y$, since the number of ones in that block is fixed. Subadditivity and concavity give average cost at most $mH(s/3)$. The last block likewise costs at most $hH(r_2/3)$. The low $e$ bits of each part are fixed. This gives the second estimate
$$
\log_2 W\le h[H(r_1)+r_1 B+H(r_2/3)]+mH(s/3)+o(L).
$$
The same estimate holds with $r_1,r_2$ interchanged. Average the two and put $r=(r_1+r_2)/2$. Take a convex combination of this average with the binomial estimate, giving the average weight $\lambda=79/125$. By concavity,
$$
\log_2 W\le
 h\big[(2-\lambda)H(r)+\lambda(Br+H(r/3))\big]
 +m\big[(1-\lambda)H(s)+\lambda H(s/3)\big]+o(L).
$$

Use $h=(1-k)L+O(1)$ and $m\le(2k-1)L$. The resulting coefficient of $L$ is strictly less than $0.9329$ on the whole unit square. Here are rational parameters to verify the comparison. Replace each occurrence of $H(t)$ by its affine majorant $-t\log_2 p-(1-t)\log_2(1-p)$, tangent at $p$. For $r,r/3,s,s/3$ use respectively $r_0,r_0/3,s_0,s_0/3$, with $r_0=62256/100000$, $s_0=67094/100000$. It suffices to test the four corners of the resulting affine expression. The [exact arithmetic verifier](check_binary_windows.py) bounds all logarithms by rational series with explicit tails. The strict margins absorb the uniform $o(L)$ and the polynomially many strata. Summing the dyadic estimates proves the first assertion.

## A bound uniform in the binary prefix

Set $\theta=19/25$ and $a=7039/10000$. Fix the high $f=L-n$ bits of the weight arbitrarily, where $n=\lfloor\theta L\rfloor$, leaving an aligned block of $2^n$ weights in the dyadic range. Let $l=\lfloor aL\rfloor$, $h=n-l$. The small-modulus estimate for $K<2^l$ bounds the count in *each* such block by $Q^{0.75975+o(1)}$, using the count for the whole dyadic range.

For $K\ge2^l$, reveal the high $h$ bits of the $n$ free bits. Write $j$ for the number of ones revealed and $u$ for the ones in the remaining $l$ bits. Denote by $F\le f$ the ones in the fixed prefix. A given part $z$ and its divisor pair now identify the weight because the unknown suffix ranges over $2^l\le K$ consecutive integers.

Give each bit of a part in a known position cost 1 and each bit in the unknown suffix cost $\log_2 5$ when it equals one. Some part has cost at most $(F+j+u\log_2 5)/3$. Counting all candidate parts with this cost bound by weighting each subset by 2 to the negative of its cost gives at most
$$
2^{(F+j+u\log_2 5)/3}(3/2)^{F+j}(6/5)^l
$$
candidates. Indeed, only known ones may be occupied above the suffix. Each candidate part allows $Q^{o(1)}$ weights for the revealed bits.

Put $d_h=1/1000$, $d_l=1/5000$. For $j\le(1/2+d_h)h$ and $u\le(1/2+d_l)l$, sum over at most $2^h$ revealed strings and use the uniform cost bound. With $C_0=H(1/6)$, the logarithm of the count divided by $L$ is at most
$$
\theta-a+B[1-\theta+(1/2+d_h)(\theta-a)]
+a[C_0+(d_l/3)\log_2 5]+o(1)<0.75997+o(1).
$$
Here $B=1/3+\log_2(3/2)$ and $C_0=(\log_2 5)/6+\log_2(6/5)$.

The number of strings in the two omitted tails is at most
$$
2^n\left(2^{-2d_h^2 h/\ln 2}+2^{-2d_l^2 l/\ln 2}\right).
$$
This follows from the binary upper-tail estimate $\sum_{j\ge tv}\binom vj\le 2^{vH(t)}$ for $t\ge1/2$, and $H(1/2+d)\le 1-2d^2/\ln2$ by the second derivative. Both savings in the exponent divided by $L$ strictly exceed $10^{-8}$ for all sufficiently large $L$. We obtain $O(Q^{\theta-10^{-8}})$ exceptions per aligned block with an absolute implied constant.

An interval in $[X,2X]$ meets at most two dyadic ranges of integers. Their scales $Q$ are comparable to $X$. Within each range, cover its intersection by all aligned blocks it touches. The number is $O(1+D/X^\theta)$, since each has length comparable to $X^\theta$. Summing and using $D\ge X^\theta$ proves the second assertion.

Run `python3 check_binary_windows.py` for the constant bounds and independent enumeration of both reconstruction statements on binary words of up to eleven bits after stripping the forced end. The existence of any exceptional weight remains open.
