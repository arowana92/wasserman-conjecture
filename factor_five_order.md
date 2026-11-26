# Multiplicative order and binary carries at a factor-five predecessor

Let $Q$ be an odd prime power, $k\ge1$, and $N=5\cdot2^kQ+1$. Each of the following hypotheses proves the full trinomial conjecture at $N$, with **no bound on $k$ or $Q$**.

* There are $d,h\ge1$ with $dQ=2^h+1$ and either $d\le2^{52}$ or $\operatorname{rad}(d)\le400000$. Here $\operatorname{rad}(d)$ is the product of its distinct prime factors.
* There are odd $d\ge1$ and positive $c,h$ with $dQ=2^h+c$, $5c\le d$, and either $\operatorname{ord}_d(2)\le h$ or $h\ge\min(15d^{3/5},164\sqrt d)$.

For every fixed finite set of primes supporting $d$, only finitely many $Q$ can be exceptional in the first construction, even allowing arbitrary powers of these primes in $d$. The proof gives an effective bound. The [independent exponent calculation](factor_five_exponents.md) deals with a range of $Q$ directly.

## The necessary binary split

The [factor-five reduction](predecessor_five_seven.md#full-factor-five) gives $Q=\ell^e$ with $\ell>5$ prime, $\ell\equiv3\bmod8$, $e$ odd, and $k\ge3$ in any counterexample. One nonendpoint row is
$$
(3QP,\ zP,\ (2Q-z)P+1),\qquad P=2^k .
$$
It avoids $2,3,5,\ell$. The other row has three equal parts $N/3$. The two positive parts after $3QP$ are multiples of 15. In particular $15\mid z$, $zP\equiv1\bmod Q$, and $PQ\equiv7\bmod15$. The bits of $z$ and $2Q-z$ partition the bits of $2Q$, and $3Q$ has no bit in common with $2Q$.

**Small inverse powers.** There can be no $j\ge0$ with $2^j<2Q$ and $2^{-k}\equiv\pm2^j\bmod Q$. The binary split makes $z$ even, so it is the unique even lift of its residue in $(0,2Q)$. For $j\ge1$, either $z=2^j$, contrary to $15\mid z$, or $2Q-z=2^j$, which would require $2^{k+j}\equiv-1\bmod15$, whereas powers of 2 occupy only $1,2,4,8$ modulo 15. For $j=0$, the complementary binary pieces are $Q-1,Q+1$. After division by two these are consecutive integers with disjoint binary digits. The larger must be a power of two, so $Q=2^s-1$. This contradicts $Q>5$, $Q\equiv3\bmod8$.

The case $d=1$ in $dQ=2^h+1$ is excluded by $Q>5$ and $Q\equiv3\bmod8$. Also $3Q=2^h+1$ is impossible. The identity $2^h\equiv-1\bmod Q$ represents every power of 2 modulo $Q$, including $2^{-k}$, as $\pm2^j$ with $0\le j<h$, and $2^{h-1}<2Q$.

## A carry at the order modulus

Let $p\ge2$ and $h\ge1$ be integers. Let $a,b,c,d,q$ be positive integers with $d>1$ satisfying
$$
dq=p^h+c,\qquad (a+b)c\le d,\qquad \gcd(p,d)=1.
$$
If $aq,bq$ add without carries in base $p$, then $\operatorname{ord}_d(p)>h$.

In fact, suppose $1\le j\le h$ and $d\mid p^j-1$. Put $m=p^j$ and $s=(m-1)/d$. Modulo $m$ we have $q\equiv-cs$, by multiplying $dq\equiv c$ by $s$. The two positive numbers $acs,bcs$ sum to at most $ds=m-1$. The least residues of $aq,bq$ therefore sum to $2m-(a+b)cs>m$. This forces a carry among the first $j$ digits. Lean 4.12.0 checks the residue inequality with the weaker conclusion $\ge m$ and the implication ruling out every $j\le h$ directly from the no-carry inequalities. The formal proof uses only Lean's standard distribution. [Lean proof](CarryOrder.lean)

For $p=2$, $q=Q$ and $(a,b)=(3,2)$ this proves the order test above. If $c=1$ and $d>3$ in a possible exception then $h$ is odd, and $d\equiv3\bmod24$, hence $d\ge27$. Indeed $2$ is a quadratic nonresidue at $\ell\equiv3\bmod8$, so its order at $\ell^e$ is even with exactly one factor of two. Its half-order gives $-1$, the unique nonidentity square root of 1 there. Every exponent giving $-1$ is odd. Now use $3\mid 2^h+1$, $\ell>5$ and $Q\equiv3\bmod8$. The cases $h=1,2$ cannot have $Q>5$. Since $2^{2h}\equiv1\bmod d$ but the order exceeds $h$, it is exactly $2h$.

## Distinct blocks of the binary expansion

Assume $0<c\le d/5$ and put $L=\lceil\log_2 d\rceil$. Then $\lfloor 2^h/d\rfloor=Q-1$. The $h$ digits of $Q-1$, padded on the left, are the first $h$ digits after the binary point of $1/d$. Since $Q$ is odd these are the digits of $Q$ except at the units position. For $0\le s<h-L$, a block of length $L$ starting after the first $s$ digits has integer value
$$\left\lfloor\frac{2^L (2^s\bmod d)}d\right\rfloor.$$
The residues are distinct by the order bound, and $2^L\ge d$ makes their integer block values distinct. All these blocks occur in an expansion of $Q$ satisfying the carry test.

Use the four states for the preceding bit and the carry in $Q+2Q$ in the order $(0,0),(1,0),(1,1),(0,1)$. The allowed successor lists are $(1,2),(3),(4),(1,3)$ with numbering from one. With $\rho=303/200$ the positive vector $(1,\rho-1,\rho(\rho-1),\rho^2(\rho-1))$ has successor sums at most $\rho$ times its entries, as in the [binary counting proof](predecessor_five_seven.md#full-factor-five). Its sum divided by its smallest entry is below eight. There are at most $8\rho^L$ readable blocks from all starting states together. The count also bounds the number of reversed blocks.

Consequently $h\le L+8\rho^L<15d^{3/5}$, including the case $h\le L$. For the last inequality use $L\le2\sqrt d$ and $\rho^5<8$. This estimate applies to any approximation $0<dQ-2^h\le d/5$ in the second hypothesis of the theorem. The [Fourier bound for two predecessor powers](paired_power_orders.md#a-gap-in-the-orbit) gives the independent estimate $h<164\sqrt d$ for these same binary rows. When $c=1$, the order is at most $2h$, and the subgroup argument in that note gives $h\ll_\epsilon d^\epsilon$ for each fixed $\epsilon>0$.

## An all-exponent check up to cofactor $2^{52}$

It remains to check $dQ=2^h+1$ for $3<d<2^{52}$, $d\equiv3\bmod24$. No upper bound on $h$ is imposed. At $h\ge52$ the lowest 52 bits $u$ of $Q$ determine $d$ uniquely through $du\equiv1\bmod2^{52}$. Enumerate all such odd prefixes passing the bit test on $3u,2u$, without imposing any terminal condition on the next bits. There are 879,862,521 prefixes, giving 293,293,896 cofactors in the required residue class after omitting 3.

The first program starts from $du=1+2^{52}v$. At a prefix of length $i$, the next inverse bit is $b=v\bmod2$, and the new residual is $(v+bd)/2$. A conflicting bit excludes every larger exponent. Reaching $v=1$ gives the only possible exponent $h=i$, since $\operatorname{ord}_d(2)=2h$ is necessary. The loop terminates in any event. At $j=\operatorname{ord}_d(2)$ and at every sufficiently large multiple of $j$, the residues of $3d^{-1},2d^{-1}$ modulo $2^j$ (or the larger power with that exponent) force a carry by the same negative-residue calculation. The six recovered pairs with $h\ge52$ have exponents $55,55,57,63,69,75$.

The independent program extends the digits of $d$ first. From $du=1+2^i v$, new bits $x$ of $d$ and $y$ of $u$ obey $y\equiv v+x\bmod2$, with new quotient $(v+xu+yd+2^ixy)/2$. It tests the overlap by direct multiplication. At 52 bits it independently inverts $d$ modulo $2^{128}$. Only 18 cofactors pass through position 127. Wider exact inversion finds a conflicting bit by position 186 for each of these. A separate divisor enumeration for every odd $h<128$, followed by exact tests of the 18 longer prefixes, reproduces the six pairs above and also covers every $h<52$. After the factor-five conditions only $Q=43{,}691$ and $174{,}763$ remain at $d>3$. At their respective $h=51,57$, direct inversion in all $k$ classes modulo $4h$ gives no binary split. A multiple of $4h$ is a period for all tests used here, since $2^{2h}\equiv1\bmod Q$.

## An unrestricted bound on smooth cofactors

For each odd $r\mid d$ prime put $e_r=v_r(2^{r-1}-1)$ and
$$C(d)=\prod_{r\mid d}r^{\min(v_r(d),e_r)}.$$
For a possible exception with $c=1$ and $d>3$, $d\mid 2^{2h}-1$. The elementary lifting formula applied to $(2^{r-1})^{2h}-1$ gives $v_r(d)\le e_r+v_r(h)$. Thus $d\le C(d)h$. Applying the two elementary bounds for $h$ gives
$$d<\min\{(15C(d))^{5/2}, (164 C(d))^2\}.$$
For a fixed set of prime divisors this bounds $d$ effectively, and then $2h=\operatorname{ord}_d(2)$ bounds $Q$.

For odd primes $r\le400000$ the only cases $2^{r-1}\equiv1\bmod r^2$ are $1093,3511$ by exact modular exponentiation. Their orders are $364,1755$, respectively. Neither can divide $d$, since $2^h\equiv-1$ for odd $h$ would force the order at $r$ to be twice an odd integer. Therefore $C(d)=\operatorname{rad}(d)$ if that radical is at most 400,000. The exact inequality $(164\cdot400000)^2<2^{52}$ puts all such candidates inside the checked cofactor range.

Run the following commands from the project directory, with python-flint installed. The Python verifier certifies its integer factorizations, compares both full searches and checks the multiplicative orders in the radical argument.

~~~sh
mkdir -p .work
c++ -std=c++17 -O3 check_factor_five_order.cpp -o .work/check_factor_five_order
c++ -std=c++17 -O3 check_factor_five_order_cross.cpp -o .work/check_factor_five_order_cross
python3 check_factor_five_order.py .work/check_factor_five_order .work/check_factor_five_order_cross
lean CarryOrder.lean
~~~
