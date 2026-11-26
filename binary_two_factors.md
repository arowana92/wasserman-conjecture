# Cofactors that work at every exponent

Let $D=PQ$ with $P,Q>1$ odd prime powers of distinct bases. For all but $O(X/\log X)$ of these $D\le X$, Wasserman's full trinomial conjecture holds at $1+2^eD$ for **every** $e\ge1$. The number of such cofactors is $\gg X\log\log X/\log X$, so the exceptional proportion tends to zero. The same conclusion holds simultaneously for any fixed finite list of prime bases in place of 2.

The list can even grow. Fix $A>0$ and put $B(X)=\exp((\log\log X)^A)$. For all but $O_A(X(1+\log\log\log X)/\log X)$ of these cofactors, the conjecture holds at $1+b^eD$ for **every prime $b\le B(X)$ and every $e\ge1$**.

There is also a finite classification with unrestricted exponents. Of the 16,805 unordered pairs of odd prime powers $P,Q\le1000$ with distinct bases, **16,796 work at every $1+2^ePQ$**. The remaining nine pairs are listed below. They too work for every $1\le e\le10^{24}$.

## Removing the exponent from a carry condition

If $b$ shares a base with $P$ or $Q$, the [two-factor predecessor result](endpoint_reduction.md) and [all-weight reduction](all_weights_endpoint_reduction.md) already prove the conjecture. Otherwise $N-1=b^e P Q$ has three distinct prime bases. Its [mandatory CRT row](endpoint_factors.md) avoids all three bases, and its part exceptional at $P$ is $kb^eQ$ with $1\le k<P$. Combining the other two parts gives $(P-k)b^eQ+1$. At $b$ the low $e$ positions hold just the digit one. Above these positions, the addition of **$kQ$ and $(P-k)Q$ must have no carry in base $b$**, independently of $e$.

For $b=2$, write $\beta=\log_2((1+\sqrt5)/2)$. The [binary affine bound](binary_affine_counts.md) gives $O_P(Y^\beta)$ possible $Q\le Y$ for any fixed $P$, after summing over $k$. In particular the proportion among prime partners also tends to zero for each fixed odd $P$.

For odd $b$, the prefix argument in the same note gives a two-summand bound as follows. Include the complemented digits of the total $PQ$. The minimum delay is zero since $b\nmid P$. If all three expressions have delay zero, at most $(2b+1)/3$ next input digits can pass, by summing the output digits over the permitted choices. Otherwise two have delay zero. Give each prefix of $Q$ weight $z^{-h}$ where $h$ is the sum of pending digits. If the current column contains assigned digit sum $S$, at most $b-S$ choices remain, and the one delayed expression permutes the new digits. The sum of relative weights is at most
$$z^S\sum_{j=0}^{b-1-S}z^{-j}.$$
For $b\ge7$ take $z=\exp(1/b)$. This bound is at most $\exp(1)/(4(1-\exp(-1/b)))\le\exp(1)(b+1)/4\le(3b+1)/4$, using $x(1-x)\le1/4$ for $x=z^{S-b}$. For $b=3,5$ take $z=4/3,6/5$ respectively. Checking $0\le S<b$ gives the same bound $(3b+1)/4$. The case of three undelayed expressions satisfies it as well.

Put $\rho_2=(1+\sqrt5)/2$ and $\rho_b=(3b+1)/4$ for odd primes. For either parity, the number of permissible length-$L$ prefixes at a fixed $P,k$ is at most $P^2\rho_b^L$. Indeed, the largest delay is at most $\log_b P$, so the total pending digit sum is at most $(b-1)\log_b P$. Each odd-base choice above satisfies $z^{b-1}<\exp(1)$, giving weight loss at most $P^{1/\log b}$. In base two the exact bound gives loss at most $P^\beta$. These bounds include all positive $Q<b^L$ and are uniform in $k$.

## Counting the cofactors

Let $B\ge2$ be the largest base under consideration, $c=\log(6/5)$ and $H=X^{c/(8\log B)}$. For $P\le H$ apply the prefix bound with $L=1+\lfloor\log_b(X/P)\rfloor$. Since $\rho_b/b\le5/6$, one choice of $b,P,k$ contributes at most $bP^2X^{1-c/\log B}$. Sum over $b\le B$, $P\le H$, $k<P$. This gives $O(B^2 H^4 X^{1-c/\log B})=O(B^2 X^{1-c/(2\log B)})$, which is $o(X/\log X)$ for fixed $B$ or the stated growing $B(X)$.

It remains to count products with both full powers greater than $H$. In a dyadic interval $[t,2t)$ there are $O(\sqrt t\log(2t))$ higher prime powers and each supplies at most $X/t$ partners. Summing gives $O(X\log(2H)/\sqrt H)$ products having a factor of exponent at least two. This contribution is again $o(X/\log X)$. For two primes, sum $\pi(X/p)$ over $H<p\le\sqrt X$. The upper prime-counting bound and partial summation give
$$\ll\frac{X}{\log X}\left(1+\log\frac{\log X}{\log H}\right)\ll\frac{X}{\log X}(1+\log\log B).$$
For the lower bound on the total family, take odd primes $p\le X^{1/4}$ and $p<q\le X/p$. The estimate $\pi(t)\asymp t/\log t$ gives $\gg X/(p\log X)$ partners per $p$. Partial summation of the same lower estimate gives $\sum_{p\le X^{1/4}}1/p\gg\log\log X$. Each product is counted once. [Rosser and Schoenfeld, Theorem 1](https://doi.org/10.1215/ijm/1255631807)

## The exhaustive binary calculation

Set $T=2^e$, $d=PQ$. The CRT row has parts
$$Tx,\quad Ty,\quad Tz+1,\qquad x=kQ,\quad y=lP,\quad z=d-x-y.$$
The bits of $x,y,z$ partition those of $d$. All three are positive, since $z=0$ would give a part equal to one on a row whose parts have a common factor of at least three. Also $Tx\equiv1\pmod P$ and $Ty\equiv1\pmod Q$. Some odd prime $u\mid N$ must be assigned to this row, with its full power in $N$ dividing all three parts. Hence $u\mid\gcd(k,l)$, $dT\equiv-1\pmod u$, and $v_u(dT+1)\le v_u(\gcd(k,l))$.

The code first tries the possible $k<P,l<Q$ and uses complete tables of powers of 2 modulo $P,Q,u$. Intersect the three exponent classes by the generalized Chinese remainder theorem. There are 1,096 binary rows with $\gcd(k,l)\ge3$ passing the two predecessor orbits, and 288 choices of row and scale prime passing all three orbits. Each gives a single class $e=e_0\bmod m$.

For each $r$ among $u$ and the bases of $P,Q$, set $F=r^{v_r(2^m-1)}$. All inputs $T$ are fixed modulo $F$. In higher positions allow arbitrary base-$r$ digits. At each modulus $w=r^j$ require
$$[xT]_w+[yT]_w+[zT+1]_w<w,$$
where brackets denote the least nonnegative residue. At $r=u$ also reject if $w\mid dT+1$ but $w\nmid\gcd(k,l)$. Emptying the set of residues proves exclusion at all exponents. The code can retain a shape when its per-level work limit is reached. Of eleven shapes left by these tests, one more has no possible terminal input in the scale-prime carry graph. For that test, fix the low digits modulo $F$, keep the three multiplication carries, and explore all reachable states. A finite input can terminate only if the three remaining carries add without overflow. This leaves ten shapes at nine pairs.

| $P$ | $Q$ | $u$ | $e_0\bmod m$ |
|---:|---:|---:|---|
| 53 | 227 | 3 | $907\bmod5876$ |
| 167 | 337 | 7 | $801,717\bmod1743$ |
| 191 | 769 | 3 | $3764\bmod36480$ |
| 223 | 577 | 3 | $619\bmod5328$ |
| 239 | 257 | 11 | $7030\bmod9520$ |
| 379 | 797 | 47 | $1149704\bmod3460212$ |
| 439 | 653 | 7 | $68989\bmod142788$ |
| 509 | 947 | 13 | $311842\bmod720852$ |
| 599 | 929 | 7 | $384366\bmod416208$ |

For the finite exponent check write $e=e_0+mj$. At increasing moduli $u^h$ the period in $j$ initially stays one, then grows by $u$ on each lift since $u$ is odd. Traverse all lifts that pass the preceding residue and valuation conditions. Once the period exceeds the maximum index, each branch contains at most one allowed exponent. Test the two predecessor bases at one layer beyond their prescribed full powers and then all three bases through 64 digits. A positive carry at any tested layer excludes that exponent.

The pair $(379,797)$ admits a second test. Here $\gcd(k,l)=47$, so $U=47$ and the last row is $V(f,g,h)$ with $f+g+h=47$ and $V=N/47$. Its exponent progression fixes $V\equiv51\pmod{81}$. This row must avoid 3. Only six sorted triples pass the first four ternary digits, and exploring their finite carry graphs leaves exactly $(8,15,24)$. The scalar 2481 verifies that this last ternary language is nonempty. Lifting the exponent on these ternary digits provides a separate exclusion at this pair through $10^{24}$.

The default calculation first imposes a predecessor layer if it admits at most two exponent classes. Using the ternary companion at $(379,797)$, the ten shapes produce 19,201,131 individual exponent candidates and 2,536,337 pass both displayed predecessor layers. All are rejected. Using the first row alone throughout gives 183,022,185 and 105,516,439 respectively, again with no surviving exponent. The second verifier generates the 288 initial classes by binary submasks, rechecks all 277 local exclusions by lifting exponents, and independently checks the 229 reachable states in the terminal exclusion. Its direct enumeration through $10^8$ verifies the individual scale-prime tests, including 41 positive controls for the joint low layers.

```sh
python3 check_binary_two_factors.py
python3 check_binary_two_factors_cross.py
python3 check_binary_two_factors.py --scale-only
```
