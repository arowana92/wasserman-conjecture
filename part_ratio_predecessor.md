# A sharp predecessor divisor for a fixed ratio of two parts

Here is an infinite collection of part ratios. Include every integer ratio $1:t$. Also take any prime power $P=p^h$ and any $m\ge1$, and include
$$
                  (1+P+\cdots+P^{m-1}):P^m .                         \tag{1}
$$
If each of two positive rows of the same weight $N$ contains a pair from this collection, **their coefficients share a prime dividing $N(N-1)$**. The two ratios, primes, exponents, scales, and positions can be chosen independently. Examples include $7:8$ against any integer ratio, and $5:16$ against $13:27$. In view of the [endpoint reduction](all_weights_endpoint_reduction.md), the two nonendpoint rows in a counterexample cannot both have such pairs.

## The exact constant

Fix coprime positive $a,b$. Write $K=N!/((ag)!(bg)!z!)$ for $g,z>0$ and $N=(a+b)g+z$. Define $C(a,b)$ to be the least positive integer with
$$
                            N-1\mid g\,C(a,b)\,K
$$
for every $g,z$. It divides $ab$. At $p\mid ab$, writing $k=v_p(ab)$, its exact exponent is
$$
 v_p(C(a,b))=\max\left(0,\ k-
       \min_{g\ge1,\ p\nmid g} v_p\binom{(a+b)g-2}{ag-1}\right).     \tag{2}
$$
Indeed factorial decrement gives
$$
 \frac{gK}{N-1}=\frac{N}{abg}
         \binom{N-2}{(a+b)g-2}\binom{(a+b)g-2}{ag-1}.
$$
At fixed $g,p$ the two factors $N$ and $\binom{N-2}{(a+b)g-2}$ are simultaneously units for $N=p^j+1>(a+b)g$, since the upper binomial argument has all digits $p-1$. This proves that the missing exponent is precisely the positive part of $v_p(abg)-v_p\binom{(a+b)g-2}{ag-1}$ at that scale. For $p\nmid ab$, adding $ag-1,bg-1$ carries at least $v_p(g)$ times in the low digits. For $p\mid ab$ and $p\mid g$, the identity
$$
 \binom{(a+b)g}{ag}=\frac{(a+b)((a+b)g-1)}{abg}
                      \binom{(a+b)g-2}{ag-1}
$$
has a unit numerator in the prefactor, so this case contributes no missing exponent either. Formula (2) follows.

The [algorithm](check_part_ratio_divisors.py) for the minimum reads the digits of a unit $g$ from low to high. The multiplication carries for $ag-1,bg-1$ start at $-1,-1$. At digit $d$ divide $x+ad$ and $y+bd$ by $p$. The quotients are the next multiplication carries and the remainders are the two digits to be added. Charge the outgoing addition carry at each step. After the first, nonzero digit, the multiplication carries are in $[0,a-1]$ and $[0,b-1]$. The initial addition carry is zero and the final state has all three carries zero. Each path gives a positive unit multiplier and every multiplier gives a path upon appending enough high zeros. A shortest path computes the exact value in (2).

For rows with pairs $(ag,bg),(ch,dh)$, let $G$ be the gcd of their coefficients with $N(N-1)$. Then
$$
 N(N-1)\mid G\,\operatorname{lcm}\big(gC(a,b),hC(c,d)\big)
          \mid G g h\,\operatorname{lcm}\big(C(a,b),C(c,d)\big).       \tag{3}
$$
At $p^e\parallel N-1$, use the defining property of $C$ on the row of smaller valuation. At $p^e\parallel N$ the retained power divides all parts of the row of smaller valuation, hence $g$ or $h$ by coprimality of the terms of the ratio. Equivalently, $N$ divides each part times the coefficient by factorial decrement. This proves (3). Positivity gives $(a+b)g,(c+d)h\le N-1$. In particular $G>1$ whenever the lcm of the two constants is at most $(a+b)(c+d)$.

## Shifting prime powers between the parts

Take arbitrary $h,k\ge1$, $D=p^{\gcd(h,k)}-1$, $a=(p^k-1)/D$ and $b=p^k B$ with $B=(p^h-1)/D$. Then
$$
     a+b=p^h a+B,\qquad
     v_p\big(C(a,b)\big)=\min(h,k).                                   \tag{4}
$$
The entries $a,b$ are coprime. The rows $(ag,p^kBg,z)$ and $(p^hag,Bg,z)$ have the same weight and the same $p$-valuation because each has the same sum of digit sums. Suppose $p^e\parallel N-1$, and remove the $p$-valuation of the coefficient, leaving $p^j$ with $j>0$. The [retained-power lemma](endpoint_pair_bound.md) gives residues $0,0,1$ in both rows modulo $p^j$. If $p\mid g$, the first two entries in each are nonunits and therefore vanish, giving $p^j\mid g$. If $p\nmid g$, the units $ag,Bg$ occupy the exceptional positions in their respective rows. Both shifted terms must vanish, so $j\le h,k$.

For sharpness let $l=\operatorname{lcm}(h,k)$ and set
$$
 g=\frac{(p^l-1)D}{(p^h-1)(p^k-1)} .
$$
This is an integer since the lcm of the two factors $p^h-1,p^k-1$ divides $p^l-1$. The digits of $ag$ are ones at $0,h,\ldots,l-h$, and those of $p^kBg$ are ones at $k,2k,\ldots,l$. These sets are disjoint. Their sum has low digit one and its next occupied digit is at $\min(h,k)$ (when $k\mid h$, only the second list contributes that digit). Taking $z$ to be a sufficiently high power of $p$ keeps the coefficient a unit and leaves $v_p(N-1)=\min(h,k)$. This proves (4).

When $k=mh$, $B=1$ and $a=1+P+\cdots+P^{m-1}$. At the other primes use $C\mid ab$. Thus $C(a,b)\mid Pa=a+b-1$. For an integer ratio $C(1,t)\mid t$. Substitution into (3) proves the stated common-prime result.

There is also an exact description of the carry-free multipliers in (1). The two entries $ag,P^mg$ add without a carry in base $p$ if and only if **all** $g,Pg,\ldots,P^mg$ add without a carry. To see the forward direction, the equal sum $g+P(ag)$ must also be carry-free by equality of digit sums. Let $F_n(x)$ be the polynomial of base-$p$ digits of $n$ with the units digit as constant coefficient. The two carry-free additions give
$$
       (x^h-1)F_{ag}(x)=(x^{mh}-1)F_g(x).
$$
Cancel $x^h-1$ in $\mathbb Z[x]$. The digit coefficients of the first $m$ shifted copies therefore coincide with those of $ag$, and adding the last copy has no carry by hypothesis. The reverse direction follows by grouping the first $m$ copies. For $p=2,h=1$ this says exactly that consecutive ones in $g$ have positions differing by at least $m+1$.

## Checks

Run `python3 check_part_ratio_divisors.py --height 512` and `python3 check_part_ratio_divisors_cross.py`. Among all 79,852 reduced ratios $1\le a\le b\le512$, exactly 35 have $C(a,b)<ab$. In every case the full reduction is accounted for by (4). There are 27 noninteger ratios with $C(a,b)\le a+b$, all from (1). This classification is finite; (2) decides any specified ratio.

The [cross-check](check_part_ratio_divisors_cross.py) matches the upper bounds against 285,934 individual prime witnesses. Each is a positive trinomial at an explicitly constructed weight, evaluated by factorial floors. It also checks the universal divisor on 5,742 exact multinomials.
