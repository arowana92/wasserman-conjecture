# A modular inverse that must split into two carry-free parts

Let $p$ be a prime, $j\ge1$, and $t\ge1$ coprime to $p$ with at most two distinct prime divisors. Consider Wasserman's three-row conjecture at $N=1+p^j t$.

* It holds for **every** $j$ if $t=(p^m-1)/d$ and $d<p$.
* The same holds whenever $d\mid p^\ell-1$ for a proper positive divisor $\ell$ of $m$, or $m=2s$ and $d\mid p^s+1$.
* It holds for every $j,m,p$ whenever $1\le d\le2500$ and $d\mid p^m-1$. For $p=2$ this computed cutoff is 30,000, and for $p=3,5,7$ it is 12,000.
* For example every $1+179(2^{89}-1)2^j$ satisfies the conjecture. Its binary expansion has 90 ones. More generally take any Mersenne prime $2^s-1$ and any full or partial odd prime power $b\mid 2^s+1$. The result at every $1+b(2^s-1)2^j$ follows by taking $m=2s$, $d=(2^s+1)/b$. No assumption on the number of Mersenne primes is needed.

## The inverse and its digits

If $t$ has at most one prime divisor, the [endpoint reduction and two-factor theorem](all_weights_endpoint_reduction.md) apply. Otherwise write $t=P Q$ for full powers of distinct primes. By the [three-factor assignment](endpoint_factors.md#when-n-1-has-three-distinct-prime-factors), one nonendpoint row $A$ of any counterexample must avoid all three predecessor primes. There is also a product $U\ge3$ of full prime powers of $N$ that divides every part of $A$ and whose prime bases $A$ avoids. Its complementary row has three positive parts summing to $U$ after normalization.

Write $W=p^j$. In the order of exceptional positions for $P,Q,W$, the actual parts of $A$ are $W x,W y,W z+1$. Here $x,y>0$, $z\ge1$ (a part equal to one would contradict $U\ge3$), and $x+y+z=t$. Since $Q\mid x$, $P\mid y$, and $Wx\equiv1\bmod P$, $Wy\equiv1\bmod Q$, the integer $r=x+y$ is the **least positive residue of $W^{-1}$ modulo $t$**. The split itself is unique:
$$
 x=Q[\,rQ^{-1}\,]_P,\qquad y=P[\,rP^{-1}\,]_Q,
$$
where $[\,\,]_v$ takes the least nonnegative residue modulo $v$. Also $U\mid\gcd(x,y)$ since $W$ is coprime to $N$.

Base-$p$ avoidance says that the digits of $x,y,z$ add without a carry. In particular the digits of $r$ are bounded by those of $t$ at every position. Its digit sum is at least four: each of $x,y$ is a positive multiple of a prime different from $p$ and cannot be a pure power of $p$. In particular $r$ cannot be a pure power of $p$. These tests apply directly to any exponent and cofactor, even without the special forms above.

There is also a pairing across different prime-power multipliers. For fixed $P,Q$ let $W,W'$ be prime powers coprime to $PQ$. If $W'\equiv-W\bmod PQ$, at least one of $1+Wt,1+W't$ satisfies the full conjecture. Indeed the corresponding positive CRT parts $x',y'$ must equal $t-x,t-y$, giving $x'+y'>t$ when $x+y<t$. More precisely, in exactly half the reduced classes of $W\bmod t$ the two CRT representatives already sum to more than $t$, proving the full conjecture throughout those classes. Their sum is never exactly $t$, and negation pairs the classes without a fixed point. For a fixed base $p$ this pairs opposite exponents whenever $p^h\equiv-1\bmod t$.

## Divisors of a power minus one

Suppose $dt=p^m-1$. Choose $0\le h<m$ congruent to $-j\bmod m$. Write $p^h=ut+r$, $0<r<t$. If $u=0$, the inverse is a pure power. This always settles $d<p$ since $t\ge p^{m-1}$ there (the trivial $t=1$ was handled separately).

More generally, if $u>0$ and some $B=p^\ell\le t$ satisfies $(u+1)(t\bmod B)<B$, the digits cannot fit. In fact $h\ge\ell$, and with $a=t\bmod B>0$ the residue of $r$ modulo $B$ is $B-u a>a$. This gives a test at individual exponents as well.

If $d\mid p^\ell-1$ and $\ell\mid m$ with $\ell<m$, put $B=p^\ell$, $a=(B-1)/d$. Then $t=a(1+B+\cdots+B^{m/\ell-1})$, so $t\bmod B=a$ and $B\le t$. Since $p^h<p^m=dt+1$, we have $u<d$ for $t>1$ (equality $p^h=dt=p^m-1$ is impossible for $p^m>2$). Thus $(u+1)a\le da=B-1$. This proves the repeated-block case.

For the other case put $H=p^s$, $b=(H+1)/d$, so $t=b(H-1)$ and $m=2s$. The cases $t=1$, $b=1$ or $d=1$ follow directly, reducing exponents modulo $s$ when $b=1$. Take $u>0$. Then $h\ge s$, and $p^h=H L$ for a power $L\le H/p$. We have $ub\le H L/(H-1)<L+1$, hence $ub\le L$. Since $b>1$ is coprime to $p$, equality is impossible. In particular $b\le ub<L<H$. The two base-$H$ blocks of the inverse and the cofactor are
$$ r=(L-ub)H+ub,\qquad t=(b-1)H+(H-b). $$
Digit containment gives $L-ub\le b-1$ in the upper block. In the lower block the digits of $H-b$ complement those of $b-1$ to $p-1$. Hence $ub$ and $b-1$ add without carries. Both are below the power $L$, so their sum is below $L$. The resulting inequalities $L\le(u+1)b-1<L$ contradict each other. Lean 4.36.0 checks the modular wrap inequality and the contradiction between the two block bounds. [Formal arithmetic lemmas](InverseResidueBlock.lean)

In particular, if the least order modulo an odd prime power $d$ is even, its halfway power is $-1$ modulo $d$ and this case applies to the least order. All larger multiples were already covered by the repeated-block argument.

## Spacing between exceptional exponents

Write $t_i$ for the base-$p$ digits of $t$ (zero beyond the last position). If both $j$ and $j+h$ are exceptional, $h\ge1$, set $s_0=s_p(t\bmod p^h)$ and
$$ D_h=\sum_{i\ge0}\min(t_i,t_{i+h}). $$
Then $s_0=1$ requires $D_h\ge3$, and $s_0=2$ requires $D_h\ge1$. For example, if the units digit is 1, all other digits are 0 or 1, and occupied positions have pairwise distance at least $k$, then at most one of every $k$ consecutive exponents can fail. This bound has no dependence on the number of nonzero digits.

To prove it, match the normalized parts $x,y,z$ by their exceptional bins in the two rows at the different weights. The CRT expressions give $x_{j+h}\equiv p^{-h}x_j\bmod t$, and likewise for $y$ and $z$. All components are in $(0,t)$. Thus every component divisible by $p^h$ at $j$ becomes its exact quotient by $p^h$ at $j+h$. Those shifted components jointly occupy at most $D_h$ digit units because both triples add without carries to $t$. At most $s_0$ of the three old components can have nonzero residues modulo $p^h$. Of the three, $x,y$ have at least two digit units each and $z$ has at least one. This proves the two inequalities.

## Verifying bounded cofactors at all exponents

For each tested $d\ge2$ it suffices to take primes $p<d$ coprime to $d$ and $m=\operatorname{ord}_d(p)$. All other positive multiples of this order have repeated blocks. Set $t=(p^m-1)/d$ and test only $h<m$ with $p^h>t$. Discard residues failing digit containment or digit sum four, cofactors without exactly two prime factors, and incorrect CRT splits. Keep only splits with $x,y,z$ adding without a carry at $p$ and $\gcd(x,y)\ge3$.

For each surviving gcd $g$ test each odd prime $v\mid g$ and also $v=4$ if $4\mid g$. At least one such $v$ divides $U$. Require $j\equiv-h\bmod m$ and $t p^j+1\equiv0\bmod v$. In each exponent progression $j\equiv j_0\bmod M$, test each prime $q$ under $P,Q,v$. At modulus $w=q^k$, $W=p^j$ is fixed if $p^M\equiv1\bmod w$. Otherwise allow all lifts of the surviving residues from the previous level. Discard a residue when the three actual part residues sum to more than $N\bmod w$. When $q\mid v$, discard also if $w\mid N$ but $w\nmid g$, since the full $q$-power of $N$ must divide all three parts. An empty level excludes the entire progression. The code checks that every progression is excluded.

The verifier recomposes factorizations and checks primality of the factors with FLINT. The scan covers 527,235 base and quotient pairs. Of 6,109 digit candidates, 46 survive the split tests and 14 exponent progressions require a further carry check. Its independent small cases construct the original three CRT parts directly at the weight and check valuations by factorial floors. Run from the repository root with Python and python-flint.

~~~sh
python3 check_inverse_residue_cofactors.py
~~~
