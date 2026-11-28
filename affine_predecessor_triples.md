# Three predecessor powers on an affine line

Fix three integral linear forms $w_i(t)=a_i t+b_i$ with $a_i>0$ and distinct rational roots. At each positive integer $t$ for which their values are powers greater than one of three distinct primes, set $M=w_1w_2w_3$ and $N=M+1$. The bases and exponents can all vary.

**Theorem.** In any arithmetic progression of $H\ge2$ positive parameters with odd common difference, at most
$$
 O_{w_1,w_2,w_3}\big(H^\alpha\log(2H)\big),\qquad
 \alpha=\log_2\lambda=0.551463\ldots,\qquad \lambda^3=\lambda^2+1,\quad\lambda>1,
$$
can give a failure of the full trinomial conjecture under the three prime-power hypotheses. The constant is effective and independent of the start and step. For each fixed line there is also a finite, computable list of possible companion row proportions, valid at every eligible parameter. The theorem bounds the absolute number of exceptions among the parameters. It requires no lower bound for the frequency of prime-power tuples.

## A bound for the common scale

By the [endpoint reduction](all_weights_endpoint_reduction.md) and its [three-factor consequence](endpoint_factors.md#when-n-1-has-three-distinct-prime-factors), a failure forces a row with parts $A_i=k_i M/w_i$, where $k_i$ is the unique integer in $[1,w_i-1]$ satisfying $k_iM/w_i\equiv1\pmod{w_i}$. This row avoids all three predecessor primes. There is a unitary factor $U\ge3$ of $N$ dividing every $k_i$. The companion has parts $(N/U)(f_1,f_2,f_3)$, with positive integers $f_i$ of sum $U$.

Put $A=\max a_i$, $B=\max_{i,j}|a_i b_j-a_j b_i|$. Whenever all $w_i>a_i^2$, even the gcd of the three CRT multipliers satisfies
$$
                 \gcd(k_1,k_2,k_3)\le B^5+2A^3B^2.                 \tag{1}
$$
This bound uses only the pairwise coprimality of the values. It applies whether or not the CRT parts sum to $N$.

Indeed write $D_i=\prod_{j\ne i}(a_i b_j-a_j b_i)$. Reduction modulo $w_i$ gives integers $m_i$ with
$$
             D_i k_i=a_i^2+m_i w_i.
$$
For $D_i>0$, $0\le m_i<D_i$. For $D_i<0$, $D_i\le m_i<0$. These follow from the bounds on $k_i$ and $a_i^2<w_i$. Order the roots $r_i=-b_i/a_i$ increasingly. Then $D_1>0>D_2$ since $a_i b_j-a_j b_i=a_i a_j(r_i-r_j)$. If $m_1=0$, the gcd divides $a_1^2$. Otherwise the affine numerators $a_i^2+m_i w_i(t)$ for $i=1,2$ have distinct roots. The first root is strictly to the left of $r_1$, and the second strictly to the right of $r_2$. The common gcd divides their nonzero determinant
$$
 m_1m_2(a_1 b_2-a_2 b_1)+m_1a_1a_2^2-m_2a_2a_1^2 .
$$
Its absolute value gives (1), since $|m_i|\le B^2$. The finitely many smaller positive parameters can be checked directly.

There are now finitely many choices of $(m_1,m_2,m_3),U,(f_1,f_2,f_3)$. On each candidate the forced parts are rational polynomials of degree two or three, and the companion parts are constant rational multiples of $N$. Keep only choices for which the forced parts sum identically to $N$, apart from the finitely many zeros of any nonzero difference. On residue classes modulo a common multiple of the denominators of the coefficients, substituting $t=r+Lj$ gives integral coefficients on every class where the values are integral. Discard constant terms that are nonintegral. Every part is nonconstant. The forced part has the two distinct real roots of the other two forms. Its possible third linear factor has a real root. Its derivative therefore has only simple roots, by Rolle's theorem with multiplicities. The total has derivative $M'$, with two distinct real roots by the same argument. The companion derivatives are nonzero constant multiples of $M'$. These properties persist under reindexing.

## A binary count for polynomial rows

Let $F_1,F_2,F_3$ be fixed integral polynomials, $F_0=F_1+F_2+F_3$. Suppose each of the four derivatives is nonzero and has only simple complex roots. Constant derivatives are allowed. Among residues modulo $2^n$, at most $O_F((n+1)\lambda^n)$ extend to an integer where the parts have nonnegative values and add without any binary carry.

At depth $h$, stop splitting a residue class if all four derivative valuations $e_i=v_2(F_i'(x))$ are strictly less than $h$ at its representative $x$. They are then constant on the class. At every depth only $O_F(1)$ classes still need splitting. In fact for each derivative $G$, the number of roots modulo $2^h$ is bounded independently of $h$. For nonconstant $G$, an integral polynomial Bézout combination of $G,G'$ equals a nonzero integer $c$. Write $b=v_2(c)$. At any root modulo $2^h$, $h>b$, one has $v_2(G'(x))\le b$. Within each class modulo $2^{b+1}$, two such roots differ by a multiple of $2^{h-b}$, by expanding $G(y)-G(x)$ at $x$. This bounds the number of roots, including the finitely many small depths. Constant derivatives satisfy the same bounded-count assertion. Consequently each depth produces only $O_F(1)$ newly stopped classes.

On a stopped class $x+2^h z$, two distinct inputs $z,z'$ have
$$
 v_2\big(F_i(x+2^h z)-F_i(x+2^h z')\big)=h+e_i+v_2(z-z').
$$
The linear term has strictly lower valuation than the higher Taylor terms. Thus the next input bit toggles the next undetermined bit of each output, at delay $e_i-\min_j e_j$ after a common shift. This is precisely the bit property used in the [affine prefix bound](binary_affine_counts.md#counting-low-digit-strings). For completeness, complement the total's bits. Carry freedom allows at most one of the four bits in each column to be one. The number of undelayed outputs is even by the derivative sum. If all four are undelayed, at most one next choice can pass. Otherwise two are undelayed. Record future columns with a determined one, reject collisions, and give a prefix with $b$ recorded columns weight $\lambda^{-b}$. An occupied current column allows at most one choice and raises weight by at most a factor $\lambda$. At a free column the two choices contribute at most $1+\lambda^{-2}=\lambda$ times the current weight, since across them the two delayed bits contribute two ones in total. The fixed initial bits can be recorded before branching. With $D=\max e_i-\min e_i$, at most $D$ columns are recorded, so at most $\lambda^{n-h+D}\le\lambda^n$ prefixes of length $n-h$ can pass. The test can use bits beyond column $n$ when their values are already determined by the input prefix. Summing over stopped classes through depth $n$ and the bounded number still open there proves the estimate.

The endpoint coefficient $N(N-1)$ is even. Hence a failure requires at least one of the two polynomial rows to add without a binary carry. On restricting the given progression to a parameter class $t=r+Lj$, the induced step in $j$ is still odd. For $n=\lceil\log_2 H\rceil$ the indices in that class have distinct residues modulo $2^n$. Applying the estimate to the finite list of rows proves the theorem.

## A pointwise family

There is a simpler test that can certify every parameter. Let $W,Z,Y$ be full powers of three distinct primes, let $a,b,c,d$ be positive integers with $|aZ-bW|=1$, and suppose
$$
 Y\equiv cW+dZ\pmod{WZ},\qquad d\mid a^2,\quad c\mid b^2,\qquad a^2/d<W,\quad b^2/c<Z .
$$
Then the full conjecture holds at $1+WZY$. The CRT multipliers at $W,Z$ are respectively $a^2/d,b^2/c$. Multiplying these numbers by $ZY,WY$ gives 1 modulo the respective powers by squaring the determinant relation. The two multipliers are coprime since $a,b$ are. This contradicts the required $U\ge3$.

For example take $W=5t+3, Z=8t+5, t\ge5$. At every eligible pair of powers, **every** third prime power $Y\equiv18t+11\pmod{WZ}$ works. Here $5Z-8W=1$, $18t+11=2W+Z$, and the two forced multipliers are 25 and 32. At $t=1{,}001{,}542$ the three numbers $5{,}007{,}713, 8{,}012{,}341, 18{,}027{,}767$ are prime. This gives the weight $723{,}337{,}184{,}691{,}130{,}355{,}012$ and infinitely many further weights by replacing the last prime by its power of exponent $1+10{,}030{,}872{,}791{,}520j$, $j\ge0$. The period is the lcm of the first two primes minus one each.

Run `python3 check_affine_predecessor_triples.py` to reproduce the checks. Direct modular inversion verifies (1) on 32,698 coprime triples. Of these, 642 have both the required row sum and a multiplier gcd of at least three. A separate enumeration tests 9,314 binary input prefixes in 201 stopped classes by comparing the determined bits directly. The verifier also checks 3,984 members of the explicit residue family, proves the three example primes by trial division and verifies the period by modular exponentiation.
