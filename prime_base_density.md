# Prime bases that work at every exponent and cofactor through 1,000

Consider $N=d p^a+1$ with $p$ prime, $a\ge1$, and $1\le d\le1000$.

- **All exponents simultaneously.** There is an explicit set of prime bases of relative natural density greater than $97.23\%$ for which the full conjecture holds at every one of these weights.
- **For every fixed prime $p$.** At most $O_p(K^{0.336})$ exponents $a\le K$ can admit a counterexample at any cofactor in this range. If $d=874$ is omitted, the union over the other 999 cofactors has at most $O_p(K^{1/3})$ exceptional exponents.
- **For individual small primes.** With $p\le10000$, the full conjecture holds at every exponent except possibly on the 13 pairs below. On those pairs it also holds through $a=10^{24}$. In particular every prime $p<1733$ works for all exponents and all 1,000 cofactors.

| $d$ | Prime bases still possible beyond the exponent bound |
|---:|---|
| 380 | 1733, 1973, 5717 |
| 629 | 4943 |
| 780 | 6379 |
| 850 | 2851, 3271, 8731, 8971 |
| 880 | 3889 |
| 950 | 4463 |
| 962 | 8861 |
| 986 | 2447 |

## From trinomial rows to residue classes of the base

By the [shifted-denominator theorem](shifted_directions.md), a counterexample must have $p\nmid d$ and rows
$$A=(a_1W,a_2W,a_3W+1),\qquad C=V(f,g,h),\qquad W=p^a,\quad UV=dW+1.$$
The theorem specifies all 21 combinations of normalized rows and scale parity. At even scale, $p=2$ and $d$ must be odd. None of the remaining odd-denominator slopes adds without carries in base 2, whereas stripping $2^a$ from the row avoiding 2 requires this. The 16 odd-scale cases each come with a congruence $W\equiv R\pmod M$ and a set of primes forced on $A$. The residues $R$ in these cases are units.

For each residue $z$ of a prime base modulo $M$, enumerate its successive powers. Retain it if some power is $R$. Then the exponent must satisfy $a\equiv k\pmod D$ where $D$ is the order of $z$ and $z^k=R$. Finitely many primes dividing the fixed moduli can be discarded when computing relative prime density.

Local restrictions on $W$ give restrictions on the prime base itself. For every prime $\ell\mid2dU$ take $w=\ell^{v_\ell(M)}$, possibly 1, and $q=w\ell^2$. Test residues $t\equiv R\pmod w$ that are units at $\ell$. The residue of $V$ is determined modulo $q/\gcd(U,q)$ from $UV=dt+1$. Require that $V$ be a unit at primes of $U$. For each row eligible to take $\ell$, read the low digits of the corresponding scalar ($W$ or $V$) in its carry graph. Retain a row only if each digit is permitted and the resulting carry state can still reach an accepting state. This last condition is checked by backward graph reachability. Require each forced row, or at least one eligible row when there is a choice. This gives a set $F_\ell$ of allowed values modulo $q$.

There are further constraints even at primes outside $2dU$. If
$$\sum_{x\in(f,g,h)} (jx\bmod\ell)\ge\ell\qquad (j=1,\ldots,\ell-1),$$
then $C$ carries for every positive $V$. Indeed the first nonzero digit of $V$ already causes a carry in $C$. Such a prime outside $2dU$ cannot divide $dW+1$, as it would then divide $V$ and be assigned to $C$. Add the corresponding allowed unit residues modulo $q=\ell$, with $w=1$. Only $\ell\le U$ can pass this test, by taking $j=1$.

## Counting prime residues

Fix $z,k,D$ as above. For the local prime $\ell$, a possible base $y\bmod q$, with $y\equiv z\pmod w$, must satisfy
$$ y^k\langle y^D\rangle\ \cap\ F_\ell\ne\varnothing\quad\text{in }(\mathbb Z/q\mathbb Z)^\times.$$
Count the fraction of its unit lifts passing this test. Conditional on $p\equiv z\pmod M$, the choices at the distinct $\ell$ are independent coordinates in the finite CRT product. Multiply their fractions, sum over $z$, and divide by $\varphi(M)$. Each local test may choose a different exponent, giving an upper bound. The prime number theorem in arithmetic progressions identifies each unit residue class in the common finite modulus with its relative natural density among primes. Summing over the 16 shapes bounds the excluded fraction by
$$\frac{14701541568851133613999}{530916553722651093597600}<0.027691.$$

The first verifier computes coset intersections by prime-power lifting. For a unit $s\equiv1\pmod\ell$ (modulo 4 at $\ell=2$), its powers modulo $q$ are exactly the units congruent to 1 modulo $\gcd(s-1,q)$. This follows from the valuation formula for $s^b-1$. The short initial cycle modulo $\ell$, or modulo 4, takes care of other units. The independent verifier constructs discrete logarithm tables at odd prime powers and directly enumerates powers at 2. It checks prefix carries by factorial floors and tests completion by searching from the residual multiplication carries. The two calculations return the same rational number.

## Counting exponents for any fixed prime

For a forced base $\ell$, let $B$ be the carry graph's nonnegative integer adjacency matrix for $A$ with input $W$. Count multiple digit edges separately. The checker forms the positive integer vector $v=(I+B)^{160}\mathbf1$. For each odd-scale case it finds a forced $\ell$ and checks componentwise
$$ (Bv)_i^b\le\ell^s v_i^b,$$
with $s/b=1/3$ except at $d=874$, where $s/b=42/125$. Since $\mathbf1\le v$, the number of admissible input blocks of length $h$ from the initial state is at most $v_0\ell^{hs/b}$. This counts low blocks even when the remaining higher digits are unknown.

Here $p\ne\ell$ since $p\nmid dU$ for a candidate, and the prime base $p$ is odd. For fixed $p,\ell$ the order of $p$ modulo $\ell^h$ is $\gg_{p,\ell}\ell^h$. Take a positive power of $p$ equal to 1 modulo $\ell$ (modulo 4 at 2), and apply the valuation formula as in the [prime-power ray proof](prime_power_ray_density.md#counting-exponents). Choose the least $h$ giving order greater than $K$. Then $\ell^h\ll_{p,\ell} K$ and the residues of $p^a$ for $1\le a\le K$ are distinct. Every surviving residue has an admissible length-$h$ block. Summing the resulting bounds over finitely many shapes gives both estimates.

## Checking all small prime bases

For $p\le10000$, first check the initial orbit modulo $M$ and the carry-free requirement on the slopes of $A$ at $p$ itself. There are 333 base/shape pairs at this point. For each prime of $2dU$, use successively the highest prime power at most 256, 4096 and 65536, increased when necessary to include its initial modulus and the unitarity digit. Include the preceding fixed-divisor conditions. At each modulus directly enumerate the exponents in a full period whose powers have admissible residues. Test joint feasibility of these exponent classes using the [generalized CRT check](successor_exponents.md#combining-prime-power-orbits-across-bases). Exactly the 13 listed pairs survive.

For a surviving pair, choose one forced prime $\ell$ dividing $M$. Starting at $a\equiv k\pmod D$, lift its powers one $\ell$-adic digit at a time. At a new layer there is one exponent class or $\ell$ classes above the preceding class, according to whether the order stays fixed or increases by $\ell$. Discard a branch only when the row carries at that layer or its least candidate exceeds the bound. When the stride exceeds $10^{24}$, test the remaining individual exponents at all forced bases through 60 digits. The traversal reaches 21,136,355 individual exponents counted across the pairs, and the modular tests reject every one. The [lifting routine and its independent small controls](check_successor_exponents_cross.py) are shared with the successor-weight calculation.

Run the following from the repository root with Python and python-flint.

```sh
python3 check_direction_prime_bases.py --controls
python3 check_direction_basis_cross.py
python3 check_direction_small_bases.py
```
