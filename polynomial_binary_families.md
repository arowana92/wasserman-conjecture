# Polynomial families evaluated at powers of two

Let $P(X)$ have nonnegative integer coefficients and positive degree. For all sufficiently large $k$, the binary digits of the coefficients occupy disjoint blocks in $P(2^k)$. In particular the binary digit sum is bounded. The [bounded-digit theorem](bounded_binary_weight.md) reduces the full conjecture there to the endpoint $E=(1,1,N-2)$.

**Polynomial criterion.** Split the bits of the coefficients of $P$ into three nonempty groups, giving $P=F_1+F_2+F_3$. For instance, a bit $2^j$ in the coefficient of $X^i$ contributes the term $2^jX^i$ to a group. Only splits with all three $F_i$ nonconstant require testing. Over $\mathbb Q[X]$ put
$$G_0=\gcd(F_1,F_2,F_3),\qquad
G_1=\operatorname{lcm}_{i<j}\gcd(P-1,F_i,F_j).$$
If for every tested split
$$3\deg G_0+2\deg G_1<2\deg P, \tag{1}$$
then the full three-row conjecture holds for all sufficiently large $P(2^k)$.

**Proof.** In an endpoint counterexample the binary-free row is one of the indicated splits. Each of its parts grows without bound by the [endpoint part-size inequality](endpoint_factors.md), so splits with a constant part cannot occur infinitely often. Write $UV=N$, $ST=N-1$ for the prime powers assigned to the binary row and its companion, respectively, with $U,S$ assigned to the binary row. We have $T^2\le(U/3)^3$; also $U$ divides all three actual parts and each full power in $S$ divides two of them.

Choose primitive integral representatives for all polynomial gcds and lcms, with positive leading coefficients. For fixed integer polynomials with gcd $G$, their common numerical gcd at an integer $x$ divides $C G(x)$ for a fixed positive integer $C$: divide out $G$, apply the polynomial Bézout identity over $\mathbb Q$, and clear denominators. Thus $U\ll |G_0(2^k)|$. Every full power in $S$ divides the numerical gcd of $P(2^k)-1,F_i(2^k),F_j(2^k)$ for some pair. Each pair's polynomial gcd divides $G_1$ in $\mathbb Z[X]$ by Gauss's lemma. A common multiple of the three Bézout constants therefore gives $S\ll |G_1(2^k)|$, using the pairwise coprimality of the assigned full powers. It follows that
$$27(P(2^k)-1)^2\le S^2 U^3\ll |G_1(2^k)|^2|G_0(2^k)|^3.$$
Degree comparison gives the contradiction in (1).

## Two consequences

1. Fix integers $d\ge1,c\ge2$ with $\gcd(d,c)\le5$ and $\gcd(d,c-1)\le5$. Then the full conjecture holds for all sufficiently large $d2^k+c$. For $P=dX+c$ any three nonconstant parts have distinct positive slopes, because disjoint nonempty subsets of the bits of $d$ have distinct sums. If $G_0$ is nonconstant, all three parts are rational multiples of $P$. If $G_1$ is nonconstant, two are rational multiples of $P-1$ and the third is 1 plus a multiple of it, by summing. In either case write the slopes as $g z_i$ with $\gcd(z_1,z_2,z_3)=1$ and $D=\sum z_i$. The distinct positive $z_i$ give $D\ge6$ and $D\mid d$. The constant terms, after subtracting the single 1 in the second case, are $z_i c/D$ or $z_i(c-1)/D$. Their integrality forces $D\mid c$ or $D\mid c-1$, respectively. The hypotheses exclude both possibilities, giving (1).
2. Fix nonnegative integers $a,b$ for which $a^2-4(b-1)$ is not an integer square. The full conjecture holds for all sufficiently large $2^{2k}+a2^k+b$. Here $P-1$ is an irreducible quadratic. The leading coefficient 1 means at most one part is quadratic. Thus $G_1$ is constant and $\deg G_0\le1$, giving $3\deg G_0+2\deg G_1\le3<4$. For example this covers the even weights $2^{2k}+31\cdot2^k+256$ and the odd weights $2^{2k}+31\cdot2^k+255$.

The symbolic checker takes the coefficients of $P$ in ascending order and uses SymPy's exact polynomial arithmetic. It reports whether (1) proves the family. A split failing the degree test gives an inconclusive result for that split.

```sh
python3 check_polynomial_family.py 256 31 1
```
