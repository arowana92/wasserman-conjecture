# Polynomial weights at powers of any prime

Let $G\in\mathbb Z[X]$ have nonnegative coefficients, degree $n\ge1$, and positive leading coefficient. Let $P=p^a$ with $p$ prime and $a\ge1$, and write $R=(1+\sum_i G_i^2)^{1/2}$, where $G_i$ are the coefficients. Irreducibility below is over $\mathbb Q$.

**Theorem.** Suppose $G$ is primitive and irreducible. If $XG-1$ is irreducible, Wasserman's three-row conjecture holds at $N=PG(P)$ for all $P>R^{(12n+3)/(2n-1)}$. If $XG+1$ is irreducible, the conjecture holds at $N=PG(P)+1$ for all $P>R^{(14n+3)/(2n)}$. In particular $P>R^{15}$ suffices in both cases, uniformly in the prime, the exponent and the polynomial degree. The threshold depends on the coefficients even if their occupied powers are far apart.

The theorem also gives a lower bound for the common divisor of the endpoint coefficient and any two coefficients $K_A,K_C$ when $p\nmid K_A$. Put $D=\gcd(N,K_A,K_C)$ and $J=\gcd(N-1,K_A,K_C)$. Whenever $P>\max_i G_i$, under the respective irreducibility hypotheses,
$$
\begin{aligned}
 PG(P):\quad &D^3J^2\ge \frac{27(N-1)^2}{P^3 R^{12n+3}},\\
 PG(P)+1:\quad &D^3J^2\ge \frac{27G(P)^2}{R^{14n+3}}.
\end{aligned}
$$
For fixed $G$, the common divisor $DJ$ therefore grows at least as fast as a positive constant times $P^{(2n-1)/3}$ or $P^{2n/3}$ respectively.

## Proof by specialization

Use the [retained prime-power assignment](endpoint_pair_bound.md). After removing the factors contributing to $D,J$, assign the remaining powers to $A,C$ as $DUV=N$ and $JST=N-1$, with $U,S$ assigned to $A$. The normalized parts of $C$ have sum $DU$, and at every retained power of $T$ two of them are divisible by that power. Hence
$$27(N-1)^2\le D^3J^2 U^3S^2. \tag{1}$$
If $p$ is avoided by both rows assign it to $A$. Its full power in $N(N-1)$ is retained.

Write $N=PG(P)+\epsilon$ with $\epsilon=0,1$. Since $P>\max_i G_i$, the coefficients of $G$ occupy disjoint base-$p$ blocks at $P$. By Kummer, the parts of $A$ have the form $P L_i(P)+\delta_i$, where the nonnegative integer polynomials $L_i$ sum coefficientwise to $G$, $\delta_i\in\{0,1\}$, and $\sum\delta_i=\epsilon$. When $\epsilon=1$ label the shifted slot 3. Every unshifted $L_i$ is nonzero. Every nonzero $L_i$ is proper because at least two slots are positive before the shift. No proper nonzero $L_i$ is a rational multiple of $G$, by primitivity and coefficientwise containment.

For integer polynomials $A,B$, every common divisor of $A(P),B(P)$ divides $\operatorname{Res}(A,B)$, by the integral Bézout identity for the resultant. Its absolute value is at most $\|A\|_2^{\deg B}\|B\|_2^{\deg A}$, by Hadamard's inequality on the Sylvester determinant. These statements include nonzero constant polynomials.

At the unshifted weight write $U=Pv$. Then $v\mid G(P)$ and $v\mid L_i(P)$ for every $i$, so $v\le R^{2n}$. Each power in $S$ divides two $L_i(P)$ and $PG(P)-1$. Thus $S^2$ divides the product of the three nonzero resultants $|\operatorname{Res}(XG-1,L_i)|$, yielding $S^2\le R^{6n+3}$. Substitution in (1) proves the first inequality. Since $N-1\ge P^{n+1}/2$, it excludes $D=J=1$ above the claimed threshold.

At the shifted weight, $U\mid L_1(P),L_2(P)$ and $U\mid PG(P)+1$, giving $U\le R^{2n+1}$. Put $S=PB$. Each retained prime power in $B$ divides $G(P)$ and at least one of $L_1(P),L_2(P)$. For a base distinct from $p$ cancel $P$ in an unshifted zero slot. At $p$ the full power from $PG(P)$ divides both unshifted parts, so the cancellation is valid here too. The two nonzero resultants with $G$ give $B\le R^{4n}$. Equation (1) now proves the second inequality. Finally, the [all-weight endpoint reduction](all_weights_endpoint_reduction.md) makes $D=J=1$ necessary for a counterexample, and some nonendpoint row must avoid $p$ at each of these two weights.

## An effective test including reducible polynomials

The computations here allow arbitrary coefficient splits, including ones that would carry in a particular base. In a counterexample $U\ge3$, so a shifted part of 1 is impossible and every $L_i$ is nonzero. For each shape we bound the assignments by evaluations of integer polynomials. This extends the [binary polynomial degree criterion](polynomial_binary_families.md) and the [quadratic specialization bounds](quadratic_powers.md).

A list of nonzero integer polynomials has an explicit polynomial $B$ such that every common evaluation divisor divides $B(P)$. First take the gcd of nonzero constants in the list and of the pairwise resultants of nonconstant members. A nonzero answer already bounds the gcd. Otherwise divide the polynomials by their common polynomial gcd and repeat on the quotients. If their pairwise resultants all vanish, take the resultant of the first quotient with an integer linear combination of the others. One can use coefficients $1,t,t^2,\ldots$ in increasing $t$ until this is nonzero. At each root of the first quotient the polynomial in $t$ is nonzero, by collective coprimality. Multiply the resulting constant by the original polynomial gcd.

For the unshifted shape this construction bounds $v$ using $(G,L_1,L_2,L_3)$. For each pair $i,j$, use $(XG-1,L_i,L_j)$, then take the integer-polynomial lcm of the three bounds for $S$. Each assigned full power divides one bound, and their bases are pairwise coprime. For a shifted shape use $(XG+1,L_1,L_2)$ for $U$. At exceptional slot $k$, the bound for its group of $B$ uses $G$ and the three congruences $P L_i(P)+\delta_i-[i=k]\equiv0$. Where $\delta_i=[i=k]$ use $L_i$ itself after cancelling $P$. The prime $p$ has exceptional slot 3, where its entire power before cancellation is available. Take the lcm over $k$ as before.

A further refinement uses full powers. Suppose a unitary divisor of $F(P)$ divides $L(P)$, with $F$ monic and $L$ a constant times a polynomial divisor of $F$. If an irreducible factor $f$ has multiplicity $k>0$ in $L$ but $e>k$ in $F$, then for every prime $r$ in the assigned divisor,
$$(e-k)v_r(f(P))\le v_r((L/f^k)(P)).$$
The assigned part of $f(P)$ therefore divides a nonzero constant bound for $\gcd(f(P),(L/f^k)(P))$. Replace $f^k$ by that constant to power $k$ in the divisor bound. This applies to $v,S,U,B$ in their respective target polynomials. It can remove a whole polynomial degree.

Use these bounds in (1) with $D=J=1$. For instance, if $v\mid f(P)$ and $S\mid g(P)$ then a counterexample requires
$$729(PG(P)-1)^4-P^6 f(P)^6 g(P)^4\le0.$$
The shifted expression is $729G(P)^4-f(P)^6g(P)^4$ for bounds on $U,B$. A polynomial with positive leading coefficient is positive on an explicitly verified half-line. The verifier doubles an integer shift until all Taylor coefficients there are nonnegative and its constant term positive. Strict degree bounds often give a smaller immediate check. At remaining prime powers above the [initial interval](full_range_1b.md), test the numerical gcds and the actual carry conditions. A symbolic bound with an inconclusive leading term leaves the shape unresolved.

## Certified coefficient boxes

For monic $G$ of degree $n$, let each of its $n$ lower coefficients independently range from 0 through $H$. The counts give polynomials for which the verifier proves the conjecture at **every** prime power $P$. The two middle columns may include different polynomials.

| $n$ | $H$ | Polynomials | $PG(P)$ | $PG(P)+1$ | Both |
|---:|---:|---:|---:|---:|---:|
| 2 | 10 | 121 | 106 | 121 | 106 |
| 3 | 3 | 64 | 58 | 64 | 58 |
| 3 | 5 | 216 | 198 | 213 | 195 |
| 4 | 3 | 256 | 243 | 247 | 234 |
| 5 | 2 | 243 | 234 | 236 | 227 |
| 6 | 1 | 64 | 61 | 62 | 59 |
| 7 | 1 | 128 | 123 | 124 | 119 |
| 8 | 1 | 256 | 250 | 248 | 242 |

In particular, **every** $P^3+bP^2+cP+1$ with $0\le b,c\le10$, and **every** $P^4+bP^3+cP^2+dP+1$ with $0\le b,c,d\le3$, satisfies Wasserman's conjecture for every prime power $P$. In all the displayed finite computations, the necessary size inequality with the actual numerical gcds already ruled out every point reaching beyond the initial interval.

Run the following with python-flint and SymPy installed. The optional flag lists precisely the polynomials and sides still unresolved in the specified box.

```sh
python3 check_polynomial_powers.py 3 5 --list-exceptions
python3 check_polynomial_powers_cross.py
```

The cross-check compares 784 Sylvester determinants with SymPy resultants, tests partial-factor removal on 60 specializations, and takes the maximal first-row assignments from 965 coefficients computed by exact binomial products. These include 99 tests involving nonconstant polynomial bounds.
