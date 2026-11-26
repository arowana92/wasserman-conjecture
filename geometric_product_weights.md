# A geometric sum times an independently shifted pair

**Theorem.** For every prime $p$ and all positive integers $a,i,k$, Wasserman's three-row conjecture holds at
$$N=1+p^a(1+p^i+p^{2i})(1+p^k).$$
The three exponents vary independently. This includes all coincident binary positions. In terms of the six labelled positive digit positions, start with a three-term arithmetic progression and take its translate by any positive distance. The theorem applies to their union with multiplicity at every outer shift.

Write $I=p^i$, $Z=p^k$ and $m=(1+I+I^2)(1+Z)$. The six terms in $m$ are indexed by the rectangle $\{0,1,2\}\times\{0,1\}$, with point $(u,v)$ representing $I^uZ^v$. The proof first treats a possible counterexample in which the row avoiding $p$ pairs its units-place 1 with one of those terms. The two binary alternatives are proved separately below.

## The shifted pair

Assume the seven labelled units add freely. The [endpoint reduction](all_weights_endpoint_reduction.md) supplies a nonendpoint row avoiding $p$. By the [seven-unit theorem for the units-place digit](seven_large_bases.md#locating-the-units-place-digit), at odd $p$ the 1 must belong to a pair. Whenever it does, the row is
$$(p^a x,\ p^a y,\ 1+p^a z),$$
where $x,y,z$ are disjoint sums of two terms, three terms and one term of the rectangle. The [prime-power assignments](seven_span_and_gaps.md#which-allocations-need-checking) give positive divisors $U,A,C$, all coprime to $p$, satisfying
$$U\mid x,y,\qquad A\mid x,y+z,\qquad C\mid y,x+z,\qquad 27m^2\le U^3 A^2 C^2. \tag{1}$$
Also $U\ge3$ and $U$ is coprime to $m$. If $x$ or $y$, after removing a coordinate monomial, divides $m$ in $\mathbb Z[I,Z]$, then $U\mid m$, which is impossible. Exact polynomial division excludes 18 of the 60 labelled allocations.

For each remaining allocation the verifier contradicts (1) using the three pairs of polynomials displayed there. The outer exponent $a$ has disappeared. The calculation below gives a finite proof for the whole positive quadrant in $(i,k)$.

## Resultants on a sector

Begin with the two coordinate directions and subdivide a sector between neighboring integer rays $v,w$ at their sum. Neighboring determinants always have absolute value one. Using both orders of the rays, every integer point in the sector has a representation
$$(i,k)=(b+\ell)v+\ell w,\qquad b,\ell\in\mathbb Z_{\ge0}.$$
In a given order use the substitution
$$I=T^{v_1}Y^{w_1},\qquad Z=T^{v_2}Y^{w_2},\qquad T=p^{b+\ell},\quad Y=p^\ell.$$
For each pair in (1), discard the common monomial factor from the terms of each individual polynomial. Its evaluated factors are all units modulo the assigned divisor. Take the resultant in $T$, or use the polynomial independent of $T$ if there is one. Remove its lowest monomial power of $Y$. The result is a nonzero integer polynomial $D(Y)$ in every sector used. The assigned divisor divides $D(p^\ell)$ by the integral resultant identity. Denote the three bounds by $D_0,D_1,D_2$ and put $(t_0,t_1,t_2)=(3,2,2)$.

Let $h(i,k)=4i+2k$ so $m^2>p^{h(i,k)}$. Keep a sector when both orientations satisfy $h(v+w)>1+\sum_j t_j\deg D_j$. Otherwise subdivide it. The verifier constructs the complete subdivision, which has 1,092 oriented final sectors across the 42 allocations. In each final sector with all three evaluations nonzero, (1) would require
$$27p^{b h(v)+\ell(h(v+w)-\sum_j t_j\deg D_j)}
<\prod_j\|D_j\|_1^{t_j}.$$
Both coefficients in the exponent are positive. Taking $p\ge2$ gives finite loops over $b,\ell$ and an integer upper bound on $p$. For $\ell=0$ use $|D_j(1)|$ in the prime ceiling. At small ceilings the code computes the actual integer divisors: $U$ is bounded by the part of $\gcd(x,y)$ coprime to $p$, and $AC$ divides $\gcd(m,\operatorname{lcm}(x,y))$ since its prime powers are pairwise coprime.

## Vanishing projections and the univariate checks

The union of all positive integer root lists is $\{1,2\}$. A root at 1 gives an entire ray $(i,k)=b v$. A root at 2 forces $p=2,\ell=1$ and gives an affine line $(i,k)=(b+1)v+w$.

For each interior ray from a root at 1, and for finite cases with prime ceiling above 64, specialize to its primitive direction $(n,m_1)$ using $I=X^n,Z=X^{m_1}$. An actual $X$ is a power of $p$. Strip powers of $X$ from each block in (1). The [polynomial gcd and integer resultant procedure](polynomial_prime_powers.md#an-effective-test-including-reducible-polynomials) supplies polynomial bounds $u(X),r(X),s(X)$ for the three divisors. The verifier checks that
$$729\big((1+I+I^2)(1+Z)\big)^4-u^6 r^4 s^4$$
has positive constant coefficient and nonnegative remaining coefficients after a translation $X=X_0+V$. Across 320 shape and ray pairs the largest translation is 8. The 148 smaller integer evaluations and 2,206 finite prime and row evaluations all contradict (1).

On a line from a root at 2, use $I=2^{w_1}X^{v_1}, Z=2^{w_2}X^{v_2}$ with $X=2^{b+1}$. In each of the 12 incidences, the two bounds for $A,C$ at $Y=2$ are nonzero. Use these constants and obtain $u(X)$ for $U$ by specializing $x,y$ and removing powers of $X$ and the even part of the bound's content. The same polynomial inequality has positive constant and nonnegative remaining coefficients after translating by 2 in every case. This finishes the shifted-pair proof.

## Binary alignment when the six positions are distinct

Take $p=2$. The six relative positions are symmetric about $H/2$, where $H=2i+k$. If 1 belongs to a triple in the row avoiding 2, [binary alignment](seven_binary_alignment.md) requires a pairing of the six positions with lengths $n,n,n\pm1$. Its size inequality gives $27m^2\le(2^n+1)^3$ and thus $H<3n/2$. Here $H\ge5$, so even $n-1>H/2$ (check $H=5$ directly). Every pair therefore goes from the bottom three positions to the top three.

If $k>2i$, the top three are a translate of the bottom three by $k$. The sum of the lengths is $3k$, contrary to $3n\pm1$. Otherwise $0<k<2i$, with $k\ne i$ by distinctness. The bottom three form the set $D=\{0,i,k\}$ and the top three are $H-D$. Every pairing has lengths $H-u-\pi(u)$, $u\in D$, for a permutation $\pi$ of $D$. The identity and both three-cycles give three distinct lengths. A transposition of $u,v$ gives repeated length $n=H-u-v$ and remaining length $H-2w$, where $w$ is the fixed position.

- Interchanging $0,i$ gives $n=i+k$ and $|i-2k|=1$. The inequality $H<3n/2$ requires $k>i$, which is impossible.
- Interchanging $0,k$ gives $n=2i$ and $|2i-k|=1$. The same inequality requires $k<i$, which is impossible.
- Interchanging $i,k$ gives $n=i$, whereas $H>2i$.

This rules out the binary alignment at every distinct profile.

## Coincident positions

At odd primes every position occurs at most twice, so the seven labelled units do add freely. In binary, coincidences require $k=i$ or $k=2i$. At $k=i$, the inner polynomial in $t=2^i$ is $1+2t+2t^2+t^3$, which gives at most five units in $N$ including the unit at 1. The [five-unit theorem](six_digit_gaps.md) applies. At $k=2i$ the inner polynomial is $1+t+2t^2+t^3+t^4$. For $i=1$ it again gives at most five units in $N$.

For $i\ge2$ the five relative positions in the latter case are $0,i,2i+1,3i,4i$. The [six-unit equal-gap theorem](equal_predecessor_gaps.md) requires two disjoint unshifted pairs of the same gap $d$, and $4i<3d/2$. The exceptional residue $2i+1$ modulo $i$ cannot lie in either pair, since the other gap is a multiple of $i$. The four remaining positions allow the equal gaps $i$ and $3i$, and the span inequality requires $d=3i$. The unshifted parts are then $2^a F$ and $2^a t F$ with $F=t^3+1$. The equal-gap theorem also gives $U\mid F$ and $27m^2\le U^3$.

Every prime $q\mid t+1$ divides the first row's coefficient. In fact after removing the full power of $q$ in $F$, the two unshifted parts have opposite nonzero residues modulo $q$, so their addition carries. Hence $U$ is coprime to $t+1$ and
$$U\le F/(t+1)=t^2-t+1<t^2.$$
But $m>t^4$, contradicting $27m^2\le U^3$. This handles the final coincidence and proves the theorem.

## Reproduction

With SymPy and python-flint installed, run the following from the repository root. The second program evaluates all three projected polynomials on the final sectors at two integers and compares them to independently formed Sylvester determinants or constant polynomials, accounting for stripped monomials. It also checks every binary pairing at 8,064 distinct profiles and the 128 coincident profiles in its control range.

~~~sh
python3 check_geometric_product.py
python3 check_geometric_product_cross.py
~~~
