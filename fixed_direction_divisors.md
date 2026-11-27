# Fixed prime divisors of every dilation of a multinomial row

For a positive vector $x$ of length $r\ge2$ and sum $s$, write $K_x(g)=(gs)!/\prod_i(gx_i)!$ and $F(x)=\gcd_{g\ge1}K_x(g)$. The exact integer $F(x)$ is computable by a shortest-path calculation on at most $2s-r$ states for each prime $p\le s$. It gives one fixed certificate across infinitely many weights. The construction also exhibits rows for which a prime always divides the coefficient even though its first carry can be arbitrarily late.

### A family at every prime

Let $p$ be prime and $a>b\ge1$. At **every** positive scale $t$, the coefficient with parts
$$
                   t(a,pb,(p-1)(a-b))
$$
is divisible by $p$. If $p^j-1$ divides $t(p-1)(a-b)$, this coefficient is divisible by $p^j$. There is no requirement that $p$ divide the weight.

Write $s_p$ for digit sum in base $p$. The parts sum to $t(pa+b)$. Legendre's digit-sum identity and $s_p(pz)=s_p(z)$ give the exact equality
$$
 v_p(K)=v_p\binom{(pa+b)t}{pat}+\frac{s_p(t(p-1)(a-b))}{p-1}.
$$
The last fraction is a positive integer. Every positive multiple of $p^j-1$ has digit sum at least $j(p-1)$. Indeed, fold the digits into $j$ bins by their exponents modulo $j$, then replace $p$ tokens in any bin with one token in the next bin, cyclically. This preserves the represented residue modulo $p^j-1$ and strictly decreases the positive number of tokens. The finished digits represent a positive multiple of $p^j-1$ at most $p^j-1$, so every bin has $p-1$ tokens.

For an odd prime take $a=2$, $b=1$. The direction is $(2,p,p-1)$ of sum $2p+1$. None of the sums of two entries divided by their gcd, or of all three primitive entries, is a multiple of $p$. For **every** $k\ge1$, the scale $t=(p^k+1)/2$ has exactly one carry, at layer $p^{k+1}$. To check this directly, the sum of the least residues at $p$ is $(p+1)/2$. At $w=p^i$ for $2\le i\le k$, the residues are $1,(w+p)/2,(p-1)/2$. Their sum is less than $w$. At the next layer each part is less than $p^{k+1}$ and their sum belongs to $[p^{k+1},2p^{k+1})$. All higher layers exceed the weight.

This gives a single common prime for two whole families of rows and the endpoint. At every $N=gp(2p+1)$, take the row in direction $(2,p,p-1)$ and **any** row of primitive sum $p$. Both coefficients and $N$ are divisible by $p$. The second coefficient contains its primitive sum by factorial decrement. In the first, the multiplier is $gp$, so at $g=(p^k+1)/2$ the sole $p$-carry occurs at layer $p^{k+2}$, even though $v_p(N)=1$ throughout.

### The exact fixed divisor

For each $p\le s$ use the slopes $(x_1,\ldots,x_r,s)$. After reading $j$ low digits of a multiplier with numerical value $h<p^j$, take as state
$$
         c_z=\lfloor zh/p^j\rfloor\qquad(z=x_1,\ldots,x_r,s).
$$
A new digit $d$ sends the state to $c'_z=\lfloor(c_z+zd)/p\rfloor$. Give this edge the nonnegative cost $c'_s-\sum_i c'_{x_i}$. This is the factorial-floor contribution at layer $p^{j+1}$ for *every* completion of the multiplier, since the unknown higher digits cancel between the weight and the parts. Once all digits have been read, append zeros until the zero state is reached. The total cost is exactly the prime valuation. Conversely any finite path starting and ending at zero represents such a multiplier. Removing initial zero digits preserves the valuation, so require the first digit to be nonzero.

The states are constant on the half-open cells made by the breakpoints $i/z$ in $[0,1)$ of the functions $\lfloor zy\rfloor$ for the slopes. Including zero, there are at most $1+(s-1)+\sum_i(x_i-1)=2s-r$ cells. After the first digit a shortest path to zero can be chosen simple, since all costs are nonnegative. Thus $v_p(F(x))$ is attained at some $1\le g<p^{2s-r}$. Primes above $s$ are already avoided at scale one. The [verifier](check_fixed_direction_divisors.py) computes shortest paths and witness multipliers, and independently obtains the minima by solving the reverse distance equations on the rational cells. This extends the [carry-free transition test](carry_intersections.md#exhaustive-search-by-intervals) to the full valuation.

The **set of prime divisors** of $F(x)$ is invariant under replacing $x$ by $m x$. To see this, suppose $K_x(t)$ avoids $p$. Remove the $p$-part from $m$ to obtain $m'$. Choose $L$ large enough to fit all digits of $st$ and so that $p^L\equiv1\pmod{m'}$, with the congruence vacuous for $m'=1$. Concatenate $m'$ copies of the carry-free block using multiplier $t(1+p^L+\cdots+p^{(m'-1)L})$. It is divisible by $m'$ and still gives a carry-free row. Appending enough low zero digits makes the multiplier divisible by $m$. The reverse inclusion follows by restricting the scales. The valuations themselves can grow after dilation, as the digit-sum bound shows.

For two base rows $A,C$ of the same weight $n$, the greatest fixed divisor of **both coefficients and the endpoint at every scale** is exactly
$$
               \gcd\big(F(A),F(C),n\gcd(2,n-1)\big).
$$
Indeed the gcd over $g$ of $gn(gn-1)$ is $n\gcd(2,n-1)$, as is seen by taking $g=1,2$ and observing $g(gn-1)$ is even for odd $n$. The formula distinguishes a fixed prime working everywhere from the scale-dependent prime in the [primitive-direction theorem](primitive_row_scales.md).

### What the first occupied digit decides

For a primitive positive triple $x$ and a prime $p$, the first occupied digit overflows for **every unit multiplier modulo $p$** if and only if $p\mid \sum x_i$ or some pair of nonzero residues sums to zero. If the total vanishes, primitivity forces the residue sum to be a positive multiple of $p$. Opposite nonzero residues already sum to $p$. When there are only two nonzero residues and their sum is nonzero, one of the multipliers $1,-1$ gives sum less than $p$. With a single nonzero residue the same conclusion holds directly.

For the converse with three nonzero residues and nonzero total, consider the standard tetrahedron with vertices $0,e_1,e_2,e_3$ in the lattice $\mathbb Z^3+\mathbb Z x/p$. In each nontrivial lattice coset the least nonnegative representative has all coordinates strictly positive, and their sum cannot equal one. Thus universal overflow at the first digit means that the tetrahedron is empty. White's theorem gives it lattice width one. A nonzero integral linear functional on this lattice then has values, after a possible sign change, in $\{0,1\}$ on the four vertices. Its coefficient vector picks out a nonempty subset of coordinates whose sum is zero modulo $p$. The hypotheses exclude singletons and the full triple, leaving an opposite pair. [Bruce Reznick, *Clean Lattice Tetrahedra*, discussion of White’s theorem](https://arxiv.org/abs/math/0606227)

Higher digits can force further primes even when every part and every nonempty subset sum is a unit. For example $5\mid F(2,7,12)$. On the three multiplication carries, the only base-five transitions that keep the addition carry-free from zero are
$$
 (0,0,0)\xrightarrow{3}(1,4,7)\xrightarrow{2}(1,3,6)
                         \xrightarrow{2}(1,3,6)
$$
and the initial zero loop. These follow by trying the five digits at each displayed state. A positive finite multiplier cannot follow this path back to zero. The infinite digits $3,2,2,\ldots$ pass every finite prefix. Both this example and the odd-prime family above show why a bounded search for an initial carry can miss a prime that divides every coefficient on a direction.

Run `python3 check_fixed_direction_divisors.py`. The verifier compares 600 shortest paths with the rational-cell computation, checks 1,050 exact multinomials, and tests the digit identities and delayed carry layers at 579 inputs. It also verifies the first-digit classification on 65,796 normalized residue triples at primes through 97.
