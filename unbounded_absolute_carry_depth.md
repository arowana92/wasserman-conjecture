# No absolute carry depth suffices at every weight

For each integer $h\ge1$, there are $\gg 6^{-h}X/(\log X)^2$ weights $N\le X$ with three positive rows whose first common carry, minimized over **all primes**, occurs at digit $h+1$. The implied constant is absolute and the threshold for $X$ may depend on $h$. Every weight counted satisfies Wasserman's conjecture for every choice of rows. Its predecessor has exactly two distinct prime factors.

Digit one is the units digit. At a given prime the depth for three rows is the maximum of their first carry positions, or infinity if one addition is carry-free. The theorem gives an unconditional obstruction to any fixed absolute depth. The 128-digit bound holds on a set of density one. Large powers at 2 and 3 enter the construction below. Another construction proves unbounded depth even when every prime occurs just once in $N(N-1)$. A 490-digit example has sharp depth 513 after maximizing over all triples at its weight. [Density-one theorem](bounded_carry_depth.md) · [Squarefree endpoints](squarefree_carry_depth.md) · [Sharp depth 513](sharp_depth_513.md)

## Choosing the rows

There is an absolute integer $b\ge1$ chosen in the proof. Write
$$
 B=3^h,\quad m=3^b,\quad c=Bm,\quad a=2^{h+1},\qquad
 N=cV=1+aq .
$$
We take $q>a$ prime and $V$ either prime or a product of two large primes. Besides the endpoint $E=(1,1,N-2)$ use
$$
 A=(2^h q,2^h q,1),\qquad
 C=(BkV,BkV,B(m-2k)V)
$$
with $k$ an odd positive integer less than $m/2$. Only primes dividing $K_E=N(N-1)$ need consideration. Row $A$ avoids $q$ completely because its upper counts sum to $a<q$ and its low residues sum to one. At 2 its first carry is at digit $h+1$, since $q$ is odd. Both $E$ and $C$ carry at 2 in the units digit, as $N,V$ and all three parts of $C$ are odd. All parts of $C$ are divisible by $3^h$. It remains to make $C$ avoid **every** prime factor of $V$ completely. These facts then give exact depth $h+1$ for the displayed triple. The [two-factor predecessor theorem and endpoint reduction](all_weights_endpoint_reduction.md) give the ordinary conjecture at the entire weight.

For $V$ prime greater than $c$, take $k=1$. If $V=st$, $s\le t$ primes and $t>cs$, the addition for $C$ is carry-free in base $t$ for every allowed $k$. After division by $t$ the three parts sum to $cs<t$. At $s$, divide out the common factor $s$ and put $z=Bt$. For each $j\ge1$, the residues of $kz,kz,(m-2k)z$ sum to the residue of $mz$ modulo $s^j$ precisely when
$$
                2\{kz/s^j\}\le\{mz/s^j\}.                       \tag{1}
$$
Braces denote fractional parts. Indeed the third residue is the nonnegative representative of the difference. All layers with $s^j>ct$ pass automatically. We will only need one or two coordinates in (1).

## How much of a digit period can fail?

For any $\epsilon>0$, we can choose $b$ independent of $h$ such that for either $J=1$ or $J=2$, outside a set of measure at most $\epsilon$ in the unit $J$-cube, some odd $1\le k<m/2$ satisfies
$$
                     2\{kx_i\}\le\{mx_i\}\quad(1\le i\le J).       \tag{2}
$$
Choose block length $w$ and number of blocks $L$, and set $b=wL$. In the ternary expansions of the $J$ coordinates, inspect the first $L$ disjoint blocks of $w$ digits. If block $i$ (numbered from zero) is zero in every coordinate, use $k=3^{wi}$. Then all $\{kx_j\}<3^{-w}$ and $2k<m$. Also require every $\{mx_j\}\ge 2\cdot3^{-w}$. The measure where either requirement fails is at most
$$
            (1-3^{-wJ})^L+2J\,3^{-w}.
$$
Choosing $w$ first and then $L$ proves the claim. Use terminating expansions on grid boundaries. The exceptional set just described is a finite union of half-open boxes for each fixed $w,L,J$.

For $s\nmid B$ prime, let $t$ run through a complete system modulo $s^J$. The residues $Bt\bmod s^J$ are uniform. The leading fractional digits of the vector $(\{Bt/s^j\})_{1\le j\le J}$ in base $s$ are the $J$ independent digits of that residue, each divided by $s$ to give the grid coordinate. Each coordinate differs from its grid coordinate by less than $1/s$. Consequently the proportion $\rho_s$ in the exceptional boxes is at most $\epsilon+o(1)$ as $s\to\infty$, uniformly for $J\in\{1,2\}$. This also follows by bounding grid cells that meet the finitely many box faces.

## Supplying the primes and preserving the digit restrictions

Apply Friedlander and Iwaniec's form of Chen's theorem to $a q+1=cV$, with $q\le(cY-1)/a$. The integers $a,c,1$ are pairwise coprime, $ac$ is even and the prime divisors of $ac$ are exactly 2 and 3. The lower bound in the theorem gives at least
$$
                    c_0\,Y/(a\log^2 Y)                              \tag{3}
$$
distinct $V\le Y$, for an absolute $c_0>0$ and all sufficiently large $Y$ depending on $a,c$. Their prime factors exceed $((cY-1)/a)^{3/11}>Y^{13/50}$ and there are at most two, counted with multiplicity. The uniform coefficient in (3) matters because $m$ will be chosen after two error tolerances. In the cited theorem the lower bound is an absolute positive constant times a factor depending only on the primes of $ac$, times $x/(c\log^2 x)$ with $x=(cY-1)/a$. [Friedlander and Iwaniec, *Opera de Cribro*, Theorem 25.11, as stated in Garcia–Luca–Shi–Udell, Theorem 6](https://arxiv.org/pdf/1906.05927#page=7)

Discard $V\le Y/(\log Y)^3$ at a cost $o(Y/(a\log^2Y))$ for fixed $a$. For composite $V$ write $V=st$, $s\le t$, so $Y^{13/50}<s\le\sqrt Y$. Fix $0<\delta<1/200$. Call $s$ a boundary prime if $Y^{1-\delta}\le s^i\le Y^{1+\delta}$ for $i=2$ or 3. By Mertens's reciprocal-prime estimate,
$$
 \sum_{s\ \mathrm{boundary}}1/s\ll\delta+o(1),\qquad
 \sum_{Y^{13/50}<s\le\sqrt Y}1/s\ll1 .                              \tag{4}
$$
For other $s$, there is a $J\in\{1,2\}$ with $s^{J+1}<Y^{1-\delta}$ and $s^{J+2}>Y^{1+\delta}$. This uses $s^4>Y^{1.04}$. Thus $ct\le cY/s<s^{J+1}$ and $t\ge Y/(s\log^3 Y)>cs$ for large $Y$. Exactly the first $J$ tests suffice in (1), with any redundant top tests harmless.

Here is the upper sieve estimate, including the dependence on the digit period. For a fixed $s$ and an arbitrary subset of residues modulo $P=s^J$ of proportion $\rho$, count $t\le Y/s$ in that subset with $t$ and $(cst-1)/a$ primes as above. Take $P=1,\rho=1$ for a boundary prime. Integrality specifies a single class of $t$ modulo $a$. Sieve at primes $5\le\ell\le R<s$. For squarefree $d$ supported on these primes, divisibility of $t(cst-1)$ by $d$ specifies $2^{\omega(d)}$ classes modulo $d$. They are independent of the classes modulo $a$ and $P$. In each class modulo $ad$, complete blocks of $P$ successive terms give the exact proportion $\rho$, and the remainder costs $O(P)$. The resulting sequence has main mass $\rho Y/(as)$, local density $2/\ell$ and remainder $O(P2^{\omega(d)})$.

Selberg's upper bound at level $R^2$ therefore bounds the count by
$$
       O\left(\frac{\rho Y}{as(\log R)^2}+P R^{2+o(1)}\right).       \tag{5}
$$
Both primes being counted exceed $R$ in the applications below. For details of the estimate, the sieve coefficients at $d\le R^2$ are bounded by $3^{\omega(d)}$, and $\sum_{d\le R^2}6^{\omega(d)}\le R^{2+o(1)}$. The Selberg denominator sums $\prod_{\ell\mid d}2/(\ell-2)$ for squarefree $d<R$. Using primes from 5 through $R^{1/10}$ gives full Euler product $\asymp(\log R)^2$. Under its normalized product weights the expectation of $\log d$ is $\sum 2\log\ell/\ell=(1/5+o(1))\log R$. Hence at least half the mass has $d<R$. All main-term constants in (5) are absolute. [Selberg upper bound, Theorem 13.1](https://kskedlaya.org/ant/chap-selberg.html)

At boundary primes take $R=Y^{1/20}$. By (4) and (5) they cost at most $O(\delta+o(1))Y/(a\log^2Y)+O(Y^{0.6+o(1)})$. At other primes set $R=Y^{\delta/8}$ and sieve only the exceptional boxes of the digit test. Their main terms cost $O((\epsilon+o(1))\delta^{-2})Y/(a\log^2Y)$. The errors sum to $O(Y^{1-3\delta/4+o(1)})$, since $P=s^J<Y^{1-\delta}/s$. These bounds include every composite case discarded from (3).

Choose $\delta$ small compared with $c_0$, then $\epsilon$ small compared with $c_0\delta^2$, then $w,L$ and $m$. All are absolute. Finally let $Y$ tend to infinity with $h$ fixed. At least half of (3) survives. Each remaining $V$ gives the asserted rows, with $k$ chosen from a successful block when $V$ is composite. Substitution of $Y=X/c$ proves the theorem.

The verifier exercises (1) against digit-by-digit addition, the digit-block count over a complete finite grid, and the CRT counts used in (5). It also checks fully factored examples of the composite case.

~~~sh
python3 check_unbounded_absolute_depth.py
~~~
