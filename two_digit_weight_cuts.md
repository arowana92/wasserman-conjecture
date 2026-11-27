# Two low digit units determine a ratio of parts

Write $N=1+DW$, where $D,W\ge1$ and $W$ has at most two distinct prime factors. For a prime $p\mid N$, put $P=p^{v_p(N)}$. Suppose that, for powers $r=p^i$ ($i\ge0$) and $B=p^m$,
$$
 N/P\equiv 1+r\pmod B,\qquad B>D(1+r).
 \tag{1}
$$
The two digit units below $B$ have sizes 1 and $r$. The case $r=1$ is possible only for odd $p$.

**Two-prime theorem.** If (1) holds at two distinct weight primes with different values of $r$, the full three-row conjecture holds at $N$. The other factors of $N$ may be arbitrary. In particular, for every $D$ coprime to six, the two conditions
$$
 N\equiv6\pmod{2B_2},\quad B_2>3D,\qquad
 N\equiv12\pmod{3B_3},\quad B_3>4D
$$
suffice, where $B_2,B_3$ are powers of 2 and 3. For each prime $q\nmid 6D$, these conditions prescribe a unit class for primes $\ell$ modulo $2B_2\cdot3B_3$ at $N=1+Dq\ell$. Hence there are infinitely many such weights for each fixed $D,q$, by Dirichlet's theorem. [Dirichlet's theorem](https://ocw.mit.edu/courses/18-785-number-theory-i-fall-2021/mit18_785f21_lec18.pdf)

**One-prime theorem.** Put $s=1+r$. Under (1), the full conjecture holds if
$$
 D\le s(4s+1),\qquad P=2\ \text{or}\ Ps-1\nmid D .
 \tag{2}
$$
In fact the size cutoff is $s(4sG+1)$ with $G$ defined below. For example, at every $N\equiv 2+2^{h+1}\bmod 2^{m+1}$, $h\ge1$, this certifies all representations above with $2^m>D(1+2^h)$ and $D\le(1+2^h)(4(1+2^h)+1)$. This uses the zero digits following the second unit. [Earlier criterion using the first unit](predecessor_weight_gap.md)

## Why the two cases suffice

Use the all-weight endpoint reduction. In a counterexample, assign each prime of $N(N-1)$ to one of the two nonendpoint rows avoiding it. Their products of assigned full weight powers are $U,V$ with $UV=N$, $\gcd(U,V)=1$, $U,V\ge3$. The actual rows are respectively multiples of $U,V$. At any assigned predecessor power their three actual residues are $0,0,1$ in some order. If $W$ has prime factors, absorb any remaining powers of those bases in $D$ into $W$ so their powers are full. This replaces $D$ by a divisor. Both size bounds persist and $Ps-1\nmid D$ persists. Use the resulting $D,W$ throughout the proof. [Endpoint reduction](all_weights_endpoint_reduction.md) · [Two predecessor powers](two_prime_powers.md)

When (1) holds, divide the parts of the row assigned $p$ by $P$. By Kummer's criterion their low digits partition the two units in $1+r$. Each least residue modulo $B$ is one of $0,1,r,1+r$, and at most two of the three residues can be nonzero. This last fact uses carry freedom at every intervening digit.

## Both predecessor bases assigned to one row

Call this row $A$, of scale $U$. Its zero positions for these two powers intersect, so one part is $kW$, with $U\mid k$ and $0<k<D$. This argument also applies to a single predecessor base or $W=1$.

A prime $p$ satisfying (1) cannot be assigned to the companion row $C=V(f_1,f_2,f_3)$. The three positive $f_i$ sum to $U\le k<D<B$. Since $V/P$ would be a unit modulo $B$, all three low residues after division by $P$ would be nonzero.

Thus $P\mid U$ and $k=Pz$. The low residue $e$ of $zW$ is positive, since $0<z<B$ and $p\nmid W$. The congruences $DW\equiv Ps-1$ and $zW\equiv e\bmod B$ give $De\equiv z(Ps-1)\bmod B$. Both sides are nonnegative and below $B$ because $e\le s$, $Pz<D$ and $Ds<B$. They coincide. In particular $e<s$, so $e=1$ or $r$. We have proved the exact identity
$$
 \frac{k}{D}=\frac{Pe}{Ps-1}. \tag{3}
$$
The fraction on the right is reduced, with numerator a positive power of $p$. Two distinct primes satisfying (1) give incompatible reduced numerators for the same $k/D$, even if their ratios agree. For just one prime, (3) requires $Ps-1\mid D$. Writing $D=d(Ps-1)$ gives $k=dPe$. Since $U\mid k$ and $\gcd(U,d)=1$, the unitary weight scale $U$ equals $P$. This is impossible if $P=2$.

## The predecessor bases assigned to different rows

Write $W=QR$ for their full powers. Two parts of $A$ are $UQx,UQy$, and two of $C$ are $VRz,VRt$. All multipliers are positive. The remaining positive parts give $V\ge Q(x+y)+1$ and $U\ge R(z+t)+1$, whence
$$
 (x+y)(z+t)<D. \tag{4}
$$
Every full power of $D$ divides one of the four multipliers, using the two zero positions in its assigned row and coprimality to $UQ,VR$.

Suppose $A$ takes a weight prime satisfying (1). The two low residues of $(U/P)Qx,(U/P)Qy$ modulo $B$ are positive since the scalar is a unit and $0<x,y<D<B$. They must be $1,r$ in some order. Name them $a,b$. Then $xb\equiv ya\bmod B$ and $xb,ya<Dr<B$. Thus $xb=ya$ and the two multipliers are $x,rx$ after interchanging slots. A prime assigned to $C$ similarly determines the unordered ratio of its displayed pair. [Lean verification of both exact lifts](TwoDigitWeightCuts.lean)

Two primes giving different ratios cannot be in the same row. If they go to different rows, write the two pairs as $x,r_p x$ and $z,r_q z$. Both ratios are coprime to $D$. Every full power of $D$ must then divide $x$ or $z$, so $D\mid xz$, contrary to (4). This completes the two-prime theorem.

For the one-prime bound assume $N$ has at least two distinct prime factors, since a prime-power weight has a universal prime. For each other prime $\ell\mid N$ put $P_\ell=\ell^{v_\ell(N)}$ and define
$$
 G=\prod_{\ell\mid N,\ \ell\ne p}\ell^{v_\ell(N/P_\ell-1)}.
$$
Each factor greater than one divides a displayed multiplier in the row taking $\ell$, since modulo this power the normalized low digits consist of just one unit. If this row is $A$, the factor divides $x$ since it is coprime to $r$. As the factors are pairwise coprime and coprime to $D$, $DG\mid xzt$. Inequality (4) now gives
$$
 sG(z+t)<zt\le (z+t)^2/4.
$$
Thus $z+t\ge4sG+1$ and $D>sx(z+t)\ge s(4sG+1)$. The grouped case was already excluded under (2).

## Numerical checks

At $N=1{,}255{,}603{,}200{,}006$, take $D=385$, $W=17\cdot191841589$. The least residues of $N/2\bmod 2048$ and $N/3\bmod 2187$ are 3 and 4, respectively. Their forced ratios are 2 and 3. A pair of assignments survives the earlier one-unit residue constraints at this weight. The two-prime argument settles it using these two prefixes.

Run `python3 check_two_digit_weight_cuts.py` and `lean TwoDigitWeightCuts.lean` (Lean 4.36). The exact coefficient check visits 396,102 rows and tests 3,234 pairs of endpoint gcd values at 94 selected weights through 360. On the 100,000 weights starting at $1{,}001{,}234{,}567{,}890$, the tests certify 19,289 weights after omitting weights with a universal weight prime or fewer than two predecessor bases. Of these, 1,476 passed none of the earlier one-unit inequalities. The full interval was independently verified by two exhaustive searches. [Verifier](check_two_digit_weight_cuts.py) · [Full interval checks](predecessor_weight_gap.md#all-100000-weights-in-the-interval)
