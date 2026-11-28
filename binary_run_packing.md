# Binary digit collisions in scaled trinomial rows

Let $h\ge3$, $1\le b\le h-2$, $A=2^{h-2}$, and put
$$
   f=(A,2A,A-2^b+1),\qquad S=4A-2^b+1,\qquad
   D=\{0\}\cup\{b,\ldots,h-1\}.
$$
For every positive $z$, the trinomial coefficient with parts $zf$ is odd **if and only if** every positive difference between two occupied binary positions of $z$ lies outside
$$
                    [1,h-b-1]\ \cup\ [b,h-1]. \tag{1}
$$
Equivalently, the translates of $D$ at these positions are pairwise disjoint, or $s_2(Sz)=(h-b+1)s_2(z)$, where $s_2$ counts ones. The previously proved [minimum separation of $h$](fermat_inverse_powers.md#a-direct-spacing-theorem) is the case $2b\le h$. For example, when $h=10,b=7$, a gap of four is allowed, but three ones at positions $0,4,8$ are forbidden by the distance eight. All pairs matter.

## Proof of the parity criterion

The supports of the three entries of $f$ partition $D$. The positive differences of elements of $D$ are exactly (1). If translates are disjoint, all the bits from all the products are distinct, and the addition is carry-free. Also, the digit sum of the product $Sz$ attains its upper bound precisely when there is no carry among the translates, since every binary carry reduces the number of units. We prove necessity using the binary carry criterion for the trinomial. [Bergman, Lemma 5](https://arxiv.org/abs/0806.0607)

Suppose the addition of the three entries of $zf$ is carry-free. The disjoint supports of $Az$ and $2Az$ imply that $z$ has no adjacent ones. Assume two translates of $D$ collide. Choose the largest $i$ that is the **smaller** of the two positions in a colliding pair. Write $z=H+L$, with $H$ formed from positions above $i$ and $L$ from those at or below $i$. All translates in $H$ are disjoint. The geometric series for alternate bits gives
$$
  2^i\le L<\tfrac43\,2^i,\qquad H\bmod 2^n<\tfrac23\,2^n\quad(n\ge1). \tag{2}
$$

First suppose $H$ has a one at $j$ with $j-i<h-b$. Since $j-i\ge2$, writing $B=2^{j+h}$ gives
$$
 2^{j+b}\le 2^{i+h-1}<SL<2^{i+h+1}\le B/2.
$$
Here $2^{h-1}<S<2^h$ and (2) give the bounds on $SL$. The positions $j+b,\ldots,j+h-1$ of $SH$ are all occupied. Thus $SH\bmod B\ge B-2^{j+b}$. Adding $SL$ crosses $B$, and $(Sz)\bmod B<SL<B/2$. This erases the bit at $j+h-1$. That bit belongs to $2Az$ and must be occupied in the carry-free total, a contradiction.

Hence every run $[j+b,j+h-1]$ belonging to $H$ starts at or above $i+h$. The collision with the translate at $i$ must therefore come from a singleton $j$ of $H$ in $[i+b,i+h-1]$. Put $T=A2^i$ and $c=A-2^b+1$. Modulo $4T=2^{i+h}$, the residue of $cH$ is $v=H\bmod 4T$. Indeed, the support of $c$ is $\{0\}\cup[b,h-3]$ (the interval may be empty), and all the translates in $H$ are disjoint. The singleton gives $v\ge2^{i+b}$. It follows from (2) and $c\le A$ that
$$
 T<2^{i+b}+c2^i\le v+cL
     <\tfrac23(4T)+\tfrac43 T=4T .
$$
Consequently $cz\bmod4T=v+cL$ has a bit at $T$ or $2T$. These two positions are already occupied by $Az$ and $2Az$ respectively, because $z$ has a one at $i$. This contradicts the hypothesis and proves (1).

## An obstruction for every weight with a specified divisor

For $1\le b<h-2$, let $U=4A-2^b+1$ and write $N=UV$. If the odd part $N_0$ of $N$ satisfies $N_0\equiv1\bmod A$, **every positive trinomial row at $N$ whose parts are multiples of $V$ has even coefficient**.

Indeed, write such a row as $V(g_1,g_2,g_3)$ and let $z$ be the odd part of $V$. Removing the common power of two preserves parity by Kummer. Now $z\sum g_i=zU=N_0$. If the coefficient is odd, only one low bit below $A$ is occupied in the total. At least two $g_i$ are therefore divisible by $A$, since $z$ is odd. They must be distinct, because equal positive parts share a bit. As $3A<U<4A$ the two are $A,2A$, leaving $A-2^b+1$. But
$$
        z\equiv(1-2^b)^{-1}\equiv 1+2^b\pmod {2^{b+1}} .
$$
Here $2^{b+1}\mid A$. Bits zero and $b$ violate (1). The excluded endpoint $b=h-2$ actually permits $z=1$ and the odd row $(A,2A,1)$.

## Most small inverse residues are impossible

Let $N=1+PW$ be even, with $P$ an odd prime power and $W\ge1$ having at most two distinct prime divisors. Suppose $A=2^m$, $m\ge2$, and the odd part of $N$ is $1\bmod A$. At any failure, $\gcd(P,W)=1$. Let $x$ be the least positive inverse of $W$ modulo $P$. If $x$ is odd and $x<4A$, a failure requires
$$
       3A<x<4A,\qquad x^{-1}\bmod A
       \text{ has no adjacent binary ones}. \tag{3}
$$
Further, all $x=4A-2^b+1$ with $1\le b<m$ are impossible.

By the [one-numerator reduction](static_square_cofactors.md#at-most-two-possible-numerators), a failure requires $3\le U\mid\gcd(x,N)$, $V=N/U$, and a companion $V(g_1,g_2,g_3)$ avoiding the primes assigned to $V$. As $x$ is odd, the power of two in $N$ belongs to $V$, so this coefficient is odd. The argument from the preceding section forces two different positive multiples of $A$ among the $g_i$, hence $U>3A$. Thus $U=x$ and the two multiples are $A,2A$. Disjointness of their scaled binary digits forbids adjacent ones in the odd part $z$ of $V$, and $z\equiv x^{-1}\bmod A$. This proves (3). The parity theorem handles the specified residues left by that test.

There are at most $F_m$ possible odd $x$ below $4A$, where $F_1=F_2=1$. In $(3A,4A)$ the odd inverses run bijectively over the odd residues modulo $A$. An odd $m$-bit word (leading zeros allowed) without adjacent ones ends in 01 and has $F_m$ choices. For $m\ge2$ the template further excludes $m-2$ of these candidates. Indeed, $x^{-1}\bmod A$ is $\sum_{0\le kb<m}2^{kb}$, so among $1\le b<m$ all $b\ge2$ survive the adjacent-ones test. Hence at most $F_m-(m-2)$ remain.

## Completing the forced low digits

For every odd $x\in(3A,4A)$, put $d=x-3A$. The table gives **all** $d$ for which some positive multiplier $z$ makes $z(A,2A,d)$ carry-free with $dz\equiv1\bmod A$. Any odd inverse outside the listed residues $x=3A+d$ is excluded under the hypotheses of (3).

| $m=\log_2 A$ | Admissible $d$ |
|---|---|
| 2 through 6 | 1 |
| 7 | 1, 57, 77 |
| 8 | 1, 25, 113 |
| 9 | 1, 113, 241 |
| 10 | 1, 141, 241, 481, 505, 653, 737, 793, 865, 945, 953 |

The input prefix is $z_0=d^{-1}\bmod A$. Only the third part can have a nonzero bit below $A$. After this prefix the three multiplication carries are $\lfloor f_i z_0/A\rfloor$, $f=(A,2A,d)$. Reading a new digit $e$ adds $e f_i$ to each carry. If the units sum to at most one, division by two gives the next carries. Each carry is less than its slope, so the graph is finite. A finite input completes the prefix precisely when one can reach a state with pairwise disjoint bits, which accepts zero as the remaining input.

A second algorithm enumerates all states $(\lfloor g_i u\rfloor)_i$ with $g=(A,2A,d,3A+d)$ for real $0\le u<1$. The left endpoints $u=j/g_i$ exhaust the cells. Keep states whose first three coordinates sum to the last. Starting at zero, invert each possible binary transition, allowing the output units to use at most one bit and requiring their sum to be the fourth unit. The initial prefix state occurs in this backwards search exactly when it can reach zero by a finite word. Actual input states are included by taking $u$ to be the prefix divided by the modulus. Conversely the inverted edges satisfy the exact carry recurrence, so any reverse path from a given actual prefix gives a valid suffix. Both algorithms give the table. For $m=7,d=57$, the multiplier $z=9$ gives the disjoint parts $1152,2304,513$.

Run **python3 check_binary_run_packing.py**. For every $3\le h\le16$ and $1\le b\le h-2$, the verifier also checks for each distance $1,\ldots,h$ whether **any** accepted path contains two ones at that distance. The 105 graphs give precisely the predicted distances, over all multiplier lengths. A separate computation of 107,520 factorial valuations and direct coefficients through weight 500 checks the parity equivalence on finite inputs.
