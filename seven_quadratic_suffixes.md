# A logarithmic exception bound at the first 10,000 cofactors

Let $P$ and $R$ be prime powers, $10P+7=DR$, and $1\le D\le10{,}000$. The full three-row conjecture holds at $(2P+1)(5P+1)$ for every $P\le 2^{65536}$. Without an upper bound, each dyadic block of inputs contains at most 24 possible exceptions across all these cofactors. Thus the number of potentially exceptional inputs up to $X$ is $O(\log(2X))$.

The [cofactor theorem](seven_quadratic_cofactors.md) settles every input at all but three cofactors. At the three possible shapes it gives $P$ odd, $T=2P+1$, and a companion of odd coefficient with parts
$$
       TRz,\qquad TRw,\qquad T\big((5P+1)-R(z+w)\big).
$$
The two remaining cofactors where $R$ must be prime and the one where it may be a higher power are all included here.

## The high parts of the binary addition

For $X=2^j\le P<2X$ put $L=16X$. Write $TR/(2L)=n+f$, $0\le f<1$. With $A=(2z,2w,D-2z-2w)$, divide each actual part by $L$ and round down. The resulting triple is $nA+s$ with
$$
 s=(\lfloor A_1 f\rfloor,\lfloor A_2 f\rfloor,\lfloor A_3 f\rfloor-e),
 \qquad e\in\{0,1\}.
$$
The last coordinate has lost $5T/(2L)<1$ before rounding. The upper triple adds without a binary carry. All allowed offsets $s$ are enumerated by taking $f=i/A_k$ at every breakpoint in $[0,1)$ and both values of $e$.

When the next bit of $n$ is $d$, put $v=s+dA$. Keep the transition if at most one coordinate of $v$ is odd, and pass to the coordinatewise floors $\lfloor v/2\rfloor$. The new offset still lies in the enumerated set: use the fractional part $(f+d)/2$ and either zero or one lost unit. After all bits have been read, accept when the remaining offset is nonnegative and adds without carry.

To count distinct words, start with this terminal set, read the leading 1 backwards through the transitions and then prepend transitions for each next bit, moving from most to least significant. Every nonempty subset so obtained admits that word. For each word length maintain the vector counting how many words end at each subset. Its update is exactly determined by the two transitions. The first repeated count vector certifies periodicity for every larger length.

| $D$ | Offsets | Reachable nonempty subsets | Initial count vectors | Period | Maximum words at one length | Words at each periodic length |
|---:|---:|---:|---:|---:|---:|---:|
| 4,587 | 8,222 | 138 | 0 | 138 | 1 | 1 |
| 7,383 | 9,192 | 283 | 28 | 1,166 | 5 | 4 |
| 9,453 | 11,272 | 68 | 0 | 68 | 1 | 1 |

Only length 28 has five words. Here length means the number of bits of the positive integer $n$.

## Returning to the input

We already know $P\le10^8$ works. In any larger block, $n=\lfloor(20P^2+24P+7)/(32DX)\rfloor$ satisfies $X/(2D)\le n<3X/D$. Indeed $X\ge2^{26}>8D$ and the numerator increases on the interval; its value at $2X-1$ is $80X^2-32X+3$. If $d=\lfloor\log_2D\rfloor$, the possible word lengths are among $j-d-1,j-d,j-d+1,j-d+2$. For a fixed word, at most one $P$ in that block has $D\mid10P+7$. Increasing $P$ by $D$ increases the numerator by $D(40P+24+20D)>32DX$. Thus beyond the length-28 exception the three cofactors together allow at most $4(1+4+1)=24$ inputs per block. The finite verification includes every shorter block.

The word and the parameter satisfy the exact square interval
$$
                Kn<(10P+6)^2\le K(n+1),\qquad K=160DX .
$$
The bin congruences in the cofactor theorem fix $R$ modulo $10ab$, hence $P$ modulo $Dab$. For every word in the four lengths we intersect the square interval with this progression and the input block. At the surviving integers we test the actual three binary expansions directly. This includes every eligible prime power, since the search allows all integers in the progression.

| $D$ | Upper words tested over all blocks | Integer inputs after the interval and congruences | Odd companions |
|---:|---:|---:|---:|
| 4,587 | 262,040 | 15 | 0 |
| 7,383 | 1,048,046 | 117 | 0 |
| 9,453 | 262,040 | 8 | 0 |

The first program maintains $s=\lfloor\sqrt{Kn}\rfloor$ and the exact remainder when appending a bit and doubling $K$. It independently recomputes integer square roots at spaced checkpoints. The second uses a common denominator for the offset breakpoints, constructs predecessor transitions by solving for each output column, and maintains the last point of the congruence progression below each square-root threshold. It confirms the counts $15,117,8$ and tests binary digit sums instead of pairwise bit intersections. Both searches cover $10^8<P<2^{65536}$. The endpoint $P=2^{65536}$ is even and is excluded by the cofactor theorem.

Run both independent Python programs with the same bit bound. The Lean file checks the strict residue inequality in the cofactor reduction, the quadratic square identity, the parameter spacing, and the two square-update identities. [First verifier](check_seven_quadratic_suffixes.py) · [Independent verifier](check_seven_quadratic_suffixes_cross.py) · [Lean lemmas](SevenQuadraticCofactor.lean)

```sh
python3 check_seven_quadratic_suffixes.py 65536
python3 check_seven_quadratic_suffixes_cross.py 65536
lean +leanprover/lean4:v4.36.0 SevenQuadraticCofactor.lean
```
