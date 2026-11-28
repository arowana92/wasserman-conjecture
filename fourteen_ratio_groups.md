# Two checks through height fourteen

Let $A$ have two parts in reduced ratio $a:b$, $1\le a\le b\le14$. Let each pair gcd in the primitive direction of $B$ be at most 14. For every common weight $N\le10^{1000}$ their coefficients share a prime dividing $N(N-1)$. Over all weights, at most $O_\epsilon(X^\epsilon)$ weights $N\le X$ could fail under these hypotheses, for every $\epsilon>0$. A full homogeneous relation on $B$ with coefficient magnitudes at most 14 implies the stated pair gcd bound.

Some unequal height bounds give the common prime at **every** weight. The following are sufficient bounds from this computation.

| Maximum ratio term in $A$ | Maximum full relation height on $B$ | Alternative maximum primitive pair gcd on $B$ |
|---:|---:|---:|
| 6 | 14 | 14 |
| 8 | 12 | 12 |
| 10 | 12 | 10 |
| 12 | 12 | 9 |
| 14 | 6 | 6 |

Apply the [exact denominator and group enumeration](mixed_relations_twelve.md) with height 14. There are 425,908 distinct tuples $(a,b,D,U,x)$ after the pair ratio bounds. After enumerating $j$ and checking carry graphs one prime at a time, only [forty families](height14_ratio_families.txt) remain. The all-weight table follows by checking the pair gcds of the primitive normalization of each remaining $x$, and testing for $rx_i+sx_j=tx_k$ in all three orientations within the indicated height.

## Verifying the weights

In each family the shapes are $A=(aUjm,bUjm,(D-(a+b)Uj)m+1)$ and $B=Vx$, with $N=Dm+1=UV$. For each prime in $DU$, and also for 2, the group conditions specify which rows can avoid the prime. In particular, if $p^e\parallel D$ is assigned to $A$, exactly one term $c$ of its reduced ratio is a unit at $p$, and $cUjm=1\bmod p^e$.

The forward check splits the 40 families into 46 choices of avoiding rows. It combines the unit congruences and $U\mid N$ by CRT and runs up each resulting progression. A positive factorial-floor difference at a power $q$ forces a carry until a denominator floor increases. The backward check uses only $N=1\bmod D$ and $N=0\bmod U$. It runs down this larger progression without choosing an assignment. If each eligible row at one prime has a positive factorial-floor difference, the entire interval on which the corresponding numerator floors stay fixed is excluded. The two runs take 308,612 and 295,425 intervals, respectively. Both find zero weights through $10^{1000}$.

## Beyond the finite interval

In each assignment a multiplier is constrained in at least two independent prime bases. For row $A$ the multiplier is $m$ and for row $B$ it is $V$. Each is affine in $N$ with fixed positive slope. As in the [block argument](six_coefficient_row_relations.md#the-count-beyond-the-finite-range), the carry graph in a required base $p$ gives a digit-block set. Using 64-digit blocks, the verifier selects numerators $k_p$ with $W_p^{1000}<p^{64k_p}$, where $W_p$ counts all paths of this length from reachable states.

In all 46 assignments a subset of $h$ required primes has $\sum k_p\le1000(h-1)-185$. For $N\le X$, normalize each multiplier by a power of its own block base comparable to $X$. The corresponding possible $N/X$ lie in dilates and translates of all the digit sets, with scales bounded above and below by positive constants independent of $X$. The uniform intersection estimate therefore gives $O_\epsilon(X^\epsilon)$ integer weights in each assignment. [Corso and Shmerkin, Corollary 1.17](https://arxiv.org/pdf/2409.04608)

The [verifier](check_height_fourteen_scales.py) also rebuilds the full finite family list with the following command.

~~~sh
python3 check_height_fourteen_scales.py --rebuild
~~~
