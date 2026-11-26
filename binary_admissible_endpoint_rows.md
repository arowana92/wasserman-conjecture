# An infinite lowest-layer family with a fully odd trinomial coefficient

There are infinitely many explicit weights and pairs of positive rows $A,B$ satisfying both of the following tests for a counterexample with the endpoint $(1,1,N-2)$.

1. Every full prime power of $N$ or $N-1$ can be assigned to a row whose residues sum to the residue of $N$ at that power.
2. The entire binary addition in $A$ is carry-free.

Nevertheless **3 divides both row coefficients and the endpoint coefficient**. The failed ternary layer in $B$ is 9. This strengthens the [lowest-layer example for rational proportions](proportional_part_pairs.md#an-infinite-family-passing-all-lowest-endpoint-layers). It identifies a limit of a search using just the full powers at the two endpoints and a complete test at 2.

## The weights and rows

Take any integer $h\ge1$ with $h\not\equiv136\bmod369$. Put
$$k=28+840h,\quad Z=2^k,\quad C=187,\quad
D=2^9+2^{47}+2^{58}+2^{89},\quad b=D/C=3310000106583215246468608.$$
Define
$$V=1+CZ,\quad U=1+DZ,\quad N=UV,\quad S=CZ,\quad T=1+bV,\quad l=T/35.$$
The rows, in the displayed coordinate order, are
$$A=U(18Z+1,33Z,136Z),\qquad B=V(23l-b,5l,7l).$$
They are positive, have normalized sums $V,U$ respectively, and $N-1=ST$. In terms of the normalized [quadratic shapes](proportional_part_pairs.md#both-nonexact-rows-in-square-root-windows), $\alpha=1$ and $\beta=-b$.

Since $2^{840}\equiv1\bmod 187\cdot1225$ and modulo 9, $V\equiv123\bmod1225$. Also $b\equiv908\bmod1225$. Thus $l$ is integral with $l\equiv6\bmod35$. We have $\gcd(S,T)=1$ since $b$ is even and $\gcd(C,b+1)=1$.

For $\gcd(U,V)=1$ use $U=bV+1-b$ and
$$b-1=3^3\cdot67\cdot389\cdot11071\cdot424866988899517.$$
For any listed prime $p$, a common divisor requires $(2^{840})^h\equiv r_p\bmod p$ where $r_p=-1/(C2^{28})\bmod p$. Put $e_p=(p-1)/\gcd(p-1,840)$. Fermat's theorem requires $r_p^{e_p}\equiv1$. The following are least nonnegative residues.

| $p$ | $e_p$ | $r_p^{e_p}\bmod p$ |
|---:|---:|---:|
| 3 | 1 | 2 |
| 67 | 11 | 37 |
| 389 | 97 | 388 |
| 11071 | 369 | 1 |
| 424866988899517 | 35405582408293 | 164056555269709 |

At $11071$ the verifier tries all 369 residues of $h$ and finds exactly 136. All prime factors in this test are recursively certified.

## The carry checks

The full powers in $U,V$ have zero residues throughout the corresponding row. In $S$, the binary power $Z$ gives residues $(1,0,0)$ in $A$. At 11 and 17 the exceptional positions of $A$ are three and two, as is seen by substituting $Z\equiv2^{28}\bmod187$. In $T$, the factor $l$ is coprime to 35. Its second and third positions in $B$ vanish, and its first is $-bV\equiv1\bmod l$. At 5 and 7 the exceptional positions in $B$ are three and two respectively. This checks every full endpoint power, including those with unknown prime base in the growing factors.

For all binary columns, let $c=(18,33,136)$ and $e_1=(1,0,0)$. Expansion gives
$$A=e_1+Z(c+De_1)+Z^2 D c.$$
The binary supports of the three $c_i$ are disjoint and lie in positions zero through seven. The four bits of $D$ are above this range, with successive gaps at least eight. Therefore the three entries in the middle block are disjoint. The three entries of $Dc$ are also disjoint, since its copies of all the $c_i$ occupy four separate eight-bit blocks. All three coefficient blocks fit below bit 100, whereas $k\ge868$. Hence $A$ is odd in coefficient.

Finally $N\equiv7\bmod9$. The three residues of $A$ modulo 3 sum to 4. The three residues of $B$ modulo 9 sum to 16. Both exceed the corresponding residues of the weight, so the coefficients are divisible by 3. These comparisons depend only on $Z\bmod9$ since $35$ is invertible there. In particular the higher ternary carry survives even in the first member, which has a weight of 1,833 bits.

Run `python3 check_binary_admissible_endpoint_rows.py` with python-flint installed. The checks include full binary digit sums for concrete weights, factorial valuations in both bases, and exact modular periods for the whole progression.
