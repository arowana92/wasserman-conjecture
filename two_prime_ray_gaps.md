# Divisibility by two or three near two rational triples

Let $W,E\ge1$ be integers. Suppose the parts of a positive trinomial row differ coordinatewise by at most $E$ from $(W,3W,6W)$ or $(2W,3W,5W)$, up to permutation. If its coefficient is coprime to six, then
$$
                   W\ll E\,(\log(2W))^{4.117}. \tag{1}
$$
The implied constant is absolute. In particular, if $E=o(W/(\log W)^{4.117})$, every sufficiently large such row has coefficient divisible by 2 or 3. Two explicit sufficient conditions for that divisibility are as follows. Write $B=\log_2(10W)$.

| Range | Sufficient condition |
|---|---|
| $B\le10^{100}$ | $W\ge 11400\,E B$ |
| $B\le10^{1000}$ | $W\ge662000\,E B$ |

At every weight divisible by three, the endpoint coefficient and the coefficient of the equal row $(N/3,N/3,N/3)$ are both multiples of six. In the equal row the first occupied digits already force carries at 2 and 3. Thus these tests give a common prime from $\{2,3\}$ whenever the third row satisfies the sufficient condition. The all-weight statement about one nearly proportional row uses no factorization of its weight. It allows errors of either sign in all three positions. For specific signed errors a single prime can already suffice at much larger errors. [Leading-digit tests](nested_predecessor_cofactors.md)

## The high parts of an addition

For a base $p=2,3$ write the slopes as $a=(a_1,a_2,a_3)$ in the stated order and choose the least $L=p^m\ge E$. Above the first $m$ places the parts have the form $a_i v+c_i$, where $v=\lfloor W/L\rfloor$ and $-1\le c_i\le a_i$. If the entire row is carry-free, so is the addition of these high parts.

On reading a digit $b$ of $v$ from the units upwards, set $t_i=c_i+b a_i$. Keep the step if $\sum_i[t_i]_p<p$, and replace $c_i$ by $\lfloor t_i/p\rfloor$. The box of states just given is closed under these steps, including for negative offsets. A valid finite word reaches $(0,0,0)$ after padding with zeros. Remove states that cannot reach this endpoint by any finite word. Next replace the retained set by the set of its successors inside itself, until the set stops changing. After $g$ such strict reductions, every accepted word has entered this last set at its $g$-th step (padding is allowed). The suffix is $y=\lfloor W/(Lp^g)\rfloor$.

The table gives all states in that last set which admit a positive suffix. Beside each state is an integer $h$ such that **every** accepted positive suffix at this state has $d y+h$ equal to a power of $p$.

| Slopes | $p$ | $g$ | $d$ | State and $h$ |
|---|---:|---:|---:|---|
| $(1,3,6)$ | 2 | 4 | 3 | $(0,1,2):1,\ (0,2,4):2$ |
| $(1,3,6)$ | 3 | 3 | 2 | $(0,1,3):1$ |
| $(2,3,5)$ | 2 | 3 | 5 | $(0,0,1):1,\ (1,1,3):3,\ (1,2,4):4,\ (0,1,2):2$ |
| $(2,3,5)$ | 3 | 3 | 5 | $(0,0,1):1,\ (0,1,2):2,\ (1,2,4):4,\ (1,1,3):3$ |

The verifier enumerates the finite graphs and checks two identities giving the arbitrary-length conclusion. For an allowed digit $b$ leading from $c$ to $c'$ with positive remaining suffix, $d b+h(c)=p h(c')$. If the remaining suffix can be zero and $b>0$, $db+h(c)$ is itself a power. Induction gives the claim for $y=b+p y'$. For $y=0$ one can take the power $1$ and offset $1$ separately. Since $1\le h\le d$ and $L<pE$,
$$
              |dW-p^i|\le d p^{g+1}E
$$
for some $i\ge0$. For the first triple this gives a power of two within $96E$ of $3W$ and a power of three within $162E$ of $2W$. For the second triple the two powers lie within $80E$ and $405E$ respectively of $5W$.

## Separating the two powers

Take $W\ge200E$. In the first case multiply the binary power by 2 and the ternary power by 3 to obtain powers $2^m,3^n$ within $192E,486E$ of $6W$. Their difference is at most $678E$, and the smaller power exceeds $3.5W$. In the second case their difference is at most $485E$ and the smaller exceeds $2.9W$. In either case the exponents are positive, $m<B$, and the mean-value theorem gives
$$
  0<|m\log 2-n\log 3|<200E/W . \tag{2}
$$
The distance of $m\log 2/\log 3$ to the nearest integer is $\gg m^{-4.117}$. This is the two-logarithm consequence of Wang and Wu used by Fan, Queffélec and Queffélec, with exponent rounded upwards. Since $\log3>1$, combining this with (2) proves (1). The case $W<200E$ is immediate. [Wang and Wu](https://doi.org/10.1016/j.jnt.2014.03.007) · [Fan, Queffélec and Queffélec, Section 2.1](https://arxiv.org/pdf/2303.06850)

For the explicit bounds, put $\alpha=\log 2/\log 3$. Given $q_j\le m<q_{j+1}$ for consecutive convergent denominators, the best-approximation property gives $\|m\alpha\|\ge\|q_j\alpha\|>1/(q_j+q_{j+1})$. The verifier encloses both logarithms using rational sums of the $\operatorname{atanh}$ series with a geometric tail bound. The upper and lower rational bounds for their ratio have identical continued-fraction quotients until the stated denominator bound is crossed. The maximum of $\lceil (q_j+q_{j+1})/q_j\rceil$ is 57 up to $10^{100}$ and 3,310 up to $10^{1000}$, including the step to the first denominator above each bound. Thus $\|m\alpha\|>1/(57m)$ or $1/(3310m)$ respectively, proving the two explicit conditions by (2).

Run **python3 check_two_prime_ray_gaps.py** for the transition and suffix certificates, direct tests of exact affine parts at small and large arguments, and the rational continued-fraction bounds. [Standard-library verifier](check_two_prime_ray_gaps.py)
