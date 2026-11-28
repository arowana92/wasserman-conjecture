# The sharp predecessor divisor on a primitive ray

Fix a positive primitive triple $y$, let $u=\sum_i y_i$, and write $K(V)$ for the coefficient of the row $Vy$ at weight $N=uV$. Let $E_y$ be the least positive integer for which
$$
                         N-1\mid E_y K(V)\qquad(V\ge1).
$$
It exists and divides $H=\operatorname{lcm}_{i<j}\gcd(y_i,y_j)$ by the [retained-power argument](endpoint_pair_bound.md). We can compute its exact value, with an individual multiplier witnessing each positive prime exponent.

At $p^e\parallel H$ take the unique unit coordinate $y_i$ at $p$. Choose another coordinate $y_j$ with $v_p(y_j)=e$. Such a choice exists since the two nonunit valuations have minimum $e$. Let $\delta_i=\delta_j=1$ and let the remaining $\delta$ be zero. If $L_p(V)$ is the $p$-valuation of the multinomial on the three nonnegative parts $y_kV-\delta_k$, then
$$
 v_p(E_y)=\max\left(0,\ e-\min_{p\nmid V} L_p(V)\right).              \tag{1}
$$
Indeed factorial decrement gives
$$
       \frac{K(V)}{uV-1}=\frac{uV}{(y_iV)(y_jV)}
                  \binom{uV-2}{y_1V-\delta_1,\ y_2V-\delta_2,\ y_3V-\delta_3}.
$$
The total $u$ is a unit at $p$. For a unit multiplier the defect $v_p(uV-1)-v_p(K(V))$ is exactly $e-L_p(V)$. When $p\mid V$ it is nonpositive. This proves (1).

In particular
$$
 u\,\frac{N-1}{\gcd(N-1,E_y)}\ \mid\ \gcd(N(N-1),K(V)).             \tag{2}
$$
Here $u\mid K(V)$ follows by factorial decrement in each coordinate and primitivity. For two rays the predecessor part of the common gcd is divisible by $(N-1)/\gcd(N-1,\operatorname{lcm}(E_y,E_z))$. Formula (2) implies that a single row with $E_y\le u$ has endpoint gcd at least $N$, since $\gcd(E_y,u)=1$. Two such rows therefore share an endpoint prime at every common weight.

## Deciding a prime power

For $1\le j\le e$ a defect of at least $j$ requires $V=r+p^j n$ where $r$ is the least positive residue of $u^{-1}$ modulo $p^j$. The low $j$ digits of the decremented parts have residues $0,p^j-1,0$: after decrement the unit coordinate becomes zero and the other decremented coordinate becomes minus one. There are no carries in those positions. On the higher digits the multiplication carries initially equal $\lfloor(y_i r-\delta_i)/p^j\rfloor$, and the addition carry is zero.

Reading an input digit $d$ changes the multiplication carries $c_i$ to $\lfloor(c_i+y_i d)/p\rfloor$. Add the three remainders and the incoming addition carry. The quotient by $p$ is the outgoing carry and the edge cost. The desired multiplier exists precisely when the all-zero state is reachable with total cost at most $e-j$. Each finite accepting path is an ordinary positive multiplier and every such multiplier gives a path after adjoining high zeros.

There are fewer than $6u$ reachable states. At a prefix $w$ of length $s$ the multiplication carries are $\lfloor(y_i w-\delta_i)/p^s\rfloor$. As functions of $w/p^s$, they are constant on each open interval between consecutive breakpoints $b/y_i$ in $[0,1]$. The subtraction can alter a floor only at an exact breakpoint because $y_iw$ is an integer. Each breakpoint itself has only one carry vector. There are fewer than $2u$ vectors and three possible addition carries. Nonnegative edge costs allow loops to be deleted from a witness. In particular, when the exponent $j>0$ is optimal, a witnessing multiplier exists below $p^{j+6u}$.

The [main computation](check_direction_predecessor_loss.py) uses this shortest-path test. The [independent computation](check_direction_predecessor_cross.py) reads the **original** row. Following the fixed low suffix, a further zero digit of the weight continues to earn one unit of predecessor valuation until the first mismatch. Each outgoing carry spends its value. The available balance starts at zero after reserving $j$ units. Before depth $e$ a matching digit has no carry. Beyond depth $e$, every further matching layer has at least one carry, since the coordinate of valuation exactly $e$ cannot vanish modulo that layer. After a mismatch the balance can only decrease. Thus a negative balance may be discarded and the balance never exceeds $e-j$. Reaching zero residuals after a mismatch proves the same defect threshold directly.

Through primitive sum 90 the two methods agree in 16,717 individual prime-power comparisons. In 8,068 directions $E_y$ is strictly less than $H$. Another check uses 9,650 exact coefficients. The smallest sum for which $E_y>u$ is 127, attained at $(10,52,65)$ with $E_y=130$. The factors 2, 5 and 13 are witnessed at multipliers 1, 63 and 4 respectively. These three witnesses need not describe the same scale.

## Density on a ray

At every fixed primitive direction, $u(N-1)\mid K(V)$ on a set of multipliers of natural density one. More precisely, at $p^e\parallel H$ the number of scales with defect at least $j>0$ in any arithmetic progression of length $X$ and positive step is
$$
 O_{p,e}\bigl(X^{b_p}(1+\log X)^{e-j}\bigr),\qquad
 b_p=\log_p\bigl(1/(2\sin(\pi/(4p+2)))\bigr)<1 .
$$
Indeed two decremented terms have slopes of valuations $0,e$. Their digit sums can overflow at at most $e-j$ positions by (1). Apply the [weighted digit count](direction_carry_counts.md#a-prime-digit-bound-for-multiples-in-a-progression) to these two affine terms. The constant offsets change the initially known digits and preserve the digit permutations and the delay $e$. At an unrestricted position the sum of successor weights is bounded by $p\rho_p$ times the parent weight, where $\rho_p=p^{b_p}$, by the displayed weight ratios in that proof. At every other position it is at most $\rho_p$ times the parent weight. Choosing at most $e-j$ unrestricted positions among $\lceil\log_p X\rceil$ gives the bound. A finite union over primes of $H$ proves the density assertion.

~~~sh
python3 check_direction_predecessor_loss.py --limit 512
python3 check_direction_predecessor_loss.py --row 10 52 65
python3 check_direction_predecessor_cross.py --limit 90
~~~
