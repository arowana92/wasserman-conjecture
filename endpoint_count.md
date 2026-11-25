# A stronger bound on the number of exceptional weights

**Theorem.** The number of weights $N\le X$ admitting a three-row counterexample is $O(X^{15/16})$. The [shorter density-one proof](almost_all_weights.md) gives $O(X^{999/1000})$ using just a sparse binary part and a single modulus. Here we also count the possible factorizations when that modulus is small, and use some of the known binary digits of $N$ when it is large.

Work in a block $X/2\le N<X$ with $X=2^L$. Selections without the endpoint contribute $O_\epsilon(X^{11/12+\epsilon})$ by the prime-gap argument in the shorter proof. For an endpoint counterexample, choose a binary carry-free row $B$, and its factor assignment $UV=N$, $ST=N-1$ as in that proof. Write $K=US$. For each individual part $z$ of $B$, $U\mid z$, $S\mid z(z-1)$ and $z\ge3$. There are $X^{o(1)}$ possible pairs $(U,S)$ for a given $z$.

## Counting assignments with small $K$

The complementary row has parts $Vy_i$ summing to $VU$. Distribute the full powers in $T$ among three products $d_1,d_2,d_3$ by the position of the exceptional residue 1. A product may be 1. Then $d_1d_2d_3=T$, and $y_i$ is a multiple of $T/d_i$. Put $y_i=(T/d_i)z_i$, where here the $z_i$ denote positive integer multipliers in the complementary row. Order the bins so $d_1\le d_2\le d_3$. We have
$$V(z_1d_2d_3+z_2d_1d_3+z_3d_1d_2)-Sd_1d_2d_3=1.\tag{1}$$

Group the assignments into dyadic ranges $U\asymp U_0$, $S\asymp S_0$. Within one box, $V\asymp X/U_0$ and $T\asymp X/S_0$. Moreover
$$z_1z_2z_3=\frac{y_1y_2y_3}{T^2}\le\frac{U^3}{27T^2}\ll \frac{U_0^3 S_0^2}{X^2}.$$
The number of ordered multiplier triples is at most $X^{o(1)}$ times this last bound, by summing the divisor estimate over the possible products. Empty boxes make no contribution.

Count the remaining parameters in two ways.

1. Choose $S,V,d_1$; we have $d_1\le T^{1/3}$. Put $D=S d_1-V z_1$. Equation (1) implies
   $$(D d_2-V d_1 z_2)(D d_3-V d_1 z_3)=V^2 d_1^2 z_2 z_3-D.$$
   Both factors are positive. For instance, the first times $d_3$ is $V y_3-1>0$. This also proves $D>0$. Hence the right side is a positive integer of polynomial size in $X$. Its factor pairs leave at most $X^{o(1)}$ possibilities for $d_2,d_3$. After including the multipliers, the box count is at most $X^{-2/3+o(1)}U_0^2 S_0^{8/3}$.

2. Choose $S,d_1,d_2$; the number of ordered pairs with $d_1d_2\le T^{2/3}$ is at most $T^{2/3}X^{o(1)}$, using the upper end of the box for $T$. Set $A=z_1d_2+z_2d_1$, $B'=z_3d_1d_2$, and $C=S d_1d_2$. Equation (1) becomes $V(A d_3+B')=C d_3+1$. Thus $A d_3+B'$ divides $CB'-A$. This last integer is positive: $(C-VA)d_3=VB'-1>0$ and $CB'-A=(C-VA)B'+A(VB'-1)$. Each divisor determines at most one $d_3$ and then $V$. This gives the second box bound $X^{-4/3+o(1)}U_0^3 S_0^{7/3}$.

The geometric mean of these upper bounds is $X^{-1+o(1)}(U_0 S_0)^{5/2}$. In particular the number of weights admitting an assignment with $K\le Y$ is at most $X^{-1+o(1)}Y^{5/2}$, on summing over $O((\log X)^2)$ boxes. This calculation only used the low residue condition on the complementary row.

## Counting assignments with large $K$

Split the $L$ binary positions into $l=\lfloor10L/13\rfloor$ low positions and $h=L-l$ high positions. The small-modulus estimate with $Y=2^l$ contributes at most $X^{12/13+o(1)}$. Take $K\ge2^l$ in the remaining count. Fix the high part $\lfloor N/2^l\rfloor$, write $j$ for its number of ones, and let $k$ be the number of low ones of $N$. There are at most $\binom lk$ weights for these data.

Give each high one weight 1 and each low one weight 2. The three parts of the binary carry-free row partition the ones of $N$. One part $z$ therefore has weighted digit count at most $(j+2k)/3$. Its high ones lie among the fixed $j$ ones. The number of candidate values for this part is at most
$$2^{(j+2k)/3}(1+1/2)^j(1+1/4)^l
=A_0^j 4^{k/3}(5/4)^l,\qquad A_0=3/2^{2/3}.$$
Indeed the generating product counts each possible subset with weight $2^{-\text{weighted digit count}}$. For each $z$ there are only $X^{o(1)}$ divisor pairs $(U,S)$. With the high part of $N$ fixed, $N$ lies in an interval of $2^l$ consecutive integers. Its residue class modulo $K$ therefore contains at most one such $N$.

Bound the minimum of the binomial count and this divisor count by their weighted geometric mean, with powers $2/5$ and $3/5$, respectively. In the sum over $k$ use
$$\binom lk^{2/5}4^{k/5}=\big(\binom lk 2^k\big)^{2/5}\le 3^{2l/5}.$$
Summing over the high bit strings gives an upper bound
$$X^{o(1)}(1+A_0^{3/5})^h(5/4)^{3l/5}3^{2l/5} = X^{\beta+o(1)},$$
where
$$2^{13\beta}=\left(1+(27/4)^{1/5}\right)^3(5/4)^6 3^4.$$
Here $\beta<15/16$. For an exact check, $(733/500)^5>27/4$ and $((1233/500)^3(5/4)^6 3^4)^{16}<2^{195}$. Each of the three exponents $11/12$, $12/13$, and $\beta$ is strictly below $15/16$. The dyadic sum proves the theorem.

Run `python3 check_endpoint_count.py` for exact integer checks of the entropy comparisons in both proofs and the divisor identities on 20,268 assignments of prime powers, including examples with all three bins occupied.
