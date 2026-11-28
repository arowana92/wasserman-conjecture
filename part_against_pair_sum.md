# A part divided by the sum of the other two

Let $B=(x,y,z)$ be a positive trinomial row of weight $N$ and write $z/(x+y)=R/S$ in lowest terms. Let $A$ be **any** other row of weight $N$. If their coefficients have no common prime dividing $N(N-1)$, then the following inequalities are necessary.

* If $R=p^k$ with $k\ge1$, then $N>\kappa_p(R/S)^3$, where $\kappa_2=729/4$, $\kappa_3=192$, and $\kappa_p=216$ for $p\ge5$.
* If $S=1$, then $N>216R^{5/2}$ for every $R$. For every prime power $P$ exactly dividing $R$, one also has $N>(729/4) R^2P$.

In particular a nonendpoint row in a counterexample whose reduced numerator is a prime power satisfies
$$
              x+y>\frac{N}{1+(N/\kappa_p)^{1/3}}.
$$
These are pointwise assertions with no restrictions on the other row. They apply at all weights by the [endpoint reduction](all_weights_endpoint_reduction.md).

## A general gcd bound

Let $u,v$ be the primitive sums of $B,A$ and $H,K$ the lcms of the three pair gcds of their respective primitive directions. Put $M=N-1$, $D=\gcd(N,K_A,K_B)$ and $J=\gcd(M,K_A,K_B)$, where $K_A,K_B$ denote the two coefficients. Then
$$
              27 M^2 u^3\le (D N)^3 (J H)^2.                           \tag{1}
$$
Indeed $u,v$ divide $N$ and each divides its own coefficient, so $\gcd(u,v)\mid D$ and $uv\le DN$. The [retained powers](endpoint_pair_bound.md) give $M\le JHK$ because every remaining power divides a pair gcd in an assigned row. The three pair gcds of a primitive triple are pairwise coprime, hence $K^2$ divides the product of its entries. AM–GM gives $27K^2\le v^3$. These three inequalities yield (1).

## Optimizing the dependence on the slope

Work with the primitive direction of $B$. Write $d,e,f$ for the pair gcds opposite the coordinates $x,y,z$, respectively. From $Rx+Ry=Sz$ we have $d,e\mid R$, $f\mid S$, $\gcd(d,e)=1$, $H=def$, and $x+y\ge f(d+e)$. Thus
$$
 \frac{u^3}{H^2}\ge\left(\frac{R+S}{S}\right)^3
          \frac{(d+e)^3}{d^2 e^2}.
$$
Let $\psi(R)$ be the minimum of the last factor for coprime positive $d,e$ with $de\mid R$. For $S=1$ this lower bound is sharp as a bound on $u^3/H^2$. For every admissible pair take the primitive triple $(e,d,R(d+e))$, whose pair-gcd lcm is $de$.

If $D=J=1$, (1) gives $N\ge27(1-1/N)^2 u^3/H^2$. Also $R+S\mid N$, so $(1-1/N)^2(R+S)^3\ge(R+S-1)^2(R+S)>R^3$. Consequently
$$
                       N>27(R/S)^3\psi(R).                             \tag{2}
$$
For every $R$, AM–GM gives $\psi(R)\ge8/\sqrt R$. Also let $m(R)$ be the largest possible $\min(d,e)$ in the same search. For $d\le e$ the polynomial identity
$$
          4(d+e)^3-27d e^2=(e-2d)^2(4e+d)
$$
gives $\psi(R)\ge27/(4m(R))$. For a full prime power $P\parallel R$, at least one of $d,e$ is coprime to $P$ and divides $R/P$, so $m(R)\le R/P$. This proves the second assertion.

For $R=p^k$ one of $d,e$ is one. The function $(1+w)^3/w^2$ decreases for $w\le2$ and increases for $w\ge2$. Among $w=1,p,\ldots,p^k$ its minimum is $27/4$ for $p=2$, $64/9$ for $p=3$, and $8$ for $p\ge5$. Equation (2) gives the first assertion.

[Lean](IntegralPartSlope.lean) verifies the composition of the three arithmetic hypotheses in (1) and the polynomial inequalities in the optimization. The factorial and primitive-row inputs are proved above. The [exact checker](check_integral_part_slope.py) tests (1) on 730,275 ordered row pairs through weight 55, checks the primitive slope bound 13,482 times, and verifies the optimizing directions for all integral slopes through 1,024.

~~~sh
lean IntegralPartSlope.lean
python3 check_integral_part_slope.py
~~~

The Lean file uses Lean 4.36 and its standard library.
