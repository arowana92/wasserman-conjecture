# Parts in different rows of a counterexample

**Theorem.** In any three-row counterexample of weight $N$, choose a part $a\ge2$ in one row and a part $b\ge2$ in a different row. Then
$$(ab)^5>27(N-1)^2.$$

In particular, at most one row can contain a part in the interval $[2,27^{1/10}(N-1)^{1/5}]$. If a row contains a 2, each of the six parts of the other rows exceeds $27^{1/5}(N-1)^{2/5}/2$.

**Proof.** Call the selected rows $A,B$ and the remaining one $C$. A counterexample has a row avoiding each prime. For primes of $N$, assign the full prime powers to avoiding rows, and call their three products $U_A,U_B,V$. For the primes of $N-1$, likewise call the assigned products $S_A,S_B,T$. Empty products are 1. Kummer’s rule on the lowest digits implies that every part of $A$ is a multiple of $U_A$, and similarly for the other rows. Thus
$$U=U_AU_B\le ab,\qquad C=(Vy_1,Vy_2,Vy_3),\qquad \sum_i y_i=U.$$

For each full power $w$ in $S_A$, the residues of the three actual parts of $A$ are $1,0,0$ in some order. Hence $w\mid a$ or $w\mid a-1$. Since $w$ is coprime to $U_A$, we obtain $S_A\mid(a/U_A)(a-1)$. Similarly $S_B\mid(b/U_B)(b-1)$. For each full power in $T$, two of the positive $y_i$ must be multiples of that power. Therefore $T^2\mid y_1y_2y_3\le(U/3)^3$. Combining the three estimates gives
$$N-1=S_AS_BT
\le\frac{a(a-1)b(b-1)}{U}\left(\frac U3\right)^{3/2}
<\frac{(ab)^{5/2}}{3\sqrt3}.$$

For the assertion when $a=2$, another row cannot contain a 1. Indeed the prime powers from $N$ avoided by either of these two rows would each divide 1 or 2. Their product would be at most 2, leaving the third row with three positive normalized parts of sum at most 2. Hence one can apply the theorem to each of the six parts.
