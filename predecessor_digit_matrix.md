# Counting carry-free cuts at a full predecessor power

Fix a prime power $P=p^a\ge5$. Among any $H$ consecutive positive values of $W$ coprime to $P$ (or all eligible $W$ in an interval of length $H$), the number of possible exceptions to Wasserman's conjecture at
$$N=1+PW,\qquad \omega(W)\le2$$
is
$$
 O_P\!\left(H^{\log_P \lambda_{P-2}}\right),\qquad
 \lambda_m=\frac1{2\sin(\pi/(4m+2))}.
$$
The same bound holds in any $H$-term arithmetic progression of $W$ whose step is coprime to $P$, uniformly in the starting point and step. Fixing the full power allows one to count the actual allowed transitions between blocks of digits.

With at most two distinct predecessor primes the [endpoint theorem](all_weights_endpoint_reduction.md) already suffices. With exactly three, a failure supplies an integer $3\le x<P$, coprime to $P$, such that $xW\equiv1\bmod P$ and the cut
$$
                     xW,\qquad (P-x)W+1
$$
adds without a carry in base $p$. These are the exceptional part at $P$ and its complement in a row avoiding all three predecessor primes. A base-$p$ carry-free sum has no carry at any boundary of a block of $a$ digits, so it also adds without a carry in base $P$. At $a=1$ the two digit tests coincide. [Three-factor assignment](endpoint_factors.md)

## An exact matrix for the block test

Write $y=P-x$ and suppose $x,y>1$. The required low digit of $W$ is $r=x^{-1}\bmod P$, taken in $1,\ldots,P-2$. It produces output digits $1,0$. For every $L\ge0$ the number of residues modulo $P^{L+1}$ with this low digit and no carry in those $L+1$ places is
$$
              (A^L\mathbf1)_{r+1},
 \qquad A_{ij}=\mathbf1\{\,xj\bmod P\in I_i\,\},\qquad
 I_i=[x-\lfloor xi/P\rfloor,\ x-\lfloor xi/P\rfloor+i-1].
$$
The indices $i,j$ range from $1$ to $P-1$. The interval $I_i$ lies in $[1,P-1]$. The intervals grow by one element at each step: the floor increases by zero or one. To evaluate a matrix-vector product, reorder the vector by the residues $xj\bmod P$ and subtract two prefix sums for each interval. This uses $O(P)$ integer operations and $O(P)$ stored integers per digit.

For the proof, take a surviving input prefix $t<P^n$, $n\ge1$, and let $k=\lfloor Pt/P^n\rfloor$ be its top digit. Since $P$ divides $P^n$, $\lfloor(Pt+1)/P^n\rfloor=k$. The multiplication carries
$$a=\lfloor xt/P^n\rfloor,\qquad b=\lfloor(yt+1)/P^n\rfloor$$
satisfy $a+b=k$. Using $t\le(k+1)P^{n-1}-1$ gives the upper bounds $\lfloor x(k+1)/P\rfloor,\lfloor y(k+1)/P\rfloor$ on $a,b$. If $k=P-1$, the strict inequalities $a<x$, $b<y$ would give $a+b\le P-2$. Otherwise the two upper bounds sum to $k$, by $\gcd(x,P)=1$, and both are attained.

Put $i=k+1$ and append a digit $d=j-1$. The two next output digits sum to $k$ modulo $P$. Their sum actually equals $k$ exactly when $(xd+a)\bmod P\le k$. Substitution gives $xj\bmod P\in I_i$, since $x-a\ge1$ and $x-a+i-1\le P-1$ (use $a+b=i-1$ and $b<y$). The digit $d=P-1$ corresponds to residue zero and fails. The new state is $j$. Appending zero always works because $x\in I_i$; once the input ends its remaining output digits therefore also add freely.

## Uniform interval bound

State 1 transitions only to itself, and every row includes state 1. Delete its row and column. The remaining matrix $B$, with $m=P-2$ rows, has nested row supports of sizes $1,\ldots,m$. A column permutation makes it the lower triangular matrix of ones. Reversing columns instead arranges that triangular matrix as the symmetric matrix $S_{ij}=\mathbf1\{i+j\ge m+1\}$, $1\le i,j\le m$. Thus $\|B\|_2=\|S\|_2=\lambda_m$. Indeed $v_i=\sin(i\pi/(2m+1))$ satisfies $Sv=\lambda_m v$ by taking successive differences. Positivity of $v$ and symmetry give the norm. A path can remain in the $m$ states, or enter state 1 once and stay. Its count through $L$ steps is bounded by a constant depending only on $P$ times $\sum_{h=0}^L\lambda_m^h=O_P(\lambda_m^L)$.

When $y=1$, the initial digit must be $P-1$. While all digits read are $P-1$, the input prefix is $P^n-1$ and the multiplication carries are $(P-2,1)$. Once a digit $d<P-1$ is read, the carries become $(d,0)$, as direct division of $(P-1)d+P-2$ and $d+1$ shows. From carries $(k,0)$ exactly the digits $0\le d\le k$ are allowed and give carries $(d,0)$. Thus across the entire word the digits are weakly decreasing from least to most significant, with at most $\binom{L+P}{P-1}$ words of length $L+1$. This polynomial is also $O_P(\lambda_m^L)$.

For $H\ge2$ choose the least $P^{L+1}\ge H$ with $L\ge0$. The progression terms have distinct residues modulo $P^{L+1}$. Sum the bounds over the fewer than $P$ possible numerators. An interval containing $H$ consecutive units can be enclosed in an interval of length $O_p(H)$, covering that formulation as well.

The [Fermat-prime application](fermat_exponent_digits.md) evaluates this matrix at $P=65537$ and checks small prime and composite bases directly by quotient sums.
