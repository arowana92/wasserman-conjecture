# Testing a quadratic factorization at the other predecessor prime

Let $N=1+DPR$, with $D\ge1$ and $P=p^k,R=q^e$ powers of distinct primes of positive exponent. Write $DR=MP+h$ where $M\ge1$, $0<h<P$. Throughout assume the [quadratic bounds](predecessor_digit_splits.md#a-second-range-for-quadratic-weights)
$$
                P>D,\qquad Dh\le P,\qquad D^2(h^2+4M)\le32P .
$$

* The full conjecture holds for $h=1,2,3,4,6$ at every cofactor. A split assignment is also impossible for $h=5,8$. At $h=7$ it requires $M=10$ and $D\equiv3,7\bmod10$.
* Write $\rho_M(v)$ for the representative of $v\bmod M$ in $\{1,\ldots,M\}$. If $\rho_M(D^2)\le4$, splitting is impossible, including when $h^2-4M$ is a square. If additionally $p\nmid D$, the full conjecture holds. At every $D$, grouping is impossible when $\gcd(D,h)\in\{1,2,3,4,6\}$ or when $R>Dh$ and $\rho_M(D^2)=1$.
* For example, $D=1001$, $P=1{,}566{,}731$ and $R=9391$ give $N=14{,}727{,}883{,}991{,}822$. Both powers are ordinary primes. Here $DR=6P+5$, the digit polynomial actually factors as $(1+2X)(1+3X)$, and $D^2\equiv1\bmod6$ settles the weight.

## An inequality at any layer of the second power

In a split assignment the quadratic theorem gives $t,u,i,j>0$, $t+u=h$, $tu=M$, $i\le j$ and $s=i+j\le t$. Its row taking $p$ is $(1+uP)(Pi,Pj,1+(t-s)P)$. Write $x/i=A/B$ in lowest terms for its first multiplier $x$. The exact ratio and the two scales give
$$
 x=Ai/B,\quad y=Aj/B,\qquad
 U=B(1+uP)/A,\quad V=A(1+tP)/B.
$$
Thus $B\mid i,j,1+tP$, $A\mid1+uP$, and $\gcd(B,t)=1$. The other row has two parts $VRz,VRw$ with $z,w>0$, and
$$
                  As(z+w)<BD,\qquad BuV=A(DR-t).
$$
If $i=j$ and $B$ is even, the two equal parts make the first coefficient even, yet $2\mid U$ is assigned to that row. This excludes that case.

Let $H=q^v$ divide $R$, with $v\ge1$. Assume $H>Dh$ and $2H\ge D(t^2-2h)$. Put $K=\rho_M(D^2 R/H)$. Then
$$
                    BK\ge2t,\qquad
       BK=2t\ \Longrightarrow\ t^2>s h.                 \tag{1}
$$
In particular $K\ge5$. Indeed $B\le s/2\le t/2$ would force equality everywhere at $K=4$, giving $s=t$ and violating the strict condition.

To see (1), put $v_z=[Vz]_H$ and $v_w=[Vw]_H$. The addition above the $R$ shift in the $q$-avoiding row is carry-free, hence $v_z+v_w\le[DP]_H$. The two positive numbers $Bu v_z+Atz$ and $Bu v_w+Atw$ are multiples of $H$. Also $M DP=D^2R-Dh$ and $0<Dh<H$, so $M[DP]_H=KH-Dh$. It follows that
$$
 (2t-BK)H\le At^2(z+w)-BDh
          < BD(t^2/s-h)\le D(t^2/2-h).
$$
The last inequality uses $2B\le s$ and $B\ge1$. The assumed bounds on $H$ imply both conclusions.

There are two exclusions before using the layer.

* If $t=u$, put $W=1+tP$. Here $N=W^2$. For each prime $\ell$ assigned to $U$, with $\ell^c\parallel W$, its full weight power $\ell^{2c}$ divides the part $WPi$. Since $\ell\nmid P$, we have $\ell^c\mid i$. Consequently $U\mid i^2$ whereas $U>2R$ and $R=(t^2P+2t)/D>t^2$.
* If $Bu=1$, then $V=A(DR-t)$. Since $At(z+w)<Dt/s\le Dt/2<R$, the residues of $Vz,Vw$ modulo $R$ sum to $2R-At(z+w)>R$. The bound $R>Dt/2$ follows from $u=1$, $M=t\ge2$ and $h^2+4t\ge17$ in the quadratic size inequality.

At $H=R$ the two layer hypotheses always hold in every remaining case. First, $R>DM(h^2+4M)/32$. With $u\ge2$ we have $M\ge4$ and $4(h^2+16)\ge32h$, hence $R>Dh$. With $u=1$, $t=2,3$ would force $B=1$ because $2B\le t$. At $t\ge4$ we have $t((t+1)^2+4t)\ge41t\ge32(t+1)$. For the second hypothesis, use $h\ge t+1$, $M\ge t$ and
$$
 t((t+1)^2+4t)-16t^2+32(t+1)=t((t-5)^2+8)+32>0.
$$
This proves $2R>D(t^2-2h)$.

The same test can use lower layers even when the upper layer passes. For example, in the surviving case at $h=7$ the only parameters are $t=2,u=5,i=j=B=1$. If $R/q^2>7D$ and $q\not\equiv1\bmod10$, splitting is impossible. The congruence $DR=10P+7$ excludes $q=2,5$. For $q\equiv3,7,9\bmod10$ one of $D^2q,D^2q^2$ is $1$ or $3\bmod10$, contradicting (1) at $R/q$ or $R/q^2$.

## Grouping and the small remainders

For grouping, the first-power argument gives $g=\gcd(D,h)$ and a part $c(N-1)/g$ with $0<c<g$ and $U\mid c$, $U\ge3$, $\gcd(U,D)=1$. This excludes $g\le3$. For $g=4,6$ the only possibility is $U=c=g-1$, a prime. Write $W=(N-1)/g$. Then $U\mid W+1$ since $U\mid N$. If $U^a\parallel W+1$, the residue of $UW$ modulo $U^{a+1}$ is $U^{a+1}-U$, and that of $W+1$ is at least $U^a$. The carry in their sum contradicts the assignment. If $p\nmid D$, then $g=\gcd(D,M)$ divides $\rho_M(D^2)$. This proves the second main statement for residues through four.

For $R>Dh$ and $\rho_M(D^2)=1$ there is also a direct grouped obstruction irrespective of $p$. The same row has the cut above $R$ given by $aP,(D-a)P$ for some $0<a<D$. Each positive number $M[aP]_R+ah$ and $M[(D-a)P]_R+(D-a)h$ is a multiple of $R$. Carry freedom would make their sum exactly $M[DP]_R+Dh=R$, which is impossible.

For each $h\le8$ enumerate $2\le t<h$, $u=h-t$, $1\le i\le j$, $s=i+j\le t$, and $B\mid\gcd(i,j)$ coprime to $t$. Exclude the equal-root, $Bu=1$, and equal-part even-denominator cases above. It suffices to test all residues $d\bmod M$ with $\gcd(d,M)\mid h$, since $D\mid MP+h$. Condition (1) leaves precisely two tuples $(t,u,i,j,B,d)$, both at seven:
$$
                   (2,5,1,1,1,3),\qquad(2,5,1,1,1,7).
$$
This enumeration uses only necessary conditions. Reproduce it with the Python checker or the Lean theorem `small_cases`. The independent layer check visits 388,705 cases, including 6,632 actual carry-free additions at the chosen layer. The grouped test then gives the asserted unconditional remainders, since every divisor of 4 or 6 belongs to the excluded list. [Python verification](check_predecessor_root_layers.py) · [Lean certificate](PredecessorDigitSplits.lean)
