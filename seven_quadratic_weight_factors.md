# The weight factors $2P+1$ and $5P+1$ in the three residual shapes

Let $10P+7=DR$, with $P,R$ prime powers and $D\le10{,}000$. The [cofactor reduction](seven_quadratic_cofactors.md) leaves just three possible shapes of a failure at $N=(2P+1)(5P+1)$. Write $T=2P+1$, $W=5P+1$. In every remaining shape the first row $A=W(P,P,1)$ must avoid **every odd prime dividing $W$**. The companion $C$ must avoid 2 and **every prime dividing $T$**. Moreover $W\equiv2\pmod4$.

Indeed the coprime products $U,V$ of the full weight powers assigned to $A,C$ satisfy $T\mid V$. Write $V=TF$. Since $F\mid W$ and the first two parts of $C$ are $TRz,TRw$, we have $F\mid\gcd(z,w)$ by $\gcd(R,W)=1$. The three gcds are 2, 6, 2. Every remaining $D$ is divisible by 3, so $P\equiv2\pmod3$ and $3\nmid W$. The repeated parts in $A$ force 2 to be assigned to $C$. As $T$ is odd, $F=2$ and $U=W/2$ is odd.

**Prime-power factors of $W$.** Let $W=FS$ with $S$ an odd prime power. If $S>F$, then $F\equiv1,4\pmod5$. In particular $W/2$ cannot be a prime power. At any prime-power layer $E\mid S$ of the same base, with $E>F$, put $k=[F^2 S/E]_5$ in $\{1,2,3,4\}$. Removing the $S$ shift from $A$ gives $(FP,FP,F)$. Since $5FP=F^2S-F$, the first part has residue $(kE-F)/5$ modulo $E$. Thus its addition requires
$$
                        (2k-5)E+3F<0.
$$
At $E=S$ only $k=1$ or 4 occurs, and 4 is impossible. This also forces the cofactor of any odd prime power exceeding $\sqrt W$ to be at least 14: it is 2 modulo 4, cannot be divisible by 3, and must be 1 or 4 modulo 5.

**Prime-power factors of $T$.** Use the row data $(D,a,b,z,w)$ of the cofactor reduction, and write $m=az/D$, $n=bw/D$. Suppose $T=BS$ for a prime power $S$. All of $B,S$ are odd, and both are coprime to $D$. At a power $E\mid S$ of the same base assume $E>B\max(3,2m,2n)$. Set
$$
       u=[5B^2 (S/E)m]_a,\qquad v=[5B^2 (S/E)n]_b .
$$
Removing the $S$ shift from $C$, its first two parts become $BRz,BRw$ and the total becomes $BW$. Since $DR=5BS+2$ and $2W=5T-3$, their residues modulo $E$ are respectively
$$
                (uE+2Bm)/a,\qquad (vE+2Bn)/b,\qquad (E-3B)/2 .
$$
These are least residues under the displayed size condition. Consequently
$$
             2E(bu+av)+B(4bm+4an+3ab)\le abE .
$$
Already the strict inequality $2(bu+av)<ab$ at $E=S$ gives the following counts among residue classes $B\bmod ab$ coprime to $ab$. Any row avoiding the base of $S$ must pass this test.

| $D$ | $(a,b)$ | $(m,n)$ | Classes passing | Unit classes | Does $B=1$ pass? |
|---:|---|---|---:|---:|---|
| 4,587 | $(11,417)$ | $(4,14)$ | 368 | 2,760 | No |
| 7,383 | $(23,107)$ | $(2,4)$ | 108 | 2,332 | No |
| 9,453 | $(23,411)$ | $(2,4)$ | 304 | 5,984 | Yes |

Thus if $T$ itself is a prime power the first two cofactors are settled for all inputs, including small $P$ by the [finite theorem](seven_quadratic_family.md). The third case retains its other carry constraints.

Run `python3 check_seven_quadratic_weight_factors.py` for the exact class counts and direct residue comparisons at 109,324 layers of $W$ and 36 layers of $T$.
