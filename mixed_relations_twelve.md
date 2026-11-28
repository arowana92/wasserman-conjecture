# A part ratio and a three-part relation of height twelve

Let $A,B$ be positive trinomial rows of common weight $N$. Suppose two parts of $A$ have reduced ratio $a:b$ with $\max(a,b)\le12$. If the parts of $B$ satisfy a homogeneous relation with three nonzero integer coefficients of magnitude at most 12, **a prime dividing $N(N-1)$ divides both coefficients at every weight**.

There is a stronger bounded-weight conclusion. Normalize $B$ by the gcd of its three parts. Suppose each of the three pair gcds of this primitive row is at most 12. Then the same assertion holds for every $N\le10^{10000}$. Over all weights, the number of possible exceptions up to $X$ under these hypotheses is $O_\epsilon(X^\epsilon)$ for every $\epsilon>0$. Replacing 12 by 10 in both height bounds gives the assertion at every weight.

For a plane $r x+s y=t z$ in primitive form, the pair gcds divide the respective coefficients. If both rows satisfy full relations of height 12, every prime power of $N-1$ in a putative failure divides a coefficient of one of the relations. Thus $N-1\mid\operatorname{lcm}(1,\ldots,12)=27720$ and the [unrestricted finite proof](full_range_1b.md) applies. Together these statements settle the two-row case with two relations of height 12 whenever at least one uses all three parts. By the [endpoint reduction](all_weights_endpoint_reduction.md), the three-row conjecture also holds when each row satisfies such a bounded relation and one chosen relation uses all three parts. If only the endpoint supplies the full relation, its parts $(1,1,N-2)$ imply $N-2\le24$.

## The denominator of the pair scale

Write $M=N-1$, $A=(ag,bg,z)$ and $C=C(a,b)$ for the [ratio constant](part_ratio_predecessor.md). Assume the endpoint gcd is one. Allocate each full prime power of $N$ to a row avoiding its prime, and let the products assigned to $A,B$ be $U,V$, respectively. Then $N=UV$, $\gcd(U,V)=1$, $U\mid g,z$, and $B=V(x_1,x_2,x_3)$ with positive entries $x_i$ summing to $U$. In particular $U\ge3$.

Define $D=M/\gcd(M,g)$ and $m=M/D$. Since $m\mid g$ and $\gcd(m,U)=1$, there is an integer $j>0$ with $g=Ujm$, $\gcd(D,Uj)=1$, and
$$
                       (a+b)Uj<D.                                      \tag{1}
$$
Indeed $z=N-(a+b)g\ge U>1$. The two rows are exactly
$$
 A=(aUjm,\ bUjm,\ (D-(a+b)Uj)m+1),\qquad B=V(x_1,x_2,x_3),\qquad UV=Dm+1. \tag{2}
$$

Consider $q=p^e\parallel D$. If $A$ avoids $p$, then $p\nmid g$. Otherwise its first two parts would both vanish modulo the full power in $M$, forcing that power to divide $g$ and contradicting $p\mid D$. Thus $q$ is the full power in $M$, $q\mid C$, one of $a,b$ vanishes modulo $q$, and the other times $g$ is $1\bmod q$. If $B$ avoids $p$, then $q$ divides two of its primitive coordinates because its scale divides $N$. It also divides those two entries of $x$. Assign such a power to the index of the remaining entry.

This gives four pairwise coprime integers $c,d_1,d_2,d_3$ of product $D$, with $c\mid C$, $d_i$ dividing the pair gcd opposite index $i$ in primitive $B$, and $d_k d_l\mid x_i$ for $\{i,k,l\}=\{1,2,3\}$. For a full relation each $d_i$ divides the corresponding coefficient. The sum bound (1) makes the list of $D,U,j,x$ finite. Against every reduced pair $i:k$ of $x$ the [two-ratio bound](part_ratio_predecessor.md) further requires $\operatorname{lcm}(C,C(i,k))>(a+b)(i+k)$.

## Two exhaustive checks

The [plane and group enumeration](check_mixed_ratio_planes.py) chooses the four coprime factors, lists multiples of the indicated steps, and in the plane case solves the homogeneous equation for the third coordinate. It writes $m=m_0+Un$, $V=V_0+Dn$ for $n\ge0$, with $m_0$ the positive residue of $-D^{-1}\bmod U$. For a specified prime, the two rows have known affine parts $s_i n+b_i$. Reading a digit $k$ of $n$ changes state $(b_i)$ to $(\lfloor(b_i+s_i k)/p\rfloor)$ only if the three remainders sum to less than $p$. The finite graph accepts exactly when the zero state can be reached, allowing high zero digits. At $p\mid U$ only $A$ can be assigned and $V$ must be a unit. At $p\mid D$ an assigned $A$ must have unit $m$. The prime 2 always needs an avoiding row.

The [CRT enumeration](check_ratio_groups_cross.py) instead chooses $D$ first. It may allocate every full power of $D$ dividing the certified upper bound for $C$ to the ratio. It partitions the remaining powers among at most three groups of sizes at most the height. Modulo their product, $x_i$ is $U$ at its group and zero at the others. These CRT classes and the exact sum enumerate the companion. The program bounds ratio constants using factorial decrement and prime-power transport. It tests $m$ or $V$ directly in the appropriate prime-power class. At $p^h\parallel U$ it requires $v_p(Dm+1)=h$. At $p^e\parallel D$ the unit ratio term fixes $m\bmod p^e$ for $A$, or $V=U^{-1}\bmod p^e$ for $B$. At 2 outside $UD$, a nonunit $m$ gives the homogeneous slope row after deleting the isolated final digit one and shifting. A unit $m$ uses the affine row in (2). The first verifier uses the second to suggest a prime, then independently checks every carry layer through the prescribed suffix and constructs the graph of reachable higher-digit states. Its optional full-parameter mode instead uses the entire common progression in (2).

At height 12 the full relations give 10,332 distinct tuples $(a,b,D,U,x)$ and 10,402 tuples after choosing $j$. All are excluded. The unrestricted group enumeration gives 27,387 tuples before $j$. Only the following two families survive the local tests. Neither companion satisfies a full relation of height 12, as is verified by trying the three possible positions on the one-term side. At height 10 there are 623 group tuples and no remaining families.

| $a:b$ | $D$ | $U$ | $j$ | $x$ | Bases forced on $B$ |
|---|---:|---:|---:|---|---|
| $4:11$ | 13,860 | 673 | 1 | $(63,70,540)$ | 3, 5, 7 |
| $7:9$ | 2,772 | 167 | 1 | $(11,44,112)$ | 2, 11 |

## The two families through ten thousand decimal powers

In the first family $A$ is forced at $2^2$ and 11, so (2) and the unit-term congruences give $V\equiv399757\bmod609840$. In the second $A$ is forced at $3^2$ and 7 and gives $V\equiv141455\bmod174636$. At every listed base the three integers $x_iV$ must add without carries.

The [scale checker](check_ratio_group_scales.py) intersects the admissible high-digit intervals in the required bases and checks the resulting progression. Independently it moves upward through the progression. If at some power $q$ the numerator floor exceeds the sum of denominator floors, the failure persists until at least one denominator floor changes. Jumping to the first possible change preserves exhaustiveness. For $UV\le10^{10000}$ the graph method visits 244,816 and 847,430 intervals in the two cases. The factorial-floor method checks 161,615 and 273,353 intervals. Both find zero admissible positive multipliers.

For the asymptotic bound it suffices to use the first two forced bases in each row of the table. At each of these four base/row choices the total number $W$ of length-256 paths from all reachable states satisfies $W^2<p^{256}$ by integer arithmetic. The set of base-$p^{256}$ numbers formed from the block labels therefore has dimension strictly below $1/2$. The [block normalization argument](six_coefficient_row_relations.md#the-count-beyond-the-finite-range) and the uniform intersection estimate for multiplicatively independent bases give $O_\epsilon(X^\epsilon)$ common multipliers through $X$. Each multiplier specifies one weight in its family. [Corso and Shmerkin, Corollary 1.17](https://arxiv.org/pdf/2409.04608)

~~~sh
python3 check_ratio_groups_cross.py --height 10
python3 check_ratio_groups_cross.py --height 12
python3 check_mixed_ratio_planes.py --height 12
python3 check_mixed_ratio_planes.py --height 12 --groups
python3 check_ratio_group_scales.py --digits 10000
~~~
