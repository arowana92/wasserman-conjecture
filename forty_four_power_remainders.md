# Remainders of the two varying powers at cofactor 44

Let $E$ be the exponent of the unit group modulo $\operatorname{lcm}(1,\ldots,45^3)$. Let $Q,R$ be powers of distinct primes with positive exponents divisible by $E$. Wasserman's three-row conjecture holds at $N=1+44QR$ if **either power leaves a remainder between 1 and 6,586,647 inclusive on division by the other**. In particular, for arbitrary distinct bases, fixing one power allows an explicit infinite progression of exponents on the other: take multiples of both $E$ and the multiplicative order of its base modulo the fixed power. The theorem gives a partial result for the first cofactor surviving the [two-power reference test](two_power_seed.md).

## Forced rows and congruences

A possible failure has the endpoint and two other rows of coprime scales $U,V$ with $UV=N$. Assign every prime of $N(N-1)$ to one of those two rows avoiding it. For $a=4,b=11$ the [static cofactor theorem](static_power_transfer.md) excludes the grouped case: its scale would be $3$, which is not a unitary divisor of $45$. If a varying base is $2$ or $11$, absorbing the shared power leaves at most one static prime power and the [four-factor criterion](two_prime_powers.md) excludes a split.

Otherwise name the powers so that $Q=q^B$ is assigned to the row taking both $4$ and $11$, and write $R=r^C$. In the split of the static theorem, the primitive positive multipliers $X,Y$ in that row satisfy $4\mid X$, $11\mid Y$, $k(X+Y)=45$ and $2\le g(z+t)<k<4$. Thus $k=3$, $X=4$, $Y=11$, $g=z=t=1$. Up to permutation the actual rows are
$$
 A=(4UQ,11UQ,N-15UQ),\qquad C_0=(VR,VR,N-2VR).
$$
Positivity gives $V>15Q$, $U>2R$. Their product $44QR+1$ then gives $U<3R$, $V<22Q$.

Put $M=2^{16}11^4=959{,}512{,}576$. At each of the two layers $2^{16},11^4$, both varying powers are one and the weight is 45 modulo that layer. Both layers exceed $11\cdot45$. The residues of the first two parts of $A$ are positive and sum to at most 45 by carry freedom. Their ratio is exactly $4:11$ by cross multiplication at that layer. In base 2 the multiplier is odd and at most three, and multiplier one fails because 11 has bit one set while 45 does not. At 11 the residue of the first part must be one modulo 11, since it is a unit there and the three parts have residues $0,0,1$. This forces multiplier three there also. Cancelling 11 at the binary layer and 4 at the other, then using $UV=N$, gives
$$
 Q\equiv R\equiv1,\qquad U\equiv3,\qquad V\equiv15\pmod M. \tag{1}
$$

For $s=UQ$, the binomial $\binom{15s}{4s}$ is divisible by 15 because $4\binom{15s}{4s}=15\binom{15s-1}{4s-1}$. It divides the coefficient of $A$. In particular $q\notin\{3,5\}$. If $r\ne3$, then $3\mid N$ since $Q,R\equiv1\bmod3$, and its full weight power is assigned to $C_0$, so $3\nmid U$. If $r\ne5$, then both varying powers are one modulo 25 and $N\equiv45\bmod25$. Thus $5\mid N$ and $C_0$ must avoid 5. In that case $U\equiv3\bmod25$ is impossible: it would give $VR\equiv15\bmod25$, and the two displayed parts of $C_0$ would already carry.

## The two inequalities

Set $T=44R\bmod Q$ and $T'=44Q\bmod R$, taking least nonnegative residues. Above the $q^B$ shift the parts of $A$ give a carry-free addition of $4U,11U,44R-15U$ at $q$. If $11T<Q$, the positive residues $i,j$ of $4U,11U$ modulo $Q$ satisfy $i+j\le T$ and $11i=4j$: both cross products are less than $Q$ and congruent. Thus $i=4k$, $j=11k$, $15k\le T$, and $U\equiv k\bmod Q$. Write $Vk=1+HQ$. Since $15Q<V<22Q$, we have $15k\le H<22k$. By (1), $H\equiv15k-1\bmod M$. Consequently
$$
 11T\ge Q\quad\text{or}\quad 15M\le7T . \tag{2}
$$

Above the $r^C$ shift the parts of $C_0$ give $V,V,44Q-2V$ without a carry at $r$. Let $k=V\bmod R>0$. Then $2k\le T'$. Write $Uk=1+HR$. Here $2k\le H<3k$, and $H\equiv3k-1\bmod M$. For $k<M$ this forces $H=3k-1$, so
$$
 k(3R-U)=R-1. \tag{3}
$$
If $r\ne3$, then $3R-U$ is a 3-adic unit and $3^{10}\mid R-1$, hence $3^{10}\mid k$. If $r\ne5$, then $R\equiv1\bmod 5^7$ while $U\not\equiv3\bmod25$, hence $5^6\mid k$. Both depths follow from $3^{10},5^7\le45^3$.

When $r=3$ or $5$, use the carry-free pair $V,V$ one more time at the low layer $r^d$, with $3^d=3^{10}$ or $5^d=5^7$ respectively. This layer divides $R$ and $44Q\equiv44\bmod r^d$. The two copies of $k\bmod r^d$ must fit digit by digit into $44$; also $r\nmid k$. At 3 this leaves $1,4$, and at 5 it leaves $1,2,6,7$. Combining these classes with the divisibilities from (3) gives the table. For a step $s$ and prescribed residue $c$ modulo $w$, the least positive candidate here is $s(c s^{-1}\bmod w)$.

| Base of $R$ | Constraints when $k<M$ | Lower bound on $k$ |
| --- | --- | ---: |
| $3$ | $5^6\mid k$, $k\bmod3^{10}\in\{1,4\}$ | $144{,}906{,}250$ |
| $5$ | $3^{10}\mid k$, $k\bmod5^7\in\{1,2,6,7\}$ | $207{,}734{,}382$ |
| Other | $3^{10}5^6\mid k$ | $922{,}640{,}625$ |

All three lower bounds also hold when $k\ge M$. In particular $T'\ge289{,}812{,}500$. Since $E$ is divisible by the order $2^{14}$ of 5 modulo $2^{16}$, both powers exceed $11\cdot289{,}812{,}500$. If either ordinary remainder is at most 6,586,647, multiplying it by 44 gives at most 289,812,468. In the $T'$ orientation this contradicts its lower bound. In the $T$ orientation it gives $11T<Q$ and $7T<15M$, contradicting (2). The two orientations cover either choice of assignment.

## Replacing the large exponent period

The remainder theorem also holds when each varying exponent is a positive multiple of **$6{,}706{,}706{,}688{,}000{,}000$**. In fact it suffices that for each of $2^{16},11^4,3^{10},5^7$ each varying power coprime to that modulus is one modulo it. The displayed period is the least common multiple of the four unit-group exponents.

Here the grouped and split arguments in the [product-congruence theorem](static_product_congruence.md) apply just using the first two layers, both larger than $44\cdot45$. At a shared base use the shifted cut in its polynomial-period proof. A grouped assignment must have $U\mid c$ or $U\mid d$, where both $c,d$ are 3; hence $U=3$. If 3 is a varying base, $3\nmid N$. Otherwise both powers are one modulo 9, so $9\mid N$ and 3 is not a unitary divisor. This rules out grouping. In a split the ratio projection gives the same rows as above. Both varying powers are one modulo $M$ and at least one of $3^{10},5^7$ as well, so each exceeds $11\cdot289{,}812{,}500$. This also guarantees $r^d\mid R$ at the layer used when $r=3,5$. The rest of the remainder proof uses exactly the four given layers.

In particular there are families with both varying exponents equal to one at this surviving cofactor. Put $J=2^{16}11^4 3^{10}5^7$. Choose any prime $R\equiv1\bmod J$ and then any prime $Q\equiv1\bmod JR$. Both exist in infinite arithmetic progressions, and the remainder on division of $Q$ by $R$ is one.

Run **lean FortyFourPowerRemainders.lean** with Lean 4.36. The formal lemmas check the integer intervals in (2) and (3) and the three Chinese remainder bounds. [Lean verification](FortyFourPowerRemainders.lean)
