# Power indices at the first thousand cofactors and at 2278

Let $q$ be prime and $e\ge1$. **Wasserman's three-row conjecture holds at every $1+dq^e$, simultaneously for $1\le d\le1000$, if $e$ is a multiple of any of the following 19 integers.**

$$
\begin{gathered}
2,\ 21,\ 51,\ 595,\ 2523,\ 17017,\ 24871,\ 40375,\ 78625,\ 83657,\\
272935,\ 382109,\ 531505,\ 2347835,\ 8921773,\ 42632447,\ 67148081,\\
166654111,\ 2573653511 .
\end{gathered}
$$

These divisibility classes have natural density $299026867069663408869/561053121385651819925$, about $53.297425\%$. Eleven of the indices extend the [earlier list of eight](common_power_exponents.md). For instance every $1+dq^{40375t}$ works for $t\ge1$ and $d\le1000$.

**At $d=2278$ the corresponding sufficient list is $3,8,49,67,289,2756$.** Its union has natural density $287087939/653714243$, about $43.9164\%$. In particular every $1+2278q^{3t}$ works, with no restrictions on $q$ or $t\ge1$.

There is a converse about the reach of low-digit tests. For each $e$ outside the first union, some $d\le1000$ admits the construction below. For each $e$ outside the second union it applies with $d=2278$. Given any prime cutoff $B$ and digit depth $k$, a reduced arithmetic progression contains arbitrarily large prime bases $q$ with specified positive rows $E,A,C$ of weight $1+dq^e$. For each prime through $B$, at least one of the three rows adds without a carry through its first $k$ digits. **The full conjecture holds at these weights.** Indeed the bases can be chosen to satisfy the conjecture at every $1+d'q^e$ with $d'\le2278$ simultaneously. With $q$ varying, both $B,k$ may grow as $\lfloor(\log q)/160\rfloor$, even while the certified cofactor range grows to $q^{1/240}$. The case of arbitrary fixed cofactor bounds is discussed in the [finite divisibility principle](cofactor_index_principle.md).

## Reducing the coefficients to local power tests

In a failure write $P=q^e$ and take a nonendpoint row $A$ avoiding $q$. The [endpoint reduction](all_weights_endpoint_reduction.md) and [primitive denominator argument](shifted_directions.md) give a positive primitive triple $(x,y,z)$ of sum $m\mid d$ and rows
$$A=(xW,yW,zW+1),\qquad C=V(f,g,h),\qquad W=(d/m)P,\qquad UV=1+mW.$$
Here $U=f+g+h\ge3$ is a unitary divisor of the weight dividing $x,y$. Each prime of $U$ is assigned to $A$ and each prime of $V$ to $C$. At every prime of $m$, one of the rows must have the full predecessor residues $0,0,1$ in some order.

For $d\le1000$ the [existing enumeration](common_power_exponents.md) settles all even exponents and all binary bases. At the remaining odd scales it gives $m=d$, $W=P$, and the fourteen row shapes in the first part of the table below. The [small-base theorem](prime_base_density.md) already includes all primes through 1000. For $d=2278=2\cdot17\cdot67$ enumerate the primitive shapes at every divisor of 2278. Denominators 1 and 2 cannot be sums of three positive slopes. At $17,34,67,134$ there are no preliminary shapes. There are 261 at 1139 and 893 at 2278, as checked both by iterating $U$ first and by iterating the two unshifted slopes first. The [finite carry tests for a shape and parity](shifted_directions.md) leave only $m=2278$, $W=P$ odd, $U=7$ and $C=V(1,2,4)$. The two possibilities are
$$
\begin{array}{c|c|c}
(x,y,z)&P\bmod15946&\text{odd bases required of }A\\\hline
(833,1407,38)&5665&7,17,67\\
(938,1309,31)&6687&7,17,67 .
\end{array}
$$
Thus $q=2,7,17,67$ are already excluded from any failure in this classification by parity or the unit residue of $P$. In general the shape search gives a progression $P\equiv R\pmod M$ and its forced odd first-row bases $r$. In this computation every such $r$ divides $M$ and $R$ is a unit at $r$.

For an odd prime $r$ and a unit target $t\bmod r$, let $o$ be its order modulo $r$. At precision $r^j$, the group of units is cyclic of order $(r-1)r^{j-1}$. An $e$th power can have residue $t$ only if $o\mid (r-1)/\gcd(e,r-1)$. Once that holds, a unit $w\equiv t\pmod r$ is an $e$th power to every precision exactly when
$$ w^{r-1}\equiv1\pmod{r^{1+v_r(e)}}. $$
This follows by taking the torsion and principal-unit factors of the cyclic groups, or directly by raising a generator to the $e$th power. There is a unique lift of $t$ satisfying $w^{r-1}=1$ modulo each $r^j$. Starting at $t$, raising the representative to the $r$th power gives this lift one digit further, by binomial expansion.

Follow those lifts at $t=R\bmod r$. At precision $w=r^j$ also impose the $r$-part of the progression modulo $M$ and the full $r$-part of $U$ in $1+dP$. Try all additional low digits needed to determine those conditions. For each permitted prefix of $P$, read the multiplication carries of $A$ and ask whether a finite high suffix can reach a carry-free termination. The [carry transition](higher_successors.md#testing-whether-a-digit-prefix-can-finish) is an exact finite-state test. If the first impossible precision is $r^j$, then $r^{j-1}\mid e$ excludes the row. At the preceding precision the verifier gives an ordinary positive scalar with an entirely carry-free row at $r$. In the two cases that continue indefinitely, the scalar 1 itself works. Together with the torsion order test, these witnesses prove local attainability for every index passing the recorded conditions.

| $d$ | First slope $x$ | Divisors of $e$ sufficient to exclude this row |
|---:|---:|:---|
| 380 | 130 | 2, 3, 5, 13, 19 |
| 589 | 133 | 2, 3, 5, 7, 19, 31 |
| 629 | 259 | 2, 3, 17, 1369 |
| 646 | 133 | 2, 3, 17, 343 |
| 703 | 133 | 2, 3, 7, 19, 37 |
| 710 | 175 | 2, 3, 5, 7, 71 |
| 780 | 224 | 2, 3, 7, 125, 169 |
| 850 | 175 | 2, 3, 7, 17, 25 |
| 874 | 209 | 2, 3, 5, 11, 19, 23 |
| 946 | 301 | 2, 3, 5, 7, 11, 43 |
| 950 | 275 | 2, 3, 5, 11, 19 |
| 962 | 143 | 2, 3, 5, 11, 37, 169 |
| 986 | 119 | 2, 7, 17, 29 |
| 986 | 203 | 2, 7, 17, 841 |
| 2278 | 833 | 3, 8, 49, 67, 289, 2756 |
| 2278 | 938 | 2, 3, 11, 17, 49, 67 |

For example, on the first 2278 row the only completable residues of $P$ above the initial digit modulo $7^3,17^3$ are 226 and 905, respectively. The lifts of order prime to 7 and 17 there are 324 and 2928, respectively, and both fail. At $67^2$ the lift 3789 fails. For the second row the corresponding failures occur at $7^3,17^2,67^2$. These give some of the prime-power indices in the table.

The extra divisor $2756=52\cdot53$ uses the weight. If $q\ne53$ and $2756\mid e$ then $P\equiv1\pmod{53^2}$. Consequently $53\mid V$ and $V/53\equiv43/7\equiv44\pmod{53}$. The column of $V,2V,4V$ at place 53 has digits $44,35,17$. Their sum forces a carry at a prime assigned to $C$. At $q=53$ both first-row slope triples carry in base 53, precluding a $q$-avoiding first row. The implication for arbitrary exponents divisible by 2756 follows by taking $x=q^{e/2756}$ in the [Lean verification](Cofactor2278Power.lean).

To combine rows, start with the divisor 1. At each new row take all least common multiples of a retained divisor with one divisor from that row and remove redundant multiples. This gives precisely the divisibility-minimal indices excluding every row processed. Inclusion–exclusion over the resulting indices gives the two densities above. A separate check counts residues by their gcd with the full period.

## Constructing the prime bases outside the lists

Choose a row whose divisors all miss $e$. At 2278 the first row suffices for every $e$ outside its list. At each forced odd first-row prime $r$ take the finite scalar witness just described. The power criterion guarantees roots modulo $r^k$ at any requested depth. Prescribe one such root for $q$; these conditions at different primes are independent.

The following choices handle the other assigned bases. Rational scalar targets with denominator prime to the base have well-defined residues at every precision.

- Through $d=1000$ the remaining indices $e$ are odd, so every odd binary unit is locally an $e$th power. At an odd $d$ choose $V=2$. At an even $d$ other than 780 choose $V=1$. In each case $P=(UV-1)/d$ is a binary unit, the listed companion adds without carries there, and the local part of the progression for $P$ holds. At $d=780$ choose $P=147\equiv3\pmod4$ locally at 2, making the first row carry-free. The companion must also avoid 3. Use $V=1$, giving $P=1/130$ locally at 3. Since $\gcd(e,6)=1$, every ternary unit is an $e$th power.
- At $d=2278$ choose $V=65$ at 2. Then $P=227/1139\equiv1\pmod{16}$ and $65,130,260$ have disjoint binary ones. For even $e$, the $e$th powers among binary units are precisely the units congruent to 1 modulo $2^{2+v_2(e)}$. This follows from the cyclic subgroup generated by 5 modulo every power of two and the sign factor. Since $8\nmid e$, the proposed scalar is attainable.
- At other odd primes $r\le B$ take $q\equiv1\pmod r$ if $r\nmid d+1$. This keeps $r$ out of the weight. At any other $r\mid d+1$ take a unit whose $e$th power is different from 1 whenever $r-1\nmid e$. That is always possible for the surviving odd indices at $d\le1000$. At $d=2278$ we have $d+1=43\cdot53$ and $3\nmid e$, so 43 causes no restriction. If $52\mid e$, the weight is necessarily divisible by 53. In this case $53\nmid e$ because $2756\nmid e$. Choose $P=(7\cdot53-1)/2278$ locally at 53. It is 1 modulo 53 and thus an $e$th power to every precision, and $C$ at $V=53$ has no carry in that base.

These choices include the modulus of the original row progression. They prescribe units for $q$ in pairwise coprime prime-power moduli. The Chinese remainder theorem and the prime number theorem in arithmetic progressions give a positive relative density of prime bases in the resulting class. Take $q>\max(d,B)$, so $A$ also avoids $q$ by its normalized digits. At each prime through $B$ that divides the endpoint, the assigned row now has the required carry-free low digits. [Prime number theorem in arithmetic progressions](https://terrytao.wordpress.com/2009/09/24/the-prime-number-theorem-in-arithmetic-progressions-and-dueling-conspiracies/)

One can place a decisive carry immediately after a prescribed long prefix. Choose an odd $r\mid d$ forced on $A$. When the cofactor has two rows use 29 at 986 and 17 at 2278, where the two initial residues differ. For its positive scalar witness $z_0$, take $j\ge k$ large enough for $r^j$ to exceed all three actual witness parts and $j>v_r(e)$. Impose $P\equiv z_0+r^j\pmod{r^{j+1}}$. This remains in the image of the $e$th powers. The row has no carry through $r^j$. In the next column its digits are the three slopes modulo $r$. Two are nonzero and opposite because $r\mid d$ and exactly one unshifted slope is an $r$-unit. Hence its first carry is at $r^{j+1}$. The other row cannot avoid this forced prime and the endpoint is divisible by it. This also excludes the chosen shape for all triples at the weight. Where there is a second shape, its different initial residue already excludes it.

For simultaneous certification throughout $d'\le D'$, refine the binary prefix of $P$ by the [binary prefix extension theorem](cofactor_binary_prefix_extensions.md) with its parameters $B=J=1$. Before refinement choose a precision at least $2+v_2(e)$, so every extension of this unit prefix remains an $e$th power locally. The resulting root class for $q$ retains all previously prescribed low bits of $P$ and is compatible with the odd-prime conditions. This proves the fixed-cutoff assertion with $D'=2278$.

For the growing assertion set $B=k=n$ and $D'=2^n$, taking $n$ sufficiently large for the fixed index. Start with at least $n$ binary digits and append at most $8+20n$ by the same theorem. All powers for the forced odd primes have precision at most $n+O(1)$, including the extra digit specifying the carry. The product of the bases needing high precision, including 2 and the optional 53, is less than $10^6$ in every case. At other primes through $n$ only one digit is prescribed. Thus the modulus $m$ of the root class has
$$\log m\le(2+\log 10^6+20\log 2)n+O(1)<30n$$
for large $n$, using $\sum_{r\le n}\log r<2n$ over primes. Its least prime $q$ is greater than $n$ since the class is a unit at all primes through $n$. Xylouris's bound $q\ll m^{5.2}$ gives $\log q<160n$ for large $n$. Hence both cutoffs exceed $\lfloor(\log q)/160\rfloor$, while $2^n\ge q^{1/240}$. Taking $n$ to infinity gives infinitely many distinct prime bases. [Xylouris, *On Linnik's constant*](https://arxiv.org/abs/0906.2749)

## An explicit exponent-52 example

The following $q$ is prime.

```text
227091210478544415535791867586401686132438401519077924179522111756785060404092160168052056211
```

At $N=1+2278q^{52}$ the conjecture holds for every triple of rows. For the displayed 833-shape with its endpoint and companion, the common primes through 67 are exactly $2,7,17,53,67$. For each of them at least one row clears the first twelve layers. The first row's carry at 17 occurs at layer 13, and the alternative 938-shape has the wrong residue modulo 17. The verifier constructs $q$ from the local roots, proves its primality with FLINT, and evaluates the actual factorial floors for this 15,965-bit weight.

Run the following with Python 3.9 or later, python-flint, and Lean 4.36.0. The second program checks the first impossible power residues by enumerating their roots and searching forward from the multiplication carries. It verifies the positive witnesses by direct addition. In total the two programs check 48 odd-base conditions.

```sh
python3 check_cofactor_power_indices.py
python3 check_cofactor_power_indices_cross.py
lean Cofactor2278Power.lean
```
