# The product of the six nonendpoint parts

Let $A,C$ be the two rows paired with the mandatory endpoint in a counterexample of weight $N$, put $M=N-1$, and write $Z_A,Z_C$ for the products of the three parts in each row. Then
$$
L=\frac{Z_A Z_C}{N^3 M^2}
$$
is an even integer and **$M<9^{40}L^{25}$**. Thus the product of the six parts has an unconditional lower bound of order $N^{5+1/25}$. The exponent 25 in the bound on $M$ is optimal if one uses only the lowest residues at primes of $N(N-1)$. An explicit family below satisfies those residues, while its two coefficients share a growing divisor of $N$.

An exhaustive search at every weight excludes $L\le32$, so $L\ge34$ at an exception and $L\ge40$ at an even-weight exception. The [limiting ratios](six_part_limits.md) further restrict families near the exponent 25.

At most $O_\epsilon(B^{5+\epsilon})$ weights can admit a counterexample with $L\le B$, independently of the size of the weight. The case of two rows each having two small parts also has the [stronger bounds in their two sums](two_row_sums.md).

## A rational product close to an integer

The argument works quantitatively for any two rows all of whose parts are at least two. Let $D=\gcd(N,K_A,K_C)$ and $J=\gcd(M,K_A,K_C)$, put $G=DJ$, and in this generality set $L=D^3J^2 Z_AZ_C/(N^3M^2)$. We will prove
$$
(M/J)G^9<9^{40}L^{25}. \tag{1}
$$

Use the [retained power assignment](endpoint_pair_bound.md) with $DUV=N$ and $JST=M$. For each power in $S$, put it into bin $s_i$ if its unique nonzero residue in $A$ is at position $i$. Write the three products as $s_1s_2s_3=S$ and similarly $t_1t_2t_3=T$ on the other row. Empty bins have product one. For positive integer multipliers $f_i,g_j$ we have
$$
 A_i=Uf_i S/s_i,\qquad C_j=Vg_j T/t_j,\qquad
 L=\prod_i f_i\prod_j g_j.
$$
Write $\alpha_i=f_i/s_i$ and $\beta_j=g_j/t_j$. Then
$$
\left(\sum_i\alpha_i\right)\left(\sum_j\beta_j\right)=G N/M>G.
$$
Every proper nonempty rectangle of this product, obtained by restricting either sum or both to a subset of its terms, is strictly below $G$. Indeed the partial product is $G\,(\sum A_i)(\sum C_j)/(NM)$, and a proper subset has sum at most $N-2$. For a counterexample the endpoint assignment gives $U,V\ge3$, so the minimum-part hypothesis is automatic.

Order both lists by decreasing fraction. Start with the two largest terms and include their bins, whose product $K_0=s_1t_1$ is less than $9 f_1g_1/G$. Put $z=s_1/f_1$. The largest cell is below $G$ and exceeds $G/9$, giving
$$g_1/t_1<Gz,\qquad 1/f_1\le z<9g_1/G.$$
At each of the next four steps extend the side whose next term is the larger proportion of its own first term. Let that proportion be $R$, and let $K$ be the product of bins already included. The known rectangle has denominator dividing $K$, so its deficit below $G$ is at least $1/K$. Each missing cell is at most $\alpha_1\beta_1 R$, and together they exceed the deficit. Thus $1/K<9\alpha_1\beta_1 R$.

If the next term is $f/s$, then $s<9K f g_1/t_1$. If it is $g/t$, then $t<9K g/z$. The new bin product is at most $9 f Gz K^2$ or $9gK^2/z$, respectively. Iterating squares the four factors with respective weights $8,4,2,1$ and gives $K_0$ weight 16. Let $a$ be the sum of the two weights on the first side and $b$ the sum of the other two. We have $a+b=15$ and $|a-b|\le9$. The factors involving $G,z$ are $G^{a-16} z^{a-b}$. For $a>b$ use $z<9g_1/G$, and for $a<b$ use $1/z\le f_1$. The resulting power of $G$ is at most $\min(a,b)-16\le-9$. Each multiplier now has exponent at most $16+9=25$, and the constant is at most $9^{31+9}$. The final bin product is $ST=M/J$, proving (1).

Whenever $G$ is odd, one of the rows avoids 2. Name it $A$ and assign 2 to it. For even $N$, $U$ contains the whole power of 2 in $N$, and the bins $s_i$ are odd. The three valuations of the entries of $A$ are distinct by binary carry freedom. Subtracting $v_2(U)$ leaves three distinct nonnegative valuations for the $f_i$, so $8\mid\prod f_i$. For odd $N$, two entries of $A$ are multiples of the full power of 2 in $S$. Their valuations differ, but the bin factor $S/s_i$ has the same 2-adic valuation at both positions. Thus at least one multiplier is even. This proves $2\mid L$ at every possible exception.

## Why the exponent is attained by the residue conditions

Take a positive integer $n\equiv2\bmod42$ and set
$$
 a=n+1,\quad b=na+1,\quad c=nab+1,\quad S=abc,\quad V=ab+ac+bc,
$$
$$
 d=V+1,\quad e=Vd-1,\quad T=de,\quad U=nde+d+e.
$$
We have $S-nV=1$ by the three successive definitions. Also $VU-ST=V(d+e)-de=1$. Use the bins $(a,b,c)$ with multipliers $(1,1,1)$, and $(1,d,e)$ with multipliers $(n,1,1)$. The actual rows at $N=UV$ are
$$A=U(bc,ac,ab),\qquad C=V(nde,e,d).$$
They have the stipulated residues $0,0,1$ at each predecessor power assigned to their side. The following coprimality checks ensure that all assigned powers are full.

The bins within either side are pairwise coprime. Modulo $S$, we have $nV\equiv-1$ and hence $nd\equiv n-1$. Thus $\gcd(d,S)=\gcd(n-1,42)=1$ by evaluating $S$ as a polynomial at 1. Also $n^2e\equiv1-n-n^2\bmod S$. Modulo $n^2+n-1$, the bins $a,b,c$ reduce to $n+1,2,3$. The modulus is coprime to $n+1$, is odd, and is never divisible by 3. Thus $\gcd(e,S)=1$. Finally $U\equiv-n\bmod V$ and $V\equiv3\bmod n$, giving $\gcd(U,V)=1$.

Here $L=n$, and $M=abcde$ is a monic polynomial of degree $1+2+4+6+12=25$ in $n$. This proves optimality for the residue conditions. The coefficient gcd in the actual construction is divisible by $V$. Indeed the first row has gcd $U$ and $N/U=V$ divides its coefficient by the one-slot factorial identity. For the second row, $\gcd(e,d)=1$ and $e+d=V(V+2)$. The binomial coefficient of the corresponding scaled pair is divisible by $e+d$, by multiplying it by $e/(e+d)$. Thus these assignments are genuine controls for the divisibility tests at the lowest layers.

## Counting and reconstructing assignments

For this section $D=J=1$. Fix six positive numerators with product $L$. Order each side by increasing bin denominator, allowing all orders of the numerators. With prefixes of lengths $i,j$ in the two lists let $S',T'$ be their bin products, $X/S',Y/T'$ the partial sums, $K=S'T'$ and $\Delta=K-XY>0$. Write $H=(\sum f_i)(\sum g_j)\le9L$. The first two bins have $K_0=s_1t_1<H$, and $\Delta_0=K_0-f_1g_1\ge1$.

Of the remaining bins, consider the smallest relative to its side's first bin. If it is $v$ on the first side, every missing cell is at most its numerator divided by $vt_1$. They total more than $\Delta/K$, so $v<HK/(\Delta t_1)$. The corresponding bound on the other side is $v<HK/(\Delta s_1)$. After adding one bin its new deficit numerator is $\Delta_1=\Delta_0 v-c>0$, where $c=f_2s_1g_1$ or $g_2t_1f_1$. In either case $c\le HK_0$.

Only two further bins need to be chosen after $s_1,t_1$. For fixed $s_1,t_1$ and fixed sides for these next two bins, their number of choices is at most
$$
 \sum_v \frac{HK_0 v}{\Delta_1}
 =\frac{HK_0}{\Delta_0}\sum_v (1+c/\Delta_1).
$$
There are at most $HK_0/\Delta_0$ possible $v$. Their positive $\Delta_1$ are spaced by at least $\Delta_0$ and are polynomially bounded in $L$. Thus $\sum 1/\Delta_1\ll1+\log(2L)/\Delta_0$. For fixed $\Delta_0$ the product $K_0$ is specified, allowing $O_\epsilon(L^\epsilon)$ factor pairs. Summing over $\Delta_0\le9L$ and using $HK_0\le81L^2$ gives $O_\epsilon(L^{4+\epsilon})$ choices for the four bins.

For the last two bins $v,w$, expansion of the product equation gives
$$
(\Delta v-c)(\Delta w-e)=ce+\Delta(h-1).
$$
Here $c,e$ are the two positive coefficients of $w,v$ in the completed numerator product after removing its known term $XYvw$. If both missing bins are on the first side, these are $f_i S'Y,f_j S'Y$ and $h=0$, and likewise on the second side. If the missing bins are on opposite sides they are $f_i S'Y,g_j T'X$ and $h=f_i g_j K$. Both factors on the left are positive because adding either one bin still gives a proper rectangle. A positive solution thus factors a fixed positive integer polynomially bounded in $L$. The divisor bound leaves $O_\epsilon(L^\epsilon)$ choices. The number of numerator sextuples of product at most $B$ is $O_\epsilon(B^{1+\epsilon})$, proving the stated count.

The first enumerator constructs these four-bin prefixes and recovers the last two bins by factoring. An independent method builds all six bins in increasing order across both rows and recovers the final bin from a linear equation. Its direct mode scans all possible fifth-bin values. For longer scans in the default mode it evaluates the cleared bilinear polynomial at three integer points, factors its elimination integer, and tests the resulting fifth-bin candidates. For its next-bin bound, if the remaining numerators on the sides sum to $A_0,B_0$, and the next bin is $v$, then
$$\Delta v^2 < (XT'B_0+YS'A_0)v+K A_0 B_0.$$
This follows by bounding every missing fraction by its numerator divided by $v$. Both searches retain only coprime scales and disjoint prime-power assignments.

Through $L=32$ the search finds 202,517 assignments up to row exchange and permutation of entries within a row. All but six give a carry in both rows at 2. The remaining six give a carry in both rows at 3, with $3\mid N(N-1)$ in each case. The largest reconstructed weight is 70,123,398,015,687,058,597,795,876,772,907,086,317. There is no cutoff on $N$ in either enumeration.

| $L$ | Assignments | $L$ | Assignments |
|---:|---:|---:|---:|
| 2 | 158 | 18 | 15,243 |
| 4 | 382 | 20 | 16,385 |
| 6 | 1,974 | 22 | 10,955 |
| 8 | 3,317 | 24 | 30,349 |
| 10 | 2,601 | 26 | 9,723 |
| 12 | 6,828 | 28 | 13,597 |
| 14 | 5,746 | 30 | 39,150 |
| 16 | 8,048 | 32 | 38,061 |

Consequently any two positive trinomial rows whose product of all six parts is at most $32N^3(N-1)^2$ share a prime divisor with the endpoint. Indeed, otherwise they would give a counterexample with the endpoint, forcing one of the enumerated even products. In the enumeration the small prime list $\{2,3\}$ applies to the assignments passing the assumed lowest-residue conditions.

With python-flint and SymPy installed, run the commands below. The first search tests the two actual row valuations by Legendre’s formula. The second checks the common primes by digit sums. Both agree on the counts and maximum weight at every product. The direct fifth-bin scan also reproduces the first 5,831 assignments. The last command checks the unrestricted inequality (1) on 278,490 pairs of exact trinomial coefficients at weights through 55, including 13,879 pairs with odd endpoint gcd. It also verifies the symbolic recurrences and residue exclusions in the [limiting-ratio proof](six_part_limits.md).

```sh
python3 check_six_part_product.py 32
python3 check_six_part_product_cross.py 32 --jobs 4
python3 check_six_part_product_cross.py 8 --direct
python3 check_product_limits.py
```
