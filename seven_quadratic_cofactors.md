# A prime-power factor of $10P+7$ with cofactor through 10,000

Let $P=p^k$ and $R=q^e$ be prime powers of positive exponent. Suppose $10P+7=DR$, where $1\le D\le10{,}000$. At weight $N=10P^2+7P+1$, Wasserman's full three-row conjecture holds for every such $P$ when $D$ is outside
$$
                            \{4587,7383,9453\}.
$$
In particular every cofactor through 4,586 works at unbounded inputs. At the three remaining cofactors, a failure would require the single row shape listed below for that cofactor. No square $P$ or square $R$ can fail. At 7,383 and 9,453, $R$ must be an ordinary prime ending in 9. At 4,587 its prime base must end in 1. The [binary continuation](seven_quadratic_suffixes.md) bounds the remaining inputs logarithmically and checks a larger finite range.

## A layer above a prime-power shift

The [finite quadratic theorem](seven_quadratic_family.md) covers $P\le10^8$ and gives the following structure for all larger inputs in a possible failure. With $\{t,u\}=\{2,5\}$, $T=1+tP$ and $W=1+uP$, one row $C$ is a multiple of $T$ and avoids every prime dividing $H=10P+7$. Its three parts divided by $T$ sum to $W$. Full prime powers of $H$ split into three nonempty groups, indexed by their unique residue-one positions in $C$.

The following argument works for any $D$ with $R>7D$. Our finite cofactor range satisfies this automatically when $P>10^8$. Permute $C$ so its first two parts are divisible by $R$. Write them as $TRz,TRw$. Since $tW=DR-u$, we have $t(z+w)<D$. Above the $R$ shift, the addition of $Tz,Tw,DP-Tz-Tw$ is carry-free in base $q$. Indeed the third actual part is $1+R(DP-Tz-Tw)$ and the two others vanish modulo $R$.

Let $E$ be a positive power of $q$ dividing $R$, with $E>7D$. Put $K$ equal to the least positive residue of $D^2 R/E$ modulo 10. Write $Z_z=[Tz]_E$ and $Z_w=[Tw]_E$ for the least nonnegative residues. Then
$$
  10[DP]_E=K E-7D,\qquad uZ_z+tz=k_z E,\qquad uZ_w+tw=k_w E
$$
for positive integers $k_z,k_w$. These follow from $10DP=D^2R-7D$ and $uT=DR-t$. The upper addition gives $Z_z+Z_w\le[DP]_E$. Thus
$$
             (t(k_z+k_w)-K)E\le t^2(z+w)-7D<-uD.
$$
In particular $K>2t$. At $E=R$, $K$ is 1 or 9 since $D$ is a unit modulo 10. This forces $t=2$, $u=5$, $K=9$ and $D\equiv3,7\pmod {10}$. The first row is $W(P,P,1)$, so $P$ must be odd and $C$ has odd coefficient. Now $0<2z,2w<D<E$ gives $1\le k_z,k_w\le5$. From $DR\equiv2\pmod5$ we get
$$
                k_z=1+((D(R/E)z-1)\bmod5),\qquad
                k_w=1+((D(R/E)w-1)\bmod5).
$$
At $E=R$ their sum is at most four.

## All shapes at a fixed cofactor

Let $a<b$ be the group products at the first two positions, so the group containing $q$ is at the third. Then $a,b>1$ are coprime unitary divisors of $D$ (a divisor is unitary if it is coprime to its complementary divisor). This also holds when $q\mid D$. We have $z=(D/a)m$, $w=(D/b)n$, with $m,n$ positive units modulo $a,b$, respectively. Since $TR$ is odd and $N$ even, the pair $z,w$ is even and has distinct 2-adic valuations by binary carry freedom. The residue-one conditions give
$$
             2Rz\equiv-5\pmod a,\qquad 2Rw\equiv-5\pmod b ,
$$
using $5T=DR-2$. All these conditions are necessary independently of the other factors of the weight.

Enumerating these shapes with $2(z+w)<D$ and $k_z+k_w\le4$ gives 31,921 possibilities. A finite obstruction at the low digits excludes 31,886. Write $R=r_0+10v$, $P=(Dr_0-7)/10+Dv$ with $Dr_0\equiv7\pmod {10}$. For each required prime $r$ (2 or a prime of $D$), lift the surviving residues of $v$ from $r^{j-1}$ to $r^j$. Discard a lift if its two actual part residues sum to more than $[N]_{r^j}$, or if the bin congruences fail at this layer. At 2 also impose odd $P$. An empty layer excludes every $v$. The longest certificate uses 781 binary digits.

Another 32 shapes are excluded at the leading digits, 31 at base 2 and one at base 3. Let $S=DTR/2$ and set
$$
    (\alpha,\beta,\gamma)=\left(2z/D,\ 2w/D,\ 1-2(z+w)/D\right).
$$
The actual parts divided by $S$ are $\alpha,\beta,\gamma-5/H$, with $H>10^9+7$. Normalize $S$ by a power of $r$ into $[1,r)$. After $j+1$ leading digits, a common scale $s\in[r^j,r^{j+1})$ gives the three truncated parts by taking the floors of these three proportions times $s$. For a possible word $(x,y,v)$ the intervals $[x/\alpha,(x+1)/\alpha)$ and $[y/\beta,(y+1)/\beta)$ must meet $[v/\gamma,(v+1)/(\gamma-5/(10^9+7)))$ and the scale interval. Enumerate next digit triples of sum less than $r$ and keep this necessary intersection. All comparisons use rational arithmetic. No high certificate needs more than 22 digits. Size and depth caps only stop unsuccessful searches for a certificate.

| $D$ | $a$ | $b$ | $z$ | $w$ |
|---:|---:|---:|---:|---:|
| 4,587 | 11 | 417 | 1,668 | 154 |
| 7,383 | 23 | 107 | 642 | 276 |
| 9,453 | 23 | 411 | 822 | 92 |

The independent verifier regenerates the parts in the opposite direction, recovering $a=D/\gcd(D,z)$ and $b=D/\gcd(D,w)$. It uses factorial-floor differences with $P$ as the binary input residue or $R$ as the odd-base residue. For the high certificates it subdivides the scale interval, directly tests all compatible integer triples and checks their digits. It independently recovers all 31,921 shapes and validates all 31,918 exclusions.

All three surviving cofactors are multiples of 3, so $P\equiv2\pmod3$. At $a$ the forced residues of $R$ are respectively 2 modulo 11, 7 modulo 23 and 10 modulo 23. None is a square. Thus $e$ is odd. For $D=4587$ we have $R\equiv1\pmod{10}$, forcing $q\equiv1\pmod{10}$ since odd powering permutes the units modulo 10. For the other two, $R\equiv9\pmod{10}$ forces $q\equiv9\pmod{10}$. If $e\ge3$, choose $E=R/q$. Here $D^2q\equiv1\pmod{10}$ contradicts $K>4$ whenever $E>7D$. This leaves only $q^{e-1}\le7D$; exhaustive integer enumeration gives respectively ten and twelve such odd powers, of which eight and ten have $P>10^8$. None satisfies both bin congruences there.

Reproduce the cofactor certificates from the repository root. The final prime-power enumeration is included in the [binary continuation verifier](check_seven_quadratic_suffixes.py).

```sh
mkdir -p .work
python3 check_seven_quadratic_cofactors.py 10000 --certificate .work/seven-quadratic-certificate.txt
python3 check_seven_quadratic_cofactors_cross.py 10000 .work/seven-quadratic-certificate.txt
```
